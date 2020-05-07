/*
 * File:   rootNEST.cpp
 * Author: mszydagis
 *
 * Created on July 7, 2018, 11:14 PM
 */

#include <cmath>
#include <fstream>
#include <iostream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TF1.h"
#include "TFeldmanCousins.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"

#include <string.h>
#include <vector>

#include "analysis.hh"
#include "TestSpectra.hh"  // contains the WIMP, earth, and escape velocities

// non-input parameters hardcoded in
#define CL 0.90  // confidence level
#define VSTEP 1e-3  // step size in keV for convolving WIMP recoil E with eff

using namespace std;

vector<vector<double> > GetBand_Gaussian(vector<vector<double> > signals);
vector<vector<double> > GetBand(vector<double> S1s, vector<double> S2s,
                                bool resol);
TRandom3 r;  // r.SetSeed(10);

double band[NUMBINS_MAX][17], band2[NUMBINS_MAX][17];
/* 
band[bin][0] = s1c_center
band[bin][1] = s1c_actual_mean
band[bin][2] = log(S2/S1) mean
band[bin][3] = log(S2/S1) stddev
band[bin][4] = log(S2/S1) alpha
band[bin][5] = log(S2/S1) mean error
band[bin][6] = log(S2/S1) stddev error
band[bin][7] = log(S2/S1) alpha error
band[bin][8] = log(S2/S1) chi2/ndf (NEST)
band[bin][9] = log(S2/S1) xi
band[bin][10] = log(S2/S1) xi error
band[bin][11] = log(S2/S1) omega
band[bin][12] = log(S2/S1) omega error
band[bin][13] = log(S2/S1) covariance xi-omega
band[bin][14] = log(S2/S1) covariance xi-alpha
band[bin][15] = log(S2/S1) covariance omega-alpha
band[bin][16] = log(S2/S1) chi2/ndf (ROOT)
*/

void GetFile(char* fileName);
vector<vector<double> > outputs, inputs;

TestSpectra myTestSpectra;
double expectedUlFc(double mub, TFeldmanCousins fc);

int modPoisRnd(double poisMean, double preFactor);  // not used currently, but
                                                    // included jic needed for
                                                    // future use, for
                                                    // stretching and squeezing
                                                    // Poissons
int modBinom(int nTot, double prob,
             double preFactor);  // just like above, but for binomial
double EstimateSkew ( double mean, double sigma, vector<double> data );
double owens_t(double h, double a);

bool loop = false;
double g1x = 1.0, g2x = 1.0;  // for looping over small changes in g1 and g2

int main(int argc, char** argv) {
  bool leak, ERis2nd;
  double NRbandX[numBins], NRbandY[numBins],
      numSigma[numBins], leakage[numBins], discrim[numBins],
      errorBars[numBins][2];
  int i = 0; if ( loopNEST ) verbosity = false;

  if (argc < 2) {
    cout << endl
         << "This program takes 1 (or 2) inputs." << endl
         << endl;
    cout << "One input means you are just doing a band of Gaussians." << endl
         << endl;
    cout << "Two inputs means you're doing both ER and NR and calculating "
            "leakage and discrimination." << endl;
    cout << "If you write NR first, you're seeking characterization of the "
            "non-Gaussian asymmetry of NR band." << endl;
    cout << "If you write ER first, you're finding non-Gaussian leakage of ER "
            "into NR band." << endl
         << endl;
    return 1;
  } else if (argc == 2)
    leak = false;
  else
    leak = true;

if ( mode == 2 ) {

  if (argc < 3) {
    cerr << "Enter 0 for SI, 1 for SD-n, and 2 for SD-p" << endl;
    return 1;
  }

  FILE* ifp = fopen(argv[1], "r");
  int ch, nLines = 0;
  while (EOF != (ch = getc(ifp))) {
    if ('\n' == ch) nLines++;
  }
  double energy[nLines], efficiency[nLines];
  rewind(ifp);
  for (i = 0; i < nLines; i++) {
    fscanf(ifp, "%lf %lf", &energy[i], &efficiency[i]);
    if ( efficiency[i] > 1.02 ) //2% margin of error is allowed
      { cerr << "eff should be frac not %" << endl; return 1; }
    if ( efficiency[i] < -.02 ) //allowing for digitization err
      { cerr << "eff must not be negative" << endl; return 1; }
  }
  fclose(ifp);
  TGraph* gr1 = new TGraph(nLines, energy, efficiency);
  TF1* fitf =
      new TF1("FracEffVkeVEnergy",
              "10.^(2.-[0]*exp(-[1]*x^[2])-[3]*exp(-[4]*x^[5]))/100.",
              energy[0], energy[nLines - 1]);  // eqn inspired by Alex Murphy
  fitf->SetParameters(10., 2., 1., 20., 1e4, -2.5);
  gr1->Fit(fitf, "q");  // remove the quiet option q if you want to see all of
                        // the parameters' values for Fractional Efficiency
                        // versus Energy in keV. Prints to screen
  double aa = fitf->GetParameter(0);
  double bb = fitf->GetParameter(1);
  double cc = fitf->GetParameter(2);
  double dd = fitf->GetParameter(3);
  double ee = fitf->GetParameter(4);
  double ff = fitf->GetParameter(5);
  delete gr1;
  delete fitf;  // it is always important to clear the memory in ROOT

  // Get input parameters for sensitivity or limit calculation
  double time, fidMass, loE, hiE, xEff, NRacc, numBGeventsExp = 0.,
                                               numBGeventsObs;
  cout << "Target Mass (kilograms): ";
  cin >> fidMass;
  cout << "Run Time (provide days): ";
  cin >> time;
  cout << "Multiplicative factor on the efficiency of NR event detection from "
          "file (usually ~>0.9): ";
  cin >> xEff;  // unitless fraction of one
  cout << "Acceptance for NR events post electron recoil background "
          "discrimination (usually ~ 50%): ";
  cin >> NRacc;  // unitless fraction of one
  cout << "Number of BG events observed: ";
  cin >> numBGeventsObs;
  if (numBGeventsObs > 0.) {
    cout << "Number of BG events expected: ";
    cin >> numBGeventsExp;  // for FC stats
  }
  cout << "Minimum energy (keV) for detection: ";
  cin >> loE;
  cout << "Maximum energy (keV) for detection: ";
  cin >> hiE;
  // Make sure inputs were valid.
  if (cin.fail() || fidMass <= 0. || time <= 0. || xEff <= 0. || NRacc <= 0. ||
      loE < 0. || hiE <= 0. || numBGeventsExp < 0. || numBGeventsObs < 0.) {
    cerr << endl
         << "Input error. Make sure all inputs were numbers (most also "
            "positive or at least 0)" << endl;
    return 1;
  }
  if ( xEff > 1. || NRacc > 1. ) {
    cerr << endl
	 << "You entered an efficiency or acceptance for NR greater than 100%"
	 << endl;
  return 1;
  }
  
  double Ul,
      v;  // Start of code block to evaluate upper limit on # of events to use
  if (numBGeventsExp == 0.) {
    for (v = 0.; v < 1e3; v += VSTEP) {
      double sum = 0.0;
      for (i = 0; i < (numBGeventsObs + 1.); i++)
        sum += exp(-v) * pow(v, i) / TMath::Factorial(i);   // using Poisson
                                                            // statistics as the
                                                            // probability
                                                            // distribution
      if (sum <= (1. - CL)) break;
    }
    Ul = 0.5 * (2. * v - VSTEP);  // should result in 2.3 events as the UL (0
                                  // BG)
  } else {
    TFeldmanCousins fc(CL);
    if (numBGeventsExp != numBGeventsObs ||
        numBGeventsObs == int(numBGeventsObs))
      Ul = fc.CalculateUpperLimit(numBGeventsObs, numBGeventsExp);
    else
      Ul = expectedUlFc(numBGeventsExp, fc);
    double powCon =
        fc.CalculateUpperLimit(0., 0.);  // can't do better than 2.44 ever!
    if (Ul < powCon) Ul = powCon;
  }

  const int masses = NUMBINS_MAX;
  double massMax = 1e5;
  double mass[masses] = {
      // list of masses to run
      2.0,  2.1,  2.2,  2.3,  2.4,  2.5,    2.6,  2.7,  2.8,  2.9,  3.0,  3.1,
      3.2,  3.3,  3.4,  3.5,  3.6,  3.7,    3.8,  3.9,  4.0,  4.1,  4.2,  4.3,
      4.4,  4.5,  4.6,  4.7,  4.8,  4.9,    5.0,  5.2,  5.4,  5.6,  5.8,  6.0,
      6.5,  7.0,  7.5,  8.0,  8.5,  9.0,    9.5,  10,   11,   12,   13,   14,
      15,   16,   17,   18,   19,   20,     22,   24,   26,   28,   30,   32,
      34,   36,   38,   40,   42,   44,     46,   48,   50,   55,   60,   65,
      70,   75,   80,   85,   90,   95,     100,  110,  120,  130,  140,  150,
      160,  170,  180,  190,  200,  220,    240,  260,  280,  300,  320,  340,
      360,  380,  400,  420,  440,  460,    480,  500,  550,  600,  650,  700,
      750,  800,  850,  900,  950,  1000,   1100, 1200, 1300, 1400, 1500, 1600,
      1700, 1800, 1900, 2000, 2200, 2400,   2600, 2800, 3000, 3200, 3400, 3600,
      3800, 4000, 4200, 4400, 4600, 4800,   5000, 5500, 6000, 6500, 7000, 7500,
      8000, 8500, 9000, 9500, 1E+4, massMax};  // in GeV
  double sigAboveThr[masses], xSect[masses];  // arrays for the fraction of WIMP
                                              // signal events above threshold
                                              // and for the cross-sections
  cout << "\nWIMP Mass [GeV/c^2]\tCross Section [cm^2]" << endl;

  i = 0;
  while (mass[i] < massMax) {  // Iterate across each sample wimp Mass
    sigAboveThr[i] = 0.;
    for (double j = VSTEP; j < hiE;
         j += VSTEP) {  // Iterate across energies within each sample wimp mass
      double eff = pow(10., 2. - aa * exp(-bb * pow(j, cc)) -
                                dd * exp(-ee * pow(j, ff))) /
                   100.;
      if ( eff > 1. || eff < 0. )
	{ cerr << "Eff cannot be greater than 100% or <0%" << endl; return 1; }
      if (j > loE)
        sigAboveThr[i] += VSTEP * myTestSpectra.WIMP_dRate(j, mass[i]) * eff * xEff *
                          NRacc;  // integrating (Riemann, left sum)
                                  // mass-dependent differential rate with
                                  // effxacc and step size
    }
    xSect[i] = 1e-36 * Ul / (sigAboveThr[i] * fidMass *
                             time);  // derive the cross-section based upon the
                                     // desired upper limit, calculated above
                                     // (aka "Ul"). Relative to 1pb
    if (atof(argv[2]) > 0.) {  // uses arXiv:1602.03489 (LUX Run03 SD) which in
                               // turn uses Klos et al.
      xSect[i] *= (1. - 2.9592 / pow(mass[i], 1.) + 352.53 / pow(mass[i], 2.) -
                   14808. / pow(mass[i], 3.) + 2.7368e+5 / pow(mass[i], 4.) -
                   2.5154e+6 / pow(mass[i], 5.) + 1.2301e+7 / pow(mass[i], 6.) -
                   3.0769e+7 / pow(mass[i], 7.) + 3.1065e+7 / pow(mass[i], 8.));
      if (atof(argv[2]) == 1.) xSect[i] *= 1.66e5;
      if (atof(argv[2]) == 2.) xSect[i] *= 5.64e6;
    }  // end spin-dep. code block
    i++;
    if (xSect[i - 1] < DBL_MAX && xSect[i - 1] > 0. &&
        !std::isnan(xSect[i - 1]))  // Print the results, skipping any weirdness
                                    // (low WIMP masses prone)
      cout << mass[i - 1] << "\t\t\t" << xSect[i - 1] << endl;  // final answer
  }
  int iMax = i;

  printf("{[");  // DM tools format (i.e. Matlab)
  for (i = 0; i < (iMax - 1); i++) {
    if (xSect[i] < DBL_MAX && xSect[i] > 0. && !std::isnan(xSect[i]))
      printf("%.1f %e; ", mass[i], xSect[i]);
  }
  printf("%.1f %e", mass[iMax - 1], xSect[iMax - 1]);  // last bin special
  printf("]}\n");

  return 0;

}

if ( mode == 1 ) {

  if (numBins == 1) {
    GetFile(argv[1]);
    return 0;
  }
  
  int DoF = numBins - freeParam;
  FILE* ifp = fopen(argv[2], "r");
  for (i = 0; i < numBins; i++) {
    fscanf(ifp, "%lf %lf %lf %lf %lf %lf", &band2[i][0], &band2[i][1],
           &band2[i][2], &band2[i][5], &band2[i][3], &band2[i][6]);
  }
  fclose(ifp);
  // comment in the next 3 lines for g1 and g2 variation loops. Use 2> /dev/null
  // when running to suppress empty data warnings from low g1 sending things out
  // of bounds
  for (g1x = 0.90; g1x <= 1.10; g1x += 0.01) {
    for (g2x = 0.90; g2x <= 1.10; g2x += 0.01) {
      if (loop)
        printf("%.2f\t%.2f\t", g1x, g2x);
      else {
        g1x = 1.;
        g2x = 1.;
      }
      GetFile(argv[1]);
      double error, chi2[4] = {0., 0., 0., 0.};
      for (i = 0; i < numBins; i++) {
	if ( fabs(band[i][0]-band2[i][0]) > 0.05 ) {
	  cerr << "Binning doesn't match for GoF calculation. Go to analysis.hh and adjust minS1, maxS1, numBins" << endl;
	  return 1;
	}
        error = sqrt(pow(band[i][5], 2.) + pow(band2[i][5], 2.));
        chi2[0] += pow((band2[i][2] - band[i][2]) / error, 2.);
        error = sqrt(pow(band[i][6], 2.) + pow(band2[i][6], 2.));
        chi2[1] += pow((band2[i][3] - band[i][3]) / error, 2.);
	chi2[2] += 100. * ( band2[i][2] - band[i][2] ) / band2[i][2];
	chi2[3] += 100. * ( band2[i][3] - band[i][3] ) / band2[i][3];
      }
      chi2[0] /= double(DoF - 1); chi2[2] /= numBins;
      chi2[1] /= double(DoF - 1); chi2[3] /= numBins;
      //cout.precision(3);
      //if ( fabs(chi2[0]) > 10. ) chi2[0] = 999.;
      //if ( fabs(chi2[1]) > 10. ) chi2[1] = 999.;
      if ( loopNEST ) { //abbreviated #only version
	cout << chi2[0] << "\t" << chi2[1] << "\t" << 0.5*(chi2[0]+chi2[1]) << "\t" << pow(chi2[0]*chi2[1],0.5) << "\t" << chi2[2] << "\t" << chi2[3] << "\t"
	     << band[0][2] << "\t" << band[0][3] << endl;
      }
      else {
	cout << "The reduced CHI^2 = " << chi2[0] << " for mean, and " << chi2[1] << " for width. ";
	cout << "Arithmetic average= " << (chi2[0]+chi2[1])/2. << " and geo. mean " << sqrt(chi2[0]*chi2[1]) << " (mean+width and mean*width)"<<endl;
	cout << "%-errors of the form (1/N)*sum{(data-NEST)/data} for mean and width are " << chi2[2] << " and " << chi2[3] << " (averages)" << endl;
      }
      if (!loop) break;
    }
    if (!loop) break;
  }  // double curly bracket goes with g1x and g2x loops above
  return 0;

}

  if (leak) {
    cout << endl << "Calculating band for second dataset" << endl;
    GetFile(argv[2]);
    for (i = 0; i < numBins; i++) {
      band2[i][0] = band[i][0];
      band2[i][1] = band[i][1];
      band2[i][2] = band[i][2];
      band2[i][3] = band[i][3];
      band2[i][4] = band[i][4];
      band2[i][5] = band[i][5];
      band2[i][6] = band[i][6];
      band2[i][7] = band[i][7];
      band2[i][8] = band[i][8];
    }
  }

  cout << endl << "Calculating band for first dataset" << endl;
  GetFile(argv[1]);
  if (leak) {
    double finalSums[3] = {0., 0., 0.};
    if (band2[0][2] > band[0][2]) {
      ERis2nd = true;
      fprintf(stderr,
              "Bin Center\tBin Actual\t#StdDev's\tLeak Frac Gaus\t+ error  - "
              "error\tDiscrim[%%]\tLower ''Half''\t+ error  - "
              "error\tUpperHf[%%]\tignore column!\n");
    } else {
      ERis2nd = false;
      fprintf(stderr,
              "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t- => Gaus more\n");
      fprintf(stderr,
              "Bin Center\tBin Actual\t#StdDev's\tLeak Frac Gaus\t+ error  - "
              "error\tDiscrim[%%]\tLeak Frac Real\t+ error  - "
              "error\tDiscrim[%%]\tReal-Gaus Leak\n");
    }
    for (i = 0; i < numBins; i++) {

      if (!skewness) {
        if (ERis2nd) {
          numSigma[i] = (band2[i][2] - band[i][2]) / band2[i][3];
          errorBars[i][0] =
            (band2[i][2] - band2[i][5] - band[i][2] - band[i][5]) /
            (band2[i][3] + band2[i][6]);  // upper (more leakage)
          errorBars[i][1] =
            (band2[i][2] + band2[i][5] - band[i][2] + band[i][5]) /
            (band2[i][3] - band2[i][6]);  // lower (less leakage)
          NRbandX[i] = band[i][0];
          NRbandY[i] = band[i][2];
        }
        else {
          numSigma[i] = (band[i][2] - band2[i][2]) / band[i][3];
          errorBars[i][0] =
            (band[i][2] - band[i][5] - band2[i][2] - band2[i][5]) /
            (band[i][3] + band[i][6]);  // upper (more leakage)
          errorBars[i][1] =
            (band[i][2] + band[i][5] - band2[i][2] + band2[i][5]) /
            (band[i][3] - band[i][6]);  // lower (less leakage)
          NRbandX[i] = band2[i][0];
          NRbandY[i] = band2[i][2];
        }
        leakage[i] = (1. - erf(numSigma[i] / sqrt(2.))) / 2.;
        errorBars[i][0] = (1. - erf(errorBars[i][0] / sqrt(2.))) / 2.;
        errorBars[i][1] = (1. - erf(errorBars[i][1] / sqrt(2.))) / 2.;
      }

      if (skewness) {
        if (ERis2nd) {
          leakage[i] = 0.5 + 0.5 * erf((band[i][2] - band2[i][9]) / band2[i][11] / sqrt(2.)) -
            2. * owens_t((band[i][2] - band2[i][9]) / band2[i][11], band2[i][4]);

          double pdf_at_NR = (1./(band2[i][11]*sqrt(2.*TMath::Pi()))) * exp(-0.5*pow(band[i][2]-band2[i][9],2.)/pow(band2[i][11],2.)) *
            (1.+erf(band2[i][4]*(band[i][2]-band2[i][9])/(band2[i][11]*sqrt(2.))));
          double dlkg_dx = pdf_at_NR;
          double dlkg_dxi = -1. * pdf_at_NR;
          double dlkg_domega = -1. * (band[i][2] - band2[i][9]) / band2[i][11] * pdf_at_NR;
          double dlkg_dalpha = -1. / TMath::Pi() / (1. + pow(band2[i][4],2.)) * 
            exp(-1. * (1. + pow(band2[i][4],2.)) * pow((band[i][2] - band2[i][9]),2.) / 2. / pow(band2[i][11],2.));

          errorBars[i][0] = sqrt(0. +
            pow(dlkg_dx,2.) * pow(band[i][5],2.) + pow(dlkg_dxi,2.) * pow(band2[i][10],2.) +
            pow(dlkg_domega,2.) * pow(band2[i][12],2.) + pow(dlkg_dalpha,2.) * pow(band2[i][7],2.) + 
            2 * dlkg_dxi * dlkg_domega * band2[i][13] + 2 * dlkg_domega * dlkg_dalpha * band2[i][15] +
            2 * dlkg_dxi * dlkg_dalpha * band2[i][14]);
          errorBars[i][1] = errorBars[i][0];
        }
        else {
          leakage[i] = 0.5 + 0.5 * erf((band2[i][2] - band[i][9]) / band[i][11] / sqrt(2.)) -
            2. * owens_t((band2[i][2] - band[i][9]) / band[i][11], band[i][4]);

          double pdf_at_NR = (1./(band[i][11]*sqrt(2.*TMath::Pi()))) * exp(-0.5*pow(band2[i][2]-band[i][9],2.)/pow(band[i][11],2.)) *
            (1.+erf(band[i][4]*(band2[i][2]-band[i][9])/(band[i][11]*sqrt(2.))));
          double dlkg_dx = pdf_at_NR;
          double dlkg_dxi = -1. * pdf_at_NR;
          double dlkg_domega = -1. * (band2[i][2] - band[i][9]) / band[i][11] * pdf_at_NR;
          double dlkg_dalpha = -1. / TMath::Pi() / (1. + pow(band[i][4],2.)) * 
            exp(-1. * (1. + pow(band[i][4],2.)) * pow((band2[i][2] - band[i][9]),2.) / 2. / pow(band[i][11],2.));

          errorBars[i][0] = sqrt(0. +
            pow(dlkg_dx,2.) * pow(band2[i][5],2.) + pow(dlkg_dxi,2.) * pow(band[i][10],2.) +
            pow(dlkg_domega,2.) * pow(band[i][12],2.) + pow(dlkg_dalpha,2.) * pow(band[i][7],2.) +               
            2 * dlkg_dxi * dlkg_domega * band[i][13] + 2 * dlkg_domega * dlkg_dalpha * band[i][15] +
            2 * dlkg_dxi * dlkg_dalpha * band[i][14]);
          errorBars[i][1] = errorBars[i][0];
        }
      }

      finalSums[2] += leakage[i];
      discrim[i] = 1. - leakage[i];
    }
    TGraph* gr1 = new TGraph(numBins, NRbandX, NRbandY);
    TF1* fitf =
        new TF1("NRbandGCentroid", "[0]/(x+[1])+[2]*x+[3]", minS1, maxS1);
    fitf->SetParameters(10., 15., -1.5e-3, 2.);
    gr1->Fit(fitf, "rq", "", minS1, maxS1);  // remove the q if you want to see
                                             // the Band Fit Parameters on
                                             // screen
    double chi2 = fitf->GetChisquare() / (double)fitf->GetNDF();
    if (chi2 > 1.5) {
      cerr << "WARNING: Poor fit to NR Gaussian band centroids i.e. means of "
              "log(S2) or log(S2/S1) histograms in S1 bins. Investigate please!"
           << endl;
      fitf = new TF1("NRbandGCentroid", "[0]+([1]-[0])/(1+(x/[2])^[3])", minS1,
                     maxS1);
      fitf->SetParameters(0.5, 3., 1e2, 0.4);
      gr1->Fit(fitf, "rq", "", minS1, maxS1);
      chi2 = fitf->GetChisquare() / (double)fitf->GetNDF();
      if (chi2 > 2.) {
        cerr << "ERROR: Even the backup plan to use sigmoid failed as well!"
             << endl;
        return 1;
      }
    }
    long below[NUMBINS_MAX] = {0};
    double NRbandGCentroid, leakTotal, poisErr[2];
    for (i = 0; i < numBins; i++) {
      below[i] = 0;
      for (long j = 0; j < inputs[i].size(); j++) {
        NRbandGCentroid =
            fitf->GetParameter(0) / (inputs[i][j] + fitf->GetParameter(1)) +
            fitf->GetParameter(2) * inputs[i][j] +
            fitf->GetParameter(3);  // use Woods function
        // NRbandGCentroid = NRbandY[i]; // use the center of the bin instead of
        // the fit, to compare to past data that did not use a smoothing spline
        // NRbandGCentroid =
        // fitf->GetParameter(0)/(NRbandX[i]+fitf->GetParameter(1))+fitf->GetParameter(2)*NRbandX[i]+fitf->GetParameter(3);
        // // compromise
        if (outputs[i][j] < NRbandGCentroid) below[i]++;
      }
      leakTotal = double(below[i]) / (double)outputs[i].size();
      poisErr[0] =
          (double(below[i]) + sqrt(double(below[i]))) /
          (double(outputs[i].size()) - sqrt(double(outputs[i].size())));
      poisErr[1] =
          (double(below[i]) - sqrt(double(below[i]))) /
          (double(outputs[i].size()) + sqrt(double(outputs[i].size())));
      fprintf(
          stderr,
          "%.2f\t\t%.6f\t%.6f\t%e\t%.2e %.2e\t%.6f\t%e\t%.2e %.2e\t%.6f\t%e\n",
          0.5 * (band[i][0] + band2[i][0]), 0.5 * (band[i][1] + band2[i][1]),
          numSigma[i], leakage[i], errorBars[i][0] - leakage[i],
          leakage[i] - errorBars[i][1], discrim[i] * 100., leakTotal,
          poisErr[0] - leakTotal, leakTotal - poisErr[1],
          (1. - leakTotal) * 100., leakTotal - leakage[i]);
      finalSums[0] += (double)below[i];
      finalSums[1] += (double)outputs[i].size();
    }
    fprintf(stderr,
            "OVERALL DISCRIMINATION or ACCEPTANCE between min and maxS1 = "
            "%.12f%%, total: Gaussian + non-Gaussian ('anomalous'). Leakage "
            "Fraction = %.12e\n",
            (1. - finalSums[0] / finalSums[1]) * 100.,
            finalSums[0] / finalSums[1]);  // Dividing by total for average
    double HighValue, LowValue;
    LowValue = (finalSums[0] + sqrt(finalSums[0])) /
               (finalSums[1] - sqrt(finalSums[1]));
    HighValue = (finalSums[0] - sqrt(finalSums[0])) /
                (finalSums[1] + sqrt(finalSums[1]));
    fprintf(stderr,
            "highest possible discrimination value (plus  1-sigma) is %.12f%%, "
            "with corresponding leakage of %.12e\n",
            (1. - HighValue) * 100., HighValue);
    fprintf(stderr,
            " lowest possible discrimination value (minus 1-sigma) is %.12f%%, "
            "with corresponding leakage of %.12e\n",
            (1. - LowValue) * 100., LowValue);
    fprintf(stderr,
            "OVERALL DISCRIMINATION or ACCEPTANCE between min and maxS1 = "
            "%.12f%%, Gaussian.                                     Leakage "
            "Fraction = %.12e\n",
            (1. - finalSums[2] / numBins) * 100., finalSums[2] / numBins);
    delete gr1;
    delete fitf;
  }

  return 0;
}

void GetFile(char* fileName) {
  FILE* ifp = fopen(fileName, "r");
  double a, b, c, d, e, f, g, h, i, j, k, l, m, n;
  double eMin = 1e100, eMax = -1e100;
  int ch, nLines = 0, o;
  vector<double> E_keV, electricField, tDrift_us, X_mm, Y_mm, Z_mm, Nph, Ne,
      S1cor_phe, S2cor_phe, S1raw_phe, S1cor_phd, S1cor_spike, Ne_Extr,
      S2raw_phe, S2cor_phd;
  
  if ( verbosity ) {
    while (EOF != (ch = getc(ifp))) {
      if ('\n' == ch && nLines)
	break;
      else
	nLines = 0;
      if (']' == ch && nLines == 0) nLines++; }
  }
  
  while (1) {
    fscanf(ifp,
           "%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
           &a, &b, &c, &d, &e, &f, &g, &h, &i, &j, &k, &l, &m, &n);
    if (feof(ifp)) break;
    // fprintf(stderr,"%.6f\t%.6f\t%.6f\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
    E_keV.push_back(a);
    if (a < eMin) eMin = a;
    if (a > eMax) eMax = a;
    electricField.push_back(b);
    tDrift_us.push_back(c);
    X_mm.push_back(d);
    Y_mm.push_back(e);
    Z_mm.push_back(f);
    Nph.push_back(g);
    Ne.push_back(h);
    S1raw_phe.push_back(i);
    if (usePD <= 0 && fabs(j * 1.2) > minS1 && (j * 1.2) < maxS1) {
      S1cor_phe.push_back(j * 1.2);  // here and down below for S2: FIX THIS by
                                     // getting P_dphe (USUALLY ~0.2) from
                                     // detector class
      S1cor_phd.push_back(0);
      S1cor_spike.push_back(0);
    } else if (usePD == 1 && fabs(j) > minS1 && j < maxS1) {
      S1cor_phe.push_back(0);
      S1cor_phd.push_back(j);
      S1cor_spike.push_back(0);
    } else if (usePD >= 2 && fabs(k) > minS1 && k < maxS1) {
      S1cor_phe.push_back(0);
      S1cor_phd.push_back(0);
      S1cor_spike.push_back(k);
    } else {
      S1cor_phe.push_back(-999.);
      S1cor_phd.push_back(-999.);
      S1cor_spike.push_back(-999.);
    }
    Ne_Extr.push_back(l);
    S2raw_phe.push_back(m);
    if (usePD <= 0 && fabs(n * 1.2) > minS2 && (n * 1.2) < maxS2) {
      S2cor_phe.push_back(n * 1.2);
      S2cor_phd.push_back(0);
    } else if (usePD >= 1 && fabs(n) > minS2 && n < maxS2) {
      S2cor_phe.push_back(0);
      S2cor_phd.push_back(n);
    } else {
      S2cor_phe.push_back(-999.);
      S2cor_phd.push_back(-999.);
    }
  }
  fclose(ifp);

  if (numBins == 1) {
    TH1F* HistogramArray = new TH1F[3];
    double minimum, maximum, average;
    unsigned long int i, numPts = E_keV.size();
    double holder[numPts];

    printf(
        "S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc Mean\t\tEc "
        "Res[%%]\n");

    for (int j = 0; j < 3; j++) {
      switch (j) {
        case 0:
          for (i = 0; i < numPts; i++) {
            if (usePD == 0)
              holder[i] = S1cor_phe[i];
            else if (usePD == 1)
              holder[i] = S1cor_phd[i];
            else
              holder[i] = S1cor_spike[i];
          }
          if (holder[numPts / 2] > 1.0) {
            minimum = holder[numPts / 2] / 2.;
            maximum = holder[numPts / 2] * 2.;
          } else {
            minimum = minS1;
            maximum = maxS1;
          }
          break;
        case 1:
          for (i = 0; i < numPts; i++) {
            if (usePD == 0)
              holder[i] = S2cor_phe[i];
            else
              holder[i] = S2cor_phd[i];
          }
          if (holder[numPts / 2] > 20.) {
            minimum = holder[numPts / 2] / 2.;
            maximum = holder[numPts / 2] * 2.;
          } else {
            minimum = minS2;
            maximum = maxS2;
          }
          break;
        default:
          minimum = eMin;
          maximum = eMax;
          for (i = 0; i < numPts; i++) holder[i] = E_keV[i];
          break;
      }
      TString HistName;
      HistName.Form("%d", j);
      HistogramArray[j].SetName(HistName.Data());
    REFIT:
      HistogramArray[j].SetBins(2 * NUMBINS_MAX, minimum, maximum);
      for (i = 0; i < numPts; i++) HistogramArray[j].Fill(holder[i]);
      // HistogramArray[j].Draw();
      TF1* f = new TF1("peak", "gaus");
      average = 0.5 * (minimum + maximum);
      f->SetParameters(average, 0.1 * average);
      HistogramArray[j].Fit(f, "q");
      average = f->GetParameter("Mean");
      if (average <= 0.) {
        minimum -= 10.;
        goto REFIT;
      }
      if (j == 2 && (average < (eMin + eMax) / 4. || average > (eMin + eMax))) {
        minimum = (eMin + eMax) / 4.;
        maximum = (eMin + eMax) * 1.;
        goto REFIT;
      }
      printf("%lf\t", average);
      printf("%lf\t", f->GetParameter("Sigma") / average * 100.);
      band[j][0] = f->GetParError(1);
      band[j][1] = f->GetParError(2) / average * 1e2;
      band[j][2] = f->GetChisquare() / double(f->GetNDF());
      delete f;
    }

    cout << endl
         << "+/-"
         << "\t\t"
         << "+/-"
         << "\t\t"
         << "+/-"
         << "\t\t"
         << "+/-"
         << "\t\t"
         << "+/-"
         << "\t\t"
         << "+/-";
    printf("\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", band[0][0], band[0][1],
           band[1][0], band[1][1], band[2][0], band[2][1]);
    cout << "X^2 / DOF"
         << "\t\t\t"
         << "X^2 / DOF"
         << "\t\t\t"
         << "X^2 / DOF";
    printf("\n%lf\t\t\t%lf\t\t\t%lf\tEnd Table\n", band[0][2], band[1][2],
           band[2][2]);

    delete[] HistogramArray;
    return;
  }

  // cout << E_keV.size() << "\t" << S2cor_phd.size() << endl;

  inputs.resize(numBins, vector<double>(1, -999.));
  outputs.resize(numBins, vector<double>(1, -999.));
  if (usePD <= 0) {
    inputs = GetBand(S1cor_phe, S2cor_phe, true);
    for (o = 0; o < numBins; o++) {
      band[o][0] = 0.;
      band[o][1] = 0.;
      band[o][2] = 0.;
      band[o][3] = 0.;
      band[o][4] = 0.;
      band[o][5] = 0.;
      band[o][6] = 0.;
      band[o][7] = 0.;
      band[o][8] = 0.;
    }
    outputs = GetBand_Gaussian(GetBand(S1cor_phe, S2cor_phe, false));
  } else if (usePD == 1) {
    inputs = GetBand(S1cor_phd, S2cor_phd, true);
    for (o = 0; o < numBins; o++) {
      band[o][0] = 0.;
      band[o][1] = 0.;
      band[o][2] = 0.;
      band[o][3] = 0.;
      band[o][4] = 0.;
      band[o][5] = 0.;
      band[o][6] = 0.;
      band[o][7] = 0.;
      band[o][8] = 0.;
    }
    outputs = GetBand_Gaussian(GetBand(S1cor_phd, S2cor_phd, false));
  } else {
    inputs = GetBand(S1cor_spike, S2cor_phd, true);
    for (o = 0; o < numBins; o++) {
      band[o][0] = 0.;
      band[o][1] = 0.;
      band[o][2] = 0.;
      band[o][3] = 0.;
      band[o][4] = 0.;
      band[o][5] = 0.;
      band[o][6] = 0.;
      band[o][7] = 0.;
      band[o][8] = 0.;
    }
    outputs = GetBand_Gaussian(GetBand(S1cor_spike, S2cor_phd, false));
  }
  if (!loop) {
    if ( verbosity ) {
      fprintf(stdout,
	      "Bin Center\tBin Actual\tBand Mean\tBand Stddev\tBand Skew\tBand Mean Err\tBand Stddev Err\tBand Skew Err\tX^2 NEST\tX^2 ROOT\n");
      for (o = 0; o < numBins; o++) {
        fprintf(stdout, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
          band[o][0], band[o][9], band[o][10], band[o][11], band[o][12], band[o][4], band[o][7], band[o][13], band[o][14], band[o][15],
          band[o][2], band[o][5], band[o][3], band[o][6]);
      }
    }
  }
  return;
}

vector<vector<double> > GetBand(vector<double> S1s, vector<double> S2s,
                                bool resol) {
  vector<vector<double> > signals;
  bool save = resol;
  resol = false;
  signals.resize(numBins, vector<double>(1, -999.));
  double binWidth, border;
  if (useS2 == 2) {
    binWidth = (maxS2 - minS2) / double(numBins);
    border = minS2;
  } else {
    binWidth = (maxS1 - minS1) / double(numBins);
    border = minS1;
  }
  int i = 0, j = 0;
  double s1c, numPts;
  unsigned long reject[NUMBINS_MAX] = {0};

  if (resol) {
    numBins = 1;
    binWidth = DBL_MAX;
  }

  if (loop) {
    for (i = 0; i < S1s.size(); i++) S1s[i] *= g1x;
    for (i = 0; i < S2s.size(); i++) S2s[i] *= g2x;
  }

  for (i = 0; i < S1s.size(); i++) {
    for (j = 0; j < numBins; j++) {
      s1c = border + binWidth / 2. + double(j) * binWidth;
      if (i == 0 && !resol) band[j][0] = s1c;
      if ((S1s[i] == 0. || S2s[i] == 0.) && j == 0) reject[j]++;
      if (fabs(S1s[i]) > (s1c - binWidth / 2.) &&
          fabs(S1s[i]) <= (s1c + binWidth / 2.)) {
        if (S1s[i] >= 0. && S2s[i] >= 0.) {
          if (save) {
            signals[j].push_back(S1s[i]);
            continue;
          } else {
            if (useS2 == 0) {
              if (S1s[i] && S2s[i])
                signals[j].push_back(log10(S2s[i] / S1s[i]));
              else
                signals[j].push_back(0.);
            } else if (useS2 == 1) {
              if (S1s[i] && S2s[i])
                signals[j].push_back(log10(S2s[i]));
              else
                signals[j].push_back(0.);
            } else {
              if (S1s[i] && S2s[i])
                signals[j].push_back(log10(S1s[i] / S2s[i]));
              else
                signals[j].push_back(0.);
            }
          }
          band[j][2] += signals[j].back();
          if (resol)
            band[j][0] += S1s[i];
          else
            band[j][1] += S1s[i];
        } else
          reject[j]++;
        break;
      }
    }
  }

  for (j = 0; j < numBins; j++) {
    if (band[j][0] <= 0. && !resol)
      band[j][0] = border + binWidth / 2. + double(j) * binWidth;
    signals[j].erase(signals[j].begin());
    numPts = (double)signals[j].size();
    if (numPts <= 0 && resol) {
      for (i = 0; i < S1s.size(); i++) band[j][0] += fabs(S1s[i]);
      numPts = S1s.size();
    }
    if (resol) band[j][0] /= numPts;
    band[j][1] /= numPts;
    band[j][2] /= numPts;
    for (i = 0; i < (int)numPts; i++) {
      if (signals[j][i] != -999.)
        band[j][3] += pow(signals[j][i] - band[j][2], 2.);  // std dev calc
    }
    for (i = 0; i < S1s.size(); i++) {
      if (resol && S1s[i] > 0.0 && S2s[i] > 0.0)
        band[j][1] += pow(S1s[i] - band[j][0], 2.);  // std dev calc
    }
    band[j][3] /= numPts - 1.;
    band[j][3] = sqrt(band[j][3]);
    if (resol) {
      band[j][1] /= numPts - 1.;
      band[j][1] = sqrt(band[j][1]);
    }
    band[j][5] = band[j][3] / sqrt(numPts);
    band[j][6] =
        numPts /
        (numPts +
         double(reject[j]));  // cerr << numPts << " " << reject[j] << endl;
  }

  return signals;
}

vector<vector<double> > GetBand_Gaussian(vector<vector<double> > signals) {
  int j = 0;
  TH1F* HistogramArray = new TH1F[numBins];

  for (j = 0; j < numBins; j++) {
    TString HistName;
    HistName.Form("%i", j);
    HistogramArray[j].SetName(HistName.Data());
    if ( skewness ) {
      logMin = band[j][2]-3.*band[j][3];
      logMax = band[j][2]+3.*band[j][3];
    }
    HistogramArray[j].SetBins(logBins, logMin, logMax);
    for (unsigned long i = 0; i < signals[j].size(); i++)
      HistogramArray[j].Fill(signals[j][i]);
    HistogramArray[j].Draw();
    if ( !skewness ) {
      TF1* f = new TF1("band", "gaus");
      f->SetParameters(signals[j].size(),band[j][2],band[j][3]);
      HistogramArray[j].Fit(f, "Q");
      band[j][2] = f->GetParameter("Mean");
      band[j][3] = f->GetParameter("Sigma");
      band[j][4] = 0.;
      band[j][5] = f->GetParError(1);
      band[j][6] = f->GetParError(2);
      band[j][7] = 0.;
      band[j][8] = f->GetChisquare() / (double)f->GetNDF();
      delete f;
    }
    else {

      TF1* f = new TF1("skewband",
          "([0]/([2]*sqrt(2.*TMath::Pi())))*exp(-0.5*(x-[1])^2/[2]^2)*(1.+TMath::Erf([3]*(x-[1])/([2]*sqrt(2.))))",
          logMin - 0.5, logMax + 0.5);
      //equation inspired by Vetri Velan
      double amplEstimate = signals[j].size();
      double alphaEstimate = EstimateSkew(band[j][2],band[j][3],signals[j]);
      double deltaEstimate = alphaEstimate / sqrt ( 1. + alphaEstimate * alphaEstimate );
      double omegaEstimate = band[j][3] / sqrt ( 1. - 2. * deltaEstimate * deltaEstimate / TMath::Pi() );
      double xiEstimate = band[j][2] - omegaEstimate * deltaEstimate * sqrt(2./TMath::Pi());
    RETRY:
      f->SetParameters(amplEstimate, xiEstimate, omegaEstimate, alphaEstimate);
      TFitResultPtr res = HistogramArray[j].Fit(f, "QS");
      double fit_xi = res->Value(1);      double fit_xi_err = res->ParError(1);
      double fit_omega = res->Value(2);   double fit_omega_err = res->ParError(2);
      double fit_alpha = res->Value(3);   double fit_alpha_err = res->ParError(3);
      double fit_delta = fit_alpha / sqrt(1 + fit_alpha * fit_alpha);
      TMatrixDSym fit_cov = res->GetCovarianceMatrix();
      double fit_cov_xo = fit_cov[1][2];
      double fit_cov_xa = fit_cov[1][3];
      double fit_cov_oa = fit_cov[2][3];

      // Calculate band mean, std dev, and skewness
      band[j][2] = fit_xi + fit_omega * fit_delta * sqrt(2. / TMath::Pi());
      band[j][3] = fit_omega * sqrt(1 - 2 * fit_delta * fit_delta / TMath::Pi());
      band[j][4] = fit_alpha;

      // Calculate error on mean
      double dmudxi = fit_omega * sqrt(2. / TMath::Pi()) * fit_delta;
      double dmudomega = sqrt(2. / TMath::Pi()) * fit_delta;
      double dmudalpha = fit_omega * sqrt(2. / TMath::Pi()) * pow(1 + fit_alpha * fit_alpha, -1.5);
      band[j][5] = sqrt(dmudxi * dmudxi * fit_xi_err * fit_xi_err + dmudomega * dmudomega * fit_omega_err * fit_omega_err +
        dmudalpha * dmudalpha * fit_alpha_err * fit_alpha_err + 2 * dmudxi * dmudomega * fit_cov_xo +
        2 * dmudomega * dmudalpha * fit_cov_oa + 2 * dmudxi * dmudalpha * fit_cov_xa);

      // Calculate error on standard deviation
      double dvardomega = 2 * fit_omega * (1. - 2. / TMath::Pi() * fit_delta * fit_delta);
      double dvardalpha = -4. / TMath::Pi() * fit_omega * fit_omega * fit_alpha / (1. + fit_alpha * fit_alpha) / (1. + fit_alpha * fit_alpha);
      double var_err = sqrt(dvardomega * dvardomega * fit_omega_err * fit_omega_err + 
        dvardalpha * dvardalpha * fit_alpha_err * fit_alpha_err + 2. * dvardomega * dvardalpha * fit_cov_oa);
      band[j][6] = 0.5 / band[j][3] * var_err;

      // Store other parameters directly from the fit
      band[j][7] = fit_alpha_err;
      band[j][9] = fit_xi;
      band[j][10] = fit_xi_err;
      band[j][11] = fit_omega;
      band[j][12] = fit_omega_err;
      band[j][13] = fit_cov_xo;
      band[j][14] = fit_cov_xa;
      band[j][15] = fit_cov_oa;

      // Store reduced chi2
      band[j][16] = res->Chi2() / res->Ndf();

      // Calculate reduced chi2 manually
      double chiSq = 0.00, modelValue, xValue, denom;
      for ( int k = 0; k < logBins; k++ ) {
        xValue = logMin + k * ( logMax - logMin ) / logBins;
        modelValue = f->GetParameter(0)*exp(-0.5*pow(xValue-f->GetParameter(1),2.)/(f->GetParameter(2)*f->GetParameter(2)))*
          (1.+erf(f->GetParameter(3)*(xValue-f->GetParameter(1))/(f->GetParameter(2)*sqrt(2.)))) / (f->GetParameter(2)*sqrt(2.*TMath::Pi()));
        double denom = max(float(modelValue+HistogramArray[j][k]),(float)1.);
        chiSq += pow ( double(HistogramArray[j][k]) - modelValue, 2. ) / denom; //combined Pearson-Neyman chi-squared (Matthew Sz.)
      }
      chiSq /= ( double(logBins) - 4. - 1. );
      band[j][8] = chiSq;

      // Retry fit if it does not converge.
      if ( chiSq > 10. || band[j][2] > 7. || band[j][5] > 2. || band[j][7] > 10. || band[j][3] > 1. || band[j][6] > 0.5 || band[j][2] <= 0.) {
          xiEstimate = fit_xi;
          omegaEstimate = fit_omega;
          alphaEstimate = 0.0;
          cerr << "Re-fitting...\n";
          goto RETRY;
      }

      delete f;
    }
  }
  
  delete[] HistogramArray;
  return signals;
}

double expectedUlFc(double mub, TFeldmanCousins fc)

{
  double cumProb = 0.0;
  double expectedUl = 0.0;
  for (int i = 0; i < (20 + 3 * mub); i++) {
    double prob = TMath::Poisson(i, mub);
    cumProb += prob;
    double lim = fc.CalculateUpperLimit(i, mub);
    expectedUl += prob * lim;
  }
  return expectedUl;
}

int modPoisRnd(double poisMean, double preFactor) {
  // Must have ROOT, to use PoissonD, the continuous, real version

  int randomNumber;

  if (preFactor >= 1.) {  // essentially, a Poisson convolved with a Gaussian
    poisMean = r.Gaus(poisMean, sqrt((preFactor - 1.) * poisMean));
    if (poisMean < 0.) poisMean = 0.;
    randomNumber = r.Poisson(poisMean);
  } else
    randomNumber =  // squeezes a Poisson thinner, opposite of above
        int(floor(r.PoissonD(poisMean / preFactor) * preFactor + r.Rndm()));

  if (randomNumber < 0) randomNumber = 0;  // just in case

  return randomNumber;
}

int modBinom(int nTot, double prob, double preFactor) {
  // based on work by Aaron Manalaysay; not used, very slow
  int randomNumber, nTotp;

  if (preFactor >= 1.) {
    prob = r.Gaus(prob, sqrt((preFactor - 1.) * prob * nTot) / double(nTot));
    if (prob < 0.) prob = 0.;
    if (prob > 1.) prob = 1.;
    randomNumber = r.Binomial(nTot, prob);
  } else {
    nTotp = int(floor(double(nTot) / preFactor + r.Rndm()));
    randomNumber = int(
        floor((double(nTot) / double(nTotp)) * double(r.Binomial(nTotp, prob)) +
              r.Rndm()));
  }

  if (randomNumber < 0) randomNumber = 0;  // just in case

  return randomNumber;
}

double EstimateSkew ( double mean, double sigma, vector<double> data ) {
  
  double xi = mean;
  double omega = sigma;
  double skew = 0.0;
  for ( int i = 0; i < (int)data.size(); i++ )
    skew += pow((data[i]-xi)/omega,3.);
  skew /= (double)data.size();
  skew = skew / fabs(skew) * std::min(0.7, fabs(skew)); // Sample skewness could theoretically be < -1 or > 1
  double delta = skew / fabs(skew) * sqrt(TMath::Pi() / 2. * pow(fabs(skew),(2./3.)) / (pow(fabs(skew),(2./3.)) + pow(((4. - TMath::Pi()) / 2.),(2./3.))));
  double alpha = delta / sqrt ( 1. - delta * delta );
  return alpha;
  
}

double owens_t (double h, double a) {

  int nIntSteps = 2500;
  double T = 0;
  double dx = a / nIntSteps;
  double x = 0;
  double integrand = 0;

  for (int k = 0; k < nIntSteps; k++) {
    x = (k + 0.5) * dx;
    integrand = exp(-0.5 * h * h * (1 + x * x)) / (1 + x * x);
    T += (integrand * dx);
  }

  T /= (2. * TMath::Pi());
  return T;
}
