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

#include <string.h>
#include <vector>

#include "../analysis.hh"
#include "../include/TestSpectra.hh"  // contains the WIMP, earth, and escape velocities

// Set the compile mode:
//#define FIT //outputs the goodness of fit for one band (Gaussian centroids of
//histogram in S1 slices)
//#define LIMIT //outputs wimp masses and cross-sections for given efficiency
// default if neither is selected is to provide the ER BG discrimination &
// leakage frac

// non-input parameters hardcoded in
#define CL 0.90  // confidence level
#define VSTEP \
  1e-3  // step size in keV for convolving WIMP recoil energy with efficiency

using namespace std;

vector<vector<double> > GetBand_Gaussian(vector<vector<double> > signals);
vector<vector<double> > GetBand(vector<double> S1s, vector<double> S2s,
                                bool resol);
TRandom3 r;  // r.SetSeed(10);
double band[NUMBINS_MAX][7];
void GetFile(char* fileName);
vector<vector<double> > outputs, inputs;

double WIMP_dRate(double ER, double mWimp);
int SelectRanXeAtom();  // to determine the isotope of Xe
long double Factorial(double x);
double expectedUlFc(double mub, TFeldmanCousins fc);

int modPoisRnd(double poisMean, double preFactor);  // not used currently, but
                                                    // included jic needed for
                                                    // future use, for
                                                    // stretching and squeezing
                                                    // Poissons
int modBinom(int nTot, double prob,
             double preFactor);  // just like above, but for binomial

// Convert all velocities from km/s into cm/s
const double v_0 = V_WIMP * 1e5;      // peak WIMP velocity
const double v_esc = V_ESCAPE * 1e5;  // escape velocity
const double v_e = V_EARTH * 1e5;     // the Earth's velocity

bool loop = false;
double g1x = 1.0, g2x = 1.0;  // for looping over small changes in g1 and g2

int main(int argc, char** argv) {
  bool leak, ERis2nd;
  double band2[numBins][7], NRbandX[numBins], NRbandY[numBins],
      numSigma[numBins], leakage[numBins], discrim[numBins],
      errorBars[numBins][2];
  int i = 0;

  if (argc < 2) {
    cout << endl << "This program takes 1 (or 2) inputs." << endl << endl;
    cout << "One input means you are just doing a band of Gaussians." << endl
         << endl;
    cout << "Two inputs means you're doing both ER and NR and calculating "
            "leakage and discrimination."
         << endl;
    cout << "If you write NR first, you're seeking characterization of the "
            "non-Gaussian asymmetry of NR band."
         << endl;
    cout << "If you write ER first, you're finding non-Gaussian leakage of ER "
            "into NR band."
         << endl
         << endl;
    return 0;
  } else if (argc == 2)
    leak = false;
  else
    leak = true;

#ifdef LIMIT

  if (argc < 3) {
    cerr << "Enter 0 for SI, 1 for SD-n, and 2 for SD-p" << endl;
    return 0;
  }

  FILE* ifp = fopen(argv[1], "r");
  int ch, nLines = 0;
  while (EOF != (ch = getc(ifp))) {
    if ('\n' == ch) nLines++;
  }
  double energy[nLines], efficiency[nLines];
  rewind(ifp);
  for (i = 0; i < nLines; i++)
    fscanf(ifp, "%lf %lf", &energy[i], &efficiency[i]);
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
            "positive or at least 0)"
         << endl;
    return 0;
  }

  double Ul,
      v;  // Start of code block to evaluate upper limit on # of events to use
  if (numBGeventsExp == 0.) {
    for (v = 0.; v < 1e3; v += VSTEP) {
      double sum = 0.0;
      for (i = 0; i < (numBGeventsObs + 1.); i++)
        sum += exp(-v) * pow(v, i) / Factorial(double(i));  // using Poisson
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
      if (j > loE)
        sigAboveThr[i] += VSTEP * WIMP_dRate(j, mass[i]) * eff * xEff *
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

  return 1;

#endif

#ifdef FIT

  if (numBins == 1) {
    GetFile(argv[1]);
    return 1;
  }

  int freeParam;
  cout << "Number of free parameters, for calculating DOF, for chi^2: ";
  cin >> freeParam;
  int DoF = numBins - freeParam;
  FILE* ifp = fopen(argv[2], "r");
  for (i = 0; i < numBins; i++) {
    fscanf(ifp, "%lf %lf %lf %lf %lf %lf", &band2[i][0], &band2[i][1],
           &band2[i][2], &band2[i][4], &band2[i][3], &band2[i][5]);
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
      double error, chi2[2] = {0., 0.};
      for (i = 0; i < numBins; i++) {
        error = sqrt(pow(band[i][4], 2.) + pow(band2[i][4], 2.));
        chi2[0] += pow((band2[i][2] - band[i][2]) / error, 2.);
        error = sqrt(pow(band[i][5], 2.) + pow(band2[i][5], 2.));
        chi2[1] += pow((band2[i][3] - band[i][3]) / error, 2.);
      }
      chi2[0] /= double(DoF - 1);
      chi2[1] /= double(DoF - 1);
      cout.precision(3);
      cout << "The reduced CHI^2 = " << chi2[0] << " for mean, and " << chi2[1]
           << " for width" << endl;
      if (!loop) break;
    }
    if (!loop) break;
  }  // double curly bracket goes with g1x and g2x loops above
  return 1;

#endif

  if (leak) {
    GetFile(argv[2]);
    for (i = 0; i < numBins; i++) {
      band2[i][0] = band[i][0];
      band2[i][1] = band[i][1];
      band2[i][2] = band[i][2];
      band2[i][3] = band[i][3];
      band2[i][4] = band[i][4];
      band2[i][5] = band[i][5];
      band2[i][6] = band[i][6];
    }
  }
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
      if (ERis2nd) {
        numSigma[i] = (band2[i][2] - band[i][2]) / band2[i][3];
        errorBars[i][0] =
            (band2[i][2] - band2[i][4] - band[i][2] - band[i][4]) /
            (band2[i][3] + band2[i][5]);  // upper (more leakage)
        errorBars[i][1] =
            (band2[i][2] + band2[i][4] - band[i][2] + band[i][4]) /
            (band2[i][3] - band2[i][5]);  // lower (less leakage)
        NRbandX[i] = band[i][0];
        NRbandY[i] = band[i][2];
      } else {
        numSigma[i] = (band[i][2] - band2[i][2]) / band[i][3];
        errorBars[i][0] =
            (band[i][2] - band[i][4] - band2[i][2] - band2[i][4]) /
            (band[i][3] + band[i][5]);  // upper (more leakage)
        errorBars[i][1] =
            (band[i][2] + band[i][4] - band2[i][2] + band2[i][4]) /
            (band[i][3] - band[i][5]);  // lower (less leakage)
        NRbandX[i] = band2[i][0];
        NRbandY[i] = band2[i][2];
      }
      leakage[i] = (1. - erf(numSigma[i] / sqrt(2.))) / 2.;
      errorBars[i][0] = (1. - erf(errorBars[i][0] / sqrt(2.))) / 2.;
      errorBars[i][1] = (1. - erf(errorBars[i][1] / sqrt(2.))) / 2.;
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
        return 0;
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

  return 1;
}

void GetFile(char* fileName) {
  FILE* ifp = fopen(fileName, "r");
  double a, b, c, d, e, f, g, h, i, j, k, l, m, n;
  double eMin = 1e100, eMax = -1e100;
  int o;
  char line[256];
  vector<double> E_keV, electricField, tDrift_us, X_mm, Y_mm, Z_mm, Nph, Ne,
      S1cor_phe, S2cor_phe, S1raw_phe, S1cor_phd, S1cor_spike, Ne_Extr,
      S2raw_phe, S2cor_phd;

  rewind(ifp);
  while (fgets(line, sizeof(line), ifp)) {
    sscanf(line,
           "%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
           &a, &b, &c, &d, &e, &f, &g, &h, &i, &j, &k, &l, &m, &n);
    if (a < DBL_MIN && b < DBL_MIN && c < DBL_MIN && d < DBL_MIN &&
        e < DBL_MIN && f < DBL_MIN && g < DBL_MIN && h < DBL_MIN &&
        i < DBL_MIN && j < DBL_MIN && k < DBL_MIN && l < DBL_MIN &&
        m < DBL_MIN && n < DBL_MIN)
      continue;
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
    }
    outputs = GetBand_Gaussian(GetBand(S1cor_spike, S2cor_phd, false));
  }
  if (!loop) {
    // fprintf(stdout,"Bin Center\tBin Actual\tHist Mean\tMean Error\tHist
    // Sigma\t\tEff[%%>thr]\n");
    fprintf(stdout,
            "Bin Center\tBin Actual\tGaus Mean\tMean Error\tGaus Sigma\tSig "
            "Error\tX^2/DOF\n");
    for (o = 0; o < numBins; o++) {
      // fprintf(stdout,"%lf\t%lf\t%lf\t%lf\t%lf\t\t%lf\n",band[o][0],band[o][1],band[o][2],band[o][4],band[o][3],band[o][5]*100.);
      fprintf(stdout, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", band[o][0],
              band[o][1], band[o][2], band[o][4], band[o][3], band[o][5],
              band[o][6]);
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
              if (S1s[i] && S2s[i] && log10(S2s[i] / S1s[i]) > logMin &&
                  log10(S2s[i] / S1s[i]) < logMax)
              signals[j].push_back(log10(S2s[i] / S1s[i]));
              else signals[j].push_back(0.);
            } else if (useS2 == 1) {
              if (S1s[i] && S2s[i] && log10(S2s[i]) > logMin &&
                  log10(S2s[i]) < logMax)
              signals[j].push_back(log10(S2s[i]));
              else signals[j].push_back(0.);
            } else {
              if (S1s[i] && S2s[i] && log10(S1s[i] / S2s[i]) > logMin &&
                  log10(S1s[i] / S2s[i]) < logMax)
              signals[j].push_back(log10(S1s[i] / S2s[i]));
              else signals[j].push_back(0.);
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
    band[j][4] = band[j][3] / sqrt(numPts);
    band[j][5] =
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
    HistogramArray[j].SetBins(50, logMin, logMax);
    for (unsigned long i = 0; i < signals[j].size(); i++)
      HistogramArray[j].Fill(signals[j][i]);
    HistogramArray[j].Draw();
    TF1* f = new TF1("band", "gaus");
    f->SetParameters(1., 0.1);
    HistogramArray[j].Fit(f, "Q");
    band[j][2] = f->GetParameter("Mean");
    band[j][3] = f->GetParameter("Sigma");
    band[j][4] = f->GetParError(1);
    band[j][5] = f->GetParError(2);
    band[j][6] = f->GetChisquare() / (double)f->GetNDF();
    delete f;
  }

  delete[] HistogramArray;
  return signals;
}

// The following was copied wholesale from TestSpectra.cpp on 06-27-2018.
// Find a way to just compile or link against this
//------++++++------++++++------++++++------++++++------++++++------++++++------
// dR() //generator written by Vic Gehman originally
//------++++++------++++++------++++++------++++++------++++++------++++++------

// This spectrum comes from Phys. Rev. D 82 (2010) 023530 (McCabe)
double WIMP_dRate(double ER, double mWimp) {
  // We are going to hard code in the astrophysical halo for now.  This may be
  // something that we make an argument later, but this is good enough to start.
  // Some constants:
  double M_N = 0.9395654;                  // Nucleon mass [GeV]
  double N_A = NEST_AVO;                   // Avogadro's number [atoms/mol]
  double c = 2.99792458e10;                // Speed of light [cm/s]
  double GeVperAMU = 0.9315;               // Conversion factor
  double SecondsPerDay = 60. * 60. * 24.;  // Conversion factor
  double KiloGramsPerGram = 0.001;         // Conversion factor
  double keVperGeV = 1.e6;                 // Conversion factor
  double SqrtPi = pow(M_PI, 0.5);
  double root2 = sqrt(2.);

  // Define the detector Z and A and the mass of the target nucleus
  double Z = ATOM_NUM;
  double A = (double)SelectRanXeAtom();
  double M_T = A * GeVperAMU;

  // Calculate the number of target nuclei per kg
  double N_T = N_A / (A * KiloGramsPerGram);

  // Rescale the recoil energy and the inelastic scattering parameter into GeV
  ER /= keVperGeV;
  double delta = 0. / keVperGeV;  // Setting this to a nonzero value will allow
  // for inelastic dark matter...
  // Set up your dummy WIMP model (this is just to make sure that the numbers
  // came out correctly for definite values of these parameters, the overall
  // normalization of this spectrum doesn't matter since we generate a definite
  // number of events from the macro).
  double m_d = mWimp;       // [GeV]
  double sigma_n = 1.e-36;  //[cm^2] 1 pb reference
  // Calculate the other factors in this expression
  double mu_ND = mWimp * M_N / (mWimp + M_N);  // WIMP-nucleON reduced mass
  double mu_TD = mWimp * M_T / (mWimp + M_T);  // WIMP-nucleUS reduced mass
  double fp =
      1.;  // Neutron and proton coupling constants for WIMP interactions.
  double fn = 1.;

  // Calculate the minimum velocity required to give a WIMP with energy ER
  double v_min = 0.;
  if (ER != 0.) {
    v_min = c * (((M_T * ER) / mu_TD) + delta) / (root2 * sqrt(M_T * ER));
  }
  double bet = 1.;

  // Start calculating the differential rate for this energy bin, starting
  // with the velocity integral:
  double x_min = v_min / v_0;  // Use v_0 to rescale the other velocities
  double x_e = v_e / v_0;
  double x_esc = v_esc / v_0;
  // Calculate overall normalization to the velocity integral
  double N = SqrtPi * SqrtPi * SqrtPi * v_0 * v_0 * v_0 *
             (erf(x_esc) -
              (4. / SqrtPi) * exp(-x_esc * x_esc) *
                  (x_esc / 2. + bet * x_esc * x_esc * x_esc / 3.));
  // Calculate the part of the velocity integral that isn't a constant
  double zeta = 0.;
  int thisCase = -1;
  if ((x_e + x_min) < x_esc) {
    thisCase = 1;
  }
  if ((x_min > fabs(x_esc - x_e)) && ((x_e + x_esc) > x_min)) {
    thisCase = 2;
  }
  if (x_e > (x_min + x_esc)) {
    thisCase = 3;
  }
  if ((x_e + x_esc) < x_min) {
    thisCase = 4;
  }
  switch (thisCase) {
    case 1:
      zeta = ((SqrtPi * SqrtPi * SqrtPi * v_0 * v_0) / (2. * N * x_e)) *
             (erf(x_min + x_e) - erf(x_min - x_e) -
              ((4. * x_e) / SqrtPi) * exp(-x_esc * x_esc) *
                  (1 + bet * (x_esc * x_esc - x_e * x_e / 3. - x_min * x_min)));
      break;
    case 2:
      zeta = ((SqrtPi * SqrtPi * SqrtPi * v_0 * v_0) / (2. * N * x_e)) *
             (erf(x_esc) + erf(x_e - x_min) -
              (2. / SqrtPi) * exp(-x_esc * x_esc) *
                  (x_esc + x_e - x_min -
                   (bet / 3.) * (x_e - 2. * x_esc - x_min) *
                       (x_esc + x_e - x_min) * (x_esc + x_e - x_min)));
      break;
    case 3:
      zeta = 1. / (x_e * v_0);
      break;
    case 4:
      zeta = 0.;
      break;
    default:
      cerr << "\tThe velocity integral in the WIMP generator broke!!!" << endl;
      exit(0);
  }

  double a = 0.52;                           // in fm
  double C = 1.23 * pow(A, 1. / 3.) - 0.60;  // fm
  double s = 0.9;  // skin depth of nucleus in fm. Originally used by Karen
                   // Gibson; XENON100 1fm; 2.30 acc. to Lewin and Smith maybe?
  double rn = sqrt(C * C + (7. / 3.) * M_PI * M_PI * a * a -
                   5. * s * s);  // alternatives: 1.14*A^1/3 given in L&S, or
                                 // rv=1.2*A^1/3 then rn =
                                 // sqrt(pow(rv,2.)-5.*pow(s,2.)); used by
                                 // XENON100 (fm)
  double q = 6.92 * sqrt(A * ER);  // in units of 1 over distance or length
  double FormFactor;
  if (q * rn > 0.)
    FormFactor =
        3. * exp(-0.5 * q * q * s * s) * (sin(q * rn) - q * rn * cos(q * rn)) /
        (q * rn * q * rn * q * rn);  // qr and qs unitless inside Bessel
                                     // function, which is dimensionless too
  else
    FormFactor = 1.;

  // Now, the differential spectrum for this bin!
  double dSpec = 0.5 * (c * c) * N_T * (RHO_NAUGHT / m_d) *
                 (M_T * sigma_n / (mu_ND * mu_ND));
  // zeta=1.069-1.4198*ER+.81058*pow(ER,2.)-.2521*pow(ER,3.)+.044466*pow(ER,4.)-0.0041148*pow(ER,5.)+0.00013957*pow(ER,6.)+2.103e-6*pow(ER,7.);
  // if ( ER > 4.36 ) squiggle = 0.; //parameterization for 7 GeV WIMP using
  // microMegas
  dSpec *= (((Z * fp) + ((A - Z) * fn)) / fn) *
           (((Z * fp) + ((A - Z) * fn)) / fn) * zeta * FormFactor * FormFactor *
           SecondsPerDay / keVperGeV;

  return dSpec;
}

int SelectRanXeAtom() {
  int A;
  double isotope = r.Uniform() * 100.;

  if (isotope > 0.000 && isotope <= 0.090)
    A = 124;
  else if (isotope > 0.090 && isotope <= 0.180)
    A = 126;
  else if (isotope > 0.180 && isotope <= 2.100)
    A = 128;
  else if (isotope > 2.100 && isotope <= 28.54)
    A = 129;
  else if (isotope > 28.54 && isotope <= 32.62)
    A = 130;
  else if (isotope > 32.62 && isotope <= 53.80)
    A = 131;
  else if (isotope > 53.80 && isotope <= 80.69)
    A = 132;
  else if (isotope > 80.69 && isotope <= 91.13)
    A = 134;
  else
    A = 136;
  return A;
}

long double Factorial(double x) { return tgammal(x + 1.); }

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
