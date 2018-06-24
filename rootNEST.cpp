#include <iostream>
#include <cmath>
#include <fstream>

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFeldmanCousins.h"
#include "TRandom3.h"
#include "TMath.h"

#include <string.h>
#include <vector>

#include "analysis.hh"

#define NUMBINS_MAX 200
//#define FIT
//#define LIMIT
#define CL 0.90 //confidence level
#define VSTEP 1e-3

using namespace std;

vector< vector<double> > GetBand_Gaussian ( vector< vector<double> > signals );
vector< vector<double> > GetBand ( vector<double> S1s, vector<double> S2s, bool resol );
TRandom3 r;
double band[NUMBINS_MAX][7];
void GetFile ( char* fileName );
vector< vector<double> > outputs, inputs;

struct WIMP_spectrum_prep {
  double base[100];
  double exponent[100];
  double integral;
  double xMax;
  double divisor;
};
WIMP_spectrum_prep wimp_spectrum_prep;

double WIMP_spectrum ( WIMP_spectrum_prep wimp_spectrum, double mass );
double WIMP_dRate ( double ER, double mWimp );
WIMP_spectrum_prep WIMP_prep_spectrum ( double mass, double eStep );
int SelectRanXeAtom ( );
vector<double> VonNeumann ( double xMin, double xMax, double yMin, double yMax,
			   double xTest, double yTest,double fValue );
long double Factorial ( double x );
double expectedUlFc ( double mub, TFeldmanCousins fc );

int main ( int argc, char** argv ) {
  
  bool leak, ERis2nd;
  double band2[NUMBINS_MAX][7], NRbandX[NUMBINS_MAX], NRbandY[NUMBINS_MAX], numSigma[NUMBINS_MAX], leakage[NUMBINS_MAX], discrim[NUMBINS_MAX], errorBars[NUMBINS_MAX][2];
  int i = 0;
  
  if ( argc < 2 ) {
    cout << endl << "This program takes 1 (or 2) inputs." << endl << endl;
    cout << "One input means you are just doing a band of Gaussians." << endl << endl;
    cout << "Two inputs means you're doing both ER and NR and calculating leakage and discrimination." << endl;
    cout << "If you write NR first, you're seeking characterization of the non-Gaussian asymmetry of NR band." << endl;
    cout << "If you write ER first, you're finding non-Gaussian leakage of ER into NR band." << endl << endl;
    return 0;
  }
  else if ( argc == 2 )
    leak =false;
  else
    leak = true;
  
#ifdef LIMIT
  
  FILE * ifp = fopen ( argv[1], "r" );
  int ch, nLines = 0;
  while ( EOF != ( ch = getc ( ifp ) ) ) {
    if ( '\n' == ch ) nLines++;
  }
  double energy[nLines], efficiency[nLines];
  rewind ( ifp );
  for ( i = 0; i < nLines; i++ )
    fscanf ( ifp, "%lf %lf", &energy[i], &efficiency[i] );
  TGraph* gr1 = new TGraph ( nLines, energy, efficiency );
  TF1* fitf = new TF1 ( "FracEffVkeVEnergy", "10.^(2.-[0]*exp(-[1]*x^[2])-[3]*exp(-[4]*x^[5]))/100.", 0.000, energy[nLines-1] ); //eqn inspired by Alex Murphy
  fitf->SetParameters(10.,2.,1.,20.,1e4,-2.5);
  gr1->Fit(fitf,"rq","",0.000,energy[nLines-1]);
  double aa = fitf->GetParameter(0);
  double bb = fitf->GetParameter(1);
  double cc = fitf->GetParameter(2);
  double dd = fitf->GetParameter(3);
  double ee = fitf->GetParameter(4);
  double ff = fitf->GetParameter(5);
  /*fprintf ( stderr, "Fractional Efficiency versus Energy in keV, Parameters: %f %f %f %f %f %f\n",
	    aa,
	    bb,
	    cc,
	    dd,
	    ee,
	    ff );*/
  delete gr1; delete fitf;
  
  // Get input parameters for sensitivity or limit calculation
  double time, fidMass, loE, hiE, xEff, NRacc, numBGeventsExp, numBGeventsObs;
  cout << "Target Mass (kilograms): ";
  cin >> fidMass;
  cout << "Run Time (provide days): ";
  cin >> time;
  cout << "Multiplicative factor on the efficiency of NR event detection from file (usually ~>0.9): ";
  cin >> xEff;
  cout << "Acceptance for NR events post electron recoil background discrimination (usually ~ 50%): ";
  cin >> NRacc;
    cout << "Number of BG events observed: ";
  cin >> numBGeventsObs;
  if ( numBGeventsObs > 0. ) {
    cout << "Number of BG events expected: ";
    cin >> numBGeventsExp;
  }
  cout << "Minimum energy (keV) for detection: ";
  cin >> loE;
  cout << "Maximum energy (keV) for detection: ";
  cin >> hiE;
  // Make sure inputs were valid.
  if ( cin.fail() ) {
    cerr << endl << "Input error. Make sure all inputs were numbers" << endl;
    return 0;
  }
  
  const int masses = 200; double massMax = 1e5;
  double mass[masses] = { 3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,
			  10,11,12,13,14,15,16,17,18,19,
			  20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,
			  50,55,60,65,70,75,80,85,90,95,
			  100,110,120,130,140,150,160,170,180,190,
			  200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,
			  500,550,600,650,700,750,800,850,900,950,
			  1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,
			  2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4200,4400,4600,4800,
			  5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,1E+4,massMax }; //in GeV
  double eStep, sigAboveThr[masses], xSect[masses]; //array for the cross-sections
  //cout << "\nWIMP Mass [GeV/c^2]\tCross Section [cm^2]" << endl;
  cout << "\nWIMP Mass [GeV/c^2]\tEff's times Acc [frac]" << endl;
  unsigned long eventsCheck, numEvents[masses]; //array for the number of WIMP signal events above threshold
  
  i = 0;
  while ( mass[i] < massMax ) { //Iterate across each sample wimp Mass
    eventsCheck = long ( floor ( 1E7 / mass[i] ) ); //number of events to run (mass-dependent)
    numEvents[i] = 0;
    eStep = 0.5 * pow ( mass[i], 0.5 );
    if ( eStep > E_step )
      eStep = E_step;
    wimp_spectrum_prep = WIMP_prep_spectrum ( mass[i], eStep );
    for ( long j = 0; j < eventsCheck; j++ ) { //Iterate across each event within each sample wimp mass
      double keV = WIMP_spectrum ( wimp_spectrum_prep, mass[i] );
      double eff = pow(10.,2.-aa*exp(-bb*pow(keV,cc))-dd*exp(-ee*pow(keV,ff)))/100.;
      if ( r.Uniform() < eff && keV > loE && keV < hiE ) numEvents[i]++;
    }
    i++;
    sigAboveThr[i-1] = ( (double)numEvents[i-1] / double(eventsCheck) ) * xEff * NRacc;
    cout << mass[i-1] << "\t\t\t" << sigAboveThr[i-1] << endl;
  }
  
  double Ul, v;
  if ( numBGeventsExp == 0. ) {
    for ( v = 0.; v < 1e3; v += VSTEP ) {
      double sum = 0.0;
      for ( i = 0; i < (numBGeventsObs+1.); i++ )
	sum += exp(-v)*pow(v,i)/Factorial(double(i));
      if ( sum <= ( 1. - CL ) ) break;
    }
    Ul = 0.5 * ( 2. * v - VSTEP );
  }
  else {
    TFeldmanCousins fc ( CL );
    if ( numBGeventsExp != numBGeventsObs || numBGeventsObs == int(numBGeventsObs) )
      Ul = fc.CalculateUpperLimit(numBGeventsObs,numBGeventsExp);
    else
      Ul = expectedUlFc ( numBGeventsExp, fc );
    double powCon = fc.CalculateUpperLimit(0.,0.);
    if ( Ul < powCon ) Ul = powCon;
  }
  
  return 1;
  
#endif
  
#ifdef FIT
  
  int freeParam;
  cout << "Number of free parameters, for calculating DOF, for chi^2: ";
  cin >> freeParam;
  int DoF = numBins - freeParam;
  FILE* ifp = fopen ( argv[2], "r" );
  for ( i = 0; i < numBins; i++ ) {
    fscanf ( ifp, "%lf %lf %lf %lf %lf %lf",
	     &band2[i][0],
	     &band2[i][1],
	     &band2[i][2],
	     &band2[i][4],
	     &band2[i][3],
	     &band2[i][5] );
  }
  GetFile ( argv[1] );
  double error, chi2[2] = { 0., 0. };
  for ( i = 0; i < numBins; i++ ) {
    error = sqrt ( pow ( band[i][4], 2. ) + pow ( band2[i][4], 2. ) );
    chi2[0] += pow ( ( band2[i][2] - band[i][2] ) / error, 2. );
    error = sqrt ( pow ( band[i][5], 2. ) + pow ( band2[i][5], 2. ) );
    chi2[1] += pow ( ( band2[i][3] - band[i][3] ) / error, 2. );
  }
  chi2[0] /= double ( DoF - 1 );
  chi2[1] /= double ( DoF - 1 );
  cout.precision ( 3 );
  cout << "The reduced CHI^2 = " << chi2[0] << " for mean, and " << chi2[1] << " for width" << endl;
  return 1;

#endif
  
  if ( leak ) {
    GetFile ( argv[2] );
    for ( i = 0; i < numBins; i++ ) {
      band2[i][0] = band[i][0];
      band2[i][1] = band[i][1];
      band2[i][2] = band[i][2];
      band2[i][3] = band[i][3];
      band2[i][4] = band[i][4];
      band2[i][5] = band[i][5];
      band2[i][6] = band[i][6];
    }
  }
  GetFile ( argv[1] );
  if ( leak ) {
    double finalSums[3] = {0.,0.,0.};
    if ( band2[0][2] > band[0][2] ) {
      ERis2nd = true;
      fprintf(stderr,"Bin Center\tBin Actual\t#StdDev's\tLeak Frac Gaus\t+ error  - error\tDiscrim[%%]\tLower ''Half''\t+ error  - error\tUpperHf[%%]\tignore column!\n");
    }
    else {
      ERis2nd = false;
      fprintf(stderr,"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t- => Gaus more\n");
      fprintf(stderr,"Bin Center\tBin Actual\t#StdDev's\tLeak Frac Gaus\t+ error  - error\tDiscrim[%%]\tLeak Frac Real\t+ error  - error\tDiscrim[%%]\tReal-Gaus Leak\n");
    }
    for ( i = 0; i < numBins; i++ ) {
      if ( ERis2nd ) {
	numSigma[i] = ( band2[i][2] - band[i][2] ) /band2[i][3];
	errorBars[i][0] = ( band2[i][2] - band2[i][4] - band[i][2] - band[i][4] ) / ( band2[i][3] + band2[i][5] ); //upper (more leakage)
	errorBars[i][1] = ( band2[i][2] + band2[i][4] - band[i][2] + band[i][4] ) / ( band2[i][3] - band2[i][5] ); //lower (less leakage)
	NRbandX[i] = band[i][0];
	NRbandY[i] = band[i][2];
      }
      else {
	numSigma[i] = ( band[i][2] - band2[i][2] ) / band[i][3];
	errorBars[i][0] = ( band[i][2] - band[i][4] - band2[i][2] - band2[i][4] ) / ( band[i][3] + band[i][5] ); //upper (more leakage)
	errorBars[i][1] = ( band[i][2] + band[i][4] - band2[i][2] + band2[i][4] ) / ( band[i][3] - band[i][5] ); //lower (less leakage)
	NRbandX[i] = band2[i][0];
	NRbandY[i] = band2[i][2];
      }
      leakage[i] = (1.-erf(numSigma[i]/sqrt(2.)))/2.;
      errorBars[i][0] = (1.-erf(errorBars[i][0]/sqrt(2.)))/2.;
      errorBars[i][1] = (1.-erf(errorBars[i][1]/sqrt(2.)))/2.;
      finalSums[2] += leakage[i];
      discrim[i] = 1. - leakage[i];
    }
    TGraph *gr1 = new TGraph ( numBins, NRbandX, NRbandY );
    TF1 *fitf = new TF1("NRbandGCentroid","[0]/(x+[1])+[2]*x+[3]",minS1,maxS1);
    fitf->SetParameters(10., 15., -1.5e-3, 2.);
    gr1->Fit(fitf,"rq","",minS1,maxS1);
    double chi2 = fitf->GetChisquare() / (double)fitf->GetNDF();
    if ( chi2 > 1.5 ) {
      cerr << "WARNING: Poor fit to NR Gaussian band centroids i.e. means of log(S2) or log(S2/S1) histograms in S1 bins. Investigate please!" << endl;
      fitf = new TF1("NRbandGCentroid","[0]+([1]-[0])/(1+(x/[2])^[3])",minS1,maxS1);
      fitf->SetParameters(0.5, 3., 1e2, 0.4);
      gr1->Fit(fitf,"rq","",minS1,maxS1);
      chi2 = fitf->GetChisquare() / (double)fitf->GetNDF();
      if ( chi2 > 2. )
	{ cerr << "ERROR: Even the backup plan to use sigmoid failed as well!" << endl; return 0; }
    }
    /*fprintf(stderr,"Band Fit Parameters: %f %f %f %f\n",
	    fitf->GetParameter(0),
	    fitf->GetParameter(1),
	    fitf->GetParameter(2),
	    fitf->GetParameter(3));*/
    long below[NUMBINS_MAX] = {0}; double NRbandGCentroid, leakTotal, poisErr[2];
    for ( i = 0; i < numBins; i++ ) {
      below[i] = 0;
      for ( long j = 0; j < inputs[i].size(); j++ ) {
	NRbandGCentroid = fitf->GetParameter(0)/(inputs[i][j]+fitf->GetParameter(1))+fitf->GetParameter(2)*inputs[i][j]+fitf->GetParameter(3); // use Woods function
	//NRbandGCentroid = NRbandY[i]; // use the center of the bin instead of the fit, to compare to past data that did not use a smoothing spline
	//NRbandGCentroid = fitf->GetParameter(0)/(NRbandX[i]+fitf->GetParameter(1))+fitf->GetParameter(2)*NRbandX[i]+fitf->GetParameter(3); // compromise
	if ( outputs[i][j] < NRbandGCentroid ) below[i]++;
      }
      leakTotal = double(below[i])/(double)outputs[i].size();
      poisErr[0] = (double(below[i])+sqrt(double(below[i]))) / (double(outputs[i].size())-sqrt(double(outputs[i].size())));
      poisErr[1] = (double(below[i])-sqrt(double(below[i]))) / (double(outputs[i].size())+sqrt(double(outputs[i].size())));
      fprintf(stderr,"%.2f\t\t%.6f\t%.6f\t%e\t%.2e %.2e\t%.6f\t%e\t%.2e %.2e\t%.6f\t%e\n",
              0.5*(band[i][0]+band2[i][0]),0.5*(band[i][1]+band2[i][1]),numSigma[i],leakage[i],errorBars[i][0]-leakage[i],leakage[i]-errorBars[i][1],discrim[i]*100.,
	      leakTotal,poisErr[0]-leakTotal,leakTotal-poisErr[1],(1.-leakTotal)*100.,leakTotal-leakage[i]);
      finalSums[0] += (double)below[i]; finalSums[1] += (double)outputs[i].size();
    }
    fprintf(stderr,"OVERALL DISCRIMINATION or ACCEPTANCE between min and maxS1 = %.12f%%, total: Gaussian + non-Gaussian ('anomalous'). Leakage Fraction = %.12e\n",
	    ( 1. - finalSums[0] / finalSums[1] ) * 100., finalSums[0] / finalSums[1]);
    double HighValue, LowValue;
    LowValue = ( finalSums[0] + sqrt(finalSums[0]) ) / ( finalSums[1] - sqrt(finalSums[1]) );
    HighValue= ( finalSums[0] - sqrt(finalSums[0]) ) / ( finalSums[1] + sqrt(finalSums[1]) );
    fprintf(stderr,"highest possible discrimination value (plus  1-sigma) is %.12f%%, with corresponding leakage of %.12e\n",(1.-HighValue)*100.,HighValue);
    fprintf(stderr," lowest possible discrimination value (minus 1-sigma) is %.12f%%, with corresponding leakage of %.12e\n",(1.-LowValue )*100., LowValue);
    fprintf(stderr,"OVERALL DISCRIMINATION or ACCEPTANCE between min and maxS1 = %.12f%%, Gaussian.                                     Leakage Fraction = %.12e\n",
            ( 1. - finalSums[2] / numBins      ) * 100., finalSums[2] / numBins     );
    delete gr1; delete fitf;
  }
  
  return 1;
  
}

void GetFile ( char* fileName ) {
  
  FILE *ifp = fopen(fileName,"r");
  double a,b,c,d,e,f,g,h,i,j,k,l,m,n;
  int ch, nLines = 0, o;
  vector<double> E_keV, electricField, tDrift_us, X_mm, Y_mm, Z_mm, Nph, Ne, S1cor_phe, S2cor_phe,
    S1raw_phe, S1cor_phd, S1cor_spike, Ne_Extr, S2raw_phe, S2cor_phd;
  
  while ( EOF != ( ch = getc ( ifp ) ) ) {
    if ( '\n' == ch )
      nLines++;
    if ( nLines == 4 ) break;
  }
  
  while ( true ) {
    fscanf(ifp,"%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
	   &a,&b,&c,&d,&e,&f,&g,&h,&i,&j,&k,&l,&m,&n);
    if ( feof(ifp) )
      break;
    //fprintf(stderr,"%.6f\t%.6f\t%.6f\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
    E_keV.push_back(a);
    electricField.push_back(b);
    tDrift_us.push_back(c);
    X_mm.push_back(d);
    Y_mm.push_back(e);
    Z_mm.push_back(f);
    Nph.push_back(g);
    Ne.push_back(h);
    S1raw_phe.push_back(i);
    if ( usePD <= 0 &&
	 fabs(j*1.2) > minS1 && (j*1.2) < maxS1 ) {
      S1cor_phe.push_back(j*1.2); // here and down below for S2: FIX THIS by getting P_dphe (USUALLY ~0.2) from detector class
      S1cor_phd.push_back(0);
      S1cor_spike.push_back(0);
    }
    else if ( usePD == 1 && fabs(j) > minS1 && j < maxS1 ) {
      S1cor_phe.push_back(0);
      S1cor_phd.push_back(j);
      S1cor_spike.push_back(0);
    }
    else if ( usePD >= 2 && fabs(k) > minS1 && k < maxS1 ) {
      S1cor_phe.push_back(0);
      S1cor_phd.push_back(0);
      S1cor_spike.push_back(k);
    }
    else {
      S1cor_phe.push_back(-999.);
      S1cor_phd.push_back(-999.);
      S1cor_spike.push_back(-999.);
    }
    Ne_Extr.push_back(l);
    S2raw_phe.push_back(m);
    if ( usePD <= 0 &&
	 fabs(n*1.2) > minS2 && (n*1.2) < maxS2 ) {
      S2cor_phe.push_back(n*1.2);
      S2cor_phd.push_back(0);
    }
    else if ( usePD >= 1 && fabs(n) > minS2 && n < maxS2 ) {
      S2cor_phe.push_back(0);
      S2cor_phd.push_back(n);
    }
    else {
      S2cor_phe.push_back(-999.);
      S2cor_phd.push_back(-999.);
    }
  }
  
  //cout << E_keV.size() << "\t" << S2cor_phd.size() << endl;
  
  inputs.resize(numBins,vector<double>(1,-999.));
  outputs.resize(numBins,vector<double>(1,-999.));
  if ( usePD <= 0 ) {
    inputs = GetBand ( S1cor_phe, S2cor_phe, true );
    for ( o = 0; o < numBins; o++ ) {
    band[o][0] = 0.; band[o][1] = 0.; band[o][2] = 0.; band[o][3] = 0.; band[o][4] = 0.; band[o][5] = 0.; band[o][6] = 0.;
    }
    outputs = GetBand_Gaussian ( GetBand ( S1cor_phe, S2cor_phe, false ) );
  }
  else if ( usePD == 1 ) {
    inputs = GetBand ( S1cor_phd, S2cor_phd, true );
    for ( o = 0; o < numBins; o++ ) {
    band[o][0] = 0.; band[o][1] = 0.; band[o][2] = 0.; band[o][3] = 0.; band[o][4] = 0.; band[o][5] = 0.; band[o][6] = 0.;
    }
    outputs = GetBand_Gaussian ( GetBand ( S1cor_phd, S2cor_phd, false ) );
  }
  else {
    inputs = GetBand(S1cor_spike, S2cor_phd, true );
    for ( o = 0; o < numBins; o++ ) {
    band[o][0] = 0.; band[o][1] = 0.; band[o][2] = 0.; band[o][3] = 0.; band[o][4] = 0.; band[o][5] = 0.; band[o][6] = 0.;
    }
    outputs = GetBand_Gaussian ( GetBand(S1cor_spike, S2cor_phd, false ) );
  }
  //fprintf(stdout,"Bin Center\tBin Actual\tHist Mean\tMean Error\tHist Sigma\t\tEff[%%>thr]\n");
  fprintf(stdout,"Bin Center\tBin Actual\tGaus Mean\tMean Error\tGaus Sigma\tSig Error\tX^2/DOF\n");
  for ( o = 0; o < numBins; o++ ) {
    //fprintf(stdout,"%lf\t%lf\t%lf\t%lf\t%lf\t\t%lf\n",band[o][0],band[o][1],band[o][2],band[o][4],band[o][3],band[o][5]*100.);
    fprintf(stdout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",band[o][0],band[o][1],band[o][2],band[o][4],band[o][3],band[o][5],band[o][6]);
  }
  
  return;
  
}

vector< vector<double> > GetBand ( vector<double> S1s,
				   vector<double> S2s, bool resol ) {
  
  vector< vector<double> > signals;
  bool save = resol; resol = false;
  signals.resize(numBins,vector<double>(1,-999.));
  double binWidth, border;
  if ( useS2 == 2 ) {
    binWidth = ( maxS2 - minS2 ) / double(numBins);
    border = minS2;
  }
  else {
    binWidth = ( maxS1 - minS1 ) / double(numBins);
    border = minS1;
  }
  int i = 0, j = 0; double s1c, numPts;
  unsigned long reject[NUMBINS_MAX] = {0};
  
  if ( resol ) {
    numBins = 1;
    binWidth = DBL_MAX;
  }
  
  for ( i = 0; i < S1s.size(); i++ ) {
    for ( j = 0; j < numBins; j++ ) {
      s1c = border + binWidth/2. + double(j) * binWidth;
      if ( i == 0 && !resol ) band[j][0] = s1c;
      if ( (S1s[i] == 0. || S2s[i] == 0.) && j == 0 ) reject[j]++;
      if ( fabs(S1s[i]) > (s1c-binWidth/2.) && fabs(S1s[i]) <= (s1c+binWidth/2.) ) {
        if ( S1s[i] >= 0. && S2s[i] >= 0. ) {
          if ( save ) {
	    signals[j].push_back(S1s[i]); continue;
          }
          else {
            if ( useS2 == 0 )
              { if ( S1s[i] && S2s[i] ) signals[j].push_back(log10(S2s[i]/S1s[i])); else signals[j].push_back(0.); }
            else if ( useS2 == 1 )
              { if ( S1s[i] && S2s[i] ) signals[j].push_back(log10(S2s[i])); else signals[j].push_back(0.); }
            else
              { if ( S1s[i] && S2s[i] ) signals[j].push_back(log10(S1s[i]/S2s[i])); else signals[j].push_back(0.); }
          }
          band[j][2] += signals[j].back();
          if ( resol )
            band[j][0] += S1s[i];
          else
	    band[j][1] += S1s[i];
	}
	else
	  reject[j]++;
        break; }
    }
  }
  
  for ( j = 0; j < numBins; j++ ) {
    if ( band[j][0] <= 0. && !resol ) band[j][0] = border + binWidth/2. + double(j) * binWidth;
    signals[j].erase(signals[j].begin());
    numPts = (double)signals[j].size();
    if (resol)
      band[j][0] /= numPts;
    band[j][1] /= numPts;
    band[j][2] /= numPts;
    for ( i = 0; i < (int)numPts; i++ ) {
      if ( signals[j][i] != -999. ) band[j][3] += pow(signals[j][i]-band[j][2],2.);
      if ( resol && S1s[i] >= 0.0 ) band[j][1] += pow(S1s[i]-band[j][0],2.); //std dev calc
    }
    band[j][3] /= numPts - 1.;
    band[j][3] = sqrt(band[j][3]);
    if ( resol ) {
      band[j][1] /= numPts - 1.;
      band[j][1] = sqrt(band[j][1]);
    }
    band[j][4] = band[j][3] / sqrt ( numPts );
    band[j][5] = numPts/
      (numPts+double(reject[j])); //cerr << numPts << " " << reject[j] << endl;
  }
  
  return signals;
  
}

vector< vector<double> > GetBand_Gaussian ( vector< vector<double> > signals ) {
  
  int j = 0;
  TH1F *HistogramArray = new TH1F[NUMBINS_MAX];
  
  for ( j = 0; j < numBins; j++ ) {
    TString HistName;
    HistName.Form("%i",j);
    HistogramArray[j].SetName(HistName.Data());
    HistogramArray[j].SetBins(50,0.6,3.6);
    for ( unsigned long i = 0; i < signals[j].size(); i++ )
      HistogramArray[j].Fill(signals[j][i]);
    HistogramArray[j].Draw();
    TF1 *f = new TF1("band","gaus");
    f->SetParameters(1.,0.1);
    HistogramArray[j].Fit(f,"Q");
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

// The following was copied wholesale from TestSpectra.cpp on 06-23-2018
//------++++++------++++++------++++++------++++++------++++++------++++++------
//dR() //generator written by Vic Gehman originally
//------++++++------++++++------++++++------++++++------++++++------++++++------

//This spectrum comes from Phys. Rev. D 82 (2010) 023530 (McCabe)
double WIMP_dRate ( double ER, double mWimp ) {
  
  // We are going to hard code in the astrophysical halo for now.  This may be 
  // something that we make an argument later, but this is good enough to start.
  // Some constants:
  double M_N = 0.9395654; //Nucleon mass [GeV]
  double N_A = 6.022e23; //Avagadro's number [atoms/mol]
  double c = 2.99792458e10; //Speed of light [cm/s]
  double GeVperAMU = 0.9315;             //Conversion factor
  double SecondsPerDay = 60. * 60. * 24.;//Conversion factor
  double KiloGramsPerGram = 0.001;       //Conversion factor
  double keVperGeV = 1.e6;               //Conversion factor
  double cmPerkm = 1.e5;                 //Conversion factor
  double SqrtPi = pow(M_PI, 0.5); double root2 = sqrt(2.);
  // Convert all velocities from km/s into cm/s
  double v_0   = 220. * cmPerkm;
  double v_esc = 544. * cmPerkm;
  double v_e   = 232. * cmPerkm;
  
  // Define the detector Z and A and the mass of the target nucleus
  double Z = 54.;
  double A = (double)SelectRanXeAtom();
  double M_T = A * GeVperAMU;
  
  // Calculate the number of target nuclei per kg
  double N_T = N_A / (A * KiloGramsPerGram);
  
  // Rescale the recoil energy and the inelastic scattering parameter into GeV
  ER /= keVperGeV;
  double delta = 0. / keVperGeV; //Setting this to a nonzero value will allow
  //for inelastic dark matter...
  // Set up your dummy WIMP model (this is just to make sure that the numbers 
  // came out correctly for definite values of these parameters, the overall 
  // normalization of this spectrum doesn't matter since we generate a definite 
  // number of events from the macro).
  double rho_D = 0.3;      // [GeV/cm^3]
  double m_d   = mWimp;      // [GeV]
  double sigma_n = 1.e-36; //[cm^2] 1 pb reference
  // Calculate the other factors in this expression
  double mu_ND = mWimp * M_N / (mWimp + M_N);// WIMP-nucleON reduced mass
  double mu_TD = mWimp * M_T / (mWimp + M_T);// WIMP-nucleUS reduced mass
  double fp = 1.;// Neutron and proton coupling constants for WIMP interactions.
  double fn = 1.;
  
  // Calculate the minimum velocity required to give a WIMP with energy ER
  double v_min = 0.;
  if(ER != 0.){
    v_min = c * (((M_T * ER) / mu_TD) + delta) / (root2*sqrt(M_T*ER));
  } double bet = 1.;
  
  // Start calculating the differential rate for this energy bin, starting 
  // with the velocity integral:
  double x_min = v_min / v_0;// Use v_0 to rescale the other velocities
  double x_e   = v_e   / v_0;
  double x_esc = v_esc / v_0;
  // Calculate overall normalization to the velocity integral
  double N = SqrtPi*SqrtPi*SqrtPi*v_0*v_0*v_0*(erf(x_esc)-(4./SqrtPi)*exp(-x_esc*x_esc)*(x_esc/2.+bet*x_esc*x_esc*x_esc/3.));
  // Calculate the part of the velocity integral that isn't a constant
  double zeta = 0.;
  int thisCase = -1;
  if((x_e + x_min) < x_esc){thisCase = 1;}
  if((x_min > fabs(x_esc - x_e)) && ((x_e + x_esc) > x_min)){thisCase = 2;}
  if(x_e > (x_min + x_esc)){thisCase = 3;}
  if((x_e + x_esc) < x_min){thisCase = 4;}
  switch(thisCase){
  case 1:
    zeta=((SqrtPi*SqrtPi*SqrtPi*v_0*v_0)/(2.*N*x_e))*(erf(x_min+x_e)-erf(x_min-x_e)
						      -((4.*x_e)/SqrtPi)*exp(-x_esc*x_esc)*(1+bet*(x_esc*x_esc-x_e*x_e/3.-x_min*x_min)));
    break;
  case 2:
    zeta=((SqrtPi*SqrtPi*SqrtPi*v_0*v_0)/(2.*N*x_e))*(erf(x_esc)+erf(x_e-x_min)
						      -(2./SqrtPi)*exp(-x_esc*x_esc)*(x_esc+x_e-x_min-(bet/3.)*(x_e-2.*x_esc-x_min)*(x_esc+x_e-x_min)*(x_esc+x_e-x_min)));
    break;
  case 3:
    zeta = 1. / (x_e * v_0);
    break;
  case 4:
    zeta = 0.;
    break;
  default:
    cerr << "\tThe velocity integral in the WIMP generator broke!!!" << endl;
  }
  
  double a = 0.52;
  double C = 1.23*pow(A,1./3.)-0.60;
  double s = 0.9;
  double rn= sqrt(C*C+(7./3.)*M_PI*M_PI*a*a-5.*s*s);
  double q = 6.92*sqrt(A*ER); double FormFactor;
  if ( q * rn > 0. ) FormFactor = 3.*exp(-0.5*q*q*s*s)*(sin(q*rn)-q*rn*cos(q*rn))/(q*rn*q*rn*q*rn);
  else FormFactor = 1.;
  
  // Now, the differential spectrum for this bin!
  double dSpec = 0.5 * (c * c) * N_T * (rho_D / m_d) * (M_T * sigma_n / (mu_ND * mu_ND));
  dSpec *= (((Z * fp) + ((A - Z) * fn)) / fn) * (((Z * fp) + ((A - Z) * fn)) / fn) * zeta * FormFactor*FormFactor * SecondsPerDay / keVperGeV;
  
  return dSpec;
  
}

WIMP_spectrum_prep WIMP_prep_spectrum ( double mass, double eStep ) {
  
  WIMP_spectrum_prep spectrum;
  double EnergySpec[10001]={0}, divisor, x1, x2;
  int numberPoints;
  
 RE_START:
  
  if ( mass < 2.0 ) { // GeV/c^2
    divisor = 100 / eStep; if ( (eStep*0.01) > 0.01 ) cerr << "WARNING, <= 0.01 keV step size recommended" << endl;
    numberPoints=int(10000./eStep);
  }
  else if ( mass < 10. ) {
    divisor = 10. / eStep;
    numberPoints = int ( 1000. / eStep );
  }
  else {
    divisor = 1.0 / eStep;
    numberPoints = int ( 100. / eStep );
  }
  
  for ( int i = 0; i < (numberPoints+1); i++ ) {
    EnergySpec[i] = WIMP_dRate( double(i)/divisor, mass );
  }
  
  spectrum.integral = 0.;
  for ( long i = 0; i < 1000000; i++ ) {
    spectrum.integral += WIMP_dRate( double(i)/1e4, mass ) / 1e4;
  }
  
  for ( int i = 0; i < numberPoints; i++ )
    {
      x1 = double(i)/divisor; x2 = double(i+1)/divisor;
      spectrum.base[i] = EnergySpec[i+1] * pow(EnergySpec[i+1] / EnergySpec[i], x2/(x1-x2));
      spectrum.exponent[i] = log(EnergySpec[i+1] / EnergySpec[i]) / ( x1 - x2 );
      if ( spectrum.base[i] > 0. && spectrum.base[i] < DBL_MAX && spectrum.exponent[i] > 0. && spectrum.exponent[i] < DBL_MAX )
	;//spectrum.integral+=spectrum.base[i]/spectrum.exponent[i]*(exp(-spectrum.exponent[i]*x1)-exp(-spectrum.exponent[i]*x2));
      else
	{
	  spectrum.xMax = double(i - 1) / divisor;
	  if ( spectrum.xMax <= 0.0 ) goto RE_START;
          break;
	}
    }
  
  spectrum.divisor = divisor; return spectrum;
  
}

double WIMP_spectrum ( WIMP_spectrum_prep wimp_spectrum, double mass ) {
  
  double xMin = 0., FuncValue = 0.00, x = 0.;
  double yMax = WIMP_dRate ( xMin, mass );
  vector<double> xyTry(3);
  xyTry[2] = 1.;
  xyTry[0] = xMin +  ( wimp_spectrum.xMax - xMin ) * r.Rndm();
  xyTry[1] = yMax * r.Rndm();

  while ( xyTry[2] > 0. )
    {
      while ( xyTry[1] > (-WIMP_dRate(0.,mass)/wimp_spectrum.xMax*xyTry[0]+WIMP_dRate(0.,mass)) ) { //triangle cut more efficient than rectangle
	xyTry[0] = (wimp_spectrum.xMax-xMin)*r.Rndm(); xyTry[1] = yMax*r.Rndm(); }
      for ( x = 0; x < wimp_spectrum.xMax; x+=(1./wimp_spectrum.divisor) )
	{
	  if ( xyTry[0] > x && xyTry[0] < (x + 1./wimp_spectrum.divisor) )
	    {
	      FuncValue = wimp_spectrum.base[int(x*wimp_spectrum.divisor)] * exp(-wimp_spectrum.exponent[int(x*wimp_spectrum.divisor)] * xyTry[0]);
	      break;
	    }
	}
      xyTry = VonNeumann ( xMin, wimp_spectrum.xMax, 0., yMax, xyTry[0], xyTry[1], FuncValue );
    }
  
  return xyTry[0];
  
}

int SelectRanXeAtom ( ) {
  
  int A;
  double isotope = r.Uniform()*100.;
  
       if ( isotope > 0.000 && isotope <= 0.090 )
    A = 124;
  else if ( isotope > 0.090 && isotope <= 0.180 )
    A = 126;
  else if ( isotope > 0.180 && isotope <= 2.100 )
    A = 128;
  else if ( isotope > 2.100 && isotope <= 28.54 )
    A = 129;
  else if ( isotope > 28.54 && isotope <= 32.62 )
    A = 130;
  else if ( isotope > 32.62 && isotope <= 53.80 )
    A = 131;
  else if ( isotope > 53.80 && isotope <= 80.69 )
    A = 132;
  else if ( isotope > 80.69 && isotope <= 91.13 )
    A = 134;
  else
    A = 136;
  return A;
  
}

vector<double> VonNeumann ( double xMin, double xMax, double yMin, double yMax,
			   double xTest, double yTest,double fValue ) {
  
  vector<double> xyTry(3);
  
  xyTry[0] = xTest;
  xyTry[1] = yTest;
  
  if ( xyTry[1] > fValue ) {
    xyTry[0] = xMin + ( xMax - xMin ) * r.Uniform ( );
    xyTry[1] = yMin + ( yMax - yMin ) * r.Uniform ( );
    xyTry[2] = 1.;
  }
  else
    xyTry[2] = 0.;
  
  return xyTry; //doing a vector means you can return 2 values at the same time
  
}

long double Factorial ( double x ) {
  
  return tgammal ( x + 1. );
  
}

double expectedUlFc ( double mub, TFeldmanCousins fc )
  
{
  
  double cumProb = 0.0;
  double expectedUl = 0.0;
  for ( int i = 0; i < (20+3*mub); i++ ) {
    double prob = TMath::Poisson(i,mub);
    cumProb += prob;
    double lim = fc.CalculateUpperLimit(i,mub);
    expectedUl += prob * lim;
  }
  return expectedUl;
  
}
