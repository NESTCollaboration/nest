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

using namespace std;

vector< vector<double> > GetBand_Gaussian ( vector< vector<double> > signals );
vector< vector<double> > GetBand ( vector<double> S1s, vector<double> S2s, bool resol );
TRandom3 r;
double band[NUMBINS_MAX][7];
void GetFile ( char* fileName );
vector< vector<double> > outputs, inputs;

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
    gr1->Fit(fitf,"rq","",minS1,maxS1);
    double chi2 = fitf->GetChisquare() / (double)fitf->GetNDF();
    if ( chi2 > 2. )
      cerr << "ERROR: Poor fit to NR Gaussian band centroids i.e. means of log(S2) or log(S2/S1) histograms in S1 bins. Investigate please!" << endl;
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
  int ch, nLines, o;
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
  fprintf(stdout,"Bin Center\tBin Actual\tHist Mean\tMean Error\tHist Sigma\tSig Error\tX^2/DOF\n");
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
    HistogramArray[j].SetBins(40,0.6,3.);
    for ( unsigned long i = 0; i < signals[j].size(); i++ )
      HistogramArray[j].Fill(signals[j][i]);
    HistogramArray[j].Draw();
    TF1 *f = new TF1("band","gaus");
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
