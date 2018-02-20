/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testNEST.cpp
 * Author: brodsky3
 *
 * Created on August 1, 2017, 1:03 PM
 */

#include <NEST.hh>
#include <TestSpectra.hh>
#include <float.h>
#include <iostream>

using namespace std;
using namespace NEST;
/*
 * 
 */
double dRate ( double ER, double mWimp );
NEST::NESTcalc n; std::default_random_engine generator;
int main ( int argc, char** argv ) {
  
  double EnergySpec[101], base[100], exponent[100];
  vector<double> xyTry, scint(8);
  double xMin, xMax, yMax, FuncValue;
  
  if (argc < 7) {
    cout << "This program takes 6 (or 7) inputs." << endl << endl;
    cout << "numEvts type_interaction E_min[keV] E_max[keV] density[g/cm^3] field_drift[V/cm] {optional:seed}" << endl;
    cout << "for 8B or WIMPs, numEvts is kg-days of exposure (still an integer)" << endl << endl;
    cout << "exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] density[g/cm^3] field_drift[V/cm] {optional:seed}" << endl;
    return 0;
  }
  unsigned long int numEvts = atoi(argv[1]);
  
    string type = argv[2];
    INTERACTION_TYPE type_num;
    if (type == "NR") type_num = NR;
    else if (type == "WIMP") {
      type_num = WIMP; int i;
      for ( i = 0; i < 101; i++ ) EnergySpec[i] = dRate ( double(i), atof(argv[3]) );
      double integral = 0.0;
      for ( i = 0; i < 100; i++ ) {
	base[i] = EnergySpec[i]*pow(EnergySpec[i]/EnergySpec[i+1],i);
	exponent[i] = log(EnergySpec[i]/EnergySpec[i+1]);
	if ( base[i] > 0. && base[i] < DBL_MAX && exponent[i] > 0. && exponent[i] < DBL_MAX )
	  integral += base[i]*(1./exponent[i]-exp(-exponent[i])/exponent[i])*exp(-exponent[i]*i);
	else { xMax = double(i-1); break; }
      }
      std::poisson_distribution<int> distribution(integral*atof(argv[1])*atof(argv[4])/1e-36);
      numEvts = distribution(generator);
    }
    else if (type == "B8") {
      type_num = B8;
      std::poisson_distribution<int> distribution(0.0026*atof(argv[1]));
      numEvts = distribution(generator);
    } else if (type == "DD") type_num = DD;
    else if (type == "AmBe")type_num = AmBe;
    else if (type == "Cf") type_num = Cf;
    else if (type == "ion") type_num = ion;
    else if (type == "gamma")type_num = gammaRay;
    else if (type == "Kr83m")type_num=Kr83m;
    else if (type == "CH3T")type_num = CH3T;
    else type_num = beta;

    double eMin = atof(argv[3]);
    double eMax = atof(argv[4]);
    double rho = atof(argv[5]);
    double field = atof(argv[6]);
    
    if ( type_num == Kr83m && eMin == 9.4 && eMax == 9.4 )
      fprintf(stdout, "t [ns]\t\tE [keV]\t\tNph\tNe-\tS1_raw [PE]\tS1_Zcorr\tS1c_spike\tNe-X\tS2_rawArea\tS2_Zcorr [phd]\n");
    else
      fprintf(stdout, "E [keV]\t\tNph\tNe-\tS1_raw [PE]\tS1_Zcorr\tS1c_spike\tNe-X\tS2_rawArea\tS2_Zcorr [phd]\n");

    if (argc >= 8) n.SetRandomSeed(atoi(argv[7]));
    
    double keV = -999;
    for (unsigned long int j = 0; j < numEvts; j++) {
      if (eMin == eMax) {
	keV = eMin;
      } else {
	switch (type_num) {
	case CH3T:
	  keV = CH3T_spectrum(eMin, eMax, n);
	  break;
	case B8: //normalize this to ~3500 / 10-ton / year, for E-threshold of 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
	  keV = B8_spectrum(eMin, eMax, n);
	  break;
	case AmBe: //for ZEPLIN-III FSR from HA (Pal '98)
	  keV = AmBe_spectrum(eMin, eMax, n);
	  break;
	case Cf:
	  keV = Cf_spectrum(eMin, eMax, n);
	  break;
	case DD:
	  keV = DD_spectrum(eMin, eMax, n);
	  break;
	case WIMP:
	  xMin = 0.;
	  yMax = dRate ( xMin, eMin );
	  xyTry = { xMin+(xMax-xMin) * n.rand_uniform(),
		    yMax * n.rand_uniform(), 1. };
	  while ( xyTry[2] > 0. ) {
	    for ( double x = 0; x < xMax; x++ ) {
	      if ( xyTry[0] > x && xyTry[0] < (x+1.) ) {
		FuncValue = base[int(x)]*exp(-exponent[int(x)]*xyTry[0]);
		break;
	      }
	    }
	    xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
	  }
	  keV = xyTry[0];
	  break;
	default:
	  keV = eMin + (eMax - eMin) * n.rand_uniform();
	  break;
	}
      }
      
      if ( type_num != WIMP ) {
	if (keV > eMax) keV = eMax;
	if (keV < eMin) keV = eMin;
      }
      
      NEST::YieldResult yields = n.GetYields ( type_num, keV, rho, field );
      scint = n.GetS1(int(floor(yields.PhotonYield+0.5)));
      printf("%.6f\t%.6f\t%.6f\t", keV, yields.PhotonYield, yields.ElectronYield);
      printf("%.6f\t%.6f\t%.6f\t", scint[2], scint[5], scint[7]);
      scint = n.GetS2(int(floor(yields.ElectronYield+0.5)));
      printf("%i\t%.6f\t%.6f\n", (int)scint[0], scint[4], scint[7]);
    }
    
    return 1;
    
}

vector<double> NESTcalc::GetS1 ( int Nph ) {
  
  vector<double> scintillation(8);
  int coinLevel = 3,numPMTs = 100;
  double g1 = 0.10, sPEres = 0.5, P_dphe = 0.2, sPEeff = 0.92, sPEthr = 0.25;
  double posDep = 0.9 + n.rand_uniform()*(1.1-0.9);
  
  int nHits=n.BinomFluct(Nph,g1*posDep), Nphe = 0;
  double pulseArea = 0., spike = 0.;
  if ( sPEeff < 1. || sPEthr > 0. ) {
    for ( int i = 0; i < nHits; i++ ) {
      double phe1 = n.rand_gauss(1.,sPEres); Nphe++;
      if ( n.rand_uniform() > sPEeff ) phe1 = 0.;
      double phe2 = 0.;
      if ( n.rand_uniform() < P_dphe ) {
	phe2 = n.rand_gauss(1.,sPEres); Nphe++;
	if ( n.rand_uniform() > sPEeff ) phe2 = 0.;
      }
      if ( (phe1+phe2) > sPEthr ) { spike++; pulseArea += phe1 + phe2; }
    }
  }
  else {
    Nphe = nHits + n.BinomFluct(nHits,P_dphe);
    pulseArea = n.rand_gauss(Nphe,sPEres*sqrt(Nphe));
    spike = (double)nHits;
  }
  if ( pulseArea < 0. ) pulseArea = 0.;
  double pulseAreaC= pulseArea / posDep;
  double Nphd = pulseArea / (1.+P_dphe);
  double NphdC= pulseAreaC/ (1.+P_dphe);
  double spikeC = spike / posDep;
  
  scintillation[0] = nHits; scintillation[1] = Nphe;
  scintillation[2] = pulseArea; scintillation[3] = pulseAreaC;
  scintillation[4] = Nphd; scintillation[5] = NphdC;
  scintillation[6] = spike; scintillation[7] = spikeC;
  
  if ( spike >= coinLevel && n.rand_uniform() < 1.-pow((double)numPMTs,1.-(double)spike) )
    { ; }
  else {
    scintillation[0] *= -1.;
    scintillation[1] *= -1.;
    scintillation[2] *= -1.;
    scintillation[3] *= -1.;
    scintillation[4] *= -1.;
    scintillation[5] *= -1.;
    scintillation[6] *= -1.;
    scintillation[7] *= -1.;
  }
  
  return scintillation;
  
}

vector<double> NESTcalc::GetS2 ( int Ne ) {
  
  vector<double> ionization(8);
  double alpha = 0.137, beta = 177., gamma = 45.7, eLife_us = 500., P_dphe = 0.2, sPEres = 0.5, Fano = 3., S2botTotRatio = 0.4;
  double driftTime = 0.0 + n.rand_uniform()*(500.-0.0);
  double g1_gas = 0.10, gasGap_cm = 0.5,p_bar = 1.5,E_gas=10.,epsilon=1.85/1.00126;
  
  double E_liq = E_gas / epsilon; //kV per cm
  double ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645; // arXiv:1710.11032
  if ( ExtEff > 1. ) ExtEff = 1.;
  if ( ExtEff < 0. ) ExtEff = 0.;
  int Nee = n.BinomFluct(Ne,ExtEff*exp(-driftTime/eLife_us));
  
  double elYield = Nee*
    (alpha*E_gas*1000.-beta*p_bar-gamma)*
    gasGap_cm; // arXiv:1207.2292
  int Nph = int(floor(rand_gauss(elYield,sqrt(Fano*elYield))+0.5));
  int nHits = n.BinomFluct(Nph,g1_gas);
  int Nphe = nHits + n.BinomFluct(nHits,P_dphe);
  double pulseArea=n.rand_gauss(Nphe,sPEres*sqrt(Nphe));
  double pulseAreaC= pulseArea/exp(-driftTime/eLife_us);
  double Nphd = pulseArea / (1.+P_dphe);
  double NphdC= pulseAreaC/ (1.+P_dphe);
  
  double S2b = n.rand_gauss(S2botTotRatio*pulseArea,sqrt(S2botTotRatio*pulseArea*(1.-S2botTotRatio)));
  double S2bc= S2b / exp(-driftTime/eLife_us);
  
  ionization[0] = Nee; ionization[1] = Nph;
  ionization[2] = nHits; ionization[3] = Nphe;
  ionization[4] = pulseArea; ionization[5] = pulseAreaC;
  ionization[6] = Nphd; ionization[7] = NphdC;
  
  return ionization;
  
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
//dR() //generator written by Vic Gehman originally
//------++++++------++++++------++++++------++++++------++++++------++++++------

//This spectrum comes from Phys. Rev. D 82 (2010) 023530 (McCabe)
double dRate ( double ER, double mWimp ) {
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
  double A = 131.293;
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
    cout << "\tThe velocity integral in the WIMP generator broke!!!" << endl;
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
