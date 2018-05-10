/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TestSpectra.cpp
 * Author: brodsky3
 * 
 * Created on December 11, 2017, 10:27 AM
 */

#include "TestSpectra.hh"

using namespace std;
//using namespace NEST;

double power = 3.7488;

double TestSpectra::CH3T_spectrum ( double xMin, double xMax ) {
  
  double m_e = 510.9989461; //e- rest mass-energy [keV]
  double aa = 0.0072973525664; //fine structure constant
  double ZZ = 2.;
  
  if(xMax>18.5898)xMax=18.5898; //tritium beta decay endpoint [keV]
  if(xMin<0.)xMin=0.;
  double yMax = 1.1e7; //top of the beta decay E histogram
  vector<double> xyTry = {xMin+(xMax-xMin)*RandomGen::rndm()->rand_uniform(),
			  yMax * RandomGen::rndm()->rand_uniform(),1.};
  while ( xyTry[2] > 0. ) {
    double B = sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e) / (xyTry[0] + m_e);
    double x = (2.*M_PI*ZZ*aa)*(xyTry[0] + m_e)/sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e);
    double FuncValue = (sqrt(2.*xyTry[0]*m_e) *
			(xyTry[0] + m_e) *
			(18.5898-xyTry[0]) * (18.5898-xyTry[0]) *
			x*(1./(1.-exp(-x)))*(1.002037-0.001427*(B)));
    xyTry = RandomGen::rndm()->VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
  }
  return xyTry[0];
  
}

double TestSpectra::B8_spectrum ( double xMin, double xMax ) {
  
  if(xMax!=4.)xMax=4.;
  if(xMin!=0.)xMin=0.;
  double yMax = pow(10.,-2.198);
  vector<double> xyTry = {xMin+(xMax-xMin)*RandomGen::rndm()->rand_uniform(),
			  yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while ( xyTry[2] > 0. ) {
    double FuncValue = 2.198 + 1.2184*xyTry[0] - 0.32849*pow(xyTry[0],2.) + 0.12441*pow(xyTry[0],3.);
    FuncValue = pow(10.,-FuncValue);
    xyTry = RandomGen::rndm()->VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
  }
  return xyTry[0];
  
}

double TestSpectra::AmBe_spectrum ( double xMin, double xMax ) {
  
  if(xMax>200.)
    xMax=200.;
  if ( xMin < DBL_MIN ) xMin = DBL_MIN;
  double yMax = pow(10.,power), yMin = 0.0;
  vector<double> xyTry = {xMin+(xMax-xMin)*RandomGen::rndm()->rand_uniform(),
			  yMax * RandomGen::rndm()->rand_uniform(),1.};
  while ( xyTry[2] > 0. ) {
    double FuncValue =
       power * pow(log10(xyTry[0]),0.)
      -0.77942*pow(log10(xyTry[0]),1.)
      +1.30300*pow(log10(xyTry[0]),2.)
      -2.75280*pow(log10(xyTry[0]),3.)
      +1.57310*pow(log10(xyTry[0]),4.)
      -0.30072*pow(log10(xyTry[0]),5.);
    FuncValue = pow(10.,FuncValue);
    xyTry = RandomGen::rndm()->VonNeumann(xMin,xMax,yMin,yMax,xyTry[0],xyTry[1],FuncValue);
  }
  
  return xyTry[0];
  
}

double TestSpectra::Cf_spectrum ( double xMin, double xMax ) {
  
  if(xMax>200.)
    xMax=200.;
  if ( xMin < DBL_MIN ) xMin = DBL_MIN;
  double yMax = 2.*pow(10.,power), yMin = 0.0;
  vector<double> xyTry = {xMin+(xMax-xMin)*RandomGen::rndm()->rand_uniform(),
			  yMax * RandomGen::rndm()->rand_uniform(),1.};
  while ( xyTry[2] > 0. ) {
    double FuncValue =
       power * pow(log10(xyTry[0]),0.)
      -0.77942*pow(log10(xyTry[0]),1.)
      +1.30300*pow(log10(xyTry[0]),2.)
      -2.75280*pow(log10(xyTry[0]),3.)
      +1.57310*pow(log10(xyTry[0]),4.)
      -0.30072*pow(log10(xyTry[0]),5.);
    FuncValue = pow(10.,FuncValue);
    FuncValue *=
      1.9929
      - .033214 * pow(xyTry[0],1.)
      +.00032857* pow(xyTry[0],2.)
      -1.000e-6 * pow(xyTry[0],3.);
    xyTry = RandomGen::rndm()->VonNeumann(xMin,xMax,yMin,yMax,xyTry[0],xyTry[1],FuncValue);
  }
  
  return xyTry[0];
  
}

double TestSpectra::DD_spectrum( double xMin, double xMax ){  //JV LUX, most closely like JENDL-4. See arXiv:1608.05381. Lower than G4/LUXSim
  
  if(xMax>80.)xMax=80.;
  if(xMin<0.000)xMin=0.000;
  double yMax = 1.1694e+6;
  vector<double> xyTry = {xMin+(xMax-xMin)*RandomGen::rndm()->rand_uniform(),
			  yMax * RandomGen::rndm()->rand_uniform(),1.};
  while ( xyTry[2] > 0. ) {
    double FuncValue = //1.*exp(-0.15*xyTry[0])+2e-3*exp(0.05*xyTry[0]); //LUXSim version (Carmen)
       1.1694e+6*pow(xyTry[0],0.)
      -1.4733e+5*pow(xyTry[0],1.)
      + 8507.0 * pow(xyTry[0],2.)
      - 273.59 * pow(xyTry[0],3.)
      + 4.3216 * pow(xyTry[0],4.)
      +0.0097428*pow(xyTry[0],5.)
      -0.0017966*pow(xyTry[0],6.)
      +3.4069e-5*pow(xyTry[0],7.)
      -2.918e-7 *pow(xyTry[0],8.)
      +9.973e-10*pow(xyTry[0],9.);
    FuncValue /= 1.+0.85*(
			  -.016698/pow(xyTry[0]-75.,1.)+
			  8.04540/pow(xyTry[0]-75.,2.)+
			  105.000/pow(xyTry[0]-75.,3.)+
			  582.400/pow(xyTry[0]-75.,4.)+
			  1218.50/pow(xyTry[0]-75.,5.)+
			  1250.90/pow(xyTry[0]-75.,6.)+
			  659.680/pow(xyTry[0]-75.,7.)+
			  161.110/pow(xyTry[0]-75.,8.)+
			  11.7710/pow(xyTry[0]-75.,9.));
    xyTry = RandomGen::rndm()->VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
  }
  return xyTry[0];
  
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
//dR() //generator written by Vic Gehman originally
//------++++++------++++++------++++++------++++++------++++++------++++++------

//This spectrum comes from Phys. Rev. D 82 (2010) 023530 (McCabe)
double TestSpectra::WIMP_dRate ( double ER, double mWimp ) {
  
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
  double A = (double)RandomGen::rndm()->SelectRanXeAtom();
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

TestSpectra::WIMP_spectrum_prep TestSpectra::WIMP_prep_spectrum ( double mass, double eStep ) {
  
  WIMP_spectrum_prep spectrum;
  double EnergySpec[10001]={0}, divisor, x1, x2;
  int numberPoints;
  
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
	  break;
	}
    }
  
  spectrum.divisor = divisor; return spectrum;
  
}

double TestSpectra::WIMP_spectrum( WIMP_spectrum_prep wimp_spectrum, double mass ){
  
  double xMin = 0., FuncValue = 0.00, x = 0.;
  double yMax = WIMP_dRate ( xMin, mass );
  vector<double> xyTry ={ xMin + (wimp_spectrum.xMax - xMin) * RandomGen::rndm()->rand_uniform(),
			  yMax * RandomGen::rndm()->rand_uniform(), 1. };
  while ( xyTry[2] > 0. )
    {
      while ( xyTry[1] > (-WIMP_dRate(0.,mass)/wimp_spectrum.xMax*xyTry[0]+WIMP_dRate(0.,mass)) ) { //triangle cut more efficient than rectangle
	xyTry[0] = (wimp_spectrum.xMax-xMin)*RandomGen::rndm()->rand_uniform(); xyTry[1] = yMax*RandomGen::rndm()->rand_uniform(); }
      for ( x = 0; x < wimp_spectrum.xMax; x+=(1./wimp_spectrum.divisor) )
	{
	  if ( xyTry[0] > x && xyTry[0] < (x + 1./wimp_spectrum.divisor) )
	    {
	      FuncValue = wimp_spectrum.base[int(x*wimp_spectrum.divisor)] * exp(-wimp_spectrum.exponent[int(x*wimp_spectrum.divisor)] * xyTry[0]);
	      break;
	    }
	}
      xyTry = RandomGen::rndm()->VonNeumann(xMin, wimp_spectrum.xMax, 0., yMax, xyTry[0], xyTry[1], FuncValue);
    }
  
  return xyTry[0];
  
}
