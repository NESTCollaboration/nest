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
using namespace NEST;
using namespace std;

double NEST::CH3T_spectrum(double xMin,double xMax, NESTcalc& n){
    double m_e = 510.9989461; //e- rest mass-energy [keV]
    double aa = 0.0072973525664; //fine structure constant
    double ZZ = 2.;
    if(xMax>18.5898)xMax=18.5898; //tritium beta decay endpoint [keV]
    if(xMin<0.)xMin=0.;
    double yMax = 1.1e7; //top of the beta decay E histogram
    vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(), xyTry[1] = yMax * n.rand_uniform(),1.};
    while ( xyTry[2] > 0. ) {
      double B = sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e) / (xyTry[0] + m_e);
      double x = (2.*M_PI*ZZ*aa)*(xyTry[0] + m_e)/sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e);
      double FuncValue = (sqrt(2.*xyTry[0]*m_e) *
              (xyTry[0] + m_e) *
              (18.5898-xyTry[0]) * (18.5898-xyTry[0]) *
              x*(1./(1.-exp(-x)))*(1.002037-0.001427*(B)));
      xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
    }
    return xyTry[0];
}

double NEST::B8_spectrum(double xMin, double xMax, NESTcalc& n){
    if(xMax>3.)xMax=3.;
      if(xMin<0.)xMin=0.;
      double yMax = 28666.6;
      vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
            xyTry[1] = yMax * n.rand_uniform(), 1.};
      double B =   5.8335;
      double x = 104.35404758803;
      while ( xyTry[2] > 0. ) {
	double FuncValue = (3484889.3845409*(pow(xyTry[0],B)+2.4276915393899*pow(xyTry[0],B-1.)+x)/pow(pow(xyTry[0],B)+x,2.)-4728.26850918)*pow(.090455252325381,xyTry[0]);
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      return xyTry[0];
}

double NEST::AmBe_spectrum(double xMin, double xMax, NESTcalc& n){
    if(xMax>200.)xMax=200.;
      if(xMin<0.00)xMin=0.00;
    double yMax = 1.;
      vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
            yMax * n.rand_uniform(),1.};
      while ( xyTry[2] > 0. ) {
	double FuncValue = exp(-sqrt(xyTry[0]))*(1.+pow(pow(xyTry[0]/31.566,5.),0.45132));
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      return xyTry[0];
}

double NEST::Cf_spectrum(double xMin, double xMax, NESTcalc& n){
    if(xMax>200.)xMax=200.;
      if(xMin<0.00)xMin=0.00;
    double yMax = 2.;
    
      vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
            yMax * n.rand_uniform(),1.};
      while ( xyTry[2] > 0. ) {
	double FuncValue = exp(-sqrt(xyTry[0]))*(1.+pow(pow(xyTry[0]/31.566,5.),0.45132));
	 FuncValue *=
				1.9929
				- .033214 * pow(xyTry[0],1.)
			        +.00032857* pow(xyTry[0],2.)
				-1.000e-6 * pow(xyTry[0],3.);
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      return xyTry[0];
}

double NEST::DD_spectrum(double xMin, double xMax, NESTcalc& n){
    if(xMax>74.15)xMax=74.15;
      if(xMin<0.000)xMin=0.000;
      double yMax = 1.2368e+6;
      vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
       yMax * n.rand_uniform(),1.};
      while ( xyTry[2] > 0. ) {
        double FuncValue =
	   1.2368e+6*pow(xyTry[0],0.)
	  -1.5285e+5*pow(xyTry[0],1.)
	  + 7531.4 * pow(xyTry[0],2.)
	  - 120.36 * pow(xyTry[0],3.)
	  - 4.1124 * pow(xyTry[0],4.)
	  +0.23758 * pow(xyTry[0],5.)
	  -0.0047675*pow(xyTry[0],6.)
	  +4.4851e-5*pow(xyTry[0],7.)
	  -1.6501e-7*pow(xyTry[0],8.);
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      return xyTry[0];
}