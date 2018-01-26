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
#include <iostream>
using namespace NEST;
using namespace std;

double NEST::CH3T_spectrum(double xMin,double xMax, NESTcalc& n){
    double m_e = 510.9989461; //e- rest mass-energy [keV]
    double aa = 0.0072973525664; //fine structure constant
    double ZZ = 2.;
    if(xMax>18.5898)xMax=18.5898; //tritium beta decay endpoint [keV]
    if(xMin<0.)xMin=0.;
    double yMax = 1.1e7; //top of the beta decay E histogram
    vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
			    yMax * n.rand_uniform(),1.};
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
    if(xMax>4.)xMax=4.;
    if(xMin<0.)xMin=0.;
    double yMax = pow(10.,-2.198);
    vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
			    yMax * n.rand_uniform(), 1.};
    while ( xyTry[2] > 0. ) {
      double FuncValue = 2.198 + 1.2184*xyTry[0] - 0.32849*pow(xyTry[0],2.) + 0.12441*pow(xyTry[0],3.);
      FuncValue = pow(10.,-FuncValue);
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

double NEST::DD_spectrum(double xMin, double xMax, NESTcalc& n){  //JV LUX, most closely like JENDL-4. See arXiv:1608.05381. Lower than G4/LUXSim
    if(xMax>75.)xMax=75.;
      if(xMin<0.000)xMin=0.000;
      double yMax = 1.1694e+6;
      vector<double> xyTry = {xMin+(xMax-xMin)*n.rand_uniform(),
       yMax * n.rand_uniform(),1.};
      while ( xyTry[2] > 0. ) {
        double FuncValue =
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
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      return xyTry[0];
}
