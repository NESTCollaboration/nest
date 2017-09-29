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
#include <iostream>

using namespace std;
using namespace NEST;

/*
 * 
 */


int main ( int argc, char** argv ) {
  
  double keV, m_e, aa, ZZ, xMax, yMax, xMin, yMin, FuncValue, B, x;
  vector<double> xyTry(3);
  
  if (argc < 7) {
    cout << "This program takes 6 (or 7) inputs." << endl;
    cout << "numEvts type_interaction E_min[keV] E_max[keV] density[g/cm^3] field_drift[V/cm] {optional:seed}" << endl;
    cout << "for 8B or WIMPs, numEvts is kg-days of exposure (still an integer)" << endl;
    return 0;
  }
  
  unsigned long int numEvts = atoi(argv[1]);
  
  string type = argv[2];
  INTERACTION_TYPE type_num;
  if (type == "NR") type_num = NR;
  else if (type == "WIMP")type_num = WIMP;
  else if (type == "B8") { type_num = B8; numEvts = (unsigned long int)(floor(atof(argv[1])*1.1/365.+0.5)); }
  else if (type == "DD") type_num = DD;
  else if (type == "AmBe")type_num = AmBe;
  else if (type == "Cf") type_num = Cf;
  else if (type == "ion") type_num = ion;
  else if (type == "gamma")type_num = gammaRay;
  else if (type == "beta")type_num = beta;
  else if (type == "CH3T")type_num = CH3T;
  else type_num = Kr83m;
  
  double eMin = atof(argv[3]); xMin = eMin;
  double eMax = atof(argv[4]); xMax = eMax;
  double rho = atof(argv[5]);
  double field = atof(argv[6]);
  
  fprintf(stdout, "E [keV]\tNph\tNe-\n");
  NEST::NESTcalc n; if ( argc >= 8 ) n.SetRandomSeed(atoi(argv[7]));
  for (unsigned long int j = 0; j < numEvts; j++) {
    xyTry[2] = 1.;
    switch (type_num) {
    case CH3T:
      m_e = 510.9989461; //e- rest mass-energy [keV]
      aa = 0.0072973525664; //fine structure constant
      ZZ = 2.;
      if(xMax>18.5898)xMax=18.5898; //tritium beta decay endpoint [keV]
      if(xMin<0.)xMin=0.;
      yMax = 1.1e7; //top of the beta decay E histogram
      xyTry[0] = xMin+(xMax-xMin)*n.rand_uniform();
      xyTry[1] = yMax * n.rand_uniform();
      while ( xyTry[2] > 0. ) {
	B = sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e) / (xyTry[0] + m_e);
	x = (2.*M_PI*ZZ*aa)*(xyTry[0] + m_e)/sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e);
	FuncValue = (sqrt(2.*xyTry[0]*m_e) *
		(xyTry[0] + m_e) *
		(18.5898-xyTry[0]) * (18.5898-xyTry[0]) *
		x*(1./(1.-exp(-x)))*(1.002037-0.001427*(B)));
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      break;
    case B8: //normalize this to ~3500 / 10-ton / year, for E-threshold of 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
      if(xMax>3.)xMax=3.;
      if(xMin<0.)xMin=0.;
      yMax = 28666.6;
      xyTry[0] = xMin+(xMax-xMin)*n.rand_uniform();
      xyTry[1] = yMax * n.rand_uniform();
      B =   5.8335;
      x = 104.35404758803;
      while ( xyTry[2] > 0. ) {
	FuncValue = (3484889.3845409*(pow(xyTry[0],B)+2.4276915393899*pow(xyTry[0],B-1.)+x)/pow(pow(xyTry[0],B)+x,2.)-4728.26850918)*pow(.090455252325381,xyTry[0]);
	xyTry = n.VonNeumann(xMin,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      break;
    case AmBe: //for ZEPLIN-III FSR from HA (Pal '98)
      if(xMax>200.)xMax=200.;
      if(xMin<0.00)xMin=0.00;
      xyTry[0] = xMin+(xMax-xMin)*n.rand_uniform();
      xyTry[1] = n.rand_uniform();
      while ( xyTry[2] > 0. ) {
	FuncValue = exp(-sqrt(xyTry[0]))*(1.+pow(pow(xyTry[0]/31.566,5.),0.45132));
	xyTry = n.VonNeumann(xMin,xMax,0.,1.,xyTry[0],xyTry[1],FuncValue);
      }
      break;
    default:
      xyTry[0] = eMin + (eMax - eMin) * n.rand_uniform();
      break;
    }
    keV = xyTry[0];
    if (keV > eMax) keV = eMax;
    if (keV < eMin) keV = eMin;
    
    
    NEST::YieldResult yields = n.GetYields(type_num, keV, rho, field);
    cout << keV << "\t" << yields.PhotonYield << "\t" << yields.ElectronYield << endl;
  }
  
  return 1;
  
}
