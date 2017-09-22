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
  
  double keV, m_e, aa, ZZ, xMax, yMax, FuncValue, B, x;
  vector<double> xyTry(3);
  
  if (argc < 7) {
    cout << "This program takes 6 inputs." << endl;
    cout << "numEvts type_interaction E_min[keV] E_max[keV] density[g/cm^3] field_drift[V/cm]" << endl;
    return 0;
  }
  
  unsigned long int numEvts = atoi(argv[1]);
  
  string type = argv[2];
  INTERACTION_TYPE type_num;
  if (type == "NR") type_num = NR;
  else if (type == "WIMP")type_num = WIMP;
  else if (type == "B8") type_num = B8;
  else if (type == "DD") type_num = DD;
  else if (type == "AmBe")type_num = AmBe;
  else if (type == "Cf") type_num = Cf;
  else if (type == "ion") type_num = ion;
  else if (type == "gamma")type_num = gammaRay;
  else if (type == "beta")type_num = beta;
  else if (type == "CH3T")type_num = CH3T;
  else type_num = Kr83m;
  
  double eMin = atof(argv[3]);
  double eMax = atof(argv[4]);
  double rho = atof(argv[5]);
  double field = atof(argv[6]);
  
  fprintf(stdout, "E [keV]\tNph\tNe-\n");
  NEST::NESTcalc n;
  for (unsigned long int j = 0; j < numEvts; j++) {
    switch (type_num) {
    case CH3T:
      m_e = 510.9989461; //e- rest mass-energy [keV]
      aa = 0.0072973525664; //fine structure constant
      ZZ = 2.;
      xMax = 18.5898; //tritium beta decay endpoint [keV]
      yMax = 1.1e7; //top of the beta decay E histogram
      xyTry[0] = xMax * n.rand_uniform();
      xyTry[1] = yMax * n.rand_uniform();
      xyTry[2] = 1.;
      while ( xyTry[2] > 0. ) {
	B = sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e) / (xyTry[0] + m_e);
	x = (2.*M_PI*ZZ*aa)*(xyTry[0] + m_e)/sqrt(xyTry[0]*xyTry[0] + 2.*xyTry[0]*m_e);
	FuncValue = (sqrt(2.*xyTry[0]*m_e) *
		(xyTry[0] + m_e) *
		(xMax-xyTry[0]) * (xMax-xyTry[0]) *
		x*(1./(1.-exp(-x)))*(1.002037-0.001427*(B)));
	xyTry = n.VonNeumann(0.,xMax,0.,yMax,xyTry[0],xyTry[1],FuncValue);
      }
      keV = xyTry[0];
      break;
    default:
      keV = eMin + (eMax - eMin) * n.rand_uniform();
      break;
    }
    if (keV > eMax) keV = eMax;
    if (keV < eMin) keV = eMin;
    
    
    NEST::YieldResult yields = n.GetYields(type_num, keV, rho, field);
    cout << keV << "\t" << yields.PhotonYield << "\t" << yields.ElectronYield << endl;
  }
  
  return 1;
  
}
