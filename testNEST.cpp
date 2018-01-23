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
#include <iostream>

using namespace std;
using namespace NEST;

/*
 * 
 */


int main(int argc, char** argv) {



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
    else if (type == "WIMP")type_num = WIMP;
    else if (type == "B8") {
        type_num = B8;
        numEvts = (unsigned long int) (floor(atof(argv[1])*1.1 / 365. + 0.5));
    } else if (type == "DD") type_num = DD;
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
    
    if ( type_num == Kr83m && eMin == 9.4 && eMax == 9.4 )
      fprintf(stdout, "t [ns]\t\tE [keV]\t\tNph\t\tNe-\n");
    else
      fprintf(stdout, "E [keV]\t\tNph\t\tNe-\n");
    NEST::NESTcalc n;
    if (argc >= 8) n.SetRandomSeed(atoi(argv[7]));
    
    double kev = -999;
    for (unsigned long int j = 0; j < numEvts; j++) {
        if (eMin == eMax) {
            kev = eMin;
        } else {
            switch (type_num) {
                case CH3T:
                    kev = CH3T_spectrum(eMin, eMax, n);
                    break;
                case B8: //normalize this to ~3500 / 10-ton / year, for E-threshold of 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
                    kev = B8_spectrum(eMin, eMax, n);
                    break;
                case AmBe: //for ZEPLIN-III FSR from HA (Pal '98)
                    kev = AmBe_spectrum(eMin, eMax, n);
                    break;
                case Cf:
                    kev = Cf_spectrum(eMin, eMax, n);
                    break;
                case DD:
                    kev = DD_spectrum(eMin, eMax, n);
                    break;
                default:
                    kev = eMin + (eMax - eMin) * n.rand_uniform();
                    break;
            }
        }

        if (kev > eMax) kev = eMax;
        if (kev < eMin) kev = eMin;


        NEST::YieldResult yields = n.GetYields(type_num, kev, rho, field);
        printf("%.6f\t%.6f\t%.6f\n", kev, yields.PhotonYield, yields.ElectronYield);
    }

    return 1;

}
