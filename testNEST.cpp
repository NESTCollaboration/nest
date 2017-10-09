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
#include <stdlib.h>

using namespace std;
using namespace NEST;

/*
 * 
 */


int main(int argc, char** argv) {

    double keV;

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
    else if (type == "Kr83m") type_num = Kr83m;
    else {cout << "Invalid particle type.\n"; return 2;}

    double eMin = atof(argv[3]);
    double eMax = atof(argv[4]);
    double rho = atof(argv[5]);
    double field = atof(argv[6]);

    fprintf(stdout, "E [keV]\tNph\tNe-\n");
    NEST::NESTcalc n;
    n.rng.seed(std::random_device()());

    for (unsigned long int j = 0; j < numEvts; j++) {
        switch (type_num) {
            default:
                keV = eMin + (eMax - eMin) * n.rand_uniform();
                break;
        }
        if (keV > eMax) keV = eMax;
        if (keV < eMin) keV = eMin;

        
        NEST::YieldResult yields = n.GetYields(type_num, keV, rho, field);
        std::cout << "Photon Yield: " << yields.PhotonYield << "\tElectron Yield: " << yields.ElectronYield << std::endl;
    }

    return 1;

}
