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

using namespace std;

/*
 * 
 */

int main(int argc, char** argv) {

    

    string type = argv[1];
    int type_num;
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

    double keV = atof(argv[2]);
    double rho = atof(argv[3]);
    double field = atof(argv[4]);
    NEST::NESTcalc n;
    NEST::YieldResult yields = n.GetYields(type_num, keV, rho, field);
    cout <<"Photon Yield: "<< yields.PhotonYield << "\tElectron Yield: " << yields.ElectronYield << endl;

    return 1;

}