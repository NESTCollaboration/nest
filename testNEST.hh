/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   testNEST.cpp
 * Author: brodsky3
 * Modified by ncarrara
 *
 * Created on August 1, 2017, 1:03 PM
 */




using namespace std;
using namespace NEST;

/*
 *
 */


vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol);

void GetEnergyRes(vector<double> Es);

void testNEST(VDetector * detector, unsigned long int numEvts, string type, double eMin, double eMax,
		double inField, std::vector<double> pos, double fPos, int seed, bool no_seed);