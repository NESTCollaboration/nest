#ifndef GAMMAHANDLER_HH
#define GAMMAHANDLER_HH

#include <assert.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <random>
#include <vector>
#include "RandomGen.hh"
#include <cmath>
#include <iostream>
#include <string>
//#include "GammaContainer.hh"
using namespace std;

class GammaHandler {
public:
	GammaHandler() {};
	/*
	The main function that combines all spectra from photoionization, compton scattering, pair production, etc...
	Takes min and max energies, a vector of monoenergetic gamma energies and a vector of their branching ratios.
	The index of the gamma energy must correspond with the index of the branch ratio
	*/
	double combineSpectra(double emin, double emax, string source);

	/*
	Get y value for given x value in xyTry. Function is just delta functions at the gammaEnergies with amplitudes given
	by yMax*branchRatio at that energy
	*/
	double photoIonization(vector<vector<double>> sourceInfo, vector<double> xyTry);

	/*Return compton spectrum from KN formula and shifted energies */
	double compton(vector<vector<double>> sourceInfo, vector<double> xyTry);

	/*Return pair production spectrum from PP energy equation */
	double pairProduction(vector<vector<double>> sourceInfo, vector<double> xyTry);

	/*return energies, branching ratios, and mass attenuation coefficients for a given source*/
	vector<vector<double>> sourceLookupTable(string source);

};
 


#endif 