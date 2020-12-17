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
	const vector<double> combineSpectra(double emin, double emax, string source);

	/*
	Get y value for given x value in xyTry. Function is just delta functions at the gammaEnergies with amplitudes given
	by yMax*branchRatio at that energy
	*/
	double photoIonization(const vector<vector<double>>& sourceInfo, const vector<double>& xyTry);

  /*Return compton spectrum from KN formula and shifted energies */
	double compton(const vector<vector<double>>& sourceInfo, const vector<double>& xyTry);

  /*Return pair production spectrum from pair production energy equation */
	double pairProduction(const vector<vector<double>>& sourceInfo, const vector<double>& xyTry);
  
  /*return energies, branching ratios, and mass attenuation coefficients for a given source*/

	const vector<vector<double>>& sourceLookupTable(const std::string& source);
};
 



#endif 
