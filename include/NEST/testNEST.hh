#ifndef __TESTNEST_H__
#define __TESTNEST_H__ 1

using namespace std;
using namespace NEST;

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol, int nFold);

void GetEnergyRes(vector<double> Es);

int testNEST(VDetector* detector, unsigned long int numEvts, string type,
             double eMin, double eMax, double inField, string position, string posiMuon,
             double fPos, int seed, bool no_seed, double dayNumber);
vector<vector<double>> runNESTvec(VDetector* detector, INTERACTION_TYPE scatterType,
	       std::vector<double> eList, std::vector<std::vector<double>> pos3dxyz);

#endif
