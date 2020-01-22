#ifndef __TESTNEST_H__
#define __TESTNEST_H__ 1

using namespace std;
using namespace NEST;

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol);

void GetEnergyRes(vector<double> Es);

int testNEST(VDetector* detector, unsigned long int numEvts, string type,
             double eMin, double eMax, double inField, string position, string posiMuon,
             double fPos, int seed, bool no_seed);

#endif
