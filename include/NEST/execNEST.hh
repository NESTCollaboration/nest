#ifndef __EXECNEST_H__
#define __EXECNEST_H__ 1

using namespace std;
using namespace NEST;

struct NESTObservableArray {
  vector<int> s1_nhits;
  vector<int> s1_nhits_thr;
  vector<int> s1_nhits_dpe;
  vector<double> s1r_phe;
  vector<double> s1c_phe;
  vector<double> s1r_phd;
  vector<double> s1c_phd;
  vector<double> s1r_spike;
  vector<double> s1c_spike;
  vector<int> Nee;
  vector<int> Nph;
  vector<int> s2_nhits;
  vector<int> s2_nhits_dpe;
  vector<double> s2r_phe;
  vector<double> s2c_phe;
  vector<double> s2r_phd;
  vector<double> s2c_phd;
};

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol, int nFold);

void GetEnergyRes(vector<double> Es);

int execNEST(VDetector* detector, uint64_t numEvts, const string& type,
             double eMin, double eMax, double inField, string position,
             const string& posiMuon, double fPos, int seed, bool no_seed,
             double dayNumber);
NESTObservableArray runNESTvec(VDetector* detector,
                               INTERACTION_TYPE scatterType,
                               std::vector<double> eList,
                               std::vector<std::vector<double>> pos3dxyz,
                               double inField = -1.0, int seed = 0);

#endif
