#ifndef __EXECNEST_H__
#define __EXECNEST_H__ 1

#include "NEST.hh"

using namespace std;
using namespace NEST;


class NESTObservableArray {
  public:
    NESTObservableArray(){};
    ~NESTObservableArray() = default;

    void store_signals(std::vector<double> s1, std::vector<double> s2);

    std::vector<int> s1_nhits;
    std::vector<int> s1_nhits_thr;
    std::vector<int> s1_nhits_dpe;
    std::vector<double> s1r_phe;
    std::vector<double> s1c_phe;
    std::vector<double> s1r_phd;
    std::vector<double> s1c_phd;
    std::vector<double> s1r_spike;
    std::vector<double> s1c_spike;
    std::vector<int> Nee;
    std::vector<int> Nph;
    std::vector<int> s2_nhits;
    std::vector<int> s2_nhits_dpe;
    std::vector<double> s2r_phe;
    std::vector<double> s2c_phe;
    std::vector<double> s2r_phd;
    std::vector<double> s2c_phd;
};

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol, int nFold);

void GetEnergyRes(vector<double> Es);

int execNEST(VDetector* detector, double numEvts, const string& type,
             double eMin, double eMax, double inField, string position,
             const string& posiMuon, double fPos, int seed, bool no_seed,
             double dayNumber);

NESTObservableArray runNESTvec(VDetector* detector,
                               INTERACTION_TYPE scatterType,
                               std::vector<double> eList,
                               std::vector<std::vector<double>> pos3dxyz,
                               double inField = -1.0, int seed = 0,   
                               std::vector<double> ERYieldsParam = default_ERYieldsParam,
                               std::vector<double> NRYieldsParam = default_NRYieldsParam,
                               std::vector<double> NRERWidthsParam = default_NRERWidthsParam,
                               S1CalculationMode s1mode = NEST::S1CalculationMode::Full,
                               S2CalculationMode s2mode = NEST::S2CalculationMode::Full);

#endif
