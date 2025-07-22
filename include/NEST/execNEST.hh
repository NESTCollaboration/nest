#ifndef __EXECNEST_H__
#define __EXECNEST_H__ 1

#include "NEST.hh"

using namespace std;
using namespace NEST;


class NESTObservableArray {
  public:
    NESTObservableArray(){};
    ~NESTObservableArray() = default;

    void store_signals(
      double energy,
      std::vector<double> pos,
      NESTresult result,
      std::vector<double> s1,
      std::vector<double> s2,
      std::vector<int64_t> s1_wf_time,
      std::vector<double> s1_wf_amp,
      std::vector<int64_t> s2_wf_time,
      std::vector<double> s2_wf_amp
    );

    // Truth 
    std::vector<double> energy_kev;
    std::vector<double> x_mm;
    std::vector<double> y_mm;
    std::vector<double> z_mm;

    // Quanta
    std::vector<int> n_electrons;
    std::vector<int> n_photons;

    // S1
    std::vector<int> s1_nhits; // MC-true integer hits in same OR different PMTs, NO double phe effect
    std::vector<int> s1_nhits_dpe; // MC-true integer hits WITH double phe effect (Nphe > nHits) 
    std::vector<int> s1_nhits_thr; // nHits post coincidence window and single PE eff! 
    std::vector<double> s1r_phe; // smeared DAQ pulse areas in phe, NO XYZ correction   
    std::vector<double> s1c_phe; // same as s1r_phe WITH XYZ correction 
    std::vector<double> s1r_phd; // same as s1r_phe, adjusted/corrected *downward* for 2-PE effect (LUX phd units)   
    std::vector<double> s1c_phd; // same as s1r_phd, but XYZ-corrected   
    std::vector<double> s1r_spike; // spike count, NO XYZ correction   
    std::vector<double> s1c_spike; // spike count, WITH XYZ correction   

    std::vector<std::vector<double>> s1_photon_times; // Time that S1 photons reach a PMT (includes undetected photons)

    // S2
    std::vector<int> s2_Nee; // integer number of electrons unabsorbed in liquid then getting extracted
    std::vector<int> s2_Nph; // raw number of photons produced in the gas gap
    std::vector<int> s2_nhits; // MC-true integer hits in same OR different PMTs, NO double phe effect
    std::vector<int> s2_nhits_dpe; // MC-true integer hits WITH double phe effect (Nphe > nHits)
    std::vector<double> s2r_phe; // smeared DAQ pulse areas in phe, NO XYZ correction
    std::vector<double> s2c_phe; // same as s2r_phe WITH XYZ correction
    std::vector<double> s2r_phd; // same s2r_phe, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
    std::vector<double> s2c_phd; // same as s2r_phd, but XYZ-corrected

    // Waveforms
    std::vector<std::vector<int64_t>> s1_waveform_time;
    std::vector<std::vector<double>> s1_waveform_amp;
    std::vector<std::vector<int64_t>> s2_waveform_time;
    std::vector<std::vector<double>> s2_waveform_amp;
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
