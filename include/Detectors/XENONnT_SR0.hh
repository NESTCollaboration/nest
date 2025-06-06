// XENONnT.hh
// Min Zhong (mizhong@ucsd.edu), Sep. 2022
// NEST detector file for XENONnT SR0 (arXiv:2207.11330)

#ifndef XENONnT_hh
#define XENONnT_hh 1

#include "VDetector.hh"

using namespace std;

class XENONnT : public VDetector {
 public:
  XENONnT() {
    //cerr << "*** Detector definition message ***" << endl;
    //cerr << "You are using the XENONnT detector file for SR0." << endl << endl;

    // Call the initialization of all the parameters
    Initialization();
  };
  ~XENONnT() override = default;
  ;

  // Do here the initialization of all the parameters that are not varying as a
  // function of time
  void Initialization() override {
    name = "XENONnT SR0";
    // From projected WIMP sensitivity paper: 2007.08796
    // Primary Scintillation (S1) parameters
    g1 = 0.12520;  // XENONnT SR0: g1 = 0.15149 PE/photon (0.12520 phd/photon)

    // 1T parameters unless otherwise noted.
    sPEres = 0.4;    // single phe resolution (Gaussian assumed)
    sPEthr = 0.;     // POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 0.935;  // Averaging top and bottom PMTs' sPE acceptance
    noiseBaseline[0] = 0.0;  // baseline noise mean in PE (Gaussian)
                             // Assumed to be small
    noiseBaseline[1] = 0.0;  // baseline noise width in PE (Gaussian)
                             // Assumed to be small
    noiseBaseline[2] = 0.0;  // baseline noise mean in e- (for grid wires)
                             // Assumed to be small
    noiseBaseline[3] = 0.0;  // baseline noise width in e- (for grid wires)
                             // Assumed to be small

    P_dphe = 0.21;  // XENONnT SR0 Value

    coinWind = 100;      // ns, unsure but seems right.
    coinLevel = 3;       // how many PMTs have to fire for an S1 to count
    numPMTs = 494 - 17;  // 17 pmts are off/abandoned for SR0

    OldW13eV = true;  // default true, which means use "classic" W instead of
                      // Baudis / EXO's
    noiseLinear[0] = 0.0e-2;  // 1910.04211 p.12, to match 1610.02076 Fig. 8.
                              // UNITS: fraction NOT %!
    noiseLinear[1] = 0.0e-2;  // 1910.04211 p.12, to match 1610.02076 Fig. 8.
                              // UNITS: fraction NOT %!

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.1533;  // phd per S2 photon in gas, used to get SE size
                      // Modified to have SE as 31.1518/1.21
    s2Fano = 1.0;     // Fano-like fudge factor for SE width				It
                   // will scale up s2 resolution. set as 1 for now.
    s2_thr = 100.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
    E_gas = 6.8903;   // field in kV/cm between liquid/gas border and anode
                      // Extraction efficiency: 52.8%
    eLife_us = 15e3;  // the drift electron mean lifetime in micro-seconds
                      // XENONnT Estimation

    // Thermodynamic Properties
    inGas = false;
    T_Kelvin = 176.7;  // for liquid drift speed calculation
    p_bar = 1.890;     // gas pressure in units of bars, it controls S2 size

    // if you are getting warnings about being in gas, lower T and/or raise p

    // Data Analysis Parameters and Geometry
    // If we're varying drift velocity, should we also make these functions?
    dt_min = 213.19;    // minimum. Top of detector fiducial volume
                        // SR0: 13.346 cm / v_drift
    dt_max = 2140.00;   // maximum. Bottom of detector fiducial volume
                        // SR0: 133.970 cm / v_drift
    dtCntr = 1176.645;  // center of detector for S1 corrections, in usec.
                        // Drift velocity: 0.626 mm/us at 23 V/cm

    radius = 607.3;  // millimeters (fiducial rad)
                     // SR0 fiducial volume cut
    radmax = 664.;   // actual physical geo. limit
                     // TPC radius

    // Electrodes positions in mm above 0 at bottom, + above
    // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
    // PMTs) In 2-phase, TopDrift=liquid/gas border. nT dimensions are from
    // https://xe1t-wiki.lngs.infn.it/doku.php?id=xenon:xenonnt:analysis:coordinate_system
    TopDrift = 1550.8;  // mm not cm or us (but, this *is* where dt=0)
                        // 5.0 mm liquid level
    anode = 1553.8;     // the level of the anode grid-wire plane in mm 3.0 mm
                        // above the interface
    gate = 1545.8;      // mm. This is where the E-field changes (higher)
    cathode = 60.2;     // mm. Defines point below which events are gamma-X

    // 2-D (X & Y) Position Reconstruction
    // Set these to zero to implement "perfect" position corrections
    PosResExp = 0.0;   // exp increase in pos recon res at hi r, 1/mm
    PosResBase = 0.0;  // baseline unc in mm, see NEST.cpp for usage
  }

  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  double FitS1(double xPos_mm, double yPos_mm, double zPos_mm,
               LCE map) override {
    return 1.0;
  }

  // Drift electric field as function of Z in mm
  // For example, use a high-order poly spline
  // We don't use this I don't think.
  double FitEF(double xPos_mm, double yPos_mm,
               double zPos_mm) override {  // in V/cm
    return 23.;
  }

  // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  double FitS2(double xPos_mm, double yPos_mm, LCE map) override {
    return 1.;  // no radial corrections for S2-bottom needed really
  }

  vector<double> FitTBA(double xPos_mm, double yPos_mm,
                        double zPos_mm) override {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.8;  // S1 bottom-to-total ratio - not that important here,
                         // rough estimation
    BotTotRat[1] =
        0.252;  // S2 bottom-to-total ratio (important when looking at cS2b!)

    return BotTotRat;
  }

  double OptTrans(double xPos_mm, double yPos_mm, double zPos_mm) override {
    double phoTravT, approxCenter = (TopDrift + cathode) / 2.,
                     relativeZ = zPos_mm - approxCenter;

    double A = 0.048467 - 7.6386e-6 * relativeZ +
               1.2016e-6 * pow(relativeZ, 2.) - 6.0833e-9 * pow(relativeZ, 3.);
    if (A < 0.) A = 0.;  // cannot have negative probability
    double B_a = 0.99373 + 0.0010309 * relativeZ -
                 2.5788e-6 * pow(relativeZ, 2.) -
                 1.2000e-8 * pow(relativeZ, 3.);
    double B_b = 1. - B_a;
    double tau_a = 11.15;  // all times in nanoseconds
    double tau_b = 4.5093 + 0.03437 * relativeZ -
                   0.00018406 * pow(relativeZ, 2.) -
                   1.6383e-6 * pow(relativeZ, 3.);
    if (tau_b < 0.) tau_b = 0.;  // cannot have negative time

    // A = 0.0574; B_a = 1.062; tau_a = 11.1; tau_b = 2.70; B_b = 1.0 - B_a;
    // //LUX D-D conditions

    if (RandomGen::rndm()->rand_uniform() < A)
      phoTravT = 0.;  // direct travel time to PMTs (low)
    else {            // using P0(t) =
      // A*delta(t)+(1-A)*[(B_a/tau_a)e^(-t/tau_a)+(B_b/tau_b)e^(-t/tau_b)]
      // LUX PSD paper, but should apply to all detectors w/ diff #'s
      if (RandomGen::rndm()->rand_uniform() < B_a)
        phoTravT = -tau_a * log(RandomGen::rndm()->rand_uniform());
      else
        phoTravT = -tau_b * log(RandomGen::rndm()->rand_uniform());
    }

    double sig = RandomGen::rndm()->rand_gauss(
        3.84, .09);  // includes stat unc but not syst
    phoTravT += RandomGen::rndm()->rand_gauss(
        0.00, sig);  // the overall width added to photon time spectra by the
                     // effects in the electronics and the data reduction
                     // pipeline

    if (phoTravT > DBL_MAX) phoTravT = tau_a;
    if (phoTravT < -DBL_MAX) phoTravT = 0.000;

    return phoTravT;  // this function follows LUX (arXiv:1802.06162) not Xe10
                      // technically but tried to make general
  }

  vector<double> SinglePEWaveForm(double area, double t0) override {
    vector<double> PEperBin;

    double threshold = PULSEHEIGHT;  // photo-electrons
    double sigma = PULSE_WIDTH;      // ns
    area *= 10. * (1. + threshold);
    double amplitude = area / (sigma * sqrt(2. * M_PI)),
           signal;  // assumes perfect Gaussian

    double tStep1 = SAMPLE_SIZE / 1e2;  // ns, make sure much smaller than
                                        // sample size; used to generate MC-true
                                        // pulses essentially
    double tStep2 =
        SAMPLE_SIZE;  // ns; 1 over digitization rate, 100 MHz assumed here

    double time = -5. * sigma;
    bool digitizeMe = false;
    while (true) {
      signal = amplitude * exp(-pow(time, 2.) / (2. * sigma * sigma));
      if (signal < threshold) {
        if (digitizeMe)
          break;
        else
          ;  // do nothing - goes down to advancing time block
      } else {
        if (digitizeMe)
          PEperBin.push_back(signal);
        else {
          if (RandomGen::rndm()->rand_uniform() < 2. * (tStep1 / tStep2)) {
            PEperBin.push_back(time + t0);
            PEperBin.push_back(signal);
            digitizeMe = true;
          } else {
          }
        }
      }
      if (digitizeMe)
        time += tStep2;
      else
        time += tStep1;
      if (time > 5. * sigma) break;
    }

    return PEperBin;
  }
};

#endif
