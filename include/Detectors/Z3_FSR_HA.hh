
#ifndef DetectorExample_ZEPLIN_hh
#define DetectorExample_ZEPLIN_hh 1

#include "VDetector.hh"

using namespace std;

class DetectorExample_ZEPLIN : public VDetector {
 public:
  DetectorExample_ZEPLIN() {
    cerr << "*** Detector definition message ***" << endl;
    cerr << "You are currently using ZEPLIN-III First Science Run." << endl
         << endl;

    // Call the initialisation of all the parameters
    Initialization();
  };
  ~DetectorExample_ZEPLIN() override = default;

  // Do here the initialisation of all the parameters that are not varying as a
  // function of time
  virtual void Initialization() override {
    // Primary Scintillation (S1) parameters
    g1 = 0.0714;    // phd per S1 phot at dtCntr (not phe). Error bar +/- 0.0050
    sPEres = 0.4;   // single phe resolution (Gaussian assumed)
    sPEthr = 0.33;  // POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 1.0;   // actual efficiency, can be used in lieu of POD threshold
    noiseBaseline[0] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noiseBaseline[1] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noiseBaseline[2] = 0.0;
    noiseBaseline[3] = 0.0;
    P_dphe = 0.11;  // chance 1 photon makes 2 phe instead of 1 in ETEL PMTs
                    // (LopezParedes, 2018)

    coinWind = 50;  // S1 coincidence window in ns (Note: Z3 used time
                    // difference method in Neves, 2011)
    coinLevel = 3;  // how many PMTs have to fire for an S1 to count
    numPMTs = 31;   // For coincidence calculation using the nFold defined in
                    // previous line

    OldW13eV = true;
    noiseLinear[0] = 0.;
    noiseLinear[1] = 0.;

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.111;  // phd per S2 photon in gas, used to get SE size. ~<g1
                     // common but ZIII unique
    s2Fano = 0.80;   // Fano-like fudge factor for SE width
                     // (causes 5.76phd/e,matching Santos 2011)
    s2_thr = 125.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
                    // (4e x 31.4 phe/e)
    E_gas =
        6.71;  // field in kV/cm between liquid surface and anode: 7-7.1 from HA
    // (7.8 kV/cm reported, later improved with electrostatic sims and lower
    // epsilon)
    eLife_us = 35.;  // the drift electron mean lifetime in micro-seconds
                     // (Santos, 2011)

    // Thermodynamic Properties
    inGas = false;
    T_Kelvin = 174.0;  // for liquid drift speed calculation
    p_bar = 1.66;      // gas pressure in units of bars, it controls S2 size
    // 174.1K, 1.65bar from Henrique: tuning needed for NEST to not do gas

    // Data Analysis Parameters and Geometry
    dtCntr = 0.5;   // center of detector for S1 corrections, in usec.
    dt_min = 0.5;   // minimum. Top of detector fiducial volume
    dt_max = 13.0;  // maximum. Bottom of detector fiducial volume
    // cathode position 14.1 us, but NEST's vD doesn't quite match Z3's (2.65
    // mm/us)

    radius = 150.;  // millimeters (fiducial rad)
    radmax = 192.;  // actual physical geo. limit

    // nominal value 36 mm not 35.75 from Henrique Araujo, ZEPLIN PI. Gas gap =
    // 4->4.25mm
    TopDrift = 35.75;  // mm not cm or us (but, this *is* where dt=0) - NOTE: Z3
                       // HAD NO GATE
    // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
    // PMTs)
    // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
    // anode!
    anode = 40.0;  // the level of the anode grid-wire plane in mm
    // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
    // above it
    gate = 35.75;  // mm. This is where the E-field changes (higher) - NOTE: Z3
                   // HAD NO GATE
    // in gas detectors, the gate is still the gate, but it's where S2 starts
    cathode = 0.;  // mm. Defines point below which events are gamma-X

    // 2-D (X & Y) Position Reconstruction (Solovov 2012)
    PosResExp = 0.015;     // exp increase in pos recon res at hi r, 1/mm
    PosResBase = 70.8364;  // baseline unc in mm, see NEST.cpp for usage
  }

  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  virtual double FitS1(double xPos_mm, double yPos_mm, double zPos_mm,
                       LCE map) override {
    return 1.;  // unitless, 1.000 at detector center
  }

  // Drift electric field as function of Z in mm
  // For example, use a high-order poly spline
  virtual double FitEF(double xPos_mm, double yPos_mm,
                       double zPos_mm) override {  // in V/cm
    return 3840.;  // compromise between 3900 (published) and 3800 (reanalysis)
  }

  // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  virtual double FitS2(double xPos_mm, double yPos_mm, LCE map) override {
    return 1.;  // unitless, 1.000 at detector center
  }

  virtual vector<double> FitTBA(double xPos_mm, double yPos_mm,
                                double zPos_mm) override {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.6;  // S1 bottom-to-total ratio
    BotTotRat[1] = 0.4;  // S2 bottom-to-total ratio, typically only used for
                         // position recon (1-this)

    return BotTotRat;
  }

  virtual double OptTrans(double xPos_mm, double yPos_mm,
                          double zPos_mm) override {
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

  virtual vector<double> SinglePEWaveForm(double area, double t0) override {
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
  // Vary VDetector parameters through custom functions
  virtual void ExampleFunction() { set_g1(0.06); }
};

#endif
