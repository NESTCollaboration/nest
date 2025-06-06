//
// LZ_SR1.hh
//
// Adapted from Quentin Riffard and Jacob Cutter by Greg Rischbieter, June 2022
//
// This file serves as a NEST input to reproduce LZ's SR1 result
//
// Please reach out to Greg Rischbieter (rischbie@umich.edu) and/or
// Matthew Szydagis (mszydagis@albany.edu) with any questions.
//
//////////////////////////////  IMPORTANT !!
/////////////////////////////////////////
//  To match our calibration data, we needed to reduce the ER and NR band
//  widths. To make these changes in your copy of NEST, make sure NEST.cpp is
//  calling an edited NRERWidthsParams vector when using the GetQuanta(...)
//  function: {0.4.,0.4.,0.04,0.5,0.19,2.25, -0.0015, 0.046452, 0.205, 0.45,
//  -0.2};
//
//  If using execNEST to simulate results, these changes can be made manually
//  to execNEST.cpp at lines 183-194 in NESTv2.3.7
//
//  Additonally, we disabled skewed-recombination to better match our tritium
//  data To make this change in v2.3.7 or earlier, replace lines 353-358 of
//  NEST.cpp with "skewness = 0.0;"
//
//  as of NESTv2.3.8, the GetQuanta(...) function takes skewness to be a
//  function input. If using this version of NEST, simply make sure GetQuanta
//  uses "skewnessER = 0.0" for ER simulations (this variable is the last
//  argument to GetQuanta).
///////////////////////////////////////////////////////////////////////////////////
//
//  To add this file into execNEST, make sure it is copied into the directory
//                        nest/include/Detectors/
//  and add the line ' #include "LZ_SR1.hh" ' , and make sure the "auto
//  detector" variable is the "LZ_Detector()" object, as opposed to the LUX
//  default. You may also want to remove the few lines of LUX-detector comments
//  between lines 42-45 of execNEST.cpp
////////////////////////////////////////////////////////////////////////////////////

#ifndef LZ_Detector_hh
#define LZ_Detector_hh 1

#include "VDetector.hh"

using namespace std;

class LZ_Detector : public VDetector {
 public:
  LZ_Detector() {
    // Call the initialisation of all the parameters
    Initialization();
  };
  ~LZ_Detector() override = default;

  // Do here the initialization of all the parameters that are not varying as a
  // function of time
  void Initialization() override {
    name = "LZ SR1";
    // Primary Scintillation (S1) parameters
    g1 = 0.113569;  // 0.114 +/- 0.002  // phd per S1 phot at dtCntr (not phe).
                    // Divide out 2-PE effect
    sPEres = 0.31;  // single phe resolution (Gaussian assumed)
    sPEthr = 0.10;  // POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 1.0;   // actual efficiency, can be used in lieu of POD threshold
    noiseBaseline[0] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noiseBaseline[1] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noiseBaseline[2] = 0.;
    noiseBaseline[3] = 0.;
    P_dphe =
        0.214;  // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

    coinWind = 300;  // S1 coincidence window in ns
    coinLevel = 3;   // how many PMTs have to fire for an S1 to count
    numPMTs = 484;  // Taking into account turned-off PMTs    // For coincidence
                    // calculation

    OldW13eV = true;
    noiseLinear[0] = 0.;
    noiseLinear[1] = 0.;

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.092103545;  // 0.092 +/- 0.002 // phd per S2 photon in gas, used
                           // to get SE size
    s2Fano = 2.0;          // Fano-like fudge factor for SE width
    s2_thr =
        600. *
        (1.0 +
         P_dphe);  // the S2 threshold in phe or PE, *not* phd. Affects NR most
    E_gas = 8.42417;   // field in kV/cm between liquid/gas border and anode
    eLife_us = 6500.;  // the drift electron mean lifetime in micro-seconds

    // Thermodynamic Properties
    inGas = false;
    T_Kelvin = 174.1;  // for liquid drift speed calculation
    p_bar = 1.79;      // gas pressure in units of bars, it controls S2 size
    // if you are getting warnings about being in gas, lower T and/or raise p

    // Data Analysis Parameters and Geometry
    dtCntr = 462.5;  // central correction bin is between 425-500us // center of
                     // detector for S1 corrections, in usec.
    dt_min = 86.;    // minimum. Top of detector fiducial volume
    dt_max = 936.5;  // maximum. Bottom of detector fiducial volume

    radius = 688.;  // millimeters (fiducial radius)
    radmax = 728.;  // actual physical geo. limit

    TopDrift = 1461.;  // mm not cm or us (but, this *is* where dt=0)
    // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
    // PMTs)
    // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
    // anode!
    anode = 1469.;  // the level of the anode grid-wire plane in mm
    // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
    // above it
    gate = 1456.;  // mm. This is where the E-field changes (higher)
    // in gas detectors, the gate is still the gate, but it's where S2 starts
    cathode = 0.;  // mm. Defines point below which events are gamma-X

    // 2-D (X & Y) Position Reconstruction
    // Set these to zero to implement "perfect" position corrections
    // Note: LZ used spatial maps to implement S1 and S2 corrections,
    // but for simplicity with this header file, we'll circumvent
    // the need for corrections entirely.
    PosResExp = 0.0;   // exp increase in pos recon res at hi r, 1/mm
    PosResBase = 0.0;  // baseline unc in mm, see NEST.cpp for usage
  }

  double FitS1(double xPos_mm, double yPos_mm, double zPos_mm,
               LCE map) override {
    return 1.0;
  }

  // Drift electric field as function of Z in mm
  double FitEF(double xPos_mm, double yPos_mm,
               double zPos_mm) override {  // in V/cm
    return 192.;
  }

  double FitS2(double xPos_mm, double yPos_mm, LCE map) override { return 1.0; }

  vector<double> FitTBA(double xPos_mm, double yPos_mm,
                        double zPos_mm) override {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.6;    // S1 bottom-to-total ratio
    BotTotRat[1] = 0.323;  // S2 bottom-to-total ratio, typically only used for
                           // position recon (1-this)

    return BotTotRat;
  }

  // The following functions were not used in LZ's SR1, and are copied from the
  // public NEST
  // file for LUX, just so NEST has them available to prevent errors.
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
