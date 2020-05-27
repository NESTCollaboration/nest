//
// VDetector.hh
//
// Adapted from Quentin Riffard by Jacob Cutter, May 8, 2018

// *********************************************************************
// THIS DEFAULT VIRTUAL DETECTOR SHOULD ONLY BE MODIFIED BY DEVELOPERS.
// PLEASE DEFINE YOUR OWN DETECTOR (see DetectorExample_XENON10.hh).
// *********************************************************************

#ifndef VDetector_hh
#define VDetector_hh 1

#include <vector>

class VDetector {
 public:
  VDetector();
  virtual ~VDetector();
  virtual void Initialization();

  // "Get Functions"
  // Primary Scintillation (S1) parameters
  const double get_g1() { return g1; }
  const double get_sPEres() { return sPEres; }
  const double get_sPEthr() { return sPEthr; }
  const double get_sPEeff() { return sPEeff; }
  const double* get_noiseB() { return &noiseB[0]; }
  const double* get_noiseL() { return &noiseL[0]; }
  const double get_P_dphe() { return P_dphe; }
  
  const bool get_extraPhot(){ return extraPhot; }
  const double get_coinWind() { return coinWind; }
  const int get_coinLevel() { return coinLevel; }
  const int get_numPMTs() { return numPMTs; }

  // Ionization and Secondary Scintillation (S2) parameters
  const double get_g1_gas() { return g1_gas; }
  const double get_s2Fano() { return s2Fano; }
  const double get_s2_thr() { return s2_thr; }
  const double get_E_gas() { return E_gas; }
  const double get_eLife_us() { return eLife_us; }

  // Thermodynamic Properties
  const bool get_inGas() { return inGas; }
  const double get_T_Kelvin() { return T_Kelvin; }
  const double get_p_bar() { return p_bar; }

  // Data Analysis Parameters and Geometry
  const double get_dtCntr() { return dtCntr; }
  const double get_dt_min() { return dt_min; }
  const double get_dt_max() { return dt_max; }
  const double get_radius() { return radius; }
  const double get_radmax() { return radmax; }
  const double get_TopDrift() { return TopDrift; }
  const double get_anode() { return anode; }
  const double get_cathode() { return cathode; }
  const double get_gate() { return gate; }

  // 2-D (X & Y) Position Reconstruction
  const double get_PosResExp() { return PosResExp; }
  const double get_PosResBase() { return PosResBase; }
  
  // Xenon properties
  const double get_molarMass() {return molarMass;}

  // "Set Functions"
  // Primary Scintillation (S1) parameters
  void set_g1(double param) { g1 = param; }
  void set_sPEres(double param) { sPEres = param; }
  void set_sPEthr(double param) { sPEthr = param; }
  void set_sPEeff(double param) { sPEeff = param; }
  void set_noiseB(double p1, double p2, double p3, double p4) {
    noiseB[0] = p1;
    noiseB[1] = p2;
    noiseB[2] = p3;
    noiseB[3] = p4;
  }
  void set_noiseL(double p1, double p2) {
    noiseL[0] = p1;
    noiseL[1] = p2;
  }
  void set_P_dphe(double param) { P_dphe = param; }

  void set_extraPhot(bool param){ extraPhot = param;}
  void set_coinWind(double param) { coinWind = param; }
  void set_coinLevel(int param) { coinLevel = param; }
  void set_numPMTs(int param) { numPMTs = param; }

  // Ionization and Secondary Scintillation (S2) parameters
  void set_g1_gas(double param) { g1_gas = param; }
  void set_s2Fano(double param) { s2Fano = param; }
  void set_s2_thr(double param) { s2_thr = param; }
  void set_E_gas(double param) { E_gas = param; }
  void set_eLife_us(double param) { eLife_us = param; }

  // Thermodynamic Properties
  void set_inGas(bool param) { inGas = param; }
  void set_T_Kelvin(double param) { T_Kelvin = param; }
  void set_p_bar(double param) { p_bar = param; }

  // Data Analysis Parameters and Geometry
  void set_dtCntr(double param) { dtCntr = param; }
  void set_dt_min(double param) { dt_min = param; }
  void set_dt_max(double param) { dt_max = param; }
  void set_radius(double param) { radius = param; }
  void set_radmax(double param) { radmax = param; }
  void set_TopDrift(double param) { TopDrift = param; }
  void set_anode(double param) { anode = param; }
  void set_cathode(double param) { cathode = param; }
  void set_gate(double param) { gate = param; }

  // 2-D (X & Y) Position Reconstruction
  void set_PosResExp(double param) { PosResExp = param; }
  void set_PosResBase(double param) { PosResBase = param; }
  
  //Xenon properties
  void set_molarMass(double param) {molarMass = param;}
  
  typedef enum {
    fold = 0,
    unfold = 1
  } LCE;
  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  virtual double FitS1(double xPos_mm, double yPos_mm, double zPos_mm, LCE map) {
    return 1.;
  }
  
  // Drift electric field as function of Z in mm
  // For example, use a high-order poly spline
  virtual double FitEF(double xPos_mm, double yPos_mm, double zPos_mm) {
    return 730.;
  }

  // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  virtual double FitS2(double xPos_mm, double yPos_mm, LCE map) { return 1.; }

  virtual std::vector<double> FitTBA(double xPos_mm, double yPos_mm,
                                     double zPos_mm) {
    std::vector<double> TopBotAsym;
    return TopBotAsym;
  }

  virtual double OptTrans(double xPos_mm, double yPos_mm, double zPos_mm) {
    return 0.;
  }
  virtual std::vector<double> SinglePEWaveForm(double area, double t0) {
    std::vector<double> PEperBin;
    return PEperBin;
  }

protected:
  // Primary Scintillation (S1) parameters
  double g1 = 0.0760;    // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
  double sPEres = 0.58;  // single phe resolution (Gaussian assumed)
  double sPEthr = 0.35;  // POD threshold in phe, usually used IN PLACE of sPEeff
  double sPEeff = 1.00;  // actual efficiency, can be used in lieu of POD threshold
  double noiseB[4] = {0.0,  // baseline noise mean and width in PE (Gaussian)
                  0.0,  // baseline noise mean and width in PE (Gaussian)
                  0.0,   //EXO noise mean
                  0.0    //EXO noise width
                  };
 
  double P_dphe = 0.2;    // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT
  
  double coinWind = 100;  // S1 coincidence window in ns
  int coinLevel = 2;   // how many PMTs have to fire for an S1 to count
  int numPMTs = 89;    // For coincidence calculation

  bool extraPhot=false;  // for matching EXO-200's W measurement
  //"Linear noise" terms as defined in Dahl thesis and by D. McK
  double noiseL[2] = {3e-2,3e-2}; // S1->S1 Gaussian-smeared w/ noiseL[0]*S1. Ditto S2

  // Ionization and Secondary Scintillation (S2) parameters
  double g1_gas = 0.06;  // phd per S2 photon in gas, used to get SE size
  double s2Fano = 3.61;  // Fano-like fudge factor for SE width
  double s2_thr = 300.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
  double E_gas = 12.;    // field in kV/cm between liquid/gas border and anode
  double eLife_us = 2200.;  // the drift electron mean lifetime in micro-seconds

  // Thermodynamic Properties
  bool inGas = false;
  double T_Kelvin = 177.;  // for liquid drift speed calculation
  double p_bar = 2.14;     // gas pressure in units of bars, it controls S2 size
  // if you are getting warnings about being in gas, lower T and/or raise p

  // Data Analysis Parameters and Geometry
  double dtCntr = 40.;  // center of detector for S1 corrections, in usec.
  double dt_min = 20.;  // minimum. Top of detector fiducial volume
  double dt_max = 60.;  // maximum. Bottom of detector fiducial volume

  double radius = 50.;  // millimeters
  double radmax = 50.;

  double TopDrift = 150.;  // mm not cm or us (but, this *is* where dt=0)
  // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
  // PMTs)
  // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
  // anode!
  double anode = 152.5;  // the level of the anode grid-wire plane in mm
  // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
  // above it
  double gate = 147.5;  // mm. This is where the E-field changes (higher)
  // in gas detectors, the gate is still the gate, but it's where S2 starts
  double cathode = 1.00;  // mm. Defines point below which events are gamma-X

  // 2-D (X & Y) Position Reconstruction
  double PosResExp = 0.015;     // exp increase in pos recon res at hi r, 1/mm
  double PosResBase = 70.8364;  // baseline unc in mm, see NEST.cpp for usage
  
  double molarMass = 131.293; //molar mass, g/mol
};
#endif
