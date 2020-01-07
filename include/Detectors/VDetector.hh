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
  double get_g1() { return g1; }
  double get_sPEres() { return sPEres; }
  double get_sPEthr() { return sPEthr; }
  double get_sPEeff() { return sPEeff; }
  double* get_noise() { return &noise[0]; }
  double get_P_dphe() { return P_dphe; }
  
  bool get_extraPhot(){ return extraPhot; }
  double get_coinWind() { return coinWind; }
  int get_coinLevel() { return coinLevel; }
  int get_numPMTs() { return numPMTs; }

  // Ionization and Secondary Scintillation (S2) parameters
  double get_g1_gas() { return g1_gas; }
  double get_s2Fano() { return s2Fano; }
  double get_s2_thr() { return s2_thr; }
  double get_E_gas() { return E_gas; }
  double get_eLife_us() { return eLife_us; }

  // Thermodynamic Properties
  bool get_inGas() { return inGas; }
  double get_T_Kelvin() { return T_Kelvin; }
  double get_p_bar() { return p_bar; }

  // Data Analysis Parameters and Geometry
  double get_dtCntr() { return dtCntr; }
  double get_dt_min() { return dt_min; }
  double get_dt_max() { return dt_max; }
  double get_radius() { return radius; }
  double get_radmax() { return radmax; }
  double get_TopDrift() { return TopDrift; }
  double get_anode() { return anode; }
  double get_cathode() { return cathode; }
  double get_gate() { return gate; }

  // 2-D (X & Y) Position Reconstruction
  double get_PosResExp() { return PosResExp; }
  double get_PosResBase() { return PosResBase; }
  
  // Xenon properties
  double get_molarMass() {return molarMass;}

  // "Set Functions"
  // Primary Scintillation (S1) parameters
  void set_g1(double param) { g1 = param; }
  void set_sPEres(double param) { sPEres = param; }
  void set_sPEthr(double param) { sPEthr = param; }
  void set_sPEeff(double param) { sPEeff = param; }
  void set_noise(double p1, double p2, double p3, double p4) {
    noise[0] = p1;
    noise[1] = p2;
    noise[2] = p3;
    noise[3] = p4;
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

  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  virtual double FitS1(double xPos_mm, double yPos_mm, double zPos_mm) {
    return 1.;
  }

  // Drift electric field as function of Z in mm
  // For example, use a high-order poly spline
  virtual double FitEF(double xPos_mm, double yPos_mm, double zPos_mm) {
    return 730.;
  }

  // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  virtual double FitS2(double xPos_mm, double yPos_mm) { return 1.; }

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
  int coinLevel, numPMTs;
  double g1, sPEres, sPEthr, sPEeff, P_dphe, coinWind;
  double noise[4];

  // Ionization and Secondary Scintillation (S2) parameters
  double g1_gas, s2Fano, s2_thr, E_gas, eLife_us;

  // Thermodynamic Properties
  bool inGas, extraPhot;
  double T_Kelvin, p_bar;

  // Data Analysis Parameters and Geometry
  double dtCntr, dt_min, dt_max, radius, radmax, TopDrift, anode, cathode, gate;

  // 2-D (X & Y) Position Reconstruction
  double PosResExp, PosResBase;
  
  //Xenon properties
  double molarMass;
};

#endif
