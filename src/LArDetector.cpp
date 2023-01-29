/**
 * @file LArDetector.cpp
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief
 * @version
 * @date 2022-04-14
 */
#include "LArDetector.hh"

LArDetector::LArDetector() { Initialization(); }

void LArDetector::Initialization() {
  // Primary Scintillation (S1) parameters
  // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
  double g1 = 0.0760;
  // single phe resolution (Gaussian assumed)
  double sPEres = 0.58;
  // POD threshold in phe, usually used IN PLACE of sPEeff
  double sPEthr = 0.35;
  // actual efficiency, can be used in lieu of POD threshold
  double sPEeff = 1.00;
  double noiseBaseline[4] = {
      0.0,  // baseline noise mean and width in PE (Gaussian)
      0.0,  // baseline noise mean and width in PE (Gaussian)
      0.0,  // EXO noise mean
      0.0   // EXO noise width
  };
  // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT
  double P_dphe = 0.2;

  double coinWind = 100;  // S1 coincidence window in ns
  int coinLevel = 2;      // how many PMTs have to fire for an S1 to count
  int numPMTs = 89;       // For coincidence calculation

  bool OldW13eV = true;  // for matching EXO-200's W measurement
  // "Linear noise" terms as defined in Dahl thesis and by D. McK
  // S1->S1 Gaussian-smeared w/ noiseL[0]*S1. Ditto S2
  double noiseLinear[2] = {3e-2, 3e-2};
  //(n)EXO quadratic noise term
  double noiseQuadratic[2] = {0.0, 0.0};

  // Ionization and Secondary Scintillation (S2) parameters
  double g1_gas = 0.06;  // phd per S2 photon in gas, used to get SE size
  double s2Fano = 3.61;  // Fano-like fudge factor for SE width
  // the S2 threshold in phe or PE, *not* phd. Affects NR most
  double s2_thr = 300.;
  double E_gas = 12.;  // field in kV/cm between liquid/gas border and anode
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
  double molarMass = 131.293;   // molar mass, g/mol
}