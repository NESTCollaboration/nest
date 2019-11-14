//
// VDetector.cpp
//
// Adapted from Quentin Riffard by Jacob Cutter, May 8, 2018

// *********************************************************************
// THIS DEFAULT VIRTUAL DETECTOR SHOULD ONLY BE MODIFIED BY DEVELOPERS.
// PLEASE DEFINE YOUR OWN DETECTOR (see DetectorExample_XENON10.hh).
// *********************************************************************

#include <math.h>

#include "VDetector.hh"

VDetector::VDetector() { Initialization(); }

VDetector::~VDetector() {}

void VDetector::Initialization() {
  // Primary Scintillation (S1) parameters
  g1 = 0.0760;    // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
  sPEres = 0.58;  // single phe resolution (Gaussian assumed)
  sPEthr = 0.35;  // POD threshold in phe, usually used IN PLACE of sPEeff
  sPEeff = 1.00;  // actual efficiency, can be used in lieu of POD threshold
  noise[0] = 0.0;  // baseline noise mean and width in PE (Gaussian)
  noise[1] = 0.0;  // baseline noise mean and width in PE (Gaussian)
  P_dphe = 0.2;    // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

  coinWind = 100;  // S1 coincidence window in ns
  coinLevel = 2;   // how many PMTs have to fire for an S1 to count
  numPMTs = 89;    // For coincidence calculation

  // Ionization and Secondary Scintillation (S2) parameters
  g1_gas = 0.06;  // phd per S2 photon in gas, used to get SE size
  s2Fano = 3.61;  // Fano-like fudge factor for SE width
  s2_thr = 300.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
  E_gas = 12.;    // field in kV/cm between liquid/gas border and anode
  eLife_us = 2200.;  // the drift electron mean lifetime in micro-seconds

  // Thermodynamic Properties
  inGas = false;
  T_Kelvin = 177.;  // for liquid drift speed calculation
  p_bar = 2.14;     // gas pressure in units of bars, it controls S2 size
  // if you are getting warnings about being in gas, lower T and/or raise p

  // Data Analysis Parameters and Geometry
  dtCntr = 40.;  // center of detector for S1 corrections, in usec.
  dt_min = 20.;  // minimum. Top of detector fiducial volume
  dt_max = 60.;  // maximum. Bottom of detector fiducial volume

  radius = 50.;  // millimeters
  radmax = 50.;

  TopDrift = 150.;  // mm not cm or us (but, this *is* where dt=0)
  // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
  // PMTs)
  // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
  // anode!
  anode = 152.5;  // the level of the anode grid-wire plane in mm
  // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
  // above it
  gate = 147.5;  // mm. This is where the E-field changes (higher)
  // in gas detectors, the gate is still the gate, but it's where S2 starts
  cathode = 1.00;  // mm. Defines point below which events are gamma-X

  // 2-D (X & Y) Position Reconstruction
  PosResExp = 0.015;     // exp increase in pos recon res at hi r, 1/mm
  PosResBase = 70.8364;  // baseline unc in mm, see NEST.cpp for usage
}
