/*
 * File:   TestSpectra.hh
 * Author: brodsky3
 *
 * Created on December 11, 2017, 10:27 AM
 */

#ifndef TESTSPECTRA_HH
#define TESTSPECTRA_HH

#include <assert.h>
#include <float.h>
#include <math.h>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <exception>

#include "ValidityTests.hh"
#include "RandomGen.hh"

#define NEST_AVO \
  6.0221409e+23       // good to keep in sync w/ NEST.hh, can't define twice
#define ATOM_NUM 54.  // ibid.

#define RHO_NAUGHT 0.3  // local DM halo density in [GeV/cm^3]. Lewin & Smith
#define V_SUN 250.2  // +/- 1.4: arXiv:2105.00599, page 12 (V_EARTH 29.8km/s)
#define V_WIMP 238.  // +/- 1.5: page 12 and Table 1
#define V_ESCAPE 544.  // M.C. Smith et al. RAVE Survey

#define NUMBINS_MAX 1000

static constexpr double ElectronRestMassEnergy = 510.9989461;

class TestSpectra {
 public:
  TestSpectra() = default;
  ;  // private so that it cannot be manually called

  struct WIMP_spectrum_prep {
    double base[100] = {1.};
    double exponent[100] = {0.};
    double integral = 0.;
    double xMax = 0.;
    double divisor = 1.;
  };
  WIMP_spectrum_prep wimp_spectrum_prep;

  static double CH3T_spectrum(double emin, double emax);
  static double C14_spectrum(double emin, double emax);
  static double B8_spectrum(double emin, double emax);
  static double AmBe_spectrum(double emin, double emax);
  static double Cf_spectrum(double emin, double emax);
  static double DD_spectrum(double xMin, double xMax, double expFall, double peakFrac, double peakMu, double peakSig);
  static double ppSolar_spectrum(double emin, double emax);
  static double atmNu_spectrum(double emin, double emax);
  static double WIMP_dRate(double ER, double mWimp, double day);
  static WIMP_spectrum_prep WIMP_prep_spectrum(double mass, double eStep,
                                               double day);
  static double WIMP_spectrum(WIMP_spectrum_prep wprep, double mass,
                              double day);
  static const vector<double> Gamma_spectrum(double xMin, double xMax,
                                             string source);

  static double
  ZeplinBackground();  // an example of how to do a better (non-flat) ER
                       // BG spectrum for a WS, from Henrique Araujo
};

#endif /* TESTSPECTRA_HH */
