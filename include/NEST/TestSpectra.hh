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
#include <exception>

#include "RandomGen.hh"

#define NEST_AVO \
  6.0221409e+23       // good to keep in sync w/ NEST.hh, can't define twice
#define ATOM_NUM 54.  // ibid.

#define RHO_NAUGHT 0.3  // local DM halo density in [GeV/cm^3]
#define V_SUN 233.
// works out to mean V_EARTH of 245 for LUX Run03; Run04, 230 km/s (arXiv:1705.03380)
#define V_WIMP 220.
#define V_ESCAPE 544.

#define NUMBINS_MAX 1000

class TestSpectra {
 public:
  TestSpectra()= default;;  // private so that it cannot be manually called

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
  static double DD_spectrum(double emin, double emax);
  static double ppSolar_spectrum(double emin, double emax);
  static double atmNu_spectrum(double emin, double emax);
  static double WIMP_dRate(double ER, double mWimp, double day);
  static WIMP_spectrum_prep WIMP_prep_spectrum(double mass, double eStep, double day);
  static double WIMP_spectrum(WIMP_spectrum_prep wprep, double mass, double day);
  double ZeplinBackground();  // an example of how to do a better (non-flat) ER
                              // BG spectrum for a WS, from Henrique Araujo
};

#endif /* TESTSPECTRA_HH */
