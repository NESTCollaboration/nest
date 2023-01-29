/**
 * @file TestLArSpectra.hh
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief
 * @version
 * @date 2022-04-14
 */
#pragma once

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

class TestLArSpectra {
 public:
  TestLArSpectra() = default;
  // struct WIMP_spectrum_prep {
  //     double base[100] = {1.};
  //     double exponent[100] = {0.};
  //     double integral = 0.;
  //     double xMax = 0.;
  //     double divisor = 1.;
  // };
  // WIMP_spectrum_prep wimp_spectrum_prep;

  // static double CH3T_spectrum(double emin, double emax);
  // static double C14_spectrum(double emin, double emax);
  // static double B8_spectrum(double emin, double emax);
  // static double AmBe_spectrum(double emin, double emax);
  // static double Cf_spectrum(double emin, double emax);
  // static double DD_spectrum(double xMin, double xMax, double expFall, double
  // peakFrac, double peakMu, double peakSig); static double
  // ppSolar_spectrum(double emin, double emax); static double
  // atmNu_spectrum(double emin, double emax); static double WIMP_dRate(double
  // ER, double mWimp, double day); static WIMP_spectrum_prep
  // WIMP_prep_spectrum(double mass, double eStep,
  //                                             double day);
  // static double WIMP_spectrum(WIMP_spectrum_prep wprep, double mass,
  //                             double day);
  // static const vector<double> Gamma_spectrum(double xMin, double xMax,
  //                                             string source);

  // static double
  // ZeplinBackground();  // an example of how to do a better (non-flat) ER
  //                    // BG spectrum for a WS, from Henrique Araujo
};