/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TestSpectra.hh
 * Author: brodsky3
 *
 * Created on December 11, 2017, 10:27 AM
 */

#ifndef TESTSPECTRA_HH
#define TESTSPECTRA_HH

#include "NEST.hh"

namespace NEST {
  
  struct WIMP_spectrum_prep {
    double base[100] = {1.};
    double exponent[100] = {0.};
    double integral = 0.;
    double xMax = 0.;
    double divisor = 1.;
  };
  
  double CH3T_spectrum(double emin, double emax, NESTcalc& n);
  double B8_spectrum(double emin, double emax, NESTcalc& n);
  double AmBe_spectrum(double emin, double emax, NESTcalc& n);
  double Cf_spectrum(double emin, double emax, NESTcalc& n);
  double DD_spectrum(double emin, double emax, NESTcalc& n);
  double WIMP_dRate(double ER, double mWimp, NESTcalc &n);
  WIMP_spectrum_prep WIMP_prep_spectrum(double mass, NESTcalc &n);
  double WIMP_spectrum(WIMP_spectrum_prep wprep, double mass, NESTcalc& n);
  
}

#endif /* TESTSPECTRA_HH */
