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

#include <math.h>
#include <vector>
#include <random>
#include <iostream>
#include <assert.h>
#include <float.h>

#include "RandomGen.hh"

#define W_DEFAULT 13.7
#define NEST_AVO 6.0221409e+23
#define ATOM_NUM 54.
#define MOLAR_MASS 131.293

class TestSpectra {
  
	public:
		TestSpectra(){};  // private so that it cannot be manually called
		
		struct WIMP_spectrum_prep {
			double base[100] = {1.};
			double exponent[100] = {0.};
			double integral = 0.;
			double xMax = 0.;
			double divisor = 1.;
		};
		WIMP_spectrum_prep wimp_spectrum_prep;
		
		double CH3T_spectrum(double emin, double emax);
		double B8_spectrum(double emin, double emax);
		double AmBe_spectrum(double emin, double emax);
		double Cf_spectrum(double emin, double emax);
		double DD_spectrum(double emin, double emax);
		double WIMP_dRate(double ER, double mWimp);
		WIMP_spectrum_prep WIMP_prep_spectrum(double mass, double eStep);
		double WIMP_spectrum(WIMP_spectrum_prep wprep, double mass);
  
};


#endif /* TESTSPECTRA_HH */
