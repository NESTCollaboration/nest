//
// DetectorExample_XENON10.hh
//
// Adapted from Quentin Riffard by Jacob Cutter, May 8, 2018
//
// This header file serves as a template for creating one's own VDetector class.
// This will be ultimately used to customize the NEST detector parameters to meet
// individual/collaboration needs.
//
// Note that the detector parameters can also be varied throughout a run, etc.

#ifndef DetectorExample_XENON10_hh
#define DetectorExample_XENON10_hh 1

#include "VDetector.hh"

using namespace std;


class DetectorExample_XENON10: public VDetector {

public:

	DetectorExample_XENON10() {
		cerr << "*** Detector definition message ***" << endl;
		cerr << "You are currently using the default XENON10 template detector." << endl << endl;

		// Call the initialisation of all the parameters
		Initialization();
	};
	virtual ~DetectorExample_XENON10() {};

	// Do here the initialization of all the parameters that are not varying as a function of time
	virtual void Initialization() {
		
		// Primary Scintillation (S1) parameters
		g1 = 0.0760; //phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
		sPEres = 0.58; //single phe resolution (Gaussian assumed)
		sPEthr = 0.35; //POD threshold in phe, usually used IN PLACE of sPEeff
		sPEeff = 1.00; //actual efficiency, can be used in lieu of POD threshold
		noise[0] = 0.0; //baseline noise mean and width in PE (Gaussian)
		noise[1] = 0.0; //baseline noise mean and width in PE (Gaussian)
		P_dphe = 0.2; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

		coinLevel= 2; //how many PMTs have to fire for an S1 to count
		numPMTs = 89; //For coincidence calculation

		// Ionization and Secondary Scintillation (S2) parameters
		g1_gas = 0.06; //phd per S2 photon in gas, used to get SE size
		s2Fano = 3.61; //Fano-like fudge factor for SE width
		s2_thr = 300.; //the S2 threshold in phe or PE, *not* phd. Affects NR most
		S2botTotRatio = 0.4; //bottom-to-total, typically only used for position recon (1-this)
		E_gas = 12.; //field in kV/cm between liquid/gas border and anode
		eLife_us = 2200.; //the drift electron mean lifetime in micro-seconds

		// Thermodynamic Properties
		inGas = false;
		T_Kelvin = 177.; //for liquid drift speed calculation
		p_bar = 2.14; //gas pressure in units of bars, it controls S2 size
		//if you are getting warnings about being in gas, lower T and/or raise p

		// Data Analysis Parameters and Geometry
		dtCntr = 40.; //center of detector for S1 corrections, in usec.
		dt_min = 20.; //minimum. Top of detector fiducial volume
		dt_max = 60.; //maximum. Bottom of detector fiducial volume

		radius = 50.; //millimeters

		TopDrift = 150.; //mm not cm or us (but, this *is* where dt=0)
		//a z-axis value of 0 means the bottom of the detector (cathode OR bottom PMTs)
		//In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not anode!
		anode = 152.5; //the level of the anode grid-wire plane in mm
		//In a gas TPC, this is not TopDrift (top of drift region), but a few mm above it
		gate = 147.5; //mm. This is where the E-field changes (higher)
		// in gas detectors, the gate is still the gate, but it's where S2 starts

		// 2-D (X & Y) Position Reconstruction
		PosResExp = 0.015; // exp increase in pos recon res at hi r, 1/mm
		PosResBase = 70.8364; // baseline unc in mm, see NEST.cpp for usage
	}

	//S1 PDE custom fit for function of z
	//s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
	virtual double FitS1 ( double xPos_mm, double yPos_mm, double zPos_mm ) {
		return 1.; // unitless, 1.000 at detector center
	}
		
	//Drift electric field as function of Z in mm
	//For example, use a high-order poly spline
	virtual double FitEF ( double xPos_mm, double yPos_mm, double zPos_mm ) { // in V/cm
		return 730.;
}

	//S2 PDE custom fit for function of r
	//s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
	virtual double FitS2 ( double xPos_mm, double yPos_mm ) {
		return 1.; // unitless, 1.000 at detector center
	}
	
	// Vary parameters as necessary based on the timestamp of the event, or any other
	// custom dependencies. Any protected parameters from VDetector can be modified here.
	virtual void SetTime(double timestamp) {
		Modify_g1(timestamp);
	}

private:

	// Parameter modification functions
	void Modify_g1(double timestamp) {
		// if timestamp falls in certain range...
		set_g1(0.0760);
	}

};


#endif
