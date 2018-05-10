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


class VDetector {

	public:
		
		VDetector();
		virtual ~VDetector();
		virtual void Initialization();
		virtual void SetTime(double timestamp) { }

		// Primary Scintillation (S1) parameters
		double get_g1(){return g1;}
		double get_sPEres(){return sPEres;}
		double get_sPEthr(){return sPEthr;}
		double get_sPEeff(){return sPEeff;}
		double* get_noise(){return &noise[0];}
		double get_P_dphe(){return P_dphe;}

		int get_coinLevel(){return coinLevel;}
		int get_numPMTs(){return numPMTs;}

		// Ionization and Secondary Scintillation (S2) parameters
		double get_g1_gas(){return g1_gas;}
		double get_s2Fano(){return s2Fano;}
		double get_s2_thr(){return s2_thr;}
		double get_S2botTotRatio(){return S2botTotRatio;}
		double get_E_gas(){return E_gas;}
		double get_eLife_us(){return eLife_us;}

		// Thermodynamic Properties
		bool get_inGas(){return inGas;}
		double get_T_Kelvin(){return T_Kelvin;}
		double get_p_bar(){return p_bar;}

		// Data Analysis Parameters and Geometry
		double get_dtCntr(){return dtCntr;}
		double get_dt_min(){return dt_min;}
		double get_dt_max(){return dt_max;}
		double get_radius(){return radius;}
		double get_TopDrift(){return TopDrift;}
		double get_anode(){return anode;}
		double get_gate(){return gate;}

		// 2-D (X & Y) Position Reconstruction
		double get_PosResExp(){return PosResExp;}
		double get_PosResBase(){return PosResBase;}

		// Functions for setting parameters which are NEST-tuned
		void set_inGas(bool InGasPhase) { inGas = InGasPhase; }

		//S1 PDE custom fit for function of z
		//s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
		virtual double FitS1 ( double xPos_mm, double yPos_mm, double zPos_mm ) { return 1.; }
		
		//Drift electric field as function of Z in mm
		//For example, use a high-order poly spline
		virtual double FitEF ( double xPos_mm, double yPos_mm, double zPos_mm ) { return 730.; }
		
		//S2 PDE custom fit for function of r
		//s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
		virtual double FitS2 ( double xPos_mm, double yPos_mm ) { return 1.; }
		
	
	protected:
		
		// Primary Scintillation (S1) parameters
		int coinLevel, numPMTs;
		double g1, sPEres, sPEthr, sPEeff, P_dphe;  
		double noise[2];

		// Ionization and Secondary Scintillation (S2) parameters
		double g1_gas, s2Fano, s2_thr, S2botTotRatio, E_gas, eLife_us;

		// Thermodynamic Properties
		bool inGas;
		double T_Kelvin, p_bar;

		// Data Analysis Parameters and Geometry
		double dtCntr, dt_min, dt_max, radius, TopDrift, anode, gate;

		// 2-D (X & Y) Position Reconstruction
		double PosResExp, PosResBase;

};


#endif
