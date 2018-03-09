
 // Primary Scintillation (S1) parameters

double g1 = 0.076348; //phd per S1 photon in liquid at dtCntr (not phe)
double sPEres = 0.58; //single phe resolution (Gaussian assumed)
double sPEthr = 0.35; //POD threshold in phe
double sPEeff = 1.00; //actual efficiency, can be used in lieu of POD threshold
double noise[2] = {0.0,0.0}; //baseline noise mean and width in PE (Gaussian)
double P_dphe = 0.2; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

int coinLevel = 2; //how many PMTs have to fire for an S1 to count
int numPMTs = 89; //For coincidence calculation

//S1 PDE quartic polynomial for function of z
//s1polA + s1polB*z[cm] + s1polC*z^2+... (QE included, for binomial distribution)
double s1poly[5] = {1.,0.,0.,0.,0.}; // unitless, 1.000 at detector center

//Drift electric field as function of Z in cm
//The coefficients for a quintic poly, in rising order
double efpoly[6] = {730.,0.,0.,0.,0.,0.}; // in V/cm

 // Ionization and Secondary Scintillation (S2) parameters

double g1_gas = 0.0246; //phd per S2 photon in gas, used to get SE size
double s2Fano = 3.61; //Fano-like fudge factor for SE width
double s2_thr = 300.; //the S2 threshold in phd. Effects NR most
double S2botTotRatio = 0.5; //S2 bottom-to-total ratio, not really used anymore
double E_gas = 12.; //field in kV/cm between liquid/gas border and anode
double eLife_us = 2175.; //the drift electron mean lifetime in micro-seconds

 // Thermodynamic Properties

double T_Kelvin = 177.; //for liquid drift speed calculation
double p_bar = 2.14; //gas pressure in units of bars, it controls S2 size

 // Data Analysis Parameters and Geometry

double dtCntr = 40.; //center of detector for S1 corrections, in usec.
double dt_min = 20.; //minimum. Top of detector fiducial volume
double dt_max = 60.; //maximum. Bottom of detector fiducial volume
double liquidBorder = 150.; // mm
double gasGap_cm = 0.25; //EL gap in cm, affecting both field and linear S2 term
