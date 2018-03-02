
 // Primary Scintillation (S1) parameters

double g1 = 0.10; //phd per S1 photon in liquid at dtCntr (not phe)
double sPEres = 0.5; //single phe resolution (Gaussian assumed)
double sPEthr = 0.25; //POD threshold in phe
double sPEeff = 0.92; //actual efficiency, can be used in lieu of POD threshold
double noise[2] = {0.0,0.1}; //baseline noise mean and width in PE (Gaussian)
double P_dphe = 0.2; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

int coinLevel = 3; //how many PMTs have to fire for an S1 to count
int numPMTs = 100; //For coincidence calculation

//S1 PDE quartic polynomial for function of z
//s1polA + s1polB*z[cm] + s1polC*z^2+... (QE included, for binomial distribution)
double s1poly[5] = {1.,0.,0.,0.,0.}; // unitless, 1.000 at detector center

//Drift electric field as function of Z in cm
//The coefficients for a quintic poly, in rising order
double efpoly[5] = {1.,0.,0.,0.,0.}; // in V/cm

 // Ionization and Secondary Scintillation (S2) parameters

double g1_gas = 0.10; //phd per S2 photon in gas, used to get SE size
double s2Fano = 3.; //Fano-like fudge factor for SE width
double s2_thr = 250.; //the S2 threshold in phd. Effects NR most
double S2botTotRatio = 0.4; //S2 bottom-to-total ratio, not really used anymore
double E_gas = 10.; //field in kV/cm between liquid/gas border and anode
double eLife_us = 500.; //the drift electron mean lifetime in micro-seconds

 // Thermodynamic Properties

double T_Kelvin = 175.; //for liquid drift speed calculation
double p_bar = 1.5; //gas pressure in units of bars, it controls S2 size

 // Data Analysis Parameters and Geometry

double dtCntr = 250.; //center of detector for S1 corrections, in usec.
double dt_min = 0.00; //minimum. Top of detector fiducial volume
double dt_max = 500.; //maximum. Bottom of detector fiducial volume
double liquidBorder = 544.2198; // mm
double gasGap_cm = 0.5; //EL gap in cm, affecting both field and linear S2 term
