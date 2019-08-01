
// Verbosity flag (for limiting output to yields; no timing)
bool verbosity = true;

// General parameters of importance changing global behavior
bool MCtruthE = true;    // false means reconstructed energy
bool MCtruthPos = true;  // false means reconstructed position

int useTiming = 0;  // photon arrival times + pulse shapes (2=eTrains)
// if 1 or 2 but verb off, then timing only saved as vectors

// 0 means PE, 1 means phd (PE/~1.2), 2 means spike count
int usePD = 0;
// band style: log(S2) with 1, while 0 means log(S2/S1)
int useS2 = 0;  // xtra feature: 2 means S2 x-axis energy scale

double minS1 = 0.;  // units are controlled by the usePE flag
// this is separate from S1 thresholds controlled by detector
double maxS1 = 165.;
int numBins = 33;

// for efficiency calculation
// minS2 need not match S2 threshold in detector.hh
// you can treat as trigger vs. analysis thresholds
double minS2 = 0.0;
double maxS2 = 1e9;

// log(S2/S1) or log(S2) admitted into analysis incl. limit
double logMax = 3.6;
double logMin = 0.6;

// some numbers for fine-tuning the speed vs. the accuracy
double z_step = 0.1;  // mm, for integrating non-uniform field
double E_step = 2.0;  // keV, for integrating WIMP spectrum

// Number of free parameters, for calculating DOF, for chi^2
int freeParam = 2;
