
// Verbosity flag (for limiting output to yields; no timing)
bool verbosity = true;

// General parameters of importance changing the global behavior
bool MCtruthE = false;    // false means reconstructed energy
bool MCtruthPos = false;  // false means reconstructed position

int useTiming = 0;  // photon arrival times + pulse shapes (2=eTrains)
// if 1 or 2 but verb off, then timing only saved as vectors
// if -1 it means a special extra-fast mode for higher energies

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
double E_step = 5.0;  // keV, for integrating WIMP spectrum

// Set the rootNEST options
int freeParam = 2; // #free param for calculating DoF in X^2
int mode = 0;
//0 default is to provide the ER BG discrimination & leakage frac
//1 outputs the goodness of fit for one band (Gaussian centroids of histogram in S1 slices)
//2 outputs wimp masses and cross-sections for given efficiency
