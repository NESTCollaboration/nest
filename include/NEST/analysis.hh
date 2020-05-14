
// Verbosity flag (for limiting output to yields; no timing)
bool verbosity = true;
// Loop for testNEST and rootNEST to find the best-fit model parameters
unsigned loopNEST = 0;
//0 for no or off, 1 for ER, 2 for NR

// General parameters of importance changing the global behavior
bool MCtruthE = false;    // false means reconstructed energy
bool MCtruthPos = false;  // false means reconstructed position

int useTiming = 0;  // photon arrival times + pulse shapes (2=eTrains)
// if 1 or 2 but verb off, then timing only saved as vectors
// if -1 it means a special extra-fast mode for higher energies

// 0 means PE, 1 means phd (PE/~1.2), 2 means spike count
int usePD = 2;
// band style: log(S2) with 1, while 0 means log(S2/S1)
int useS2 = 0;  // xtra feature: 2 means S2 x-axis energy scale

double minS1 = 1.5; //units are controlled by the usePE flag
// this is separate from S1 thresholds controlled by detector
double maxS1 = 99.5;
int numBins = 98; //for DD, change these to 1.7,110.6,99

// for efficiency calculation
// minS2 need not match S2 threshold in detector.hh
// you can treat as trigger vs. analysis thresholds
double minS2 = 42.;
double maxS2 = 1e4; //5e3 for DD. At least 2e5 for post-Run04 14C

// log(S2/S1) or log(S2) admitted into analysis incl. limit
double logMax = 3.6;
double logMin = 0.6;
int logBins = 50; //#bins in between logMin & logMax for fits

// some numbers for fine-tuning the speed vs. the accuracy
double z_step = 0.1;  // mm, for integrating non-uniform field
double E_step = 5.0;  // keV, for integrating WIMP spectrum
// Rec >~20GeV 6keV, <~5GeV 0.5keV

// Set the rootNEST options
int freeParam= 2; // #free param for calculating DoF in X^2; 2 for Ly and Qy, or g1 and g2
int skewness = 2; // true means skew-Gaussian fits
int mode = 0;
//0 default is to provide the ER BG discrimination & leakage frac
//1 outputs the goodness of fit for one band (Gaussian centroids of histogram in S1 slices)
//2 outputs wimp masses and cross-sections for given efficiency
