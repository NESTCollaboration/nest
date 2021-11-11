#include "NEST.hh"

// Verbosity flag (for limiting output to yields; no timing)
bool verbosity = true;
// Loop for execNEST and rootNEST to find the best-fit model parameters
unsigned loopNEST = 0;
// 0 for no or off, 1 for ER, 2 for NR
// an executable bash script was created in the build directory
bool PrintSubThr = true;
// include or exclude the infamous negative S1 & S2 pulses (and 0's)

// General parameters of importance changing the global behavior
bool MCtruthE = false;    // false means reconstructed energy
bool MCtruthPos = false;  // false means reconstructed position

// Setting the S1 and S2 calculation modes
NEST::S1CalculationMode s1CalculationMode = NEST::S1CalculationMode::Hybrid;
// S1 calculation mode options are:
// Full [Default]: calculating the pulse area by looping over all the pmt hits
// Parametric: calculating the pulse area by using a parametric equation
// Hybrid: Using Full and Parametric with a transition point: n_pmt_hits >
// n_pmts Waveform: calculating the pulse area with the Full calculation mode
// and the waveform

NEST::S2CalculationMode s2CalculationMode = NEST::S2CalculationMode::Full;
// S2 calculation mode options are:
// Full [Default]: calculate only the pulse area
// Waveform: calculate the pulse area and the waveform
// WaveformWithEtrain: calculate the pulse area and the waveform with etrain

// 0 means PE, 1 means phd (PE/~1.2), 2 means spike count
int usePD = 2;
// band style: log(S2) with 1, while 0 means log(S2/S1)
int useS2 = 0;  // xtra feature: 2 means S2 x-axis energy scale

double minS1 = 1.5;  // units are controlled by the usePE flag
// this is separate from S1 thresholds controlled by detector
double maxS1 = 99.5;
int numBins = 98;  // for LUXRun03 DD, change these to 1.7,110.6,99

// for efficiency calculation
// minS2 need not match S2 threshold in detector.hh
// you can treat as trigger vs. analysis thresholds
double minS2 = 42.;
double maxS2 = 1e4;  // 5e3 for DD. At least 2e5 for post-Run04 14C

// log(S2/S1) or log(S2) admitted into analysis incl. limit
double logMax = 3.6;  // when skewness=1 or 2 ROOT ignores these and does raw
                      // mean +/- 3-sigma
double logMin = 0.6;
int logBins = 30;  //#bins in between logMin & logMax for fits

// some numbers for fine-tuning the speed vs. the accuracy
double z_step =
    0.1;  // mm, for integrating non-uni EF. Larger detectors require larger
          // z_step ~0.5-1mm for tonne-scale TPCs like LZ and XENONnT
double E_step = 5.0;  // keV, for integrating WIMP spectrum. NEST will warn you
                      // if choice poor
// Rec >~20GeV 6keV, <~5GeV 0.5keV

// Set the rootNEST options
int freeParam =
    2;  // #free param for calculating DoF in X^2; 2 for Ly and Qy, or g1 and g2
int skewness =
    1;  // 1 means skew-Gaussian fits (2 more rigorous fit, more detail output)
int mode = 0;
// 0 default is to provide 1 band (no data comp) or if 2 args ER BG discrim &
// leakage frac 1 outputs GoF for sim band 1st cf. data band 2nd (Gauss centroids
// of hist in S1 slices) 2 outputs wimp masses and cross-sections for given
// efficiency
int NRbandCenter =
    3;  // 1 is Gauss mean, 3 is fit to means, 2 is compromise. Neg->medians
