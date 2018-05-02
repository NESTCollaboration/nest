
//General parameters of importance changing global behavior
bool MCtruthE = true; //false means reconstructed energy
bool MCtruthPos=true; //false means reconstructed position

bool useTiming = false; //photon arrival times + pulse shapes

//0 means PE, 1 means phd (PE/~1.2), 2 means spike count
int usePE = 0;
//band style: log(S2) with 1, while 0 means log(S2/S1)
int useS2 = 0; // xtra feature: 2 means S2 x-axis energy scale

double minS1 = 0.; //units are controlled by the usePE flag
double maxS1 = 165.;
int numBins = 33;

//for efficiency calculation
//minS2 need not match S2 threshold in detector.hh
//you can treat as trigger vs. analysis thresholds
double minS2 = 0.0;
double maxS2 = 1e9;

//some numbers for fine-tuning the speed vs. the accuracy
double z_step = 0.1; //mm, for integrating non-uniform field
double E_step = 2.0; //keV, for integrating WIMP spectrum
