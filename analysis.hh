
bool MCtruthE = true; //false means reconstructed energy

//0 means PE, 1 means phd (PE/~1.2), 2 means spike count
int usePE = 0;
//band style: log(S2) with true, while false means log(S2/S1)
bool useS2 = false;

double minS1 = 0.; //units are controlled by the usePE flag
double maxS1 = 165.;
int numBins = 33;

//for efficiency calculation
//minS2 need not match S2 threshold in detector.hh
//you can treat as trigger vs. analysis thresholds
double minS2 = 0.0;
double maxS2 = 1e9;
