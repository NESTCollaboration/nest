using namespace std;
using namespace NEST;

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol);

void GetEnergyRes(vector<double> Es);

int testNEST(VDetector* detector, unsigned long int numEvts, string type,
             double eMin, double eMax, double inField, string position, string posiMuon,
             double fPos, int seed, bool no_seed);

double band[NUMBINS_MAX][6];
double energies[3];
