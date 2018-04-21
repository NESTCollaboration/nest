#ifndef __NEST_H__
#define __NEST_H__ 1

#include <math.h>
#include <vector>
#include <random>
#include <iostream>
#include <assert.h>
#include <float.h>

#define W_DEFAULT 13.7
#define NEST_AVO 6.0221409e+23
#define ATOM_NUM 54.
#define MOLAR_MASS 131.293

namespace NEST {
  
  typedef enum {
    
    //nuclear recoil
    NR = 0,
    WIMP = 1,
    B8 = 2,
    DD = 3,
    AmBe = 4,
    Cf = 5,
    ion = 6, //includes alphas, Pb-206
    //electron recoil
    gammaRay = 7,
    beta = 8,
    CH3T = 9,
    Kr83m = 10,
    
  } INTERACTION_TYPE;
  
  struct YieldResult {
    double PhotonYield;
    double ElectronYield;
    double ExcitonRatio;
    double Lindhard;
    double ElectricField;
  };
  
  struct QuantaResult{
    int photons;
    int electrons;
  };
  
  typedef std::vector<double> photonstream;
  
  struct NESTresult {
    YieldResult yields;
    QuantaResult quanta;
    photonstream photon_times;
  };
  
  struct DetectorParameters {
    double temperature;
    double pressure;
    double GXeInterface;
    double efFit;
    double rad;
    double dtExtrema[2];
  };
  
  class NESTcalc {
    
  private:

    long double Factorial ( double x );
    double nCr ( double n, double r );
    
  protected:

    std::ranlux24 rng;
    
  public:

    NESTcalc();
    double rand_uniform();
    int poisson_draw(double mean);
    double rand_exponential(double half_life);
    std::vector<double> VonNeumann(double xMin, double xMax, double yMin,double yMax,
				   double xTest,double yTest,double fValue);
    double rand_gauss( double mean, double sigma );
    int BinomFluct(int, double);
    NESTresult FullCalculation(INTERACTION_TYPE species,double energy,double density,double dfield,double A,double Z);
    double PhotonTime(INTERACTION_TYPE species,bool exciton);
    photonstream GetPhotonTimes(/*inputs*/);
    YieldResult GetYields(INTERACTION_TYPE species,double energy,double density,double dfield,double A,double Z);
    QuantaResult GetQuanta(YieldResult yields, double density);
    DetectorParameters GetDetector ( double x, double y, double z, bool IsInGasPhase );
    void DriftRangeOverride ( double drift_low, double drift_high, DetectorParameters &detParam );
    std::vector<double> GetS1 ( int Nph,double dx, double dy, double dz, double driftSpeed );
    std::vector<double> GetSpike(int Nph,double dx,double dy, double dz, double driftSpeed, std::vector<double> origScint );
    std::vector<double> GetS2 ( int Ne, double dx, double dy, double dt, bool IsInGasPhase );
    int SelectRanXeAtom (double isotope);
    void SetRandomSeed(unsigned long int);
    
  };

}

#endif
