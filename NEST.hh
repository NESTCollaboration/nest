#ifndef __NEST_H__
#define __NEST_H__ 1

#include <math.h>
#include <vector>
#include <random>
#include <iostream>
#include <assert.h>
#include <float.h>

#include "RandomGen.hh"
#include "Detectors/VDetector.hh"

#define W_DEFAULT 13.7
#define NEST_AVO 6.0221409e+23
#define ATOM_NUM 54.
#define MOLAR_MASS 131.293
#define PHE_MIN 1e-6

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
    NoneType=11
    
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
    int ions;
    int excitons;
  };
  
  typedef std::vector<double> photonstream;
  
  struct NESTresult {
    YieldResult yields;
    QuantaResult quanta;
    photonstream photon_times;
  };
  
  
	class NESTcalc {
    
		private:

			long double Factorial ( double x );
			double nCr ( double n, double r );
			
		public:

			NESTcalc();
			NESTcalc(VDetector* detector);
			long BinomFluct(long, double);
			NESTresult FullCalculation(INTERACTION_TYPE species,double energy,double density,double dfield,double A,double Z,std::vector<double> NuisParam);
	  double PhotonTime(INTERACTION_TYPE species,bool exciton, double dfield, double energy);
	  photonstream GetPhotonTimes(INTERACTION_TYPE species, QuantaResult result, double dfield, double energy);
			YieldResult GetYields ( INTERACTION_TYPE species, double energy, double density, double dfield,double A,double Z,std::vector<double> NuisParam);
			QuantaResult GetQuanta(YieldResult yields, double density);
			std::vector<double> GetS1 ( int Nph,double dx, double dy, double dz, double driftSpeed, double dS_mid, INTERACTION_TYPE species );
			std::vector<double> GetSpike(int Nph,double dx,double dy, double dz, double driftSpeed, double dS_mid, std::vector<double> origScint );
			std::vector<double> GetS2 ( int Ne, double dx, double dy, double dt, double driftSpeed );
			double SetDriftVelocity ( double T, double D, double F );
			double SetDriftVelocity_MagBoltz ( double D, double F );
			std::vector<double> SetDriftVelocity_NonUniform ( double rho, double zStep );
			double SetDensity ( double T, double P );
			std::vector<double> xyResolution ( double xPos_mm, double yPos_mm, double A_top );
			
		protected:
			
			VDetector* fdetector;
	  
  };

}


#endif
