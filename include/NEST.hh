#ifndef __NEST_H__
#define __NEST_H__ 1

#include <math.h>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <algorithm>

#include "RandomGen.hh"
#include "VDetector.hh"

#define W_DEFAULT 13.7 //default work function, in eV
#define NEST_AVO 6.0221409e+23
#define ATOM_NUM 54. //period to make float
#define MOLAR_MASS 131.293 //grams per mole
#define PHE_MIN 1e-6 //area
#define ELEC_MASS 9.109e-31 //kg
#define FIELD_MIN 1. //min elec field to make S2 (in V/cm)

#define SAMPLE_SIZE 10 //nano-seconds
#define PULSE_WIDTH 10 //nano-seconds
#define PULSEHEIGHT 0.005 //threshold height, in PE

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
    C14 = 10,    
    Kr83m = 11,
    NoneType=12
    
  } INTERACTION_TYPE;
  
  struct YieldResult {
    double PhotonYield;
    double ElectronYield;
    double ExcitonRatio;
    double Lindhard;
    double ElectricField;
  	double DeltaT_Scint;
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
    ~NESTcalc();
    
    long BinomFluct(long, double);
    NESTresult FullCalculation(INTERACTION_TYPE species,double energy,double density,double dfield,double A,double Z,std::vector<double> NuisParam={1,1});
    double PhotonTime(INTERACTION_TYPE species,bool exciton, double dfield, double energy);
    photonstream AddPhotonTransportTime (photonstream emitted_times, double x, double y, double z);
    photonstream GetPhotonTimes(INTERACTION_TYPE species, int total_photons, int excitons, double dfield, double energy);
    YieldResult GetYields ( INTERACTION_TYPE species, double energy, double density, double dfield,double A,double Z,std::vector<double> NuisParam);
    QuantaResult GetQuanta(YieldResult yields, double density);
    std::vector<double> GetS1 ( QuantaResult quanta, double dx, double dy, double dz, double driftSpeed, double dS_mid, INTERACTION_TYPE species, long evtNum, double dfield, double energy, bool useTiming,  bool outputTiming, vector<long int>& wf_time, vector<double>& wf_amp );
    std::vector<double> GetSpike(int Nph,double dx,double dy, double dz, double driftSpeed, double dS_mid, std::vector<double> origScint );
    std::vector<double> GetS2 ( int Ne, double dx, double dy, double dt, double driftSpeed, long evtNum, double dfield, bool useTiming, bool outputTiming, vector<long int>& wf_time, vector<double>& wf_amp,vector<double> &g2_params );
    std::vector<double> CalculateG2( bool verbosity = true );
    double SetDriftVelocity ( double T, double D, double F );
    double SetDriftVelocity_MagBoltz ( double D, double F );
    std::vector<double> SetDriftVelocity_NonUniform ( double rho, double zStep );
    double SetDensity ( double T, double P );
    std::vector<double> xyResolution ( double xPos_mm, double yPos_mm, double A_top );

  private:
    ofstream pulseFile;
    
  protected:
			
    VDetector* fdetector;
	  
  };

}

#endif
