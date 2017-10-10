#ifndef __NEST_H__
#define __NEST_H__ 1

#include <math.h>
#include <vector>
#include <random>



#define NEST_AVO 6.022e23


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
    

    struct YieldResult{
        double PhotonYield;
        double ElectronYield;
    };
    
    struct QuantaResult{
        int photons;
        int electrons;
    };
    
    typedef std::vector<double> photonstream;
    
    struct NESTresult{
        YieldResult yields;
        QuantaResult quanta;
        photonstream photon_times;
        
    };

    
    
    class NESTcalc {
    private:


    protected:
        std::ranlux24 rng;
        
        double rand_gauss( double mean, double sigma );
        int BinomFluct(int, double);

    public:
        NESTcalc();
        double rand_uniform();
      std::vector<double> VonNeumann(double xMin, double xMax, double yMin,double yMax,
				     double xTest,double yTest,double fValue);
        NESTresult FullCalculation(INTERACTION_TYPE species, double energy, double density, double dfield);
        double PhotonTime(INTERACTION_TYPE species,bool exciton);
        photonstream GetPhotonTimes(/*inputs*/);
        YieldResult GetYields( INTERACTION_TYPE species, double energy, double density, double dfield );
        QuantaResult GetQuanta(YieldResult yields);
        void SetRandomSeed(unsigned long int);

    };
}


#endif
