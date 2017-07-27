#ifndef __NEST_H__
#define __NEST_H__ 1

#include <math.h>
#include <vector>
#include <random>
//#include <time.h>
//#include <cmath>
//#include <stdlib.h>
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <assert.h>
//#include <array>
//#include <functional>



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

    } INTERACTION_TYPES;
    
//    struct Vertex {
//    public:
//
//        Vertex(double E, double A_in, std::array<double, 3>pos_in, double t_in) : energy(E), A(A_in), time(t_in), pos(pos_in), size(0) { ; }
//        Vertex(const Vertex& orig) : energy(orig.energy), A(orig.A), time(orig.time), pos(orig.pos), size(orig.size) { ; }
//
//        static const Vertex merge(Vertex va, Vertex vb);
//
//        std::array<double, 3> getPos() const {
//            return pos;
//        }
//
//        double getA() const {
//            return A;
//        }
//        double gett() const{return time;}
//        double getEnergy() const {
//            return energy;
//        }
//    private:
//        double energy; // keV
//        double A;      // -1 for electron/positron/gamma, A for (neutral) nucleus
//        double time;   // ns
//        std::array<double, 3> pos; // mm
//        double size;   // mm
//
//    };
//
//    std::vector<Vertex> cluster(std::vector<Vertex>);
//
//    struct ExitonIonRes {
//        int    nExcitons;
//        int    nIons;
//        double fExcitons;   // like nExcitons but without fluctuations
//        double fIons;       // .. same for nIons
//    };
    struct YieldResult{
        double PhotonYield;
        double ElectronYield;
    };
//    struct PhotonElectronRes {
//        double recombProby; // recombination probability for ions
//        int    nPhotons;
//        int    nElectrons;
//        double fPhotons;    // like nPhotons but without fluctuations
//        double fElectrons;  // .. same for nElectrons
//    };
//
//    struct FullRes {
//        ExitonIonRes EIRes;
//        PhotonElectronRes PERes;
//        FullRes(ExitonIonRes e, PhotonElectronRes p): EIRes(e), PERes(p){}
//    };
//    enum ParticleType {
//        NR,
//        ER
//    };
//
//    struct Efield_dependent_params {
//        double TIB_lo, TIB100;
//        double m1, m2, m3, m4, m5;
//        double pol0_lowE, pol1_lowE, pol2_lowE, pol0_medE, pol1_medE, pol2_medE;
//        double t_a,t_b,t_c,t_n;
//    };
//
//    class Detector {
//    public:
//        Detector();
//        void XeSettings();
//
//        // to be honest, code knows only about xenon
//        int    Z;         // Z
//        double molarMass; // g/mol
//        double density;   // g/cm3
//    };

    class NESTcalc {
    private:
//        const Detector det;
//        unsigned long int randomSeed;
//        TRandom3 rng;
//        const double W; //derived from detector parameters
//        const double alf_special; //derived from detector parameters
//        const double detfactor,zfactor; //derived from detector parameters

    protected:
        std::ranlux24 rng;
        double rand_uniform();
        double rand_gauss( double mean, double sigma );
//        double LindhardFactor(int A, double lambda);
//        int modPoisRnd(double, double);
//        int modBinom(int, double, double);
        int BinomFluct(int, double);

    public:
        NESTcalc();
        
        YieldResult GetYields( int species, double energy, double density, double dfield );

//        void SetRandomSeed(unsigned long int);
//        unsigned long int GetRandomSeed();
//
//        Detector GetDetectorParameters();
//
//        double TIB(const Vertex &vertex, double electricField);
//        double RecombProby(int numIons, const Vertex &vertex, double electricField);
//        static std::function<Efield_dependent_params(double)> mem_calc_Efield_params;
//        static std::function<double (int, double) > LindhardKFactor;
//
//        ExitonIonRes ExcitonsIons(const Vertex &vertex, double electricField);
//        PhotonElectronRes PhotonsElectrons(const Vertex &vertex, double electricField);
//        PhotonElectronRes PhotonsElectrons(const Vertex &vertex, double electricField, const ExitonIonRes &ei_result);
//        FullRes FullCalculation(const Vertex &vertex, double electricField);

    };
}

//std::ostream& operator<<(std::ostream &os, const NEST::ExitonIonRes &res);
//std::ostream& operator<<(std::ostream &os, const NEST::PhotonElectronRes &res);
//std::ostream& operator<<(std::ostream &os, const NEST::FullRes &res);

#endif
