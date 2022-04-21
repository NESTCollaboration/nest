/**
 * @file LArNEST.hh
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Justin Mueller [Justin.Mueller@colostate.edu]
 * @author Michael Mooney [mrmooney@colostate.edu]
 * @brief 
 * @version 
 * @date 2022-04-13
 */
#pragma once

#include "NEST.hh"
#include "RandomGen.hh"

// #define HIGH_E_NR 330.
// #define W_DEFAULT \
//   13.4  // default work func, in eV. arXiv:1611.10322. +/- 0.35. 19.5-19.6 eV
//         // for LAr
// #define W_SCINT \
//   8.5e-3  // the *max* possible energy of 1 scint phot, keV. Make this at least
//           // 10 eV for LAr
// #define NEST_AVO 6.0221409e+23
// #define ATOM_NUM \
//   18.  // period to make float. 18 for LAr. If changed here go to TestSpectra.hh
//        // too

// #define PHE_MIN 1e-6  // area
// #define FIELD_MIN 1.  // min elec field to make S2 (in V/cm)
// #define DENSITY 2.90  // g/cm^3, ref density for dependent effects. ~1.4 for LAr

// #define EPS_GAS \
//   1.00126  // poly-morphic: make negative to use Aprile/PandaX instead of
//            // LLNL/PIXeY's e- EE
// // for GAr it is 1.000574 at least at room T (doi.org/10.1103/PhysRev.34.615)
// #define EPS_LIQ \
//   1.85  // LXe dielectric constant explicitly NOT 1.96 (old). Update thx to Dan
//         // M. LAr 1.325

// #define SAMPLE_SIZE 10  // nano-seconds
// #define PULSE_WIDTH 10  // nano-seconds
// #define PULSEHEIGHT \
//   0.005                  // threshold height, in PE, for writing to photon_times
// #define SPIKES_MAXM 120  // above this switch to pulse area (70 phd in 1 array)
// #define PHE_MAX 180      // saturation threshold, in PE per bin i.e. sample

#define ZurichEXOW 1.1716263232
// Qy Boost Factor to fit EXO-200 Data (error on this quantity is +/- 0.03)
#define ZurichEXOQ 1.08  

namespace NEST
{
    static constexpr double LAr_Z{18};
    static constexpr double legacy_density_LAr{1.393};
    static constexpr double legacy_scint_yield{1.0 / (19.5 * 1.e-6)};
    static constexpr double legacy_resolution_scale{0.107}; // Doke 1976
    static constexpr double two_PI = 2. * M_PI;
    static constexpr double sqrt2 = gcem::sqrt(2.);
    static constexpr double sqrt2_PI = gcem::sqrt(2. * M_PI);
    static constexpr double inv_sqrt2_PI = 1. / gcem::sqrt(2. * M_PI);

    struct LArNRYieldsParameters
    {
        double alpha =  {11.10};
        double beta =   {1.087};
        double gamma =  {0.1};
        double delta =  {-0.0932};
        double epsilon ={2.998};
        double zeta =   {0.3};
        double eta =    {2.94};
    };
    struct LArERExcitonYieldsAlphaParameters
    {
        double A = {32.988};
        double B = {-552.988};
        double C = {17.2346};
        double D = {-4.7};
        double E = {0.025115};
        double F = {3.768456};
        double G = {-0.242671};
    };
    struct LArERExcitonYieldsBetaParameters
    {
        double A = {0.778482};
        double B = {25.9};
        double C = {1.105};
        double D = {0.4};
        double E = {4.55};
        double F = {-7.502};
    };
    struct LArERExcitonYieldsGammaParameters
    {
        double A = {0.659509};
        double B = {1000};
        double C = {6.5};
        double D = {5.0};
        double E = {-0.5};
        double F = {1047.408};
        double G = {0.01851};
    };
    struct LArERExcitonYieldsDokeBirksParameters
    {
        double A = {1052.264};
        double B = {14159350000 - 1652.264};
        double C = {-5.0};
        double D = {0.157933};
        double E = {1.83894};
    };
    struct LArERYieldsParameters
    {
        LArERExcitonYieldsAlphaParameters alpha;
        LArERExcitonYieldsBetaParameters beta;
        LArERExcitonYieldsGammaParameters gamma;
        LArERExcitonYieldsDokeBirksParameters doke_birks;
        double p1 = {1.0};
        double p2 = {10.304};
        double p3 = {13.0654};
        double p4 = {0.10535};
        double p5 = {0.7};
        double delta = {15.7489};
        double let = {-2.07763};
    };

    struct LArYieldResult 
    {
        double TotalYield;
        double QuantaYield;
        double LightYield;
        double Nph;
        double Ne;
        double Nex;
        double Nion;
        double Lindhard;
        double ElectricField;
    };

    struct LArNESTResult
    {   
        LArYieldResult yields;
        QuantaResult quanta;
        photonstream photon_times;
    };

    /**
     * @brief 
     * 
     */
    class LArNEST : public NESTcalc
    {
    public:
        explicit LArNEST(VDetector *detector);

        void setDensity(double density)     { fDensity = density; }
        void setWDefault(double wDefault)   { fWDefault = wDefault; }
        void setRIdealGas(double RIdealGas) { fRIdealGas = RIdealGas; }
        void setRealGasA(double RealGasA)   { fRealGasA = RealGasA; }
        void setRealGasB(double RealGasB)   { fRealGasB = RealGasB; }
        void setNuisanceParameters(std::vector<double> nuisanceParameters);
        void setTemperature(std::vector<double> temperature);
        
        /// setters for various parameters
        void setNRYieldsParameters(LArNRYieldsParameters NRYieldsParameters);
        void setERYieldsParameters(LArERYieldsParameters ERYieldsParameters);
        void setERExcitonYieldsAlphaParameters(
            LArERExcitonYieldsAlphaParameters ERExcitonYieldsAlphaParameters
        );
        void setERExcitonYieldsBetaParameters(
            LArERExcitonYieldsBetaParameters ERExcitonYieldsBetaParameters
        );
        void setERExcitonYieldsGammaParameters(
            LArERExcitonYieldsGammaParameters ERExcitonYieldsGammaParameters
        );
        void setERExcitonYieldsDokeBirksParameters(
            LArERExcitonYieldsDokeBirksParameters ERExcitonYieldsDokeBirksParameters
        );

        /**
         * @brief NR Total Yields function:
         *      Ty = alpha * E^beta
         * 
         * @param energy 
         * @return double 
         */
        double GetNRTotalYields(double energy);
        /**
         * @brief NR Exciton Yields function:
         *      Qy = (gamma * E^delta)*(sqrt(E + eps)^-1)
         *              * (1 - (1 + (E/zeta)^eta))^-1)
         * 
         * @param energy 
         * @param efield 
         * @return double 
         */
        double GetNRExcitonYields(double energy, double efield);
        /**
         * @brief NR Photon Yields function:
         *      Ly = 
         * 
         * @param energy 
         * @param efield 
         * @return double 
         */
        double GetNRPhotonYields(double energy, double efield);
        /**
         * @brief NR Photon Yields function (conserved)
         *      Ly = Ty/E - Qy
         * 
         * @param energy 
         * @param efield 
         * @return double 
         */
        double GetNRPhotonYieldsConserved(double energy, double efield);

        /**
         * @brief ER Total Yields function:
         *      Ty
         * 
         * @param energy 
         * @return double 
         */
        double GetERTotalYields(double energy);
        /**
         * @brief ER Exciton Yields Alpha function:
         * 
         * @param efield 
         * @param density 
         * @return double 
         */
        double GetERExcitonYieldsAlpha(double efield, double density);
        /**
         * @brief ER Exciton Yields Beta function:
         * 
         * @param efield 
         * @return double 
         */
        double GetERExcitonYieldsBeta(double efield);
        /**
         * @brief ER Exciton Yields Gamma function:
         * 
         * @param efield 
         * @return double 
         */
        double GetERExcitonYieldsGamma(double efield);
        /**
         * @brief ER Exciton Yields Doke-Birks function:
         * 
         * @param efield 
         * @return double 
         */
        double GetERExcitonYieldsDokeBirks(double efield);
        /**
         * @brief ER Exciton Yields function:
         * 
         * @param energy 
         * @param efield 
         * @param density 
         * @return double 
         */
        double GetERExcitonYields(double energy, double efield, double density);
        
        /// various defined constants
        inline double GetWorkQuantaFunction();
        inline double GetWorkIonFunction();
        inline double GetWorkPhotonFunction();
        inline double GetFanoER();
        inline double GetNexOverNion(); 

        double GetLinearEnergyTransfer(
            double energy, bool CSDA=false
        );
        double GetDensity(
            double Kelvin, double bara, bool &inGas,
            uint64_t evtNum, double molarMass
        );
        double GetDriftVelocity_Liquid(
            double Kelvin, double eField
        );
        double GetDriftVelocity_MagBoltz(
            double density, double efieldinput,
            double molarMass
        );
        double GetPhotonEnergy(bool state);
        double GetPhotonTime(
            INTERACTION_TYPE species, bool exciton,
            double dfield, double energy
        );
        
        /**
         * @brief 
         * 
         */
        QuantaResult GetQuanta(
            const YieldResult &yields, double density
        );

        /**
         * @brief Calculate yields for nuclear recoils.  The formulas
         * for these are given by:
         *      Ty = alpha * E^beta
         *      Qy = (gamma * F^delta)^-1 * (sqrt(E + epsilon)) 
         *              * (1 - (1 + (E/zeta)^eta)^-1)
         * 
         *      Ly = alpha * E^(beta - 1) - (gamma * F^delta)^-1 
         *              * (sqrt(E + epsilon))^-1
         * First, we fit the total yield model and take it as fixed.
         * The model is a simple power law in the deposited energy.  There
         * is no field dependence.  (J. Mueller and E. Kozlova)
         * @param energy 
         * @param density 
         * @param dfield 
         * @param NuisParam 
         * @return YieldResult 
         */
        LArYieldResult GetNRYields(
            double energy, double dfield, double density
        ); 

        /**
         * @brief Calculate yields for beta particles.  The formulas
         * for these are given by:
         *      Qy = alpha * beta + (gamma - alpha * beta)/(p1 + p2 * E^p3)^p4
         *                        + delta / (p5 + Doke * E^LET)
         *      
         *      Ly = Nq - Qy
         * 
         *      The p1-p5 parameters are stored in ERQuantaParameters.
         * @param energy 
         * @param density 
         * @param dfield 
         * @return LArYieldResult 
         */
        LArYieldResult GetERYields(
            double energy, double density, double dfield
        );
        LArYieldResult GetAlphaYields(
            double energy, double density, double dfield
        );
        LArYieldResult GetYields(
            INTERACTION_TYPE species, double energy, 
            double density, double dfield,
            bool oldModelER=false
        );
        /**
         * @brief 
         * 
         */
        LArNESTResult FullCalculation(
            INTERACTION_TYPE species, double energy, 
            double density, double efield,
            bool do_times
        );


        std::vector<double> CalculateG2(bool verbosity=false);

        /**
         * @brief Below are legacy LAr calculation 
         * functions which are almost copied verbadim 
         * from LArSoft, see - 
         * 
         * The default value for the dx=track_length is
         * 1/3 mm, which is the standard from Geant4.
         * @return QuantaResult 
         */
        QuantaResult LegacyCalculation(
            int pdgcode, double energy,
            double density, double eField, 
            double track_length=0.0003
        );
        double LegacyGetLinearEnergyTransfer(double E);

    private:
        double fDensity = {1.4};
        double fWDefault = {19.5};
        std::vector<double> fNuisanceParameters = {
            11.0, 1.1, 0.0480, 
            -0.0533, 12.6, 0.3, 
            2., 0.3, 2., 0.5, 
            1., 1.
        };
        double fRIdealGas = {8.31446261815324};
        double fRealGasA = {0.1355};  // m^6*Pa/mol^2 or m^4*N/mol^2.
        double fRealGasB = {3.201e-5};  // m^3/mol.

        std::vector<double> fTemperature = {
           84., 86., 88., 92., 96., 110., 125., 140.
        };
        
        LArNRYieldsParameters fLArNRYieldsParameters;
        LArERYieldsParameters fLArERYieldsParameters;
    };
}