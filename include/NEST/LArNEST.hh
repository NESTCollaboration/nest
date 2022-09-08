/**
 * @file LArNEST.hh
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Justin Mueller [Justin.Mueller@colostate.edu]
 * @author Ekaterina Kozlova [aspelene@gmail.com]
 * @author Michael Mooney [mrmooney@colostate.edu]
 * @brief 
 * @version 
 * @date 2022-04-13
 */
#pragma once

#include "NEST.hh"
#include "RandomGen.hh"

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

    enum class LArInteraction
    {
        NR = 0,
        ER = 1,
        Alpha = 2,
    };

    enum class LArFluctuationModel
    {
        Default = 0,
    };

    struct LArNRYieldsParameters
    {
        double alpha =  {11.10};
        double beta =   {0.087};
        double gamma =  {0.1};
        double delta =  {-0.0932};
        double epsilon ={2.998};
        double zeta =   {0.3};
        double eta =    {2.94};
    };
    struct LArERElectronYieldsAlphaParameters
    {
        double A = {32.988};
        double B = {-552.988};
        double C = {17.2346};
        double D = {-4.7};
        double E = {0.025115};
        double F = {0.265360653};
        double G = {0.242671};
    };
    struct LArERElectronYieldsBetaParameters
    {
        double A = {0.778482};
        double B = {25.9};
        double C = {1.105};
        double D = {0.4};
        double E = {4.55};
        double F = {-7.502};
    };
    struct LArERElectronYieldsGammaParameters
    {
        double A = {0.659509};
        double B = {1000};
        double C = {6.5};
        double D = {5.0};
        double E = {-0.5};
        double F = {1047.408};
        double G = {0.01851};
    };
    struct LArERElectronYieldsDokeBirksParameters
    {
        double A = {1052.264};
        double B = {14159350000 - 1652.264};
        double C = {-5.0};
        double D = {0.157933};
        double E = {1.83894};
    };
    struct LArERYieldsParameters
    {
        LArERElectronYieldsAlphaParameters alpha;
        LArERElectronYieldsBetaParameters beta;
        LArERElectronYieldsGammaParameters gamma;
        LArERElectronYieldsDokeBirksParameters doke_birks;
        double p1 = {1.0};
        double p2 = {10.304};
        double p3 = {13.0654};
        double p4 = {0.10535};
        double p5 = {0.7};
        double delta = {15.7489};
        double let = {-2.07763};
    };

    struct LArAlphaElectronYieldsParameters
    {
        double A = {1.0/6200.0};
        double B = {64478398.7663};
        double C = {0.173553719};
        double D = {1.21};
        double E = {0.02852};
        double F = {0.01};
        double G = {4.71598};
        double H = {7.72848};
        double I = {-0.109802};
        double J = {3.0};
    };

    struct LArAlphaPhotonYieldsParameters
    {
        double A = {1.5};
        double B = {-0.012};
        double C = {1.0/6500.0};
        double D = {278037.250283};
        double E = {0.173553719};
        double F = {1.21};
        double G = {2};
        double H = {0.653503};
        double I = {4.98483};
        double J = {10.0822};
        double K = {1.2076};
        double L = {-0.97977};
        double M = {3.0};
    };

    struct LArAlphaYieldsParameters
    {
        LArAlphaElectronYieldsParameters Ye;
        LArAlphaPhotonYieldsParameters Yph;
    };

    struct ThomasImelParameters
    {
        double A = {0.1};
        double B = {-0.0932};
    };

    struct DriftParameters
    {
        std::vector<double> A = {
            0.937729, 0.80302379, 0.7795972, 0.6911897,
            0.76551511, 0.502022794, 0.24207633
        };
        std::vector<double> B = {
            -0.0734108, -0.06694564, -0.0990952, -0.092997,
            -0.0731659, -0.06644517, -0.03558428
        };
        std::vector<double> C = {
            0.315338, 0.331798, 0.320876, 0.3295202,
            0.317972, 0.3290246, 0.33645519
        };

        std::vector<double> TempLow = {
            84., 86., 88., 92., 96., 110., 125.
        };
        std::vector<double> TempHigh = {
            86., 88., 92., 96., 110., 125., 140.
        };
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

    struct LArYieldFluctuationResult
    {
        double NphFluctuation;
        double NeFluctuation;
        double NexFluctuation;
        double NionFluctuation;
    };
    
    struct LArNESTResult
    {   
        LArYieldResult yields;
        LArYieldFluctuationResult fluctuations;
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

        //-------------------------Parameters-------------------------//
        /// set LAr parameters
        void setDensity(double density)     { fDensity = density; }
        void setRIdealGas(double RIdealGas) { fRIdealGas = RIdealGas; }
        void setRealGasA(double RealGasA)   { fRealGasA = RealGasA; }
        void setRealGasB(double RealGasB)   { fRealGasB = RealGasB; }
        void setWorkQuantaFunction(double workQuantaFunction) 
        { fWorkQuantaFunction = workQuantaFunction; }
        void setWorkIonFunction(double workIonFunction) 
        { fWorkIonFunction = workIonFunction; }
        void setWorkPhotonFunction(double workPhotonFunction) 
        { fWorkPhotonFunction = workPhotonFunction; }
        void setFanoER(double FanoER) { fFanoER = FanoER; }
        void setNexOverNion(double NexOverNion) { fNexOverNion = NexOverNion; }

        // TODO: do we need nuisance parameters?
        void setNuisanceParameters(std::vector<double> nuisanceParameters);
        void setTemperature(std::vector<double> temperature);
        
        /// setters for various parameters
        void setNRYieldsParameters(LArNRYieldsParameters NRYieldsParameters);
        void setERYieldsParameters(LArERYieldsParameters ERYieldsParameters);
        void setERElectronYieldsAlphaParameters(
            LArERElectronYieldsAlphaParameters ERElectronYieldsAlphaParameters
        );
        void setERElectronYieldsBetaParameters(
            LArERElectronYieldsBetaParameters ERElectronYieldsBetaParameters
        );
        void setERElectronYieldsGammaParameters(
            LArERElectronYieldsGammaParameters ERElectronYieldsGammaParameters
        );
        void setERElectronYieldsDokeBirksParameters(
            LArERElectronYieldsDokeBirksParameters ERElectronYieldsDokeBirksParameters
        );
        void setThomasImelParameters(ThomasImelParameters thomasImelParameters);
        void setDriftParameters(DriftParameters driftParameters);

        /// get LAr parameters
        double getDensity()             const { return fDensity; }
        double getRIdealGas()           const { return fRIdealGas; }
        double getRealGasA()            const { return fRealGasA; }
        double getRealGasB()            const { return fRealGasB; }
        double getWorkQuantaFunction()  const { return fWorkQuantaFunction; }
        double getWorkIonFunction()     const { return fWorkIonFunction; }
        double getWorkPhotonFunction()  const { return fWorkPhotonFunction; }
        double getFanoER()              const { return fFanoER; }
        double getNexOverNion()         const { return fNexOverNion; }

        LArNRYieldsParameters getNRYieldsParameters() { return fNR; }
        LArERYieldsParameters getERYieldsParameters() { return fER; }
        LArERElectronYieldsAlphaParameters getERElectronYieldsAlphaParameters()
            { return fER.alpha; }
        LArERElectronYieldsBetaParameters getERElectronYieldsBetaParameters()
            { return fER.beta; }
        LArERElectronYieldsGammaParameters getERElectronYieldsGammaParameters()
            { return fER.gamma; }
        LArERElectronYieldsDokeBirksParameters getERElectronYieldsDokeBirksParameters()
            { return fER.doke_birks; }
        ThomasImelParameters getThomasImelParameters() {return fThomasImelParameters; }
        DriftParameters getDriftParameters() {return fDriftParameters; }

        //-------------------------All Yields-------------------------//
        LArYieldResult GetRecombinationYields(
            double TotalYields, double ElectronYields, double PhotonYields,
            double energy, double efield
        );
        LArYieldResult GetYields(
            LArInteraction species, double energy, 
            double efield, double density
        );
        /**
         * @brief Calculate fluctions on the mean yields
         * 
         */
        LArYieldFluctuationResult GetYieldFluctuations(
            const LArYieldResult &yields, 
            double density
        );
        /**
         * @brief 
         * 
         */
        LArNESTResult FullCalculation(
            LArInteraction species, double energy, 
            double efield, double density, 
            bool do_times
        );

        //-------------------------NR Yields-------------------------//
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
        double GetNRElectronYields(double energy, double efield);
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
         * @param efield 
         * @param NuisParam 
         * @return YieldResult 
         */
        LArYieldResult GetNRYields(
            double energy, double efield, double density
        ); 

        //-------------------------ER Yields-------------------------//
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
        double GetERElectronYieldsAlpha(double efield, double density);
        /**
         * @brief ER Exciton Yields Beta function:
         * 
         * @param efield 
         * @return double 
         */
        double GetERElectronYieldsBeta(double efield);
        /**
         * @brief ER Exciton Yields Gamma function:
         * 
         * @param efield 
         * @return double 
         */
        double GetERElectronYieldsGamma(double efield);
        /**
         * @brief ER Exciton Yields Doke-Birks function:
         * 
         * @param efield 
         * @return double 
         */
        double GetERElectronYieldsDokeBirks(double efield);
        /**
         * @brief ER Exciton Yields function:
         * 
         * @param energy 
         * @param efield 
         * @param density 
         * @return double 
         */
        double GetERElectronYields(double energy, double efield, double density);
        /**
         * @brief Calculate yields for charged particles.  The formulas
         * for these are given by:
         *      Qy = alpha * beta + (gamma - alpha * beta)/(p1 + p2 * E^p3)^p4
         *                        + delta / (p5 + Doke * E^LET)
         *      
         *      Ly = Nq - Qy
         * 
         *      The p1-p5 parameters are stored in ERQuantaParameters.
         * @param energy 
         * @param density 
         * @param efield 
         * @return LArYieldResult 
         */
        LArYieldResult GetERYields(
            double energy, double efield, double density
        );
        //-------------------------Alpha Yields-------------------------//
        /**
         * @brief Get the Alpha Total Yields object
         * 
         * @param energy 
         * @return double 
         */
        double GetAlphaTotalYields(double energy);
        /**
         * @brief Get the Alpha Electron Yields object
         * 
         * @param efield 
         * @return double 
         */
        double GetAlphaElectronYields(double efield);
        /**
         * @brief Get the Alpha Photon Yields object
         * 
         * @param efield 
         * @return double 
         */
        double GetAlphaPhotonYields(double efield);
        /**
         * @brief Calculate yields for alpha particles
         * 
         * @param energy 
         * @param efield 
         * @param density 
         * @return LArYieldResult 
         */
        LArYieldResult GetAlphaYields(
            double energy, double efield, double density
        );
        //-------------------------Fluctuation Yields-------------------------//
        /**
         * @brief Get the Default Fluctuations object
         * 
         * @param yields 
         * @param density 
         * @return LArYieldFluctuationResult 
         */
        LArYieldFluctuationResult GetDefaultFluctuations(
            const LArYieldResult &yields, double density
        );
        //-------------------------Photon Times-------------------------//
        double GetPhotonTime(
            LArInteraction species, bool exciton,
            double energy
        );
        double GetPhotonEnergy(bool state);

        //-------------------------Drift Velocity-------------------------//
        double GetDriftVelocity_Liquid(
            double Kelvin, double eField
        );
        double GetDriftVelocity_MagBoltz(
            double density, double efieldinput,
            double molarMass
        );

        //-------------------------Utilities-------------------------//
        double GetLinearEnergyTransfer(
            double energy, bool CSDA=false
        );
        double GetDensity(
            double Kelvin, double bara, bool &inGas,
            uint64_t evtNum, double molarMass
        );
        std::vector<double> CalculateG2(int verbosity=-1);

        //-------------------------Legacy LArNEST-------------------------//
        /**
         * @brief 
         * 
         * @param energy 
         * @param efield 
         * @param yieldFactor 
         * @param excitationRatio 
         * @param epsilon 
         * @param recombProb 
         * @return LArYieldResult 
         */
        LArYieldResult LegacyGetYields(
            double energy, double efield, 
            double yieldFactor, double excitationRatio, 
            double epsilon, double recombProb
        );
        /**
         * @brief Below are legacy LAr calculation 
         * functions which are almost copied verbatim 
         * from LArSoft, see - 
         * 
         * The default value for the dx=track_length is
         * 1/3 mm, which is the standard from Geant4.
         * @return LArYieldResult 
         */
        LArYieldResult LegacyCalculation(
            int pdgcode, double energy,
            double efield, double density,
            double track_length=0.0003
        );
        double LegacyGetRecombinationProbability(
            double energy, double efield, double density,
            int pdgcode, double track_length
        );
        double LegacyGetLinearEnergyTransfer(double E);

    private:
        double fDensity =   {1.393};        
        double fRIdealGas = {8.31446261815324};
        double fRealGasA =  {0.1355};  // m^6*Pa/mol^2 or m^4*N/mol^2.
        double fRealGasB =  {3.201e-5};  // m^3/mol.

        double fWorkQuantaFunction = {19.5};
        double fWorkIonFunction =    {23.6};
        double fWorkPhotonFunction = {14.544};

        double fNexOverNion = {0.21};
        double fALF = {1. / (1. + 0.21)};
        double fFanoER = {0.1115};
    
        LArNRYieldsParameters fNR;
        LArERYieldsParameters fER;
        LArAlphaYieldsParameters fAlpha;

        ThomasImelParameters fThomasImelParameters;
        DriftParameters fDriftParameters;

        enum LArFluctuationModel fLArFluctuationModel;
    };
}
