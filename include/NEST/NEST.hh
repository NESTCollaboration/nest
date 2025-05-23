#ifndef __NEST_H__
#define __NEST_H__ 1

/*

  NEST.hh and NEST.cpp are the core, the heart, the most important parts of the
  code and model.
  Here the mean light and charge yields are produced but also the fluctuations.
  The yields and fluctuated quanta are functions of the particle type, energy,
  field, and density/phase.
  The raw numbers of photons and electrons are converting into S1 and S2 as
  observed in real life.
  S1 and S2 totals have arrival timing and basic waveform digitization sim added
  on top of them well.
  Lastly, the NEST functions handle a few more advanced detector property
  calculations.
  You'll find the drift v, density as a function of T and P, and XY position
  resolution smearing.

*/

/* *** A word on the units for variables! ***

   distances/lengths are typically mm although a few may be in cm. Velocity
   always uses mm and field always uses cm
   energies are in keV. Work functions are eV sometimes and converted internally
   inside of functions
   fields can be either V/cm or kV/cm, depending on what is most natural for a
   given problem
   times generally end up in ns (and NOT 100 MHz "samples" e.g.) If us or rarely
   ms they're each ultimately converted
   - Microseconds are used as the initial units for drift time, before final
   conversions into ns for digitized pulses
   The magnitude of the drift velocity is given in mm per us, another exception
   to the ns rule
   Temperatures must be in Kelvin not Celsius and pressures are always in bar
   (absolute)
   S1 and S2 units vary from photons and electrons to photo-electrons (labeled
   PE or phe since there are 2 conventions)
   - phd or dph or "photons detected" or "detected photons" refers to PE after
   2pe effect is approximately averaged out (reduction by ~20%)
   - Spikes refer to an attempt at photon counting, hits in PMTs, improving S1
   resolution compared to traditional smeared area
   Even before XYZ position corrections, a spike does not have to be a whole
   number due to things like probabilistic overlap bias correction
   Expect mass densities to be in grams per cm^3 or mL. If atomic density is
   needed, then we derive it
   Particles like photons, electrons, excitons, ions are unit-less, as are the
   "Nuisance Parameters" (multiplicative factors on yields)
   - Event number also dimensionless, and follows C++ convention of starting at
   0, not at 1, for the first event in your simulation
   Probabilities are all initially fractional even if converted to percentages
   for display purposes
   A=atomic mass, Z=mass number (#protons) for simulating the yields of random
   non-noble ions/atoms traveling through medium.
   - These are round numbers even if input as double's since go into
   calculations with other doubles as a mixture
   The interaction type ("species") is a new type defined for use in NEST that
   attempts to cover all possible kinds of particle/scatter
   For waveforms (with output timing set to >=1) the unit is PE per bin
   (sample, usually 10 ns) but has option of total area in PE

*/

#include "RandomGen.hh"
#include "VDetector.hh"
#include "ValidityTests.hh"
#include "gcem.hpp"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <array>
#include <vector>
#include <string>

#define HIGH_E_NR 330.
#define W_DEFAULT \
  13.4  // default work func, in eV. arXiv:1611.10322. +/- 0.35. 19.5-19.6 eV
        // for LAr
#define W_SCINT \
  8.5e-3  // the *max* possible energy of 1 scint phot, keV. Make this at least
          // 10 eV for LAr
#define NEST_AVO 6.0221409e+23
#define ATOM_NUM \
  54.  // period to make float. 18 for LAr. If changed here go to TestSpectra.hh
       // too

#define PHE_MIN 1e-6  // area
#define FIELD_MIN 1.  // min elec field to make S2 (in V/cm)
#define DENSITY 2.90  // g/cm^3, ref density for dependent effects. ~1.4 for LAr

#define EPS_GAS \
  1.00126  // poly-morphic: make negative to use Aprile/PandaX instead of
           // LLNL/PIXeY's e- EE
// for GAr it is 1.000574 at least at room T (doi.org/10.1103/PhysRev.34.615)
#define EPS_LIQ \
  1.85  // LXe dielectric constant explicitly NOT 1.96 (old). Update thx to Dan
        // M. LAr 1.5

#define SAMPLE_SIZE 10  // nano-seconds, 5 for LAr
#define PULSE_WIDTH 10  // nano-seconds, 5 for LAr
#define PULSEHEIGHT \
  0.005                  // threshold height, in PE, for writing to photon_times
#define SPIKES_MAXM 120  // above this switch to pulse area (70 phd in 1 array)
#define PHE_MAX 140      // saturation threshold, in PE/sample (LZ). LUX val 180

static constexpr int XYcorr =
    3;  // 0 means no corrections, 1 is for S1, 2 for S2, 3 for both

static constexpr double RidealGas = 8.31446261815324;  // Joules/mole/Kelvin
static constexpr double RealGasA =
    0.4250;  // m^6*Pa/mol^2 or m^4*N/mol^2. For Ar: 0.1355
static constexpr double RealGasB = 5.105e-5;  // m^3/mol. For Ar: 3.201e-5

const std::vector<double> default_NRYieldsParam = {
    11., 1.1, 0.0480, -0.0533, 12.6, 0.3, 2., 0.3, 2., 0.5, 1., 1.};//{11.1, 1.087, 0.1, -0.0932, 2.998, 0.3, 2.94, 0.3, 2., 0.5, 1., 1.}; for LAr
const std::vector<double> default_NRERWidthsParam = {
    0.4, 0.4, 0.04, 0.5, 0.19, 2.25, 1., 0.046452, 0.205, 0.45, -0.2};
// 1 => -0.0015 for old-style ER Fq function instead of a constant value
const std::vector<double> default_ERYieldsParam = {-1., -1., -1., -1., -1.,
                                                   -1., -1., -1., -1., -1.};
// Fano factor of ~3 at least for ionization in NRERWidthsParam if using
// OldW13eV (look at first 2 values). Also, 0.046452 used to be 0.05(53).
const std::vector<double> default_EnergyParams = {0.23, 0.77, 2.95, -1.44};
const std::vector<double> default_FieldParams = {421.15, 3.27};

namespace NEST {

typedef enum {

  // nuclear recoil
  NR = 0,
  WIMP = 1,
  B8 = 2,
  DD = 3,
  AmBe = 4,
  Cf = 5,
  ion = 6,  // includes alphas, Pb-206
  // electron recoil
  gammaRay = 7,  // only photoelectric
  beta = 8,      // includes comptons
  CH3T = 9,
  C14 = 10,
  Kr83m = 11,
  ppSolar = 12,
  atmNu = 13,
  fullGamma = 14,
  fullGamma_PE = 15,
  fullGamma_Compton_PP = 16,
  NoneType = 17

} INTERACTION_TYPE;

enum class S1CalculationMode { Full, Parametric, Hybrid, Waveform };

enum class S2CalculationMode { Full, Waveform, WaveformWithEtrain };

struct YieldResult {
  double PhotonYield;
  double ElectronYield;
  double ExcitonRatio;
  double Lindhard;
  double ElectricField;
  double DeltaT_Scint;
};

struct QuantaResult {
  int photons;
  int electrons;
  int ions;
  int excitons;
  double recombProb;
  double Variance;
};

typedef std::vector<double> photonstream;

struct NESTresult {
  YieldResult yields;
  QuantaResult quanta;
  photonstream photon_times;
};

class NESTcalc {
 protected:
  VDetector *fdetector;

 private:
  ofstream pulseFile;

  static long double Factorial(double x);

  static double nCr(double n, double r);

  vector<vector<double>> photon_areas;
  vector<double> scintillation;  // return vector
  vector<double> newSpike;       // for re-doing spike counting more precisely
  vector<double> ionization;

  static constexpr double two_PI = 2. * M_PI;
  static constexpr double sqrt2 = gcem::sqrt(2.);
  static constexpr double sqrt2_PI = gcem::sqrt(2. * M_PI);
  static constexpr double inv_sqrt2_PI = 1. / gcem::sqrt(2. * M_PI);

 public:
  NESTcalc(const NESTcalc &) = delete;

  NESTcalc &operator=(const NESTcalc &) = delete;

  explicit NESTcalc(VDetector *detector);

  virtual ~NESTcalc();

  NESTresult FullCalculation(
      INTERACTION_TYPE species, double energy, double density, double dfield,
      double A, double Z,
      const std::vector<double> &NRYieldsParam = default_NRYieldsParam,
      const std::vector<double> &NRERWidthsParam = default_NRERWidthsParam,
      const std::vector<double> &ERYieldsParam = default_ERYieldsParam,
      bool do_times =
          true);  // the so-called full NEST calculation puts together all the
                  // individual functions/calculations below

  double PhotonTime(
      INTERACTION_TYPE species, bool exciton, double dfield,
      double energy);  // gives you the birth times of S1 as well as S2
                       // scintillation photons, taking singlet, triplet, and
                       // recombination times into account, depending on
                       // particle, energy, field

  photonstream AddPhotonTransportTime(
      const photonstream &emitted_times, double x, double y,
      double
          z);  // adds an approximately ray-traced (analytical function) photon
               // travel time in the detector to each photon birth time

  photonstream GetPhotonTimes(
      INTERACTION_TYPE species, int total_photons, int excitons, double dfield,
      double energy);  // this function loops over the photon times above to add
                       // times to multiple photons in a loop

  YieldResult GetYields(
      INTERACTION_TYPE species, double energy, double density, double dfield,
      double A, double Z,
      const std::vector<double> &NRYieldsParam = default_NRYieldsParam,
      const std::vector<double> &ERYieldsParam = default_ERYieldsParam);
  // the innermost heart of NEST, this provides floating-point average values
  // for photons and electrons per keV. Nuis(ance)Param included for varying the
  // NR Ly & Qy up and down

  virtual YieldResult GetYieldGamma(double energy, double density,
                                    double dfield, double multFact = 1.);
  // Called by GetYields in the Gamma/x-ray/Photoabsorption Case

  virtual YieldResult GetYieldERWeighted(
      double energy, double density, double dfield,
      const std::vector<double> &ERYieldsParam = default_ERYieldsParam,
      const std::vector<double> &EnergyParams = default_EnergyParams,
      const std::vector<double> &FieldParams = default_FieldParams);
  // Weights beta/gamma models to account for ER sources with differing
  // recombination profiles (such as L-shell electron-capture interactions)

  virtual NESTresult GetYieldERdEOdxBasis(
      const std::vector<double> &dEOdxParam, string muonInitPos,
      vector<double> eDriftVelTable,
      const std::vector<double> &NRERWidthsParam);
  // Use dE/dx-based yield models instead of energy-based, as everywhere else

  virtual YieldResult GetYieldNR(
      double energy, double density, double dfield, double massNum,
      const std::vector<double> &NRYieldsParam = default_NRYieldsParam);
  // Called by GetYields in the NR (and related) cases

  virtual YieldResult GetYieldNROld(double energy, int alt);
  // Quick and dirty simple analytical approximations saved for earlier NEST
  // versions that were first principles: power laws, ln, sigmoid, exponentials

  virtual YieldResult GetYieldIon(
      double energy, double density, double dfield, double massNum,
      double atomNum,
      const std::vector<double> &NRYieldsParam = default_NRYieldsParam);
  // Called by GetYields in the ion case

  virtual YieldResult GetYieldKr83m(double energy, double density,
                                    double dfield,
                                    double maxTimeSeparation = 1000.,
                                    double minTimeSeparation = 300.);
  // Called by GetYields in the Kr83m case
  
  virtual YieldResult GetYieldsAndQuanta ( double energy, double density, double dfield,
  INTERACTION_TYPE species, const vector<double> &betaMeansPara, const vector<double> &nuclMeansPara );
  
  virtual YieldResult GetYieldBetaGR(
      double energy, double density, double dfield,
      const std::vector<double> &ERYieldsParam = default_ERYieldsParam,
				     double multFact = 1.);
  // Called by GetYields in the Beta/Compton/etc.(IC,Auger,EC) Case

  virtual YieldResult YieldResultValidity(YieldResult &res, const double energy,
                                          const double Wq_eV);
  // Confirms and sometimes adjusts YieldResult to make physical sense

  virtual QuantaResult GetQuanta(
      const YieldResult &yields, double density,
      const std::vector<double> &NRERWidthsParam = default_NRERWidthsParam,
      double SkewnessER = -999.);
  // GetQuanta takes the yields from above and fluctuates them, both the total
  // quanta (photons+electrons) with a Fano-like factor, and the "slosh" between
  // photons and electrons
  // Namely, the recombination fluctuations

  virtual double RecombOmegaNR(
      double elecFrac,
      const std::vector<double> &NRERWidthsParam = default_NRERWidthsParam);
  // Calculates the Omega parameter governing non-binomial recombination
  // fluctuations for nuclear recoils and ions (Lindhard<1)

  virtual double RecombOmegaER(
      double efield, double elecFrac,
      const std::vector<double> &NRERWidthsParam = default_NRERWidthsParam);
  // Calculates the Omega parameter governing non-binomial recombination
  // fluctuations for gammas and betas (Lindhard==1)

  virtual double FanoER(
      double density, double Nq_mean, double efield,
      const std::vector<double> &NRERWidthsParam = default_NRERWidthsParam);
  // Fano-factor (and Fano-like additional energy resolution model) for gammas
  // and betas (Lindhard==1)

  const std::vector<double> &GetS1(
      const QuantaResult &quanta, double truthPosX, double truthPosY,
      double truthPosZ, double smearPosX, double smearPosY, double smearPosZ,
      double driftSpeed, double dS_mid, INTERACTION_TYPE species,
      uint64_t evtNum, double dfield, double energy, S1CalculationMode mode,
      int outputTiming, vector<int64_t> &wf_time, vector<double> &wf_amp);
  // Very comprehensive conversion of the "original" intrinsic scintillation
  // photons into the many possible definitions of S1 as measured by
  // photo-sensors

  const std::vector<double> &GetSpike(int Nph, double dx, double dy, double dz,
                                      double driftSpeed, double dS_mid,
                                      const std::vector<double> &origScint);
  // GetSpike takes the extremely basic digital/integer number of spike counts
  // provided by GetS1 and does more realistic smearing

  const std::vector<double> &GetS2(
      int Ne, double truthPosX, double truthPosY, double truthPosZ,
      double smearPosX, double smearPosY, double smearPosZ, double dt,
      double driftSpeed, uint64_t evtNum, double dfield, S2CalculationMode mode,
      int outputTiming, vector<int64_t> &wf_time, vector<double> &wf_amp,
      const vector<double> &g2_params);
  // Exhaustive conversion of the intrinsic ionization electrons into the many
  // possible definitions of S2 pulse areas as observed in the photo-tubes
  // This function also applies the extraction efficiency (binomial) and finite
  // electron mean free path or life time caused by electronegative impurities
  // (exponential)

  std::vector<double> CalculateG2(int verbosity = 1);
  // Calculates "g2" by combining the single electron size with the extraction
  // efficiency. Called by GetS2 above. Includes helper variables like gas gap
  // and SE width.

  double SetDriftVelocity(double T, double D, double F, double P);
  // Gives one the drift velocity as a function of temperature and electric
  // field in liquid or solid. If density implies gas, kicks calculation down to
  // the next function below

  static double GetDriftVelocity(double T, double D, double F, bool inGas, double P);
  // Gives one the drift velocity as a function of temperature and electric
  // field in liquid or solid. If density implies gas, kicks calculation down to
  // the next function below

  static double GetDriftVelocity_Liquid(double T, double F, double D, double P = 1.5, short SD = -1);
  // Gives one the drift velocity as a function of temperature and electric
  // field in liquid or solid. If density implies gas, kicks calculation down to
  // the next function below. NOTE: Density default implies liquid

  static double GetDriftVelocity_MagBoltz(double T, double D, double F, double P,
                                          double molarMass = 131.293);
  // Gas electron drift speed for S2 gas gap in 2-phase TPCs or the whole
  // detector for all gas. Based on simple fits to complicated MagBoltz software
  // output.

  std::vector<double> SetDriftVelocity_NonUniform(double rho, double zStep, double T, double P,
                                                  double dx, double dy);
  // Special handling for the case of non-uniform electric fields in a detector,
  // this integrates over position to find the correct total drift time from any
  // starting point

  double SetDensity(double T, double P);
  // A simple, approximate but good, density is returned for solid, liquid, or
  // gaseous xenon, as a function of temperature and pressure

  static double GetDensity(double T, double P, bool &inGas, uint64_t evtNum = 0,
                           double molarMass = 131.293);
  // A simple, approximate but good, density is returned for solid, liquid, or
  // gaseous xenon, as a function of temperature and pressure

  std::vector<double> xyResolution(double xPos_mm, double yPos_mm,
                                   double A_top);
  // Utilizing a dependence on radius and the size of the S2 signal, takes MC
  // truth X and Y and outputs smeared values as if you did position
  // reconstruction like in real data

  virtual double PhotonEnergy(bool s2Flag, bool state, double tempK);
  // Determines the birth energies in electron-Volts of scintillation photons,
  // for either S1 or S2, including fluctuations in them, so that you can apply
  // proper QE in G4 for ex.

  double CalcElectronLET(double E, int Z, bool CSDA = true);
  // Linear Energy Transfer in units of MeV*cm^2/gram which when combined with
  // density can provide the dE/dx, as a function of energy in keV. Will be more
  // useful in the future

  struct Wvalue {
    double Wq_eV;
    double alpha;
  };

  static Wvalue WorkFunction(double rho, double MolarMass,
                             bool OldW13eV = true);
  // the W-value as a func of density in g/cm^3

  virtual double NexONi(double energy, double density);
  // calculate exciton/ion

  VDetector *GetDetector() { return fdetector; }

  void SetDetector(VDetector *detector) { fdetector = detector; }

  // Access the diffusion coefficient for transverse diffusion in liquid
  static double GetDiffTran_Liquid(double dfield, bool highFieldModel = false,
                                   double T = 175., double P = 2., double rho = DENSITY, int Z = 54);

  // Access the diffusion coefficient for longitudinal diffusion in liquid
  static double GetDiffLong_Liquid(double dfield, bool highFieldModel = false,
                                   double T = 175., double P = 2., double rho = DENSITY, int Z = 54, short SD = -1);

  // Function helpful for interpolation of the new diffusion coefficient model
  // (Boyle)
  static double interpolateFunction(
      const std::vector<std::pair<double, double>> &func, double x,
      bool isLogLog);

  // Read in the Boyle model data for DT
  static std::vector<std::pair<double, double>> GetBoyleModelDT();

  // Read in the Boyle model data for DL
  static std::vector<std::pair<double, double>> GetBoyleModelDL();

  static constexpr int clamp(int v, const int lo, const int hi);
};
}  // namespace NEST

#endif
