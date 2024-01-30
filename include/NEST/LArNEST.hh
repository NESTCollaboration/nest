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
#include "LArParameters.hh"

namespace NEST {
struct LArYieldResult {
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

struct LArYieldFluctuationResult {
  double NphFluctuation;
  double NeFluctuation;
  double NexFluctuation;
  double NionFluctuation;
};

struct LArNESTResult {
  LArYieldResult yields;
  LArYieldFluctuationResult fluctuations;
  photonstream photon_times;
};

/**
 * @brief LArNEST class for simulating recombination physics
 * in LAr.  There are two modes of operation, (a) data-driven and (b)
 * first-principles driven.  The data-driven method takes in an
 * energy deposition and computes
 *    <Nq>, <Ne> and <Ngamma> + fluctuations,
 * whereas the first-principles method computes
 *    <Nq>, <Nion> and <Nex> + recombination.
 *
 * Currently the data-driven method computes according to
 * NR, ER or Alpha types, but not yet configured to work with dE/dx.
 *
 * The first-principles method is only invoked for dE/dx and
 * LET models.  It uses either Doke-Birks or Thomas-Imel to
 * compute recombination.
 */
class LArNEST : public NESTcalc {
 public:
  explicit LArNEST(VDetector *detector);

  //-------------------------Parameters-------------------------//
  /// Set LAr parameters
  void SetDensity(double density) { fDensity = density; }
  void SetRIdealGas(double RIdealGas) { fRIdealGas = RIdealGas; }
  void SetRealGasA(double RealGasA) { fRealGasA = RealGasA; }
  void SetRealGasB(double RealGasB) { fRealGasB = RealGasB; }
  void SetWorkQuantaFunction(double workQuantaFunction) {
    fWorkQuantaFunction = workQuantaFunction;
  }
  void SetWorkIonFunction(double workIonFunction) {
    fWorkIonFunction = workIonFunction;
  }
  void SetWorkPhotonFunction(double workPhotonFunction) {
    fWorkPhotonFunction = workPhotonFunction;
  }
  void SetFanoER(double FanoER) { fFanoER = FanoER; }
  void SetNexOverNion(double NexOverNion) { fNexOverNion = NexOverNion; }
  void SetTemperature(double temperature) { fTemperature = temperature; }

  /// Setters for various parameters
  void SetNRYieldsParameters(LArNRYieldsParameters NRYieldsParameters) {
    fNR = NRYieldsParameters;
  }
  void SetERYieldsParameters(LArERYieldsParameters ERYieldsParameters) {
    fER = ERYieldsParameters;
  }
  void SetERElectronYieldsAlphaParameters(
      LArERElectronYieldsAlphaParameters ERElectronYieldsAlphaParameters) {
    fER.alpha = ERElectronYieldsAlphaParameters;
  }
  void SetERElectronYieldsBetaParameters(
      LArERElectronYieldsBetaParameters ERElectronYieldsBetaParameters) {
    fER.beta = ERElectronYieldsBetaParameters;
  }
  void SetERElectronYieldsGammaParameters(
      LArERElectronYieldsGammaParameters ERElectronYieldsGammaParameters) {
    fER.gamma = ERElectronYieldsGammaParameters;
  }
  void SetERElectronYieldsDokeBirksParameters(
      LArERElectronYieldsDokeBirksParameters
          ERElectronYieldsDokeBirksParameters) {
    fER.doke_birks = ERElectronYieldsDokeBirksParameters;
  }
  void SetThomasImelParameters(ThomasImelParameters thomasImelParameters) {
    fThomasImelParameters = thomasImelParameters;
  }
  void SetDriftParameters(DriftParameters driftParameters) {
    fDriftParameters = driftParameters;
  }

  /// Get LAr parameters
  double GetDensity() const { return fDensity; }
  double GetRIdealGas() const { return fRIdealGas; }
  double GetRealGasA() const { return fRealGasA; }
  double GetRealGasB() const { return fRealGasB; }
  double GetWorkQuantaFunction() const { return fWorkQuantaFunction; }
  double GetWorkIonFunction() const { return fWorkIonFunction; }
  double GetWorkPhotonFunction() const { return fWorkPhotonFunction; }
  double GetEffectiveWorkIonFunction() const {
    return fWorkIonFunction + fNexOverNion * fWorkPhotonFunction;
  }
  double GetFanoER() const { return fFanoER; }
  double GetNexOverNion() const { return fNexOverNion; }

  LArNRYieldsParameters GetNRYieldsParameters() { return fNR; }
  LArERYieldsParameters GetERYieldsParameters() { return fER; }
  LArERElectronYieldsAlphaParameters GetERElectronYieldsAlphaParameters() {
    return fER.alpha;
  }
  LArERElectronYieldsBetaParameters GetERElectronYieldsBetaParameters() {
    return fER.beta;
  }
  LArERElectronYieldsGammaParameters GetERElectronYieldsGammaParameters() {
    return fER.gamma;
  }
  LArERElectronYieldsDokeBirksParameters
  GetERElectronYieldsDokeBirksParameters() {
    return fER.doke_birks;
  }
  ThomasImelParameters GetThomasImelParameters() {
    return fThomasImelParameters;
  }
  DriftParameters GetDriftParameters() { return fDriftParameters; }

  //-------------------------All Yields-------------------------//
  LArYieldResult GetRecombinationYields(double TotalYields,
                                        double ElectronYields,
                                        double PhotonYields, double energy,
                                        double efield);
  /**
   * @brief GetYields computes the quanta N_ex, N_ion, N_gamma,
   * and N_e- using one of two methods, either the first principles
   * approach (which calculates N_ex and N_ion first) or the
   * data-driven approach (which gives N_gamma and N_e- directly).
   *
   * @param species
   * @param energy
   * @param dx
   * @param efield
   * @param density
   * @return LArYieldResult
   */
  LArYieldResult GetYields(LArInteraction species, double energy, double dx,
                           double efield, double density);
  /**
   * @brief Calculate fluctions on the mean yields
   *
   */
  LArYieldFluctuationResult GetYieldFluctuations(LArInteraction species,
                                                 const LArYieldResult &yields,
                                                 double density);
  /**
   * @brief
   *
   */
  LArNESTResult FullCalculation(LArInteraction species, double energy,
                                double dx, double efield, double density,
                                bool do_times);

  //---------------------Data-driven methods--------------------//
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
  LArYieldResult GetNRYields(double energy, double efield, double density);

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
  LArYieldResult GetERYields(double energy, double efield, double density);
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
  LArYieldResult GetAlphaYields(double energy, double efield, double density);

  //---------------------First-principles methods--------------------//
  double GetCanonicalTotalYields(double energy);
  double GetCanonicalIonizationYields(double energy);

  //-----------------------------LET Yields-----------------------------//
  LArYieldResult GetLeptonLETYields(double energy, double dx, double efield,
                                    double density);
  LArYieldResult GetLETYields(double energy, double dx, double efield,
                              double density);
  double GetLETRecombinationProbability(double LET, double efield);
  LArYieldResult GetLETRecombinationYields(double ionization_yields,
                                           double exciton_yields, double energy,
                                           double LET, double efield);

  //----------------------------dE/dx Yields----------------------------//
  LArYieldResult GetdEdxRecombinationYields(double total_yields,
                                            double ionization_yields,
                                            double energy, double dx,
                                            double efield);
  double GetdEdxRecombinationProbability(double dEdx, double efield);
  /**
   * @brief Calculate yields based on the Thomas-Imel model using dE/dx.
   * This should only be used for particles which have a dE/dx, other
   * particle types should use the point-like models NR/ER/Alpha.
   *
   * @return LArYieldResult
   */
  LArYieldResult GetdEdxYields(double energy, double dx, double efield,
                               double density);
  //-------------------------Fluctuation Yields-------------------------//
  LArYieldFluctuationResult GetDataDrivenFluctuations(
      const LArYieldResult &yields, double density);
  /**
   * @brief Get the Default Fluctuations object
   *
   * @param yields
   * @param density
   * @return LArYieldFluctuationResult
   */
  LArYieldFluctuationResult GetDefaultFluctuations(const LArYieldResult &yields,
                                                   double density);
  //-------------------------Photon Times-------------------------//
  double GetPhotonTime(LArInteraction species, bool exciton, double energy);
  double GetPhotonEnergy(bool state);

  //-------------------------Drift Velocity-------------------------//
  double GetDriftVelocity_Liquid(double Kelvin, double eField);
  double GetDriftVelocity_MagBoltz(double density, double efieldinput,
                                   double molarMass);

  //-------------------------Utilities-------------------------//
  double GetLinearEnergyTransfer(double energy, bool CSDA = false);
  double GetDensity(double Kelvin, double bara, bool &inGas, uint64_t evtNum,
                    double molarMass);
  std::vector<double> CalculateG2(int verbosity = -1);

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
  LArYieldResult LegacyGetYields(double energy, double efield,
                                 double yieldFactor, double excitationRatio,
                                 double epsilon, double recombProb);
  /**
   * @brief Below are legacy LAr calculation
   * functions which are almost copied verbatim
   * from LArSoft, see -
   *
   * The default value for the dx=track_length is
   * 1/3 mm, which is the standard from Geant4.
   * @return LArYieldResult
   */
  LArYieldResult LegacyCalculation(int pdgcode, double energy, double efield,
                                   double density,
                                   double track_length = 0.0003);
  double LegacyGetRecombinationProbability(double energy, double efield,
                                           double density, int pdgcode,
                                           double track_length);
  double LegacyGetLinearEnergyTransfer(double E);

 private:
  bool fUseDokeBirks = {false};

  double fDensity = {1.393};
  double fRIdealGas = {8.31446261815324};
  double fRealGasA = {0.1355};    // m^6*Pa/mol^2 or m^4*N/mol^2.
  double fRealGasB = {3.201e-5};  // m^3/mol.
  double fTemperature = {85.0};   // K

  double fWorkQuantaFunction = {19.5};
  double fWorkIonFunction = {23.6};
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
}  // namespace NEST
