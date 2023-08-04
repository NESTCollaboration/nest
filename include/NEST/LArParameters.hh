/**
 * @file LArParameters.hh
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
#include <vector>

#include "RandomGen.hh"

namespace NEST {
static constexpr double LAr_Z{18};
static constexpr double legacy_density_LAr{1.393};
static constexpr double legacy_scint_yield{1.0 / (19.5 * 1.e-6)};
static constexpr double legacy_resolution_scale{0.107};  // Doke 1976
static constexpr double two_PI = 2. * M_PI;
static constexpr double sqrt2 = gcem::sqrt(2.);
static constexpr double sqrt2_PI = gcem::sqrt(2. * M_PI);
static constexpr double inv_sqrt2_PI = 1. / gcem::sqrt(2. * M_PI);

enum class LArInteraction {
  NR = 0,
  ER = 1,
  Alpha = 2,
  dEdx = 3,
  LeptonLET = 4,
  LET = 5
};

enum class LArFluctuationModel {
  Default = 0,
};

struct LArNRYieldsParameters {
  double alpha = {11.10};
  double beta = {0.087};
  double gamma = {0.1};
  double delta = {-0.0932};
  double epsilon = {2.998};
  double zeta = {0.3};
  double eta = {2.94};
};
struct LArERElectronYieldsAlphaParameters {
  double A = {32.988};
  double B = {-552.988};
  double C = {17.2346};
  double D = {-4.7};
  double E = {0.025115};
  double F = {0.265360653};
  double G = {0.242671};
};
struct LArERElectronYieldsBetaParameters {
  double A = {0.778482};
  double B = {25.9};
  double C = {1.105};
  double D = {0.4};
  double E = {4.55};
  double F = {-7.502};
};
struct LArERElectronYieldsGammaParameters {
  double A = {0.659509};
  double B = {1000};
  double C = {6.5};
  double D = {5.0};
  double E = {-0.5};
  double F = {1047.408};
  double G = {0.01851};
};
struct LArERElectronYieldsDokeBirksParameters {
  double A = {1052.264};
  double B = {14159350000 - 1652.264};
  double C = {-5.0};
  double D = {0.157933};
  double E = {1.83894};
};
struct LArERYieldsParameters {
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

struct LArAlphaElectronYieldsParameters {
  double A = {1.0 / 6200.0};
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

struct LArAlphaPhotonYieldsParameters {
  double A = {1.16};
  double B = {-0.012};
  double C = {1.0 / 6500.0};
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

struct LArAlphaYieldsParameters {
  LArAlphaElectronYieldsParameters Ye;
  LArAlphaPhotonYieldsParameters Yph;
};

struct ThomasImelParameters {
  double A = {0.1};
  double B = {-0.0932};
};

struct DriftParameters {
  std::vector<double> A = {0.937729,   0.80302379,  0.7795972, 0.6911897,
                           0.76551511, 0.502022794, 0.24207633};
  std::vector<double> B = {-0.0734108, -0.06694564, -0.0990952, -0.092997,
                           -0.0731659, -0.06644517, -0.03558428};
  std::vector<double> C = {0.315338, 0.331798,  0.320876,  0.3295202,
                           0.317972, 0.3290246, 0.33645519};

  std::vector<double> TempLow = {84., 86., 88., 92., 96., 110., 125.};
  std::vector<double> TempHigh = {86., 88., 92., 96., 110., 125., 140.};
};
}  // namespace NEST