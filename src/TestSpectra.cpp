/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   TestSpectra.cpp
 * Author: brodsky3
 *
 * Created on December 11, 2017, 10:27 AM
 */

#include <exception>
#include "TestSpectra.hh"
#include "GammaHandler.hh"
#include <stdexcept>

using namespace std;

double power =
    3.7488;  // this is a global variable because it is for both AmBe and 252Cf

const vector<double> TestSpectra::Gamma_spectrum(double xMin, double xMax,
                                                 string source) {
  GammaHandler gh;

  return gh.combineSpectra(xMin, xMax, source);
}

double TestSpectra::CH3T_spectrum(double xMin, double xMax) {
  double m_e = ElectronRestMassEnergy;  // e- rest mass-energy [keV]
  double aa = 0.0072973525664;          // fine structure constant
  double ZZ = 2.;
  double qValue = 18.5898;  // tritium beta decay endpoint [keV]

  if (xMax > qValue) xMax = qValue;
  if (xMin < 0.) xMin = 0.;
  if (xMin != 0. || xMax != qValue)
    cerr << "WARNING: Recommended energy range is 0 to " << qValue << " keV"
         << endl;
  double yMax = 1.1e7;  // top of the beta decay E histogram
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double B =
        sqrt(xyTry[0] * xyTry[0] + 2. * xyTry[0] * m_e) / (xyTry[0] + m_e);
    double x = (2. * M_PI * ZZ * aa) * (xyTry[0] + m_e) /
               sqrt(xyTry[0] * xyTry[0] + 2. * xyTry[0] * m_e);
    double FuncValue = (sqrt(2. * xyTry[0] * m_e) * (xyTry[0] + m_e) *
                        (qValue - xyTry[0]) * (qValue - xyTry[0]) * x *
                        (1. / (1. - exp(-x))) * (1.002037 - 0.001427 * (B)));
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }
  return xyTry[0];
}

double TestSpectra::C14_spectrum(double xMin, double xMax) {
  double m_e = ElectronRestMassEnergy;  // e- rest mass-energy [keV]
  double aa = 0.0072973525664;          // fine structure constant
  double ZZ = 7.;
  double V0 = 0.495;  // effective offset in T due to screening of the nucleus
                      // by electrons
  double qValue = 156.;  // C14 beta decay endpoint [keV]

  if (xMax > qValue) xMax = qValue;
  if (xMin < 0.) xMin = 0.;
  if (xMin != 0. || xMax != qValue)
    cerr << "WARNING: Recommended energy range is 0 to " << qValue << " keV"
         << endl;
  double yMax = 2.5e9;  // top of the beta decay E histogram
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double Ee = xyTry[0] + m_e;             // Total energy of electron
    double pe = sqrt(Ee * Ee - m_e * m_e);  // momentum of the electron
    // phase space part of spectrum
    double dNdE_phasespace =
        pe * Ee * (qValue - xyTry[0]) * (qValue - xyTry[0]);

    // Fermi function (Bethe-Bacher approximation)
    double Ee_screen = Ee - V0;
    double W_screen = (Ee_screen) / m_e;
    double p_screen = sqrt(W_screen * W_screen - 1);
    double WW = (Ee) / m_e;
    double pp = sqrt(WW * WW - 1);
    double G_screen = (Ee_screen) / (m_e);  // Gamma, Total energy(KE+M) over M
    double B_screen = sqrt((G_screen * G_screen - 1) /
                           (G_screen * G_screen));  // v/c of electron. Ratio of
                                                    // velocity to speed of
                                                    // light in vacuum.
    double x_screen = (2 * M_PI * ZZ * aa) / B_screen;
    double F_nr_screen =
        W_screen * p_screen / (WW * pp) * x_screen * (1 / (1 - exp(-x_screen)));
    double F_bb_screen =
        F_nr_screen *
        pow(W_screen * W_screen * (1 + 4 * (aa * ZZ) * (aa * ZZ)) - 1,
            sqrt(1 - aa * aa * ZZ * ZZ) - 1);

    double FuncValue = dNdE_phasespace * F_bb_screen;
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }
  return xyTry[0];
}

double TestSpectra::B8_spectrum(double xMin, double xMax) {
  if (xMax != 4.) xMax = 4.;
  if (xMin != 0.) xMin = 0.;
  double yMax = pow(10., -2.198);
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double FuncValue = 2.198 + 1.2184 * xyTry[0] - 0.32849 * pow(xyTry[0], 2.) +
                       0.12441 * pow(xyTry[0], 3.);
    FuncValue = pow(10., -FuncValue);
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }
  return xyTry[0];
}

double TestSpectra::AmBe_spectrum(double xMin, double xMax) {
  if (xMax > 200.) xMax = 200.;
  if (xMin < DBL_MIN) xMin = DBL_MIN;
  double yMax = pow(10., power), yMin = 0.0;
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double FuncValue =
        power * pow(log10(xyTry[0]), 0.) - 0.77942 * pow(log10(xyTry[0]), 1.) +
        1.30300 * pow(log10(xyTry[0]), 2.) -
        2.75280 * pow(log10(xyTry[0]), 3.) +
        1.57310 * pow(log10(xyTry[0]), 4.) - 0.30072 * pow(log10(xyTry[0]), 5.);
    FuncValue = pow(10., FuncValue);
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, yMin, yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }

  return xyTry[0];
}

double TestSpectra::Cf_spectrum(double xMin, double xMax) {
  if (xMax > 200.) xMax = 200.;
  if (xMin < DBL_MIN) xMin = DBL_MIN;
  double yMax = 2. * pow(10., power), yMin = 0.0;
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double FuncValue =
        power * pow(log10(xyTry[0]), 0.) - 0.77942 * pow(log10(xyTry[0]), 1.) +
        1.30300 * pow(log10(xyTry[0]), 2.) -
        2.75280 * pow(log10(xyTry[0]), 3.) +
        1.57310 * pow(log10(xyTry[0]), 4.) - 0.30072 * pow(log10(xyTry[0]), 5.);
    FuncValue = pow(10., FuncValue);
    FuncValue *= 1.9929 - .033214 * pow(xyTry[0], 1.) +
                 .00032857 * pow(xyTry[0], 2.) - 1.000e-6 * pow(xyTry[0], 3.);
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, yMin, yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }

  return xyTry[0];
}

double TestSpectra::DD_spectrum(
    double xMin, double xMax) {  // JV LUX, most closely like JENDL-4. See
                                 // arXiv:1608.05381. Lower than G4/LUXSim
  if (xMax > 80.) xMax = 80.;
  if (xMin < 0.000) xMin = 0.000;
  double yMax = 1.1;
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double FuncValue =
        exp(-xyTry[0] / 10.) + 0.1 * exp(-pow((xyTry[0] - 60.) / 25., 2.));
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }
  return xyTry[0];
}

double TestSpectra::ppSolar_spectrum(double xMin, double xMax) {
  if (xMax > 250.) xMax = 250.;
  if (xMin < 0.00) xMin = 0.00;
  double yMax = 0.000594;
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double FuncValue =
        9.2759e-6 - 3.5556e-8 * xyTry[0] - 5.4608e-12 * xyTry[0] * xyTry[0];
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., yMax, xyTry[0],
                                          xyTry[1], FuncValue);
  }
  return xyTry[0];
}

double TestSpectra::atmNu_spectrum(double xMin, double xMax) {
  if (xMax > 85.) xMax = 85.;
  if (xMin < 0.0) xMin = 0.0;
  vector<double> xyTry = {
      xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {
    double FuncValue =
        (1. + 0.0041482 * xyTry[0] + 0.00079972 * xyTry[0] * xyTry[0] -
         1.0201e-5 * xyTry[0] * xyTry[0] * xyTry[0]) *
        exp(-xyTry[0] / 12.355);
    xyTry = RandomGen::rndm()->VonNeumann(xMin, xMax, 0., 1., xyTry[0],
                                          xyTry[1], FuncValue);
  }
  return xyTry[0];
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
// dR() //generator written by Vic Gehman originally
//------++++++------++++++------++++++------++++++------++++++------++++++------

// This spectrum comes from Phys. Rev. D 82 (2010) 023530 (McCabe). It's PreGAIA
double TestSpectra::WIMP_dRate(double ER, double mWimp, double dayNum) {
  // We are going to hard code in the astrophysical halo for now.  This may be
  // something that we make an argument later, but this is good enough to start.
  // Some constants:
  double M_N = 0.9395654;                  // Nucleon mass [GeV]
  double N_A = NEST_AVO;                   // Avogadro's number [atoms/mol]
  double c = 2.99792458e10;                // Speed of light [cm/s]
  double GeVperAMU = 0.9315;               // Conversion factor
  double SecondsPerDay = 60. * 60. * 24.;  // Conversion factor
  double KiloGramsPerGram = 0.001;         // Conversion factor
  double keVperGeV = 1.e6;                 // Conversion factor
  double cmPerkm = 1.e5;                   // Conversion factor
  double SqrtPi = pow(M_PI, 0.5);
  double root2 = sqrt(2.);
  // Convert all velocities from km/s into cm/s
  double v_0 = V_WIMP * cmPerkm;      // peak WIMP velocity
  double v_esc = V_ESCAPE * cmPerkm;  // escape velocity
  double v_e =
      (V_SUN + (0.49 * 29.8 *
                cos((dayNum * 2. * M_PI / 365.24) - (0.415 * 2. * M_PI)))) *
      cmPerkm;  // the Earth's velocity
  // used Eq. 18 for SHM w/ June 1 as reference date (MAX!) from arXiv 0607121
  // [Savage, Freese, Gondolo 2006] - Juergen Reichenbacher 09/17/2020

  // Define the detector Z and A and the mass of the target nucleus
  double Z = ATOM_NUM;
  double A = (double)RandomGen::rndm()->SelectRanXeAtom();
  double M_T = A * GeVperAMU;

  // Calculate the number of target nuclei per kg
  double N_T = N_A / (A * KiloGramsPerGram);

  // Rescale the recoil energy and the inelastic scattering parameter into GeV
  ER /= keVperGeV;
  double delta = 0. / keVperGeV;  // Setting this to a nonzero value will allow
  // for inelastic dark matter...
  // Set up your dummy WIMP model (this is just to make sure that the numbers
  // came out correctly for definite values of these parameters, the overall
  // normalization of this spectrum doesn't matter since we generate a definite
  // number of events from the macro).
  double m_d = mWimp;       // [GeV]
  double sigma_n = 1.e-36;  //[cm^2] 1 pb reference
  // Calculate the other factors in this expression
  double mu_ND = mWimp * M_N / (mWimp + M_N);  // WIMP-nucleON reduced mass
  double mu_TD = mWimp * M_T / (mWimp + M_T);  // WIMP-nucleUS reduced mass
  double fp =
      1.;  // Neutron and proton coupling constants for WIMP interactions.
  double fn = 1.;

  // Calculate the minimum velocity required to give a WIMP with energy ER
  double v_min = 0.;
  if (ER != 0.) {
    v_min = c * (((M_T * ER) / mu_TD) + delta) / (root2 * sqrt(M_T * ER));
  }
  double bet = 1.;

  // Start calculating the differential rate for this energy bin, starting
  // with the velocity integral:
  double x_min = v_min / v_0;  // Use v_0 to rescale the other velocities
  double x_e = v_e / v_0;
  double x_esc = v_esc / v_0;
  // Calculate overall normalization to the velocity integral
  double N = SqrtPi * SqrtPi * SqrtPi * v_0 * v_0 * v_0 *
             (erf(x_esc) - (4. / SqrtPi) * exp(-x_esc * x_esc) *
                               (x_esc / 2. + bet * x_esc * x_esc * x_esc / 3.));
  // Calculate the part of the velocity integral that isn't a constant
  double zeta = 0.;
  int thisCase = -1;
  if ((x_e + x_min) < x_esc) {
    thisCase = 1;
  }
  if ((x_min > std::abs(x_esc - x_e)) && ((x_e + x_esc) > x_min)) {
    thisCase = 2;
  }
  if (x_e > (x_min + x_esc)) {
    thisCase = 3;
  }
  if ((x_e + x_esc) < x_min) {
    thisCase = 4;
  }
  switch (thisCase) {
    case 1:
      zeta = ((SqrtPi * SqrtPi * SqrtPi * v_0 * v_0) / (2. * N * x_e)) *
             (erf(x_min + x_e) - erf(x_min - x_e) -
              ((4. * x_e) / SqrtPi) * exp(-x_esc * x_esc) *
                  (1 + bet * (x_esc * x_esc - x_e * x_e / 3. - x_min * x_min)));
      break;
    case 2:
      zeta = ((SqrtPi * SqrtPi * SqrtPi * v_0 * v_0) / (2. * N * x_e)) *
             (erf(x_esc) + erf(x_e - x_min) -
              (2. / SqrtPi) * exp(-x_esc * x_esc) *
                  (x_esc + x_e - x_min -
                   (bet / 3.) * (x_e - 2. * x_esc - x_min) *
                       (x_esc + x_e - x_min) * (x_esc + x_e - x_min)));
      break;
    case 3:
      zeta = 1. / (x_e * v_0);
      break;
    case 4:
      zeta = 0.;
      break;
    default:
      throw std::runtime_error(
          "\tThe velocity integral in the WIMP generator broke!!!");
  }

  double a = 0.52;                           // in fm
  double C = 1.23 * pow(A, 1. / 3.) - 0.60;  // fm
  double s = 0.9;  // skin depth of nucleus in fm. Originally used by Karen
                   // Gibson; XENON100 1fm; 2.30 acc. to Lewin and Smith maybe?
  double rn = sqrt(C * C + (7. / 3.) * M_PI * M_PI * a * a -
                   5. * s * s);    // alternatives: 1.14*A^1/3 given in L&S, or
                                   // rv=1.2*A^1/3 then rn =
                                   // sqrt(pow(rv,2.)-5.*pow(s,2.)); used by
                                   // XENON100 (fm)
  double q = 6.92 * sqrt(A * ER);  // in units of 1 over distance or length
  double FormFactor;
  if (q * rn > 0.)
    FormFactor =
        3. * exp(-0.5 * q * q * s * s) * (sin(q * rn) - q * rn * cos(q * rn)) /
        (q * rn * q * rn * q * rn);  // qr and qs unitless inside Bessel
                                     // function, which is dimensionless too
  else
    FormFactor = 1.;

  // Now, the differential spectrum for this bin!
  double dSpec = 0.5 * (c * c) * N_T * (RHO_NAUGHT / m_d) *
                 (M_T * sigma_n / (mu_ND * mu_ND));
  // zeta=1.069-1.4198*ER+.81058*pow(ER,2.)-.2521*pow(ER,3.)+.044466*pow(ER,4.)-0.0041148*pow(ER,5.)+0.00013957*pow(ER,6.)+2.103e-6*pow(ER,7.);
  // if ( ER > 4.36 ) squiggle = 0.; //parameterization for 7 GeV WIMP using
  // microMegas
  dSpec *= (((Z * fp) + ((A - Z) * fn)) / fn) *
           (((Z * fp) + ((A - Z) * fn)) / fn) * zeta * FormFactor * FormFactor *
           SecondsPerDay / keVperGeV;

  return dSpec;
}

TestSpectra::WIMP_spectrum_prep TestSpectra::WIMP_prep_spectrum(double mass,
                                                                double eStep,
                                                                double dayNum) {
  WIMP_spectrum_prep spectrum;
  double divisor, x1, x2;
  vector<double> EnergySpec;
  int numberPoints;
  if (mass < 2.0) {  // GeV/c^2
    divisor = 10. / eStep;
    numberPoints = int(1000. / eStep);
  } else if (mass < 10.) {
    divisor = 10. / eStep;
    numberPoints = int(1000. / eStep);
  } else {
    divisor = 1.0 / eStep;
    numberPoints = int(100. / eStep);
  }
  int nZeros = 0;  // keep track of the number of zeros in a row
  for (int i = 0; i < (numberPoints + 1); ++i) {
    EnergySpec.push_back(WIMP_dRate(double(i) / divisor, mass, dayNum));
    if (ValidityTests::nearlyEqual(EnergySpec[i], 0.))
      ++nZeros;
    else
      nZeros = 0;  // reset the count if EnergySpec[i] != zero
    if (nZeros == 100)
      break;  // quit the for-loop once we're sure we're only getting zeros
  }

  for (uint64_t i = 0; i < 1000000; ++i) {
    spectrum.integral += WIMP_dRate(double(i) / 1e4, mass, dayNum) / 1e4;
  }
  spectrum.xMax = ((double)EnergySpec.size() - 1.) / divisor;
  // defualt value -- will be overwritten if
  // xMax is acutally smaller
  for (int i = 0; i < (int)EnergySpec.size() - 1; ++i) {
    x1 = double(i) / divisor;
    x2 = double(i + 1) / divisor;
    spectrum.base[i] = EnergySpec[i + 1] *
                       pow(EnergySpec[i + 1] / EnergySpec[i], x2 / (x1 - x2));
    spectrum.exponent[i] = log(EnergySpec[i + 1] / EnergySpec[i]) / (x1 - x2);
    if (spectrum.base[i] > 0. && spectrum.base[i] < DBL_MAX &&
        spectrum.exponent[i] > 0. && spectrum.exponent[i] < DBL_MAX)
      ;  // spectrum.integral+=spectrum.base[i]/spectrum.exponent[i]*(exp(-spectrum.exponent[i]*x1)-exp(-spectrum.exponent[i]*x2));
    else {
      if (EnergySpec[i + 1] >
          10.) {  // i.e. the calculation stopped before event rate was low
        throw std::runtime_error(
            "ERROR: WIMP E_step is too small (or large)! Increase(decrease) it "
            "slightly to avoid noise in the calculation.");
      }
      spectrum.xMax = double(i - 1) / divisor;
      if (spectrum.xMax <= 0.0) {
        throw std::runtime_error(
            "ERROR: The maximum possible WIMP recoil is not +-ive, which "
            "usually means your E_step is too small (OR it is too large).");
      }
      break;
    }
  }

  spectrum.divisor = divisor;
  return spectrum;
}

double TestSpectra::WIMP_spectrum(WIMP_spectrum_prep wimp_spectrum, double mass,
                                  double dayNum) {
  int count = 0;  // added by Jack Genovesi of LZ
  double xMin = 0., FuncValue = 0.00, x = 0.;
  double yMax = WIMP_dRate(xMin, mass, dayNum);
  vector<double> xyTry = {
      xMin + (wimp_spectrum.xMax - xMin) * RandomGen::rndm()->rand_uniform(),
      yMax * RandomGen::rndm()->rand_uniform(), 1.};
  while (xyTry[2] > 0.) {  // start outer while loop
    while (
        xyTry[1] >
        (-WIMP_dRate(0., mass, dayNum) / wimp_spectrum.xMax * xyTry[0] +
         WIMP_dRate(0., mass,
                    dayNum))) {  // triangle cut more efficient than rectangle
      xyTry[0] =
          (wimp_spectrum.xMax - xMin) * RandomGen::rndm()->rand_uniform();
      xyTry[1] =
          yMax *
          RandomGen::rndm()->rand_uniform();  // needs to be a lower value
    }                                         // end of the inner while loop
    for (x = 0; x < wimp_spectrum.xMax;
         x += (1. / wimp_spectrum.divisor)) {  // start inner for loop
      if (xyTry[0] > x &&
          xyTry[0] <
              (x + 1. / wimp_spectrum.divisor)) {  // start inner if statement
        FuncValue =
            wimp_spectrum.base[int(x * wimp_spectrum.divisor)] *
            exp(-wimp_spectrum.exponent[int(x * wimp_spectrum.divisor)] *
                xyTry[0]);
        break;
      }  // end inner if statement
    }    // end inner for loop
    xyTry = RandomGen::rndm()->VonNeumann(xMin, wimp_spectrum.xMax, 0., yMax,
                                          xyTry[0], xyTry[1], FuncValue);

    ++count;  // for avoiding an infinite loop
    if (count >= 100) {
      xyTry[0] = 0.;
      break;
    }
  }

  return xyTry[0];
}

double TestSpectra::ZeplinBackground() {  // Z3 FSR ex.

  double selector = RandomGen::rndm()->rand_uniform();
  double selEnerg;

  if (selector > 0.000000 && selector <= 0.038602)
    selEnerg =
        1.0482 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.038602 && selector <= 0.081630)
    selEnerg =
        1.1494 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.081630 && selector <= 0.085197)
    selEnerg =
        1.2603 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.085197 && selector <= 0.098211)
    selEnerg =
        1.3820 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.098211 && selector <= 0.116010)
    selEnerg =
        1.5153 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.116010 && selector <= 0.134960)
    selEnerg =
        1.6616 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.134960 && selector <= 0.181840)
    selEnerg =
        1.8219 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.181840 && selector <= 0.215600)
    selEnerg =
        1.9977 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.215600 && selector <= 0.250500)
    selEnerg =
        2.1905 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.250500 && selector <= 0.280450)
    selEnerg =
        2.4019 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.280450 && selector <= 0.307760)
    selEnerg =
        2.6337 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.307760 && selector <= 0.335780)
    selEnerg =
        2.8879 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.335780 && selector <= 0.362760)
    selEnerg =
        3.1665 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.362760 && selector <= 0.404200)
    selEnerg =
        3.4721 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.404200 && selector <= 0.437260)
    selEnerg =
        3.8072 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.437260 && selector <= 0.459880)
    selEnerg =
        4.1746 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.459880 && selector <= 0.493280)
    selEnerg =
        4.5775 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.493280 && selector <= 0.527320)
    selEnerg =
        5.0192 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.527320 && selector <= 0.548560)
    selEnerg =
        5.5036 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.548560 && selector <= 0.577610)
    selEnerg =
        6.0347 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.577610 && selector <= 0.609550)
    selEnerg =
        6.6171 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.609550 && selector <= 0.635570)
    selEnerg =
        7.2556 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.635570 && selector <= 0.656480)
    selEnerg =
        7.9558 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.656480 && selector <= 0.689470)
    selEnerg =
        8.7236 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.689470 && selector <= 0.720960)
    selEnerg =
        9.5654 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.720960 && selector <= 0.749250)
    selEnerg =
        10.489 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.749250 && selector <= 0.779750)
    selEnerg =
        11.501 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.779750 && selector <= 0.814330)
    selEnerg =
        12.611 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.814330 && selector <= 0.842290)
    selEnerg =
        13.828 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.842290 && selector <= 0.878470)
    selEnerg =
        15.162 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.878470 && selector <= 0.908490)
    selEnerg =
        16.625 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.908490 && selector <= 0.939570)
    selEnerg =
        18.230 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.939570 && selector <= 0.971280)
    selEnerg =
        19.989 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else if (selector > 0.971280 && selector < 1.0000000)
    selEnerg =
        21.918 *
        (1. + 0.08801 / 2. * (2. * RandomGen::rndm()->rand_uniform() - 1.));
  else
    selEnerg = RandomGen::rndm()->rand_uniform() * 20.;

  return selEnerg;  // selection under the curve is made
}
