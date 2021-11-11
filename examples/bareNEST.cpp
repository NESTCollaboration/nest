/*
 * File:   bareNEST.cpp
 * Author: Jacob Cutter
 *
 * Created on June 22, 2018
 *
 * This skeleton program is meant to demonstrate the bare-bones usage of NEST
 * functions.
 * For the sake of simplicity, this code does not implement any of execNEST's
 * statistical
 * functionality (bands, energy resolution, etc.). This will generate a SINGLE
 * event.
 */

// Include necessary NEST classes
#include "NEST.hh"
#include "TestSpectra.hh"
#include "analysis.hh"

// Include user-specific detector header, located in Detectors/
#include "DetectorExample_XENON10.hh"

using namespace std;
using namespace NEST;

// Skeleton function for doing a simple NEST calculation
int main(int argc, char** argv) {
  double dayNum = 0.;
  // Give a message establishing the use of this code
  cerr << "NOTE: This is a skeleton code meant to be a starting point for "
          "custom uses of NEST. But you should really look at execNEST."
       << endl
       << endl;

  // Instantiate your own VDetector class here, then load into NEST class
  // constructor
  DetectorExample_XENON10* detector = new DetectorExample_XENON10();

  // Construct NEST class using detector object
  NEST::NESTcalc n(detector);

  // Set random seed of RandomGen class
  RandomGen::rndm()->SetSeed(0);

  // Construct NEST objects for storing calculation results
  NEST::YieldResult yields;
  NEST::QuantaResult quanta;

  // Declare needed temporary variables
  vector<double> vTable, NuisParam = {11., 1.1, 0.0480, -0.0533, 12.6, 0.3,
                                      2.,  0.3, 2.,     0.5,     1.,   1.};
  int index;  // index for Z step (for getting pre-calculated drift field)
  double g2, pos_x, pos_y, pos_z, r, phi, driftTime, field, vD, vD_middle;
  // Energy min and max for source spectrum
  double eMin = 10;
  double eMax = 100;
  // For ion calculations
  double atomNum = 0, massNum = 0;
  // Only necessary for type_num == WIMP
  double wimp_mass_GeV = 1.;

  // Choose a possible interaction type
  INTERACTION_TYPE type_num = NR;
  // INTERACTION_TYPE type_num = WIMP;
  // INTERACTION_TYPE type_num = B8;
  // INTERACTION_TYPE type_num = DD;
  // INTERACTION_TYPE type_num = AmBe;
  // INTERACTION_TYPE type_num = Cf;
  // NOTE: For heavy nuclei or ions, ion requires specification of "atomNum" and
  // "massNum"
  // INTERACTION_TYPE type_num = ion;
  // INTERACTION_TYPE type_num = gammaRay;
  // INTERACTION_TYPE type_num = Kr83m;
  // INTERACTION_TYPE type_num = CH3T;
  // INTERACTION_TYPE type_num = C14;
  // INTERACTION_TYPE type_num = NEST::beta;

  // Do checks for special-case source spectra
  if (type_num == Kr83m) {
    if (eMin == 9.4 && eMax == 9.4) {
    } else if (eMin == 32.1 && eMax == 32.1) {
    } else {
      if (verbosity)
        cerr << "ERROR: For Kr83m, put both energies as 9.4 or both as 32.1 "
                "keV please."
             << endl;
      return 1;
    }
  } else if (type_num == gammaRay && verbosity) {
    if (eMin < 10. || eMax < 10.) {
      cerr << "WARNING: Typically beta model works better for ER BG at low "
              "energies as in a WS."
           << endl;
      cerr << "ER data is often best matched by a weighted average of the beta "
              "& gamma models."
           << endl;
    }
  }

  // Draw an event energy from the appropriate spectrum
  TestSpectra spec;
  double keV = -999.;
  if (eMin == eMax && eMin >= 0. && eMax > 0.)
    keV = eMin;
  else {
    switch (type_num) {
      case CH3T:
        keV = spec.CH3T_spectrum(eMin, eMax);
        break;
      case C14:
        keV = spec.C14_spectrum(eMin, eMax);
        break;
      case B8:  // normalize this to ~3500 / 10-ton / year, for E-threshold of
                // 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
        keV = spec.B8_spectrum(eMin, eMax);
        break;
      case AmBe:  // for ZEPLIN-III FSR from HA (Pal '98)
        keV = spec.AmBe_spectrum(eMin, eMax);
        break;
      case Cf:
        keV = spec.Cf_spectrum(eMin, eMax);
        break;
      case DD:
        keV = spec.DD_spectrum(eMin, eMax);
        break;
      case WIMP:
        spec.wimp_spectrum_prep =
            spec.WIMP_prep_spectrum(wimp_mass_GeV, E_step, dayNum);
        keV =
            spec.WIMP_spectrum(spec.wimp_spectrum_prep, wimp_mass_GeV, dayNum);
        break;
      default:
        keV = eMin + (eMax - eMin) * RandomGen::rndm()->rand_uniform();
        break;
    }
  }

  // For most sources, ensure that the energy is in bounds
  if (type_num != WIMP && type_num != B8) {
    if (keV > eMax) keV = eMax;
    if (keV < eMin) keV = eMin;
  }

  // Calculate noble density based on temperature and pressure.
  // Use this to determine whether we are in gas phase
  double rho = n.SetDensity(detector->get_T_Kelvin(),
                            detector->get_p_bar());  // cout.precision(12);
  if (rho < 1.) detector->set_inGas(true);

  // Calculate and print g1, g2 parameters (once per detector)
  vector<double> g2_params = n.CalculateG2();
  g2 = g2_params.back();

  // Calculate a drift velocity table for non-uniform fields,
  // and calculate the drift velocity at detector center for normalization
  // purposes
  vTable = n.SetDriftVelocity_NonUniform(rho, z_step, pos_x, pos_y);
  vD_middle = vTable[int(
      floor(.5 * (detector->get_gate() - 100. + detector->get_cathode() + 1.5) /
                z_step +
            0.5))];

  // Regenerate a random position for the event until the corresponding drift
  // time is within bounds (drift time uses the cumulative Z-dependent drift
  // velocity)
  driftTime = -999;
  while (driftTime < detector->get_dt_min() ||
         driftTime > detector->get_dt_max() || pos_z <= 0 ||
         pos_z > detector->get_TopDrift()) {
    pos_z = 0. +
            (detector->get_TopDrift() - 0.) * RandomGen::rndm()->rand_uniform();
    r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
    phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
    pos_x = r * cos(phi);
    pos_y = r * sin(phi);
    field = detector->FitEF(pos_x, pos_y, pos_z);
    index = int(floor(pos_z / z_step + 0.5));
    vD = vTable[index];
    driftTime =
        (detector->get_TopDrift() - pos_z) / vD;  // (mm - mm) / (mm / us) = us
  }

  // Get yields from NEST calculator, along with number of quanta
  yields = n.GetYields(type_num, keV, rho, field, double(massNum),
                       double(atomNum), NuisParam);
  vector<double> FreeParam = {1, 1, .1, .5, .19, 2.25};
  quanta = n.GetQuanta(yields, rho, FreeParam);

  // Calculate S2 photons using electron lifetime correction
  double Nphd_S2 =
      g2 * quanta.electrons * exp(-driftTime / detector->get_eLife_us());

  // Vectors for saving times and amplitudes of waveforms (with calculationMode
  // and verbosity boolean flags both set to true in analysis.hh)
  vector<double> wf_amp;
  vector<int64_t> wf_time;

  double truthPos[3] = {pos_x, pos_y, pos_z};
  double smearPos[3] = {pos_x, pos_y, pos_z};

  // Calculate the S1 based on the quanta generated
  vector<double> scint =
      n.GetS1(quanta, truthPos[0], truthPos[1], truthPos[2], smearPos[0],
              smearPos[1], smearPos[2], vD, vD_middle, type_num, 0, field, keV,
              s1CalculationMode, verbosity, wf_time, wf_amp);

  // Take care of gamma-X case for positions below cathode
  if (truthPos[2] < detector->get_cathode()) quanta.electrons = 0;
  vector<double> scint2 =
      n.GetS2(quanta.electrons, truthPos[0], truthPos[1], truthPos[2],
              smearPos[0], smearPos[1], smearPos[2], driftTime, vD, 0, field,
              s2CalculationMode, verbosity, wf_time, wf_amp, g2_params);

  // If using the reconstructed energy, back-calculate energy as measured from
  // yields
  if (!MCtruthE) {
    double Nph, g1 = detector->get_g1(), Ne;
    if (usePD == 0)
      Nph = std::abs(scint[3]) / (g1 * (1. + detector->get_P_dphe()));
    else if (usePD == 1)
      Nph = std::abs(scint[5]) / g1;
    else
      Nph = std::abs(scint[7]) / g1;

    if (usePD == 0)
      Ne = std::abs(scint2[5]) / (g2 * (1. + detector->get_P_dphe()));
    else
      Ne = std::abs(scint2[7]) / g2;

    if (yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0.) {
      keV = (Nph + Ne) * W_DEFAULT * 1e-3 / yields.Lindhard;
    } else
      keV = 0.;
  }

  // Possible outputs from "scint" vector
  // scint[0] = nHits; // MC-true integer hits in same OR different PMTs, NO
  // double phe effect
  // scint[1] = Nphe; // MC-true integer hits WITH double phe effect (Nphe >
  // nHits)
  // scint[2] = pulseArea; // floating real# smeared DAQ pulse areas in phe, NO
  // XYZ correction
  // scint[3] = pulseAreaC; // smeared DAQ pulse areas in phe, WITH XYZ
  // correction
  // scint[4] = Nphd; // same as pulse area, adjusted/corrected *downward* for
  // 2-PE effect (LUX phd units)
  // scint[5] = NphdC; // same as Nphd, but XYZ-corrected
  // scint[6] = spike; // floating real# spike count, NO XYZ correction
  // scint[7] = spikeC; // floating real# spike count, WITH XYZ correction
  // scint[8] = fdetector->get_g1(); // g1 (light collection efficiency in
  // liquid)

  // Possible outputs from "scint2" vector
  // scint2[0] = Nee; // integer number of electrons unabsorbed in liquid then
  // getting extracted
  // scint2[1] = Nph; // raw number of photons produced in the gas gap
  // scint2[2] = nHits; // MC-true integer hits in same OR different PMTs, NO
  // double phe effect
  // scint2[3] = Nphe; // MC-true integer hits WITH double phe effect (Nphe >
  // nHits). S2 has more steps than S1 (e's 1st)
  //
  // If S2 threshold is set to positive (normal mode)
  // scint2[4] = pulseArea; // floating real# smeared DAQ pulse areas in phe, NO
  // XYZ correction
  // scint2[5] = pulseAreaC; // smeared DAQ pulse areas in phe, WITH XYZ
  // correction
  // scint2[6] = Nphd; // same as pulse area, adjusted/corrected *downward* for
  // 2-PE effect (LUX phd units)
  // scint2[7] = NphdC; // same as Nphd, but XYZ-corrected
  //
  // If S2 threshold is set to negative (switches from S2 -> S2 bottom, NOT
  // literally negative)
  // scint2[4] = S2b; // floating real# smeared pulse areas in phe ONLY
  // including bottom PMTs, NO XYZ correction
  // scint2[5] = S2bc; // floating real# smeared pulse areas in phe ONLY
  // including bottom PMTs, WITH XYZ correction
  // scint2[6] = S2b / (1.+fdetector->get_P_dphe()); // same as S2b, but
  // adjusted for 2-PE effect (LUX phd units)
  // scint2[7] = S2bc / (1.+fdetector->get_P_dphe()); // same as S2bc, but
  // adjusted for 2-PE effect (LUX phd units)
  // scint2[8] = g2; // g2 = ExtEff * SE, light collection efficiency of EL in
  // gas gap (from CalculateG2)

  // Print selected outputs a la execNEST
  fprintf(
      stdout,
      "\nE [keV]\t\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\tNe-\tS1 [PE "
      "or phe]\tS1_3Dcor [phd]\tS1c_spike\tNe-Extr\tS2_rawArea "
      "[PE]\tS2_3Dcorr [phd]\n");
  printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%d\t%d\t", keV, field, driftTime,
         truthPos[0], truthPos[1], truthPos[2], quanta.photons,
         quanta.electrons);  // comment this out when below line in
  // If energy > 1 MeV, switch output notation
  if (keV > 1000. || scint[5] > maxS1 || scint2[7] > maxS2) {
    printf("%e\t%e\t%e\t", scint[2], scint[5], scint[7]);
    printf("%lli\t%e\t%e\n", (int64_t)scint2[0], scint2[4], scint2[7]);
  } else {
    printf("%.6f\t%.6f\t%.6f\t", scint[2], scint[5],
           scint[7]);  // see GetS1 inside of NEST.cpp for full explanation of
                       // all 8 scint return vector elements. Sample 3 most
                       // common
    printf("%i\t%.6f\t%.6f\n", (int)scint2[0], scint2[4],
           scint2[7]);  // see GetS2 inside of NEST.cpp for full explanation of
                        // all 8 scint2 vector elements. Change as you desire
  }

  return 0;
}
