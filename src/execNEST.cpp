/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor. CHANGE ..Param @LINES ~180-240 (or 30)
 */

/*
 * File:   execNEST.cpp
 * Author: brodsky3
 *
 * Created on August 1, 2017, 1:03 PM
 */

#include <iterator>
#include "NEST.hh"
#include "TestSpectra.hh"
#include "analysis_G2.hh"
#include "execNEST.hh"

#include "LZ_SR1.hh"

#define tZero 0.00  // day{of the year, 0 is ~Jan. 1}
#define tStep 0.03
#define tMax 365
#define hiEregime 1E+2  // keV

using namespace std;
using namespace NEST;

vector<double> NRERWidthsParam, NRYieldsParam, ERWeightParam,
    ERYieldsParam = default_ERYieldsParam;
double band[NUMBINS_MAX][7], energies[3],
    AnnModERange[2] = {1.5, 6.5};  // keVee or nr (recon)
bool BeenHere = false;
unsigned int SaveTheDates[tMax] = {0};
bool dEOdxBasis = false;
double multFact = 1.;
double minTimeSeparation = 0.;  // ns (Kr83m)

int main(int argc, char** argv) {
  // Instantiate your own VDetector class here, then load into NEST class
  // constructor
  auto* detector = new LZ_Detector();
  if (verbosity > 0) cerr << "*** Detector definition message ***" << endl;
  if (verbosity > 0)
    cerr << "You are currently using the " << detector->getName() << " detector." << endl << endl;
  // Custom parameter modification functions
  // detector->ExampleFunction();

  if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) {
    detector->set_molarMass(39.948);
    if (verbosity > 0)
      cerr << "\nWARNING: Argon is currently only in alpha testing mode!! Many "
              "features copied over from Xenon wholesale still. Use models at "
              "your own risk.\n"
           << endl;
  }

  /* vector<double> eList = { 1., 2., 3. }; // fast example--for PLR, ML train
  vector<vector<double>> pos3dxyz = { {0.,-1.,60.},{-1.,0.,70.},{1.,0.,80.} };
  runNESTvec ( detector, NEST::beta, eList, pos3dxyz );
  delete detector;
  return EXIT_SUCCESS; */

  if (argc < 7) {
    cout << "This program takes 6 (or 7) inputs, with Z position in mm from "
            "bottom of detector:"
         << endl;
    cout << "\t./execNEST numEvts type_interaction E_min[keV] E_max[keV] "
            "field_drift[V/cm] x,y,z-position[mm] {optional:seed}"
         << endl
         << endl;
    cout << "For Kr83m time-dependent 9.4, 32.1, or 41.5 keV yields: " << endl;
    cout << "\t ./execNEST numEvents Kr83m Energy[keV] maxTimeDiff[ns] "
            "field_drift[V/cm] x,y,z-position[mm] {optional:seed}"
         << endl
         << endl;
    cout << "For 8B or pp or atmNu, numEvts is kg-days of exposure with "
            "everything else same. "
            "For WIMPs:"
         << endl;
    cout << "\t./execNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] "
            "field_drift[V/cm] x,y,z-position[mm] {optional:seed}"
         << endl
         << endl;
    cout << "For cosmic-ray muons or other similar particles with elongated "
            "track lengths:"
         << endl;
    cout << "\t./execNEST numEvts {MIP} LET[MeV*cm^2/gram] "
            "x,y,z-position[mm](Initial) field_drift[V/cm] "
            "x,y,z-position[mm](Final) {optional:seed}"
         << endl
         << endl;
    return 1;
  }

  double numEvts;
  string type, position, posiMuon;
  double eMin, eMax, inField, fPos;
  int seed;
  bool no_seed = false;

  if (loopNEST) {
    numEvts = 200000; //10000 for speed
    type = "MIP";
    dEOdxBasis = true;
    eMin = -1.;  // this means do the default, which is CH3T
    s1CalculationMode = NEST::S1CalculationMode::Hybrid;  // Param for speed!
    eMax = -1.;  // for monoE peak: just use eMin (keep neg)
    inField = 180.;
    position = "-1";
    posiMuon = "-1";
    fPos = -1.;
    seed = 0;
    no_seed = false;
    NRERWidthsParam.clear();
    NRYieldsParam.clear();
    ERWeightParam.clear();
    verbosity = -1;

    NRYieldsParam.push_back(0.0);
    NRYieldsParam.push_back(0.);
    NRYieldsParam.push_back(300.);
    NRYieldsParam.push_back(eMax);
    NRYieldsParam.push_back(eMin);
    NRYieldsParam.push_back(atof(argv[1]));
    NRYieldsParam.push_back(inField);
    NRYieldsParam.push_back(1.51);
    NRYieldsParam.push_back(2.888);
    NRYieldsParam.push_back(1e3);
    NRYieldsParam.push_back(atof(argv[2]));
    NRYieldsParam.push_back(-0.5);
    NRYieldsParam.push_back(0);
    detector->set_g1(atof(argv[3]));
    detector->set_noiseLinear(atof(argv[4]),atof(argv[11]));
    
    ERWeightParam.push_back(-atof(argv[5]));
    ERWeightParam.push_back(atof(argv[6]));
    ERWeightParam.push_back(atof(argv[7]));
    ERWeightParam.push_back(atof(argv[8]));
    ERWeightParam.push_back(0.5);
    detector->set_E_gas(atof(argv[9]));
    ERWeightParam.push_back(atof(argv[10]));
    ERWeightParam.push_back(1e3);
    
    NRERWidthsParam.push_back(0.4);  // Leave fluctuations parameters
    NRERWidthsParam.push_back(0.4);  // fixed for now
    NRERWidthsParam.push_back(0.04);
    NRERWidthsParam.push_back(0.50);
    NRERWidthsParam.push_back(0.19);
    NRERWidthsParam.push_back(2.25);
    NRERWidthsParam.push_back(1.);
    NRERWidthsParam.push_back(0.046452);
    NRERWidthsParam.push_back(0.205);
    NRERWidthsParam.push_back(0.45);
    NRERWidthsParam.push_back(-0.2);

  } else {
    numEvts = (double)atof(argv[1]);
    if (numEvts <= 0.) {
      if (verbosity > 0)
        cerr << "ERROR, you must simulate at least one event, or non-zero kg*days"
             << endl;
      return 1;
    }
    type = argv[2];
    eMin = atof(argv[3]);
    eMax = atof(argv[4]);
    inField = atof(argv[5]);
    position = argv[6];
    posiMuon = argv[4];
    fPos = atof(argv[6]);

    seed = 1;  // if not given make 1
    if (argc == 8) {
      multFact = atof(argv[7]);
      seed = atoi(argv[7]);
    } else {
      RandomGen::rndm()->SetSeed(1);
      no_seed = true;
    }
    
    NRERWidthsParam.clear();
    NRYieldsParam.clear();
    ERWeightParam.clear();
    NRERWidthsParam.push_back(0.4);  // NR Fi (Fano factor for ionization)
    NRERWidthsParam.push_back(0.4);  // NR Fex
    NRERWidthsParam.push_back(
        0.04);  // amplitude for non-binomial NR recombination fluctuations
    NRERWidthsParam.push_back(0.50);  // center in e-Frac (NR)
    NRERWidthsParam.push_back(0.19);  // width parameter (Gaussian 1-sigma)
    NRERWidthsParam.push_back(2.25);  // raw skewness, for NR
    NRERWidthsParam.push_back(
        1.);  // ER Fano normalization for non-density dependence
    // negative 0.0015 restores https://arxiv.org/abs/2211.10726v3 Eq. 8
    NRERWidthsParam.push_back(
        0.046452);  // Minimum amplitude for ER non-binom recomb flucts
    NRERWidthsParam.push_back(0.205);  // center in e-frac (ER)
    NRERWidthsParam.push_back(0.45);   // width parameter
    NRERWidthsParam.push_back(-0.2);   // ER non-binom skewness in e-frac

    // if (type == "ER") {  // Based on XELDA L-shell 5.2 keV yields
    // https://arxiv.org/abs/2109.11487
    ERWeightParam.push_back(0.23);   // 0.5 for LUX Run03
    ERWeightParam.push_back(0.77);   // 0.5
    ERWeightParam.push_back(2.95);   // 1.1
    ERWeightParam.push_back(-1.44);  //-5.
    ERWeightParam.push_back(1.0);    // 1.01
    ERWeightParam.push_back(1.0);    // 0.95
    ERWeightParam.push_back(0.);     // 1.4e-2
    ERWeightParam.push_back(0.);     // 1.8e-2

    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) {  // liquid Ar
      NRYieldsParam.push_back(
          11.1025);  // +/-1.10 Everything from
                     // https://docs.google.com/document/d/1vLg8vvY5bcdl4Ah4fzyE182DGWt0Wr7_FJ12_B10ujU
      NRYieldsParam.push_back(1.087399);  // +/-0.025
      NRYieldsParam.push_back(0.1);       // +/-0.005
      NRYieldsParam.push_back(-0.0932);   // +/-0.0095
      NRYieldsParam.push_back(2.998);     // +/-1.026
      NRYieldsParam.push_back(0.3);       // Fixed
      NRYieldsParam.push_back(2.94);      // +/-0.12
      NRYieldsParam.push_back(W_DEFAULT / 1000.);
      NRYieldsParam.push_back(DBL_MAX);
      NRYieldsParam.push_back(0.5);  // square root
      NRYieldsParam.push_back(1.0);
      NRYieldsParam.push_back(1.0);
    } else {
      NRYieldsParam.push_back(
          11.);  // alpha, for NR model. See http://nest.physics.ucdavis.edu
      NRYieldsParam.push_back(1.1);      // beta
      NRYieldsParam.push_back(0.0480);   // gamma
      NRYieldsParam.push_back(-0.0533);  // delta
      NRYieldsParam.push_back(12.6);     // epsilon
      NRYieldsParam.push_back(0.3);      // zeta
      NRYieldsParam.push_back(2.);       // eta
      NRYieldsParam.push_back(0.3);      // theta
      NRYieldsParam.push_back(2.);       // iota
      // last 3 are the secret extra parameters for additional flexibility
      NRYieldsParam.push_back(0.5);  // changes sqrt in Qy equation
      NRYieldsParam.push_back(
          1.0);  // makes low-E sigmoid an asymmetric one, for charge
      NRYieldsParam.push_back(
          1.0);  // makes low-E sigmoid an asymmetric one, for light
    }
  }

  auto exec = execNEST(detector, numEvts, type, eMin, eMax, inField, position,
                       posiMuon, fPos, seed, no_seed, tZero);
  delete detector;
  return exec;
}

NESTObservableArray runNESTvec(
    VDetector* detector,
    INTERACTION_TYPE
        particleType,  // func suggested by Xin Xiang, PD Brown U. for RG, LZ
    vector<double> eList, vector<vector<double>> pos3dxyz, double inField,
    int seed) {
  verbosity = -1;
  NESTcalc n(detector);
  NESTresult result;
  QuantaResult quanta;
  double x, y, z, driftTime, vD;
  RandomGen::rndm()->SetSeed(seed);
  ERYieldsParam = default_ERYieldsParam;
  NRYieldsParam = default_NRYieldsParam;
  NRERWidthsParam = default_NRERWidthsParam;
  vector<double> scint, scint2, wf_amp;
  vector<int64_t> wf_time;
  NESTObservableArray OutputResults;
  double useField;
  vector<double> g2_params = n.CalculateG2(verbosity);
  double rho = n.SetDensity(detector->get_T_Kelvin(), detector->get_p_bar());

  for (uint64_t i = 0; i < eList.size(); ++i) {
    x = pos3dxyz[i][0];
    y = pos3dxyz[i][1];
    z = pos3dxyz[i][2];
    if (inField <= 0.)
      useField = detector->FitEF(x, y, z);
    else
      useField = inField;
    double truthPos[3] = {x, y, z};
    double smearPos[3] = {x, y, z};  // ignoring the difference in this quick
                                     // function caused by smearing
    result = n.FullCalculation(
        particleType, eList[i], rho, useField, detector->get_molarMass(),
        ATOM_NUM, NRYieldsParam, NRERWidthsParam, ERYieldsParam, verbosity);
    quanta = result.quanta;
    vD = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, useField, detector->get_p_bar());
    scint = n.GetS1(quanta, truthPos[0], truthPos[1], truthPos[2], smearPos[0],
                    smearPos[1], smearPos[2], vD, vD, particleType, i, useField,
                    eList[i], NEST::S1CalculationMode::Full, verbosity, wf_time,
                    wf_amp);
    driftTime = (detector->get_TopDrift() - z) /
                vD;  // vD,vDmiddle assumed same (uniform field)
    scint2 = n.GetS2(quanta.electrons, truthPos[0], truthPos[1], truthPos[2],
                     smearPos[0], smearPos[1], smearPos[2], driftTime, vD, i,
                     useField, S2CalculationMode::Full, verbosity, wf_time,
                     wf_amp, g2_params);
    if (scint[7] > PHE_MIN &&
        scint2[7] > PHE_MIN) {  // unlike usual, kill (don't skip, just -> 0)
                                // sub-thr evts
      OutputResults.s1_nhits.push_back(std::abs(int(scint[0])));
      OutputResults.s1_nhits_thr.push_back(std::abs(int(scint[8])));
      OutputResults.s1_nhits_dpe.push_back(std::abs(int(scint[1])));
      OutputResults.s1r_phe.push_back(std::abs(scint[2]));
      OutputResults.s1c_phe.push_back(std::abs(scint[3]));
      OutputResults.s1r_phd.push_back(std::abs(scint[4]));
      OutputResults.s1c_phd.push_back(std::abs(scint[5]));
      OutputResults.s1r_spike.push_back(std::abs(scint[6]));
      OutputResults.s1c_spike.push_back(std::abs(
          scint[7]));  // default is S1c in units of spikes, 3-D XYZ corr
      OutputResults.Nee.push_back(std::abs(int(scint2[0])));
      OutputResults.Nph.push_back(std::abs(int(scint2[1])));
      OutputResults.s2_nhits.push_back(std::abs(int(scint2[2])));
      OutputResults.s2_nhits_dpe.push_back(std::abs(int(scint2[3])));
      OutputResults.s2r_phe.push_back(std::abs(scint2[4]));
      OutputResults.s2c_phe.push_back(std::abs(scint2[5]));
      OutputResults.s2r_phd.push_back(std::abs(scint2[6]));
      OutputResults.s2c_phd.push_back(std::abs(
          scint2[7]));  // default is S2c in terms of phd, not phe a.k.a. PE
    } else {
      OutputResults.s1_nhits.push_back(0);
      OutputResults.s1_nhits_thr.push_back(0);
      OutputResults.s1_nhits_dpe.push_back(0);
      OutputResults.s1r_phe.push_back(0.0);
      OutputResults.s1c_phe.push_back(0.0);
      OutputResults.s1r_phd.push_back(0.0);
      OutputResults.s1c_phd.push_back(0.0);
      OutputResults.s1r_spike.push_back(0.0);
      OutputResults.s1c_spike.push_back(0.0);
      OutputResults.Nee.push_back(0);
      OutputResults.Nph.push_back(0);
      OutputResults.s2_nhits.push_back(0);
      OutputResults.s2_nhits_dpe.push_back(0);
      OutputResults.s2r_phe.push_back(0.);
      OutputResults.s2c_phe.push_back(0.);
      OutputResults.s2r_phd.push_back(0.);
      OutputResults.s2c_phd.push_back(0.);
    }
  }

  return OutputResults;
}

int execNEST(VDetector* detector, double numEvts, const string& type,
             double eMin, double eMax, double inField, string position,
             const string& posiMuon, double fPos, int seed, bool no_seed,
             double dayNumber) {
  // Construct NEST class using detector object
  NESTcalc n(detector);

  if (detector->get_TopDrift() <= 0. || detector->get_anode() <= 0. ||
      detector->get_gate() <= 0.) {
    if (verbosity > 0)
      cerr << "ERROR, unphysical value(s) of position within the detector "
              "geometry.";  // negative or 0 for cathode position is OK (e.g.,
                            // LZ)
    return 1;
  }
  
  if ( !PrintSubThr &&
       ( s1CalculationMode == NEST::S1CalculationMode::Waveform ||
	 s2CalculationMode == NEST::S2CalculationMode::Waveform ) )
    PrintSubThr = true;  // to ensure proper event alignment for PSD work
  
  vector<double> signal1, signal2, signalE, vTable;
  string delimiter, token;
  size_t loc;
  int index;
  double g2, pos_x, pos_y, pos_z, r, phi, driftTime, field, vD,
      vD_middle = 0., atomNum = 0, keVee = 0.0;
  double massNum = detector->get_molarMass();
  YieldResult yieldsMax{};  // for warnings about S1 range

  if (!no_seed) {
    if (seed == -1) {
      RandomGen::rndm()->SetSeed(time(nullptr));
    }

    else {
      RandomGen::rndm()->SetSeed(seed);
    }
  }

  if (ValidityTests::nearlyEqual(eMin, -1.) && type != "muon" &&
      type != "MIP" && type != "LIP" && type != "mu" && type != "mu-")
    eMin = 0.;

  if (ValidityTests::nearlyEqual(eMax, -1.) &&
      ValidityTests::nearlyEqual(eMin, 0.))
    eMax = hiEregime;  // the default energy max
  if (ValidityTests::nearlyEqual(eMax, 0.) && type != "muon" && type != "MIP" &&
      type != "LIP" && type != "mu" && type != "mu-") {
    if (verbosity > 0)
      cerr << "ERROR: The maximum energy (or Kr time sep) cannot be 0 keV (or "
              "0 ns)!"
           << endl;
    return 1;
  }

  INTERACTION_TYPE type_num;
  string gamma_source;

  TestSpectra spec;
  if (type == "NR" || type == "neutron" || type == "-1")
    type_num = NR;  //-1: default particle type is also NR
  else if (type == "WIMP") {
    if (eMin < 0.44) {  // here eMin is WIMP mass
      if (verbosity > 0) cerr << "WIMP mass too low, you're crazy!" << endl;
      return 1;
    }
    type_num = WIMP;
    spec.wimp_spectrum_prep =
        TestSpectra::WIMP_prep_spectrum(eMin, E_step, dayNumber);
    numEvts = (double)RandomGen::rndm()->poisson_draw(
        spec.wimp_spectrum_prep.integral *  // here eMax is cross-section
        1.0 * numEvts * eMax / 1e-36);
  } else if (type == "B8" || type == "Boron8" || type == "8Boron" ||
             type == "8B" || type == "Boron-8") {
    type_num = B8;
    numEvts = (double)RandomGen::rndm()->poisson_draw(0.0026 * numEvts);
  } else if (type == "DD" || type == "D-D")
    type_num = DD;
  else if (type == "AmBe")
    type_num = AmBe;
  else if (type == "Cf" || type == "Cf252" || type == "252Cf" ||
           type == "Cf-252")
    type_num = Cf;
  else if (type == "ion" || type == "nucleus" || type == "alpha") {
    type_num = ion;
    if (type == "alpha") {
      atomNum = 2;
      massNum = 4;
    } else {
      cerr << "Atomic Number: ";
      cin >> atomNum;
      cerr << "Mass Number: ";
      cin >> massNum;
    }
    if (atomNum == ATOM_NUM) type_num = NR;
  } else if (type == "gamma" || type == "gammaRay" || type == "x-ray" ||
             type == "xray" || type == "xRay" || type == "X-ray" ||
             type == "Xray" || type == "XRay")
    type_num = gammaRay;  // includes photo-absorption and electron capture
  else if (type == "Kr83m" || type == "83mKr" || type == "Kr83")
    type_num = Kr83m;
  else if (type == "CH3T" || type == "tritium")
    type_num = CH3T;
  else if (type == "C14" || type == "Carbon14" || type == "14C" ||
           type == "C-14" || type == "Carbon-14")
    type_num = C14;
  else if (type == "beta" || type == "ER" || type == "Compton" ||
           type == "compton" || type == "electron" || type == "e-")
    type_num = NEST::beta;  // default electron recoil model
  else if (type == "muon" || type == "MIP" || type == "LIP" || type == "mu" ||
           type == "mu-") {
    dEOdxBasis = true;
    type_num = NEST::beta;
  } else if (type == "pp" || type == "ppsolar" || type == "ppSolar" ||
             type == "pp_Solar" || type == "pp_solar" || type == "pp-Solar" ||
             type == "pp-solar") {
    type_num = ppSolar;
    numEvts = (double)RandomGen::rndm()->poisson_draw(
        0.0011794 * numEvts);  // normalization: counts per kg-day from 0-250 keV(ee)
  } else if (type == "atmNu" || type == "AtmNu" || type == "atm_Nu" ||
             type == "Atm_Nu" || type == "atm-Nu" || type == "Atm-Nu" ||
             type == "atm_nu" || type == "atm-nu") {
    type_num = atmNu;
    numEvts = (double)RandomGen::rndm()->poisson_draw(1.5764e-7 * numEvts);
  } else if (type == "fullGamma") {
    type_num = fullGamma;
    if (verbosity > 0)
      cerr << "Please choose gamma source. The allowed sources "
              "are:\n\"Co57\"\n\"Co60\"\n\"Cs137\"\nSource: ";
    cin >> gamma_source;
    if (gamma_source == "Co60" && verbosity > 0) {
      cerr << "WARNING: This source is in the pair production range. "
              "Electron/positron pairs are not accounted for after initial "
              "interaction, and some "
           << "photons and electrons may go unaccounted." << endl;
    }
  } else {
    if (verbosity > 0) {
      string particleTypes =
          "UNRECOGNIZED PARTICLE TYPE!! VALID OPTIONS ARE:\n"
          "NR or neutron,\n"
          "WIMP,\n"
          "B8 or Boron8 or 8Boron or 8B or Boron-8,\n"
          "DD or D-D,\n"
          "AmBe,\n"
          "Cf or Cf252 or 252Cf or Cf-252,\n"
          "ion or nucleus,\n"
          "alpha,\n"
          "gamma or gammaRay,\n"
          "x-ray or xray or xRay or X-ray or Xray or XRay,\n"
          "Kr83m or 83mKr or Kr83,\n"
          "CH3T or tritium,\n"
          "Carbon14 or 14C or C14 or C-14 or Carbon-14,\n"
          "beta or ER or Compton or compton or electron or e-,\n"
          "pp or ppSolar with many various underscore, hyphen and "
          "capitalization permutations permitted,\n"
          "atmNu,\n"
          "muon or MIP or LIP or mu or mu-, and\n"
          "fullGamma\n";
      copy(particleTypes.begin(), particleTypes.end(),
           std::ostream_iterator<char>(cerr, ""));
    }
    return 1;
  }
  
  if (numEvts < 1. && type_num != WIMP && type_num != B8 && type_num != ppSolar && type_num != atmNu) {
    if (verbosity > 0) cerr << "ERROR, you must simulate at least one event, or non-zero kg*days" << endl;
    return 1;
  }
  
  double maxTimeSep = DBL_MAX;
  if (type_num == Kr83m) {
    if (((eMin > 9.35 && eMin < 9.45) || (eMin >= 32. && eMin < 32.2) ||
         (eMin > 41. && eMin <= 41.6)) &&
        eMin != eMax) {
      maxTimeSep = eMax;
      if (eMax <= 0.) {
        if (verbosity > 0) cerr << "Max t sep must be +." << endl;
        return 1;
      }
    } else {
      if (verbosity > 0)
        cerr << "ERROR: For Kr83m, put E_min as 9.4, 32.1, or 41.5 keV and "
                "E_max as the max time-separation [ns] between the two decays "
                "please."
             << endl;
      return 1;
    }
  }

  if ((eMin < 10. || eMax < 10.) && type_num == gammaRay && verbosity > 0) {
    cerr << "WARNING: Typically beta model works better for ER BG at low "
            "energies as in a WS."
         << endl;
    cerr << "ER data is often best matched by a weighted average of the beta & "
            "gamma models."
         << endl;
  }

  double rho = n.SetDensity(detector->get_T_Kelvin(), detector->get_p_bar());
  if (rho <= 0. || detector->get_T_Kelvin() <= 0. ||
      detector->get_p_bar() <= 0.) {
    if (verbosity > 0) cerr << "ERR: Unphysical thermodynamic property!";
    return 1;
  }
  if (rho < 1.75 && ValidityTests::nearlyEqual(ATOM_NUM, 54.))
    detector->set_inGas(true);

  double Wq_eV = NESTcalc::WorkFunction(rho, detector->get_molarMass(),
                                        detector->get_OldW13eV())
                     .Wq_eV;
  // if ( rho > 3. ) detector->set_OldW13eV(false); //solid OR enriched. Units
  // of g/mL
  if (!detector->get_OldW13eV())
    Wq_eV = 11.5;  // 11.5±0.5(syst.)±0.1(stat.) from EXO
  if (loopNEST && dEOdxBasis) Wq_eV = -1000. / ERWeightParam[0];

  // Calculate and print g1, g2 parameters (once per detector)
  vector<double> g2_params = n.CalculateG2(verbosity);
  g2 = std::abs(g2_params[3]);
  double g1 = detector->get_g1();

  double centralZ =
      (detector->get_gate() * 0.8 + detector->get_cathode() * 1.03) /
      2.;  // fid vol def usually shave more off the top, because of gas
  // interactions (100->10cm)
  double centralField = detector->FitEF(0.0, 0.0, centralZ);
  if (inField != -1.) centralField = inField;

  if (type_num == WIMP) {
    yieldsMax =
        n.GetYields(NR, 25.0, rho, centralField, detector->get_molarMass(),
                    double(atomNum), NRYieldsParam, ERYieldsParam);
  } else if (type_num == B8) {
    yieldsMax =
        n.GetYields(NR, 4.00, rho, centralField, detector->get_molarMass(),
                    double(atomNum), NRYieldsParam, ERYieldsParam);
  } else {
    double energyMaximum;
    if (eMax < 0.)
      energyMaximum = 1. / std::abs(eMax);
    else
      energyMaximum = eMax;
    if (!energyMaximum || std::isnan(energyMaximum)) energyMaximum = 1.;
    if (type_num == Kr83m)
      yieldsMax = n.GetYields(Kr83m, eMin, rho, centralField, 400., 100.,
                              NRYieldsParam, ERYieldsParam);
    else
      yieldsMax = n.GetYields(type_num, energyMaximum, rho, centralField,
                              double(massNum), double(atomNum), NRYieldsParam,
                              ERYieldsParam);
  }
  if ((g1 * yieldsMax.PhotonYield) > (2. * maxS1) && eMin != eMax &&
      type_num != Kr83m && verbosity > 0 && !dEOdxBasis)
    cerr
        << "\nWARNING: Your energy maximum may be too high given your maxS1.\n";

  if (type_num < 6) massNum = 0;
  if (type_num == Kr83m) massNum = maxTimeSep;
  // use massNum to input maxTimeSep into GetYields(...)
  double keV = -999.;
  double timeStamp = dayNumber;
  vector<double> keV_vec;
  for (uint64_t j = 0; j < (uint64_t)numEvts; ++j) {
    // timeStamp = tMax * RandomGen::rndm()->rand_uniform(); /*OR: += tStep*/
    // detector->set_eLife_us(5e1+1e3*(timeStamp/3e2)); for E-recon when you've
    // changed g1,g2-related stuff, redo g1 and g2 calculations in line 477+
    try {
      if ((ValidityTests::nearlyEqual(eMin, eMax) && eMin >= 0. && eMax > 0.) ||
          type_num == Kr83m) {
        keV = eMin;
      } else {
        switch (type_num) {
          case CH3T:
            keV = TestSpectra::CH3T_spectrum(eMin, eMax);
            break;
          case C14:
            keV = TestSpectra::C14_spectrum(eMin, eMax);
            break;
          case B8:  // normalize this to ~3500 / 10-ton / year, for E-threshold
                    // of
            // 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
            keV = TestSpectra::B8_spectrum(eMin, eMax);
            break;
          case AmBe:
            //keV = TestSpectra::AmBe_spectrum(eMin, eMax); //for ZEPLIN-III FSR from HA (Pal '98)
	    keV = TestSpectra::PowLawFit_spectrum(eMin, eMax, -0.95351); //0.5,300 recommended Es
            break;
          case Cf:
            keV = TestSpectra::Cf_spectrum(eMin, eMax);
            break;
          case DD:
            //keV = TestSpectra::DD_spectrum(eMin, eMax, 10., 0.1, 60., 25., 0.); //LUX Run03
            keV = TestSpectra::DD_spectrum(eMin, eMax, 13., 0.12, 71.2, 20., -20.5); //LZ SR1 Ryan McM
            break;
          case WIMP:
            keV = TestSpectra::WIMP_spectrum(spec.wimp_spectrum_prep, eMin,
                                             timeStamp);
            break;
          case ppSolar:
            keV = TestSpectra::ppSolar_spectrum(eMin, eMax);
            break;
          case atmNu:
            keV = TestSpectra::atmNu_spectrum(eMin, eMax);
            break;
          case fullGamma:
            keV_vec = TestSpectra::Gamma_spectrum(eMin, eMax, gamma_source);
            if (keV_vec[0] > -1) {
              type_num = fullGamma_PE;  // here PE means photo-electric (effect)
                                        // and not photo-electrons :-)
              keV = keV_vec[0];
            } else {
              type_num = fullGamma_Compton_PP;  // PP = pair production
              if (keV_vec[1] > -1) {
                keV = keV_vec[1];
              } else {
                keV = keV_vec[2];
              }
            }
            break;
          case fullGamma_PE:
            keV_vec = TestSpectra::Gamma_spectrum(eMin, eMax, gamma_source);
            if (keV_vec[0] > -1) {
              type_num = fullGamma_PE;
              keV = keV_vec[0];
            } else {
              type_num = fullGamma_Compton_PP;
              if (keV_vec[1] > -1) {
                keV = keV_vec[1];
              } else {
                keV = keV_vec[2];
              }
            }
            break;
          case fullGamma_Compton_PP:
            keV_vec = TestSpectra::Gamma_spectrum(eMin, eMax, gamma_source);
            if (keV_vec[0] > -1) {
              type_num = fullGamma_PE;
              keV = keV_vec[0];
            } else {
              type_num = fullGamma_Compton_PP;
              if (keV_vec[1] > -1) {
                keV = keV_vec[1];
              } else {
                keV = keV_vec[2];
              }
            }
            break;
          default:
            // keV = TestSpectra::ZeplinBackground(); //example of continuous
            // ER/beta/gamma BG spec break;
            if (eMin < 0.) {
              if (!dEOdxBasis)
                return 1;
              else
                break;
            }
            if (eMax > 0.) {
              if (eMax > eMin)
                keV = eMin + (eMax - eMin) * RandomGen::rndm()->rand_uniform();
              else  // polymorphic feature: Gauss E distribution
                keV = RandomGen::rndm()->rand_zero_trunc_gauss(eMin, eMax);
            } else {  // negative eMax signals to NEST that you want to use an
              // exponential energy spectrum profile
              if (ValidityTests::nearlyEqual(eMin, 0.)) return 1;
              // eMin will be used in place of eMax as the maximum
              // energy in exponential scenario
              keV = RandomGen::rndm()->rand_exponential(-eMax * log(2.), 0.,
                                                        eMin);
            }
            break;
        }
      }

      if (type_num != WIMP && type_num != B8 && type_num != ppSolar &&
          type_num != atmNu && eMax > 0. && eMin < eMax) {
        if (keV > eMax) keV = eMax;
        if (keV < eMin) keV = eMin;
      }

      if (!dEOdxBasis) {
        if (keV < 0) {
          if (verbosity > 0)
            cerr << "ERROR: Get ready for time travel or FTL--negative energy!"
                 << endl;
          return 1;
        }
        if (ValidityTests::nearlyEqual(keV, 0.) && verbosity > 0) {
          cerr << "WARNING: Zero energy has occurred, and this is not normal"
               << endl;
        }
      }

      double FudgeFactor[2] = {1., 1.};
    Z_NEW:
      if (ValidityTests::nearlyEqual(
              fPos, -1.)) {  // -1 means default, random location mode
        pos_z = 0. + (detector->get_TopDrift() - 0.) *
                         RandomGen::rndm()->rand_uniform();  // initial guess
        r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
        phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
        pos_x = r * cos(phi);
        pos_y = r * sin(phi);
      } else {
        delimiter = ",";
        loc = 0;
        int i = 0;
        while ((loc = position.find(delimiter)) != string::npos) {
          token = position.substr(0, loc);
          if (i == 0)
            pos_x = stof(token);
          else
            pos_y = stof(token);
          position.erase(0, loc + delimiter.length());
          ++i;
        }
        if (!dEOdxBasis) {
          if (pos_x != -999. && pos_y != -999. &&
              (pos_x * pos_x + pos_y * pos_y) >
                  detector->get_radius() * detector->get_radius() &&
              j == 0 && verbosity > 0)
            cerr << "WARNING: outside fiducial radius." << endl;
          if (pos_x != -999. && pos_y != -999. &&
              (pos_x * pos_x + pos_y * pos_y) >
                  detector->get_radmax() * detector->get_radmax()) {
            if (verbosity > 0)
              cerr << "\nERROR: outside physical radius!!!" << endl;
            return EXIT_FAILURE;
          }
          pos_z = stof(position);
        } else {
          pos_z = stof(position);
        }
        if (ValidityTests::nearlyEqual(stof(position), -1.))
          pos_z = 0. + (detector->get_TopDrift() - 0.) *
                           RandomGen::rndm()->rand_uniform();
        if (ValidityTests::nearlyEqual(stof(token), -999.)) {
          r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
          phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
          pos_x = r * cos(phi);
          pos_y = r * sin(phi);
        }
        // if ( j == 0 ) { origX = pos_x; origY = pos_y; }
      }

      if (ValidityTests::nearlyEqual(
              inField, -1.)) {  // -1 means use poly position dependence
        field = detector->FitEF(pos_x, pos_y, pos_z);
      } else
        field = inField;  // no fringing

      if (field < 0. || detector->get_E_gas() < 0.) {
        if (verbosity > 0)
          cerr << "\nERROR: Neg field is not permitted. We don't simulate "
                  "field dir (yet). Put in magnitude.\n";
        return 1;
      }
      if ((ValidityTests::nearlyEqual(field, 0.) || std::isnan(field)) &&
          verbosity > 0)
        cerr
            << "\nWARNING: A LITERAL ZERO (or undefined) FIELD MAY YIELD WEIRD "
               "RESULTS. USE A SMALL VALUE INSTEAD.\n";
      if ((field > 12e3 || detector->get_E_gas() > 17e3) && verbosity > 0)
        cerr << "\nWARNING: Your field is >12,000 V/cm. No data out here. Are "
                "you sure about this?\n";

      if (j == 0 && ValidityTests::nearlyEqual(vD_middle, 0.)) {
        if (ValidityTests::nearlyEqual(inField, -1.)) {
          // build a vD table for non-uniform field, but if field varies in XY
          // not just Z you need to do more coding
          vTable = n.SetDriftVelocity_NonUniform(rho, z_step, detector->get_T_Kelvin(), detector->get_p_bar(), pos_x, pos_y);
          vD_middle = vTable[int(floor(centralZ / z_step + 0.5))];
          // for ( int jj = 0; jj < vTable.size(); ++jj ) //DEBUG
          // cerr << double(jj)*z_step << "\t" << vTable[jj] << endl;
        } else {
          vD_middle =
	    n.SetDriftVelocity(detector->get_T_Kelvin(), rho, inField, detector->get_p_bar());
          vD = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, field, detector->get_p_bar());
        }
        if (verbosity > 0) {
          cout << "Density = " << rho << " g/mL"
               << "\t";
          cout << "central vDrift = " << vD_middle << " mm/us\n";
          cout
              << "\t\t\tPRO TIP: neg s2_thr in detector->S2bot!\tW = " << Wq_eV
              << " eV\tNegative numbers are flagging things below threshold!   "
                 "phe=(1+P_dphe)*phd & phd=phe/(1+P_dphe)\n";

          if (type_num == Kr83m && (eMin < 10. || eMin > 40.))
            fprintf(stdout, "t [ns]\t\t");
          if (timeStamp > tZero) fprintf(stdout, "dayNum\t\t");
          if ((ValidityTests::nearlyEqual(eMax, eMin) || type_num == Kr83m) &&
              numBins == 1 && MCtruthE) {
            MCtruthE = false;
            if (verbosity > 0)
              fprintf(stderr,
                      "Simulating a mono-E peak; setting MCtruthE false.\n");
          }
          if (eMax > hiEregime)
            fprintf(stdout, "Energy [keV]");
          else {
            if (MCtruthE)
              fprintf(stdout, "E_truth [keV]");
            else
              fprintf(stdout, "E_recon [keV]");
          }
          if (seed == -2) {
            if (dEOdxBasis) {
              printf(
                  "\tfield [V/cm]\ttDrift [us]\tX,Y,Z "
                  "[mm]\tNph/keV\t\tNe-/keV\t\tS1 [PE "
                  "or phe]\tS1_3Dcor "
                  "[phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea "
                  "[PE]\tS2_3Dcorr [phd]\n");
	    }
            else {
              printf(
                  "\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\t\tNe-\t\tS1 "
                  "[PE "
                  "or phe]\tS1_3Dcor "
                  "[phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea "
                  "[PE]\tS2_3Dcorr [phd]\n");
	    }
          } else {
	    if ( verbosity >= 2 )
	      printf(
		       "\tfield [V/cm]\ttDrift [us]\tZ [mm]\tNph\tNe-\tS1 [PE or "
		     "phe]\tS1_3Dcor [phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea "
		     "[PE]\tS2_3Dcorr [phd]\n");
	    else
	      printf(
		   "\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\tNe-\tS1 [PE or "
		     "phe]\tS1_3Dcor [phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea "
		     "[PE]\tS2_3Dcorr [phd]\n");
	  }
        }
      }
      if (ValidityTests::nearlyEqual(inField, -1.)) {
        index = int(floor(pos_z / z_step + 0.5));
        vD = vTable[index];
      }
      driftTime = (detector->get_TopDrift() - pos_z) /
                  vD;  // (mm - mm) / (mm / us) = us
      if (inField != -1. &&
          detector->get_dt_min() > (detector->get_TopDrift() - 0.) / vD &&
          field >= FIELD_MIN) {
        if (verbosity > 0)
          cerr << "ERROR: dt_min is too restrictive (too large)" << endl;
        return 1;
      }
      if ((driftTime > detector->get_dt_max() ||
           driftTime < detector->get_dt_min()) &&
          (ValidityTests::nearlyEqual(fPos, -1.) ||
           ValidityTests::nearlyEqual(stof(position), -1.)) &&
          field >= FIELD_MIN)
        goto Z_NEW;
      if (detector->get_dt_max() > (detector->get_TopDrift() - 0.) / vD && !j &&
          field >= FIELD_MIN && verbosity > 0) {
        cerr << "WARNING: dt_max is greater than max possible" << endl;
      }

      // The following should never happen: this is simply a just-in-case
      // code-block dealing with user error (rounding another possible culprit)
      if (!dEOdxBasis) {
        if (pos_z <= 0) {
          if (verbosity > 0)
            cerr << "WARNING: unphysically low Z coordinate (vertical axis of "
                    "detector) of "
                 << pos_z << " mm" << endl;  // warn user on screen
          pos_z = z_step;
        }
        if ((pos_z > (detector->get_TopDrift() + z_step) || driftTime < 0.0) &&
            field >= FIELD_MIN) {
          if (verbosity > 0)
            cerr << "WARNING: unphysically big Z coordinate (vertical axis of "
                    "detector) of "
                 << pos_z << " mm" << endl;  // give the specifics
          driftTime = 0.0;
          pos_z = detector->get_TopDrift() - z_step;  // just fix it and move on
        }
      }

      YieldResult yields;
      QuantaResult quanta;
      NESTresult result;
      if (dEOdxBasis) {
        NRYieldsParam[0] = pos_x;
        NRYieldsParam[1] = pos_y;
        NRYieldsParam[2] = pos_z;
        NRYieldsParam[7] = vD_middle;
        NRYieldsParam[8] = rho;
        if (loopNEST && eMin == -1.)
          NRYieldsParam[4] = -TestSpectra::CH3T_spectrum(
              0., 18.6);  // replace with different rad calib source continuous
                          // in energy, as desired
        if (j == 0 && !loopNEST) {
          NRYieldsParam[3] = eMax;
          NRYieldsParam[4] = eMin;
          NRYieldsParam[5] =
              z_step;  // 5mm for fast, 6um for accurate; .01mm LUXRun3 WS
          NRYieldsParam[6] = inField;
          NRYieldsParam[9] = 1e3;  // 1MeV
          NRYieldsParam[10] =
              0.;  //+/-1=true, means use continuous slowing-down approx
          NRYieldsParam[11] =
              0.;  // use w/[10] for dE/dx = [10]*keV^[11] e.g. 50 & -0.5
          NRYieldsParam[12] = 0.;  // fractional variation: e.g .15 = 15%
        }
        result = n.GetYieldERdEOdxBasis(NRYieldsParam, posiMuon, vTable,
                                        ERWeightParam);
        yields = result.yields;
        quanta = result.quanta;
        if (ERWeightParam[0] >= 0.)
          driftTime = 0.00;
        else
          driftTime = (detector->get_TopDrift() - pos_z) / vD_middle;
        field = yields.ElectricField;
        pos_z = yields.DeltaT_Scint;
      } else {
        if (keV > 0.001 * Wq_eV) {
          if (type == "ER") {
            if (verbosity > 0 && j == 0) {
              cerr << "CAUTION: Are you sure you don't want beta model instead "
                      "of ER? This is a recomb-enhanced (Qy-reduced) version of"
                      " it for non-betas e.g. EC, DEC. Seed is now Q reduction."
		   << endl <<
		      ">1 means increase, and negative number means gamma model"
		      "; use 0 for a weighted average of the gamma, beta models"
                   << endl;
            }
	    if ( multFact > 0. )
	      yields = n.GetYieldBetaGR(keV, rho, field, ERYieldsParam, multFact);
	    else if ( multFact < 0. )
	      yields = n.GetYieldGamma(keV, rho, field, -multFact);
	    else
	      yields = n.GetYieldERWeighted(keV, rho, field, ERYieldsParam,
			        default_EnergyParams, default_FieldParams);
          } else {
            if (seed < 0 && seed != -1 && type_num <= 5)
              massNum = detector->get_molarMass();
            if (type_num == 11) atomNum = minTimeSeparation;  // Kr83m
            yields = n.GetYields(type_num, keV, rho, field, double(massNum),
                                 double(atomNum), NRYieldsParam, ERYieldsParam);
          }
          if (type_num == ion) {  // alphas +other nuclei, lighter/heavier than
                                  // medium's default nucleus
            NRERWidthsParam.clear();
            NRERWidthsParam = {
                1.00, 1.00,     0.,    0.50, 0.19, 0.,
                1.,   0.046452, 0.205, 0.45, -0.2};  // zero out non-binom
                                                     // recomb fluct & skew (NR)
          }
          if (!dEOdxBasis)  // last argument -999 for default VV skew mod
            quanta = n.GetQuanta(yields, rho, NRERWidthsParam, 0.);
        } else {
          yields.PhotonYield = 0.;
          yields.ElectronYield = 0.;
          yields.ExcitonRatio = 0.;
          yields.Lindhard = 0.;
          yields.ElectricField = 0.;
          yields.DeltaT_Scint = 0.;
          quanta.photons = 0;
          quanta.electrons = 0;
          quanta.ions = 0;
          quanta.excitons = 0;
        }
      }

      if (detector->get_noiseBaseline()[2] != 0. ||
          detector->get_noiseBaseline()[3] != 0.)
        quanta.electrons +=
            int(floor(RandomGen::rndm()->rand_gauss(
                          detector->get_noiseBaseline()[2],
                          detector->get_noiseBaseline()[3], false) +
                      0.5));

      // If we want the smeared positions (non-MC truth), then implement
      // resolution function
      double truthPos[3] = {pos_x, pos_y, pos_z};
      double smearPos[3] = {pos_x, pos_y, pos_z};
      double Nphd_S2 =
          g2 * quanta.electrons * exp(-driftTime / detector->get_eLife_us());
      if (!MCtruthPos && Nphd_S2 > PHE_MIN) {
        vector<double> xySmeared(2);
        xySmeared = n.xyResolution(pos_x, pos_y, Nphd_S2);
        smearPos[0] = xySmeared[0];
        smearPos[1] = xySmeared[1];
      }

      vector<int64_t> wf_time;
      vector<double> wf_amp;
      vector<double> scint =
          n.GetS1(quanta, truthPos[0], truthPos[1], truthPos[2], smearPos[0],
                  smearPos[1], smearPos[2], vD, vD_middle, type_num, j, field,
                  keV, s1CalculationMode, verbosity, wf_time, wf_amp);
      if (truthPos[2] < detector->get_cathode() && !dEOdxBasis)
        quanta.electrons = 0;
      vector<double> scint2;

      scint2 = n.GetS2(quanta.electrons, truthPos[0], truthPos[1], truthPos[2],
                       smearPos[0], smearPos[1], smearPos[2], driftTime, vD, j,
                       field, s2CalculationMode, verbosity, wf_time, wf_amp,
                       g2_params);
      if (dEOdxBasis && !loopNEST) {
        driftTime = (detector->get_TopDrift() - pos_z) / vD_middle;
        if (scint2[7] != -PHE_MIN && NRERWidthsParam[0] >= 0.)
          scint2[7] *= exp(driftTime / detector->get_eLife_us());
      }

    NEW_RANGES:
      if (usePD == 0 && std::abs(scint[3]) > minS1 && scint[3] < maxS1)
        signal1.push_back(scint[3]);
      else if (usePD == 1 && std::abs(scint[5]) > minS1 && scint[5] < maxS1)
        signal1.push_back(scint[5]);
      else if ((usePD >= 2 && std::abs(scint[7]) > minS1 && scint[7] < maxS1) ||
               maxS1 >= 998.5)  // xtra | handles bizarre bug of ~0eff, S1=999
        signal1.push_back(scint[7]);
      else
        signal1.push_back(-999.);

      double lowest1 = (double)detector->get_coinLevel();
      double lowest2 = minS2;
      if (lowest1 <= 0.) lowest1 = 1.;
      if (lowest2 <= 0.) lowest2 = 1.;
      if ((useS2 == 0 && logMax <= log10(maxS2 / maxS1)) ||
          (useS2 == 1 && logMax <= log10(maxS2)) ||
          (useS2 == 2 && logMax <= log10(maxS2 / maxS1))) {
        if (j == 0 && verbosity > 0)
          cerr << "err: You may be chopping off the upper half of your (ER?) "
                  "band; increase logMax and/or maxS2"
               << endl;
      }
      if ((useS2 == 0 && logMin >= log10(lowest2 / lowest1)) ||
          (useS2 == 1 && logMin >= log10(lowest2)) ||
          (useS2 == 2 && logMin >= log10(lowest2 / lowest1))) {
        if (j == 0 && verbosity > 0)
          cerr << "err: You may be chopping off the lower half of your (NR?) "
                  "band; decrease logMin and/or minS2"
               << endl;
      }

      if (usePD == 0 && std::abs(scint2[5]) > minS2 && scint2[5] < maxS2)
        signal2.push_back(scint2[5]);
      else if (usePD >= 1 && std::abs(scint2[7]) > minS2 && scint2[7] < maxS2)
        signal2.push_back(scint2[7]);  // no spike option for S2
      else
        signal2.push_back(-999.);

      if (ValidityTests::nearlyEqual(eMin, eMax) || type_num == Kr83m) {
        if ((scint[3] > maxS1 || scint[5] > maxS1 || scint[7] > maxS1) &&
            j < 10 && verbosity > 0)
          cerr << "WARNING: Some S1 pulse areas are greater than maxS1" << endl;
        if ((scint2[5] > maxS2 || scint2[7] > maxS2) && j < 10 &&
            verbosity > 0)  // don't repeat too much: only if within first 10
                            // events then show (+above)
          cerr << "WARNING: Some S2 pulse areas are greater than maxS2" << endl;
      }

      double Nph = 0.0, Ne = 0.0;
      if (!MCtruthE) {
        double MultFact = 1., eff = detector->get_sPEeff();

        if (s1CalculationMode == S1CalculationMode::Full ||
            s1CalculationMode == S1CalculationMode::Hybrid ||
            s1CalculationMode == S1CalculationMode::Waveform) {
          if (detector->get_sPEthr() >= 0. && detector->get_sPEres() > 0. &&
              eff > 0.) {
            MultFact = 0.5 *
                       (1. + erf((detector->get_sPEthr() - 1.) /
                                 (detector->get_sPEres() * sqrt(2.)))) /
                       (detector->get_sPEres() * sqrt(2. * M_PI));
            if (usePD == 2 && scint[7] < SPIKES_MAXM)
              MultFact =
                  1.04;  // see p. 5 Phys. Rev. D 97, 112002 (2018). Everything
                         // else here is just made up (empirical)
            else
              MultFact =
                  1. +
                  MultFact *
                      detector->get_P_dphe();  // correction factor discovered
                                               // by UA student Emily Mangus
            if (eff < 1.)
              eff += ((1. - eff) / (2. * double(detector->get_numPMTs()))) *
                     scint[0];
            eff = max(0., min(eff, 1.));
            MultFact /= eff;  // for smaller detectors leave it as 1.00
          } else
            MultFact = 1.;
        }
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
        Nph *= MultFact;
        if (signal1.back() <= 0. && timeStamp == tZero) Nph = 0.;
        if (signal2.back() <= 0. && timeStamp == tZero) Ne = 0.;
        if (detector->get_coinLevel() <= 0 && Nph <= PHE_MIN) Nph = DBL_MIN;
        if (yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0.) {
          // if(!detector->get_OldW13eV()) yields.Lindhard = 1.;
          if (ValidityTests::nearlyEqual(yields.Lindhard, 1.))
            keV = (Nph / FudgeFactor[0] + Ne / FudgeFactor[1]) * Wq_eV * 1e-3;
          else {
            if (type_num <= NEST::INTERACTION_TYPE::Cf) {
              keV = pow((Ne + Nph) / NRYieldsParam[0], 1. / NRYieldsParam[1]);
              Ne *= 1. - 1. / pow(1. + pow((keV / NRYieldsParam[5]),
                                           NRYieldsParam[6]),
                                  NRYieldsParam[10]);
              Nph *= 1. - 1. / pow(1. + pow((keV / NRYieldsParam[7]),
                                            NRYieldsParam[8]),
                                   NRYieldsParam[11]);
              keV = pow((Ne + Nph) / NRYieldsParam[0], 1. / NRYieldsParam[1]);
            } else {
              const double logden = log10(rho);
              const double Wq_eV_alpha =
                  28.259 + 25.667 * logden - 33.611 * pow(logden, 2.) -
                  123.73 * pow(logden, 3.) - 136.47 * pow(logden, 4.) -
                  74.194 * pow(logden, 5.) - 20.276 * pow(logden, 6.) -
                  2.2352 * pow(logden, 7.);
              keV = (Nph + Ne) * Wq_eV_alpha * 1e-3 / yields.Lindhard;  // cheat
            }  // nuclear recoil other than "NR" (Xe-Xe)
          }
          keVee +=
              (Nph + Ne) * Wq_eV * 1e-3;  // as alternative, use W_DEFAULT in
          // both places, but will not account for density dependence
        } else
          keV = 0.;
      }
      if ((signal1.back() <= 0. || signal2.back() <= 0.) &&
          field >= FIELD_MIN && detector->get_coinLevel() != 0 &&
          timeStamp == tZero)
        signalE.push_back(0.);
      else
        signalE.push_back(keV);

      if (type == "MIP" && eMin >= 1. && eMin <= 2. && MCtruthE) {
        keV = -1.;
      }

      if ((ValidityTests::nearlyEqual(Ne + Nph, 0.00) || std::isnan(keVee)) &&
          eMin == eMax && eMin > hiEregime &&
          !BeenHere) {  // efficiency is zero?
        minS1 = -999.;
        minS2 = -999.;
        detector->set_s2_thr(0.0);  // since needs GetS2() re-run, only
                                    // "catches" after the first caught event
        if (maxS2 > 1e10) BeenHere = true;
        maxS1 = 1e9;
        maxS2 = 1e11;
        numBins = 1;
        MCtruthE = false;
        signal1.pop_back();
        signal2.pop_back();
        signalE.pop_back();
        if (verbosity > 0) {
          cerr << endl
               << "CAUTION: Efficiency seems to have been zero, so trying "
                  "again with full S1 and S2 ranges."
               << endl;
          cerr << "OR, you tried to simulate a mono-energetic peak with MC "
                  "truth E turned on. Silly! Setting MCtruthE to false."
               << endl;
        }
        goto NEW_RANGES;
      }

      // Possible outputs from "scint" vector
      // scint[0] = nHits; // MC-true integer hits in same OR different PMTs, NO
      // double phe effect
      // scint[1] = Nphe; // MC-true integer hits WITH double phe effect (Nphe >
      // nHits)
      // scint[2] = pulseArea; // floating real# smeared DAQ pulse areas in phe,
      // NO XYZ correction
      // scint[3] = pulseAreaC; // smeared DAQ pulse areas in phe, WITH XYZ
      // correction
      // scint[4] = Nphd; // same as pulse area, adjusted/corrected *downward*
      // for 2-PE effect (LUX phd units) scint[5] = NphdC; // same as Nphd, but
      // XYZ-corrected scint[6] = spike; // floating real# spike count, NO XYZ
      // correction scint[7] = spikeC; // floating real# spike count, WITH XYZ
      // correction scint[8] = nHits post coincidence window and single PE eff!

      // Possible outputs from "scint2" vector
      // scint2[0] = Nee; // integer number of electrons unabsorbed in liquid
      // then getting extracted scint2[1] = Nph; // raw number of photons
      // produced in the gas gap scint2[2] = nHits; // MC-true integer hits in
      // same OR different PMTs, NO double phe effect scint2[3] = Nphe; //
      // MC-true integer hits WITH double phe effect (Nphe > nHits). S2 has more
      // steps than S1 (e's 1st)
      //
      // If S2 threshold is set to positive (normal mode)
      // scint2[4] = pulseArea; // floating real# smeared DAQ pulse areas in
      // phe, NO XYZ correction scint2[5] = pulseAreaC; // smeared DAQ pulse
      // areas in phe, WITH XYZ correction scint2[6] = Nphd; // same as pulse
      // area, adjusted/corrected *downward* for 2-PE effect (LUX phd units)
      // scint2[7] = NphdC; // same as Nphd, but XYZ-corrected
      //
      // If S2 threshold is set to negative (switches from S2 -> S2 bottom, NOT
      // literally negative)
      // scint2[4] = S2b; // floating real# smeared pulse areas in phe ONLY
      // including bottom PMTs, NO XYZ correction
      // scint2[5] = S2bc; // floating real# smeared pulse areas in phe ONLY
      // including bottom PMTs, WITH XYZ correction
      // scint2[6] = S2b / (1.+detector->get_P_dphe()); // same as S2b, but
      // adjusted for 2-PE effect (LUX phd units)
      // scint2[7] = S2bc / (1.+detector->get_P_dphe()); // same as S2bc, but
      // adjusted for 2-PE effect (LUX phd units)
      // scint2[8] = g2; // g2 = ExtEff * SE, light collection efficiency of EL
      // in gas gap (from CalculateG2)

      if (truthPos[2] < detector->get_cathode() && verbosity > 0 && !BeenHere &&
          !dEOdxBasis) {
        BeenHere = true;
        if (verbosity > 0)
          fprintf(
              stderr,
              "gamma-X i.e. MSSI may be happening. This may be why even high-E "
              "eff is <100%%. Check your cathode position definition.\n\n");
      }

      if (PrintSubThr ||
          (scint[0] > PHE_MIN && scint[1] > PHE_MIN && scint[2] > PHE_MIN &&
           scint[3] > PHE_MIN && scint[4] > PHE_MIN && scint[5] > PHE_MIN &&
           scint[6] > PHE_MIN && scint[7] > PHE_MIN && scint[8] > PHE_MIN &&
           scint2[0] > PHE_MIN && scint2[1] > PHE_MIN && scint2[2] > PHE_MIN &&
           scint2[3] > PHE_MIN && scint2[4] > PHE_MIN && scint2[5] > PHE_MIN &&
           scint2[6] > PHE_MIN && scint2[7] > PHE_MIN && scint2[8] > PHE_MIN)) {
        // for skipping specific sub-threshold events (save screen/disk space)
        // other suggestions: minS1, minS2 (or s2_thr) for tighter cuts
        // depending on analysis.hh settings (think of as analysis v. trigger
        // thresholds) and using max's too, pinching both ends
        if (type_num == Kr83m && verbosity > 0 && (eMin < 10. || eMin > 40.))
          printf("%.6f\t", yields.DeltaT_Scint);
        if (timeStamp > tZero && verbosity > 0) {
          printf("%.6f\t", timeStamp);
          if (keV > AnnModERange[0] && keV < AnnModERange[1])
            SaveTheDates[int(timeStamp) % tMax]++;
        }
        if (seed < 0 && seed != -1) {  // for when you want means
          if (dEOdxBasis) {
            if (eMin < 0.) {
              yields.PhotonYield = -quanta.photons / NRYieldsParam[4];
              yields.ElectronYield = -quanta.electrons / NRYieldsParam[4];
            } else {
              yields.PhotonYield /= NRYieldsParam[9];
              yields.ElectronYield /= NRYieldsParam[9];
            }
          }
          printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%lf\t%lf\t", keV, field,
                 driftTime, smearPos[0], smearPos[1], smearPos[2],
                 yields.PhotonYield, yields.ElectronYield);
        } else {
	  if ( verbosity >= 2 )
	    printf("%.6f\t%.6f\t%.6f\t%.0f\t%d\t%d\t", keV, field,
		   driftTime, smearPos[2],
		   quanta.photons, quanta.electrons);
	  else
	    printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%d\t%d\t", keV, field,
		             driftTime, smearPos[0], smearPos[1], smearPos[2],
		   quanta.photons, quanta.electrons);
	}
        if (keV > 10. * hiEregime || scint[5] > maxS1 || scint2[7] > maxS2 ||
            // switch to engineering notation to make output more readable, if
            // energy is too high (>1 MeV)
            (dEOdxBasis && eMin > 0.)) {
          printf("%e\t%e\t%e\t", scint[2], scint[5], scint[7]);
          printf("%lli\t%e\t%e\n", (long long int)scint2[0], scint2[4],
                 scint2[7]);
        } else {
          printf(
              "%.6f\t%.6f\t%.6f\t", scint[2], scint[5],
              scint[7]);  // see GetS1 inside of NEST.cpp for full explanation
          // of all 8 scint return vector elements. Sample 3
          // most common
          printf(
              "%i\t%.6f\t%.6f\n", (int)scint2[0], scint2[4],
              scint2[7]);  // see GetS2 inside of NEST.cpp for full explanation
          // of all 8 scint2 vector elements. Change as you
          // desire
        }
      }  // always execute statement, if(1) above, because if is just
         // place-holder
      // in case you want to drop all sub-threshold data
    } catch (exception& e) {
      if (verbosity > 0) cerr << e.what() << endl;
      return 1;
    }
  }  // end of the gigantic primary event number loop

  if (timeStamp > tZero) {
    for (unsigned int j = 0; j < tMax; ++j) {
      cerr << j << "\t" << SaveTheDates[j] << endl;
    }
  }

  if (verbosity >= 0) {
    if (eMin != eMax && type_num != Kr83m) {
      if (useS2 == 2)
        GetBand(signal2, signal1, false, detector->get_coinLevel());
      else
        GetBand(signal1, signal2, false, detector->get_coinLevel());
      fprintf(stderr,
              "Bin Center\tBin Actual\tHist Mean\tMean Error\tHist "
              "Sigma\tSkewness\t\tEff[%%>thr]\n");
      for (int j = 0; j < numBins; ++j) {
        fprintf(stderr, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\t%lf\n", band[j][0],
                band[j][1], band[j][2], band[j][4], band[j][3], band[j][6],
                band[j][5] * 100.);
        if (band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 ||
            band[j][3] <= 0.0 || band[j][4] <= 0.0 || band[j][5] <= 0.0 ||
            std::isnan(band[j][0]) || std::isnan(band[j][1]) ||
            std::isnan(band[j][2]) || std::isnan(band[j][3]) ||
            std::isnan(band[j][4]) || std::isnan(band[j][5])) {
          if (eMax != -999. && verbosity > 0) {
            if (((g1 * yieldsMax.PhotonYield) < maxS1 ||
                 (g2 * yieldsMax.ElectronYield) < maxS2) &&
                j != 0)
              cerr << "WARNING: Insufficient number of high-energy events to "
                      "populate highest bins is likely. Decrease maxS1 and/or "
                      "increase maxS2. Also: up the stats, change E range\n";
            else
              cerr << "WARNING: Insufficient number of low-energy events to "
                      "populate lowest bins is likely. Increase minS1 and/or "
                      "decrease minS2,s2_thr. Up the stats, change E range\n";
          }
          eMax = -999.;
        }
      }
    } else {
      GetBand(signal1, signal2, true, detector->get_coinLevel());
      GetEnergyRes(signalE);
      if (type_num == NR || type_num == WIMP || type_num == B8 ||
          type_num == DD || type_num == AmBe || type_num == Cf ||
          type_num == ion) {
        if (verbosity > 0)
          fprintf(stderr,
                  "S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc "
                  "[keVnr]\tEc Res[%%]\tEff[%%>thr]\tEc [keVee]\n");
        keVee /= numEvts;
      } else
        fprintf(stderr,
                "S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc Mean\t\tEc "
                "Res[%%]\tEff[%%>thr]\n");  // the C here refers to the combined
      // (S1+S2) energy scale
      for (int j = 0; j < numBins; ++j) {
        fprintf(stderr, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", band[j][0],
                band[j][1] / band[j][0] * 100., band[j][2],
                band[j][3] / band[j][2] * 100., energies[0],
                energies[1] / energies[0] * 100., energies[2] * 100.);
        if (type_num < 7)  // 0-6=NR/related (WIMPs,etc.)
          fprintf(stderr, "%lf\n", keVee / energies[2]);
        else
          fprintf(stderr, "\n");
        if ((band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 ||
             band[j][3] <= 0.0 || std::isnan(band[j][0]) ||
             std::isnan(band[j][1]) || std::isnan(band[j][2]) ||
             std::isnan(band[j][3])) &&
            field >= FIELD_MIN) {
          if (verbosity > 0) {
            if (numEvts > 1)
              cerr << "CAUTION: YOUR S1 and/or S2 MIN and/or MAX may be set to "
                      "be too restrictive, please check.\n";
            else
              cerr << "CAUTION: Poor stats. You must have at least 2 events to "
                      "calculate S1 and S2 and E resolutions.\n";
          }
        } else if ((ValidityTests::nearlyEqual(energies[0], eMin) ||
                    ValidityTests::nearlyEqual(energies[0], eMax) ||
                    energies[1] <= 1E-6) &&
                   field >= FIELD_MIN)
          if (verbosity > 0)
            cerr << "If your energy resolution is 0% then you probably still "
                    "have MC truth energy on."
                 << endl;
      }
    }
  }
  return 0;
}

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol, int nFold) {
  if (numBins > NUMBINS_MAX) {
    if (verbosity > 0)
      cerr << "ERROR: Too many bins. Decrease numBins (analysis.hh) or "
              "increase NUMBINS_MAX (TestSpectra.hh)"
           << endl;
    exit(EXIT_FAILURE);
  }

  vector<vector<double>> signals;
  signals.resize(numBins, vector<double>(1, -999.));
  double binWidth, border;
  if (useS2 == 2) {
    binWidth = (maxS2 - minS2) / double(numBins);
    border = minS2;
  } else {
    binWidth = (maxS1 - minS1) / double(numBins);
    border = minS1;
  }
  int i = 0, j = 0;
  double s1c, numPts;
  uint64_t reject[NUMBINS_MAX] = {0};

  if (resol) {
    numBins = 1;
    binWidth = DBL_MAX;
  }

  double ReqS1 = 0.;
  if (nFold == 0) ReqS1 = -DBL_MAX;
  for (i = 0; i < S1s.size(); ++i) {
    for (j = 0; j < numBins; ++j) {
      s1c = border + binWidth / 2. + double(j) * binWidth;
      if (i == 0 && !resol) band[j][0] = s1c;
      if ((std::abs(S1s[i]) > (s1c - binWidth / 2.) &&
           std::abs(S1s[i]) <= (s1c + binWidth / 2.)) ||
          !nFold) {
        if (S1s[i] >= ReqS1 && S2s[i] >= 0.) {
          if (resol) {
            signals[j].push_back(S2s[i]);
          } else {
            if (useS2 == 0) {
              if (S1s[i] && S2s[i] && log10(S2s[i] / S1s[i]) > logMin &&
                  log10(S2s[i] / S1s[i]) < logMax)
                signals[j].push_back(log10(S2s[i] / S1s[i]));
              else
                signals[j].push_back(0.);
            } else if (useS2 == 1) {
              if (S1s[i] && S2s[i] && log10(S2s[i]) > logMin &&
                  log10(S2s[i]) < logMax)
                signals[j].push_back(log10(S2s[i]));
              else
                signals[j].push_back(0.);
            } else {
              if (S1s[i] && S2s[i] && log10(S1s[i] / S2s[i]) > logMin &&
                  log10(S1s[i] / S2s[i]) < logMax)
                signals[j].push_back(log10(S1s[i] / S2s[i]));
              else
                signals[j].push_back(0.);
            }
          }
          band[j][2] += signals[j].back();
          if (resol)
            band[j][0] += S1s[i];
          else
            band[j][1] += S1s[i];
        } else
          ++reject[j];
        break;
      }
    }
  }

  for (j = 0; j < numBins; ++j) {
    if (band[j][0] <= 0. && !resol)
      band[j][0] = border + binWidth / 2. + double(j) * binWidth;
    signals[j].erase(signals[j].begin());
    numPts = (double)signals[j].size();
    if (numPts <= 0 && resol) {
      for (i = 0; i < S1s.size(); ++i) band[j][0] += std::abs(S1s[i]);
      numPts = S1s.size();
    }
    if (resol) band[j][0] /= numPts;
    band[j][1] /= numPts;
    band[j][2] /= numPts;
    if (numPts > signals[j].size())
      numPts = signals[j].size();  // seg fault prevention line
    for (i = 0; i < (int)numPts; ++i) {
      if (signals[j][i] != -999.)
        band[j][3] += pow(signals[j][i] - band[j][2], 2.);  // std dev calc
    }
    for (i = 0; i < S1s.size(); ++i) {
      if (resol && S1s[i] > ReqS1 && S2s[i] > 0.0)
        band[j][1] += pow(S1s[i] - band[j][0], 2.);  // std dev calc
    }
    band[j][3] /= numPts - 1.;
    band[j][3] = sqrt(band[j][3]);
    if (resol) {
      band[j][1] /= numPts - 1.;
      band[j][1] = sqrt(band[j][1]);
    }
    band[j][4] = band[j][3] / sqrt(numPts);
    band[j][5] = numPts / (numPts + double(reject[j]));
    for (i = 0; i < (int)numPts; ++i) {
      if (signals[j][i] != -999.)
        band[j][6] +=
            pow((signals[j][i] - band[j][2]) / band[j][3], 3.);  // skew calc
    }
    band[j][6] /= (numPts - 2.);
  }

  return signals;
}

void GetEnergyRes(vector<double> Es) {
  int i, numPts = Es.size();
  double numerator = 0.;

  for (i = 0; i < numPts; ++i) {
    if (Es[i] > 0.) {
      energies[0] += Es[i];
      ++numerator;
    }
  }

  energies[0] /= numerator;

  for (i = 0; i < numPts; ++i) {
    if (Es[i] > 0.) energies[1] += pow(energies[0] - Es[i], 2.);
  }

  energies[1] /= numerator - 1.;
  energies[1] = sqrt(energies[1]);

  energies[2] = numerator / double(numPts);
}
