/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   testNEST.cpp
 * Author: brodsky3
 *
 * Created on August 1, 2017, 1:03 PM
 */

#include "NEST.hh"
#include "TestSpectra.hh"
#include "analysis.hh"
#include "testNEST.hh"

#include "LUX_Run03.hh"

#define tZero 0.00 //day{of the year, 0 is average WIMP velocity}
#define tStep 0.03
#define hiEregime 1E+2 //keV

using namespace std;
using namespace NEST;

vector<double> FreeParam,  NuisParam;
double band[NUMBINS_MAX][7];
double energies[3];
bool BeenHere = false;

int main(int argc, char** argv) {
  // Instantiate your own VDetector class here, then load into NEST class
  // constructor
  DetectorExample_LUX_RUN03* detector = new DetectorExample_LUX_RUN03();

  // Custom parameter modification functions
  // detector->ExampleFunction();
  
  /* vector<double> eList = { 1., 2., 3. }; // fast example--for PLR, ML train
  vector<vector<double>> pos3dxyz = { {0.,-1.,60.},{-1.,0.,70.},{1.,0.,80.} };
  runNESTvec ( detector, NEST::beta, eList, pos3dxyz );
  return EXIT_SUCCESS; */
  
  if (argc < 7) {
    cout << "This program takes 6 (or 7) inputs, with Z position in mm from "
            "bottom of detector:" << endl;
    cout << "\t./testNEST numEvts type_interaction E_min[keV] E_max[keV] "
            "field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl
         << endl;
    cout << "For Kr83m time-dependent 9.4, 32.1, or 41.5 keV yields: " << endl;
    cout << "\t ./testNEST numEvents Kr83m Energy[keV] maxTimeDiff[ns] "
	    "field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl 
         << endl;
    cout << "For 8B, numEvts is kg-days of exposure with everything else same. "
            "For WIMPs:" << endl;
    cout << "\t./testNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] "
            "field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl
         << endl;
    cout << "For cosmic-ray muons or other similar particles with elongated "
            "track lengths:" << endl;
    cout << "\t./testNEST numEvts {MIP} LET[MeV*cm^2/gram] "
            "x,y-position[mm](Initial) field_drift[V/cm] "
            "x,y,z-position[mm](Final) {optional:seed}" << endl
         << endl;
    return 1;
  }
  
  unsigned long int numEvts;
  string type, position, posiMuon;
  double eMin, eMax, inField, fPos;
  int seed; bool no_seed = false;
  if ( loopNEST ) {
    
    numEvts = 100000; //10,000 faster but of course less precise
    if ( loopNEST == 1 )
      type = "ER";
    else
      type = "NR";
    eMin = 0.00;
    eMax = 156.; //Carbon-14
    inField = -1.;
    position = "-1";
    posiMuon = "-1";
    fPos = -1.;
    seed = 0;
    no_seed = false;
    FreeParam.clear(); NuisParam.clear();
    verbosity = false;
    useTiming = 0; //-1 for faster (less accurate)
    
    if ( type == "ER" ) {
      
      FreeParam.push_back(atof(argv[1])); //-0.1 for LUX C-14 ~200V/cm
      FreeParam.push_back(atof(argv[2])); //0.5
      FreeParam.push_back(atof(argv[3])); //0.06
      FreeParam.push_back(atof(argv[4])); //-0.6
      FreeParam.push_back(atof(argv[5])); //1.11
      FreeParam.push_back(atof(argv[6])); //0.95
      FreeParam.push_back(atof(argv[7])); //8e-2
      FreeParam.push_back(atof(argv[7])); //repeat
      
      NuisParam.push_back(11.);
      NuisParam.push_back(1.1);
      NuisParam.push_back(0.0480);
      NuisParam.push_back(-0.0533);
      NuisParam.push_back(12.6);
      NuisParam.push_back(0.3);
      NuisParam.push_back(2.);
      NuisParam.push_back(0.3);
      NuisParam.push_back(2.);
      NuisParam.push_back(0.5);
      NuisParam.push_back(1.0);
      NuisParam.push_back(1.0);
      
    }
    
    else {
      
      NuisParam.push_back(atof(argv[1])); //11.0 XENON10
      NuisParam.push_back(atof(argv[2])); //1.09
      NuisParam.push_back(0.0480);
      NuisParam.push_back(-0.0533);
      NuisParam.push_back(12.6);
      NuisParam.push_back(0.3);
      NuisParam.push_back(2.);
      NuisParam.push_back(0.3);
      NuisParam.push_back(atof(argv[3])); //2.00
      NuisParam.push_back(0.5);
      NuisParam.push_back(1.0);
      NuisParam.push_back(1.0);
      detector->set_g1(atof(argv[5])); //0.0725
      detector->set_g1_gas(atof(argv[6])); //0.0622
      detector->set_noiseL(atof(argv[7]),atof(argv[7])); //0,0
      FreeParam.push_back(1.00);
      FreeParam.push_back(1.00);
      FreeParam.push_back(atof(argv[4])); //0.070
      FreeParam.push_back(0.50);
      FreeParam.push_back(0.19);
      
    }
    
  }
  else {
    
    numEvts = (unsigned long)atof(argv[1]);
    if ( numEvts <= 0 ) {
      cerr << "ERROR, you must simulate at least 1 event" << endl;
      return 1;
    }
    type = argv[2];
    eMin = atof(argv[3]);
    eMax = atof(argv[4]);
    inField = atof(argv[5]);
    position = argv[6];
    posiMuon = argv[4];
    fPos = atof(argv[6]);
    
    seed = 0; //if not given make 0
    if (argc == 8) {
      seed = atoi(argv[7]);
    }
    else {
      RandomGen::rndm()->SetSeed(0);
      no_seed = true;
    }
    
    FreeParam.clear();
    NuisParam.clear();
    if ( type == "ER" ) {
      FreeParam.push_back(0.500);//LUX Run03
      FreeParam.push_back(0.5);
      FreeParam.push_back(1.1);
      FreeParam.push_back(-5.);
      FreeParam.push_back(1.01);
      FreeParam.push_back(0.95);
      FreeParam.push_back(1.4e-2);
      FreeParam.push_back(1.8e-2);
    }
    else {
      FreeParam.push_back(1.00); // Fi (Fano factor for ionization)
      FreeParam.push_back(1.00); // Fex
      FreeParam.push_back(0.10); // amplitude for non-binomial recombination fluctuations
      FreeParam.push_back(0.50); // center in e-Frac
      FreeParam.push_back(0.19); // width parameter (Gaussian 1-sigma)
    }
    NuisParam.push_back(11.); //alpha, for NR model. See http://nest.physics.ucdavis.edu
    NuisParam.push_back(1.1); //beta
    NuisParam.push_back(0.0480); //gamma
    NuisParam.push_back(-0.0533); //delta
    NuisParam.push_back(12.6); //epsilon
    NuisParam.push_back(0.3); //zeta
    NuisParam.push_back(2.); //eta
    NuisParam.push_back(0.3); //theta
    NuisParam.push_back(2.); //iota
    // last 3 are the secret extra parameters for additional flexibility
    NuisParam.push_back(0.5); //changes sqrt in Qy equation
    NuisParam.push_back(1.0); //makes low-E sigmoid an asymmetric one, for charge
    NuisParam.push_back(1.0); //makes low-E sigmoid an asymmetric one, for light
    
  }
  
  return testNEST(detector, numEvts, type, eMin, eMax, inField, position, posiMuon, fPos,
                  seed, no_seed, tZero);
}

NESTObservableArray runNESTvec ( VDetector* detector, INTERACTION_TYPE particleType, //func suggested by Xin Xiang, PD Brown U. for RG, LZ
				 vector<double> eList, vector<vector<double>> pos3dxyz, double inField, int seed ) {
  
  verbosity = false; NESTcalc n(detector);
  NESTresult result; QuantaResult quanta;
  double x, y, z, driftTime, vD; RandomGen::rndm()->SetSeed(seed);
  NuisParam = {11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1., 1.};
  FreeParam = {1.,1.,0.10,0.5,0.19};
  vector<double> scint, scint2, wf_amp; vector<long int> wf_time;
  NESTObservableArray OutputResults; double useField;
  vector<double> g2_params = n.CalculateG2(verbosity);
  double rho = n.SetDensity(detector->get_T_Kelvin(),detector->get_p_bar());
  
  for ( long i = 0; i < eList.size(); i++ ) {
    x = pos3dxyz[i][0];
    y = pos3dxyz[i][1];
    z = pos3dxyz[i][2];
    if ( inField <= 0. ) useField = detector->FitEF(x,y,z);
    else useField = inField;
    double truthPos[3] = { x, y, z };
    double smearPos[3] = { x, y, z }; //ignoring the difference in this quick function caused by smearing
    result = n.FullCalculation(particleType,eList[i],rho,useField,detector->get_molarMass(),ATOM_NUM,NuisParam,FreeParam,verbosity);
    quanta = result.quanta;
    vD = n.SetDriftVelocity(detector->get_T_Kelvin(),rho,useField);
    scint = n.GetS1(quanta,truthPos,smearPos,vD,vD,particleType,i,useField,eList[i],0,verbosity,wf_time,wf_amp); //0 means useTiming = 0
    driftTime = (detector->get_TopDrift()-z)/vD; //vD,vDmiddle assumed same (uniform field)
    scint2= n.GetS2(quanta.electrons,truthPos,smearPos,driftTime,vD,i,useField,0,verbosity,wf_time,wf_amp,g2_params);
    if ( scint[7] > PHE_MIN && scint2[7] > PHE_MIN ) { //unlike usual, kill (don't skip, just -> 0) sub-thr evts
      OutputResults.s1_nhits.push_back(abs(int(scint[0])));
      OutputResults.s1_nhits_thr.push_back(abs(int(scint[8])));
      OutputResults.s1_nhits_dpe.push_back(abs(int(scint[1])));
      OutputResults.s1r_phe.push_back(fabs(scint[2]));
      OutputResults.s1c_phe.push_back(fabs(scint[3]));
      OutputResults.s1r_phd.push_back(fabs(scint[4]));
      OutputResults.s1c_phd.push_back(fabs(scint[5]));
      OutputResults.s1r_spike.push_back(fabs(scint[6]));
      OutputResults.s1c_spike.push_back(fabs(scint[7])); //default is S1c in units of spikes, 3-D XYZ corr
      OutputResults.Nee.push_back(abs(int(scint2[0])));
      OutputResults.Nph.push_back(abs(int(scint2[1])));
      OutputResults.s2_nhits.push_back(abs(int(scint2[2])));
      OutputResults.s2_nhits_dpe.push_back(abs(int(scint2[3])));
      OutputResults.s2r_phe.push_back(fabs(scint2[4]));
      OutputResults.s2c_phe.push_back(fabs(scint2[5]));
      OutputResults.s2r_phd.push_back(fabs(scint2[6]));
      OutputResults.s2c_phd.push_back(fabs(scint2[7])); //default is S2c in terms of phd, not phe a.k.a. PE
    }
    else {
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

int testNEST(VDetector* detector, unsigned long int numEvts, string type,
             double eMin, double eMax, double inField, string position, string posiMuon,
             double fPos, int seed, bool no_seed, double dayNumber ) {
  // Construct NEST class using detector object
  NESTcalc n(detector);

  if (detector->get_TopDrift() <= 0. || detector->get_anode() <= 0. ||
      detector->get_gate() <= 0.) {
    cerr << "ERROR, unphysical value(s) of position within the detector "
            "geometry.";  // negative or 0 for cathode position is OK (e.g., LZ)
    return 1;
  }

  vector<double> signal1, signal2, signalE, vTable;
  string delimiter, token;
  size_t loc;
  int index;
  double g2, pos_x, pos_y, pos_z, r, phi, driftTime, field, vD,
    vD_middle = 0., atomNum = 0, keVee = 0.0;
  double massNum = detector->get_molarMass();
  YieldResult yieldsMax; //for warnings about S1 range

  if (no_seed != true) {
    if (seed == -1) {
      RandomGen::rndm()->SetSeed(time(NULL));
    }

    else {
      RandomGen::rndm()->SetSeed(seed);
    }
  }

  if (eMin == -1.) eMin = 0.;

  if (eMax == -1. && eMin == 0.)
    eMax = hiEregime;  // the default energy max
  if (eMax == 0.) {
    cerr << "ERROR: The maximum energy cannot be 0 keV!" << endl;
    return 1;
  }

  INTERACTION_TYPE type_num;
  TestSpectra spec;
  if (type == "NR" || type == "neutron" || type == "-1")
    type_num = NR;  //-1: default particle type is also NR
  else if (type == "WIMP") {
    if (eMin < 0.44) { //here eMin is WIMP mass
      cerr << "WIMP mass too low, you're crazy!" << endl;
      return 1;
    }
    type_num = WIMP;
    spec.wimp_spectrum_prep = spec.WIMP_prep_spectrum(eMin, E_step, dayNumber);
    numEvts = RandomGen::rndm()->poisson_draw(spec.wimp_spectrum_prep.integral * //here eMax is cross-section
                                              1.0 * double(numEvts) * eMax / 1e-36);
  } else if (type == "B8" || type == "Boron8" || type == "8Boron" ||
             type == "8B" || type == "Boron-8") {
    type_num = B8;
    numEvts = RandomGen::rndm()->poisson_draw(0.0026 * double(numEvts));
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
  else if (type == "C14" || type == "Carbon14" || type == "14C" || type == "C-14" || type == "Carbon-14")
    type_num = C14;
  else if (type == "beta" || type == "ER" || type == "Compton" ||
           type == "compton" || type == "electron" || type == "e-" ||
           type == "muon" || type == "MIP" || type == "LIP" || type == "mu" ||
           type == "mu-")
    type_num = NEST::beta;  // default electron recoil model
  else {
    cerr << "UNRECOGNIZED PARTICLE TYPE!! VALID OPTIONS ARE:" << endl;
    cerr << "NR or neutron," << endl;
    cerr << "WIMP," << endl;
    cerr << "B8 or Boron8 or 8Boron or 8B or Boron-8," << endl;
    cerr << "DD or D-D," << endl;
    cerr << "AmBe," << endl;
    cerr << "Cf or Cf252 or 252Cf or Cf-252," << endl;
    cerr << "ion or nucleus," << endl;
    cerr << "alpha," << endl;
    cerr << "gamma or gammaRay," << endl;
    cerr << "x-ray or xray or xRay or X-ray or Xray or XRay," << endl;
    cerr << "Kr83m or 83mKr or Kr83," << endl;
    cerr << "CH3T or tritium," << endl;
    cerr << "Carbon14 or 14C or C14 or C-14 or Carbon-14," << endl;
    cerr << "beta or ER or Compton or compton or electron or e-, and" << endl;
    cerr << "muon or MIP or LIP or mu or mu-" << endl;
    return 1;
  }
  
  double maxTimeSep = DBL_MAX;
  if (type_num == Kr83m) {
    if ( (eMin == 9.4 || eMin == 32.1 || eMin == 41.5) && eMin != eMax) {
      maxTimeSep = eMax;
    } else {
      cerr << "ERROR: For Kr83m, put E_min as 9.4, 32.1, or 41.5 keV "
              "and E_max as the max time-separation [ns] between the two decays "
	      "please." << endl;
      return 1;
    }
  }

  if ((eMin < 10. || eMax < 10.) && type_num == gammaRay) {
    cerr << "WARNING: Typically beta model works better for ER BG at low "
            "energies as in a WS." << endl;
    cerr << "ER data is often best matched by a weighted average of the beta & "
            "gamma models." << endl;
  }

  double rho = n.SetDensity(detector->get_T_Kelvin(), detector->get_p_bar());
  if (rho <= 0. || detector->get_T_Kelvin() <= 0. ||
      detector->get_p_bar() <= 0.) {
    cerr << "ERR: Unphysical thermodynamic property!";
    return 1;
  }
  if (rho < 1.75) detector->set_inGas(true);

  double Wq_eV = n.WorkFunction(rho).Wq_eV;
  //if ( rho > 3. ) detector->set_extraPhot(true); //solid OR enriched. Units of g/mL
  if ( detector->get_extraPhot() )
    Wq_eV = 11.5; //11.5±0.5(syst.)±0.1(stat.) from EXO
  
  // Calculate and print g1, g2 parameters (once per detector)
  vector<double> g2_params = n.CalculateG2(verbosity);
  g2 = fabs(g2_params[3]);
  double g1 = detector->get_g1();
  
  double centralZ =
      (detector->get_gate() * 0.8 + detector->get_cathode() * 1.03) /
      2.;  // fid vol def usually shave more off the top, because of gas
  // interactions (100->10cm)
  double centralField = detector->FitEF(0.0, 0.0, centralZ);
  if ( inField != -1. ) centralField = inField;

  if (type_num == WIMP) {
    yieldsMax = n.GetYields(NR, 25.0, rho, centralField, detector->get_molarMass(),
                            double(atomNum), NuisParam);
  } else if (type_num == B8) {
    yieldsMax = n.GetYields(NR, 4.00, rho, centralField, detector->get_molarMass(),
                            double(atomNum), NuisParam);
  } else {
    double energyMaximum;
    if (eMax < 0.)
      energyMaximum = 1. / fabs(eMax);
    else
      energyMaximum = eMax;
    if (type_num == Kr83m)
      yieldsMax = n.GetYields(NEST::beta, eMin, rho, centralField,
                              double(massNum), double(atomNum),
                              NuisParam);  // the reason for this: don't do the
    // special Kr stuff when just
    // checking max
    else
      yieldsMax = n.GetYields(type_num, energyMaximum, rho, centralField,
                              double(massNum), double(atomNum), NuisParam);
  }
  if ((g1 * yieldsMax.PhotonYield) > (2. * maxS1) && eMin != eMax)
    cerr
        << "\nWARNING: Your energy maximum may be too high given your maxS1.\n";
  
  if ( type_num < 6 ) massNum = 0;
  if ( type_num == Kr83m ) massNum = maxTimeSep; 
      //use massNum to input maxTimeSep into GetYields(...)
  double keV = -999.; double timeStamp = dayNumber;
  for (unsigned long int j = 0; j < numEvts; j++) {
    //timeStamp += tStep; //detector->set_eLife_us(5e1+1e3*(timeStamp/3e2));
    //for E-recon when you've changed g1,g2-related stuff, redo line 341+
    if ( (eMin == eMax && eMin >= 0. && eMax > 0.) || type_num == Kr83m ) {
      keV = eMin;
    } else {
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
          keV = spec.WIMP_spectrum(spec.wimp_spectrum_prep, eMin, timeStamp);
	  break;
        default:
          if (eMin < 0.) return 1;
          if (eMax > 0.)
            keV = eMin + (eMax - eMin) * RandomGen::rndm()->rand_uniform();
          else {  // negative eMax signals to NEST that you want to use an
            // exponential energy spectrum profile
            if (eMin == 0.) return 1;
            keV = 1e100;  // eMin will be used in place of eMax as the maximum
            // energy in exponential scenario
            while (keV > eMin)
              keV = eMax * log(RandomGen::rndm()->rand_uniform());
          }
          break;
      }
    }

    if (type_num != WIMP && type_num != B8 && eMax > 0.) {
      if (keV > eMax) keV = eMax;
      if (keV < eMin) keV = eMin;
    }
    
    if ( keV < 0. ) {
      cerr << "ERROR: Get ready for time travel or FTL--negative energy!" << endl;
      return 1;
    }
    if( keV == 0. ) {
      cerr << "WARNING: Zero energy has occurred, and this is not normal" << endl;
    }
    
    double FudgeFactor[2] = { 1., 1. };
  Z_NEW:
    if (fPos == -1.) {  // -1 means default, random location mode
      pos_z = 0. +
              (detector->get_TopDrift() - 0.) *
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
        i++;
      }
      if ( sqrt ( pos_x*pos_x+pos_y*pos_y ) > detector->get_radius() && j == 0 )
	cerr << "WARNING: outside fiducial radius." << endl;
      if ( sqrt ( pos_x*pos_x+pos_y*pos_y ) > detector->get_radmax() ) {
	cerr << "\nERROR: outside physical radius!!!" << endl; return EXIT_FAILURE; }
      pos_z = stof(position);
      if (stof(position) == -1.)
        pos_z =
            0. +
            (detector->get_TopDrift() - 0.) * RandomGen::rndm()->rand_uniform();
      if (stof(token) == -999.) {
        r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
        phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
        pos_x = r * cos(phi);
        pos_y = r * sin(phi);
      }
      // if ( j == 0 ) { origX = pos_x; origY = pos_y; }
    }

    if (inField == -1.) {  // -1 means use poly position dependence
      field = detector->FitEF(pos_x, pos_y, pos_z);
    } else
      field = inField;  // no fringing

    if (field < 0. || detector->get_E_gas() < 0.) {
      cerr << "\nERROR: Neg field is not permitted. We don't simulate field "
              "dir (yet). Put in magnitude.\n";
      return 1;
    }
    if (field == 0. || std::isnan(field))
      cerr << "\nWARNING: A LITERAL ZERO (or undefined) FIELD MAY YIELD WEIRD "
              "RESULTS. USE A SMALL VALUE INSTEAD.\n";
    if (field > 12e3 || detector->get_E_gas() > 17e3)
      cerr << "\nWARNING: Your field is >12,000 V/cm. No data out here. Are "
              "you sure about this?\n";

    if (j == 0 && vD_middle == 0.) {
      if (inField == -1.) {
        // build a vD table for non-uniform field, but if field varies in XY not
        // just Z you need to do more coding
        vTable = n.SetDriftVelocity_NonUniform(rho, z_step, pos_x, pos_y);
        vD_middle = vTable[int(floor(centralZ / z_step + 0.5))];
        // for ( int jj = 0; jj < vTable.size(); jj++ ) //DEBUG
        // cerr << double(jj)*z_step << "\t" << vTable[jj] << endl;
      } else {
        vD_middle = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, inField);
        vD = n.SetDriftVelocity(detector->get_T_Kelvin(), rho, field);
      }
      if (verbosity) {
        cout << "Density = " << rho << " g/mL"
             << "\t";
        cout << "central vDrift = " << vD_middle << " mm/us\n";
        cout << "\t\t\tPRO TIP: neg s2_thr in detector->S2bot!\tW = " << Wq_eV
             << " eV\tNegative numbers are flagging things below threshold!   "
                "phe=(1+P_dphe)*phd & phd=phe/(1+P_dphe)\n";

        if (type_num == Kr83m && eMin != 32.1)
          fprintf(stdout,
                  "t [ns]\t\t");
	if (type_num == WIMP && timeStamp > (tZero+tStep))
	  fprintf(stdout,
		  "dayNum\t");
	if ( eMax == eMin && numBins == 1 && MCtruthE == true ) {
	  MCtruthE = false;
	  fprintf(stderr,"Simulating a mono-E peak; setting MCtruthE false.\n");
	}
	if ( eMax > hiEregime )
	  fprintf(stdout,"Energy [keV]");
	else {
	  if ( MCtruthE ) fprintf(stdout,"E_truth [keV]");
	  else fprintf(stdout,"E_recon [keV]");
	}
	fprintf(stdout,
		"\tfield [V/cm]\ttDrift [us]\tX,Y,Z "
		"[mm]\tNph\tNe-\tS1 [PE or phe]\tS1_3Dcor "
		"[phd]\tspikeC(NON-INT)\tNe-Extr\tS2_rawArea [PE]\tS2_3Dcorr "
		"[phd]\n");
      }
    }
    if (inField == -1.) {
      index = int(floor(pos_z / z_step + 0.5));
      vD = vTable[index];
    }
    driftTime =
        (detector->get_TopDrift() - pos_z) / vD;  // (mm - mm) / (mm / us) = us
    if (inField != -1. &&
        detector->get_dt_min() > (detector->get_TopDrift() - 0.) / vD &&
        field >= FIELD_MIN) {
      cerr << "ERROR: dt_min is too restrictive (too large)" << endl;
      return 1;
    }
    if ((driftTime > detector->get_dt_max() ||
         driftTime < detector->get_dt_min()) &&
        (fPos == -1. || stof(position) == -1.) && field >= FIELD_MIN)
      goto Z_NEW;
    if (detector->get_dt_max() > (detector->get_TopDrift() - 0.) / vD && !j &&
        field >= FIELD_MIN) {
      cerr << "WARNING: dt_max is greater than max possible" << endl;
    }

    // The following should never happen: this is simply a just-in-case
    // code-block dealing with user error (rounding another possible culprit)
    if (pos_z <= 0.) {
      cerr << "WARNING: unphysically low Z coordinate (vertical axis of "
	"detector) of " << pos_z << " mm" << endl; //warn user on screen
      pos_z = z_step;
    }
    if ((pos_z > (detector->get_TopDrift() + z_step) || driftTime < 0.0) &&
        field >= FIELD_MIN) {
      cerr << "WARNING: unphysically big Z coordinate (vertical axis of "
	"detector) of " << pos_z << " mm" << endl; // give the specifics
      driftTime = 0.0;
      pos_z = detector->get_TopDrift() - z_step; //just fix it and move on
    }
    
    YieldResult yields;
    QuantaResult quanta;
    if (type == "muon" || type == "MIP" || type == "LIP" || type == "mu" ||
        type == "mu-") {
      double xi = -999., yi = -999.;
      if (eMax == -1.) {
        r = detector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
        phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
        xi = r * cos(phi);
        yi = r * sin(phi);
      } else {
        position = posiMuon;
        delimiter = ",";
        loc = 0;
        int ii = 0;
        while ((loc = position.find(delimiter)) != string::npos) {
          token = position.substr(0, loc);
          if (ii == 0)
            xi = stof(token);
          else
            yi = stof(token);
          position.erase(0, loc + delimiter.length());
          ii++;
        }
        yi = stof(position);
      }
      double dEOdx = eMin, eStep = dEOdx * rho * z_step * 1e2, refEnergy = 1e6;
      keV = 0.;
      int Nph = 0, Ne = 0;
      double xx = xi, yy = yi, zz = detector->get_TopDrift();
      double xf = pos_x;
      double yf = pos_y;
      double distance = sqrt(pow(xf - xi, 2.) + pow(yf - yi, 2.) +
                             pow(detector->get_TopDrift(), 2.));
      double norm[3];
      norm[0] = (xf - xi) / distance;
      norm[1] = (yf - yi) / distance;
      norm[2] = -detector->get_TopDrift() /
                distance;  // have not yet tested muons which leave before
      // hitting Z=0, would have to modify code here
      while (zz > 0. &&
             sqrt(pow(xx, 2.) + pow(yy, 2.)) <
                 detector->get_radmax()) {  // stop making S1 and S2 if particle
        // exits Xe vol
        yields =
            n.GetYields(NEST::beta, refEnergy, rho, detector->FitEF(xx, yy, zz),
                        double(massNum), double(atomNum), NuisParam);
        quanta = n.GetQuanta(yields, rho, FreeParam);
        Nph += quanta.photons * (eStep / refEnergy);
        index = int(floor(zz / z_step + 0.5));
        if (index >= vTable.size()) index = vTable.size() - 1;
        vD = vTable[index];
        driftTime = (detector->get_TopDrift() - zz) / vD;
        if (pos_z >= detector->get_cathode())
          Ne += quanta.electrons * (eStep / refEnergy) *
                exp(-driftTime / detector->get_eLife_us());
        keV += eStep;
        xx += norm[0] * z_step;
        yy += norm[1] * z_step;
        zz += norm[2] * z_step;  // cout << xx << " " << yy << " " << zz <<
                                 // endl;
      }
      quanta.photons = Nph;
      quanta.electrons = Ne;
      driftTime = 0.00;
      vD = vD_middle;  // approximate things not already done right in loop as
      // middle of detector since muon traverses whole length
      pos_x = .5 * (xi + xf);
      pos_y = .5 * (yi + yf);
      field = detector->FitEF(pos_x, pos_y, centralZ);
    } else {
      if (keV > .001 * Wq_eV) {
	if ( type == "ER" ) {
	  YieldResult yieldsB = n.GetYields(NEST::beta, keV, rho, field,
					    double(massNum), double(atomNum), NuisParam);
	  YieldResult yieldsG = n.GetYields(gammaRay, keV, rho, field,
					    double(massNum), double(atomNum), NuisParam);
	  double weightG = FreeParam[0] + FreeParam[1] * erf ( FreeParam[2] * ( log ( keV ) + FreeParam[3] ) ); // Xe10:1,.55,-1.6,-1.0
	  double weightB = 1. - weightG;
	  yields.PhotonYield = weightG * yieldsG.PhotonYield + weightB * yieldsB.PhotonYield;
	  yields.ElectronYield = weightG * yieldsG.ElectronYield + weightB * yieldsB.ElectronYield;
	  yields.ExcitonRatio = weightG * yieldsG.ExcitonRatio + weightB * yieldsB.ExcitonRatio;
	  yields.Lindhard = weightG * yieldsG.Lindhard + weightB * yieldsB.Lindhard;
	  yields.ElectricField = weightG * yieldsG.ElectricField + weightB * yieldsB.ElectricField;
	  yields.DeltaT_Scint = weightG * yieldsG.DeltaT_Scint + weightB * yieldsB.DeltaT_Scint;
	  FudgeFactor[0] = FreeParam[4];//0.99;
	  FudgeFactor[1] = FreeParam[5];//1.04;
	  yields.PhotonYield *= FudgeFactor[0];
	  yields.ElectronYield*=FudgeFactor[1];
	  detector->set_noiseL(FreeParam[6], FreeParam[7]); // XENON10: 1.0, 1.0. Hi-E gam: ~0-2%,6-5%
	}
	else {
	  yields = n.GetYields(type_num, keV, rho, field, double(massNum), double(atomNum), NuisParam);
	}
	//FreeParam.clear();
	//FreeParam = { 1.00, 1.00, 0.100, 0.50, 0.19 };
        quanta = n.GetQuanta(yields, rho, FreeParam);
      }
      else {
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
    
    if ( detector->get_noiseB()[2] != 0. || detector->get_noiseB()[3] != 0. )
      quanta.electrons += int(floor(RandomGen::rndm()->rand_gauss(
		   detector->get_noiseB()[2],detector->get_noiseB()[3])+0.5));
    
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
    
    vector<long int> wf_time;
    vector<double> wf_amp;
    vector<double> scint =
      n.GetS1(quanta, truthPos, smearPos, vD, vD_middle, type_num, j, field,
                keV, useTiming, verbosity, wf_time, wf_amp);
    if (truthPos[2] < detector->get_cathode()) quanta.electrons = 0;
    vector<double> scint2 =
      n.GetS2(quanta.electrons, truthPos, smearPos, driftTime, vD, j, field,
	      useTiming, verbosity, wf_time, wf_amp, g2_params);
  NEW_RANGES:
    if (usePD == 0 && fabs(scint[3]) > minS1 && scint[3] < maxS1)
      signal1.push_back(scint[3]);
    else if (usePD == 1 && fabs(scint[5]) > minS1 && scint[5] < maxS1)
      signal1.push_back(scint[5]);
    else if (usePD >= 2 && fabs(scint[7]) > minS1 && scint[7] < maxS1)
      signal1.push_back(scint[7]);
    else
      signal1.push_back(-999.);
    
    if (usePD == 0 && fabs(scint2[5]) > minS2 && scint2[5] < maxS2)
      signal2.push_back(scint2[5]);
    else if (usePD >= 1 && fabs(scint2[7]) > minS2 && scint2[7] < maxS2)
      signal2.push_back(scint2[7]);  // no spike option for S2
    else
      signal2.push_back(-999.);
    
    if ( eMin == eMax ) {
      if ( (scint[3] > maxS1 || scint[5] > maxS1 || scint[7] > maxS1) && j < 10 )
	cerr << "WARNING: Some S1 pulse areas are greater than maxS1" << endl;
      if ( (scint2[5] > maxS2 || scint2[7] > maxS2) && j < 10 ) //don't repeat too much: only if within first 10 events then show (+above)
	cerr << "WARNING: Some S2 pulse areas are greater than maxS2" << endl;
    }
    
    double Nph = 0.0, Ne = 0.0;
    if (!MCtruthE) {
      double MultFact = 1., eff = detector->get_sPEeff();
      if ( useTiming >= 0 ) {
	if ( detector->get_sPEthr() >= 0. && detector->get_sPEres() > 0. && eff > 0. ) {
	  MultFact = 0.5*(1.+erf((detector->get_sPEthr()-1.)/(detector->get_sPEres()*sqrt(2.)))) / (detector->get_sPEres()*sqrt(2.*M_PI));
          MultFact = 1./(1.-MultFact); MultFact = 1.02+(MultFact-1.02)/(1.+pow(scint[0]/detector->get_numPMTs(),2.)); //2% by Emily Mangus
	  if ( eff < 1. )
	    eff += ((1.-eff)/(1.*double(detector->get_numPMTs())))*scint[0];
	  if ( eff > 1. ) eff = 1.;
	  if ( eff < 0. ) eff = 0.;
	  MultFact /= eff; //for smaller detectors leave it as 1.00
	}
	else MultFact = 1.;
      }
      if (usePD == 0)
        Nph = fabs(scint[3]) / (g1 * (1. + detector->get_P_dphe()));
      else if (usePD == 1)
        Nph = fabs(scint[5]) / g1;
      else
        Nph = fabs(scint[7]) / g1;
      if (usePD == 0)
        Ne = fabs(scint2[5]) / (g2 * (1. + detector->get_P_dphe()));
      else
        Ne = fabs(scint2[7]) / g2;
      Nph *= MultFact;
      if (signal1.back() <= 0.) Nph = 0.;
      if (signal2.back() <= 0.) Ne = 0.;
      if ( detector->get_coinLevel() <= 0 && Nph <= PHE_MIN ) Nph = DBL_MIN;
      if (yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0.) {
	if ( detector->get_extraPhot() ) yields.Lindhard = 1.;
	if ( yields.Lindhard == 1. )
	  keV = (Nph/FudgeFactor[0] + Ne/FudgeFactor[1]) * Wq_eV * 1e-3;
	else {
	  if ( type_num <= NEST::INTERACTION_TYPE::Cf ) {
	    keV = pow((Ne+Nph)/NuisParam[0],1./NuisParam[1]);
	    Ne *= 1. - 1. / pow(1. + pow((keV / NuisParam[5]), NuisParam[6]), NuisParam[10]);
	    Nph *=1. - 1. / pow(1. + pow((keV / NuisParam[7]), NuisParam[8]), NuisParam[11]);
	    keV = pow((Ne+Nph)/NuisParam[0],1./NuisParam[1]);
	  }
	  else
	    keV =(Nph + Ne) * Wq_eV * 1e-3 / yields.Lindhard; // cheating :-)
	}
        keVee += (Nph + Ne) * Wq_eV * 1e-3;  // as alternative, use W_DEFAULT in
                                             // both places, but will not account
                                             // for density dependence
      } else
        keV = 0.;
    }
    if ((signal1.back() <= 0. || signal2.back() <= 0.) && field >= FIELD_MIN && detector->get_coinLevel() != 0)
      signalE.push_back(0.);
    else
      signalE.push_back(keV);
    
    if ( (Ne+Nph == 0.00 || std::isnan(keVee)) && eMin == eMax && eMin > hiEregime && !BeenHere ) { //efficiency is zero?
      minS1 = -999.;
      minS2 = -999.;
      detector->set_s2_thr(0.0); //since needs GetS2() re-run, only "catches" after the first caught event
      if ( maxS2 > 1e10 )
	BeenHere = true;
      maxS1 = 1e9;
      maxS2 = 1e11;
      numBins = 1;
      MCtruthE = false;
      signal1.pop_back();
      signal2.pop_back();
      signalE.pop_back();
      cerr << endl << "CAUTION: Efficiency seems to have been zero, so trying again with full S1 and S2 ranges." << endl;
      cerr << "OR, you tried to simulate a mono-energetic peak with MC truth E turned on. Silly! Setting MCtruthE to false." << endl;
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
    // scint[4] = Nphd; // same as pulse area, adjusted/corrected *downward* for
    // 2-PE effect (LUX phd units)
    // scint[5] = NphdC; // same as Nphd, but XYZ-corrected
    // scint[6] = spike; // floating real# spike count, NO XYZ correction
    // scint[7] = spikeC; // floating real# spike count, WITH XYZ correction
    // scint[8] = nHits post coincidence window and N-fold coincidence cuts

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
    // scint2[4] = pulseArea; // floating real# smeared DAQ pulse areas in phe,
    // NO XYZ correction
    // scint2[5] = pulseAreaC; // smeared DAQ pulse areas in phe, WITH XYZ
    // correction
    // scint2[6] = Nphd; // same as pulse area, adjusted/corrected *downward*
    // for 2-PE effect (LUX phd units)
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
    // scint2[8] = g2; // g2 = ExtEff * SE, light collection efficiency of EL in
    // gas gap (from CalculateG2)

    if (1) {  // fabs(scint[7]) > PHE_MIN && fabs(scint2[7]) > PHE_MIN ) { //if
      // you want to skip specific sub-threshold events, comment in this
      // if statement (save screen/disk space)
      // other suggestions: minS1, minS2 (or s2_thr) for tighter cuts depending
      // on analysis.hh settings (think of as analysis v. trigger thresholds)
      // and using max's too, pinching both ends
      if (type_num == Kr83m && eMin != 32.1 )
        printf("%.6f\t", yields.DeltaT_Scint);
      if (type_num == WIMP && timeStamp > (tZero+tStep))
	printf("%.0f\t", timeStamp);
      printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%d\t%d\t", keV, field,
             driftTime, smearPos[0], smearPos[1], smearPos[2], quanta.photons,
             quanta.electrons);  // comment this out when below line in
      // printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f,%.0f\t%lf\t%lf\t",keV,field,driftTime,smearPos[0],smearPos[1],smearPos[2],yields.PhotonYield,yields.ElectronYield);
      //for when you want means
      // if (truthPos[2] < detector->get_cathode() && verbosity) printf("g-X ");
      if (keV > 10.*hiEregime || scint[5] > maxS1 || scint2[7] > maxS2 ||
          // switch to exponential notation to make output more readable, if
          // energy is too high (>1 MeV)
          type == "muon" || type == "MIP" || type == "LIP" || type == "mu" ||
          type == "mu-") {
        printf("%e\t%e\t%e\t", scint[2], scint[5], scint[7]);
        printf("%li\t%e\t%e\n", (long)scint2[0], scint2[4], scint2[7]);
      } else {
        printf("%.6f\t%.6f\t%.6f\t", scint[2], scint[5],
               scint[7]);  // see GetS1 inside of NEST.cpp for full explanation
        // of all 8 scint return vector elements. Sample 3
        // most common
        printf("%i\t%.6f\t%.6f\n", (int)scint2[0], scint2[4],
               scint2[7]);  // see GetS2 inside of NEST.cpp for full explanation
                            // of all 8 scint2 vector elements. Change as you
                            // desire
      }
    }  // always execute statement, if(1) above, because if is just place-holder
       // in case you want to drop all sub-threshold data
  }
  
  if (verbosity) {
    if (eMin != eMax && type_num != Kr83m) {
      if (useS2 == 2)
        GetBand(signal2, signal1, false, detector->get_coinLevel());
      else
        GetBand(signal1, signal2, false, detector->get_coinLevel());
      fprintf(stderr,
              "Bin Center\tBin Actual\tHist Mean\tMean Error\tHist "
              "Sigma\tSkewness\t\tEff[%%>thr]\n");
      for (int j = 0; j < numBins; j++) {
        fprintf(stderr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t\t%lf\n", band[j][0],
                band[j][1], band[j][2], band[j][4], band[j][3], band[j][6],
                band[j][5] * 100.);
        if (band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 ||
            band[j][3] <= 0.0 || band[j][4] <= 0.0 || band[j][5] <= 0.0 ||
            std::isnan(band[j][0]) || std::isnan(band[j][1]) ||
            std::isnan(band[j][2]) || std::isnan(band[j][3]) ||
            std::isnan(band[j][4]) || std::isnan(band[j][5])) {
          if (eMax != -999.) {
            if (((g1 * yieldsMax.PhotonYield) < maxS1 ||
                 (g2 * yieldsMax.ElectronYield) < maxS2) &&
                j != 0)
              cerr << "WARNING: Insufficient number of high-energy events to "
                      "populate highest bins is likely.\n";
            else
              cerr << "WARNING: Insufficient number of low-energy events to "
                      "populate lowest bins is likely. Increase minS1 and/or "
                      "minS2.\n";
          }
          eMax = -999.;
        }
      }
    } else {
      GetBand(signal1, signal2, true, detector->get_coinLevel());
      GetEnergyRes(signalE);
      if (type_num == NR || type_num == WIMP || type_num == B8 || type_num == DD || type_num == AmBe || type_num == Cf || type_num == ion) {
        fprintf(stderr,
                "S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc "
                "[keVnr]\tEc Res[%%]\tEff[%%>thr]\tEc [keVee]\n");
        keVee /= numEvts;
      } else
        fprintf(stderr,
                "S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc Mean\t\tEc "
                "Res[%%]\tEff[%%>thr]\n");  // the C here refers to the combined
      // (S1+S2) energy scale
      for (int j = 0; j < numBins; j++) {
        fprintf(stderr, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", band[j][0],
                band[j][1] / band[j][0] * 100., band[j][2],
                band[j][3] / band[j][2] * 100., energies[0],
                energies[1] / energies[0] * 100., energies[2] * 100.);
        if (type_num < 7) //0-6=NR/related (WIMPs,etc.)
          fprintf(stderr, "%lf\n", keVee / energies[2]);
        else
          fprintf(stderr, "\n");
        if ((band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 ||
             band[j][3] <= 0.0 || std::isnan(band[j][0]) ||
             std::isnan(band[j][1]) || std::isnan(band[j][2]) ||
             std::isnan(band[j][3])) &&
            field >= FIELD_MIN) {
          if (numEvts > 1)
            cerr << "CAUTION: YOUR S1 and/or S2 MIN and/or MAX may be set to "
                    "be too restrictive, please check.\n";
          else
            cerr << "CAUTION: Poor stats. You must have at least 2 events to "
                    "calculate S1 and S2 and E resolutions.\n";
        } else if ((energies[0] == eMin || energies[0] == eMax ||
                    energies[1] <= 1E-6) &&
                   field >= FIELD_MIN)
          cerr << "If your energy resolution is 0% then you probably still "
                  "have MC truth energy on." << endl;
        else
          ;
      }
    }
  }
  return 0;
}

vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol, int nFold) {
  
  if ( numBins > NUMBINS_MAX ) {
    cerr << "ERROR: Too many bins. Decrease numBins (analysis.hh) or increase NUMBINS_MAX (TestSpectra.hh)" << endl;
    exit ( EXIT_FAILURE );
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
  unsigned long reject[NUMBINS_MAX] = {0};

  if (resol) {
    numBins = 1;
    binWidth = DBL_MAX;
  }
  
  double ReqS1 = 0.;
  if ( nFold == 0 )
    ReqS1 = -DBL_MAX;
  for (i = 0; i < S1s.size(); i++) {
    for (j = 0; j < numBins; j++) {
      s1c = border + binWidth / 2. + double(j) * binWidth;
      if (i == 0 && !resol) band[j][0] = s1c;
      if ( (fabs(S1s[i]) > (s1c - binWidth / 2.) &&
	    fabs(S1s[i]) <=(s1c + binWidth / 2.)) || !nFold ) {
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
          reject[j]++;
        break;
      }
    }
  }

  for (j = 0; j < numBins; j++) {
    if (band[j][0] <= 0. && !resol)
      band[j][0] = border + binWidth / 2. + double(j) * binWidth;
    signals[j].erase(signals[j].begin());
    numPts = (double)signals[j].size();
    if (numPts <= 0 && resol) {
      for (i = 0; i < S1s.size(); i++) band[j][0] += fabs(S1s[i]);
      numPts = S1s.size();
    }
    if (resol) band[j][0] /= numPts;
    band[j][1] /= numPts;
    band[j][2] /= numPts;
    if ( numPts > signals[j].size() )
      numPts = signals[j].size();  // seg fault prevention line
    for (i = 0; i < (int)numPts; i++) {
      if (signals[j][i] != -999.)
        band[j][3] += pow(signals[j][i] - band[j][2], 2.);  // std dev calc
    }
    for (i = 0; i < S1s.size(); i++) {
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
    for ( i = 0; i < (int)numPts; i++ ) {
      if ( signals[j][i] != -999. )
	band[j][6] += pow ( ( signals[j][i] - band[j][2] ) / band[j][3], 3. ); // skew calc
    }
    band[j][6] /= ( numPts - 2. );
  }

  return signals;
}

void GetEnergyRes(vector<double> Es) {
  int i, numPts = Es.size();
  double numerator = 0.;

  for (i = 0; i < numPts; i++) {
    if (Es[i] > 0.) {
      energies[0] += Es[i];
      numerator++;
    }
  }

  energies[0] /= numerator;

  for (i = 0; i < numPts; i++) {
    if (Es[i] > 0.) energies[1] += pow(energies[0] - Es[i], 2.);
  }

  energies[1] /= numerator - 1.;
  energies[1] = sqrt(energies[1]);

  energies[2] = numerator / double(numPts);
  return;
}
