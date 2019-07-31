
#include "NEST.hh"

using namespace std;
using namespace NEST;

const std::vector<double> NESTcalc::default_NuisParam = {11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.};
const std::vector<double> NESTcalc::default_FreeParam = {1.,1.,0.1,0.5,0.07};

long NESTcalc::BinomFluct(long N0, double prob) {
  double mean = N0 * prob;
  double sigma = sqrt(N0 * prob * (1. - prob));
  int N1 = 0;

  if (prob <= 0.00) return N1;
  if (prob >= 1.00) return N0;

  if (N0 < 10) {
    for (int i = 0; i < N0; i++) {
      if (RandomGen::rndm()->rand_uniform() < prob) N1++;
    }
  } else {
    N1 = int(floor(RandomGen::rndm()->rand_gauss(mean, sigma) + 0.5));
  }

  if (N1 > N0) N1 = N0;
  if (N1 < 0) N1 = 0;

  return N1;
}

NESTresult NESTcalc::FullCalculation(INTERACTION_TYPE species, double energy,
                                     double density, double dfield, double A,
                                     double Z,
                                     vector<double> NuisParam /*={11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.}*/,
				     vector<double> FreeParam /*={1.,1.,0.1,0.5,0.07}*/,
                                     bool do_times /*=true*/) {
  NESTresult result;
  result.yields = GetYields(species, energy, density, dfield, A, Z, NuisParam);
  result.quanta = GetQuanta(result.yields, density, FreeParam);
  if (do_times)
    result.photon_times = GetPhotonTimes(
        species, result.quanta.photons, result.quanta.excitons, dfield, energy);
  else
    result.photon_times = photonstream(result.quanta.photons, 0.0);
  return result;
}

double NESTcalc::PhotonTime(INTERACTION_TYPE species, bool exciton,
                            double dfield, double energy) {
  double time_ns = 0., SingTripRatio, tauR = 0., tau3 = 23.97,
         tau1 = 3.27;  // arXiv:1802.06162
  if (fdetector->get_inGas() ||
      energy < W_DEFAULT * 0.001) {  // from G4S1Light.cc in old NEST
    tau1 = 5.18;                     // uncertainty of 1.55 ns from G4S2Light
    tau3 = 100.1;                    // uncertainty of 7.90 ns from G4S2Light
  }
  // tau1 = 3.5*ns; tau3 = 20.*ns; tauR = 40.*ns for solid Xe from old NEST.
  // Here assuming same as in liquid

  if (species <= Cf)                           // NR
    SingTripRatio = 0.15 * pow(energy, 0.15);  // arXiv:1802.06162
  else if (species == ion)                     // e.g., alphas
    SingTripRatio =
        0.065 *
        pow(energy,
            0.416);  // spans 2.3 (alpha) and 7.8 (Cf in Xe) from NEST v1
  else {             // ER
    if (!exciton) {
      tauR = 0.5 * exp(-0.00900 * dfield) *
             (7.3138 + 3.8431 * log10(energy));    // arXiv:1310.1117
      SingTripRatio = 0.069 * pow(energy, -0.12);  // see comment below
    } else
      SingTripRatio =
          0.015 *
          pow(energy, -0.12);  // mixing arXiv:1802.06162 with Kubota 1979
  }

  if (fdetector->get_inGas() || energy < W_DEFAULT * 0.001) {
    SingTripRatio = 0.1;
    tauR = 0.;
  }

  // the recombination time is non-exponential, but approximates
  // to exp at long timescales (see Kubota '79)
  time_ns += tauR * (1.0 / RandomGen::rndm()->rand_uniform() - 1.);

  if (RandomGen::rndm()->rand_uniform() < SingTripRatio / (1. + SingTripRatio))
    time_ns -= tau1 * log(RandomGen::rndm()->rand_uniform());
  else
    time_ns -= tau3 * log(RandomGen::rndm()->rand_uniform());

  return time_ns;
}

// in analytical model, an empirical distribution (Brian Lenardo and Dev A.K.)
photonstream NESTcalc::AddPhotonTransportTime(photonstream emitted_times,
                                              double x, double y, double z) {
  photonstream return_photons;
  for (auto t : emitted_times) {
    double newtime = t + fdetector->OptTrans(x, y, z);
    return_photons.push_back(newtime);
  }
  return return_photons;
}

photonstream NESTcalc::GetPhotonTimes(INTERACTION_TYPE species,
                                      int total_photons, int excitons,
                                      double dfield, double energy) {
  photonstream return_photons;
  int total, excit;

  for (int ip = 0; ip < total_photons; ++ip) {
    bool isExciton = false;
    if (ip < excitons) isExciton = true;
    return_photons.push_back(PhotonTime(species, isExciton, dfield, energy));
  }

  return return_photons;
}

QuantaResult NESTcalc::GetQuanta(YieldResult yields, double density,
				 vector<double> FreeParam/*={1.,1.,0.1,0.5,0.07}*/) {
  QuantaResult result;
  bool HighE;
  int Nq_actual, Ne, Nph, Ni, Nex;

  double NexONi = yields.ExcitonRatio, Fano = 1.;
  double Nq_mean = yields.PhotonYield + yields.ElectronYield;

  double elecFrac = yields.ElectronYield / Nq_mean;
  if (elecFrac > 1.) elecFrac = 1.;
  if (elecFrac < 0.) elecFrac = 0.;

  if (NexONi < 0.) {
    NexONi = 0.;
    HighE = true;
  } else
    HighE = false;
  double alf = 1. / (1. + NexONi);
  double recombProb = 1. - (NexONi + 1.) * elecFrac;
  if (recombProb < 0.) NexONi = 1. / elecFrac - 1.;

  if (yields.Lindhard == 1.) {
    Fano = 0.12707 - 0.029623 * density -  // Fano factor is  << 1
           0.0057042 *
               pow(density,
                   2.) +  //~0.1 for GXe w/ formula from Bolotnikov et al. 1995
           0.0015957 *
               pow(density,
                   3.);  // to get it to be ~0.03 for LXe (E Dahl Ph.D. thesis)
    if (!fdetector->get_inGas())
      Fano += 0.0015 * sqrt(Nq_mean) * pow(yields.ElectricField, 0.5);
    Nq_actual = int(floor(
        RandomGen::rndm()->rand_gauss(Nq_mean, sqrt(Fano * Nq_mean)) + 0.5));
    if (Nq_actual < 0 || Nq_mean == 0.) Nq_actual = 0;

    Ni = BinomFluct(Nq_actual, alf);
    Nex = Nq_actual - Ni;

  }

  else {
    Fano = FreeParam[0];
    Ni = int(floor(RandomGen::rndm()->rand_gauss(Nq_mean * alf,
                                                 sqrt(Fano * Nq_mean * alf)) +
                   0.5));
    if (Ni < 0) Ni = 0;
    Fano = FreeParam[1];
    Nex = int(
        floor(RandomGen::rndm()->rand_gauss(
                  Nq_mean * NexONi * alf, sqrt(Fano * Nq_mean * NexONi * alf)) +
              0.5));
    if (Nex < 0) Nex = 0;
    Nq_actual = Nex + Ni;
  }

  if (Nq_actual == 0) {
    result.ions = 0;
    result.excitons = 0;
    result.photons = 0;
    result.electrons = 0;
    return result;
  }

  if (Nex < 0) Nex = 0;
  if (Ni < 0) Ni = 0;
  if (Nex > Nq_actual) Nex = Nq_actual;
  if (Ni > Nq_actual) Ni = Nq_actual;

  result.ions = Ni;
  result.excitons = Nex;

  if (Nex <= 0 && HighE) recombProb = yields.PhotonYield / double(Ni);
  if (recombProb < 0.) recombProb = 0.;
  if (recombProb > 1.) recombProb = 1.;
  if (std::isnan(recombProb) || std::isnan(elecFrac) || Ni == 0 ||
      recombProb == 0.0) {
    recombProb = 0.0;
    elecFrac = 1.0;
    result.photons = Nex;
    result.electrons = Ni;
    return result;
  }

  double ef = yields.ElectricField;
  double cc =
      0.075351 +
      (0.050461 - 0.075351) / pow(1. + pow(ef / 30057., 3.0008), 2.9832e5);
  if (cc < 0.) cc = 0.;
  double bb = 0.54;
  double aa = cc / pow(1. - bb, 2.);
  double omega = -aa * pow(recombProb - bb, 2.) + cc;
  if (omega < 0.0) omega = 0.0;

  if (yields.Lindhard < 1.)
    omega = FreeParam[2] * exp(-pow(elecFrac - FreeParam[3], 2.) / FreeParam[4]);
  double Variance =
      recombProb * (1. - recombProb) * Ni + omega * omega * Ni * Ni;
  Ne = int(floor(
      RandomGen::rndm()->rand_gauss((1. - recombProb) * Ni, sqrt(Variance)) +
      0.5));
  if (Ne < 0) Ne = 0;
  if (Ne > Ni) Ne = Ni;

  Nph = Nq_actual - Ne;
  if (Nph > Nq_actual) Nph = Nq_actual;
  if (Nph < Nex) Nph = Nex;

  if ((Nph + Ne) != (Nex + Ni)) {
    cerr << "\nERROR: Quanta not conserved. Tell Matthew Immediately!\n";
    exit(1);
  }

  result.photons = Nph;
  result.electrons = Ne;

  return result;  // quanta returned with recomb fluctuations
}

YieldResult NESTcalc::GetYields(INTERACTION_TYPE species, double energy,
                                double density, double dfield, double massNum,
                                double atomNum, vector<double> NuisParam
				/*={11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.}*/) {
  // For temporary variables for storing results
  double Ne = -999;
  double Nph = -999;
  double NexONi = -999, deltaT_ns = -999;
  double m8 = 2., L = 1.;
  const double deltaT_ns_halflife = 154.4;

  double Wq_eV =
      1.9896 + (20.8 - 1.9896) / (1. + pow(density / 4.0434, 1.4407));
  double alpha = 0.067366 + density * 0.039693;
  switch (species) {
    case NR:
    case WIMP:
    case B8:
    case DD:
    case AmBe:
    case Cf:  // this doesn't mean all NR is Cf, this is like a giant if
              // statement. Same intrinsic yields, but different energy spectra
              // (TestSpectra)
      {
        int massNumber;
        double ScaleFactor[2] = {1., 1.};
        if (massNum != 0.)
          massNumber = int(massNum);
        else
          massNumber = RandomGen::rndm()->SelectRanXeAtom();
        ScaleFactor[0] = sqrt(MOLAR_MASS / (double)massNumber);
        ScaleFactor[1] = ScaleFactor[0];
        double Nq = NuisParam[0] * pow(energy, NuisParam[1]);
        double ThomasImel =
            NuisParam[2] * pow(dfield, NuisParam[3]) * pow(density / DENSITY, 0.3);
        double Qy = 1. / (ThomasImel*pow(energy+NuisParam[4],NuisParam[9]));
        Qy *= 1. - 1. / pow(1. + pow((energy / NuisParam[5]), NuisParam[6]),NuisParam[10]);
        double Ly = Nq / energy - Qy;
        if (Qy < 0.0) Qy = 0.0;
        if (Ly < 0.0) Ly = 0.0;
        Ne = Qy * energy * ScaleFactor[1];
        Nph = Ly * energy * ScaleFactor[0] *
              (1. - 1. / (1. + pow((energy / NuisParam[7]), NuisParam[8])));
        Nq = Nph + Ne;
        double Ni = (4. / ThomasImel) * (exp(Ne * ThomasImel / 4.) - 1.);
        double Nex = (-1. / ThomasImel) * (4. * exp(Ne * ThomasImel / 4.) -
                                           (Ne + Nph) * ThomasImel - 4.);
        if (fabs(Nex - (Nq - Ni)) > PHE_MIN ||
            fabs(Ni - (Nq - Nex)) > PHE_MIN) {
          cerr << "\nERROR: Quanta not conserved. Tell Matthew Immediately!\n";
          exit(1);
        }
        NexONi = Nex / Ni;
        L = (Nq / energy) * Wq_eV * 1e-3;
      }
      break;
    case ion: {
      double A1 = massNum, A2 = RandomGen::rndm()->SelectRanXeAtom();
      double Z1 = atomNum, Z2 = ATOM_NUM;
      double Z_mean = pow(pow(Z1, (2. / 3.)) + pow(Z2, (2. / 3.)), 1.5);
      double E1c = pow(A1, 3.) * pow(A1 + A2, -2.) * pow(Z_mean, (4. / 3.)) *
                   pow(Z1, (-1. / 3.)) * 500.;
      double E2c = pow(A1 + A2, 2.) * pow(A1, -1.) * Z2 * 125.;
      double gamma = 4. * A1 * A2 / pow(A1 + A2, 2.);
      double Ec_eV = gamma * E2c;
      double Constant =
          (2. / 3.) * (1. / sqrt(E1c) + 0.5 * sqrt(gamma / Ec_eV));
      L = Constant * sqrt(energy * 1e3);
      double L_max = 0.96446 / (1. + pow(massNum * massNum / 19227., 0.99199));
      if (atomNum == 2. && massNum == 4.) L = 0.56136 * pow(energy, 0.056972);
      if (L > L_max) L = L_max;
      double densDep = pow(density / 0.2679, -2.3245);
      double massDep =
          0.02966094 * exp(0.17687876 * (massNum / 4. - 1.)) + 1. - 0.02966094;
      double fieldDep = pow(1. + pow(dfield / 95., 8.7), 0.0592);
      if (fdetector->get_inGas()) fieldDep = sqrt(dfield);
      double ThomasImel = 0.00625 * massDep / (1. + densDep) / fieldDep;
      const double logden = log10(density);
      Wq_eV = 28.259 + 25.667 * logden - 33.611 * pow(logden, 2.) -
              123.73 * pow(logden, 3.) - 136.47 * pow(logden, 4.) -
              74.194 * pow(logden, 5.) - 20.276 * pow(logden, 6.) -
              2.2352 * pow(logden, 7.);
      alpha = 0.64 / pow(1. + pow(density / 10., 2.), 449.61);
      NexONi = alpha + 0.00178 * pow(atomNum, 1.587);
      double Nq = 1e3 * L * energy / Wq_eV;
      double Ni = Nq / (1. + NexONi);
      double recombProb;
      if (Ni > 0. && ThomasImel > 0.)
        recombProb =
            1. - log(1. + (ThomasImel / 4.) * Ni) / ((ThomasImel / 4.) * Ni);
      else
        recombProb = 0.0;
      Nph = Nq * NexONi / (1. + NexONi) + recombProb * Ni;
      Ne = Nq - Nph;
    } break;
    case gammaRay: {
      const double m3 = 2., m4 = 2., m6 = 0.;
      double m1 =
          33.951 + (3.3284 - 33.951) / (1. + pow(dfield / 165.34, .72665));
      double m2 = 1000 / Wq_eV;
      double m5 =
          23.156 + (10.737 - 23.156) / (1. + pow(dfield / 34.195, .87459));
      double densCorr = 240720. / pow(density, 8.2076);
      double m7 =
          66.825 + (829.25 - 66.825) / (1. + pow(dfield / densCorr, .83344));
      double Nq = energy * 1000. / Wq_eV;
      if (fdetector->get_inGas()) m8 = -2.;
      double Qy = m1 + (m2 - m1) / (1. + pow(energy / m3, m4)) + m5 +
                  (m6 - m5) / (1. + pow(energy / m7, m8));
      double Ly = Nq / energy - Qy;
      Ne = Qy * energy;
      Nph = Ly * energy;
      NexONi = alpha * erf(0.05 * energy);
    } break;
    case Kr83m: {
      double Nq = 0.;
      if (energy == 9.4) {
        deltaT_ns = RandomGen::rndm()->rand_exponential(deltaT_ns_halflife);
        Nq = energy * (1e3 / Wq_eV + 6.5);
        double medTlevel =
            47.8 + (69.201 - 47.8) / pow(1. + pow(dfield / 250.13, 0.9), 1.);
        double highTrise =
            1.15 + (1. - 1.15) / (1. + pow(deltaT_ns / 1200., 18.));
        double lowTdrop = 14. * pow(dfield, 0.19277);
        Nph = energy * highTrise *
              (5.1e4 * pow(2. * deltaT_ns + 10., -1.5) + medTlevel) /
              (1. + pow(deltaT_ns / lowTdrop, -3.));
        alpha = 0.;
      } else {
        Nq = energy * 1000. / Wq_eV;
        Nph =
            energy *
            (6. + (69.742 - 6.) / pow(1. + pow(dfield / 9.515, 1.9), 0.063032));
      }
      Ne = Nq - Nph;
      NexONi = alpha * erf(0.05 * energy);
    } break;
    default:  // beta, CH3T
    {
      double QyLvllowE =
          1e3 / Wq_eV + 6.5 * (1. - 1. / (1. + pow(dfield / 47.408, 1.9851)));
      double HiFieldQy =
          1. + 0.4607 / pow(1. + pow(dfield / 621.74, -2.2717), 53.502);
      double QyLvlmedE =
          32.988 -
          32.988 /
              (1. + pow(dfield / (0.026715 * exp(density / 0.33926)), 0.6705));
      QyLvlmedE *= HiFieldQy;
      double DokeBirks =
          1652.264 +
          (1.415935e10 - 1652.264) / (1. + pow(dfield / 0.02673144, 1.564691));
      double Nq = energy * 1e3 /
                  Wq_eV;  //( Wq_eV+(12.578-Wq_eV)/(1.+pow(energy/1.6,3.5)) );
      double LET_power = -2.;
      if (fdetector->get_inGas()) LET_power = 2.;
      double QyLvlhighE = 28.;
      //      if (density > 3.) QyLvlhighE = 49.; Solid Xe effect from Yoo. But,
      //      beware of enabling this line: enriched liquid Xe for neutrinoless
      //      double beta decay has density higher than 3g/cc;
      double Qy = QyLvlmedE +
                  (QyLvllowE - QyLvlmedE) /
                      pow(1. + 1.304 * pow(energy, 2.1393), 0.35535) +
                  QyLvlhighE / (1. + DokeBirks * pow(energy, LET_power));
      if (Qy > QyLvllowE && energy > 1. && dfield > 1e4) Qy = QyLvllowE;
      double Ly = Nq / energy - Qy;
      Ne = Qy * energy;
      Nph = Ly * energy;
      NexONi = alpha * erf(0.05 * energy);
    } break;
  }

  assert(Ne != -999 && Nph != -999 && NexONi != -999);
  if (Nph > energy / W_SCINT)
    Nph = energy / W_SCINT;  // yields can never exceed 1 / [ W ~ few eV ]
  if (Ne > energy / W_SCINT) Ne = energy / W_SCINT;
  if (Nph < 0.) Nph = 0.;
  if (Ne < 0.) Ne = 0.;
  // if (NexONi < 0.) NexONi = 0.;
  if (L < 0.) L = 0.;
  if (L > 1.) L = 1.;  // Lindhard Factor
  if (energy < 0.001 * Wq_eV / L) {
    Nph = 0.;
    Ne = 0.;
  }

  YieldResult result;
  result.PhotonYield = Nph;
  result.ElectronYield = Ne;
  result.ExcitonRatio = NexONi;
  result.Lindhard = L;
  result.ElectricField = dfield;
  result.DeltaT_Scint = deltaT_ns;
  return result;  // everything needed to calculate fluctuations
}

NESTcalc::NESTcalc() { fdetector = NULL; }

NESTcalc::NESTcalc(VDetector* detector) {
  NESTcalc();
  fdetector = detector;
}

NESTcalc::~NESTcalc() {
  if (pulseFile) pulseFile.close();
  if (fdetector) delete fdetector;
}

vector<double> NESTcalc::GetS1(QuantaResult quanta, double truthPos[3],
                               double smearPos[3], double driftVelocity,
                               double dV_mid, INTERACTION_TYPE type_num,
                               long evtNum, double dfield, double energy,
                               int useTiming, bool outputTiming,
                               vector<long int>& wf_time,
                               vector<double>& wf_amp) {
  int Nph = quanta.photons;

  wf_time.clear();
  wf_amp.clear();

  vector<double> photon_areas[2];
  vector<double> scintillation(9);  // return vector
  vector<double> newSpike(2);  // for re-doing spike counting more precisely

  // Add some variability in g1 drawn from a polynomial spline fit
  double posDep = fdetector->FitS1(truthPos[0], truthPos[1], truthPos[2]);
  double posDepSm = fdetector->FitS1(smearPos[0], smearPos[1], smearPos[2]);
  double dt = (fdetector->get_TopDrift() - truthPos[2]) / driftVelocity;
  double dz_center = fdetector->get_TopDrift() -
                     dV_mid * fdetector->get_dtCntr();  // go from t to z
  posDep /=
      fdetector->FitS1(0., 0., dz_center);  // XYZ always in mm now never cm
  posDepSm /= fdetector->FitS1(0., 0., dz_center);

  // generate a number of PMT hits drawn from a binomial distribution.
  // Initialize number of photo-electrons
  int nHits = BinomFluct(Nph, fdetector->get_g1() * posDep), Nphe = 0;

  // Initialize the pulse area and spike count variables
  double pulseArea = 0., spike = 0., prob;

  // If single photo-electron efficiency is under 1 and the threshold is above 0
  // (some phe will be below threshold)
  if (useTiming ||
      (fdetector->get_sPEthr() > 0. &&
       nHits < fdetector->get_numPMTs())) {  // digital nHits eventually becomes
                                             // spikes (spike++) based upon
                                             // threshold
    // Step through the pmt hits
    for (int i = 0; i < nHits; i++) {
      // generate photo electron, integer count and area
      double phe1 = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres()) +
                    RandomGen::rndm()->rand_gauss(fdetector->get_noise()[0],
                                                  fdetector->get_noise()[1]);
      Nphe++;
      if (phe1 > DBL_MAX) phe1 = 1.;
      if (phe1 < -DBL_MAX) phe1 = 0.;
      prob = RandomGen::rndm()->rand_uniform();
      // zero the area if random draw determines it wouldn't have been observed.
      if (prob > fdetector->get_sPEeff()) {
        phe1 = 0.;
      }  // add an else with Nphe++ if not doing mc truth
      // Generate a double photo electron if random draw allows it
      double phe2 = 0.;
      if (RandomGen::rndm()->rand_uniform() < fdetector->get_P_dphe()) {
        // generate area and increment the photo-electron counter
        phe2 = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres()) +
               RandomGen::rndm()->rand_gauss(fdetector->get_noise()[0],
                                             fdetector->get_noise()[1]);
        Nphe++;
        if (phe2 > DBL_MAX) phe2 = 1.;
        if (phe2 < -DBL_MAX) phe2 = 0.;
        // zero the area if phe wouldn't have been observed
        if (RandomGen::rndm()->rand_uniform() > fdetector->get_sPEeff() &&
            prob > fdetector->get_sPEeff()) {
          phe2 = 0.;
        }  // add an else with Nphe++ if not doing mc truth
        // The dphe occurs simultaneously to the first one from the same source
        // photon. If the first one is seen, so should be the second one
      }
      // Save the phe area and increment the spike count (very perfect spike
      // count) if area is above threshold
      if (useTiming) {
        if ((phe1 + phe2) > fdetector->get_sPEthr()) {
          pulseArea += phe1 + phe2;
          spike++;
          photon_areas[0].push_back(phe1);
          photon_areas[1].push_back(phe2);
        }
      } else {  // use approximation to find timing
        if ((phe1 + phe2) > fdetector->get_sPEthr() &&
            (-20. * log(RandomGen::rndm()->rand_uniform()) <
                 fdetector->get_coinWind() ||
             nHits > fdetector->get_coinLevel())) {
          spike++;
          pulseArea += phe1 + phe2;
        } else
          ;
      }
    }
  } else {  // apply just an empirical efficiency by itself, without direct area
            // threshold
    Nphe = nHits + BinomFluct(nHits, fdetector->get_P_dphe());
    double eff = fdetector->get_sPEeff();
    if (nHits >= fdetector->get_numPMTs()) eff = 1.;
    pulseArea = RandomGen::rndm()->rand_gauss(
        BinomFluct(Nphe, 1. - (1. - eff) / (1. + fdetector->get_P_dphe())),
        fdetector->get_sPEres() * sqrt(Nphe));
    spike = (double)nHits;
  }

  if (useTiming) {
    vector<double> PEperBin, AreaTable[2], TimeTable[2];
    int numPts = 1100 - 100 * SAMPLE_SIZE;
    AreaTable[0].resize(numPts, 0.);
    AreaTable[1].resize(numPts, 0.);

    int total_photons = (int)fabs(spike);
    int excitons =
        int((double(nHits) / double(quanta.photons)) * double(quanta.excitons) +
            0.5);
    photonstream photon_emission_times =
        GetPhotonTimes(type_num, total_photons, excitons, dfield, energy);
    photonstream photon_times = AddPhotonTransportTime(
        photon_emission_times, truthPos[0], truthPos[1], truthPos[2]);

    if (outputTiming && !pulseFile.is_open()) {
      pulseFile.open("photon_times.txt");
      // pulseFile << "Event#\tt [ns]\tA [PE]" << endl;
      pulseFile << "Event#\tt [ns]\tPEb/bin\tPEt/bin\tin win" << endl;
    }

    int ii, index;
    double min = 1e100, pTime;
    for (ii = 0; ii < (int)fabs(spike); ++ii) {
      PEperBin.clear();
      PEperBin = fdetector->SinglePEWaveForm(
          photon_areas[0][ii] + photon_areas[1][ii], photon_times[ii]);
      int total = (unsigned int)PEperBin.size() - 1;
      int whichArray;
      if (RandomGen::rndm()->rand_uniform() <
          fdetector->FitTBA(truthPos[0], truthPos[1], truthPos[2])[0])
        whichArray = 0;
      else
        whichArray = 1;
      for (int kk = 0; kk < total; ++kk) {
        pTime = PEperBin[0] + kk * SAMPLE_SIZE;
        index = int(floor(pTime / SAMPLE_SIZE + 0.5)) + numPts / 2;
        if (index < 0) index = 0;
        if (index >= numPts) index = numPts - 1;
        AreaTable[whichArray][index] +=
            10. * (1. + PULSEHEIGHT) *
            (photon_areas[0][ii] + photon_areas[1][ii]) /
            (PULSE_WIDTH * sqrt(2. * M_PI)) *
            exp(-pow(pTime - photon_times[ii], 2.) /
                (2. * PULSE_WIDTH * PULSE_WIDTH));
      }
      if (total >= 0) {
        if (PEperBin[0] < min) min = PEperBin[0];
        TimeTable[0].push_back(PEperBin[0]);
      }
      // else
      // TimeTable[0].push_back(-999.);
      // TimeTable[1].push_back(photon_areas[0][ii]+photon_areas[1][ii]);
    }
    for (ii = 0; ii < numPts; ++ii) {
      if ((AreaTable[0][ii] + AreaTable[1][ii]) <= PULSEHEIGHT) continue;

      wf_time.push_back((ii - numPts / 2) * SAMPLE_SIZE);
      wf_amp.push_back(AreaTable[0][ii] + AreaTable[1][ii]);

      if (outputTiming) {
        char line[80];
        sprintf(line, "%lu\t%ld\t%.2f\t%.2f", evtNum, wf_time.back(),
                AreaTable[0][ii], AreaTable[1][ii]);
        pulseFile << line << flush;
      }

      if (((ii - numPts / 2) * SAMPLE_SIZE - (int)min) >
              fdetector->get_coinWind() &&
          nHits <= fdetector->get_coinLevel()) {
        pulseArea -= (AreaTable[0][ii] + AreaTable[1][ii]);
        pulseFile << "\t0" << endl;
      } else
        pulseFile << "\t1" << endl;
    }
    for (ii = 0; ii < TimeTable[0].size(); ++ii) {
      if ((TimeTable[0][ii] - min) > fdetector->get_coinWind() &&
          nHits <= fdetector->get_coinLevel())
        --spike;
      // char line[80]; sprintf ( line, "%lu\t%.1f\t%.2f", evtNum,
      // TimeTable[0][ii], TimeTable[1][ii] ); pulseFile << line << endl;
    }
  }

  pulseArea = RandomGen::rndm()->rand_gauss(
      pulseArea, fdetector->get_noise()[2] * pulseArea);
  if (pulseArea < fdetector->get_sPEthr()) pulseArea = 0.;
  if (spike < 0) spike = 0;
  double pulseAreaC = pulseArea / posDepSm;
  double Nphd = pulseArea / (1. + fdetector->get_P_dphe());
  double NphdC = pulseAreaC / (1. + fdetector->get_P_dphe());
  double spikeC = spike / posDepSm;

  scintillation[0] = nHits;  // MC-true integer hits in same OR different PMTs,
                             // NO double phe effect
  scintillation[1] =
      Nphe;  // MC-true integer hits WITH double phe effect (Nphe > nHits)
  scintillation[2] = pulseArea;  // floating real# smeared DAQ pulse areas in
                                 // phe, NO XYZ correction
  scintillation[3] =
      pulseAreaC;  // smeared DAQ pulse areas in phe, WITH XYZ correction
  scintillation[4] = Nphd;  // same as pulse area, adjusted/corrected *downward*
                            // for 2-PE effect (LUX phd units)
  scintillation[5] = NphdC;   // same as Nphd, but XYZ-corrected
  scintillation[6] = spike;   // floating real# spike count, NO XYZ correction
  scintillation[7] = spikeC;  // floating real# spike count, WITH XYZ correction
  scintillation[8] = spike;

  if (spike < fdetector->get_coinLevel())  // no chance of meeting coincidence
                                           // requirement. Here, spike is still
                                           // "perfect" (integer)
    prob = 0.;
  else {
    if (spike > 10.)
      prob = 1.;
    else {
      if (fdetector->get_coinLevel() == 0)
        prob = 1.;
      else if (fdetector->get_coinLevel() == 1) {
        if (spike >= 1.)
          prob = 1.;
        else
          prob = 0.;
      } else if (fdetector->get_coinLevel() == 2)
        prob = 1. - pow((double)fdetector->get_numPMTs(), 1. - spike);
      else {
        double numer = 0., denom = 0.;
        for (int i = spike; i > 0; i--) {
          denom += nCr(fdetector->get_numPMTs(), i);
          if (i >= fdetector->get_coinLevel())
            numer += nCr(fdetector->get_numPMTs(), i);
        }
        prob = numer / denom;
      }  // end of case of coinLevel of 3 or higher
    }    // end of case of spike is equal to 9 or lower
  }      // the end of case of spike >= coinLevel

  if (spike >= 1. && spikeC > 0.) {
    // over-write spike with smeared version, ~real data. Last chance to read
    // out digital spike and spikeC in scintillation[6] and [7]
    newSpike = GetSpike(Nph, smearPos[0], smearPos[1], smearPos[2],
                        driftVelocity, dV_mid, scintillation);
    scintillation[6] = newSpike[0];  // S1 spike count (NOT adjusted for double
                                     // phe effect), IF sufficiently small nHits
                                     // (otherwise = Nphd)
    scintillation[7] =
        newSpike[1];  // same as newSpike[0], but WITH XYZ correction
  }

  if (RandomGen::rndm()->rand_uniform() < prob ||
      prob >= 1.)  // coincidence has to happen in different PMTs
  {
    ;
  } else {  // some of these are set to -1 to flag them as having been below
            // threshold
    if (scintillation[0] == 0.) scintillation[0] = PHE_MIN;
    scintillation[0] *= -1.;
    if (scintillation[1] == 0.) scintillation[1] = PHE_MIN;
    scintillation[1] *= -1.;
    if (scintillation[2] == 0.) scintillation[2] = PHE_MIN;
    scintillation[2] *= -1.;
    if (scintillation[3] == 0.) scintillation[3] = PHE_MIN;
    scintillation[3] *= -1.;
    if (scintillation[4] == 0.) scintillation[4] = PHE_MIN;
    scintillation[4] *= -1.;
    if (scintillation[5] == 0.) scintillation[5] = PHE_MIN;
    scintillation[5] *= -1.;
    if (scintillation[6] == 0.) scintillation[6] = PHE_MIN;
    scintillation[6] *= -1.;
    if (scintillation[7] == 0.) scintillation[7] = PHE_MIN;
    scintillation[7] *= -1.;
  }

  // scintillation[8] =
  //  fdetector->get_g1();  // g1 (light collection efficiency in liquid)

  return scintillation;
}

vector<double> NESTcalc::GetS2(int Ne, double truthPos[3], double smearPos[3],
                               double dt, double driftVelocity, long evtNum,
                               double dfield, int useTiming, bool outputTiming,
                               vector<long int>& wf_time,
                               vector<double>& wf_amp,
                               vector<double>& g2_params) {
  double elYield = g2_params[0];
  double ExtEff = g2_params[1];
  double SE = g2_params[2];
  double g2 = g2_params[3];
  double gasGap = g2_params[4];

  vector<double> ionization(9);
  int i;
  bool eTrain = false;
  if (useTiming >= 2 && !fdetector->get_inGas()) eTrain = true;

  if (dfield < FIELD_MIN  //"zero"-drift-field detector has no S2
      || elYield <= 0. || ExtEff <= 0. || SE <= 0. || g2 <= 0. ||
      gasGap <= 0.) {
    for (i = 0; i < 8; i++) ionization[i] = 0.;
    return ionization;
  }

  // Add some variability in g1_gas drawn from a polynomial spline fit
  double posDep = fdetector->FitS2(
      truthPos[0], truthPos[1]);  // XY is always in mm now, never cm
  double posDepSm = fdetector->FitS2(smearPos[0], smearPos[1]);
  posDep /= fdetector->FitS2(0., 0.);
  posDepSm /= fdetector->FitS2(0., 0.);
  double dz = fdetector->get_TopDrift() - dt * driftVelocity;

  int Nee = BinomFluct(
      Ne, ExtEff * exp(-dt / fdetector->get_eLife_us()));  // MAKE this 1 for
                                                           // SINGLE e-
                                                           // DEBUGGING
  long Nph = 0, nHits = 0, Nphe = 0;
  double pulseArea = 0.;

  if (useTiming) {
    long k;
    int stopPoint;
    double tau1, tau2, E_liq, amp2;
    vector<double> electronstream, AreaTableBot[2], AreaTableTop[2], TimeTable;
    if (eTrain)
      stopPoint = BinomFluct(Ne, exp(-dt / fdetector->get_eLife_us()));
    else
      stopPoint = Nee;
    electronstream.resize(stopPoint, dt);
    double elecTravT = 0., DL, DL_time, DT, phi, sigX, sigY, newX, newY;
    double Diff_Tran = 37.368 * pow(dfield, .093452) *
                       exp(-8.1651e-5 * dfield);  // arXiv:1609.04467 (EXO-200)
    double Diff_Long = 345.92 * pow(dfield, -0.47880) *
                       exp(-81.3230 / dfield);  // fit to Aprile & Doke review
                                                // paper and to arXiv:1102.2865;
                                                // plus, LUX Run02+03
    // a good rule of thumb but only for liquids, as gas kind of opposite:
    // Diff_Long ~ 0.15 * Diff_Tran, as it is in LAr, at least as field goes to
    // infinity
    double Diff_Long_Gas =
        4.265 + 19097. / (1e3 * fdetector->get_E_gas()) -
        1.7397e6 / pow(1e3 * fdetector->get_E_gas(), 2.) +
        1.2477e8 / pow(1e3 * fdetector->get_E_gas(), 3.);  // Nygren, NEXT
    double Diff_Tran_Gas = Diff_Tran * 0.01;
    if (fdetector->get_inGas()) {
      Diff_Tran *= 0.01;
      Diff_Long = 4.265 + 19097. / dfield - 1.7397e6 / pow(dfield, 2.) +
                  1.2477e8 / pow(dfield, 3.);
    }
    double sigmaDL =
        10. * sqrt(2. * Diff_Long * dt *
                   1e-6);  // sqrt of cm^2/s * s = cm; times 10 for mm.
    double sigmaDT = 10. * sqrt(2. * Diff_Tran * dt * 1e-6);
    double rho = fdetector->get_p_bar() * 1e5 /
                 (fdetector->get_T_Kelvin() * 8.314) * MOLAR_MASS * 1e-6;
    double driftVelocity_gas =
        SetDriftVelocity_MagBoltz(rho, fdetector->get_E_gas() * 1000.);
    double dt_gas = gasGap / driftVelocity_gas;
    double sigmaDLg = 10. * sqrt(2. * Diff_Long_Gas * dt_gas * 1e-6);
    double sigmaDTg = 10. * sqrt(2. * Diff_Tran_Gas * dt_gas * 1e-6);
    double tauTrap = 0.185;  // microseconds from arXiv:1310.1117, modified to
                             // better fit XENON10 and LUX data at same time
    double min = 1e100;
    for (i = 0; i < stopPoint; ++i) {
      elecTravT = 0.;  // resetting for the current electron
      DL = RandomGen::rndm()->rand_gauss(0., sigmaDL);
      DT = fabs(RandomGen::rndm()->rand_gauss(0., sigmaDT));
      phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
      sigX = DT * cos(phi);
      sigY = DT * sin(phi);
      newX = truthPos[0] + sigX;
      newY = truthPos[1] + sigY;
      if (newX * newX + newY * newY >
          fdetector->get_radmax() * fdetector->get_radmax())
        continue;  // remove electrons from outside the maximum possible radius
                   // (alternative is squeeze them, depending on your E-fields)
      DL_time = DL / driftVelocity;
      elecTravT += DL_time;
      if (!fdetector->get_inGas() && fdetector->get_E_gas() != 0.)
        elecTravT -= tauTrap * log(RandomGen::rndm()->rand_uniform());
      electronstream[i] += elecTravT;
      // char line[80]; sprintf ( line, "%lu\t%.0f\t%.3f\t%.3f", evtNum,
      // electronstream[i]*1e+3, newX, newY ); pulseFile << line << endl;
      SE = floor(RandomGen::rndm()->rand_gauss(
                     elYield, sqrt(fdetector->get_s2Fano() * elYield)) +
                 0.5);
      Nph += long(SE);
      SE = (double)BinomFluct(long(SE), fdetector->get_g1_gas() * posDep);
      nHits += long(SE);
      double KE = 0.5 * ELEC_MASS * driftVelocity_gas * driftVelocity_gas *
                  1e6 / 1.602e-16;
      double origin = fdetector->get_TopDrift() + gasGap / 2.;
      QuantaResult quanta;
      quanta.photons = int(SE);
      quanta.electrons = 0;
      quanta.ions = 0;
      quanta.excitons = int(floor(0.0566 * SE + 0.5));
      photonstream photon_emission_times = GetPhotonTimes(
          INTERACTION_TYPE::beta, quanta.photons, quanta.excitons, dfield, KE);
      photonstream photon_times =
          AddPhotonTransportTime(photon_emission_times, newX, newY, origin);
      SE += (double)BinomFluct(long(SE), fdetector->get_P_dphe());
      Nphe += long(SE);
      DL = RandomGen::rndm()->rand_gauss(0., sigmaDLg);
      DT = RandomGen::rndm()->rand_gauss(0., sigmaDTg);
      phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
      sigX = DT * cos(phi);
      sigY = DT * sin(phi);
      newX += sigX;
      newY += sigY;
      DL_time = DL / driftVelocity_gas;
      electronstream[i] += DL_time;
      if (i >= Nee && eTrain) {  // all numbers are based on arXiv:1711.07025
        E_liq = fdetector->get_E_gas() / (EPS_LIQ / EPS_GAS);
        tau2 = 0.58313 * exp(0.20929 * E_liq) * 1e3;
        tau1 = 1.40540 * exp(0.15578 * E_liq) * 1e3 * 1e-2;
        amp2 = 0.38157 * exp(0.21177 * E_liq) * 1e-2;
        if (RandomGen::rndm()->rand_uniform() < (amp2 / 0.035) * ExtEff)
          electronstream[i] -= tau2 * log(RandomGen::rndm()->rand_uniform());
        else
          electronstream[i] -= tau1 * log(RandomGen::rndm()->rand_uniform());
      }
      for (double j = 0.; j < SE; j += 1.) {
        double phe = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres());
        if (i < Nee || !eTrain) pulseArea += phe;
        origin = fdetector->get_TopDrift() +
                 gasGap * RandomGen::rndm()->rand_uniform();
        k = long(j);
        if (k >= photon_times.size()) k -= photon_times.size();
        double offset = ((fdetector->get_anode() - origin) / driftVelocity_gas +
                         electronstream[i]) *
                            1e3 +
                        photon_times[k];
        if (offset < min) min = offset;
        if (RandomGen::rndm()->rand_uniform() <
            fdetector->FitTBA(truthPos[0], truthPos[1], truthPos[2])[1]) {
          AreaTableBot[0].push_back(phe);
          AreaTableTop[0].push_back(0.0);
        } else {
          AreaTableBot[0].push_back(0.0);
          AreaTableTop[0].push_back(phe);
        }
        TimeTable.push_back(offset);
        // char line[80]; sprintf ( line, "%lu\t%.0f\t%.2f\t%i", evtNum, offset,
        // phe, (int)(i<Nee) ); pulseFile << line << endl;
      }
    }
    int numPts = Nphe * 1000;
    if (numPts < 1000 / SAMPLE_SIZE) numPts = 1000 / SAMPLE_SIZE;
    AreaTableBot[1].resize(numPts, 0.);
    AreaTableTop[1].resize(numPts, 0.);
    min -= 5. * SAMPLE_SIZE;
    for (k = 0; k < Nphe; ++k) {
      vector<double> PEperBin;
      if (AreaTableBot[0][k] == 0.0)
        PEperBin =
            fdetector->SinglePEWaveForm(AreaTableTop[0][k], TimeTable[k] - min);
      else
        PEperBin =
            fdetector->SinglePEWaveForm(AreaTableBot[0][k], TimeTable[k] - min);
      for (i = 0; i < int(PEperBin.size()) - 1; ++i) {
        double eTime = PEperBin[0] + i * SAMPLE_SIZE;
        int index = int(floor(eTime / SAMPLE_SIZE + 0.5));
        if (index < 0) index = 0;
        if (index >= numPts) index = numPts - 1;
        if (AreaTableBot[0][k] == 0.0)
          AreaTableTop[1][index] += 10. * AreaTableTop[0][k] /
                                    (PULSE_WIDTH * sqrt(2. * M_PI)) *
                                    exp(-pow(eTime - TimeTable[k] + min, 2.) /
                                        (2. * PULSE_WIDTH * PULSE_WIDTH));
        else
          AreaTableBot[1][index] += 10. * AreaTableBot[0][k] /
                                    (PULSE_WIDTH * sqrt(2. * M_PI)) *
                                    exp(-pow(eTime - TimeTable[k] + min, 2.) /
                                        (2. * PULSE_WIDTH * PULSE_WIDTH));
      }
    }
    for (k = 0; k < numPts; ++k) {
      if ((AreaTableBot[1][k] + AreaTableTop[1][k]) <= PULSEHEIGHT) continue;
      wf_time.push_back(k * SAMPLE_SIZE + long(min + SAMPLE_SIZE / 2.));
      wf_amp.push_back(AreaTableBot[1][k] + AreaTableTop[1][k]);

      if (outputTiming) {
        char line[80];
        sprintf(line, "%lu\t%ld\t%.2f\t%.2f", evtNum, wf_time.back(),
                AreaTableBot[1][k], AreaTableTop[1][k]);
        pulseFile << line << endl;
      }
    }
  } else {
    Nph = long(floor(RandomGen::rndm()->rand_gauss(
                         elYield * double(Nee), sqrt(fdetector->get_s2Fano() *
                                                     elYield * double(Nee))) +
                     0.5));
    nHits = BinomFluct(Nph, fdetector->get_g1_gas() * posDep);
    Nphe = nHits + BinomFluct(nHits, fdetector->get_P_dphe());
    pulseArea = RandomGen::rndm()->rand_gauss(
        Nphe, fdetector->get_sPEres() * sqrt(Nphe));
  }

  pulseArea = RandomGen::rndm()->rand_gauss(
      pulseArea, fdetector->get_noise()[3] * pulseArea);
  double pulseAreaC =
      pulseArea / exp(-dt / fdetector->get_eLife_us()) / posDepSm;
  double Nphd = pulseArea / (1. + fdetector->get_P_dphe());
  double NphdC = pulseAreaC / (1. + fdetector->get_P_dphe());

  double S2b = RandomGen::rndm()->rand_gauss(
      fdetector->FitTBA(truthPos[0], truthPos[1], truthPos[2])[1] * pulseArea,
      sqrt(fdetector->FitTBA(truthPos[0], truthPos[1], truthPos[2])[1] *
           pulseArea *
           (1. - fdetector->FitTBA(truthPos[0], truthPos[1], truthPos[2])[1])));
  double S2bc =
      S2b / exp(-dt / fdetector->get_eLife_us()) /
      posDepSm;  // for detectors using S2 bottom-only in their analyses

  ionization[0] = Nee;  // integer number of electrons unabsorbed in liquid then
                        // getting extracted
  ionization[1] = Nph;  // raw number of photons produced in the gas gap
  ionization[2] = nHits;  // MC-true integer hits in same OR different PMTs, NO
                          // double phe effect
  ionization[3] = Nphe;   // MC-true integer hits WITH double phe effect (Nphe >
                          // nHits). S2 has more steps than S1 (e's 1st)
  if (fdetector->get_s2_thr() >= 0) {
    ionization[4] = pulseArea;  // floating real# smeared DAQ pulse areas in
                                // phe, NO XYZ correction
    ionization[5] =
        pulseAreaC;  // smeared DAQ pulse areas in phe, WITH XYZ correction
    ionization[6] = Nphd;   // same as pulse area, adjusted/corrected *downward*
                            // for 2-PE effect (LUX phd units)
    ionization[7] = NphdC;  // same as Nphd, but XYZ-corrected
  }
  // Negative S2 threshold is a hidden feature, which switches from S2 -> S2
  // bottom (NOT literally negative)
  else {
    ionization[4] = S2b;   // floating real# smeared pulse areas in phe ONLY
                           // including bottom PMTs, NO XYZ correction
    ionization[5] = S2bc;  // floating real# smeared pulse areas in phe ONLY
                           // including bottom PMTs, WITH XYZ correction
    ionization[6] =
        S2b / (1. + fdetector->get_P_dphe());  // same as S2b, but adjusted for
                                               // 2-PE effect (LUX phd units)
    ionization[7] = S2bc / (1. + fdetector->get_P_dphe());  // same as S2bc, but
                                                            // adjusted for 2-PE
                                                            // effect (LUX phd
                                                            // units)
  }

  if (pulseArea < fabs(fdetector->get_s2_thr()))
    for (i = 0; i < 8; i++) {
      if (ionization[i] == 0.) ionization[i] = PHE_MIN;
      ionization[i] *= -1.;
    }

  ionization[8] =
      g2;  // g2 = ExtEff * SE, gain of EL in gas gap (from CalculateG2)

  return ionization;
}

vector<double> NESTcalc::CalculateG2(bool verbosity) {
  vector<double> g2_params(5);

  // Set parameters for calculating EL yield and extraction
  double alpha = 0.137, beta = 177.,
         gamma = 45.7;  // note the value of alpha is similar to ~1/7eV. Not
                        // coincidence. Noted in Mock et al.
  double epsilonRatio = EPS_LIQ / EPS_GAS;
  if (fdetector->get_inGas())
    epsilonRatio = 1.;  // in an all-gas detector, E_liq variable below simply
                        // becomes the field value between anode and gate

  // Convert gas extraction field to liquid field
  double E_liq = fdetector->get_E_gas() / epsilonRatio;  // kV per cm
  double ExtEff = -0.03754 * pow(E_liq, 2.) + 0.52660 * E_liq -
                  0.84645;  // arXiv:1710.11032 (PIXeY)
  if (ExtEff > 1. || fdetector->get_inGas() || E_liq > 7.) ExtEff = 1.;
  if (ExtEff < 0. || E_liq <= 0.) ExtEff = 0.;

  double gasGap =
      fdetector->get_anode() -
      fdetector
          ->get_TopDrift();  // EL gap in mm -> cm, affecting S2 size linearly
  if (gasGap <= 0. && E_liq > 0.) {
    cerr << "\tERR: The gas gap in the S2 calculation broke!!!!" << endl;
    exit(1);
  }

  // Calculate EL yield based on gas gap, extraction field, and pressure
  double elYield = (alpha * fdetector->get_E_gas() * 1e3 -
                    beta * fdetector->get_p_bar() - gamma) *
                   gasGap * 0.1;  // arXiv:1207.2292 (HA, Vitaly C.)
  if (elYield <= 0.0 && E_liq != 0.) {
    cerr << "\tWARNING, the field in gas must be at least "
         << 1e-3 * (beta * fdetector->get_p_bar() + gamma) / alpha
         << " kV/cm, for S2 to work," << endl;
    cerr << "\tOR: your pressure for gas must be less than "
         << (alpha * fdetector->get_E_gas() * 1e3 - gamma) / beta << " bar."
         << endl;
  }
  // Calculate single electron size and then g2
  double SE = elYield * fdetector->get_g1_gas();  // multiplying by light
                                                  // collection efficiency in
                                                  // the gas gap
  if (fdetector->get_s2_thr() < 0)
    SE *= fdetector->FitTBA(0., 0., fdetector->get_TopDrift() / 2.)[1];
  double g2 = ExtEff * SE;
  double StdDev = sqrt((1. - fdetector->get_g1_gas()) * SE +
                       fdetector->get_s2Fano() * fdetector->get_s2Fano() +
                       fdetector->get_sPEres());

  if (verbosity) {
    cout << endl
         << "g1 = " << fdetector->get_g1() << " phd per photon\tg2 = " << g2
         << " phd per electron (e-EE = ";
    cout << ExtEff * 100. << "%, SE_mean,width = " << SE << "," << StdDev
         << ")\t";
  }

  // Store the g2 parameters in a vector for later (calculated once per
  // detector)
  g2_params[0] = elYield;
  g2_params[1] = ExtEff;
  g2_params[2] = SE;
  g2_params[3] = g2;
  g2_params[4] = gasGap;

  return g2_params;
}

long double NESTcalc::Factorial(double x) { return tgammal(x + 1.); }

double NESTcalc::nCr(double n, double r) {
  return Factorial(n) / (Factorial(r) * Factorial(n - r));
}

vector<double> NESTcalc::GetSpike(int Nph, double dx, double dy, double dz,
                                  double driftSpeed, double dS_mid,
                                  vector<double> oldScint) {
  vector<double> newSpike(2);

  if (oldScint[7] > SPIKES_MAXM) {
    newSpike[0] = oldScint[4];
    newSpike[1] = oldScint[5];
    return newSpike;
  }
  newSpike[0] = fabs(oldScint[6]);
  newSpike[0] = RandomGen::rndm()->rand_gauss(
      newSpike[0], (fdetector->get_sPEres() / 4.) * sqrt(newSpike[0]));
  if (newSpike[0] < 0.0) newSpike[0] = 0.0;
  newSpike[1] = newSpike[0] / fdetector->FitS1(dx, dy, dz) *
                fdetector->FitS1(0., 0., fdetector->get_TopDrift() -
                                             dS_mid * fdetector->get_dtCntr());

  return newSpike;  // regular and position-corrected spike counts returned
}

double NESTcalc::SetDensity(double Kelvin,
                            double bara) {  // currently only for fixed pressure
                                            // (saturated vapor pressure); will
                                            // add pressure dependence later

  if (Kelvin < 161.40) {  // solid Xenon
    cerr << "\nWARNING: SOLID PHASE. IS THAT WHAT YOU WANTED?\n";
    return 3.41;  // from Yoo at 157K
    // other sources say 3.100 (Wikipedia, 'maximum') and 3.64g/mL at an unknown
    // temperature
  }

  double VaporP_bar;  // we will calculate using NIST
  if (Kelvin < 289.7)
    VaporP_bar = pow(10., 4.0519 - 667.16 / Kelvin);
  else
    VaporP_bar = DBL_MAX;
  if (bara < VaporP_bar) {
    double density =
        bara * 1e5 / (Kelvin * 8.314);  // ideal gas law approximation, mol/m^3
    density *= MOLAR_MASS * 1e-6;
    cerr << "\nWARNING: GAS PHASE. IS THAT WHAT YOU WANTED?\n";
    fdetector->set_inGas(true);
    return density;  // in g/cm^3
  }

  return 2.9970938084691329E+02 * exp(-8.2598864714323525E-02 * Kelvin) -
         1.8801286589442915E+06 * exp(-pow((Kelvin - 4.0820251276172212E+02) /
                                               2.7863170223154846E+01,
                                           2.)) -
         5.4964506351743057E+03 * exp(-pow((Kelvin - 6.3688597345042672E+02) /
                                               1.1225818853661815E+02,
                                           2.)) +
         8.3450538370682614E+02 * exp(-pow((Kelvin + 4.8840568924597342E+01) /
                                               7.3804147172071107E+03,
                                           2.)) -
         8.3086310405942265E+02;  // in grams per cubic centimeter based on
                                  // zunzun fit to NIST data; will add gas later
}

double NESTcalc::SetDriftVelocity(double Kelvin, double Density,
                                  double eField) {  // for liquid and solid only

  if (fdetector->get_inGas()) return SetDriftVelocity_MagBoltz(Density, eField);

  double speed =
      0.0;  // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
  int i, j;
  double vi, vf, slope, Ti, Tf, offset;

  double polyExp[11][7] = {
      {-3.1046, 27.037, -2.1668, 193.27, -4.8024, 646.04, 9.2471},  // 100K
      {-2.7394, 22.760, -1.7775, 222.72, -5.0836, 724.98, 8.7189},  // 120
      {-2.3646, 164.91, -1.6984, 21.473, -4.4752, 1202.2, 7.9744},  // 140
      {-1.8097, 235.65, -1.7621, 36.855, -3.5925, 1356.2, 6.7865},  // 155
      {-1.5000, 37.021, -1.1430, 6.4590, -4.0337, 855.43,
       5.4238},  // 157, merging Miller with Yoo
      {-1.4939, 47.879, 0.12608, 8.9095, -1.3480, 1310.9,
       2.7598},  // 163, merging Miller with Yoo
      {-1.5389, 26.602, -.44589, 196.08, -1.1516, 1810.8, 2.8912},  // 165
      {-1.5000, 28.510, -.21948, 183.49, -1.4320, 1652.9, 2.884},   // 167
      {-1.1781, 49.072, -1.3008, 3438.4, -.14817, 312.12, 2.8049},  // 184
      {1.2466, 85.975, -.88005, 918.57, -3.0085, 27.568, 2.3823},   // 200
      {334.60, 37.556, 0.92211, 345.27, -338.00, 37.346, 1.9834}};  // 230

  double Temperatures[11] = {100., 120., 140., 155., 157., 163.,
                             165., 167., 184., 200., 230.};

  if (Kelvin >= Temperatures[0] && Kelvin < Temperatures[1])
    i = 0;
  else if (Kelvin >= Temperatures[1] && Kelvin < Temperatures[2])
    i = 1;
  else if (Kelvin >= Temperatures[2] && Kelvin < Temperatures[3])
    i = 2;
  else if (Kelvin >= Temperatures[3] && Kelvin < Temperatures[4])
    i = 3;
  else if (Kelvin >= Temperatures[4] && Kelvin < Temperatures[5])
    i = 4;
  else if (Kelvin >= Temperatures[5] && Kelvin < Temperatures[6])
    i = 5;
  else if (Kelvin >= Temperatures[6] && Kelvin < Temperatures[7])
    i = 6;
  else if (Kelvin >= Temperatures[7] && Kelvin < Temperatures[8])
    i = 7;
  else if (Kelvin >= Temperatures[8] && Kelvin < Temperatures[9])
    i = 8;
  else if (Kelvin >= Temperatures[9] && Kelvin <= Temperatures[10])
    i = 9;
  else {
    cerr << "\nERROR: TEMPERATURE OUT OF RANGE (100-230 K)\n";
    exit(1);
  }

  j = i + 1;
  Ti = Temperatures[i];
  Tf = Temperatures[j];
  // functional form from http://zunzun.com
  vi = polyExp[i][0] * exp(-eField / polyExp[i][1]) +
       polyExp[i][2] * exp(-eField / polyExp[i][3]) +
       polyExp[i][4] * exp(-eField / polyExp[i][5]) + polyExp[i][6];
  vf = polyExp[j][0] * exp(-eField / polyExp[j][1]) +
       polyExp[j][2] * exp(-eField / polyExp[j][3]) +
       polyExp[j][4] * exp(-eField / polyExp[j][5]) + polyExp[j][6];
  if (Kelvin == Ti) return vi;
  if (Kelvin == Tf) return vf;
  if (vf < vi) {
    offset = (sqrt((Tf * (vf - vi) - Ti * (vf - vi) - 4.) * (vf - vi)) +
              sqrt(Tf - Ti) * (vf + vi)) /
             (2. * sqrt(Tf - Ti));
    slope = -(sqrt(Tf - Ti) *
                  sqrt((Tf * (vf - vi) - Ti * (vf - vi) - 4.) * (vf - vi)) -
              (Tf + Ti) * (vf - vi)) /
            (2. * (vf - vi));
    speed = 1. / (Kelvin - slope) + offset;
  } else {
    slope = (vf - vi) / (Tf - Ti);
    speed = slope * (Kelvin - Ti) + vi;
  }

  if (speed <= 0.) {
    if (eField < 1e2 && eField >= FIELD_MIN) {
      cerr << "\nERROR: DRIFT SPEED NON-POSITIVE -- FIELD TOO LOW\n";
      exit(1);
    }
    if (eField > 1e4) {
      cerr << "\nERROR: DRIFT SPEED NON-POSITIVE -- FIELD TOO HIGH\n";
      exit(1);
    }
  }
  return speed;
}

double NESTcalc::SetDriftVelocity_MagBoltz(
    double density, double efieldinput)  // Nichole Barry UCD 2011
{
  density *= NEST_AVO / MOLAR_MASS;
  // Gas equation one coefficients (E/N of 1.2E-19 to 3.5E-19)
  double gas1a = 395.50266631436, gas1b = -357384143.004642,
         gas1c = 0.518110447340587;
  // Gas equation two coefficients (E/N of 3.5E-19 to 3.8E-17)
  double gas2a = -592981.611357632, gas2b = -90261.9643716643,
         gas2c = -4911.83213989609, gas2d = -115.157545835228,
         gas2f = -0.990440443390298, gas2g = 1008.30998933704,
         gas2h = 223.711221224885;
  double edrift = 0., gasdep = efieldinput / density, gas1fix = 0.,
         gas2fix = 0.;

  if (gasdep < 1.2e-19 && gasdep >= 0.) edrift = 4e22 * gasdep;
  if (gasdep < 3.5e-19 && gasdep >= 1.2e-19) {
    gas1fix = gas1b * pow(gasdep, gas1c);
    edrift = gas1a * pow(gasdep, gas1fix);
  }
  if (gasdep < 3.8e-17 && gasdep >= 3.5e-19) {
    gas2fix = log(gas2g * gasdep);
    edrift = (gas2a + gas2b * gas2fix + gas2c * pow(gas2fix, 2.) +
              gas2d * pow(gas2fix, 3.) + gas2f * pow(gas2fix, 4.)) *
             (gas2h * exp(gasdep));
  }
  if (gasdep >= 3.8e-17) edrift = 6e21 * gasdep - 32279.;

  return edrift * 1e-5;  // from cm/s into mm per microsecond
}

vector<double> NESTcalc::SetDriftVelocity_NonUniform(double rho, double zStep,
                                                     double dx, double dy) {
  vector<double> speedTable;
  double driftTime, zz;

  for (double pos_z = 0.0; pos_z < fdetector->get_TopDrift(); pos_z += zStep) {
    driftTime = 0.0;
    for (zz = pos_z; zz < fdetector->get_TopDrift(); zz += zStep) {
      if (pos_z > fdetector->get_gate()) {
        if (!fdetector->get_inGas())
          driftTime += zStep / SetDriftVelocity(fdetector->get_T_Kelvin(), rho,
                                                fdetector->get_E_gas() /
                                                    (EPS_LIQ / EPS_GAS) * 1e3);
        else  // if gate == TopDrift properly set, shouldn't happen
          driftTime += zStep / SetDriftVelocity_MagBoltz(
                                   rho, fdetector->get_E_gas() * 1e3);
      } else
        driftTime +=
            zStep /
            SetDriftVelocity(
                fdetector->get_T_Kelvin(), rho,
                fdetector->FitEF(dx, dy,
                                 zz));  // update x and y if you want 3-D fields
    }

    speedTable.push_back((zz - pos_z) / driftTime);  // uses highest zz
  }

  return speedTable;
}

vector<double> NESTcalc::xyResolution(double xPos_mm, double yPos_mm,
                                      double A_top) {
  vector<double> xySmeared(2);
  A_top *=
      1. -
      fdetector->FitTBA(xPos_mm, yPos_mm, fdetector->get_TopDrift() / 2.)[1];

  double rad = sqrt(pow(xPos_mm, 2.) + pow(yPos_mm, 2.));
  double kappa = fdetector->get_PosResBase() +
                 exp(fdetector->get_PosResExp() * rad);  // arXiv:1710.02752
  double sigmaR = kappa / sqrt(A_top);                   // ibid.

  double phi = 2. * M_PI * RandomGen::rndm()->rand_uniform();
  sigmaR = fabs(RandomGen::rndm()->rand_gauss(0.0, sigmaR));
  double sigmaX = sigmaR * cos(phi);
  double sigmaY = sigmaR * sin(phi);

  xySmeared[0] = xPos_mm + sigmaX;
  xySmeared[1] = yPos_mm + sigmaY;

  return xySmeared;  // new X and Y position in mm with empirical smearing. LUX
                     // Run03 example
}

double NESTcalc::PhotonEnergy(bool s2Flag, bool state, double tempK) {
  double wavelength, E_keV;  // wavelength is in nanometers

  if (!state)  // liquid or solid
    wavelength = RandomGen::rndm()->rand_gauss(
        178., 14. / 2.355);  // taken from Jortner JchPh 42 '65
  else                       // gas
    wavelength =
        RandomGen::rndm()->rand_gauss(175., 5.);  // G4S1Light, probably Doke

  if (s2Flag) {  // S2 different from ordinary gas (or just measurement error?)
    if (tempK < 200.)  // cold gas
      wavelength = RandomGen::rndm()->rand_gauss(
          179., 5.);  // source is G4S2Light.cc from the old NEST
    else
      wavelength = RandomGen::rndm()->rand_gauss(174., 5.);  // ditto
  }

  E_keV = 1240e-3 / wavelength;  // h*c in keV-nm divided by lambda in nm
  if (E_keV > W_SCINT)
    E_keV = W_SCINT;  // don't go so high breaks G4. Gauss in lambda -> non-G in
                      // E, high tail

  return E_keV * 1000.;  // convert from keV into eV. Eventually add full T, P
                         // dependence
}

double NESTcalc::CalcElectronLET(double E) {
  double LET;

  // use a spline fit to online ESTAR data
  if (E >= 1.)
    LET = 58.482 - 61.183 * log10(E) + 19.749 * pow(log10(E), 2) +
          2.3101 * pow(log10(E), 3) - 3.3469 * pow(log10(E), 4) +
          0.96788 * pow(log10(E), 5) - 0.12619 * pow(log10(E), 6) +
          0.0065108 * pow(log10(E), 7);
  // at energies <1 keV, use a different spline, determined manually by
  // generating sub-keV electrons in Geant4 and looking at their ranges, since
  // ESTAR does not go this low (4.9.4)
  else if (E > 0. && E < 1.)
    LET = 6.9463 + 815.98 * E - 4828 * pow(E, 2) + 17079 * pow(E, 3) -
          36394 * pow(E, 4) + 44553 * pow(E, 5) - 28659 * pow(E, 6) +
          7483.8 * pow(E, 7);
  else
    LET = 0.;

  return LET;
}
