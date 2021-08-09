
#include <exception>
#include "NEST.hh"
#include <stdexcept>

#define InfraredER 1.1716263232
//#define InfraredNR 7.

#define ChargeLoss 0.20

using namespace std;
using namespace NEST;

static constexpr int podLength = 1100; //roughly 100-1,000 ns for S1

bool kr83m_reported_low_deltaT = false; //to aid in verbosity 

const std::vector<double> NESTcalc::default_NuisParam = {11., 1.1, 0.0480, -0.0533, 12.6, 0.3, 2., 0.3, 2., 0.5, 1.,
                                                         1.};
const std::vector<double> NESTcalc::default_FreeParam = {1., 1., 0.1, 0.5, 0.19, 2.25};

int64_t NESTcalc::BinomFluct(int64_t N0, double prob) {
    double mean = N0 * prob;
    double sigma = sqrt(N0 * prob * (1. - prob));
    int N1 = 0;

    if (prob <= 0.00) return N1;
    if (prob >= 1.00) return N0;

    if (N0 <= 9. * (1. - prob) / prob || N0 <= 9. * prob / (1. - prob)) {
        //https://en.wikipedia.org/wiki/Binomial_distribution#Normal_approximation
        for (int i = 0; i < N0; ++i) {
            if (RandomGen::rndm()->rand_uniform() < prob) ++N1;
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
                                     const std::vector<double> &NuisParam /*={11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/,
                                     const std::vector<double> &FreeParam /*={1.,1.,0.1,0.5,0.19,2.25}*/,
                                     bool do_times /*=true*/) {

    if (density < 1.) fdetector->set_inGas(true);
    NESTresult result;
    result.yields = GetYields(species, energy, density, dfield, A, Z, NuisParam);
    result.quanta = GetQuanta(result.yields, density, FreeParam);
    if (do_times)
        result.photon_times = GetPhotonTimes(
                species, result.quanta.photons, result.quanta.excitons, dfield, energy);
    else {
        result.photon_times = photonstream(result.quanta.photons, 0.0);
    }
    return result;

}

double NESTcalc::PhotonTime(INTERACTION_TYPE species, bool exciton,
                            double dfield, double energy) {
    double time_ns = 0., SingTripRatio, tauR = 0., tau3 = 23.97,
            tau1 = 3.27; //arXiv:1802.06162. NR may need tauR ~0.5-1ns instead of 0
    if (fdetector->get_inGas() ||
        energy < W_DEFAULT * 0.001) {  
	// from G4S1Light.cc in old NEST
        tau1 = 5.18;                     // uncertainty of 1.55 ns from G4S2Light
        tau3 = 100.1;                    // uncertainty of 7.90 ns from G4S2Light
    }
    // tau1 = 3.5*ns; tau3 = 20.*ns; tauR = 40.*ns for solid Xe from old NEST.
    // Here assuming same as in liquid

    // if ( dfield > 60. ) dfield = 60. // makes Xed work. 200 for LUX Run04 instead. A mystery! Why no field dep?

    double LET = CalcElectronLET(energy, ATOM_NUM);
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) { 
	//copied from 2013 NEST version for LAr on LBNE
        tau1 = RandomGen::rndm()->rand_gauss(6.5, 0.8); // error from weighted average
        tau3 = RandomGen::rndm()->rand_gauss(1300, 50); // ibid.
        tauR = RandomGen::rndm()->rand_gauss(0.8, 0.2); // Kubota 1979
        if (species <= Cf) SingTripRatio = 0.22218 * pow(energy, 0.48211);
        else if (species == ion) { //really only alphas here
            SingTripRatio = (-0.065492 + 1.9996 * exp(-energy / 1e3)) / (1. + 0.082154 / pow(energy / 1e3, 2.)) +
                            2.1811; //uses energy in MeV not keV
	}
        else {
            SingTripRatio = 0.2701 + 0.003379 * LET - 4.7338e-5 * pow(LET, 2.) + 8.1449e-6 * pow(LET, 3.);
            if (LET < 3. && !exciton) SingTripRatio = RandomGen::rndm()->rand_gauss(0.5, 0.2);
            if (LET < 3. && exciton) SingTripRatio = RandomGen::rndm()->rand_gauss(0.36, 0.06);
        } //lastly is ER for LAr
    } else {
        if (species <= Cf) // NR
            SingTripRatio =
                    (0.21 - 0.0001 * dfield) * pow(energy, 0.21 - 0.0001 * dfield); //arXiv:1803.07935. LUX:0.15*E^0.15
        else if (species == ion) { 
            // e.g., alphas
            SingTripRatio = 0.065 * pow(energy, 0.416); // spans 2.3 (alpha) and 7.8 (Cf in Xe) from NEST v1
	}
        else { 
            // ER
            if (!exciton) {
                if (energy > 1e3) energy = 1e3; // MIP above ~1 MeV. Fix thanks to Austin de St. Croix
                tauR = exp(-0.00900 * dfield) * (7.3138 + 3.8431 * log10(energy)); // arXiv:1310.1117
                if (tauR < 3.5) tauR = 3.5; //used to be for gammas only but helpful for matching beta data better
                if (dfield > 8e2) dfield = 8e2; //to match Kubota's 4,000 V/cm
                SingTripRatio = 1.00 * pow(energy, -0.45 + 0.0005 *
                                                           dfield); // see comment below; also, dfield may need to be fixed at ~100-200 V/cm (for NR too)
            } else {
                SingTripRatio = 0.20 * pow(energy, -0.45 + 0.0005 * dfield); // mixing arXiv:1807.07121 with Kubota 1979
	    }
        }
    }
    if (fdetector->get_inGas() || energy < W_DEFAULT * 0.001) {
        SingTripRatio = 0.1;
        if (fdetector->get_inGas() && !exciton)
            tauR = 28e3;
        else {
	  tauR = 0.; //28 microseconds comes from Henrique: https://doi.org/10.1016/j.astropartphys.2018.04.006
	}
        if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) {
            tau3 = 1600.;
            tau1 = 6.;
        } // from old G4S2Light
    }
    if (tauR < 0.) tauR = 0.; //in case varied with Gaussian earlier

    // the recombination time is non-exponential, but approximates
    // to exp at long timescales (see Kubota '79)
    time_ns += tauR * (1.0 / RandomGen::rndm()->rand_uniform() - 1.);

    if (RandomGen::rndm()->rand_uniform() < SingTripRatio / (1. + SingTripRatio))
        time_ns -= tau1 * log(RandomGen::rndm()->rand_uniform());
    else {
        time_ns -= tau3 * log(RandomGen::rndm()->rand_uniform());
    }

    return time_ns;
}

// in analytical model, an empirical distribution (Brian Lenardo and Dev A.K.)
photonstream NESTcalc::AddPhotonTransportTime(const photonstream &emitted_times,
                                              double x, double y, double z) {
    photonstream return_photons;
    for (auto t : emitted_times) {
        double newtime = t + fdetector->OptTrans(x, y, z);
        return_photons.emplace_back(newtime);
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
        return_photons.emplace_back(PhotonTime(species, isExciton, dfield, energy));
    }

    return return_photons;
}

double NESTcalc::RecombOmegaNR(double elecFrac, const std::vector<double> &FreeParam/*={1.,1.,0.1,0.5,0.19,2.25}*/) {
    double omega = FreeParam[2] * exp(-0.5 * pow(elecFrac - FreeParam[3], 2.) / (FreeParam[4] * FreeParam[4]));
    if (omega < 0.)
        omega = 0;
    return omega;
}

double NESTcalc::RecombOmegaER(double efield, double elecFrac, const std::vector<double> &FreeParam) {
    double ampl = 0.14 + (0.043 - 0.14) / (1. + pow(efield / 1210.,
                                                    1.25)); //0.086036+(0.0553-0.086036)/pow(1.+pow(efield/295.2,251.6),0.0069114); //pair with GregR mean yields model
    if (ampl < 0.)
        ampl = 0.;
    double wide = 0.205; //or: FreeParam #2, like amplitude (#1)
    double cntr = 0.5; //0.41-45 agrees better with Dahl thesis. Odd! Reduces fluctuations for high e-Frac (high EF,low E). Also works with GregR LUX Run04 model. FreeParam #3
    //for gamma-rays larger than 100 keV at least in XENON10 use 0.43 as the best fit. 0.62-0.37 for LUX Run03
    double skew = -0.2; //FreeParam #4
    double mode = cntr + 2.*inv_sqrt2_PI * skew * wide / sqrt(1. + skew * skew);
    double norm = 1. / (exp(-0.5 * pow(mode - cntr, 2.) / (wide * wide)) *
                        (1. + erf(skew * (mode - cntr) / (wide * sqrt2)))); //makes sure omega never exceeds ampl
    double omega = norm * ampl * exp(-0.5 * pow(elecFrac - cntr, 2.) / (wide * wide)) *
                   (1. + erf(skew * (elecFrac - cntr) / (wide * sqrt2)));
    if (omega < 0.)
        omega = 0;
    return omega;
}

double NESTcalc::FanoER(double density, double Nq_mean, double efield) {
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) // liquid argon
        return 0.1115; // T&I. conflicting reports of .107 (Doke) & ~0.1 elsewhere
    double Fano = 0.12707 - 0.029623 * density -  // Fano factor is  << 1
                  0.0057042 *
                  pow(density,
                      2.) +  //~0.1 for GXe w/ formula from Bolotnikov et al. 1995
                  0.0015957 *
                  pow(density,
                      3.);  // to get it to be ~0.03 for LXe (E Dahl Ph.D. thesis)
    if (!fdetector->get_inGas())
        Fano += 0.0015 * sqrt(Nq_mean) * pow(efield, 0.5);
    return Fano;
}


QuantaResult NESTcalc::GetQuanta(const YieldResult &yields, double density,
                                 const std::vector<double> &FreeParam/*={1.,1.,0.1,0.5,0.19,2.25}*/) {
    QuantaResult result{};
    bool HighE;
    int Nq_actual, Ne, Nph, Ni, Nex;

    if (FreeParam.size() < 6) {
        throw std::runtime_error("ERROR: You need a minimum of 6 free parameters for the resolution model.");
    }

    double excitonRatio = yields.ExcitonRatio;
    double Nq_mean = yields.PhotonYield + yields.ElectronYield;

    double elecFrac = max(0., min(yields.ElectronYield / Nq_mean, 1.));

    if (excitonRatio < 0.) {
        excitonRatio = 0.;
        HighE = true;
    } else {
        HighE = false;
    }

    double alf = 1. / (1. + excitonRatio);
    double recombProb = 1. - (excitonRatio + 1.) * elecFrac;
    if (recombProb < 0.) {
        excitonRatio = 1. / elecFrac - 1.;
    }

    if (ValidityTests::nearlyEqual(yields.Lindhard, 1.)) {
        double Fano = FanoER(density, Nq_mean, yields.ElectricField);
        Nq_actual = int(floor(
                RandomGen::rndm()->rand_gauss(Nq_mean, sqrt(Fano * Nq_mean)) + 0.5));
        if (Nq_actual < 0 || ValidityTests::nearlyEqual(Nq_mean, 0.)) Nq_actual = 0;

        Ni = BinomFluct(Nq_actual, alf);
        Nex = Nq_actual - Ni;

    } else {
        double Fano = FreeParam[0];
        Ni = int(floor(RandomGen::rndm()->rand_gauss(Nq_mean * alf,
                                                     sqrt(Fano * Nq_mean * alf)) +
                       0.5));
        if (Ni < 0) Ni = 0;
        Fano = FreeParam[1];
        Nex = int(
                floor(RandomGen::rndm()->rand_gauss(
                        Nq_mean * excitonRatio * alf, sqrt(Fano * Nq_mean * excitonRatio * alf)) +
                      0.5));
        if (Nex < 0) Nex = 0;
        Nq_actual = Nex + Ni;
    }

    if (Nq_actual == 0) {
        result.ions = 0;
        result.excitons = 0;
        result.photons = 0;
        result.electrons = 0;
        result.Variance = 0;
        result.recombProb = 0;
        return result;
    }

    if (Nex < 0) Nex = 0;
    if (Ni < 0) Ni = 0;
    if (Nex > Nq_actual) Nex = Nq_actual;
    if (Ni > Nq_actual) Ni = Nq_actual;

    result.ions = Ni;
    result.excitons = Nex;

    if (Nex <= 0 && HighE)
        recombProb = yields.PhotonYield / double(Ni);
    recombProb = max(0., min(recombProb, 1.));
    if (std::isnan(recombProb) || std::isnan(elecFrac) || Ni == 0 || ValidityTests::nearlyEqual(recombProb, 0.0)) {
        result.photons = Nex;
        result.electrons = Ni;
        elecFrac = 1.0;
        result.recombProb = 0.;
        result.Variance = 0.;
        return result;
    }

    //set omega (non-binomial recombination fluctuations parameter) according to whether the Lindhard <1, i.e. this is NR.
    double omega =
            yields.Lindhard < 1 ? RecombOmegaNR(elecFrac, FreeParam) : RecombOmegaER(yields.ElectricField, elecFrac,
                                                                                     FreeParam);
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) omega = 0.0; // Ar has no non-binom sauce
    double Variance =
            recombProb * (1. - recombProb) * Ni + omega * omega * Ni * Ni;

    double skewness;
    if ((yields.PhotonYield + yields.ElectronYield) > 1e4 || yields.ElectricField > 4e3 || yields.ElectricField < 50.) {
        skewness = 0.00; //make it a constant 0 when outside the range of Vetri Velan's Run04 models.
    } else { 
	// LUX Skewness Model
        Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
        double Wq_eV = wvalue.Wq_eV;
        double engy = 1e-3 * Wq_eV * (yields.PhotonYield + yields.ElectronYield);
        double fld = yields.ElectricField;

        double alpha0 = 1.39;
        double cc0 = 4.0, cc1 = 22.1;
        double E0 = 7.7, E1 = 54., E2 = 26.7, E3 = 6.4;
        double F0 = 225., F1 = 71.;

        skewness = 0.;

        if (ValidityTests::nearlyEqual(yields.Lindhard, 1.)) {
            skewness = 1. / (1. + exp((engy - E2) / E3)) *
                       (alpha0 + cc0 * exp(-1. * fld / F0) * (1. - exp(-1. * engy / E0))) +
                       1. / (1. + exp(-1. * (engy - E2) / E3)) * cc1 * exp(-1. * engy / E1) *
                       exp(-1. * sqrt(fld) / sqrt(F1));
            //if ( std::abs(skewness) <= DBL_MIN ) skewness = DBL_MIN;
        } else {
            skewness = FreeParam[5]; //2.25 but ~5-20 also good (for NR). All better than zero, but 0 is OK too
        } //note to self: find way to make 0 for ion (wall BG) incl. alphas?
    }

    double widthCorrection = sqrt(1. - (2. / M_PI) * skewness * skewness / (1. + skewness * skewness));
    double muCorrection =
            (sqrt(Variance) / widthCorrection) * (skewness / sqrt(1. + skewness * skewness)) * 2.*inv_sqrt2_PI;
    if (std::abs(skewness) > DBL_MIN && ValidityTests::nearlyEqual(ATOM_NUM, 54.)) //skewness model only for Xenon!
        Ne = int(floor(RandomGen::rndm()->rand_skewGauss((1. - recombProb) * Ni - muCorrection,
                                                         sqrt(Variance) / widthCorrection, skewness) + 0.5));
    else {
        Ne = int(floor(RandomGen::rndm()->rand_gauss((1. - recombProb) * Ni, sqrt(Variance)) + 0.5));
    }

    if (Ne < 0) Ne = 0;
    if (Ne > Ni) Ne = Ni;

    Nph = Nq_actual - Ne;
    if (Nph > Nq_actual) Nph = Nq_actual;
    if (Nph < Nex) Nph = Nex;

    if ((Nph + Ne) != (Nex + Ni)) {
        throw std::runtime_error("ERROR: Quanta not conserved. Tell Matthew Immediately!");
    }

    result.Variance = Variance;
    result.recombProb = recombProb;
    result.photons = Nph;
    result.electrons = Ne;

    return result;  // quanta returned with recomb fluctuations
}

YieldResult NESTcalc::GetYieldGamma(double energy, double density, double dfield) {
    Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
    double Wq_eV = wvalue.Wq_eV;
    double alpha = wvalue.alpha;
    constexpr double m3 = 2., m4 = 2., m6 = 0.;

    const double m1 =
            33.951 + (3.3284 - 33.951) / (1. + pow(dfield / 165.34, .72665));
    double m2 = 1000 / Wq_eV;
    double m5 =
            23.156 + (10.737 - 23.156) / (1. + pow(dfield / 34.195, .87459));
    double densCorr = 240720. / pow(density, 8.2076);
    double m7 =
            66.825 + (829.25 - 66.825) / (1. + pow(dfield / densCorr, .83344));

    double Nq = energy * 1000. / Wq_eV;
    double m8 = 2.;
    if (fdetector->get_inGas()) m8 = -2.;
    double Qy = m1 + (m2 - m1) / (1. + pow(energy / m3, m4)) + m5 +
                (m6 - m5) / (1. + pow(energy / m7, m8));
    double Ly = Nq / energy - Qy;

    YieldResult result{};
    result.PhotonYield = Ly * energy;
    result.ElectronYield = Qy * energy;
    result.ExcitonRatio = NexONi(energy, density);
    result.Lindhard = 1;
    result.ElectricField = dfield;
    result.DeltaT_Scint = -999;
    return YieldResultValidity(result, energy, Wq_eV);
}

YieldResult NESTcalc::GetYieldNROld(double energy, int option) { // possible anti-correlation in NR ignored totally

    double Ne, Nph, FakeField;

    if (option == 0) { 
	// with old-school L_eff*S_nr conversion, and more explicit Thomas-Imel box model formulae
	// but approximation where Qy can sail off as energy to 0
        FakeField = 1.00;
        Nph = 0.95 * 64. * energy * (0.050295 * pow(energy, 1.3483) * (log(1. + (0.84074 * pow(energy, 1.3875))) /
                                                                       (0.84074 * pow(energy,
                                                                                      1.3875)))); // NESTv0.97beta (2011) ~0 V/cm arXiv:1106.1613
        Ne = energy * (10.661 * pow(energy, 0.16199) *
                       (log(1. + (0.745330 * pow(energy, 1.0880))) / (0.745330 * pow(energy, 1.0880))) +
                       0.93739); // essentially Lindhard with k of 0.166
    } else if (option == 1) {
        FakeField = 200.;
        Ne = energy * (4.1395 * pow(energy, 0.13816) * (log(1. + (0.040945 * pow(energy, 1.1388))) / (0.040945 *
                                                                                                      pow(energy,
                                                                                                          1.1388)))); // conservative: k < or ~ to 0.110 (Hitachi-like)
        Nph = energy * 3.35 * pow(energy,
                                  0.29222); // NESTv0.98 (2013) model for legacy comparisons (arXiv:1307.6601) using simple power law (optional constant; optional for Nph)
    } else if (option == 2) {
        FakeField = 730.;
        Nph = 7582.3 + (-1.3728 - 7582.3) / (1. + pow(energy / 385.46,
                                                      1.2669)); // NEST v1.0 (2015) arXiv:1412.4417, with ability to cut off yield at either low or high energy. Fit to Nph
        Ne = (60.914 * pow(energy, 0.32220)) * (1. - exp(-0.12684 * pow(energy,
                                                                        0.65729))); // middle of road: k = 0.14. By Brian Lenardo, still pre-LUX-DD. With explicit bi-Exc quenching
    } else { 
	// simplest by far, but not 1st principles: based on work of M. Wyman UAlbany
        FakeField = 180.; // fit to LUX Run03 D-D alone
        Ne = (-3.8780 + 12.372 * pow(energy, 0.73502)) *
             exp(-0.0034329 * energy); // slide 68 of Matthew's private Google Slides mega-deck
        Nph = 0.81704 + 3.8584 * pow(energy, 1.30180); // ditto
    }

    YieldResult result{};
    if (Nph < 0.) Nph = 0.;
    if (Ne < 0.) Ne = 0.;
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = 1.;
    result.Lindhard = ((Nph + Ne) / energy) * W_DEFAULT * 1e-3;
    result.ElectricField = FakeField;
    result.DeltaT_Scint = -999.;
    return YieldResultValidity(result, energy, W_DEFAULT);

}

YieldResult NESTcalc::GetYieldNR(double energy, double density, double dfield, double massNum,
                                 const std::vector<double> &NuisParam/*{11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/) {

    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) massNum = fdetector->get_molarMass();

    if (NuisParam.size() < 12) {
        throw std::runtime_error("ERROR: You need a minimum of 12 nuisance parameters for the mean yields.");
    }
    if (energy > HIGH_E_NR)
        cerr << "\nWARNING: No data out here, you are beyond the AmBe endpoint of about 300 keV.\n";
    int massNumber;
    double ScaleFactor[2] = {1., 1.};
    if (ValidityTests::nearlyEqual(massNum, 0.))
        massNumber = RandomGen::rndm()->SelectRanXeAtom();
    else {
        massNumber = int(massNum);
    }
    ScaleFactor[0] = sqrt(fdetector->get_molarMass() / (double) massNumber);
    ScaleFactor[1] = ScaleFactor[0];
    double Nq = NuisParam[0] * pow(energy, NuisParam[1]);
    double ThomasImel =
            NuisParam[2] * pow(dfield, NuisParam[3]) * pow(density / DENSITY, 0.3);
    double Qy = 1. / (ThomasImel * pow(energy + NuisParam[4], NuisParam[9]));
    Qy *= 1. - 1. / pow(1. + pow((energy / NuisParam[5]), NuisParam[6]), NuisParam[10]);
    double Ly = Nq / energy - Qy;
    if (Qy < 0.0) Qy = 0.0;
    if (Ly < 0.0) Ly = 0.0;
    double Ne = Qy * energy * ScaleFactor[1];
    double Nph = Ly * energy * ScaleFactor[0] *
                 (1. - 1. / pow(1. + pow((energy / NuisParam[7]), NuisParam[8]), NuisParam[11]));
    Nq = Nph + Ne;
    double Ni = (4. / ThomasImel) * (exp(Ne * ThomasImel / 4.) - 1.);
    double Nex = (-1. / ThomasImel) * (4. * exp(Ne * ThomasImel / 4.) - (Ne + Nph) * ThomasImel - 4.);

    double NexONi = Nex / Ni;
    Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
    if (NexONi < wvalue.alpha && energy > 1e2) {
        NexONi = wvalue.alpha;
        Ni = Nq / (1. + NexONi);
        Nex = Nq - Ni;
    }
    if (NexONi > 1.0 && energy < 1.) {
        NexONi = 1.00;
        Ni = Nq / (1. + NexONi);
        Nex = Nq - Ni;
    }

    if (Nex <= 0.)
        cerr
                << "\nCAUTION: You are approaching the border of NEST's validity for high-energy (OR, for LOW) NR, or are beyond it, at "
                << energy << " keV." << endl;
    if (std::abs(Nex + Ni - Nq) > 2. * PHE_MIN) {
        throw std::runtime_error("ERROR: Quanta not conserved. Tell Matthew Immediately!");
    }

    double Wq_eV = wvalue.Wq_eV;
    double L = (Nq / energy) * Wq_eV * 1e-3;

    YieldResult result{};
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = NexONi;
    result.Lindhard = L;
    result.ElectricField = dfield;
    result.DeltaT_Scint = -999;
    return YieldResultValidity(result, energy, Wq_eV);  // everything needed to calculate fluctuations
}

YieldResult NESTcalc::GetYieldIon(double energy, double density, double dfield, double massNum, double atomNum,
                                  const std::vector<double> &NuisParam/*{11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/) { //simply uses the original Lindhard model, but as cited by Hitachi in:
    //https://indico.cern.ch/event/573069/sessions/230063/attachments/1439101/2214448/Hitachi_XeSAT2017_DM.pdf

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
    double L = Constant * sqrt(energy * 1e3);
    double L_max = 0.96446 / (1. + pow(massNum * massNum / 19227., 0.99199));
    if (ValidityTests::nearlyEqual(atomNum, 2.) && ValidityTests::nearlyEqual(massNum, 4.))
        L = 0.56136 * pow(energy, 0.056972);
    if (L > L_max) L = L_max;
    double densDep = pow(density / 0.2679, -2.3245);
    double massDep =
            0.02966094 * exp(0.17687876 * (massNum / 4. - 1.)) + 1. - 0.02966094;
    double fieldDep = pow(1. + pow(dfield / 95., 8.7), 0.0592);
    if (fdetector->get_inGas()) fieldDep = sqrt(dfield);
    double ThomasImel = 0.00625 * massDep / (1. + densDep) / fieldDep;
    if (ValidityTests::nearlyEqual(A1, 206.) &&
        ValidityTests::nearlyEqual(Z1, 82.)) { 
	//Pb-206 (from Po-210 alpha decay).
        ThomasImel = 79.9 * pow(dfield, -0.868); //Nishat Parveen
    }
    const double logden = log10(density);
    double Wq_eV = 28.259 + 25.667 * logden - 33.611 * pow(logden, 2.) -
                   123.73 * pow(logden, 3.) - 136.47 * pow(logden, 4.) -
                   74.194 * pow(logden, 5.) - 20.276 * pow(logden, 6.) -
                   2.2352 * pow(logden, 7.);
    double alpha = 0.64 / pow(1. + pow(density / 10., 2.), 449.61);
    double NexONi = alpha + 0.00178 * pow(atomNum, 1.587); //Wq_eV=13.7;ThomasImel=0.05;NexONi=0.05;
    double Nq = 1e3 * L * energy / Wq_eV;
    double Ni = Nq / (1. + NexONi);
    double recombProb;
    if (Ni > 0. && ThomasImel > 0.)
        recombProb =
                1. - log(1. + (ThomasImel / 4.) * Ni) / ((ThomasImel / 4.) * Ni);
    else {
        recombProb = 0.0;
    }
    double Nph = Nq * NexONi / (1. + NexONi) + recombProb * Ni;
    double Ne = Nq - Nph;
    if (ValidityTests::nearlyEqual(A1, 206.) && ValidityTests::nearlyEqual(Z1, 82.))
        Ne = RandomGen::rndm()->rand_gauss(Ne / ChargeLoss, 2. * sqrt(Ne / (ChargeLoss * ChargeLoss)));
    // to compensate for accidentally including Q-loss in fits to Xed data
    if (ValidityTests::nearlyEqual(Z2, 18.) && ValidityTests::nearlyEqual(Z1, 2.) &&
        ValidityTests::nearlyEqual(A1, 4.)) { //alphas in argon
        double factorE = pow((4.71598 + pow(dfield, 7.72848)), 0.109802);
        double Qy = (1. / 6200.) *
	  (64478398.7663 - (64478398.7663 * 0.173553719 + (64478398.7663 / 1.21) * (1. - (0.02852 * log(1. + (64478398.7663 / 1.21) * ((0.01 / factorE) / 3.))) /
										    ((64478398.7663 / 1.21) * (0.01 / factorE)))));
        double qu = 1. / (1.5 * pow(dfield, -0.012));
        factorE = pow(4.98483 + pow(dfield / 10.0822, 1.2076), 0.97977);
        double Ly= qu * (1. / 6500.) * (278037.250283 * 0.173553719 + (278037.250283 / 1.21) * (1. - (2. * log(1. + (278037.250283 / 1.21) * ((0.653503 / factorE) / 3.))) /
												((278037.250283 / 1.21) * (0.653503 / factorE))));
        Ne = Qy * energy;
        Nph = Ly * energy;
        L = 0.0;
        NexONi = 0.21;
        Wq_eV = 19.5;
    }

    YieldResult result{};
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = NexONi;
    result.Lindhard = L;
    result.ElectricField = dfield;
    result.DeltaT_Scint = -999;
    return YieldResultValidity(result, energy, Wq_eV);  // everything needed to calculate fluctuations
}

YieldResult NESTcalc::GetYieldKr83m(double energy, double density, double dfield, double maxTimeSeparation,
                                    double deltaT_ns = -999) {
    double Nq = -999;
    double Nph = -999;
    double Ne = -999;
    Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
    double Wq_eV = wvalue.Wq_eV;
    double alpha = wvalue.alpha;
    constexpr double deltaT_ns_halflife = 154.4;
    if (ValidityTests::nearlyEqual(energy, 9.4)) {
        if (deltaT_ns < 0) {
            while (deltaT_ns > maxTimeSeparation || deltaT_ns < 0) {
                deltaT_ns = RandomGen::rndm()->rand_exponential(deltaT_ns_halflife);
            }
        }
        if (deltaT_ns < 100. && energy < 41.5 && kr83m_reported_low_deltaT == false) {
            kr83m_reported_low_deltaT = true;
            cerr << "\tWARNING! Past Kr83m model fit validity region. Details: "
                 << " deltaT_ns is <100 ns and your input energy is either 32.1 or 9.4 keV. "
                 << " Data for separated Kr83m decays does not yet exist for deltaT_ns <100 ns. "
                 << " 9.4 & 32.1 keV yields are still summed to physically accurate result, but individually will be nonsensical."
                 << endl;
        }
        Nq = energy * 1e3 / Wq_eV;
        double medTlevel = 57.462 + (69.201 - 57.462) / pow(1. + pow(dfield / 250.13, 0.9), 1.);
        double lowTdrop = 35. + (75. - 35.) / pow(1. + pow(dfield / 60, 1), 1);
        double lowTpeak = 6.2831e4 - (6.2831e4 - 5.949e4) / pow(1. + pow(dfield / 60, 1.), 1);
        Nph = energy * (lowTpeak * pow(2. * deltaT_ns + 10., -1.5) + medTlevel) /
              (1. + pow(deltaT_ns / lowTdrop, -1. * lowTdrop / 5.));
        Ne = Nq - Nph;
        if (Ne < 0)
            Ne = 0.;
        alpha = 0.;
    } else {
        if (ValidityTests::nearlyEqual(energy, 32.1)) {
            Nq = energy * 1e3 / Wq_eV;
            Nph = energy * (6. + (66.742 - 6.) / pow(1. + pow(dfield / 115.037, 0.6409), 0.3215));
            Ne = Nq - Nph;
            if (Ne < 0.)
                Ne = 0.;
        } else {   //merged 41.5 keV decay
            if (deltaT_ns < 0) {
                while (deltaT_ns > maxTimeSeparation || deltaT_ns < 0) {
                    deltaT_ns = RandomGen::rndm()->rand_exponential(deltaT_ns_halflife);
                }
            }
            double medTlevel = 57.462 + (69.201 - 57.462) / pow(1. + pow(dfield / 250.13, 0.9), 1.);
            double lowTdrop = 35. + (75. - 35.) / pow(1. + pow(dfield / 60, 1), 1);
            double lowTpeak = 6.2831e4 - (6.2831e4 - 5.949e4) / pow(1. + pow(dfield / 60, 1.), 1);
            //9.4 keV model
            Nph = 9.4 * (lowTpeak * pow(2. * deltaT_ns + 10., -1.5) + medTlevel) /
                  (1. + pow(deltaT_ns / lowTdrop, -1. * lowTdrop / 5.));
            Ne = (9.4 * 1e3 / Wq_eV) - Nph;
            if (Ne < 0.) Ne = 0.;  //can't have negative charge from either decay
            //Add in 32.1 keV yields

            double Nph_32 = 32.1 * (6. + (66.742 - 6.) / pow(1. + pow(dfield / 115.037, 0.6409), 0.3215));
            double Ne_32 = (32.1 * 1e3 / Wq_eV) - Nph_32;
            Nph += Nph_32;
            Ne += Ne_32;
            if (Ne < 0.) Ne = 0.;
        }
    }

    YieldResult result{};
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = NexONi(energy, density);
    result.Lindhard = 1;
    result.ElectricField = dfield;
    result.DeltaT_Scint = deltaT_ns;
    return YieldResultValidity(result, energy, Wq_eV);  // everything needed to calculate fluctuations
}

YieldResult NESTcalc::GetYieldBeta(double energy, double density, double dfield) {
    Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
    double Qy, Nq;
    double Wq_eV = wvalue.Wq_eV;
    //double alpha = wvalue.alpha; // duplicate definition below. We don't even need this here (it is Nex/Ni)

    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) { 
	// Liquid Argon
        double alpha =
                32.988 - 552.988 / (17.2346 + pow(dfield / (-4.7 + 0.025115 * exp(1.3954 / 0.265360653)), 0.242671));
        double beta = 0.778482 + 25.9 / pow(1.105 + pow(dfield / 0.4, 4.55), 7.502);
        double gamma = 0.659509 * (1000 / 19.5 + 6.5 * (5 - 0.5 / pow(dfield / 1047.408, 0.01851)));
        double delta = 15.7489;
        double DB = 1052.264 + (14159350000 - 1652.264) / (-5 + pow(dfield / 0.157933, 1.83894));
        double p1 = 1;
        double p2 = 10.304;
        double p3 = 13.0654;
        double p4 = 0.10535;
        double p5 = 0.7;
        double LET = -2.07763;
        Nq = energy * 1e3 / Wq_eV;
        Qy = alpha * beta + (gamma - alpha * beta) / pow(p1 + p2 * pow(energy + 0.5, p3), p4) +
             delta / (p5 + DB * pow(energy, LET));
    } else {
        double QyLvllowE = 1e3 / Wq_eV + 6.5 * (1. - 1. / (1. + pow(dfield / 47.408, 1.9851)));
        double HiFieldQy = 1. + 0.4607 / pow(1. + pow(dfield / 621.74, -2.2717), 53.502);
        double QyLvlmedE = 32.988 - 32.988 / (1. + pow(dfield / (0.026715 * exp(density / 0.33926)), 0.6705));
        QyLvlmedE *= HiFieldQy;
        double DokeBirks = 1652.264 + (1.415935e10 - 1652.264) / (1. + pow(dfield / 0.02673144, 1.564691));
        Nq = energy * 1e3 / Wq_eV;  //( Wq_eV+(12.578-Wq_eV)/(1.+pow(energy/1.6,3.5)) );
        double LET_power = -2.;
        if (fdetector->get_inGas()) LET_power = 2.;
        double QyLvlhighE = 28.;
        if (density > 3.100) QyLvlhighE = 49.; //SXe effect from Yoo.
        Qy = QyLvlmedE + (QyLvllowE - QyLvlmedE) / pow(1. + 1.304 * pow(energy, 2.1393), 0.35535) +
             QyLvlhighE / (1. + DokeBirks * pow(energy, LET_power));
        if (Qy > QyLvllowE && energy > 1. && dfield > 1e4) Qy = QyLvllowE;
    }

    double Ly = Nq / energy - Qy;
    double Ne = Qy * energy;
    double Nph = Ly * energy;

    YieldResult result{};
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = NexONi(energy, density);
    result.Lindhard = 1;
    result.ElectricField = dfield;
    result.DeltaT_Scint = -999;
    return YieldResultValidity(result, energy, Wq_eV);  // everything needed to calculate fluctuations;
}

YieldResult
NESTcalc::GetYieldBetaGR(double energy, double density, double dfield, const std::vector<double> &NuisParam) {

    if (RecombOmegaER(0.0, 0.5, NuisParam) < 0.05)
        cerr << "WARNING! You need to change RecombOmegaER to go along with GetYieldBetaGR" << endl;

    Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
    double Wq_eV = wvalue.Wq_eV;
    double alpha = wvalue.alpha;

    double Nq = energy * 1e3 / Wq_eV;
    double m1 = 30.66 + (6.1978 - 30.66) / pow(1. + pow(dfield / 73.855, 2.0318), 0.41883); //NuisParam[0];
    double m5 = Nq / energy / (1 + alpha * erf(0.05 * energy)) - m1;
    double m10 = (0.0508273937 + (0.1166087199 - 0.0508273937) / (1 + pow(dfield / 1.39260460e+02, -0.65763592)));

    double Qy = m1 + (77.2931084 - m1) /
                     pow((1. + pow(energy / (log10(dfield) * 0.13946236 + 0.52561312), 1.82217496 +
                                                                                       (2.82528809 - 1.82217496) / (1 +
                                                                                                                    pow(dfield /
                                                                                                                        144.65029656,
                                                                                                                        -2.80532006)))),
                         0.3344049589)
                + m5 + (0. - m5) / pow((1. + pow(energy / (7.02921301 + (98.27936794 - 7.02921301) /
                                                                        (1. + pow(dfield / 256.48156448, 1.29119251))),
                                                 4.285781736)), m10);

    double coeff_TI = pow(1. / DENSITY, 0.3);
    double coeff_Ni = pow(1. / DENSITY, 1.4);
    double coeff_OL = pow(1. / DENSITY, -1.7) / log(1. + coeff_TI * coeff_Ni * pow(DENSITY, 1.7));
    Qy *= coeff_OL * log(1. + coeff_TI * coeff_Ni * pow(density, 1.7)) * pow(density, -1.7);
    double Ly = Nq / energy - Qy;
    double Ne = Qy * energy;
    double Nph = Ly * energy;
    double lux_NexONi = alpha * erf(0.05 * energy);

    YieldResult result{};
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = lux_NexONi;
    result.Lindhard = 1;
    result.ElectricField = dfield;
    result.DeltaT_Scint = -999;
    return YieldResultValidity(result, energy, Wq_eV);  // everything needed to calculate fluctuations;

}

YieldResult NESTcalc::GetYields(INTERACTION_TYPE species, double energy, double density, double dfield, double massNum,
                                double atomNum, const std::vector<double> &NuisParam
        /*={11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/) {
    switch (species) {
        case NR:
        case WIMP:
        case B8:
        case atmNu:
        case DD:
        case AmBe:
        case Cf:  // this doesn't mean all NR is Cf, this is like a giant if
            // statement. Same intrinsic yields, but different energy spectra
            // (TestSpectra)

            return GetYieldNR(energy, density, dfield, massNum, NuisParam);

            //return GetYieldNROld ( energy, 1 );
            break;
        case ion:
            return GetYieldIon(energy, density, dfield, massNum, atomNum, NuisParam);
            break;
        case gammaRay:
            return GetYieldGamma(energy, density, dfield);
            break;
        case Kr83m:
            return GetYieldKr83m(energy, density, dfield, massNum);
            //not actually massNumber, but a place holder for maxTime
            break;
        case fullGamma_PE:
            return GetYieldGamma(energy, density, dfield); //PE of the full gamma spectrum
            break;
        default:  // beta, CH3T, 14C, the pp solar neutrino background, and Compton/PP spectra of fullGamma
            return GetYieldBeta(energy, density, dfield);
            //return GetYieldBetaGR(energy,density,dfield,NuisParam);
            break;
    }


}

YieldResult NESTcalc::YieldResultValidity(YieldResult &res, const double energy, const double Wq_eV) {
    assert(res.ElectronYield != -999 && res.PhotonYield != -999 && res.ExcitonRatio != -999);
    if (res.PhotonYield > energy / W_SCINT)
        res.PhotonYield = energy / W_SCINT;  // yields can never exceed 1 / [ W ~ few eV ]
    if (res.ElectronYield > energy / W_SCINT) res.ElectronYield = energy / W_SCINT;
    if (res.PhotonYield < 0.) res.PhotonYield = 0.;
    if (res.ElectronYield < 0.) res.ElectronYield = 0.;
    res.Lindhard = max(0., min(res.Lindhard, 1.));  // Lindhard Factor
    if (energy < 0.001 * Wq_eV / res.Lindhard) {
        res.PhotonYield = 0.;
        res.ElectronYield = 0.;
    }
    return res;
}


NESTcalc::NESTcalc(VDetector *detector) {
    assert(detector);
    fdetector = detector;

    photon_areas.reserve(2);
    photon_areas.resize(2);
    scintillation.resize(9);
    newSpike.resize(2);
    ionization.resize(9);

}

NESTcalc::~NESTcalc() {
    if (pulseFile) pulseFile.close();
}

const vector<double> &NESTcalc::GetS1(const QuantaResult &quanta, double truthPosX, double truthPosY, double truthPosZ,
                                      double smearPosX, double smearPosY, double smearPosZ, double driftVelocity,
                                      double dV_mid, INTERACTION_TYPE type_num,
                                      uint64_t evtNum, double dfield, double energy,
                                      S1CalculationMode mode, bool outputTiming,
                                      vector<int64_t> &wf_time,
                                      vector<double> &wf_amp) {

    int Nph = quanta.photons;
    double subtract[2] = {0., 0.};

    photon_areas.clear();
    photon_areas.resize(2);

    // This will clear and reset the vector values
    std::fill(scintillation.begin(), scintillation.end(), 0);

    wf_time.clear();
    wf_amp.clear();

    // Add some variability in g1 drawn from a polynomial spline fit
    double posDep = fdetector->FitS1(truthPosX, truthPosY, truthPosZ, VDetector::fold);
    double posDepSm;
    if (XYcorr == 1 || XYcorr == 3) {
        posDepSm = fdetector->FitS1(smearPosX, smearPosY, smearPosZ, VDetector::unfold);
    } else {
        posDepSm = fdetector->FitS1(0, 0, smearPosZ, VDetector::unfold);
    }

    double dz_center = fdetector->get_TopDrift() -
                       dV_mid * fdetector->get_dtCntr();  // go from t to z
    posDep /=
            fdetector->FitS1(0., 0., dz_center, VDetector::fold);  // XYZ always in mm now never cm
    posDepSm /= fdetector->FitS1(0., 0., dz_center, VDetector::unfold);
    double g1_XYZ = max(0., min(fdetector->get_g1() * posDep, 1.));
    if (ValidityTests::nearlyEqual(fdetector->get_g1(), 1.)) {
        g1_XYZ = 1.;
        posDep = 1.;
        posDepSm = 1.;
    }
    if (ValidityTests::nearlyEqual(fdetector->get_g1(), 0.)) {
        g1_XYZ = 0.;
        posDep = 0.;
        posDepSm = 0.;
    }

    // generate a number of PMT hits drawn from a binomial distribution.
    // Initialize number of photo-electrons
    int nHits = static_cast<int>(BinomFluct(Nph, g1_XYZ));
    int Nphe = 0;

    double eff = fdetector->get_sPEeff();
    if (eff < 1.) {
        eff += ((1. - eff) / (2. * static_cast<double>(fdetector->get_numPMTs()))) * static_cast<double>(nHits);
    }
    //this functional form is just linear approximation taking us from sPEeff at ~few PMTs firing to 100% when there are enough photons detected for there to be *2* in each PMT
    eff = max(0., min(eff, 1.));

    // Initialize the pulse area and spike count variables
    double pulseArea = 0., spike = 0., prob;

    // If single photo-electron efficiency is under 1 and the threshold is above 0
    // (some phe will be below threshold)
    
    S1CalculationMode current_mode = mode;
    if (current_mode == S1CalculationMode::Hybrid)
      {
	//100 keV_ee at 100 V/cm gives approximately 500 photon hits at g1=10%
	current_mode = energy < 100. ? S1CalculationMode::Full : S1CalculationMode::Parametric;
      }
    
    // Follow https://en.wikipedia.org/wiki/Truncated_normal_distribution
    double TruncGaussAlpha = -1. / fdetector->get_sPEres();
    double LittlePhi_Alpha = inv_sqrt2_PI * exp(-0.5 * TruncGaussAlpha * TruncGaussAlpha);
    double BigPhi_Alpha = 0.5 * (1. + erf(TruncGaussAlpha / sqrt2));
    double NewMean = 1. + (LittlePhi_Alpha / (1. - BigPhi_Alpha)) * fdetector->get_sPEres();
    switch (current_mode) {
        case S1CalculationMode::Full:
        case S1CalculationMode::Waveform: {
            // Step through the pmt hits
            for (int i = 0; i < nHits; ++i) {
                // generate photo electron, integer count and area
                // real photo cathodes can't make negative phe
                double TruncGauss = RandomGen::rndm()->rand_zero_trunc_gauss(1. / NewMean, fdetector->get_sPEres());

                double phe1 = TruncGauss +
                              RandomGen::rndm()->rand_gauss(fdetector->get_noiseB()[0],
                                                            fdetector->get_noiseB()[1]);
                ++Nphe;
                if (phe1 > DBL_MAX) {
                    phe1 = 1.; // for Box-Mueller fuck-ups
                }
                prob = RandomGen::rndm()->rand_uniform();
                // zero the area if random draw determines it wouldn't have been observed.
                if (prob > eff) {
                    phe1 = 0.;
                }  // add an else with Nphe++ if not doing mc truth
                // Generate a double photo electron if random draw allows it

                double phe2 = 0.;
                if (RandomGen::rndm()->rand_uniform() < fdetector->get_P_dphe()) {
                    // generate area and increment the photo-electron counter
                    TruncGauss = RandomGen::rndm()->rand_zero_trunc_gauss(1. / NewMean, fdetector->get_sPEres());

                    phe2 = TruncGauss +
                           RandomGen::rndm()->rand_gauss(fdetector->get_noiseB()[0],
                                                         fdetector->get_noiseB()[1]);
                    ++Nphe;
                    if (phe2 > DBL_MAX) {
                        phe2 = 1.;
                    }
                    // zero the area if phe wouldn't have been observed
                    if (RandomGen::rndm()->rand_uniform() > eff &&
                        prob > eff) {
                        phe2 = 0.;
                    }  // add an else with Nphe++ if not doing mc truth
                    // The dphe occurs simultaneously to the first one from the same source
                    // photon. If the first one is seen, so should be the second one
                }
                // Save the phe area and increment the spike count (very perfect spike
                // count) if area is above threshold
                if (current_mode == S1CalculationMode::Waveform) {
                    if ((phe1 + phe2) > fdetector->get_sPEthr()) {
                        pulseArea += phe1 + phe2;
                        ++spike;
                        photon_areas[0].emplace_back(phe1);
                        photon_areas[1].emplace_back(phe2);
                    }
                } else {  
		    // use approximation to find timing
                    if ((phe1 + phe2) > fdetector->get_sPEthr()
                        && (-20. * log(RandomGen::rndm()->rand_uniform()) < fdetector->get_coinWind()
                            || nHits > fdetector->get_coinLevel())) {
                        ++spike;
                        pulseArea += phe1 + phe2;
                    }
                }
            }
            break;
        }
        case S1CalculationMode::Parametric: {
           //Take into account truncated Gaussian shape of SPE and DPE distributions
           //Use the CDFs of the SPE and DPE dists to get the fraction of events below the detectors sPE threshold
           // This will be applied to the digital spike count variable
           double BigPhi_alpha_SPE = 0.5*(1. + erf( -1./fdetector->get_sPEres()/sqrt2 ) );
           double BigPhi_xi_SPE = 0.5*(1. + erf( (fdetector->get_sPEthr()-1.)/fdetector->get_sPEres()/sqrt2 ) ); //where we'll evaluate the CDF
           double sPE_belowThresh_percentile = ( BigPhi_xi_SPE - BigPhi_alpha_SPE ) / (1. - BigPhi_alpha_SPE );
        
           //get the same parameters for the DPE distribution
           double BigPhi_alpha_DPE = 0.5*(1. + erf( -2./(sqrt2*fdetector->get_sPEres())/sqrt2 ) );
           double BigPhi_xi_DPE = 0.5*(1. + erf( (fdetector->get_sPEthr() - 2.)/(sqrt2*fdetector->get_sPEres())/sqrt2) );
           double dPE_belowThresh_percentile = ( BigPhi_xi_DPE - BigPhi_alpha_DPE ) / (1. - BigPhi_alpha_DPE);
        
           double belowThresh_percentile = sPE_belowThresh_percentile * (1. - fdetector->get_P_dphe() )
                                         + dPE_belowThresh_percentile * fdetector->get_P_dphe();
        
           Nphe = nHits + static_cast<int>(BinomFluct(nHits, fdetector->get_P_dphe()));
           eff = fdetector->get_sPEeff();
           if ( eff < 1. ) {
             eff += ((1.-eff)/(2.*static_cast<double>(fdetector->get_numPMTs())))*static_cast<double>(nHits); //same as Full S1CalculationMode case
	   }
	   // Above efficiency adjustment needs to be 'NPhe' instead of 'nHits' for Flamedisx implementation 
           eff = max ( 0., min ( eff, 1. ) );
           auto Nphe_det = static_cast<double>(BinomFluct( Nphe, 1. - ( 1. - eff ) / ( 1. + fdetector->get_P_dphe())));
           //take into account the truncation of the PE distributions
           pulseArea = RandomGen::rndm()->rand_gauss ( Nphe_det*(1./NewMean), fdetector->get_sPEres() * sqrt(Nphe_det) );
           spike = static_cast<double>(BinomFluct(nHits, eff*(1. - belowThresh_percentile)));

	   break;
        }

        case S1CalculationMode::Hybrid: {
        default:

            break;
        }
    }


    if (current_mode == S1CalculationMode::Waveform) {

        vector<double> PEperBin, AreaTable[2], TimeTable[2];
        int numPts = podLength - 100 * SAMPLE_SIZE;
        AreaTable[0].resize(numPts, 0.);
        AreaTable[1].resize(numPts, 0.);

        int total_photons = (int) std::abs(spike);
        int excitons = int(double(total_photons) * double(quanta.excitons) / double(quanta.photons));
        photonstream photon_emission_times =
                GetPhotonTimes(type_num, total_photons, excitons, dfield, energy);
        photonstream photon_times = AddPhotonTransportTime(
                photon_emission_times, truthPosX, truthPosY, truthPosZ);

        if (outputTiming && !pulseFile.is_open()) {
            pulseFile.open("photon_times.txt");
            // pulseFile << "Event#\tt [ns]\tA [PE]" << endl;
            pulseFile << "Event#\tt [ns]\tPEb/bin\tPEt/bin\tin win" << endl;
        }

        int ii, index;
        double min = 1e100, pTime;
        for (ii = 0; ii < (int) std::abs(spike); ++ii) {
            PEperBin.clear();
            PEperBin = fdetector->SinglePEWaveForm(
                    photon_areas[0][ii] + photon_areas[1][ii], photon_times[ii]);
            int total = static_cast<int>(PEperBin.size()) - 1;
            int whichArray;
            if (RandomGen::rndm()->rand_uniform() <
                fdetector->FitTBA(truthPosX, truthPosY, truthPosZ)[0]) {
                whichArray = 0;
            }
            else {
                whichArray = 1;
            }
            for (int kk = 0; kk < total; ++kk) {
                pTime = PEperBin[0] + kk * SAMPLE_SIZE;
                index = int(floor(pTime / SAMPLE_SIZE + 0.5)) + numPts / 2;
                if (index < 0) {
                    index = 0;
                }
                if (index >= numPts) {
                    index = numPts - 1;
                }
                AreaTable[whichArray][index] += 10. * (1. + PULSEHEIGHT)
                                                * (photon_areas[0][ii] + photon_areas[1][ii])
                                                / (PULSE_WIDTH * sqrt2_PI)
                                                * exp(-pow(pTime - photon_times[ii], 2.)
                                                      / (2. * PULSE_WIDTH * PULSE_WIDTH));
            }
            if (total >= 0) {
                if (PEperBin[0] < min) {
                    min = PEperBin[0];
                }
                TimeTable[0].emplace_back(PEperBin[0]);
            }
        }
        double tRandOffset = (SAMPLE_SIZE / 2.) * (2. * RandomGen::rndm()->rand_uniform() -
                                                   1.); //-16,20 was good for LUX, but made weird skew in fP
        for (ii = 0; ii < numPts; ++ii) {
            if ((AreaTable[0][ii] + AreaTable[1][ii]) <= PULSEHEIGHT) {
                continue;
            }

            wf_time.emplace_back((ii - numPts / 2) * SAMPLE_SIZE);
            wf_amp.emplace_back(AreaTable[0][ii] + AreaTable[1][ii]);

            if (outputTiming) {
                char line[80];
                if (AreaTable[0][ii] > PHE_MAX) {
                    subtract[0] = AreaTable[0][ii] - PHE_MAX;
                } else {
                    subtract[0] = 0.0;
                }
                if (AreaTable[1][ii] > PHE_MAX) {
                    subtract[1] = AreaTable[1][ii] - PHE_MAX;
                } else {
                    subtract[1] = 0.0;
                }
                sprintf(line, "%llu\t%lld\t%.3f\t%.3f", evtNum, wf_time.back() + (int64_t) tRandOffset,
                        AreaTable[0][ii] - subtract[0], AreaTable[1][ii] - subtract[1]);
                pulseFile << line << flush;
            }

            if (((ii - numPts / 2) * SAMPLE_SIZE - (int) min) > fdetector->get_coinWind()
                && nHits <= fdetector->get_coinLevel()) {
                pulseArea -= (AreaTable[0][ii] + AreaTable[1][ii]);
                pulseFile << "\t0" << endl;
            } else {
                pulseFile << "\t1" << endl;
            }
        }
        for (ii = 0; ii < TimeTable[0].size(); ++ii) {
            if ((TimeTable[0][ii] - min) > fdetector->get_coinWind() &&
                nHits <= fdetector->get_coinLevel()) {
                --spike;
            }
        }
    }

    if (fdetector->get_noiseL()[0] >= 0.1)
        cerr
                << " !!WARNING!! S1 linear noise term is greater than or equal to 10% (i.e. 0.1) Did you mistake fraction for percent??"
                << endl;

    if (fdetector->get_noiseQ()[0] != 0) {
        pulseArea = RandomGen::rndm()->rand_gauss(pulseArea,
                                                  sqrt(pow(fdetector->get_noiseQ()[0] * pulseArea * pulseArea, 2) +
                                                       pow(fdetector->get_noiseL()[0] * pulseArea, 2)));
    } else {
        pulseArea = RandomGen::rndm()->rand_gauss(pulseArea, fdetector->get_noiseL()[0] * pulseArea);
    }
    pulseArea -= (subtract[0] + subtract[1]);
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
                else {
                    prob = 0.;
		}
            } else if (fdetector->get_coinLevel() == 2) {
                prob = 1. - pow((double) fdetector->get_numPMTs(), 1. - spike);
	    }
            else {
                double numer = 0., denom = 0.;
                for (int i = spike; i > 0; --i) {
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
        GetSpike(Nph, smearPosX, smearPosY, smearPosZ,
                            driftVelocity, dV_mid, scintillation);
        scintillation[6] = newSpike[0];  // S1 spike count (NOT adjusted for double
        // phe effect), IF sufficiently small nHits
        // (otherwise = Nphd)
        scintillation[7] =
                newSpike[1];  // same as newSpike[0], but WITH XYZ correction
    }

    if (RandomGen::rndm()->rand_uniform() > prob && prob < 1.) {  // coincidence has to happen in different PMTs
        // some of these are set to -1 to flag them as having been below threshold
        if (ValidityTests::nearlyEqual(scintillation[0], 0.)) scintillation[0] = PHE_MIN;
        scintillation[0] = -1. * std::abs(scintillation[0]);
        if (ValidityTests::nearlyEqual(scintillation[1], 0.)) scintillation[1] = PHE_MIN;
        scintillation[1] = -1. * std::abs(scintillation[1]);
        if (ValidityTests::nearlyEqual(scintillation[2], 0.)) scintillation[2] = PHE_MIN;
        scintillation[2] = -1. * std::abs(scintillation[2]);
        if (ValidityTests::nearlyEqual(scintillation[3], 0.)) scintillation[3] = PHE_MIN;
        scintillation[3] = -1. * std::abs(scintillation[3]);
        if (ValidityTests::nearlyEqual(scintillation[4], 0.)) scintillation[4] = PHE_MIN;
        scintillation[4] = -1. * std::abs(scintillation[4]);
        if (ValidityTests::nearlyEqual(scintillation[5], 0.)) scintillation[5] = PHE_MIN;
        scintillation[5] = -1. * std::abs(scintillation[5]);
        if (ValidityTests::nearlyEqual(scintillation[6], 0.)) scintillation[6] = PHE_MIN;
        scintillation[6] = -1. * std::abs(scintillation[6]);
        if (ValidityTests::nearlyEqual(scintillation[7], 0.)) scintillation[7] = PHE_MIN;
        scintillation[7] = -1. * std::abs(scintillation[7]);
    }

    // scintillation[8] =
    //fdetector->get_g1();  // g1 (light collection efficiency in liquid)

    return scintillation;
}

const vector<double> &
NESTcalc::GetS2(int Ne, double truthPosX, double truthPosY, double truthPosZ, double smearPosX, double smearPosY,
                double smearPosZ,
                double dt, double driftVelocity, uint64_t evtNum,
                double dfield, S2CalculationMode mode, bool outputTiming,
                vector<int64_t> &wf_time,
                vector<double> &wf_amp,
                const vector<double> &g2_params) {

    double elYield = g2_params[0];
    double ExtEff = g2_params[1];
    double SE = g2_params[2];
    double g2 = g2_params[3];
    double gasGap = g2_params[4];

    std::fill(ionization.begin(), ionization.end(), 0);
    
    double subtract[2] = {0., 0.};
    int i;
    bool eTrain = false;
    if (mode == S2CalculationMode::WaveformWithEtrain && !fdetector->get_inGas()) {
        eTrain = true;
    }

    if (dfield < FIELD_MIN  //"zero"-drift-field detector has no S2
        || elYield <= 0. || ExtEff <= 0. || SE <= 0. || g2 <= 0. ||
        gasGap <= 0.) {

        return ionization;
    }

    // Add some variability in g1_gas drawn from a polynomial spline fit
    double posDep = fdetector->FitS2(
            truthPosX, truthPosY, VDetector::fold);  // XY is always in mm now, never cm
    double posDepSm;
    if (XYcorr == 2 || XYcorr == 3) {
        posDepSm = fdetector->FitS2(smearPosX, smearPosY, VDetector::unfold);
    } else {
        posDepSm = fdetector->FitS2(0, 0, VDetector::unfold);
    }
    posDep /= fdetector->FitS2(0., 0., VDetector::fold);
    posDepSm /= fdetector->FitS2(0., 0., VDetector::unfold);
    double dz = fdetector->get_TopDrift() - dt * driftVelocity;

    if (ValidityTests::nearlyEqual(fdetector->get_g1_gas(), 1.)) {
        posDep = 1.;
        posDepSm = 1.;
    }
    if (ValidityTests::nearlyEqual(fdetector->get_g1_gas(), 0.)) {
        posDep = 0.;
        posDepSm = 0.;
    }

    int Nee = BinomFluct(Ne, ExtEff * exp(-dt / fdetector->get_eLife_us()));
    //MAKE this 1 for SINGLE e- DEBUG

    uint64_t Nph = 0, nHits = 0, Nphe = 0;
    double pulseArea = 0.;

    if (mode == S2CalculationMode::Waveform || mode == S2CalculationMode::WaveformWithEtrain)
    {
        uint64_t k;
        int stopPoint;
        double tau1, tau2, E_liq, amp2;
        vector<double> electronstream, AreaTableBot[2], AreaTableTop[2], TimeTable;
        if (eTrain) {
            stopPoint = BinomFluct(Ne, exp(-dt / fdetector->get_eLife_us()));
        } else {
            stopPoint = Nee;
        }
        electronstream.resize(stopPoint, dt);
        double elecTravT = 0., DL, DL_time, DT, phi, sigX, sigY, newX, newY;
        double Diff_Tran = GetDiffTran_Liquid(dfield, false, fdetector->get_T_Kelvin(), ATOM_NUM);
        double Diff_Long = GetDiffLong_Liquid(dfield, false, fdetector->get_T_Kelvin(), ATOM_NUM);

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
        bool YesGas = true;
        double rho = GetDensity(fdetector->get_T_Kelvin(), fdetector->get_p_bar(), YesGas, 1,
                                fdetector->get_molarMass());
        double driftVelocity_gas =
                GetDriftVelocity_MagBoltz(rho, fdetector->get_E_gas() * 1000.);
        double dt_gas = gasGap / driftVelocity_gas;
        double sigmaDLg = 10. * sqrt(2. * Diff_Long_Gas * dt_gas * 1e-6);
        double sigmaDTg = 10. * sqrt(2. * Diff_Tran_Gas * dt_gas * 1e-6);
        double tauTrap = 0.185;  // microseconds from arXiv:1310.1117, modified to
        // better fit XENON10 and LUX data at same time
        double min = 1e100;
        for (i = 0; i < stopPoint; ++i) {
            elecTravT = 0.;  // resetting for the current electron
            DL = RandomGen::rndm()->rand_gauss(0., sigmaDL);
            DT = std::abs(RandomGen::rndm()->rand_gauss(0., sigmaDT));
            phi = two_PI * RandomGen::rndm()->rand_uniform();
            sigX = DT * cos(phi);
            sigY = DT * sin(phi);
            newX = truthPosX + sigX;
            newY = truthPosY + sigY;
            if (newX * newX + newY * newY >
                fdetector->get_radmax() * fdetector->get_radmax()) {
                continue;  // remove electrons from outside the maximum possible radius
                // (alternative is squeeze them, depending on your E-fields)
            }
            DL_time = DL / driftVelocity;
            elecTravT += DL_time;
            if (!fdetector->get_inGas() && fdetector->get_E_gas() != 0.) {
                elecTravT -= tauTrap * log(RandomGen::rndm()->rand_uniform());
            }
            electronstream[i] += elecTravT;
            // char line[80]; sprintf ( line, "%lu\t%.0f\t%.3f\t%.3f", evtNum,
            // electronstream[i]*1e+3, newX, newY ); pulseFile << line << endl;
            SE = floor(RandomGen::rndm()->rand_gauss(
                    elYield, sqrt(fdetector->get_s2Fano() * elYield)) +
                       0.5);
            Nph += uint64_t(SE);
            SE = (double) BinomFluct(uint64_t(SE), fdetector->get_g1_gas() * posDep);
            nHits += uint64_t(SE);
            double KE = 0.5 * 9.109e-31 * driftVelocity_gas * driftVelocity_gas *
                        1e6 / 1.602e-16;
            double origin = fdetector->get_TopDrift() + gasGap / 2.;
            QuantaResult quanta{};
            quanta.photons = int(SE);
            quanta.electrons = 0;
            quanta.ions = 0;
            quanta.excitons = int(floor(0.0566 * SE + 0.5));
            photonstream photon_emission_times = GetPhotonTimes(NEST::beta, quanta.photons,
                                                                quanta.excitons, dfield, KE);
            photonstream photon_times =
                    AddPhotonTransportTime(photon_emission_times, newX, newY, origin);
            SE += (double) BinomFluct(uint64_t(SE), fdetector->get_P_dphe());
            Nphe += uint64_t(SE);
            DL = RandomGen::rndm()->rand_gauss(0., sigmaDLg);
            DT = RandomGen::rndm()->rand_gauss(0., sigmaDTg);
            phi = two_PI * RandomGen::rndm()->rand_uniform();
            sigX = DT * cos(phi);
            sigY = DT * sin(phi);
            newX += sigX;
            newY += sigY;
            DL_time = DL / driftVelocity_gas;
            electronstream[i] += DL_time;
            if (i >= Nee && eTrain) {  // exponential based on arXiv:1711.07025, power on 2004.07791
                E_liq = fdetector->get_E_gas() / (EPS_LIQ / std::abs(EPS_GAS));
                tau2 = (fdetector->get_TopDrift() / driftVelocity);//0.58313 * exp(0.20929 * E_liq) * 1e3;
                tau1 = 1.40540 * exp(0.15578 * E_liq) * 1e3 * 1e-2;
                amp2 = 0.38157 * exp(0.21177 * E_liq) * 1e-2;
                if (RandomGen::rndm()->rand_uniform() < (amp2 / 0.035) * ExtEff) {
                    tau2 /= RandomGen::rndm()->rand_uniform();
                    if (tau2 > 1e6)
                        tau2 = 1e6; // 1 sec. cap
                    electronstream[i] += tau2;
                } else {
                    electronstream[i] -= tau1 * log(RandomGen::rndm()->rand_uniform());
		}
            }
            for (double j = 0.; j < SE; j += 1.) {
                double phe = RandomGen::rndm()->rand_gauss(1., fdetector->get_sPEres());
                if (i < Nee || !eTrain) pulseArea += phe;
                origin = fdetector->get_TopDrift() +
                         gasGap * RandomGen::rndm()->rand_uniform();
                k = uint64_t(j);
                if (k >= photon_times.size()) k -= photon_times.size();
                double offset = ((fdetector->get_anode() - origin) / driftVelocity_gas +
                                 electronstream[i]) *
                                1e3 +
                                photon_times[k];
                if (offset < min) min = offset;
                if (RandomGen::rndm()->rand_uniform() <
                    fdetector->FitTBA(truthPosX, truthPosY, truthPosZ)[1]) {
                    AreaTableBot[0].emplace_back(phe);
                    AreaTableTop[0].emplace_back(0.0);
                } else {
                    AreaTableBot[0].emplace_back(0.0);
                    AreaTableTop[0].emplace_back(phe);
                }
                TimeTable.emplace_back(offset);
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
            if (ValidityTests::nearlyEqual(AreaTableBot[0][k], 0.0))
                PEperBin =
                        fdetector->SinglePEWaveForm(AreaTableTop[0][k], TimeTable[k] - min);
            else {
                PEperBin =
                        fdetector->SinglePEWaveForm(AreaTableBot[0][k], TimeTable[k] - min);
	    }
            for (i = 0; i < int(PEperBin.size()) - 1; ++i) {
                double eTime = PEperBin[0] + i * SAMPLE_SIZE;
                int index = int(floor(eTime / SAMPLE_SIZE + 0.5));
                if (index < 0) index = 0;
                if (index >= numPts) index = numPts - 1;
                if (ValidityTests::nearlyEqual(AreaTableBot[0][k], 0.0))
                    AreaTableTop[1][index] += 10. * AreaTableTop[0][k] /
                                              (PULSE_WIDTH * sqrt2_PI) *
                                              exp(-pow(eTime - TimeTable[k] + min, 2.) /
                                                  (2. * PULSE_WIDTH * PULSE_WIDTH));
                else {
                    AreaTableBot[1][index] += 10. * AreaTableBot[0][k] /
                                              (PULSE_WIDTH * sqrt2_PI) *
                                              exp(-pow(eTime - TimeTable[k] + min, 2.) /
                                                  (2. * PULSE_WIDTH * PULSE_WIDTH));
		}
            }
        }
        for (k = 0; k < numPts; ++k) {
            if ((AreaTableBot[1][k] + AreaTableTop[1][k]) <= PULSEHEIGHT) continue;
            wf_time.emplace_back(k * SAMPLE_SIZE + int64_t(min + SAMPLE_SIZE / 2.));
            wf_amp.emplace_back(AreaTableBot[1][k] + AreaTableTop[1][k]);

            if (outputTiming) {
                char line[80];
                if (AreaTableBot[1][k] > PHE_MAX) subtract[0] = AreaTableBot[1][k] - PHE_MAX;
                else {
		  subtract[0] = 0.0;
		}
                if (AreaTableTop[1][k] > PHE_MAX) subtract[1] = AreaTableTop[1][k] - PHE_MAX;
                else {
		  subtract[1] = 0.0;
		}
                sprintf(line, "%llu\t%lld\t%.3f\t%.3f", evtNum, wf_time.back(),
                        AreaTableBot[1][k] - subtract[0], AreaTableTop[1][k] - subtract[1]);
                pulseFile << line << endl;
            }
        }
    } else {
        Nph = uint64_t(floor(RandomGen::rndm()->rand_gauss(
                elYield * double(Nee), sqrt(fdetector->get_s2Fano() *
                                            elYield * double(Nee))) +
                             0.5));
        nHits = BinomFluct(Nph, fdetector->get_g1_gas() * posDep);
        Nphe = nHits + BinomFluct(nHits, fdetector->get_P_dphe());
        pulseArea = RandomGen::rndm()->rand_gauss(
                Nphe, fdetector->get_sPEres() * sqrt(Nphe));
    }

    if (fdetector->get_noiseL()[1] >= 0.1)
        cerr
                << " !!WARNING!! S2 linear noise term is greater than or equal to 10% (i.e. 0.1) Did you mistake fraction for percent??"
                << endl;
    if (fdetector->get_noiseQ()[1] != 0) {
        pulseArea = RandomGen::rndm()->rand_gauss(
                pulseArea, sqrt(pow(fdetector->get_noiseQ()[1] * pow(pulseArea, 2), 2) +
                                pow(fdetector->get_noiseL()[1] * pulseArea, 2)));
    } else {
        pulseArea = RandomGen::rndm()->rand_gauss(
                pulseArea, fdetector->get_noiseL()[1] * pulseArea);
    }
    pulseArea -= (subtract[0] + subtract[1]);
    double pulseAreaC =
            pulseArea / exp(-dt / fdetector->get_eLife_us()) / posDepSm;
    double Nphd = pulseArea / (1. + fdetector->get_P_dphe());
    double NphdC = pulseAreaC / (1. + fdetector->get_P_dphe());

    double S2b = RandomGen::rndm()->rand_gauss(
            fdetector->FitTBA(truthPosX, truthPosY, truthPosZ)[1] * pulseArea,
            sqrt(fdetector->FitTBA(truthPosX, truthPosY, truthPosZ)[1] *
                 pulseArea *
                 (1. - fdetector->FitTBA(truthPosX, truthPosY, truthPosZ)[1])));
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

    if (pulseArea < std::abs(fdetector->get_s2_thr())) {
        for (i = 0; i < 8; ++i) {
            if (ValidityTests::nearlyEqual(ionization[i], 0.)) ionization[i] = PHE_MIN;
            ionization[i] = -1. * std::abs(ionization[i]);
        }
    }

    ionization[8] =
            g2;  // g2 = ExtEff * SE, gain of EL in gas gap (from CalculateG2)

    return ionization;
}

vector<double> NESTcalc::CalculateG2(bool verbosity) {
    vector<double> g2_params(5);

    // Set parameters for calculating EL yield and extraction
    double alpha = 0.137, beta = 4.70e-18, gamma = 0;  // note the value of alpha is similar to ~1/7eV. Not coincidence. Noted in Mock et al.
    // actually listed as 'a' and 'b' in ref (below). Units 1/V, cm^2
    double epsilonRatio = EPS_LIQ / std::abs(EPS_GAS);
    if (fdetector->get_inGas())
        epsilonRatio = 1.;  // in an all-gas detector, E_liq variable below simply
    // becomes the field value between anode and gate

    // Convert gas extraction field to liquid field
    double E_liq = fdetector->get_E_gas() / epsilonRatio;  // kV per cm
    double ExtEff;
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) // Argon
        ExtEff = 1. - 1.1974 * exp(-1.003 * pow(E_liq, 1.3849));  // Guschin 1978-79 and 1981-82
    else {
        if (EPS_GAS > 0.) {
            ExtEff = -0.03754 * pow(E_liq, 2.) + 0.52660 * E_liq - 0.84645;  // arXiv:1710.11032 (PIXeY)
            if (E_liq > 7.) ExtEff = 1.;
        } else {
            ExtEff = 1. - 1. / (1. + pow(E_liq / 3.4832, 4.9443));  // arXiv:1904.02885 (Livermore)
            //ExtEff = -0.0012052 + 0.1638 * fdetector->get_E_gas() - 0.0063782 * pow ( fdetector->get_E_gas(), 2. );  // Aprile 2004 IEEE No. 5
            //ExtEff = 1. - 1.3558 * exp ( -0.074312 * pow ( E_liq, 2.4259 ) );//Gus, favored by RED-100
        } //the alternative options
    }
    if (ExtEff > 1. || fdetector->get_inGas()) ExtEff = 1.;
    if (ExtEff < 0. || E_liq <= 0.) ExtEff = 0.;

    double gasGap =
            fdetector->get_anode() -
            fdetector
                    ->get_TopDrift();  // EL gap in mm -> cm, affecting S2 size linearly
    if (gasGap <= 0. && E_liq > 0.) {
        throw std::runtime_error("\tERR: The gas gap in the S2 calculation broke!!!!");
    }

    // Calculate EL yield based on gas gap, extraction field, and pressure
    //double elYield = (alpha * fdetector->get_E_gas() * 1e3 -
    //                beta * fdetector->get_p_bar() - gamma) *
    //               gasGap * 0.1;  // arXiv:1207.2292 (HA, Vitaly C.)
    bool YesGas = true;
    double rho = GetDensity(fdetector->get_T_Kelvin(), fdetector->get_p_bar(), YesGas, 1, fdetector->get_molarMass());
    double elYield;
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) { // Henrique Araujo and Vitaly Chepel again
        alpha = 0.0813;
        beta = 1.90e-18;
    }
    elYield = (alpha * fdetector->get_E_gas() * 1e3 - beta * (NEST_AVO * rho / fdetector->get_molarMass())) * gasGap *
              0.1;
    // replaced with more accurate version also from 1207.2292, but works for room temperature gas
    if (elYield <= 0.0 && E_liq != 0.) {
        cerr << "\tWARNING, the field in gas must be at least "
             << 1e-3 * beta * NEST_AVO * rho / (alpha * fdetector->get_molarMass())
             << " kV/cm, for S2 to work," << endl;
        cerr << "\tOR: your density for gas must be less than "
             << fdetector->get_molarMass() * alpha * fdetector->get_E_gas() * 1e3 / (beta * NEST_AVO) << " g/cm^3."
             << endl;
    }
    // Calculate single electron size and then g2
    double SE = elYield * fdetector->get_g1_gas();  // multiplying by light
    // collection efficiency in
    // the gas gap
    if (fdetector->get_s2_thr() < 0)
        SE *= fdetector->FitTBA(0., 0., fdetector->get_TopDrift() / 2.)[1];
    double g2 = ExtEff * SE;
    double StdDev = 0., Nphe, pulseArea, pulseAreaC, NphdC, phi, posDep, r, x, y;
    int Nph, nHits;
    if (verbosity) {
        for (int i = 0; i < 10000; ++i) { 
            // calculate properly the width (1-sigma std dev) in the SE size
            Nph = int(floor(RandomGen::rndm()->rand_gauss(elYield, sqrt(fdetector->get_s2Fano() * elYield)) + 0.5));
            phi = two_PI * RandomGen::rndm()->rand_uniform();
            r = fdetector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
            x = r * cos(phi);
            y = r * sin(phi);
            posDep = fdetector->FitS2(x, y, VDetector::fold) /
                     fdetector->FitS2(0., 0., VDetector::fold); //future upgrade: smeared pos
            nHits = BinomFluct(Nph, fdetector->get_g1_gas() * posDep);
            Nphe = nHits + BinomFluct(nHits, fdetector->get_P_dphe());
            pulseArea = RandomGen::rndm()->rand_gauss(Nphe, fdetector->get_sPEres() * sqrt(Nphe));
            if (fdetector->get_noiseQ()[1] != 0) {
                pulseArea = RandomGen::rndm()->rand_gauss(
                        pulseArea, sqrt(pow(fdetector->get_noiseQ()[1] * pow(pulseArea, 2), 2) +
                                        pow(fdetector->get_noiseL()[1] * pulseArea, 2)));
            } else {
                pulseArea = RandomGen::rndm()->rand_gauss(
                        pulseArea, fdetector->get_noiseL()[1] * pulseArea);
            }
            if (fdetector->get_s2_thr() < 0.)
                pulseArea = RandomGen::rndm()->rand_gauss(
                        fdetector->FitTBA(0.0, 0.0, fdetector->get_TopDrift() / 2.)[1] * pulseArea, sqrt
                                (fdetector->FitTBA(0.0, 0.0, fdetector->get_TopDrift() / 2.)[1] *
                                 pulseArea * (1. - fdetector->FitTBA(0.0, 0.0, fdetector->get_TopDrift() / 2.)[1])));
            pulseAreaC = pulseArea / posDep;
            NphdC = pulseAreaC / (1. + fdetector->get_P_dphe());
            StdDev += (SE - NphdC) * (SE - NphdC);
        }
        StdDev = sqrt(StdDev) / sqrt(9999.); // N-1 from above (10,000)

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

const vector<double>& NESTcalc::GetSpike(int Nph, double dx, double dy, double dz,
                                  double driftSpeed, double dS_mid,
                                  const vector<double> &oldScint)
{
    std::fill(newSpike.begin(), newSpike.end(), 0);

    if (oldScint[7] > SPIKES_MAXM) {
        newSpike[0] = oldScint[4];
        newSpike[1] = oldScint[5];
        return newSpike;
    }
    newSpike[0] = std::abs(oldScint[6]);
    double TruncGauss = 0.;
    while (TruncGauss <= 0.)
        TruncGauss = RandomGen::rndm()->rand_gauss(newSpike[0], (fdetector->get_sPEres() / 4.) * sqrt(newSpike[0]));
    newSpike[0] = TruncGauss;
    if (XYcorr == 0 || XYcorr == 2) {
        dx = 0.;
        dy = 0.;
    }
    newSpike[1] = newSpike[0] / fdetector->FitS1(dx, dy, dz, VDetector::unfold) *
                  fdetector->FitS1(0., 0., fdetector->get_TopDrift() -
                                           dS_mid * fdetector->get_dtCntr(), VDetector::unfold);

    return newSpike;  // regular and position-corrected spike counts returned
}

double NESTcalc::SetDensity(double Kelvin,
                            double bara) {
    bool inGas = false;
    double density = GetDensity(Kelvin, bara, inGas, 0);
    fdetector->set_inGas(inGas);

    return density;
}

double NESTcalc::GetDensity(double Kelvin, double bara, bool &inGas, uint64_t evtNum, double molarMass) {
    // currently only for fixed pressure (the saturated vapor pressure); will add pressure dependence later

    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.) && !inGas) {
        inGas = false;
        if (DENSITY > 2.) return 1.4;
        else {
            return DENSITY;
	}
    }

    //if (MOLAR_MASS > 134.5) //enrichment for 0vBB expt (~0.8 Xe-136)
    //return 3.0305; // ±0.0077 g/cm^3, EXO-200 @167K: arXiv:1908.04128

    if (Kelvin < 161.40 && ValidityTests::nearlyEqual(ATOM_NUM, 54.)) {  // solid Xenon
        cerr << "\nWARNING: SOLID PHASE. IS THAT WHAT YOU WANTED?\n";
        return 3.41;  // from Yoo at 157K
        // other sources say 3.1 (Wikipedia, 'minimum') and 3.640g/mL at an unknown temperature
    }

    double VaporP_bar;  // we will calculate using NIST
    if (Kelvin < 289.7) {
        if (ValidityTests::nearlyEqual(ATOM_NUM, 54.)) VaporP_bar = pow(10., 4.0519 - 667.16 / Kelvin);
        else {
	  VaporP_bar = 1.;
	}
    } else {
        VaporP_bar = DBL_MAX;
    }
    if (bara < VaporP_bar || inGas) {
        double p_Pa = bara * 1e5;
        double density = 1. /
                         (pow(RidealGas * Kelvin, 3.) / (p_Pa * pow(RidealGas * Kelvin, 2.) + RealGasA * p_Pa * p_Pa) +
                          RealGasB);  // Van der Waals equation, mol/m^3
        density *= molarMass * 1e-6;
        if (bara < VaporP_bar && evtNum == 0) cerr << "\nWARNING: GAS PHASE. IS THAT WHAT YOU WANTED?\n";
        inGas = true;
        return density;  // in g/cm^3
    }
    inGas = false;
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

double NESTcalc::SetDriftVelocity(double Kelvin, double Density, double eField) {
    return GetDriftVelocity(Kelvin, Density, eField, fdetector->get_inGas());
}

double NESTcalc::GetDriftVelocity(double Kelvin, double Density, double eField, bool inGas) {
    if (inGas) return GetDriftVelocity_MagBoltz(Density, eField);
    else {
        return GetDriftVelocity_Liquid(Kelvin, eField, Density);
    }
}

double NESTcalc::GetDriftVelocity_Liquid(double Kelvin, double eField,
                                         double Density) {  // for liquid and solid only
    double speed =
            0.0;  // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
    int i, j;
    double vi, vf, slope, Ti, Tf, offset;

    if (ValidityTests::nearlyEqual(ATOM_NUM,
                                   18.)) { 
      //replace eventually with Kate's work, once T's splined as done for LXe and SXe
      /*x_nest is field in kV/cm:
      Liquid:
      y_nest_85 = 2.553969*(0.929235**(1/x_nest))*x_nest**0.315673 // temperatures in Kelvin
      y_nest_87 = 2.189388*(0.961701**(1/x_nest))*x_nest**0.339414
      y_nest_89 = 2.201489*(0.931080**(1/x_nest))*x_nest**0.320762
      y_nest_94 = 1.984534*(0.936081**(1/x_nest))*x_nest**0.332491
      y_nest_100 = 2.150101*(0.929447**(1/x_nest))*x_nest**0.317973
      y_nest_120 = 1.652059*(0.935714**(1/x_nest))*x_nest**0.329025
      y_nest_130 = 1.273891*(0.965041**(1/x_nest))*x_nest**0.336455
      Solid:
      y_nest_80 = 7.063288*(0.929753**(1/x_nest))*x_nest**0.314640
      y_nest_82 = 5.093097*(0.920459**(1/x_nest))*x_nest**0.458724*/
        speed = 0.097384 * pow(log10(eField), 3.0622) - 0.018614 * sqrt(eField);
        if (speed < 0.) speed = 0.;
        return speed;
    }

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
            {1.2466,  85.975, -.88005, 918.57, -3.0085, 27.568, 2.3823},   // 200
            {334.60,  37.556, 0.92211, 345.27, -338.00, 37.346, 1.9834}};  // 230

    double Temperatures[11] = {100., 120., 140., 155., 157., 163.,
                               165., 167., 184., 200., 230.};

    if (Kelvin < 100. || Kelvin > 230.) {
        cerr << "\nWARNING: TEMPERATURE OUT OF RANGE (100-230 K) for vD\n";
        if (Kelvin < 100.) Kelvin = 100.;
        if (Kelvin > 230.) Kelvin = 230.;
        cerr << "Using value at closest temp for a drift speed estimate\n";
    }

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
        throw std::runtime_error("ERROR: TEMPERATURE OUT OF RANGE (100-230 K)");
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
    if (ValidityTests::nearlyEqual(Kelvin, Ti)) return vi;
    if (ValidityTests::nearlyEqual(Kelvin, Tf)) return vf;
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
        cerr << "\nWARNING: DRIFT SPEED NON-POSITIVE. Setting to 0.1 mm/us\t" <<
             "Line Number ~1700 of NEST.cpp, in function NESTcalc::GetDriftVelocity_Liquid\t" <<
             "Stop bothering Matthew about this, and fix underlying cause in your code!\n";
        if (eField < 1e2 && eField >= FIELD_MIN) {
            cerr << "FIELD MAY BE TOO LOW. ";
        } else if (eField > 1e4) {
            cerr << "FIELD MAYBE TOO HIGH. ";
        } else {
            cerr << "Something unknown went wrong: are you in a noble element?? ";
        }
        cerr << "EF = " << eField << " V/cm. T = " << Kelvin << " Kelvin" << endl;
        speed = 0.1;
    }
    return speed; //mm per microsecond
}

double NESTcalc::GetDriftVelocity_MagBoltz(
        double density, double efieldinput, double molarMass)  // Nichole Barry UCD 2011
{
    density *= NEST_AVO / molarMass;
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

    return std::abs(edrift) * 1e-5;  // from cm/s into mm per microsecond
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
                                                          (EPS_LIQ / std::abs(EPS_GAS)) * 1e3);
                else { 
		    // if gate == TopDrift properly set, shouldn't happen
                    driftTime += zStep / GetDriftVelocity_MagBoltz(
                            rho, fdetector->get_E_gas() * 1e3);
		}
            } else {
                driftTime +=
                        zStep /
                        SetDriftVelocity(
                                fdetector->get_T_Kelvin(), rho,
                                fdetector->FitEF(dx, dy,
                                                 zz));  // update x and y if you want 3-D fields
	    }
        }

        speedTable.emplace_back((zz - pos_z) / driftTime);  // uses highest zz
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

    double phi = two_PI * RandomGen::rndm()->rand_uniform();
    sigmaR = std::abs(RandomGen::rndm()->rand_gauss(0.0, sigmaR));
    double sigmaX = sigmaR * cos(phi);
    double sigmaY = sigmaR * sin(phi);

    if (sigmaR > 1e2 || std::isnan(sigmaR) || sigmaR <= 0. ||
        std::abs(sigmaX) > 1e2 || std::abs(sigmaY) > 1e2) {
        if (A_top > 20.) {
            cerr << "WARNING: your position resolution is worse than 10 cm. Is that correct?!" << endl;
            cerr << "Setting resolution to perfect." << endl;
            sigmaX = 0.;
            sigmaY = 0.;
        } // this is only a problem if the area is large and the res is still bad
    }

    xySmeared[0] = xPos_mm + sigmaX;
    xySmeared[1] = yPos_mm + sigmaY;

    return xySmeared;  // new X and Y position in mm with empirical smearing. LUX
    // Run03 example
}

double NESTcalc::PhotonEnergy(bool s2Flag, bool state, double tempK) {
    double wavelength, E_keV;  // wavelength is in nanometers

    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) { 
	// liquid Argon
        if (state) return RandomGen::rndm()->rand_gauss(9.7, 0.2);
        return RandomGen::rndm()->rand_gauss(9.69, 0.22);
    }

    if (!state)  // liquid or solid
        wavelength = RandomGen::rndm()->rand_gauss(
                178., 14. / 2.355);  // taken from Jortner JchPh 42 '65
    else {                     // gas
        wavelength =
                RandomGen::rndm()->rand_gauss(175., 5.);  // G4S1Light, probably Doke
    }
    if (s2Flag) {  
	// S2 different from ordinary gas (or just measurement error?)
        if (tempK < 200.)  // cold gas
            wavelength = RandomGen::rndm()->rand_gauss(
                    179., 5.);  // source is G4S2Light.cc from the old NEST
        else {
            wavelength = RandomGen::rndm()->rand_gauss(174., 5.);  // ditto
	}
    }

    E_keV = 1240e-3 / wavelength;  // h*c in keV-nm divided by lambda in nm
    if (E_keV > W_SCINT)
        E_keV = W_SCINT;  // don't go so high breaks G4. Gauss in lambda -> non-G in
    // E, high tail

    return E_keV * 1000.;  // convert from keV into eV. Eventually add full T, P
    // dependence
}

double NESTcalc::CalcElectronLET(double E, int Z) {
    double LET;

    // use a spline fit to online ESTAR data
    if (ValidityTests::nearlyEqual(ATOM_NUM, 54.)) {
        if (E >= 1.)
            LET = 58.482 - 61.183 * log10(E) + 19.749 * pow(log10(E), 2) +
                  2.3101 * pow(log10(E), 3) - 3.3469 * pow(log10(E), 4) +
                  0.96788 * pow(log10(E), 5) - 0.12619 * pow(log10(E), 6) +
                  0.0065108 * pow(log10(E), 7);
            // at energies <1 keV, use a different spline, determined manually by
            // generating sub-keV electrons in Geant4 and looking at their ranges, since
            // ESTAR does not go this low (4.9.4)
        else if (E > 0. && E < 1.) {
            LET = 6.9463 + 815.98 * E - 4828 * pow(E, 2) + 17079 * pow(E, 3) -
                  36394 * pow(E, 4) + 44553 * pow(E, 5) - 28659 * pow(E, 6) +
                  7483.8 * pow(E, 7);
	}
        else {
            LET = 0.;
	}
    } else { 
	//replace with Justin and Prof. Mooney's work
        if (E >= 1.)
            LET = 116.70 - 162.97 * log10(E) + 99.361 * pow(log10(E), 2) -
                  33.405 * pow(log10(E), 3) + 6.5069 * pow(log10(E), 4) -
                  0.69334 * pow(log10(E), 5) + .031563 * pow(log10(E), 6);
        else if (E > 0. && E < 1.) {
	  LET = 100.;
	}
        else {
            LET = 0.;
	}
    }

    return LET;
}

NESTcalc::Wvalue NESTcalc::WorkFunction(double density, double MolarMass, bool rmQuanta) {
    double alpha, Wq_eV;
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) { 
	// Liquid argon
        alpha = 0.21;
        Wq_eV = 1000. / 51.9; //23.6/1.21; // ~19.2-5 eV
        return Wvalue{.Wq_eV=Wq_eV, .alpha=alpha};
    }
    alpha = 0.067366 + density * 0.039693;
    /*double xi_se = 9./(1.+pow(density/2.,2.));
    double I_ion = 9.+(12.13-9.)/(1.+pow(density/2.953,65.));
    double I_exc = I_ion / 1.46;
    double Wq_eV = I_exc*(alpha/(1.+alpha))+I_ion/(1.+alpha)
    +xi_se/(1.+alpha);*/
    double eDensity = (density / MolarMass) * NEST_AVO * ATOM_NUM;
    Wq_eV = 18.7263 - 1.01e-23 * eDensity;
    if (rmQuanta) Wq_eV *= InfraredER; //EXO
    return Wvalue{.Wq_eV=Wq_eV, .alpha=alpha}; //W and Nex/Ni together
}

double NESTcalc::NexONi(double energy, double density) {
    if (ValidityTests::nearlyEqual(ATOM_NUM, 18.)) return 0.21; // https://arxiv.org/pdf/1903.05815.pdf
    Wvalue wvalue = WorkFunction(density, fdetector->get_molarMass(), fdetector->get_rmQuanta());
    double alpha = wvalue.alpha;
    return alpha * erf(0.05 * energy);
}

//This function returns the transverse diffusion coefficient in liquid. It allows a user
//to select whether they use the canonical NEST model, or a model modified to accommodate higher
//field values (from Boyle et al., 2016, arXiv:1603.04157v1)
double NESTcalc::GetDiffTran_Liquid(double dfield, bool highFieldModel, double Kelvin,
                                    int Z) // for gas: look for Diff_Tran_Gas above
{
    double output;

    if (Z == 18) { 
	// G4S2Light.cc from LUXSim
        double nDensity = NEST_AVO * DENSITY / 40.; // 1 over cm^3
        return 93.342 * pow(dfield / nDensity, 0.041322);
    }

    //Use the standard NEST parametrization
    if (!highFieldModel) {
        output = 37.368 * pow(dfield, .093452) *
                 exp(-8.1651e-5 * dfield);  // arXiv:1609.04467 (EXO-200)
    }
        //Use the Boyle model, which is drastically different at high (>5kV/cm) fields. Note here that
        //the Boyle model is only at one temperature
        //First double in pair is field, second is DT
    else {
        const std::vector<std::pair<double, double> > BoyleModelDT = GetBoyleModelDT();
        output = interpolateFunction(BoyleModelDT, dfield, true);
        if (output == 0)
            cerr << "Looks like your desired drift field, " << dfield
                 << ", may be either too low or too high to interpolate a DT. Returning DT=0.\n";
    }
    return output;
}

//This function returns the longitudinal diffusion coefficient in liquid. It allows a user
//to select whether they use the canonical NEST model, or a model modified to accommodate higher
//field values (from Boyle et al., 2016, arXiv:1603.04157v1)
double NESTcalc::GetDiffLong_Liquid(double dfield, bool highFieldModel, double Kelvin,
                                    int Z) // for gas: look for Diff_Long_Gas above
{
    double output;

    if (Z == 18) {
        output = GetDiffTran_Liquid(dfield, false, Kelvin, 18);
        return 0.15 * output; // lacking data, just assume that D_L = 0.15 * D_T
    }

    //Use the standard NEST parametrization
    if (!highFieldModel) {
        output = 345.92 * pow(dfield, -0.47880) *
                 exp(-81.3230 / dfield);  // fit to Aprile & Doke review
        // paper and to arXiv:1102.2865;
        // plus, LUX Run02+03
    }
        //Use the Boyle model, which is drastically different at high (>5kV/cm) fields. Note here that
        //the Boyle model is only at one temperature.
        //First double in pair is field, second is DL
    else {
        const std::vector<std::pair<double, double> > BoyleModelDL = GetBoyleModelDL();
        output = interpolateFunction(BoyleModelDL, dfield, true);
        if (output == 0)
            cerr
                    << "Looks like your desired drift field may be either too low or too high to interpolate a DL Returning DL=0.\n";
    }
    return output;
}

//Simple function for interpolating
double NESTcalc::interpolateFunction(const std::vector<std::pair<double, double> > &func, double x, bool isLogLog) {
    //Linear interpolation
    if (!isLogLog) {
        double y_desired = 0;
        for (int iP = 0; iP < func.size() - 1; ++iP) {
            double x1 = func[iP].first;
            double x2 = func[iP + 1].first;
            if (x1 < x && x2 >= x) {
                double y1 = func[iP].second;
                double y2 = func[iP + 1].second;
                double slope = (y2 - y1) / (x2 - x1);
                y_desired = y1 + (slope * (x - x1));
                break;
            }
        }
        return y_desired;
    }

        //Linear interpolation on a log-log scale (more accurate at high fields, where points are far apart on a linear scale)
    else {
        double y_desired = 0;
        for (int iP = 0; iP < func.size() - 1; ++iP) {
            double logx1 = log10(func[iP].first);
            double logx2 = log10(func[iP + 1].first);
            if (logx1 < log10(x) && logx2 > log10(x)) {
                double logy1 = log10(func[iP].second);
                double logy2 = log10(func[iP + 1].second);
                double slope = (logy2 - logy1) / (logx2 - logx1);
                double logy_desired = logy1 + slope * (log10(x) - logx1);
                y_desired = pow(10, logy_desired);
                break;
            }
        }
        return y_desired;
    }
}

// Organized nicely to just input the Boyle Model's curve data. Using vectors and pairs so
// all of the size accounting is done nicely downstream without having to pass container sizes around.
std::vector<std::pair<double, double> > NESTcalc::GetBoyleModelDT() {
    std::vector<std::pair<double, double> > output;
    const int nPts = 25;
    double modelDT[nPts][2] = {{14.6236, 24.674},
                               {24.5646, 26.4954},
                               {36.675,  29.2043},
                               {53.4802, 32.6845},
                               {74.3939, 37.9944},
                               {111.071, 45.3424},
                               {173.835, 58.2226},
                               {306.106, 74.6236},
                               {467.921, 89.2124},
                               {749.809, 96.9913},
                               {1093.38, 98.4549},
                               {1711.25, 97.3852},
                               {2615.85, 92.5493},
                               {3998.65, 85.4419},
                               {6258.24, 78.2405},
                               {9794.71, 67.7294},
                               {16453.1, 58.544},
                               {27637.8, 50.3599},
                               {38445.7, 48.5424},
                               {49828.3, 43.4695},
                               {70967.5, 36.8038},
                               {98719.8, 32.7631},
                               {127948,  29.5377},
                               {169785,  27.4412},
                               {220053,  27.9686}};
    for (auto &iP : modelDT) {
        std::pair<double, double> thePair(iP[0], iP[1]);
        output.emplace_back(iP[0], iP[1]);
    }
    return output;
}

// Organized nicely to just input the Boyle Model's curve data. Using vectors and pairs so
// all of the size accounting is done nicely downstream without having to pass container sizes around. Returns cm^2/s
std::vector<std::pair<double, double> > NESTcalc::GetBoyleModelDL() {
    std::vector<std::pair<double, double> > output;
    const int nPts = 25;
    double modelDL[nPts][2] = {{14.6236, 22.2977},
                               {24.5646, 22.4238},
                               {36.675,  22.933},
                               {53.4802, 23.4595},
                               {74.3939, 25.4994},
                               {111.071, 29.4632},
                               {173.835, 33.9251},
                               {306.106, 38.1019},
                               {467.921, 37.2352},
                               {749.809, 31.8425},
                               {1093.38, 25.495},
                               {1711.25, 18.3343},
                               {2615.85, 12.6461},
                               {3998.65, 8.51594},
                               {6258.24, 5.23814},
                               {9794.71, 3.01815},
                               {16453.1, 1.58294},
                               {27637.8, 0.779995},
                               {38445.7, 0.525689},
                               {49828.3, 0.364016},
                               {70967.5, 0.228694},
                               {98719.8, 0.148206},
                               {127948,  0.116139},
                               {169785,  0.0957632},
                               {220053,  0.0900259}};
    for (auto &iP : modelDL) {
        std::pair<double, double> thePair(iP[0], iP[1]);
        output.emplace_back(iP[0], iP[1]);
    }
    return output;
}
