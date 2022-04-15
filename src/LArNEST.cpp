/**
 * @file LArNEST.cpp
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Justin Mueller [Justin.Mueller@colostate.edu]
 * @author Michael Mooney [mrmooney@colostate.edu]
 * @brief 
 * @version
 * @date 2022-04-13 
 */
#include "LArNEST.hh"

namespace NEST
{
    LArNEST::LArNEST(VDetector *detector)
    : NESTcalc(detector)
    {

    }

    NESTresult LArNEST::FullCalculation(
        INTERACTION_TYPE species, double energy, 
        double density, double dfield,
        double A, double Z,
        const std::vector<double>
        &NuisParam /*={11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/,
        const std::vector<double> &FreeParam /*={1.,1.,0.1,0.5,0.19,2.25}*/,
        bool do_times /*=true*/
    ) 
    {
        if (density < 1.) fdetector->set_inGas(true);
        NESTresult result;
        result.yields = GetYields(species, energy, density, dfield, A, Z, NuisParam);
        result.quanta = GetQuanta(result.yields, density, FreeParam);
        if (do_times)
            result.photon_times = GetPhotonTimes(
                species, result.quanta.photons, result.quanta.excitons, dfield, energy
            );
        else {
            result.photon_times = photonstream(result.quanta.photons, 0.0);
        }
        return result;
    }

    inline double LArNEST::FanoER() 
    {
        // T&I. conflicting reports of .107 (Doke) & ~0.1 elsewhere
        return 0.1115;  
    }

    double LArNEST::GetDensity(
        double Kelvin, double bara, bool &inGas,
        uint64_t evtNum, double molarMass
    ) 
    {
        // currently only for fixed pressure (the saturated vapor pressure); will add
        // pressure dependence later
        double density = 0.;
        double VaporP_bar = 0.;
        VaporP_bar = exp(45.973940 - 1464.718291 / Kelvin - 6.539938 * log(Kelvin));
        if (bara < VaporP_bar || inGas) 
        {
            double p_Pa = bara * 1e5;
            density =
                1. /
                (pow(RidealGas * Kelvin, 3.) /
                    (p_Pa * pow(RidealGas * Kelvin, 2.) + RealGasA * p_Pa * p_Pa) +
                RealGasB);  // Van der Waals equation, mol/m^3
            density *= molarMass * 1e-6;
            if (bara < VaporP_bar && evtNum == 0)
                cerr << "\nWARNING: ARGON GAS PHASE. IS THAT WHAT YOU WANTED?\n";
            inGas = true;
            return density;
        }
        inGas = false;
        return 1.4;
    }

    double LArNEST::GetDriftVelocity_Liquid(
        double Kelvin, double eField,
        double Density
    ) 
    {  
        // for liquid and solid only
        // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
        double speed = 0.0;  
        int i, j;
        double vi, vf, slope, Ti, Tf, offset;

        // replace eventually with Kate's work, once T's splined as done for LXe and
        // SXe
        /*x_nest is field in kV/cm:
        Liquid:
        y_nest_85 = 2.553969*(0.929235**(1/x_nest))*x_nest**0.315673 // temperatures
        in Kelvin y_nest_87 = 2.189388*(0.961701**(1/x_nest))*x_nest**0.339414
        y_nest_89 = 2.201489*(0.931080**(1/x_nest))*x_nest**0.320762
        y_nest_94 = 1.984534*(0.936081**(1/x_nest))*x_nest**0.332491
        y_nest_100 = 2.150101*(0.929447**(1/x_nest))*x_nest**0.317973
        y_nest_120 = 1.652059*(0.935714**(1/x_nest))*x_nest**0.329025
        y_nest_130 = 1.273891*(0.965041**(1/x_nest))*x_nest**0.336455
        Solid:
        y_nest_80 = 7.063288*(0.929753**(1/x_nest))*x_nest**0.314640
        y_nest_82 = 5.093097*(0.920459**(1/x_nest))*x_nest**0.458724*/
        double Temperature[8] = {84., 86., 88., 92., 96., 110., 125., 140.};

        if (Kelvin < 84. || Kelvin > 140.) 
        {
            cerr << "\nWARNING: TEMPERATURE OUT OF RANGE (84-140 K) for vD\n";
            cerr << "Using value at closest temp for a drift speed estimate\n";
            Kelvin = (double)NESTcalc::clamp(int(Kelvin), 84, 140);
        }

        if (Kelvin >= Temperature[0] && Kelvin < Temperature[1]) {
            speed = exp(0.937729 - 0.0734108 / (eField / 1000) +
                        0.315338 * log(eField / 1000));
        }
        else if (Kelvin >= Temperature[1] && Kelvin < Temperature[2]) {
            speed = exp(0.80302379 - 0.06694564 / (eField / 1000) +
                        0.331798 * log(eField / 1000));
        }
        else if (Kelvin >= Temperature[2] && Kelvin < Temperature[3]) {
            speed = exp(0.7795972 - 0.0990952 / (eField / 1000) +
                        0.320876 * log(eField / 1000));
        }
        else if (Kelvin >= Temperature[3] && Kelvin < Temperature[4]) {
            speed = exp(0.6911897 - 0.092997 / (eField / 1000) +
                        0.3295202 * log(eField / 1000));
        }
        else if (Kelvin >= Temperature[4] && Kelvin < Temperature[5]) {
            speed = exp(0.76551511 - 0.0731659 / (eField / 1000) +
                        0.317972 * log(eField / 1000));
        }
        else if (Kelvin >= Temperature[5] && Kelvin < Temperature[6]) {
            speed = exp(0.502022794 - 0.06644517 / (eField / 1000) +
                        0.3290246 * log(eField / 1000));
        }
        else if (Kelvin >= Temperature[6] && Kelvin <= Temperature[7]) {
            speed = exp(0.24207633 - 0.03558428 / (eField / 1000) +
                        0.33645519 * log(eField / 1000));
        }

        if (speed < 0.) {
            speed = 0.;
        }
        return speed;
    }

    double LArNEST::GetDriftVelocity_MagBoltz(
        double density, double efieldinput,
        double molarMass
    )  
    {
        // Nichole Barry UCD 2011
        molarMass = 39.948;
        density *= NEST_AVO / molarMass;
        double edrift = 0., gasdep = efieldinput / density;
        edrift = 2.991205 * pow(1.00113657, 1. / (gasdep * 1.e17)) *
                pow(gasdep * 1.e17, 0.2570253);
        return edrift;
    }

    double LArNEST::PhotonEnergy(bool state) 
    {
        // liquid Argon
        if (state) {
            return RandomGen::rndm()->rand_gauss(9.7, 0.2);
        }
        else {
            return RandomGen::rndm()->rand_gauss(9.69, 0.22);
        }
    }

    inline LArNEST::Wvalue LArNEST::WorkFunction() 
    {
        // Liquid argon
        double alpha = 0.21;
        double Wq_eV = 1000. / 51.9;  // 23.6/1.21; // ~19.2-5 eV
        return Wvalue{.Wq_eV = Wq_eV, .alpha = alpha};
    }

    inline double LArNEST::NexONi() 
    {
        // https://arxiv.org/pdf/1903.05815.pdf
        return 0.21;  
    }
    double LArNEST::PhotonTime(
        INTERACTION_TYPE species, bool exciton,
        double dfield, double energy
    ) 
    {
        // arXiv:1802.06162. NR may need tauR ~0.5-1ns instead of 0
        double time_ns = 0.;
        double SingTripRatio = 0.;
        double tauR = 0.;
        double tau3 = 23.97;
        double tau1 = 3.27;  
        if (fdetector->get_inGas() || energy < W_DEFAULT * 0.001) 
        {
            // from G4S1Light.cc in old NEST
            tau1 = 5.18;   // uncertainty of 1.55 ns from G4S2Light
            tau3 = 100.1;  // uncertainty of 7.90 ns from G4S2Light
        }
        // tau1 = 3.5*ns; tau3 = 20.*ns; tauR = 40.*ns for solid Xe from old NEST.
        // Here assuming same as in liquid

        // if ( dfield > 60. ) dfield = 60. // makes Xed work. 200 for LUX Run04
        // instead. A mystery! Why no field dep?

        double LET = CalcElectronLET(energy, ATOM_NUM);
        // copied from 2013 NEST version for LAr on LBNE
        tau1 = RandomGen::rndm()->rand_gauss(6.5, 0.8);  // error from weighted average
        tau3 = RandomGen::rndm()->rand_gauss(1300, 50);  // ibid.
        tauR = RandomGen::rndm()->rand_gauss(0.8, 0.2);  // Kubota 1979
        if (species <= Cf) {
            SingTripRatio = 0.22218 * pow(energy, 0.48211);
        }
        else if (species == ion) 
        {  // really only alphas here
            SingTripRatio = (-0.065492 + 1.9996 * exp(-energy / 1e3)) /
                            (1. + 0.082154 / pow(energy / 1e3, 2.)) +
                        2.1811;  // uses energy in MeV not keV
        } 
        else 
        {
            SingTripRatio = 0.2701 + 0.003379 * LET - 4.7338e-5 * pow(LET, 2.) +
                            8.1449e-6 * pow(LET, 3.);
            if (LET < 3. && !exciton) {
                SingTripRatio = RandomGen::rndm()->rand_gauss(0.5, 0.2);
            }
            if (LET < 3. && exciton) {
                SingTripRatio = RandomGen::rndm()->rand_gauss(0.36, 0.06);
            }
        }  // lastly is ER for LAr
        if (fdetector->get_inGas() || energy < W_DEFAULT * 0.001) 
        {
            SingTripRatio = 0.1;
            if (fdetector->get_inGas() && !exciton) {
                tauR = 28e3;
            }
            else {
                tauR = 0.;  // 28 microseconds comes from Henrique:
                        // https://doi.org/10.1016/j.astropartphys.2018.04.006
            }
            // from old G4S2Light
            tau3 = 1600.;
            tau1 = 6.; 
        }
        if (tauR < 0.) {
            tauR = 0.;  
        } // in case varied with Gaussian earlier
        // the recombination time is non-exponential, but approximates
        // to exp at long timescales (see Kubota '79)
        time_ns += tauR * (1.0 / RandomGen::rndm()->rand_uniform() - 1.);

        if (RandomGen::rndm()->rand_uniform() < SingTripRatio / (1. + SingTripRatio)) {
            time_ns -= tau1 * log(RandomGen::rndm()->rand_uniform());
        }
        else {
            time_ns -= tau3 * log(RandomGen::rndm()->rand_uniform());
        }
        return time_ns;
    }

    QuantaResult LArNEST::GetQuanta(
        const YieldResult &yields, double density,
        const std::vector<double> &FreeParam /*={1.,1.,0.1,0.5,0.19,2.25}*/
    ) 
    {
        QuantaResult result{};
        bool HighE;
        int Nq_actual, Ne, Nph, Ni, Nex;

        if (FreeParam.size() < 6) 
        {
            throw std::runtime_error(
                "ERROR: You need a minimum of 6 free parameters for the resolution "
                "model.");
        }

        double excitonRatio = yields.ExcitonRatio;
        double Nq_mean = yields.PhotonYield + yields.ElectronYield;

        double elecFrac = max(0., min(yields.ElectronYield / Nq_mean, 1.));

        if (excitonRatio < 0.) 
        {
            excitonRatio = 0.;
            HighE = true;
        } 
        else {
            HighE = false;
        }

        double alf = 1. / (1. + excitonRatio);
        double recombProb = 1. - (excitonRatio + 1.) * elecFrac;
        if (recombProb < 0.) {
            excitonRatio = 1. / elecFrac - 1.;
        }

        if (ValidityTests::nearlyEqual(yields.Lindhard, 1.)) 
        {
            double Fano = FanoER();
            Nq_actual = int(floor(
                RandomGen::rndm()->rand_gauss(Nq_mean, sqrt(Fano * Nq_mean)) + 0.5)
            );
            if (Nq_actual < 0 || ValidityTests::nearlyEqual(Nq_mean, 0.)) {
                Nq_actual = 0;
            }
            Ni = RandomGen::rndm()->binom_draw(Nq_actual, alf);
            Nex = Nq_actual - Ni;

        } 
        else {
            double Fano = FreeParam[0];
            Ni = int(floor(RandomGen::rndm()->rand_gauss(Nq_mean * alf,
                                                        sqrt(Fano * Nq_mean * alf)) +
                        0.5));
            if (Ni < 0) {
                Ni = 0;
            }
            Fano = FreeParam[1];
            Nex = int(floor(RandomGen::rndm()->rand_gauss(
                Nq_mean * excitonRatio * alf,
                sqrt(Fano * Nq_mean * excitonRatio * alf)) + 0.5)
            );
            if (Nex < 0) {
                Nex = 0;
            }
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

        if (Nex <= 0 && HighE) recombProb = yields.PhotonYield / double(Ni);
        recombProb = max(0., min(recombProb, 1.));
        if (std::isnan(recombProb) || std::isnan(elecFrac) || Ni == 0 ||
            ValidityTests::nearlyEqual(recombProb, 0.0)) 
        {
            result.photons = Nex;
            result.electrons = Ni;
            elecFrac = 1.0;
            result.recombProb = 0.;
            result.Variance = 0.;
            return result;
        }

        // set omega (non-binomial recombination fluctuations parameter) according to
        // whether the Lindhard <1, i.e. this is NR.
        double omega = 0.0;  // Ar has no non-binom sauce
        double Variance = recombProb * (1. - recombProb) * Ni + omega * omega * Ni * Ni;
        // if ( !fdetector->get_OldW13eV() ) Variance /= sqrt ( ZurichEXOQ );

        double skewness;
        if ((yields.PhotonYield + yields.ElectronYield) > 1e4 ||
            yields.ElectricField > 4e3 || yields.ElectricField < 50.) {
            skewness = 0.00;  // make it a constant 0 when outside the range of Vetri
                            // Velan's Run04 models.
        } 
        else 
        {
            // LUX Skewness Model
            Wvalue wvalue = WorkFunction();
            double Wq_eV = wvalue.Wq_eV;
            double engy = 1e-3 * Wq_eV * (yields.PhotonYield + yields.ElectronYield);
            double fld = yields.ElectricField;

            double alpha0 = 1.39;
            double cc0 = 4.0, cc1 = 22.1;
            double E0 = 7.7, E1 = 54., E2 = 26.7, E3 = 6.4;
            double F0 = 225., F1 = 71.;

            skewness = 0.;

            if (ValidityTests::nearlyEqual(yields.Lindhard, 1.)) 
            {
                skewness = 1. / (1. + exp((engy - E2) / E3)) *
                        (alpha0 +
                        cc0 * exp(-1. * fld / F0) * (1. - exp(-1. * engy / E0))) +
                    1. / (1. + exp(-1. * (engy - E2) / E3)) * cc1 *
                        exp(-1. * engy / E1) * exp(-1. * sqrt(fld) / sqrt(F1));
            // if ( std::abs(skewness) <= DBL_MIN ) skewness = DBL_MIN;
            } 
            else {
                skewness = FreeParam[5];  // 2.25 but ~5-20 also good (for NR). All better
                                        // than zero, but 0 is OK too
            }  // note to self: find way to make 0 for ion (wall BG) incl. alphas?
        }

        double widthCorrection =
            sqrt(1. - (2. / M_PI) * skewness * skewness / (1. + skewness * skewness));
        double muCorrection = (
            sqrt(Variance) / widthCorrection) *
            (skewness / sqrt(1. + skewness * skewness)) * 2. *
            NEST::inv_sqrt2_PI;
        Ne = int(floor(
            RandomGen::rndm()->rand_gauss((1. - recombProb) * Ni, sqrt(Variance)) +
            0.5));
        
        Ne = NESTcalc::clamp(Ne, 0, Ni);

        Nph = Nq_actual - Ne;
        if (Nph > Nq_actual) Nph = Nq_actual;
        if (Nph < Nex) Nph = Nex;

        if ((Nph + Ne) != (Nex + Ni)) {
            throw std::runtime_error(
                "ERROR: Quanta not conserved. Tell Matthew Immediately!");
        }

        result.Variance = Variance;
        result.recombProb = recombProb;
        result.photons = Nph;
        result.electrons = Ne;

        return result;  // quanta returned with recomb fluctuations
    }

    YieldResult LArNEST::GetYieldNR(
        double energy, double density, double dfield,
        const std::vector<double>
            &NuisParam /*{11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/
    ) 
    {
        double massNum = fdetector->get_molarMass();

        if (NuisParam.size() < 12) 
        {
            throw std::runtime_error(
                "ERROR: You need a minimum of 12 nuisance parameters for the mean "
                "yields.");
        }
        if (energy > HIGH_E_NR)
            cerr << "\nWARNING: No data out here, you are beyond the AmBe endpoint of "
                    "about 300 keV.\n";
        int massNumber;
        double ScaleFactor[2] = {1., 1.};
        if (ValidityTests::nearlyEqual(massNum, 0.)) {
            massNumber = RandomGen::rndm()->SelectRanXeAtom();
        }
        else {
            massNumber = int(massNum);
        }
        ScaleFactor[0] = sqrt(fdetector->get_molarMass() / (double)massNumber);
        ScaleFactor[1] = ScaleFactor[0];
        double Nq = NuisParam[0] * pow(energy, NuisParam[1]);
        if (!fdetector->get_OldW13eV()) {
            Nq *= ZurichEXOW;
        }
        double ThomasImel =
            NuisParam[2] * pow(dfield, NuisParam[3]) * pow(density / DENSITY, 0.3);
        double Qy = 1. / (ThomasImel * pow(energy + NuisParam[4], NuisParam[9]));
        Qy *= 1. - 1. / pow(1. + pow((energy / NuisParam[5]), NuisParam[6]),
                            NuisParam[10]);
        if (!fdetector->get_OldW13eV()) {
            Qy *= ZurichEXOQ;
        }
        double Ly = Nq / energy - Qy;
        if (Qy < 0.0) {
            Qy = 0.0;
        }
        if (Ly < 0.0) {
            Ly = 0.0;
        }
        double Ne = Qy * energy * ScaleFactor[1];
        double Nph = Ly * energy * ScaleFactor[0] *
                    (1. - 1. / pow(1. + pow((energy / NuisParam[7]), NuisParam[8]),
                                    NuisParam[11]));
        Nq = Nph + Ne;
        double Ni = (4. / ThomasImel) * (exp(Ne * ThomasImel / 4.) - 1.);
        double Nex = (-1. / ThomasImel) *
                    (4. * exp(Ne * ThomasImel / 4.) - (Ne + Nph) * ThomasImel - 4.);

        double NexONi = Nex / Ni;
        Wvalue wvalue = WorkFunction();
        if (NexONi < wvalue.alpha && energy > 1e2) 
        {
            NexONi = wvalue.alpha;
            Ni = Nq / (1. + NexONi);
            Nex = Nq - Ni;
        }
        if (NexONi > 1.0 && energy < 1.) 
        {
            NexONi = 1.00;
            Ni = Nq / (1. + NexONi);
            Nex = Nq - Ni;
        }

        if (Nex < 0.) {
            cerr << "\nCAUTION: You are approaching the border of NEST's validity for "
                    "high-energy (OR, for LOW) NR, or are beyond it, at "
                << energy << " keV." << endl;
        }
        if (std::abs(Nex + Ni - Nq) > 2. * PHE_MIN) 
        {
            throw std::runtime_error(
                "ERROR: Quanta not conserved. Tell Matthew Immediately!");
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
        return YieldResultValidity(result, energy, Wq_eV);  
    }

    YieldResult LArNEST::GetYieldBeta(
        double energy, double density, double dfield
    ) 
    {  
        Wvalue wvalue = WorkFunction();
        double Qy, Nq;
        double Wq_eV = wvalue.Wq_eV;
        // double alpha = wvalue.alpha; // duplicate definition below. We don't even
        // need this here (it is Nex/Ni)

        // Liquid Argon
        double alpha =
            32.988 -
            552.988 / (17.2346 +
                    pow(dfield / (-4.7 + 0.025115 * exp(1.3954 / 0.265360653)),
                        0.242671));
        double beta = 0.778482 + 25.9 / pow(1.105 + pow(dfield / 0.4, 4.55), 7.502);
        double gamma =
            0.659509 *
            (1000 / 19.5 + 6.5 * (5 - 0.5 / pow(dfield / 1047.408, 0.01851)));
        double delta = 15.7489;
        double DB = 1052.264 + (14159350000 - 1652.264) /
                                (-5 + pow(dfield / 0.157933, 1.83894));
        double p1 = 1;
        double p2 = 10.304;
        double p3 = 13.0654;
        double p4 = 0.10535;
        double p5 = 0.7;
        double LET = -2.07763;
        Nq = energy * 1e3 / Wq_eV;
        Qy = alpha * beta +
            (gamma - alpha * beta) / pow(p1 + p2 * pow(energy + 0.5, p3), p4) +
            delta / (p5 + DB * pow(energy, LET));

        if (!fdetector->get_OldW13eV()) {
            Qy *= ZurichEXOQ;
        }
        double Ly = Nq / energy - Qy;
        double Ne = Qy * energy;
        double Nph = Ly * energy;

        YieldResult result{};
        result.PhotonYield = Nph;
        result.ElectronYield = Ne;
        result.ExcitonRatio = NexONi();
        result.Lindhard = 1;
        result.ElectricField = dfield;
        result.DeltaT_Scint = -999;
        return YieldResultValidity(result, energy, Wq_eV);  
    }

    YieldResult LArNEST::GetYields(
        INTERACTION_TYPE species, double energy, 
        double density, double dfield,
        double massNum, double atomNum,
        const std::vector<double> &NuisParam
        /*={11.,1.1,0.0480,-0.0533,12.6,0.3,2.,0.3,2.,0.5,1.,1.}*/,
        bool oldModelER
    ) 
    {
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
                return GetYieldNR(energy, density, dfield, NuisParam);
                // return GetYieldNROld ( energy, 1 );
                break;
            case ion:
                return GetYieldIon(energy, density, dfield, massNum, atomNum, NuisParam);
                break;
            case gammaRay:
                return GetYieldGamma(energy, density, dfield);
                break;
            case Kr83m:
                return GetYieldKr83m(energy, density, dfield, massNum, atomNum);
                // not actually massNumber, but a place holder for maxTime
                break;
            case fullGamma_PE:
                return GetYieldGamma(energy, density, dfield);  // PE of the full gamma spectrum
                break;
            default:  // beta, CH3T, 14C, the pp solar neutrino background, and
                    // Compton/PP spectra of fullGamma
                if (oldModelER) {
                    return GetYieldBeta(energy, density, dfield);  // OLD
                }
                break;
        }
    }

    std::vector<double> LArNEST::CalculateG2(bool verbosity) 
    {
        std::vector<double> g2_params(5);

        // Set parameters for calculating EL yield and extraction
        double alpha = 0.137, beta = 4.70e-18,
                gamma = 0;  // note the value of alpha is similar to ~1/7eV. Not
                            // coincidence. Noted in Mock et al.
        // actually listed as 'a' and 'b' in ref (below). Units 1/V, cm^2
        double epsilonRatio = EPS_LIQ / std::abs(EPS_GAS);
        if (fdetector->get_inGas()) {
            epsilonRatio = 1.;  // in an all-gas detector, E_liq variable below simply
        }
        // becomes the field value between anode and gate
        // Convert gas extraction field to liquid field
        double E_liq = fdetector->get_E_gas() / epsilonRatio;  // kV per cm
        double ExtEff =
                1. - 1.1974 * exp(-1.003 *
                                pow(E_liq, 1.3849));  // Guschin 1978-79 and 1981-82
        if (ExtEff > 1. || fdetector->get_inGas()) ExtEff = 1.;
        if (ExtEff < 0. || E_liq <= 0.) ExtEff = 0.;

        double gasGap =
            fdetector->get_anode() -
            fdetector
                ->get_TopDrift();  // EL gap in mm -> cm, affecting S2 size linearly
        if (gasGap <= 0. && E_liq > 0.) {
            throw std::runtime_error(
                "\tERR: The gas gap in the S2 calculation broke!!!!");
        }

        // Calculate EL yield based on gas gap, extraction field, and pressure
        // double elYield = (alpha * fdetector->get_E_gas() * 1e3 -
        //                beta * fdetector->get_p_bar() - gamma) *
        //               gasGap * 0.1;  // arXiv:1207.2292 (HA, Vitaly C.)
        bool YesGas = true;
        double rho = GetDensity(
            fdetector->get_T_Kelvin(), fdetector->get_p_bar(),
            YesGas, 1, fdetector->get_molarMass()
        );
        double elYield;
        // Henrique Araujo and Vitaly Chepel again
        alpha = 0.0813;
        beta = 1.90e-18;
        elYield = (alpha * fdetector->get_E_gas() * 1e3 -
                    beta * (NEST_AVO * rho / fdetector->get_molarMass())) *
                    gasGap * 0.1;
        // replaced with more accurate version also from 1207.2292, but works for room
        // temperature gas
        if (elYield <= 0.0 && E_liq != 0.) 
        {
            cerr << "\tWARNING, the field in gas must be at least "
                << 1e-3 * beta * NEST_AVO * rho / (alpha * fdetector->get_molarMass())
                << " kV/cm, for S2 to work," << endl;
            cerr << "\tOR: your density for gas must be less than "
                << fdetector->get_molarMass() * alpha * fdetector->get_E_gas() * 1e3 /
                        (beta * NEST_AVO)
                << " g/cm^3." << endl;
        }
        // Calculate single electron size and then g2
        double SE = elYield * fdetector->get_g1_gas();  // multiplying by light
        // collection efficiency in
        // the gas gap
        if (fdetector->get_s2_thr() < 0) {
            SE *= fdetector->FitTBA(0., 0., fdetector->get_TopDrift() / 2.)[1];
        }
        double g2 = ExtEff * SE;
        double StdDev = 0., Nphe, pulseArea, pulseAreaC, NphdC, phi, posDep, r, x, y;
        int Nph, nHits;
        if (verbosity) 
        {
            for (int i = 0; i < 10000; ++i) {
            // calculate properly the width (1-sigma std dev) in the SE size
            Nph = int(floor(RandomGen::rndm()->rand_gauss(
                                elYield, sqrt(fdetector->get_s2Fano() * elYield), true) +
                            0.5));
            phi = NEST::two_PI * RandomGen::rndm()->rand_uniform();
            r = fdetector->get_radius() * sqrt(RandomGen::rndm()->rand_uniform());
            x = r * cos(phi);
            y = r * sin(phi);
            posDep = fdetector->FitS2(x, y, VDetector::fold) /
                    fdetector->FitS2(
                        0., 0., VDetector::fold);  // future upgrade: smeared pos
            nHits =
                RandomGen::rndm()->binom_draw(Nph, fdetector->get_g1_gas() * posDep);
            Nphe =
                nHits + RandomGen::rndm()->binom_draw(nHits, fdetector->get_P_dphe());
            pulseArea = RandomGen::rndm()->rand_gauss(
                Nphe, fdetector->get_sPEres() * sqrt(Nphe));
            if (fdetector->get_noiseQuadratic()[1] != 0) {
                pulseArea = RandomGen::rndm()->rand_gauss(
                    pulseArea,
                    sqrt(
                        pow(fdetector->get_noiseQuadratic()[1] * pow(pulseArea, 2), 2) +
                        pow(fdetector->get_noiseLinear()[1] * pulseArea, 2)));
            } else {
                pulseArea = RandomGen::rndm()->rand_gauss(
                    pulseArea, fdetector->get_noiseLinear()[1] * pulseArea);
            }
            if (fdetector->get_s2_thr() < 0.)
                pulseArea = RandomGen::rndm()->rand_gauss(
                    fdetector->FitTBA(0.0, 0.0, fdetector->get_TopDrift() / 2.)[1] *
                        pulseArea,
                    sqrt(
                        fdetector->FitTBA(0.0, 0.0, fdetector->get_TopDrift() / 2.)[1] *
                        pulseArea *
                        (1. - fdetector->FitTBA(0.0, 0.0,
                                                fdetector->get_TopDrift() / 2.)[1])));
            pulseAreaC = pulseArea / posDep;
            NphdC = pulseAreaC / (1. + fdetector->get_P_dphe());
            StdDev += (SE - NphdC) * (SE - NphdC);
            }
            StdDev = sqrt(StdDev) / sqrt(9999.);  // N-1 from above (10,000)

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

    QuantaResult LArNEST::LegacyCalculation(
        int pdgcode, double energy,
        double density, double eField, 
        double track_length
    )
    {
        // default quenching factor, for electronic recoils
        double yieldFactor = 1.0; 
        // ratio for light particle in LAr, such as e-, mu-, Aprile et. al book
        double excitationRatio = 0.21; 
        // set up DokeBirks coefficients
        double DokeBirks[3];
        if (eField) {
            DokeBirks[0] = 0.07 * pow((eField / 1.0e3), -0.85);
            DokeBirks[2] = 0.00;
        }
        else {
            DokeBirks[0] = 0.0003;
            DokeBirks[2] = 0.75;
        }
        DokeBirks[1] = DokeBirks[0] / (1 - DokeBirks[2]); //B=A/(1-C) (see paper)
        // nuclear recoil quenching "L" factor: total yield is
        // reduced for nuclear recoil as per Lindhard theory
        double epsilon = 11.5 * (energy / 1.e-3 * pow(LAr_Z, (-7. / 3.)));

        if (abs(pdgcode) == 2112) //nuclear recoil
        {
            yieldFactor = 0.23 * (1 + exp(-5 * epsilon)); //liquid argon L_eff
            excitationRatio = 0.69337 + 0.3065 * exp(-0.008806 * pow(eField, 0.76313));
        }

        // determine ultimate number of quanta from current E-deposition (ph+e-) 
        // total mean number of exc/ions the total number of either quanta produced 
        // is equal to product of the work function, the energy deposited, 
        // and yield reduction, for NR
        double MeanNumQuanta = legacy_scint_yield * energy;
        double sigma = sqrt(legacy_resolution_scale * MeanNumQuanta); //Fano
        int NumQuanta = int(floor(RandomGen::rndm()->rand_gauss(MeanNumQuanta, sigma) + 0.5));
        double LeffVar = RandomGen::rndm()->rand_gauss(yieldFactor, 0.25 * yieldFactor);
        LeffVar = std::clamp(LeffVar, 0., 1.);

        if (yieldFactor < 1) // nuclear recoil
        {
            NumQuanta = RandomGen::rndm()->binom_draw(NumQuanta, LeffVar);
        }

        //if Edep below work function, can't make any quanta, and if NumQuanta
        //less than zero because Gaussian fluctuated low, update to zero
        if (energy < 1 / legacy_scint_yield || NumQuanta < 0) { NumQuanta = 0; }

        // next section binomially assigns quanta to excitons and ions
        int NumExcitons = RandomGen::rndm()->binom_draw(
            NumQuanta, excitationRatio / (1 + excitationRatio)
        );
        int NumIons = NumQuanta - NumExcitons;

        // this section calculates recombination following the modified 
        // Birks'Law of Doke, deposition by deposition, may be overridden 
        // later in code if a low enough energy necessitates switching to the
        // Thomas-Imel box model for recombination instead (determined by site)
        double dE = energy;
        double dx = 0.0;
        double LET = 0.0;
        double recombProb;

        if (abs(pdgcode) != 11 && abs(pdgcode) != 13) 
        {
            // e-: 11, e+: -11, mu-: 13, mu+: -13
            // in other words, if it's a gamma,ion,proton,alpha,pion,et al. do not
            // use the step length provided by Geant4 because it's not relevant,
            // instead calculate an estimated LET and range of the electrons that
            // would have been produced if Geant4 could track them
            LET = LegacyCalcElectronLET(1000 * dE);

            if (LET) {
                //find the range based on the LET
                dx = dE / (density * LET); 
            }
            if (abs(pdgcode) == 2112) // nuclear recoils
            {
                dx = 0;
            }
        }
        else //normal case of an e-/+ energy deposition recorded by Geant
        {
            dx = track_length / 10.;
            if (dx) 
            {
                LET = (dE / dx) * (1 / density); //lin. energy xfer (prop. to dE/dx)
            }
            if (LET > 0 && dE > 0 && dx > 0) 
            {
                double ratio = LegacyCalcElectronLET(dE * 1e3) / LET;
                if (ratio < 0.7 && pdgcode == 11) 
                {
                    dx /= ratio;
                    LET *= ratio;
                }
            }
        }
        recombProb = (DokeBirks[0] * LET) / (1 + DokeBirks[1] * LET) +
                    DokeBirks[2]; //Doke/Birks' Law as spelled out in the NEST pape
        recombProb *= (density / legacy_density_LAr);

        //check against unphysicality resulting from rounding errors
        recombProb = std::clamp(recombProb, 0., 1.);

        //use binomial distribution to assign photons, electrons, where photons
        //are excitons plus recombined ionization electrons, while final
        //collected electrons are the "escape" (non-recombined) electrons
        int const NumPhotons = NumExcitons + RandomGen::rndm()->binom_draw(NumIons, recombProb);
        int const NumElectrons = NumQuanta - NumPhotons;

        // create the quanta results
        QuantaResult result;
        result.photons = NumPhotons;
        result.electrons = NumElectrons;
        
        return result;

    }
    double LArNEST::LegacyCalcElectronLET(double E)
    {
        double LET;
        if (E >= 1) 
        {
            LET = 116.70 - 
                  162.97 * log10(E) + 
                  99.361 * pow(log10(E), 2) - 
                  33.405 * pow(log10(E), 3) +
                  6.5069 * pow(log10(E), 4) - 
                  0.69334 * pow(log10(E), 5) + 
                  0.031563 * pow(log10(E), 6);
        }
        else if (E > 0 && E < 1) 
        {
            LET = 100;
        }
        else 
        {
            LET = 0;
        }
        return LET;
    }
}