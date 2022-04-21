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
        double density, double efield,
        bool do_times
    ) 
    {
        if (density < 1.) fdetector->set_inGas(true);
        NESTresult result;
        result.yields = GetYields(species, energy, density, efield);
        result.quanta = GetQuanta(result.yields, density);
        if (do_times)
            result.photon_times = GetPhotonTimes(
                species, result.quanta.photons, 
                result.quanta.excitons, efield, energy
            );
        else {
            result.photon_times = photonstream(result.quanta.photons, 0.0);
        }
        return result;
    }
    double LArNEST::GetNRTotalYields(double energy)
    {
        return fLArNRYieldsParameters.alpha * pow(energy, fLArNRYieldsParameters.beta);
    }
    double LArNEST::GetNRExcitonYields(double energy, double efield)
    {
        return (1.0 / (fLArNRYieldsParameters.gamma * pow(efield, fLArNRYieldsParameters.delta))) *
               (1.0 / sqrt(energy + fLArNRYieldsParameters.epsilon)) *
               (1.0 - 1.0 / (1.0 + pow(energy/fLArNRYieldsParameters.zeta, fLArNRYieldsParameters.eta)));
    }
    double LArNEST::GetNRPhotonYields(double energy, double efield)
    {
        return fLArNRYieldsParameters.alpha * pow(energy, fLArNRYieldsParameters.beta - 1) - 
               (1.0 / (fLArNRYieldsParameters.gamma * pow(efield, fLArNRYieldsParameters.delta))) *
               (1.0 / sqrt(energy + fLArNRYieldsParameters.epsilon));
    }
    double LArNEST::GetNRPhotonYieldsConserved(double energy, double efield)
    {
        return fLArNRYieldsParameters.alpha * pow(energy, fLArNRYieldsParameters.beta) - 
               (1.0 / (fLArNRYieldsParameters.gamma * pow(efield, fLArNRYieldsParameters.delta))) *
               (1.0 / sqrt(energy + fLArNRYieldsParameters.epsilon)) *
               (1.0 - 1.0 / (1.0 + pow(energy/fLArNRYieldsParameters.zeta, fLArNRYieldsParameters.eta)));
    }
    // TODO: This is likely temporary, otherwise it says that the total
    // yield is linear with energy, but need to incorporate thermal quenching.
    double LArNEST::GetERTotalYields(double energy)
    {
        return energy *1e3 / GetWorkFunction();
    }
    double LArNEST::GetERExcitonYieldsAlpha(double efield, double density)
    {
        return fLArERYieldsParameters.alpha.A + 
               fLArERYieldsParameters.alpha.B *
               pow(fLArERYieldsParameters.alpha.C + 
                   efield / (fLArERYieldsParameters.alpha.D + fLArERYieldsParameters.alpha.E * 
                             exp(fLArERYieldsParameters.alpha.F * density)), fLArERYieldsParameters.alpha.G);
    }
    double LArNEST::GetERExcitonYieldsBeta(double efield)
    {
        return fLArERYieldsParameters.beta.A +         
               fLArERYieldsParameters.beta.B *
               pow(fLArERYieldsParameters.beta.C + 
                   pow(efield / fLArERYieldsParameters.beta.D, fLArERYieldsParameters.beta.E), fLArERYieldsParameters.beta.F);
    }
    double LArNEST::GetERExcitonYieldsGamma(double efield)
    {
        return fLArERYieldsParameters.gamma.A *
               (fLArERYieldsParameters.gamma.B / GetWorkFunction() + 
                fLArERYieldsParameters.gamma.C *
                (fLArERYieldsParameters.gamma.D + fLArERYieldsParameters.gamma.E / 
                 pow(efield / fLArERYieldsParameters.gamma.F, fLArERYieldsParameters.gamma.G)));
    }
    double LArNEST::GetERExcitonYieldsDokeBirks(double efield)
    {
        return fLArERYieldsParameters.doke_birks.A + 
               fLArERYieldsParameters.doke_birks.B / 
               (fLArERYieldsParameters.doke_birks.C + 
                pow(efield / fLArERYieldsParameters.doke_birks.D, fLArERYieldsParameters.doke_birks.E));
    }
    double LArNEST::GetERExcitonYields(double energy, double efield, double density)
    {
        double alpha = GetERExcitonYieldsAlpha(efield, density);
        double beta = GetERExcitonYieldsBeta(efield);
        double gamma = GetERExcitonYieldsGamma(efield);
        double doke_birks = GetERExcitonYieldsDokeBirks(efield);
        return alpha * beta + 
               (gamma - alpha * beta) / 
               pow(fLArERYieldsParameters.p1 + fLArERYieldsParameters.p2 * pow(energy + 0.5, fLArERYieldsParameters.p3), fLArERYieldsParameters.p4) + 
               fLArERYieldsParameters.delta / (fLArERYieldsParameters.p5 + pow(doke_birks, fLArERYieldsParameters.let));
    }
    double LArNEST::GetWorkFunction()
    {
        return 1000. / 51.9;
    }
    inline double LArNEST::GetFanoER() 
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
        double VaporP_bar = exp(45.973940 - 1464.718291 / Kelvin - 6.539938 * log(Kelvin));
        if (bara < VaporP_bar || inGas) 
        {
            double p_Pa = bara * 1e5;
            double density =
                1. /
                (pow(fRIdealGas * Kelvin, 3.) /
                    (p_Pa * pow(fRIdealGas * Kelvin, 2.) + fRealGasA * p_Pa * p_Pa) +
                fRealGasB) * molarMass * 1e-6;  // Van der Waals equation, mol/m^3
            if (bara < VaporP_bar && evtNum == 0)
            // TODO: get rid of these output statements
                std::cerr << "\nWARNING: ARGON GAS PHASE. IS THAT WHAT YOU WANTED?\n";
            inGas = true;
            return density;
        }
        else 
        {
            inGas = false;
            return 1.4;
        }
    }
    double LArNEST::GetDriftVelocity_Liquid(
        double Kelvin, double eField
    ) 
    {  
        // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
        double speed = 0.0;  
        if (Kelvin < 84. || Kelvin > 140.) 
        {
            cerr << "\nWARNING: TEMPERATURE OUT OF RANGE (84-140 K) for vD\n";
            cerr << "Using value at closest temp for a drift speed estimate\n";
            Kelvin = (double)NESTcalc::clamp(int(Kelvin), 84, 140);
        }

        if (Kelvin >= fTemperature[0] && Kelvin < fTemperature[1]) {
            speed = exp(0.937729 - 0.0734108 / (eField / 1000) +
                        0.315338 * log(eField / 1000));
        }
        else if (Kelvin >= fTemperature[1] && Kelvin < fTemperature[2]) {
            speed = exp(0.80302379 - 0.06694564 / (eField / 1000) +
                        0.331798 * log(eField / 1000));
        }
        else if (Kelvin >= fTemperature[2] && Kelvin < fTemperature[3]) {
            speed = exp(0.7795972 - 0.0990952 / (eField / 1000) +
                        0.320876 * log(eField / 1000));
        }
        else if (Kelvin >= fTemperature[3] && Kelvin < fTemperature[4]) {
            speed = exp(0.6911897 - 0.092997 / (eField / 1000) +
                        0.3295202 * log(eField / 1000));
        }
        else if (Kelvin >= fTemperature[4] && Kelvin < fTemperature[5]) {
            speed = exp(0.76551511 - 0.0731659 / (eField / 1000) +
                        0.317972 * log(eField / 1000));
        }
        else if (Kelvin >= fTemperature[5] && Kelvin < fTemperature[6]) {
            speed = exp(0.502022794 - 0.06644517 / (eField / 1000) +
                        0.3290246 * log(eField / 1000));
        }
        else if (Kelvin >= fTemperature[6] && Kelvin <= fTemperature[7]) {
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
        // TODO: What is going on here?  molarMass is an input but
        // is being overwritten.
        molarMass = 39.948;
        density *= NEST_AVO / molarMass;
        double edrift = 0., gasdep = efieldinput / density;
        edrift = 2.991205 * pow(1.00113657, 1. / (gasdep * 1.e17)) *
                pow(gasdep * 1.e17, 0.2570253);
        return edrift;
    }
    double LArNEST::GetPhotonEnergy(bool state) 
    {
        // liquid Argon
        // TODO: What is this function about?
        if (state) {
            return RandomGen::rndm()->rand_gauss(9.7, 0.2);
        }
        else {
            return RandomGen::rndm()->rand_gauss(9.69, 0.22);
        }
    }
    inline double LArNEST::GetNexONi() 
    {
        // https://arxiv.org/pdf/1903.05815.pdf
        // TODO: What is this?
        return 0.21;  
    }
    double LArNEST::GetLinearEnergyTransfer(double energy, bool CSDA) 
    {
        double LET;
        if (!CSDA) 
        { //total stopping power directly from ESTAR (radiative + collision)
            LET = 1.8106 - 
                  0.45086 * log10(energy) - 
                  0.33151 * pow(log10(energy),2.) + 
                  0.25916 * pow(log10(energy),3.) - 
                  0.2051  * pow(log10(energy),4.) + 
                  0.15279 * pow(log10(energy),5.) - 
                  0.084659 * pow(log10(energy),6.) + 
                  0.030441 * pow(log10(energy),7.) - 
                  0.0058953 * pow(log10(energy),8.) + 
                  0.00045633 * pow(log10(energy),9.);
            LET = pow(10.,LET);
            if ( std::isnan(LET) || LET <= 0. ) {
                LET = 1e2;
            }
        }
        else 
        {
            // the "continuous slowing down approximation" (CSDA)
            // replace with Justin and Prof. Mooney's work
            if (energy >= 1.) 
            {
                LET = 116.70 - 
                      162.97 * log10(energy) + 
                      99.361 * pow(log10(energy), 2) -
                      33.405 * pow(log10(energy), 3) + 
                      6.5069 * pow(log10(energy), 4) -
                      0.69334 * pow(log10(energy), 5) + 
                      0.031563 * pow(log10(energy), 6);
            }
            else if (energy > 0. && energy < 1.) {
                LET = 100.;
            } 
            else {
                LET = 0.;
            }
        }
        return LET;
    }

    double LArNEST::GetPhotonTime(
        INTERACTION_TYPE species, bool exciton,
        double efield, double energy
    ) 
    {
        // arXiv:1802.06162. NR may need tauR ~0.5-1ns instead of 0
        double SingTripRatio = 0.; 
        // if ( efield > 60. ) efield = 60. // makes Xed work. 200 for LUX Run04
        // instead. A mystery! Why no field dep?

        double LET = GetLinearEnergyTransfer(energy);
        // copied from 2013 NEST version for LAr on LBNE
        double tau1 = RandomGen::rndm()->rand_gauss(6.5, 0.8);  // error from weighted average
        double tau3 = RandomGen::rndm()->rand_gauss(1300, 50);  // ibid.
        double tauR = RandomGen::rndm()->rand_gauss(0.8, 0.2);  // Kubota 1979
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
        else
        {
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
                if (LET < 3. && !exciton) {
                    SingTripRatio = RandomGen::rndm()->rand_gauss(0.5, 0.2);
                }
                else if (LET < 3. && exciton) {
                    SingTripRatio = RandomGen::rndm()->rand_gauss(0.36, 0.06);
                }
                else 
                {
                    SingTripRatio = 0.2701 + 
                                    0.003379 * LET - 
                                    4.7338e-5 * pow(LET, 2.) +
                                    8.1449e-6 * pow(LET, 3.);
                }
            }  // lastly is ER for LAr
        }
        if (tauR < 0.) {
            tauR = 0.;  
        } // in case varied with Gaussian earlier
        // the recombination time is non-exponential, but approximates
        // to exp at long timescales (see Kubota '79)
        double time_ns = tauR * (1.0 / RandomGen::rndm()->rand_uniform() - 1.);

        if (RandomGen::rndm()->rand_uniform() < SingTripRatio / (1. + SingTripRatio)) {
            time_ns -= tau1 * log(RandomGen::rndm()->rand_uniform());
        }
        else {
            time_ns -= tau3 * log(RandomGen::rndm()->rand_uniform());
        }
        return time_ns;
    }

    YieldResult LArNEST::GetNRYields(
        double energy, double efield, double density    
    ) 
    {
        // TODO: Fix these error messages
        // if (energy > HIGH_E_NR)
        //     cerr << "\nWARNING: No data out here, you are beyond the AmBe endpoint of "
        //             "about 300 keV.\n";

        double NRTotalYields = GetNRTotalYields(energy);
        double NRExcitonYields = GetNRExcitonYields(energy, efield);
        double NRPhotonYields = GetNRPhotonYields(energy, efield);

        if (NRExcitonYields < 0.0) {
            NRExcitonYields = 0.0;
        }
        // Perhaps don't need this unless energy is negative
        // if (NRPhotonYields < 0.0) {
        //     NRPhotonYields = 0.0;
        // }
        double Ne = NRExcitonYields * energy;
        double Nph = NRPhotonYields * energy;

        double Wq_eV = GetWorkFunction();
        double L = (NRTotalYields / energy) * Wq_eV * 1e-3;
        YieldResult result{
            Nph, Ne, GetNexONi(), L, efield, -999
        };
        return YieldResultValidity(result, energy, Wq_eV);  
    }

    YieldResult LArNEST::GetERYields(
        double energy, double density, double efield
    ) 
    {  
        double ERTotalYields = GetERTotalYields(energy);
        double ERExcitonYields = GetERExcitonYields(energy, efield, density);
        double ERPhotonYields = ERTotalYields - ERExcitonYields;

        if (ERExcitonYields < 0.0) {
            ERExcitonYields = 0.0;
        }
        // Perhaps don't need this unless energy is negative
        // if (ERPhotonYields < 0.0) {
        //     ERPhotonYields = 0.0;
        // }
        double Ne = ERExcitonYields * energy;
        double Nph = ERPhotonYields * energy;

        double Wq_eV = GetWorkFunction();
        double L = (ERTotalYields / energy) * Wq_eV * 1e-3;

        YieldResult result{
            Nph, Ne, GetNexONi(), L, efield, -999
        };
        return YieldResultValidity(result, energy, Wq_eV);  
    }

    YieldResult LArNEST::GetYields(
        INTERACTION_TYPE species, double energy, 
        double density, double efield,
        bool oldModelER
    ) 
    {
        switch (species) {
            case NR:
                return GetNRYields(energy, efield, density);
                break;
            case WIMP:
            case B8:
            case atmNu:
            case DD:
            case AmBe:
            case Cf:  
                return GetNRYields(energy, efield, density);
                break;
            case ion:
                return GetERYields(energy, density, efield);
                break;
            case gammaRay:
            case Kr83m:
            case fullGamma_PE:
            default:  // beta, CH3T, 14C, the pp solar neutrino background, and
                    // Compton/PP spectra of fullGamma
                if (oldModelER) {
                    return GetYieldBeta(energy, density, efield);  // OLD
                }
                break;
        }
    }
    //  inline LArNEST::Wvalue LArNEST::WorkFunction() 
    // {
    //     // Liquid argon
    //     // TODO: Is this work function independent of
    //     // other parameters?
    //     double alpha = 0.21;
    //     double Wq_eV = 1000. / 51.9;  // 23.6/1.21; // ~19.2-5 eV
    //     return Wvalue{.Wq_eV = Wq_eV, .alpha = alpha};
    // }
    QuantaResult LArNEST::GetQuanta(
        const YieldResult &yields, double density
        /*={1.,1.,0.1,0.5,0.19,2.25}*/
    ) 
    {
        // TODO: Understand and break up this function.
        QuantaResult result{};
        return result;  // quanta returned with recomb fluctuations
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
            LET = LegacyGetLinearEnergyTransfer(1000 * dE);

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
                double ratio = LegacyGetLinearEnergyTransfer(dE * 1e3) / LET;
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
    double LArNEST::LegacyGetLinearEnergyTransfer(double E)
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