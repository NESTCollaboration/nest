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

    //------------------------------------NR Yields------------------------------------//
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
        return fLArNRYieldsParameters.alpha * pow(energy, fLArNRYieldsParameters.beta - 1.0) - 
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
    //------------------------------------ER Yields------------------------------------//
    // TODO: This is likely temporary, otherwise it says that the total
    // yield is linear with energy, but need to incorporate thermal quenching.
    double LArNEST::GetERTotalYields(double energy)
    {
        return (energy * 1.0e3 / fWorkQuantaFunction);
    }
    double LArNEST::GetERExcitonYieldsAlpha(double efield, double density)
    {
        return fLArERYieldsParameters.alpha.A + 
               fLArERYieldsParameters.alpha.B /
               (fLArERYieldsParameters.alpha.C + 
                   pow(efield / (fLArERYieldsParameters.alpha.D + fLArERYieldsParameters.alpha.E * 
                             exp(density / fLArERYieldsParameters.alpha.F)), fLArERYieldsParameters.alpha.G));
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
               (fLArERYieldsParameters.gamma.B / fWorkQuantaFunction + 
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
               fLArERYieldsParameters.delta / (fLArERYieldsParameters.p5 + doke_birks * pow(energy, fLArERYieldsParameters.let));
    }
    LArYieldResult LArNEST::GetNRYields(
        double energy, double efield, double density    
    ) 
    {
        double NRTotalYields = GetNRTotalYields(energy);
        double NRExcitonYields = GetNRExcitonYields(energy, efield);
        if (NRExcitonYields < 0.0) {
            NRExcitonYields = 0.0;
        }
        //double NRPhotonYields = NRTotalYields / energy - NRExcitonYields;
        double NRPhotonYields = GetNRPhotonYields(energy, efield);
        if (NRPhotonYields < 0.0) {
            NRPhotonYields = 0.0;
        }

        double Ne = NRExcitonYields * energy;
        double Nph = NRPhotonYields * energy;
        double Nq = Ne + Nph;

        /// recombination
        double ThomasImel = fThomasImelParameters.A * pow(efield, fThomasImelParameters.B);
        double Nion = (4. / ThomasImel) * (exp(Ne * ThomasImel / 4.) - 1.);
        double Nex = Nq - Nion;

        // double Nex = (Ne + Nph)/(1.0 + fNexOverNion);
        // double Nion = (Ne + Nph) - Nex;

        double WQ_eV = fWorkQuantaFunction;
        double Lindhard = (NRTotalYields / energy) * WQ_eV * 1e-3;
        LArYieldResult result{
            NRTotalYields, NRExcitonYields, NRPhotonYields,
            Nph, Ne, Nex, Nion, Lindhard, efield
        };
        return result;
        //return YieldResultValidity(result, energy, WQ_eV);  
    }

    LArYieldResult LArNEST::GetERYields(
        double energy, double density, double efield
    ) 
    {  
        double ERTotalYields = GetERTotalYields(energy);
        double ERExcitonYields = GetERExcitonYields(energy, efield, density);
        if (ERExcitonYields < 0.0) {
            ERExcitonYields = 0.0;
        }
        double ERPhotonYields = ERTotalYields / energy - ERExcitonYields;
        if (ERPhotonYields < 0.0) {
            ERPhotonYields = 0.0;
        }

        double Ne = ERExcitonYields * energy;
        double Nph = ERPhotonYields * energy;
        double Nq = Ne + Nph;

        /// recombination
        double ThomasImel = fThomasImelParameters.A * pow(efield, fThomasImelParameters.B);
        double Nion = (4. / ThomasImel) * (exp(Ne * ThomasImel / 4.) - 1.);
        double Nex = Nq - Nion;

        double WQ_eV = fWorkQuantaFunction;
        double Lindhard = (ERTotalYields / energy) * WQ_eV * 1e-3;
        LArYieldResult result{
            ERTotalYields, ERExcitonYields, ERPhotonYields,
            Nph, Ne, Nex, Nion, Lindhard, efield
        };
        return result;
        //return YieldResultValidity(result, energy, WQ_eV);
    }
    //------------------------------------Full Calculations------------------------------------//
    LArYieldResult LArNEST::GetYields(
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
                    //return GetYieldBeta(energy, density, efield);  // OLD
                }
                break;
        }
    }

    LArYieldFluctuationResult LArNEST::GetYieldFluctuations(
        const YieldResult &yields, double density
    ) 
    {
        // TODO: Understand and break up this function.
        LArYieldFluctuationResult result{};
        return result;  // quanta returned with recomb fluctuations
    }

    LArNESTResult LArNEST::FullCalculation(
        INTERACTION_TYPE species, double energy, 
        double density, double efield,
        bool do_times
    ) 
    {
        if (density < 1.) fdetector->set_inGas(true);
        LArNESTResult result;
        result.yields = GetYields(species, energy, density, efield);
        // result.quanta = GetQuanta(result.yields, density);
        // if (do_times)
        //     result.photon_times = GetPhotonTimes(
        //         species, result.quanta.photons, 
        //         result.quanta.excitons, efield, energy
        //     );
        // else {
        //     result.photon_times = photonstream(result.quanta.photons, 0.0);
        // }
        return result;
    }
    
    //------------------------------------Various Functions------------------------------------//
    std::vector<double> LArNEST::CalculateG2(bool verbosity) 
    {
    }

    //------------------------------------Legacy NEST------------------------------------//
    LArYieldResult LArNEST::LegacyCalculation(
        int pdgcode, double energy,
        double density, double efield, 
        double track_length
    )
    {
        // default quenching factor, for electronic recoils
        double yieldFactor = 1.0; 
        // ratio for light particle in LAr, such as e-, mu-, Aprile et. al book
        double excitationRatio = 0.21; 
        // set up DokeBirks coefficients
        double DokeBirks[3];
        if (efield) {
            DokeBirks[0] = 0.07 * pow((efield / 1.0e3), -0.85);
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
            excitationRatio = 0.69337 + 0.3065 * exp(-0.008806 * pow(efield, 0.76313));
        }

        // determine ultimate number of quanta from current E-deposition (ph+e-) 
        // total mean number of exc/ions the total number of either quanta produced 
        // is equal to product of the work function, the energy deposited, 
        // and yield reduction, for NR
        double MeanNq = legacy_scint_yield * energy;
        double sigma = sqrt(legacy_resolution_scale * MeanNq); //Fano
        int Nq = int(floor(RandomGen::rndm()->rand_gauss(MeanNq, sigma) + 0.5));
        double LeffVar = RandomGen::rndm()->rand_gauss(yieldFactor, 0.25 * yieldFactor);
        LeffVar = clamp(LeffVar, 0., 1.);

        if (yieldFactor < 1) // nuclear recoil
        {
            Nq = RandomGen::rndm()->binom_draw(Nq, LeffVar);
        }

        //if Edep below work function, can't make any quanta, and if Nq
        //less than zero because Gaussian fluctuated low, update to zero
        if (energy < 1 / legacy_scint_yield || Nq < 0) { Nq = 0; }

        // next section binomially assigns quanta to excitons and ions
        double Nex = RandomGen::rndm()->binom_draw(
            Nq, excitationRatio / (1 + excitationRatio)
        );
        double Nion = Nq - Nex;

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
        recombProb = clamp(recombProb, 0., 1.);

        //use binomial distribution to assign photons, electrons, where photons
        //are excitons plus recombined ionization electrons, while final
        //collected electrons are the "escape" (non-recombined) electrons
        double Nph = Nex + RandomGen::rndm()->binom_draw(Nion, recombProb);
        double Ne = Nq - Nph;

        // create the quanta results
        LArYieldResult result {
            (Ne + Nph)/energy, Ne/energy, Nph/energy,
            Nph, Ne, Nex, Nion, 0.0, efield
        };
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