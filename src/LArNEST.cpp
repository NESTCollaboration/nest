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

        DokeBirks[1] = DokeBirks[0] / (1 - DokeBirks[2]); //B=A/(1-C) (see paper)
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