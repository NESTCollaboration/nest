/**
 * @file LArNEST.cpp
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Justin Mueller [Justin.Mueller@colostate.edu]
 * @author Ekaterina Kozlova [aspelene@gmail.com]
 * @author Michael Mooney [mrmooney@colostate.edu]
 * @brief
 * @version
 * @date 2022-04-13
 */
#include "LArNEST.hh"

namespace NEST {
LArNEST::LArNEST(VDetector *detector) : NESTcalc(detector) {}
//-------------------------Full Calculation-------------------------//
LArNESTResult LArNEST::FullCalculation(LArInteraction species, double energy,
                                       double dx, double efield, double density,
                                       bool do_times) {
  if (density < 1.) fdetector->set_inGas(true);
  LArNESTResult result;
  result.yields = GetYields(species, energy, dx, efield, density);
  result.fluctuations = GetYieldFluctuations(species, result.yields, density);
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
//-------------------------All Yields-------------------------//
LArYieldResult LArNEST::GetYields(LArInteraction species, double energy,
                                  double dx, double efield, double density) {
  if (species == LArInteraction::NR) {
    return GetNRYields(energy, efield, density);
  } else if (species == LArInteraction::ER) {
    return GetERYields(energy, efield, density);
  } else if (species == LArInteraction::Alpha) {
    return GetAlphaYields(energy, efield, density);
  } else if (species == LArInteraction::LeptonLET) {
    return GetLeptonLETYields(energy, dx, efield, density);
  } else if (species == LArInteraction::LET) {
    return GetLETYields(energy, dx, efield, density);
  } else {
    return GetdEdxYields(energy, dx, efield, density);
  }
}
//-------------------------Yields Fluctuations-------------------------//
LArYieldFluctuationResult LArNEST::GetYieldFluctuations(
    LArInteraction species, const LArYieldResult &yields, double density) {
  if (species == LArInteraction::NR || species == LArInteraction::ER ||
      species == LArInteraction::Alpha) {
    return GetDefaultFluctuations(yields, density);
  } else {
    return GetDefaultFluctuations(yields, density);
  }
}
//-------------------------Recombination Yields-------------------------//
LArYieldResult LArNEST::GetRecombinationYields(double TotalYields,
                                               double ElectronYields,
                                               double PhotonYields,
                                               double energy, double efield) {
  double Ne = ElectronYields * energy;
  double Nph = PhotonYields * energy;
  double Nq = Ne + Nph;

  /// recombination
  // Thomas-Imel is not used here, since alpha = Nex/Nion
  // is assumed to be a fixed quantity, and so the
  // recombination probability is constrained through
  // the yields functions:
  //      Nion = Nq/(1 + alpha) = (Ye + Yph)E/(1 + alpha)
  double Nion = Nq / (1 + fNexOverNion);
  double Nex = Nq - Nion;

  double WQ_eV = fWorkQuantaFunction;
  double Lindhard = (TotalYields / energy) * WQ_eV * 1e-3;
  LArYieldResult result{TotalYields, ElectronYields, PhotonYields, Nph,   Ne,
                        Nex,         Nion,           Lindhard,     efield};
  return result;
}
//-------------------------NR Yields-------------------------//
double LArNEST::GetNRTotalYields(double energy) {
  return fNR.alpha * pow(energy, fNR.beta);
}
double LArNEST::GetNRElectronYields(double energy, double efield) {
  return (1.0 / (fNR.gamma * pow(efield, fNR.delta))) *
         (1.0 / sqrt(energy + fNR.epsilon)) *
         (1.0 - 1.0 / (1.0 + pow(energy / fNR.zeta, fNR.eta)));
}
double LArNEST::GetNRPhotonYields(double energy, double efield) {
  return fNR.alpha * pow(energy, fNR.beta) -
         (1.0 / (fNR.gamma * pow(efield, fNR.delta))) *
             (1.0 / sqrt(energy + fNR.epsilon));
}
double LArNEST::GetNRPhotonYieldsConserved(double energy, double efield) {
  return fNR.alpha * pow(energy, fNR.beta) -
         (1.0 / (fNR.gamma * pow(efield, fNR.delta))) *
             (1.0 / sqrt(energy + fNR.epsilon)) *
             (1.0 - 1.0 / (1.0 + pow(energy / fNR.zeta, fNR.eta)));
}
LArYieldResult LArNEST::GetNRYields(double energy, double efield,
                                    double density) {
  double NRTotalYields = GetNRTotalYields(energy);
  double NRElectronYields = GetNRElectronYields(energy, efield);
  if (NRElectronYields < 0.0) {
    NRElectronYields = 0.0;
  }
  double NRPhotonYields = GetNRPhotonYields(energy, efield);
  if (NRPhotonYields < 0.0) {
    NRPhotonYields = 0.0;
  }
  return GetRecombinationYields(NRTotalYields, NRElectronYields, NRPhotonYields,
                                energy, efield);
}
//-------------------------ER Yields-------------------------//
double LArNEST::GetERTotalYields(double energy) {
  return (energy * 1.0e3 / fWorkQuantaFunction);
}
double LArNEST::GetERElectronYieldsAlpha(double efield, double density) {
  return fER.alpha.A +
         fER.alpha.B / (fER.alpha.C +
                        pow(efield / (fER.alpha.D +
                                      fER.alpha.E * exp(density / fER.alpha.F)),
                            fER.alpha.G));
}
double LArNEST::GetERElectronYieldsBeta(double efield) {
  return fER.beta.A +
         fER.beta.B *
             pow(fER.beta.C + pow(efield / fER.beta.D, fER.beta.E), fER.beta.F);
}
double LArNEST::GetERElectronYieldsGamma(double efield) {
  return fER.gamma.A *
         (fER.gamma.B / fWorkQuantaFunction +
          fER.gamma.C * (fER.gamma.D +
                         fER.gamma.E / pow(efield / fER.gamma.F, fER.gamma.G)));
}
double LArNEST::GetERElectronYieldsDokeBirks(double efield) {
  return fER.doke_birks.A +
         fER.doke_birks.B / (fER.doke_birks.C +
                             pow(efield / fER.doke_birks.D, fER.doke_birks.E));
}
double LArNEST::GetERElectronYields(double energy, double efield,
                                    double density) {
  double alpha = GetERElectronYieldsAlpha(efield, density);
  double beta = GetERElectronYieldsBeta(efield);
  double gamma = GetERElectronYieldsGamma(efield);
  double doke_birks = GetERElectronYieldsDokeBirks(efield);
  return alpha * beta +
         (gamma - alpha * beta) /
             pow(fER.p1 + fER.p2 * pow(energy + 0.5, fER.p3), fER.p4) +
         fER.delta / (fER.p5 + doke_birks * pow(energy, fER.let));
}
LArYieldResult LArNEST::GetERYields(double energy, double efield,
                                    double density) {
  double ERTotalYields = GetERTotalYields(energy);
  double ERElectronYields = GetERElectronYields(energy, efield, density);
  if (ERElectronYields < 0.0) {
    ERElectronYields = 0.0;
  }
  double ERPhotonYields = ERTotalYields / energy - ERElectronYields;
  if (ERPhotonYields < 0.0) {
    ERPhotonYields = 0.0;
  }
  return GetRecombinationYields(ERTotalYields, ERElectronYields, ERPhotonYields,
                                energy, efield);
}
//-------------------------Alpha Yields-------------------------//
double LArNEST::GetAlphaTotalYields(double energy) { return 0; }
double LArNEST::GetAlphaElectronYields(double efield) {
  double fieldTerm =
      fAlpha.Ye.F * pow((fAlpha.Ye.G + pow(efield, fAlpha.Ye.H)), fAlpha.Ye.I);
  return fAlpha.Ye.A *
         (fAlpha.Ye.B -
          (fAlpha.Ye.B * fAlpha.Ye.C +
           (fAlpha.Ye.B / fAlpha.Ye.D) *
               (1 - (fAlpha.Ye.E *
                     log(1 + (fAlpha.Ye.B / fAlpha.Ye.D) *
                                 (fieldTerm / fAlpha.Ye.J)) /
                     ((fAlpha.Ye.B / fAlpha.Ye.D) * fieldTerm)))));
}
double LArNEST::GetAlphaPhotonYields(double efield) {
  double quench = 1.0 / (fAlpha.Yph.A * pow(efield, fAlpha.Yph.B));
  double fieldTerm =
      fAlpha.Yph.H *
      pow((fAlpha.Yph.I + pow((efield / fAlpha.Yph.J), fAlpha.Yph.K)),
          fAlpha.Yph.L);
  return quench * fAlpha.Yph.C *
         (fAlpha.Yph.D * fAlpha.Yph.E +
          (fAlpha.Yph.D / fAlpha.Yph.F) *
              (1 - (fAlpha.Yph.G *
                    log(1 + (fAlpha.Yph.D / fAlpha.Yph.F) *
                                (fieldTerm / fAlpha.Yph.M)) /
                    ((fAlpha.Yph.D / fAlpha.Yph.F) * fieldTerm))));
}
LArYieldResult LArNEST::GetAlphaYields(double energy, double efield,
                                       double density) {
  double AlphaElectronYields = GetAlphaElectronYields(efield);
  if (AlphaElectronYields < 0.0) {
    AlphaElectronYields = 0.0;
  }
  double AlphaPhotonYields = GetAlphaPhotonYields(efield);
  if (AlphaPhotonYields < 0.0) {
    AlphaPhotonYields = 0.0;
  }
  double AlphaTotalYields = (AlphaElectronYields + AlphaPhotonYields) * energy;
  return GetRecombinationYields(AlphaTotalYields, AlphaElectronYields,
                                AlphaPhotonYields, energy, efield);
}
//---------------------------Lepton LET Yields--------------------------//
LArYieldResult LArNEST::GetLeptonLETYields(double energy, double dx,
                                           double efield, double density) {
  LArYieldResult result;
  double ionization_yields = GetCanonicalIonizationYields(energy);
  if (ionization_yields < 0.0) {
    ionization_yields = 0.0;
  }
  double exciton_yields = fNexOverNion * ionization_yields;
  double LET = 0.0;
  if (dx) {
    LET = (energy / dx) * (1 / density);  // lin. energy xfer (prop. to dE/dx)
  }
  if (LET > 0 && energy > 0 && dx > 0) {
    double ratio = GetLinearEnergyTransfer(energy * 1e3) / LET;
    if (ratio < 0.7) {
      dx /= ratio;
      LET *= ratio;
    }
  }
  return GetLETRecombinationYields(ionization_yields, exciton_yields, energy,
                                   LET, efield);
}
//------------------------------LET Yields-----------------------------//
LArYieldResult LArNEST::GetLETYields(double energy, double dx, double efield,
                                     double density) {
  LArYieldResult result;
  double ionization_yields = GetCanonicalIonizationYields(energy);
  if (ionization_yields < 0.0) {
    ionization_yields = 0.0;
  }
  double exciton_yields = fNexOverNion * ionization_yields;
  double LET = GetLinearEnergyTransfer(1000 * energy);

  if (LET) {
    dx = energy / (density * LET);  // find the range based on the LET
  }
  return GetLETRecombinationYields(ionization_yields, exciton_yields, energy,
                                   LET, efield);
}
double LArNEST::GetLETRecombinationProbability(double LET, double efield) {
  if (fUseDokeBirks) {
    // set up DokeBirks coefficients
    double DokeBirksA = 0.07 * pow((efield / 1.0e3), -0.85);
    double DokeBirksC = 0.00;
    if (efield == 0.0) {
      DokeBirksA = 0.0003;
      DokeBirksC = 0.75;
    }
    // B=A/(1-C) (see paper)
    double DokeBirksB = DokeBirksA / (1 - DokeBirksC);
    double recombProb =
        (DokeBirksA * LET) / (1 + DokeBirksB * LET) + DokeBirksC;

    // check against unphysicality resulting from rounding errors
    if (recombProb < 0.0) {
      recombProb = 0.0;
    }
    if (recombProb > 1.0) {
      recombProb = 1.0;
    }
    return recombProb;
  } else {
    return 0;
  }
}
LArYieldResult LArNEST::GetLETRecombinationYields(double ionization_yields,
                                                  double exciton_yields,
                                                  double energy, double LET,
                                                  double efield) {
  LArYieldResult result;
  double recombination_probability =
      GetLETRecombinationProbability(LET, efield);
  result.TotalYield = (ionization_yields + exciton_yields) / energy;
  result.Nex = exciton_yields;
  result.Nion = ionization_yields;
  result.Nph = exciton_yields + ionization_yields * recombination_probability;
  result.Ne = ionization_yields * (1.0 - recombination_probability);
  result.ElectricField = efield;
  return result;
}
//-----------------------------dEdx Yields-----------------------------//
double LArNEST::GetCanonicalTotalYields(double energy) {
  double mean_quanta = energy / (fWorkQuantaFunction * 1e-6);
  double Fano = sqrt(fFanoER * mean_quanta);
  return RandomGen::rndm()->rand_gauss(mean_quanta, Fano, true);
}
double LArNEST::GetCanonicalIonizationYields(double energy) {
  double mean_quanta = energy / (GetEffectiveWorkIonFunction() * 1e-6);
  double Fano = sqrt(fFanoER * mean_quanta);
  return RandomGen::rndm()->rand_gauss(mean_quanta, Fano, true);
}
LArYieldResult LArNEST::GetdEdxRecombinationYields(double ionization_yields,
                                                   double exciton_yields,
                                                   double energy, double dx,
                                                   double efield) {
  LArYieldResult result;
  double recombination_probability =
      GetdEdxRecombinationProbability(energy / dx, efield);
  result.TotalYield = (ionization_yields + exciton_yields) / energy;
  result.Nex = exciton_yields;
  result.Nion = ionization_yields;
  result.Nph = exciton_yields + ionization_yields * recombination_probability;
  result.Ne = ionization_yields * (1.0 - recombination_probability);
  result.ElectricField = efield;
  return result;
}
double LArNEST::GetdEdxRecombinationProbability(double dEdx, double efield) {
  if (fUseDokeBirks) {
    // set up DokeBirks coefficients
    double DokeBirksA = 0.07 * pow((efield / 1.0e3), -0.85);
    double DokeBirksC = 0.00;
    if (efield == 0.0) {
      DokeBirksA = 0.0003;
      DokeBirksC = 0.75;
    }
    // B=A/(1-C) (see paper)
    double DokeBirksB = DokeBirksA / (1 - DokeBirksC);
    double recombProb =
        (DokeBirksA * dEdx) / (1 + DokeBirksB * dEdx) + DokeBirksC;

    // check against unphysicality resulting from rounding errors
    if (recombProb < 0.0) {
      recombProb = 0.0;
    }
    if (recombProb > 1.0) {
      recombProb = 1.0;
    }
    return recombProb;
  } else {
    return 0;
  }
}
LArYieldResult LArNEST::GetdEdxYields(double energy, double dx, double efield,
                                      double density) {
  LArYieldResult result;
  double ionization_yields = GetCanonicalIonizationYields(energy);
  if (ionization_yields < 0.0) {
    ionization_yields = 0.0;
  }
  double exciton_yields = fNexOverNion * ionization_yields;
  return GetdEdxRecombinationYields(ionization_yields, exciton_yields, energy,
                                    dx, efield);
}
//-------------------------Fluctuation Yields-------------------------//
LArYieldFluctuationResult LArNEST::GetDataDrivenFluctuations(
    const LArYieldResult &yields, double density) {
  LArYieldFluctuationResult result{0, 0, 0, 0};
  return result;
}
LArYieldFluctuationResult LArNEST::GetDefaultFluctuations(
    const LArYieldResult &yields, double density) {
  double Fano = 0.1115;
  double NexOverNion = yields.Nex / yields.Nion;
  double p_r = (yields.Nph - yields.Nex) / yields.Nion;
  double alf = 1. / (1. + NexOverNion);
  double Nq = 0;
  double Nex = 0;
  double Nion = 0;
  double Nq_mean = yields.Ne + yields.Nph;
  // If nuclear recoils
  if (yields.Lindhard == 1.) {
    Nq = floor(
        RandomGen::rndm()->rand_gauss(Nq_mean, sqrt(Fano * Nq_mean), true) +
        0.5);
    if (Nq < 0.) {
      Nq = 0.;
    }
    Nion = double(RandomGen::rndm()->binom_draw(Nq, alf));
    if (Nion > Nq) {
      Nion = Nq;
    }
    Nex = Nq - Nion;
    // otherwise
  } else {
    Nion = floor(RandomGen::rndm()->rand_gauss(
                     Nq_mean * alf, sqrt(Fano * Nq_mean * alf), true) +
                 0.5);
    Nex = floor(RandomGen::rndm()->rand_gauss(
                    Nq_mean * alf * NexOverNion,
                    sqrt(Fano * Nq_mean * alf * NexOverNion), true) +
                0.5);
    if (Nex < 0.) {
      Nex = 0.;
    }
    if (Nion < 0.) {
      Nion = 0.;
    }
    Nq = Nion + Nex;
  }
  LArYieldFluctuationResult result{Nex + p_r * Nion, (1. - p_r) * Nion, Nex,
                                   Nion};
  return result;
}
//-------------------------Photon Times-------------------------//
double LArNEST::GetPhotonTime(LArInteraction species, bool exciton,
                              double energy) {
  // arXiv:1802.06162. NR may need tauR ~0.5-1ns instead of 0
  double SingTripRatio = 0.;
  // if ( efield > 60. ) efield = 60. // makes Xed work. 200 for LUX Run04
  // instead. A mystery! Why no field dep?

  double LET = GetLinearEnergyTransfer(energy);
  // copied from 2013 NEST version for LAr on LBNE
  double tau1 = RandomGen::rndm()->rand_gauss(
      6.5, 0.8, true);  // error from weighted average
  double tau3 = RandomGen::rndm()->rand_gauss(1300, 50, true);  // ibid.
  double tauR = RandomGen::rndm()->rand_gauss(0.8, 0.2, true);  // Kubota 1979
  if (energy < fWorkQuantaFunction * 0.001) {
    // from old G4S2Light
    tau1 = 6.;
    tau3 = 1600.;
    SingTripRatio = 0.1;
    if (!exciton) {
      tauR = 28e3;
    } else {
      tauR = 0.;  // 28 microseconds comes from Henrique:
                  // https://doi.org/10.1016/j.astropartphys.2018.04.006
    }
  } else {
    if (species == LArInteraction::NR) {
      SingTripRatio = 0.22218 * pow(energy, 0.48211);
    } else if (species == LArInteraction::Alpha) {  // really only alphas here
      SingTripRatio = (-0.065492 + 1.9996 * exp(-energy / 1e3)) /
                          (1. + 0.082154 / pow(energy / 1e3, 2.)) +
                      2.1811;  // uses energy in MeV not keV
    } else {
      if (LET < 3.) {
        if (!exciton) {
          SingTripRatio = RandomGen::rndm()->rand_gauss(0.5, 0.2, true);
        } else if (exciton) {
          SingTripRatio = RandomGen::rndm()->rand_gauss(0.36, 0.06, true);
        }
      } else {
        SingTripRatio = 0.2701 + 0.003379 * LET - 4.7338e-5 * pow(LET, 2.) +
                        8.1449e-6 * pow(LET, 3.);
      }
    }  // lastly is ER for LAr
  }
  // in case varied with Gaussian earlier
  // the recombination time is non-exponential, but approximates
  // to exp at long timescales (see Kubota '79)
  double time_ns = tauR * (1.0 / RandomGen::rndm()->rand_uniform() - 1.);

  if (RandomGen::rndm()->rand_uniform() <
      SingTripRatio / (1. + SingTripRatio)) {
    time_ns -= tau1 * log(RandomGen::rndm()->rand_uniform());
  } else {
    time_ns -= tau3 * log(RandomGen::rndm()->rand_uniform());
  }
  return time_ns;
}
double LArNEST::GetPhotonEnergy(bool state) {
  // liquid Argon
  // TODO: What is this function about?
  if (state) {
    return RandomGen::rndm()->rand_zero_trunc_gauss(9.7, 0.2);
  } else {
    return RandomGen::rndm()->rand_zero_trunc_gauss(9.69, 0.22);
  }
}
//-------------------------Drift Velocity-------------------------//
double LArNEST::GetDriftVelocity_Liquid(double Kelvin, double efield) {
  // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
  if (Kelvin < 84. || Kelvin > 140.) {
    cerr << "\nWARNING: TEMPERATURE OUT OF RANGE (84-140 K) for vD\n";
    cerr << "Using value at closest temp for a drift speed estimate\n";
    if (Kelvin < 84.) {
      Kelvin = 84.;
    } else {
      Kelvin = 140.;
    }
  }
  for (size_t i = 0; i < fDriftParameters.A.size() - 1; i++) {
    if (Kelvin >= fDriftParameters.TempLow[i] and
        Kelvin < fDriftParameters.TempHigh[i + 1]) {
      double speed =
          exp(fDriftParameters.A[i] + fDriftParameters.B[i] / (efield / 1000) +
              fDriftParameters.C[i] * log(efield / 1000));
      if (speed > 0.0) {
        return speed;
      } else {
        return 0.0;
      }
    }
  }
  return 0.0;
}
//-------------------------Utilities-------------------------//
double LArNEST::GetDensity(double Kelvin, double bara, bool &inGas,
                           uint64_t evtNum, double molarMass) {
  // currently only for fixed pressure
  // (the saturated vapor pressure); will add
  // pressure dependence later
  double VaporP_bar =
      exp(45.973940 - (1464.718291 / Kelvin) - (6.539938 * log(Kelvin)));
  if (bara < VaporP_bar || inGas) {
    double p_Pa = bara * 1e5;
    double density =
        1. /
        (pow(fRIdealGas * Kelvin, 3.) /
             (p_Pa * pow(fRIdealGas * Kelvin, 2.) + fRealGasA * p_Pa * p_Pa) +
         fRealGasB) *
        molarMass * 1e-6;  // Van der Waals equation, mol/m^3
    if (bara < VaporP_bar && evtNum == 0)
      // TODO: get rid of these output statements
      std::cerr << "\nWARNING: ARGON GAS PHASE. IS THAT WHAT YOU WANTED?\n";
    inGas = true;
    return density;
  } else {
    inGas = false;
    return 1.4;
  }
}
double LArNEST::GetLinearEnergyTransfer(double energy, bool CSDA) {
  double LET;
  if (!CSDA) {  // total stopping power directly from ESTAR (radiative +
                // collision)
    LET = 1.8106 - 0.45086 * log10(energy) - 0.33151 * pow(log10(energy), 2.) +
          0.25916 * pow(log10(energy), 3.) - 0.2051 * pow(log10(energy), 4.) +
          0.15279 * pow(log10(energy), 5.) - 0.084659 * pow(log10(energy), 6.) +
          0.030441 * pow(log10(energy), 7.) -
          0.0058953 * pow(log10(energy), 8.) +
          0.00045633 * pow(log10(energy), 9.);
    LET = pow(10., LET);
    if (std::isnan(LET) || LET <= 0.) {
      LET = 1e2;
    }
  } else {
    // the "continuous slowing down approximation" (CSDA)
    // replace with Justin and Prof. Mooney's work
    if (energy >= 1.) {
      LET = 116.70 - 162.97 * log10(energy) + 99.361 * pow(log10(energy), 2) -
            33.405 * pow(log10(energy), 3) + 6.5069 * pow(log10(energy), 4) -
            0.69334 * pow(log10(energy), 5) + 0.031563 * pow(log10(energy), 6);
    } else if (energy > 0. && energy < 1.) {
      LET = 100.;
    } else {
      LET = 0.;
    }
  }
  return LET;
}
std::vector<double> LArNEST::CalculateG2(int verbosity) {
  return std::vector<double>();
}
//------------------------------------Legacy
// LArNEST------------------------------------//
LArYieldResult LArNEST::LegacyGetYields(double energy, double efield,
                                        double yieldFactor,
                                        double excitationRatio, double epsilon,
                                        double recombProb) {
  // determine ultimate number of quanta from current E-deposition (ph+e-)
  // total mean number of exc/ions the total number of either quanta produced
  // is equal to product of the work function, the energy deposited,
  // and yield reduction, for NR
  double MeanNq = legacy_scint_yield * energy;
  double sigma = sqrt(legacy_resolution_scale * MeanNq);  // Fano
  double leftvar =
      RandomGen::rndm()->rand_gauss(yieldFactor, 0.25 * yieldFactor, true);
  if (leftvar > 1.0) {
    leftvar = 1.0;
  }
  int Nq = int(floor(RandomGen::rndm()->rand_gauss(MeanNq, sigma, true) + 0.5));
  if (yieldFactor < 1) {
    Nq = RandomGen::rndm()->binom_draw(Nq, leftvar);
  }

  // if Edep below work function, can't make any quanta, and if Nq
  // less than zero because Gaussian fluctuated low, update to zero
  if (energy < 1 / legacy_scint_yield || Nq < 0) {
    Nq = 0;
  }

  // next section binomially assigns quanta to excitons and ions
  double Nex = RandomGen::rndm()->binom_draw(
      Nq, excitationRatio / (1 + excitationRatio));
  double Nion = Nq - Nex;

  // use binomial distribution to assign photons, electrons, where photons
  // are excitons plus recombined ionization electrons, while final
  // collected electrons are the "escape" (non-recombined) electrons
  double Nph = Nex + RandomGen::rndm()->binom_draw(Nion, recombProb);
  double Ne = Nq - Nph;

  // create the quanta results
  LArYieldResult result{(Ne + Nph) / energy,
                        Ne / energy,
                        Nph / energy,
                        Nph,
                        Ne,
                        Nex,
                        Nion,
                        0.0,
                        efield};
  return result;
}
LArYieldResult LArNEST::LegacyCalculation(int pdgcode, double energy,
                                          double efield, double density,
                                          double track_length) {
  // determine various parameter values, such as the
  // yieldfactor, excitationratio, dokebirks parameters,
  // epislon and the recombination probability
  // default quenching factor, for electronic recoils
  double yieldFactor = 1.0;
  // ratio for light particle in LAr, such as e-, mu-, Aprile et. al book
  double excitationRatio = 0.21;

  // nuclear recoil quenching "L" factor: total yield is
  // reduced for nuclear recoil as per Lindhard theory
  double epsilon = 11.5 * (energy * pow(LAr_Z, (-7. / 3.)));
  if (abs(pdgcode) == 2112)  // nuclear recoil
  {
    yieldFactor = 0.23 * (1 + exp(-5 * epsilon));  // liquid argon L_eff
    excitationRatio = 0.69337 + 0.3065 * exp(-0.008806 * pow(efield, 0.76313));
  }
  // get the recombination probability
  double recombProb = LegacyGetRecombinationProbability(energy, efield, density,
                                                        pdgcode, track_length);
  return LegacyGetYields(energy, efield, yieldFactor, excitationRatio, epsilon,
                         recombProb);
}
double LArNEST::LegacyGetRecombinationProbability(double energy, double efield,
                                                  double density, int pdgcode,
                                                  double track_length) {
  // this section calculates recombination following the modified
  // Birks' Law of Doke, deposition by deposition, may be overridden
  // later in code if a low enough energy necessitates switching to the
  // Thomas-Imel box model for recombination instead (determined by site)
  double dE = energy;
  double dx = 0.0;
  double LET = 0.0;
  double recombProb;
  if (abs(pdgcode) != 11 && abs(pdgcode) != 13) {
    // e-: 11, e+: -11, mu-: 13, mu+: -13
    // in other words, if it's a gamma,ion,proton,alpha,pion,et al. do not
    // use the step length provided by Geant4 because it's not relevant,
    // instead calculate an estimated LET and range of the electrons that
    // would have been produced if Geant4 could track them
    LET = LegacyGetLinearEnergyTransfer(dE);
    if (LET) {
      // find the range based on the LET
      dx = (dE / 1e3) / (density * LET);
    }
    if (abs(pdgcode) == 2112)  // nuclear recoils
    {
      dx = 0;
    }
  } else  // normal case of an e-/+ energy deposition recorded by Geant
  {
    dx = track_length / 10.;
    if (dx) {
      LET = ((dE / 1e3) / dx) *
            (1 / density);  // lin. energy xfer (prop. to dE/dx)
    }
    if (LET > 0 && dE > 0 && dx > 0) {
      double ratio = LegacyGetLinearEnergyTransfer(dE) / LET;
      if (ratio < 0.7 && pdgcode == 11) {
        dx /= ratio;
        LET *= ratio;
      }
    }
  }
  // set up DokeBirks coefficients
  double DokeBirksA = 0.07 * pow((efield / 1.0e3), -0.85);
  double DokeBirksC = 0.00;
  if (efield == 0.0) {
    DokeBirksA = 0.0003;
    DokeBirksC = 0.75;
  }
  // B=A/(1-C) (see paper)
  double DokeBirksB = DokeBirksA / (1 - DokeBirksC);
  recombProb = ((DokeBirksA * LET) / (1 + DokeBirksB * LET) + DokeBirksC) *
               (density / legacy_density_LAr);

  // check against unphysicality resulting from rounding errors
  if (recombProb < 0.0) {
    recombProb = 0.0;
  }
  if (recombProb > 1.0) {
    recombProb = 1.0;
  }
  return recombProb;
}
double LArNEST::LegacyGetLinearEnergyTransfer(double E) {
  double LET;
  if (E >= 1) {
    LET = 116.70 - 162.97 * log10(E) + 99.361 * pow(log10(E), 2) -
          33.405 * pow(log10(E), 3) + 6.5069 * pow(log10(E), 4) -
          0.69334 * pow(log10(E), 5) + 0.031563 * pow(log10(E), 6);
  } else if (E > 0 && E < 1) {
    LET = 100;
  } else {
    LET = 0;
  }
  return LET;
}
}  // namespace NEST