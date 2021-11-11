
#ifndef DetectorExample_LUX_RUN03_hh
#define DetectorExample_LUX_RUN03_hh 1

#include "VDetector.hh"
using namespace std;

// NOTES: best g1 for DD 0.115, but for tritium 0.119; g1_gas = 0.1 for both,
// s2Fano & s2_thr equal E_gas 6.17, 6.35; e- lifet 650, 750 us; fid vol 80-130,
// 38-305 us; gasGap and noiseL0,1 same (both ~0) DISCLAIMER: Slight differences
// from official published values due to private LUX algorithms

class DetectorExample_LUX_RUN03 : public VDetector {
 public:
  DetectorExample_LUX_RUN03() {
    // cerr << "*** Detector definition message ***" << endl;
    // cerr << "You are currently using the LUX Run03 template detector." <<
    // endl << endl;

    // Call the initialization of all the parameters
    Initialization();
  };
  ~DetectorExample_LUX_RUN03() override = default;

  // Do here the initialization of all the parameters that are not varying as a
  // function of time
  void Initialization() override {
    // Primary Scintillation (S1) parameters
    g1 = 0.1170;  // 0.117+/-0.003 WS,0.115+/-0.005 D-D,0.115+/-0.005
                  // CH3T,0.119+/-0.001 LUXSim. UNITS: phd per photon (NOT
                  // photoelectrons!!)
    sPEres = 0.37;  // arXiv:1910.04211. UNITS: phe a.k.a. PE or photoelectrons
    sPEthr = (0.3 * 1.173) / 0.915;  // arXiv:1910.04211. UNITS: phe
    sPEeff = 1.00;                   // arXiv:1910.04211. UNITS: fractional
    noiseBaseline[0] = 0.00;         // arXiv:1910.04211 says -0.01. UNITS: phe
    noiseBaseline[1] = 0.08;         // arXiv:1910.04211. UNITS: phe
    noiseBaseline[2] = 0.;           // UNITS: e-'s
    noiseBaseline[3] = 0.;           // UNITS: e-'s
    P_dphe = 0.173;                  // arXiv:1910.04211. UNITS: fractional

    coinWind = 100;  // 1310.8214. UNITS: ns
    coinLevel =
        2;  // 1512.03506. UNITS: number of PMTs for coincidence requirement
    numPMTs = 119;  // 122 minus 3 off. UNITS: number of PMTs

    OldW13eV = true;  // default true, which means use "classic" W instead of
                      // Baudis / EXO's
    noiseLinear[0] = 0.0e-2;  // 1910.04211 p.12, to match 1610.02076 Fig. 8.
                              // UNITS: fraction NOT %!
    noiseLinear[1] = 0.0e-2;  // 1910.04211 p.12, to match 1610.02076 Fig. 8.
                              // UNITS: fraction NOT %!

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.1;  // 0.1 in 1910.04211. UNITS: phd per e-
    s2Fano = 3.6;  // 3.7 in 1910.04211; this matches 1608.05381 better.
                   // Dimensionless
    s2_thr = (150. * 1.173) / 0.915;  // 65-194 pe in 1608.05381. UNITS: phe
    E_gas = 6.25;                     // 6.55 in 1910.04211. UNITS: kV/cm
    eLife_us =
        800.;  // p.44 of James Verbus PhD thesis Brown. UNIT: microseconds (us)

    // Thermodynamic Properties
    // inGas = false; //duh
    T_Kelvin = 173.;  // 1910.04211. UNITS: Temperature in Kelvin
    p_bar = 1.57;     // 1910.04211. UNITS: pressure in bar

    // Data Analysis Parameters and Geometry
    dtCntr =
        160.;  // p.61 Dobi thesis UMD, 159 in 1708.02566. UNITS: microseconds
    dt_min = 38.;   // 1608.05381. UNITS: microseconds
    dt_max = 305.;  // 1608.05381. UNITS: microseconds

    radius = 200.;  // 1512.03506. UNITS: mm
    radmax = 235.;  // 1910.04211. UNITS: mm

    TopDrift = 544.95;  // 544.95 in 1910.04211. UNITS: mm
    anode = 549.2;      // 1910.04211 and 549 in 1708.02566. UNITS: mm
    gate = 539.2;       // 1910.04211 and 539 in 1708.02566. UNITS: mm
    cathode = 55.90;    // 55.9-56 in 1910.04211,1708.02566. UNITS: mm

    // 2-D (X & Y) Position Reconstruction
    PosResExp = 0.015;     // arXiv:1710.02752 indirectly. UNITS: mm^-1
    PosResBase = 70.8364;  // 1710.02752 indirectly. UNITS: mm
  }

  // S1 PDE custom fit for function of xyz
  // 1712.05696 indirectly, 1708.02566 Figure 10 color map
  double FitS1(double xPos_mm, double yPos_mm, double zPos_mm,
               LCE map) override {
    double radius = sqrt(pow(xPos_mm, 2.) + pow(yPos_mm, 2.));
    double amplitude = 307.9 - 0.3071 * zPos_mm + 0.0002257 * pow(zPos_mm, 2.);
    double shape = 1.1525e-7 * sqrt(std::abs(zPos_mm - 318.84));
    double finalCorr = -shape * pow(radius, 3.) + amplitude;
    finalCorr /= 307.9;
    if ((finalCorr < 0.5 || finalCorr > 1.5 || std::isnan(finalCorr)) &&
        radius < radmax) {
      cerr << "ERR: S1 corrections exceed a 50% difference. Are you sure you "
              "didn't forget to change LUX numbers for your own detector??"
           << endl;
      return 1.;
    } else
      return finalCorr;
  }

  // Drift electric field as function of Z in mm
  // 1709.00095, 1904.08979, 1708.02566 Fig. 13
  double FitEF(double xPos_mm, double yPos_mm,
               double zPos_mm) override {  // in V/cm

    double finalEF =
        158.92  // NOTE: DO NOT JUST RETURN A CONSTANT, THAT IS A SILLY USE of
                // FitEF
        - 0.2209000 * pow(zPos_mm, 1.) + 0.0024485 * pow(zPos_mm, 2.) -
        8.7098e-6 * pow(zPos_mm, 3.) + 1.5049e-8 * pow(zPos_mm, 4.) -
        1.0110e-11 * pow(zPos_mm, 5.);
    if (finalEF <= FIELD_MIN || finalEF > 1e5 || std::isnan(finalEF)) {
      cerr << "ERR: Very weird drift electric field value!! Are you sure you "
              "didn't forget to change LUX numbers for your own detector??"
           << endl;
      return FIELD_MIN;
    } else
      return finalEF;
  }

  // S2 PDE custom fit for function of r
  // 1712.05696 & 1710.02752 indirectly. Fig. 13 1708.02566
  double FitS2(double xPos_mm, double yPos_mm, LCE map) override {
    double radius = sqrt(pow(xPos_mm, 2.) + pow(yPos_mm, 2.));

    double finalCorr =  // unitless, 1.000 at detector center
        9156.3 + 6.22750 * pow(radius, 1.) + 0.38126 * pow(radius, 2.) -
        0.017144 * pow(radius, 3.) + 0.0002474 * pow(radius, 4.) -
        1.6953e-6 * pow(radius, 5.) + 5.6513e-9 * pow(radius, 6.) -
        7.3989e-12 * pow(radius, 7.);
    finalCorr /= 9156.3;
    if ((finalCorr < 0.5 || finalCorr > 1.5 || std::isnan(finalCorr)) &&
        radius < radmax) {
      cerr << "ERR: S2 corrections exceed a 50% difference. Are you sure you "
              "didn't forget to change LUX numbers for your own detector??"
           << endl;
      return 1.;
    } else
      return finalCorr;
  }

  vector<double> FitTBA(double xPos_mm, double yPos_mm,
                        double zPos_mm) override {
    vector<double> BotTotRat(2);

    double radSq = (pow(xPos_mm, 2.) + pow(yPos_mm, 2.)) / 1e2;
    double TBAzS1 = -0.853 + 0.00925 * (zPos_mm / 10.);
    double TBArS2 = 0.126 + 0.000545 * radSq - 1.90e-6 * radSq * radSq +
                    1.20e-9 * radSq * radSq * radSq;

    if (TBAzS1 < -1.) TBAzS1 = -1.;
    if (TBAzS1 > 1.0) TBAzS1 = 1.0;
    if (TBArS2 < -1.) TBArS2 = -1.;
    if (TBArS2 > 1.0) TBArS2 = 1.0;

    BotTotRat[0] = (1. - TBAzS1) / 2.;  // 1712.05696
    BotTotRat[1] = 0.449;  //(1.-TBArS2)/2.;  // 1712.05696 and 1710.02752
                           // position recon (1-this)

    return BotTotRat;
  }

  double OptTrans(double xPos_mm, double yPos_mm, double zPos_mm) override {
    double phoTravT, approxCenter = (TopDrift + cathode) / 2.,
                     relativeZ = zPos_mm - approxCenter;

    double A = 0.048467 - 7.6386e-6 * relativeZ +
               1.2016e-6 * pow(relativeZ, 2.) - 6.0833e-9 * pow(relativeZ, 3.);
    if (A < 0.) A = 0.;  // cannot have negative probability
    double B_a = 0.99373 + 0.0010309 * relativeZ -
                 2.5788e-6 * pow(relativeZ, 2.) -
                 1.2000e-8 * pow(relativeZ, 3.);
    double B_b = 1. - B_a;
    double tau_a = 11.15;  // all times in nanoseconds
    double tau_b = 4.5093 + 0.03437 * relativeZ -
                   0.00018406 * pow(relativeZ, 2.) -
                   1.6383e-6 * pow(relativeZ, 3.);
    if (tau_b < 0.) tau_b = 0.;  // cannot have negative time

    // A = 0.0574; B_a = 1.062; tau_a = 11.1; tau_b = 2.70; B_b = 1.0 - B_a;
    // //LUX D-D conditions

    if (RandomGen::rndm()->rand_uniform() < A)
      phoTravT = 0.;  // direct travel time to PMTs (low)
    else {            // using P0(t) =
            // A*delta(t)+(1-A)*[(B_a/tau_a)e^(-t/tau_a)+(B_b/tau_b)e^(-t/tau_b)]
            // LUX PSD paper, but should apply to all detectors w/ diff #'s
      if (RandomGen::rndm()->rand_uniform() < B_a)
        phoTravT = -tau_a * log(RandomGen::rndm()->rand_uniform());
      else
        phoTravT = -tau_b * log(RandomGen::rndm()->rand_uniform());
    }

    double sig = RandomGen::rndm()->rand_gauss(
        3.84, .09);  // includes stat unc but not syst
    phoTravT += RandomGen::rndm()->rand_gauss(
        0.00,
        sig);  // the overall width added to photon time spectra by the effects
               // in the electronics and the data reduction pipeline

    if (phoTravT > DBL_MAX) phoTravT = tau_a;
    if (phoTravT < -DBL_MAX) phoTravT = 0.000;

    return phoTravT;  // this function follows LUX (arXiv:1802.06162)
  }

  vector<double> SinglePEWaveForm(double area, double t0) override {
    vector<double> PEperBin;

    double threshold = PULSEHEIGHT;  // photo-electrons
    double sigma = PULSE_WIDTH;      // ns
    area *= 10. * (1. + threshold);
    double amplitude = area / (sigma * sqrt(2. * M_PI)),
           signal;  // assumes perfect Gaussian

    double tStep1 =
        SAMPLE_SIZE / 1e2;  // ns, make sure much smaller than sample size; used
                            // to generate MC-true pulses essentially
    double tStep2 =
        SAMPLE_SIZE;  // ns; 1 over digitization rate, 100 MHz assumed here

    double time = -5. * sigma;
    bool digitizeMe = false;
    while (true) {
      signal = amplitude * exp(-pow(time, 2.) / (2. * sigma * sigma));
      if (signal < threshold) {
        if (digitizeMe)
          break;
        else
          ;  // do nothing - goes down to advancing time block
      } else {
        if (digitizeMe)
          PEperBin.push_back(signal);
        else {
          if (RandomGen::rndm()->rand_uniform() < 2. * (tStep1 / tStep2)) {
            PEperBin.push_back(time + t0);
            PEperBin.push_back(signal);
            digitizeMe = true;
          } else {
          }
        }
      }
      if (digitizeMe)
        time += tStep2;
      else
        time += tStep1;
      if (time > 5. * sigma) break;
    }

    return PEperBin;
  }

  // Vary VDetector parameters through custom functions
  virtual void ExampleFunction() { set_g1(0.1167); }
  virtual void ExampleFunction2() { set_molarMass(132.); }
};

#endif
