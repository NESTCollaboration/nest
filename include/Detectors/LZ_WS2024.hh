//
// LZ_WS2024.hh
//
// Adapted from Quentin Riffard and Jacob Cutter by Greg Rischbieter, July 2024
//
// This file serves as a NEST input to reproduce LZ's WS2024 result
//
// This version of the file is very similar to https://www.hepdata.net/record/155182
// as linked from https://arxiv.org/abs/2410.17036
// Additional methods to return the specified yield parameters below have been added
//
// Please reach out to
// Matthew Szydagis (mszydagis@albany.edu) with any questions.
//
//
//////////////////////////////  IMPORTANT !! //////////////////////////////////////
//  
//  It is recommended to use NEST v2.4.0 with this detector file. 
//
//  To match our calibration data (3H and D-D), the NEST models parameters were adjusted.
//  To make these changes in your copy of NEST, make sure NEST.cpp is calling the
//  edited NRYieldsParam, ERYieldsParam, and NRERWidthsParams vectors
//  when using the GetYields(...) and GetQuanta(...) functions:
//
//    NRYieldsParam = {10.19, 1.11, 0.0498, -0.0533, 12.46, 0.2942, 1.899, 0.3197, 2.066, 0.509, 0.996, 0.999};
//    ERYieldsParam = {12.4886, 85.0, 0.6050, 2.14687, 25.721, -1.0, 59.651, 3.6869, 0.2872, 0.1121};
//    NRERWidthsParam = {0.404, 0.393, 0.0383, 0.497, 0.1906, 2.220, 0.3, 0.04311, 0.46894, 0.15505, -0.26564, 0., 0.};
//
//  The easiest way to implement this in NEST is to change the "default" vectors in 
//  lines 120-124 of nest/include/NEST/NEST.hh.
//  Else, If using execNEST to simulate results, these changes can be made manually 
//  to execNEST.cpp at lines 186-200 in NESTv2.4.0
//  Or the updated "default" vectors from NEST.hh can be added to the execNEST function
//  around line 30 of execNEST (i.e. NRYieldsParam = default_NRYieldsParam; etc.)
//
//  To add this file into execNEST, make sure it is copied into the directory 
//                        nest/include/Detectors/
//  and add the line ' #include "LZ_WS2024.hh" ' , and make sure the "auto detector"
//  variable is the "LZ_Detector_2024()" object, as opposed to any earlier NEST default.
////////////////////////////////////////////////////////////////////////////////////
//
//  To use this detector with nestpy, you'll need to (1) compile a local install of 
//  nestpy from github.com:NESTCollaboration/nestpy, and (2) then add this file to 
//  nestpy/lib/nest/include/Detectors
//  Then (3) add "#include LZ_WS2024.hh" in nestpy/src/nestpy/bindings.cpp, and 
//  finally (4) create the LZ-specific bindings, similar to the XENON10 and LUX 
//  bindings examples.
//
//  With nestpy, you can pass the NRYieldsParam and ERYieldsParam vectors as lists 
//  to NESTCalc.GetYields(...) with the kwargs "nuisance_parameters = NRYieldsParams"
//  and "ERYieldsParams = ERYieldsParams".
//  For the width parameters, you can add "free_parameters = NRERWidthsParam" in 
//  NESTcalc.GetQuanta(...) 
//
////////////////////////////////////////////////////////////////////////////////////

#ifndef LZ_Detector_2024_hh
#define LZ_Detector_2024_hh 1

#include "VDetector.hh"

using namespace std;

class LZ_Detector_2024 : public VDetector {
 public:
  LZ_Detector_2024() {
    // Call the initialization of all the parameters
    Initialization();
  };
  ~LZ_Detector_2024() override = default;

  // Do here the initialization of all the parameters that are not varying as a
  // function of time
  void Initialization() override {
    name = "LZ SR3 or WS2024";
    // Primary Scintillation (S1) parameters
    g1 = 0.1122;// +/- 0.002  // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
    sPEres = 0.338;   // single phe resolution (Gaussian assumed)
    sPEthr = 0.10;   // POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 1.0;   // actual efficiency, can be used in lieu of POD threshold
    noiseBaseline[0] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noiseBaseline[1] = 0.0;  // baseline noise mean and width in PE (Gaussian)
    noiseBaseline[2] = 0.; 
    noiseBaseline[3] = 0.;
    P_dphe = 0.214;  // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

    coinWind = 150;  // S1 coincidence window in ns
    coinLevel = 3;   // how many PMTs have to fire for an S1 to count
    numPMTs = 481; //Taking into account turned-off PMTs    // For coincidence calculation

    OldW13eV = true;
    noiseLinear[0] = 0.; //(.083 to increase ER band widths at higher Es: arXiv:2312.02030)
    noiseLinear[1] = 0.;

    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.076404;// +/- 0.002 // phd per S2 photon in gas, used to get SE size
    s2Fano = 4.0;  // Fano-like fudge factor for SE width
    s2_thr = 783.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
    E_gas = 8.301178;    // effective field in kV/cm between liquid/gas border and anode
			 // gives extraction efficiency of ~75%
    eLife_us = 9000.; //the drift electron mean lifetime in micro-seconds

    // Thermodynamic Properties
    inGas = false;
    T_Kelvin = 175.9;  // for liquid drift speed calculation
    p_bar = 1.86;     // gas pressure in units of bars, it controls S2 size
    // if you are getting warnings about being in gas, lower T and/or raise p

    // Data Analysis Parameters and Geometry
    dtCntr = 525.; // center of detector for S1 corrections, in usec.
    dt_min = 71.; //minimum. Top of detector fiducial volume
    dt_max = 1029.;  // maximum. Bottom of detector fiducial volume

    radius = 688.;  // millimeters (fiducial radius)
    radmax = 728.;  // actual physical geo. limit

    TopDrift = 1461.;  // mm not cm or us (but, this *is* where dt=0)
    // a z-axis value of 0 means the bottom of the detector (cathode OR bottom
    // PMTs)
    // In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not
    // anode!
    anode = 1469.;  // the level of the anode grid-wire plane in mm
    // In a gas TPC, this is not TopDrift (top of drift region), but a few mm
    // above it
    gate = 1456.;  // mm. This is where the E-field changes (higher)
    // in gas detectors, the gate is still the gate, but it's where S2 starts
    cathode = 0.;  // mm. Defines point below which events are gamma-X

    // 2-D (X & Y) Position Reconstruction
    // Set these to zero to implement "perfect" position corrections
    // Note: LZ used spatial maps to implement S1 and S2 corrections
    PosResFlat = 3.000;  // constant syst unc in position recon res, mm
    PosResBase = 120.6;  // amplitude in mm, divided by sqrt of S2 area
  }
  
  // S1 PDE custom fit for function of z
  // s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  double FitS1(double xPos_mm, double yPos_mm, double zPos_mm, LCE map) override {
    //Based on S1 Map from C. Nedlik using MDC3 83m-Kr data
    double zPos_cm = zPos_mm/10;
    double Rsq_cm = xPos_mm*xPos_mm/10./10. + yPos_mm*yPos_mm/10/10.;
    double Rsq_max = (get_radmax()/10.)*(get_radmax()/10.);
    if (Rsq_cm > Rsq_max)
      Rsq_cm = Rsq_max;
    double p0 = 1.3198 + -5.2742e-05*Rsq_cm + 2.2966e-08*pow(Rsq_cm, 2.) + -6.8098e-12*pow(Rsq_cm, 3.) + 6.4871e-16*pow(Rsq_cm, 4.);
    double p1 = -6.2930e-03 + 2.9826e-06*Rsq_cm + -1.8987e-09*pow(Rsq_cm, 2.) + 5.2800e-13*pow(Rsq_cm, 3.) + -4.4126e-17*pow(Rsq_cm, 4.);
    double p2 = 3.0634e-05 + -8.1693e-08*Rsq_cm + 6.1091e-11*pow(Rsq_cm, 2.) + -1.7352e-14*pow(Rsq_cm, 3.) + 1.4823e-18*pow(Rsq_cm, 4.);
    double p3 = -8.0736e-08 + 8.5948e-10*Rsq_cm + -6.8764e-13*pow(Rsq_cm, 2.) + 1.9889e-16*pow(Rsq_cm, 3.) + -1.7312e-20*pow(Rsq_cm, 4.);
    double p4 = 1.3926e-10 + -2.9924e-12*Rsq_cm + 2.4693e-15*pow(Rsq_cm, 2.) + -7.2454e-19*pow(Rsq_cm, 3.) + 6.3858e-23*pow(Rsq_cm, 4.);
    return p0 + p1*zPos_cm + p2*pow(zPos_cm, 2.) + p3*pow(zPos_cm, 3.) + p4*pow(zPos_cm, 4.);
  }

  // Drift electric field as function of Z in mm
  double FitEF(double xPos_mm, double yPos_mm,
                       double zPos_mm) override {  // in V/cm
		return 96.5; // V/cm
  }

  // S2 PDE custom fit for function of r
  // s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  double FitS2(double xPos_mm, double yPos_mm, LCE map) override {
    double Rsq_cm = xPos_mm*xPos_mm/10./10. + yPos_mm*yPos_mm/10./10.;
    double Rsq_max = (get_radmax()/10.)*(get_radmax()/10.);
    if (Rsq_cm > Rsq_max)
      Rsq_cm = Rsq_max;
    double p0 =  1.01728; //   +/-   0.00951576
    double p1 = -0.000148818;//   +/-   3.73192e-05
    double p2 =  1.11934e-07;//   +/-   4.70456e-08
    double p3 = -4.27149e-11;//  +/-   2.50216e-11
    double p4 =  6.89494e-15;//   +/-   5.85633e-15
    double p5 = -4.29735e-19;//  +/-   4.97409e-19
    //polynomial fit by A. Stevens, using April MDC3 Xe-131m data
    return p0 + Rsq_cm*p1 + Rsq_cm*Rsq_cm*p2 + pow(Rsq_cm, 3.)*p3 + pow(Rsq_cm, 4.)*p4 + pow(Rsq_cm, 5.)*p5; //unitless: `S2(x,y)/S2(0,0)
  }

  vector<double> FitTBA(double xPos_mm, double yPos_mm,
                                double zPos_mm) override {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.60;  // S1 bottom-to-total ratio
    BotTotRat[1] = 0.356; // S2 bottom-to-total ratio, typically only used for
                          // position recon (1-this)

    return BotTotRat;
  }
  
  double OptTrans(double xPos_mm, double yPos_mm, double zPos_mm) override {
    double phoTravT, z_m = zPos_mm/1000.;
    
    double A = 0.113 - 0.19 * z_m + 0.185 * pow(z_m, 2.);
    if (A < 0.) A = 0.;  // cannot have negative probability
    double B_a = 0.983 + 0.137 * z_m + 0.033 * pow(z_m, 2.);
    double B_b = 1. - B_a;
    double tau_a = 47.304 + 3.26 * z_m - 5.077 * pow(z_m, 2.) + 1.468 * pow(z_m, 3.);  // all times in nanoseconds
    double tau_b = 14.138 - 30.180 * z_m + 52.954 * pow(z_m, 2.) - 19.330 * pow(z_m, 3.);
    if (tau_b < 0.) tau_b = 0.;  // cannot have negative time
    
    double Zmax = max(153.9-z_m*100., 14.8+z_m*100.);
    double Tdirect = sqrt(pow(82.73, 2.)+pow(Zmax, 2.))/17.75; //in ns
    
    vector<double> t_array;
    for (int ii = 0; ii < 250; ++ii) t_array.push_back(double(ii));
    vector<double> pdf_array;
    for (int ip = 0; ip < 250; ++ip) {
      double pdf = B_a/tau_a*exp(-t_array[ip]/tau_a)+B_b/tau_b*exp(-t_array[ip]/tau_b);
      pdf_array.push_back(pdf);
    }
    double total = accumulate(pdf_array.begin(), pdf_array.end(), 0.0);
    vector<double> pdfNorm_array;
    for (int pp = 0; pp < 250; ++pp) {
      double newpdf = pdf_array[pp]/total;
      pdfNorm_array.push_back(newpdf);
    }
    vector<double> cdf_array;
    cdf_array.push_back(0.);
    for (int c = 1; c < 250; ++c) {
      double tem = 0.;
      for (int cc = 0; cc < c; ++cc) tem += pdfNorm_array[cc];
      cdf_array.push_back(tem);
    }
    
    if ( RandomGen::rndm()->rand_uniform() < A ) {
      double Zmin = min(153.9-z_m*100., 14.8+z_m*100.);
      double TdirectMin = Zmin/17.75; // in ns
      phoTravT = TdirectMin + (Tdirect-TdirectMin)*RandomGen::rndm()->rand_uniform();  // direct travel time to PMTs (low)
    }
    else {
      double y = RandomGen::rndm()->rand_uniform();
      vector<double> diff_array;
      for (int dd = 0; dd < 250; ++dd) {
	double diff = fabs(cdf_array[dd]-y);
	diff_array.push_back(diff);
      }
      auto it = min_element(begin(diff_array), end(diff_array));
      int ind = distance(begin(diff_array), it);
      phoTravT = t_array[ind];
      phoTravT += Tdirect; //ns
    }
    
    double k_add = z_m*100./2.;  // if you want more realistic photon transit times (cf. BACCARAT) keep this smearing
    phoTravT += k_add * RandomGen::rndm()->rand_uniform();
    //double sig = RandomGen::rndm()->rand_gauss(3.84, .09, true); // includes stat unc but not syst
    //phoTravT += RandomGen::rndm()->rand_gauss(0.00, sig,false);  // the overall width added to photon time spectra by the effects in the electronics and the data reduction pipeline
    
    if (phoTravT > DBL_MAX) phoTravT = tau_a;
    if (phoTravT < -DBL_MAX) phoTravT = 0.00;
    
    return phoTravT;  // this function follows arXiv:2603.26877
    
  }
  
  // The following function was not used in LZ's WS2024. It has been copied from the public NEST file for LUX, just so that NEST has it available, in order to prevent an error.
  vector<double> SinglePEWaveForm(double area, double t0) override {
    vector<double> PEperBin;

    double threshold = PULSEHEIGHT;  // photo-electrons
    double sigma = PULSE_WIDTH;      // ns
    area *= 10. * (1. + threshold);
    double amplitude = area / (sigma * sqrt(2. * M_PI)),
           signal;  // assumes perfect Gaussian

    double tStep1 = SAMPLE_SIZE / 1e2;  // ns, make sure much smaller than
                                        // sample size; used to generate MC-true
                                        // pulses essentially
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

    std::vector<double> get_nr_yield_params(){
        return m_NRYieldsParam;
    }
    std::vector<double> get_er_yield_params(){
        return m_ERYieldsParam;
    }
    std::vector<double> get_nr_er_width_params(){
        return m_NRERWidthsParam;
    }

  private:
    std::vector<double> m_NRYieldsParam = {10.19, 1.11, 0.0498, -0.0533, 12.46, 0.2942, 1.899, 0.3197, 2.066, 0.509, 0.996, 0.999};
    std::vector<double> m_ERYieldsParam = {12.4886, 85.0, 0.6050, 2.14687, 25.721, -1.0, 59.651, 3.6869, 0.2872, 0.1121};
  std::vector<double> m_NRERWidthsParam = {0.404, 0.393, 0.0383, 0.497, 0.1906, 2.220, 0.3, 0.04311, 0.46894, 0.15505, -0.26564, 0., 0.};
};

#endif
