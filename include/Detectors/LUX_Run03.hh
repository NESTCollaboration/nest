
#ifndef DetectorExample_LUX_RUN03_hh
#define DetectorExample_LUX_RUN03_hh 1

#include "VDetector.hh"
using namespace std;

class DetectorExample_LUX_RUN03: public VDetector {
  
public:
  
  DetectorExample_LUX_RUN03() {
    cerr << "*** Detector definition message ***" << endl;
    cerr << "You are currently using the LUX Run03 template detector." << endl << endl;
    
    // Call the initialization of all the parameters
    Initialization();
  };
  virtual ~DetectorExample_LUX_RUN03() {};
  
  // Do here the initialization of all the parameters that are not varying as a function of time
  virtual void Initialization() {
    
    // Primary Scintillation (S1) parameters
    g1 = 0.115; //phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
    sPEres = 0.37; //single phe resolution (Gaussian assumed)
    sPEthr = (0.3*1.173)/0.915; //POD threshold in phe, usually used IN PLACE of sPEeff
    sPEeff = 1.00; //actual efficiency, can be used in lieu of POD threshold
    noiseB[0] =-0.01; //baseline noise mean and width in PE (Gaussian)
    noiseB[1] = 0.08; //baseline noise mean and width in PE (Gaussian)
    noiseB[2] = 0.;
    noiseB[3] = 0.;
    P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT
    
    coinWind= 100;// S1 coincidence window in ns
    coinLevel=2;  //how many PMTs have to fire for an S1 to count
    numPMTs = 119; //For coincidence calculation
    
    extraPhot =false;
    noiseL[0]=1.4e-2;
    noiseL[1]=1.4e-2;
    
    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.09775; //phd per S2 photon in gas, used to get SE size
    s2Fano = 3.0; //Fano-like fudge factor for SE width
    s2_thr = 164.;//(150.*1.173)/0.915; //the S2 threshold in phe or PE, *not* phd. Affects NR most
    E_gas = 6.2; //field in kV/cm between liquid/gas border and anode
    eLife_us = 650.; //the drift electron mean lifetime in micro-seconds
    
    // Thermodynamic Properties
    inGas = false;
    T_Kelvin = 173.; //for liquid drift speed calculation
    p_bar = 1.57; //gas pressure in units of bars, it controls S2 size
    //if you are getting warnings about being in gas, lower T and/or raise p
    
    // Data Analysis Parameters and Geometry
    dtCntr = 160.; //center of detector for S1 corrections, in usec.
    dt_min = 80.; //minimum. Top of detector fiducial volume
    dt_max = 130.; //maximum. Bottom of detector fiducial volume
    
    radius = 200.; //millimeters
    radmax = 235.;
    
    TopDrift = 544.7; //mm not cm or us (but, this *is* where dt=0)
    //a z-axis value of 0 means the bottom of the detector (cathode OR bottom PMTs)
    //In 2-phase, TopDrift=liquid/gas border. In gas detector it's GATE, not anode!
    anode = 549.2; //the level of the anode grid-wire plane in mm
    //In a gas TPC, this is not TopDrift (top of drift region), but a few mm above it
    gate = 539.2; //mm. This is where the E-field changes (higher)
    // in gas detectors, the gate is still the gate, but it's where S2 starts
    cathode = 55.90; //mm. Defines point below which events are gamma-X
    
    // 2-D (X & Y) Position Reconstruction
    PosResExp = 0.015; // exp increase in pos recon res at hi r, 1/mm
    PosResBase = 70.8364; // baseline unc in mm, see NEST.cpp for usage
  }
  
  //S1 PDE custom fit for function of z
  //s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binom dist) e.g.
  virtual double FitS1 ( double xPos_mm, double yPos_mm, double zPos_mm ) {
    
    double radius = sqrt(pow(xPos_mm,2.)+pow(yPos_mm,2.));
    double amplitude = 307.9-0.3071*zPos_mm+0.0002257*pow(zPos_mm,2.);
    double shape = 1.1525e-7*sqrt(fabs(zPos_mm-318.84));
    return -shape * pow ( radius, 3. ) + amplitude;
    
  }
  
  //Drift electric field as function of Z in mm
  //For example, use a high-order poly spline
  virtual double FitEF ( double xPos_mm, double yPos_mm, double zPos_mm ) { // in V/cm
    
    return 158.92
      -0.2209000 *pow(zPos_mm,1.)
      +0.0024485 *pow(zPos_mm,2.)
      -8.7098e-6 *pow(zPos_mm,3.)
      +1.5049e-8 *pow(zPos_mm,4.)
      -1.0110e-11*pow(zPos_mm,5.);
    
  }
  
  //S2 PDE custom fit for function of r
  //s2polA + s2polB*r[mm] + s2polC*r^2+... (QE included, for binom dist) e.g.
  virtual double FitS2 ( double xPos_mm, double yPos_mm ) {
    
    double radius = sqrt(pow(xPos_mm,2.)+pow(yPos_mm,2.));
    
    return // unitless, 1.000 at detector center
      9156.3
      +6.22750*pow(radius,1.)
      +0.38126*pow(radius,2.)
      -0.017144*pow(radius,3.)+
      0.0002474*pow(radius,4.)-
      1.6953e-6*pow(radius,5.)+
      5.6513e-9*pow(radius,6.)
      -7.3989e-12*pow(radius,7.);
    
  }
  
  virtual vector<double> FitTBA(double xPos_mm, double yPos_mm,
				double zPos_mm) {
    vector<double> BotTotRat(2);
    
    BotTotRat[0] = 0.650;  // S1 bottom-to-total ratio
    BotTotRat[1] = 0.449;  // S2 bottom-to-total ratio, typically only used for
                           // position recon (1-this)
    
    return BotTotRat;
  }
  
  virtual double OptTrans ( double xPos_mm, double yPos_mm, double zPos_mm ) {
    
    double phoTravT, approxCenter = ( TopDrift + cathode ) / 2., relativeZ = zPos_mm - approxCenter;
    
    double A = 0.048467 - 7.6386e-6 * relativeZ + 1.2016e-6 * pow ( relativeZ, 2. ) - 6.0833e-9 * pow ( relativeZ, 3. );
    if ( A < 0. ) A = 0.; //cannot have negative probability
    double B_a =0.99373 + 0.0010309 * relativeZ - 2.5788e-6 * pow ( relativeZ, 2. ) - 1.2000e-8 * pow ( relativeZ, 3. );
    double B_b = 1. - B_a;
    double tau_a = 11.15; //all times in nanoseconds
    double tau_b = 4.5093 + 0.03437 * relativeZ -0.00018406 * pow ( relativeZ, 2. ) - 1.6383e-6 * pow ( relativeZ, 3. );
    if ( tau_b < 0. ) tau_b = 0.; //cannot have negative time
    
    A = 0.0574; B_a = 1.062; tau_a = 11.1; tau_b = 2.70; B_b = 1.0 - B_a; //LUX D-D conditions
    
    if ( RandomGen::rndm()->rand_uniform() < A )
      phoTravT = 0.; //direct travel time to PMTs (low)
    else { //using P0(t) = A*delta(t)+(1-A)*[(B_a/tau_a)e^(-t/tau_a)+(B_b/tau_b)e^(-t/tau_b)] LUX PSD paper, but should apply to all detectors w/ diff #'s
      if ( RandomGen::rndm()->rand_uniform() < B_a )
	phoTravT = -tau_a * log ( RandomGen::rndm()->rand_uniform() );
      else
	phoTravT = -tau_b * log ( RandomGen::rndm()->rand_uniform() );
    }
    
    double sig= RandomGen::rndm()->rand_gauss(3.84,.09); //includes stat unc but not syst
    phoTravT += RandomGen::rndm()->rand_gauss(0.00,sig); //the overall width added to photon time spectra by the effects in the electronics and the data reduction pipeline
    
    if ( phoTravT > DBL_MAX ) phoTravT = tau_a;
    if ( phoTravT <-DBL_MAX ) phoTravT = 0.000;
    
    return phoTravT; //this function follows LUX (arXiv:1802.06162) not Xe10 technically but tried to make general
  }
  
  virtual vector<double> SinglePEWaveForm ( double area, double t0 ) {
    
    vector<double> PEperBin;
    
    double threshold = PULSEHEIGHT; //photo-electrons
    double sigma = PULSE_WIDTH; //ns
    area *= 10. * ( 1. + threshold );
    double amplitude = area / ( sigma * sqrt ( 2. * M_PI ) ), signal; //assumes perfect Gaussian
    
    double tStep1 = SAMPLE_SIZE/1e2; //ns, make sure much smaller than sample size; used to generate MC-true pulses essentially
    double tStep2 = SAMPLE_SIZE; //ns; 1 over digitization rate, 100 MHz assumed here
    
    double time = -5.*sigma;
    bool digitizeMe = false;
    while ( true ) {
      signal = amplitude * exp(-pow(time,2.)/(2.*sigma*sigma));
      if ( signal < threshold ) {
	if ( digitizeMe ) break;
	else ; //do nothing - goes down to advancing time block
      }
      else {
	if ( digitizeMe )
	  PEperBin.push_back(signal);
	else {
	  if ( RandomGen::rndm()->rand_uniform() < 2.*(tStep1/tStep2) ) {
	    PEperBin.push_back(time+t0);
	     PEperBin.push_back(signal);
	    digitizeMe = true;
	  }
	  else {}
	}
      }
      if ( digitizeMe ) time += tStep2;
      else time += tStep1;
      if ( time > 5.*sigma ) break;
    }
    
    return PEperBin;
    
  }
  
  // Vary VDetector parameters through custom functions
  virtual void ExampleFunction() {
    set_g1 ( 0.1167 );
  }
  
};

#endif
