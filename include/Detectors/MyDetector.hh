
#ifndef DetectorExample_LUX_RUN03_hh
#define DetectorExample_LUX_RUN03_hh 1

#include "VDetector.hh"
using namespace std;

//NOTES: best g1 for DD 0.1193, but for tritium 0.1146; S2 noise 1.9, 7.5%; g1_gas 0.1019, 0.1012
//s2fano 3.6, 0.9; eField in gas 6.25, 6.2; e- life 650, 750 us; fid vol 80-130, 38-305 us; gasGap 4.25, 4.5 mm
//DISCLAIMER: Slight differences from official published values due to private LUX algorithms

class DetectorExample_LUX_RUN03: public VDetector {
  
public:
  
  DetectorExample_LUX_RUN03() {
    if ( verbosity ) cerr << "*** Detector definition message ***" << endl;
    if ( verbosity ) cerr << "You are currently using the LUX Run03 template detector." << endl << endl;
    
    // Call the initialization of all the parameters
    Initialization();
  };
  virtual ~DetectorExample_LUX_RUN03() {};
  
  // Do here the initialization of all the parameters that are not varying as a function of time
  virtual void Initialization() {
    
    // Primary Scintillation (S1) parameters
    g1 = 0.1170; //0.117+/-0.003 WS,0.115+/-0.005 D-D,0.115+/-0.005 CH3T,0.119+/-0.001 LUXSim
    sPEres = 0.37; //arXiv:1910.04211
    sPEthr = (0.3*1.173)/0.915; //arXiv:1910.04211
    sPEeff = 1.00; //arXiv:1910.04211
    noiseB[0] =-0.01; //arXiv:1910.04211
    noiseB[1] = 0.08; //arXiv:1910.04211
    noiseB[2] = 0.;
    noiseB[3] = 0.;
    P_dphe = 0.173; //arXiv:1910.04211
    
    coinWind= 100;// 1310.8214
    coinLevel=2;  //1512.03506
    numPMTs = 119;// 122 minus 3 off
    
    extraPhot =false; //default
    noiseL[0]=1.4e-2; //1910.04211 p.12, to match 1610.02076 Fig. 8
    noiseL[1]=5.0e-2; //1910.04211 p.12, to match 1610.02076 Fig. 8
    
    // Ionization and Secondary Scintillation (S2) parameters
    g1_gas = 0.1016; //0.1 in 1910.04211
    s2Fano = 2.2; //3.7 in 1910.04211; this matches 1608.05381 better
    s2_thr = 165.;//(150.*1.173)/0.915; //65-194 pe in 1608.05381
    E_gas = 6.23; //6.55 in 1910.04211
    eLife_us = 800.; //p.44 of James Verbus PhD thesis Brown
    
    // Thermodynamic Properties
    inGas = false; //duh
    T_Kelvin = 173.; //1910.04211
    p_bar = 1.57; //1910.04211
    
    // Data Analysis Parameters and Geometry
    dtCntr = 160.; //p.61 Dobi thesis UMD, 159 in 1708.02566
    dt_min = 38.; //1608.05381
    dt_max = 305.; //1608.05381
    
    radius = 200.; //1512.03506
    radmax = 235.; //1910.04211
    
    TopDrift = 544.8; //544.95 in 1910.04211
    anode = 549.2; //1910.04211 and 549 in 1708.02566
    gate = 539.2; //1910.04211 and 539 in 1708.02566
    cathode = 55.90; //55.9-56 in 1910.04211,1708.02566
    
    // 2-D (X & Y) Position Reconstruction
    PosResExp = 0.015; //arXiv:1710.02752 indirectly
    PosResBase = 70.8364; //1710.02752 indirectly
  }
  
  // S1 PDE custom fit for function of xyz
  // 1712.05696 indirectly, 1708.02566 Figure 10 color map
  virtual double FitS1 ( double xPos_mm, double yPos_mm, double zPos_mm, LCE map ) {
    
    double radius = sqrt(pow(xPos_mm,2.)+pow(yPos_mm,2.));
    double amplitude = 307.9-0.3071*zPos_mm+0.0002257*pow(zPos_mm,2.);
    double shape = 1.1525e-7*sqrt(fabs(zPos_mm-318.84));
    return -shape * pow ( radius, 3. ) + amplitude;
    
  }
  
  // Drift electric field as function of Z in mm
  // 1709.00095, 1904.08979, 1708.02566 Fig. 13
  virtual double FitEF ( double xPos_mm, double yPos_mm, double zPos_mm ) { // in V/cm
    
    return 158.92 // NOTE: DO NOT JUST RETURN A CONSTANT, THAT IS A SILLY USE of FitEF
      -0.2209000 *pow(zPos_mm,1.)
      +0.0024485 *pow(zPos_mm,2.)
      -8.7098e-6 *pow(zPos_mm,3.)
      +1.5049e-8 *pow(zPos_mm,4.)
      -1.0110e-11*pow(zPos_mm,5.);
    
  }
  
  // S2 PDE custom fit for function of r
  // 1712.05696 & 1710.02752 indirectly. Fig. 13 1708.02566
  virtual double FitS2 ( double xPos_mm, double yPos_mm, LCE map ) {
    
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
    
    double radSq = ( pow(xPos_mm,2.) + pow(yPos_mm,2.) ) / 1e2;
    double TBAzS1 = -0.853 + 0.00925 * ( zPos_mm / 10. );
    double TBArS2 = 0.126+0.000545*radSq-1.90e-6*radSq*radSq+1.20e-9*radSq*radSq*radSq;
    
    if ( TBAzS1 < -1. ) TBAzS1 = -1.;
    if ( TBAzS1 > 1.0 ) TBAzS1 = 1.0;
    if ( TBArS2 < -1. ) TBArS2 = -1.;
    if ( TBArS2 > 1.0 ) TBArS2 = 1.0;
    
    BotTotRat[0] = (1.-TBAzS1)/2.;  // 1712.05696
    BotTotRat[1] = 0.449;//(1.-TBArS2)/2.;  // 1712.05696 and 1710.02752
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
    
    //A = 0.0574; B_a = 1.062; tau_a = 11.1; tau_b = 2.70; B_b = 1.0 - B_a; //LUX D-D conditions
    
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
    
    return phoTravT; //this function follows LUX (arXiv:1802.06162)
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
  virtual void ExampleFunction() { set_g1(0.1167); }
  virtual void ExampleFunction2() { set_molarMass(132.); }
};

#endif
