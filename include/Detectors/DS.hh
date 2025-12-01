
#ifndef DetectorExample_LUX_RUN03_hh
#define DetectorExample_LUX_RUN03_hh 1

#include "VDetector.hh"
using namespace std;



class DS_Detector : public VDetector {
 public:
  DS_Detector() {
    Initialization();
  };
  ~DS_Detector() override = default;
  
  
 void Initialization() override {
    name = "DS-50";
   	g1 = 0.16;   
    sPEres = 0.05;   // Personal communication with Shawn
    sPEthr = 0.6;   // 
    sPEeff = 1;   //I take 1 just for the sake of lowering free parameters
    noiseBaseline[0] = 0.00;        
    noiseBaseline[1] = 0.00;         
    noiseBaseline[2] = 0.;          
    noiseBaseline[3] = 0.;  
    P_dphe = 0.02;  

    //also from arXiv:2311.18647
    coinWind = 100;  
    coinLevel = 2;  
    numPMTs = 38;  

    OldW13eV = true;  
    noiseLinear[0] = 0; 
    noiseLinear[1] = 0;  
    
    // Ionization and Secondary Scintillation (S2) parameters
    //these parameters are incorrect - S2 still in pre-alpha mode
    g1_gas = 0.1;  
    s2Fano = 3.61;  
    s2_thr = 300.; 
    //these parameters are incorrect  - S2 still in pre-alpha mode
    
    E_gas = 4.2;    
    eLife_us = 10000.; 
    
    // Thermodynamic Properties
   //right now they're incorrect and just used for reference taking into account the LAr drift velocities haven't been updated yet (you can uncomment them manually, but)
    T_Kelvin = 89.; 
    p_bar = 1.2; 
    
    // Data Analysis Parameters and Geometry
    dtCntr = 177.5;//a bit random
    dt_min = 20; //a bit random
    dt_max = 375.; //calculated from arXiv:2311.18647
    
    //from arXiv:2311.18647
    radius = 230.; 
    radmax = 230.; 
    
    //arXiv:1410.0653
    TopDrift = 306; 
    anode = 316; 
    gate = 256; 
    cathode = 0; 
    
    // 2-D (X & Y) Position Reconstruction - still from LUX
    PosResExp = 0.015; 
    PosResBase = 70.8364; 
  }
  
 double FitS1(double xPos_mm, double yPos_mm, double zPos_mm,
               LCE map) override {
    return 1.0;
  }

  // Drift electric field as function of Z in mm
  double FitEF(double xPos_mm, double yPos_mm,
               double zPos_mm) override {  // in V/cm
    return 200.;
  }

  double FitS2(double xPos_mm, double yPos_mm, LCE map) override { return 1.0; }

  vector<double> FitTBA(double xPos_mm, double yPos_mm,
                        double zPos_mm) override {
    vector<double> BotTotRat(2);

    BotTotRat[0] = 0.6;    // S1 bottom-to-total ratio
    BotTotRat[1] = 0.323;  // S2 bottom-to-total ratio, typically only used for
                           // position recon (1-this)

    return BotTotRat;
  }

  //very preliminary opttrans function for LAr
  double OptTrans(double xPos_mm, double yPos_mm, double zPos_mm) override {
    double phoTravT, approxCenter = (TopDrift + cathode) / 2.,
                     relativeZ = zPos_mm - approxCenter;
    return - 1.001*pow(abs(relativeZ),1.6142)+1.0008*pow(abs(relativeZ),1.6142)+5.04080;
  }
//unchanged from LUX
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
};

#endif
