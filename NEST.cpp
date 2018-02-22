
#include "NEST.hh"
#include <iostream>
#include <assert.h>

using namespace NEST;
using namespace std;

double NESTcalc::rand_uniform() {

    return (double) (rng() - rng.min()) / (double) (rng.max() - rng.min());

}

double NESTcalc::rand_gauss(double mean, double sigma) {

    double u = rand_uniform(), v = rand_uniform();
    return mean + sigma * sqrt(-2. * log(u)) * cos(2. * M_PI * v);

}

int NESTcalc::poisson_draw(double mean)
{
  std::poisson_distribution<int> distribution(mean);
  return distribution(rng);
}


double NESTcalc::rand_exponential(double half_life) {

    double r = rand_uniform();
    return log(1 - r) * -1 * half_life / log(2.);

}

vector<double> NESTcalc::VonNeumann(double xMin, double xMax, double yMin,double yMax,
				    double xTest,double yTest,double fValue){
  
  vector<double> xyTry(3);
  
  xyTry[0]= xTest ;
  xyTry[1]= yTest ;
  
  if ( xyTry[1] > fValue ) {
    xyTry[0] = xMin+(xMax-xMin)*rand_uniform() ;
    xyTry[1] = yMin+(yMax-yMin)*rand_uniform() ;
    xyTry[2] = 1. ;
  }
  else
    xyTry[2] = 0. ;
  
  return xyTry; //doing a vector means you can return 2 values at the same time
  
}

int NESTcalc::BinomFluct(int N0, double prob) {

    double mean = N0*prob;
    double sigma = sqrt(N0 * prob * (1. - prob));
    int N1 = 0;

    if (prob <= 0.00) return N1;
    if (prob >= 1.00) return N0;

    if (N0 < 10) {
        for (int i = 0; i < N0; i++) {
            if (rand_uniform() < prob) N1++;
        }
    } else {
        N1 = int(floor(rand_gauss(mean, sigma) + 0.5));
    }

    if (N1 > N0) N1 = N0;
    if (N1 < 0) N1 = 0;

    return N1;

}

NESTresult NESTcalc::FullCalculation(INTERACTION_TYPE species, double energy, double density, double dfield){
    NESTresult result;
    result.yields = GetYields(species,energy,density,dfield);
    result.quanta = GetQuanta(result.yields);
    result.photon_times = GetPhotonTimes(/*stuff*/);
    return result;
}

double NESTcalc::PhotonTime(INTERACTION_TYPE species, bool exciton){
    //old code put here by Jason
    //times in ns
    double return_time=0;
    double tau1 = rand_gauss(3.1,.7); //err from wgted avg.
    double tau3 = rand_gauss(24.,1.); //ibid.
    //these singlet and triplet times may not be the ones you're
    //used to, but are the world average: Kubota 79, Hitachi 83 (2
    //data sets), Teymourian 11, Morikawa 89, and Akimov '02
    double SingTripRatioX, SingTripRatioR;
    if(species==beta || species==Kr83m || species == gammaRay){ //these classes are questionable
        //disregard tauR from original model--it's very small for any electric field.
        SingTripRatioX = rand_gauss(0.17,0.05);
        SingTripRatioR = rand_gauss(0.8, 0.2);
    }
    else if(species==ion){ //these classes are questionable
        SingTripRatioR = rand_gauss(2.3,0.51);
        SingTripRatioX = SingTripRatioR;
    }
    else{//NR //these classes are questionable
        SingTripRatioR = rand_gauss(7.8,1.5);
        SingTripRatioX = SingTripRatioR;
    }
    if(exciton){
        if (rand_uniform() < SingTripRatioR / (1 + SingTripRatioR))
            return_time = tau1 * -log(rand_uniform());
        else return_time = tau3 * -log(rand_uniform());
    } else {
        if (rand_uniform() < SingTripRatioX / (1 + SingTripRatioX))
            return_time = tau1 * -log(rand_uniform());
        else return_time = tau3 * -log(rand_uniform());
    }
    
    return return_time;
}

photonstream NESTcalc::GetPhotonTimes(/*inputs*/){
    //TODO by MATTHEW
    photonstream return_photons;
    return_photons.push_back(PhotonTime(beta,true));//example line, modify
    return return_photons;
    
}


QuantaResult NESTcalc::GetQuanta(YieldResult yields){
    //TODO by MATTHEW
    QuantaResult result;
    result.photons = yields.PhotonYield; //fake, fix this
    result.electrons = yields.ElectronYield;
    return result;
}

YieldResult NESTcalc::GetYields ( INTERACTION_TYPE species, double energy, double density, double dfield ) {
  
  double massNum = 4.; //TODO: make this flexible
  const double m3 = 2., m4 = 2., m6 = 0.;
  double Ne = -999; double Nph = -999, m8 = 2.;
  const double deltaT_ns_halflife = 154.4;
  
  double Wq_eV = 1.9896 + (20.8 - 1.9896) / (1. + pow(density / 4.0434, 1.4407));
  const double alpha = 0.067366 + density * 0.039693;
  switch ( species ) {
  case NR:
  case WIMP:
  case B8:
  case DD:
  case AmBe:
  case Cf:
    {
      double Nq = 12.6*pow(energy,1.05);
      double TIB = 0.0522*pow(dfield,-0.0694)*pow(density/2.9,0.3);
      double Qy = 1. / (TIB*sqrt(energy+9.75));
      double Ly = Nq / energy - Qy;
      Ne = Qy * energy;
      Nph= Ly * energy;
      double NexONi = 1.00*erf(0.01*energy);
      double Ni = Nq/(1.+NexONi); double Nex = Nq - Ni;
      double elecFrac = Ne / Nq;
      double recombProb = 1.-(NexONi+1.)*elecFrac;
    } break;
        case ion:
        {
            double L = 0.96446 / (1. + pow(massNum * massNum / 19227., 0.99199));
            if (massNum == 4) L = 0.56136 * pow(energy, 0.056972);
            double ThomasImel = 0.0067 / pow(1. + pow(dfield / 95.768, 8.5673), 0.060318) * pow(density / 2.857, 0.3);
            double totQ = 1e3 * L * energy / Wq_eV;
            double Ni = totQ / (1. + alpha);
            double recombProb = 1. - log(1. + (ThomasImel / 4.) * Ni) / ((ThomasImel / 4.) * Ni);
            Nph = totQ * alpha / (1. + alpha) + recombProb*Ni;
            Ne = totQ - Nph;
        } break;
  case gammaRay:
    {
      double m1 = 33.951 + (3.3284 - 33.951) / (1. + pow(dfield / 165.34, .72665));
      double m2 = 1000 / Wq_eV;
      double m5 = 23.156 + (10.737 - 23.156) / (1. + pow(dfield / 34.195, .87459));
      double densCorr = 240720. / pow ( density, 8.2076 );
      double m7 = 66.825 + (829.25 - 66.825) / (1. + pow(dfield /densCorr,.83344));
      double totQ = energy * 1000. / Wq_eV;
      if ( density < 1. ) m8 = -2.;
      double Qy = m1 + (m2 - m1) / (1. + pow(energy / m3, m4)) + m5 + (m6 - m5) / (1. + pow(energy / m7, m8));
      double Ly = totQ/energy-Qy;
      Ne = Qy * energy;
      Nph =Ly * energy;
      double NexONi = alpha*erf(0.05*energy);
      double Ni = totQ/(1.+NexONi); double Nex = totQ - Ni;
      double elecFrac = Ne / totQ;
      double recombProb = 1.-(NexONi+1.)*elecFrac;
    } break;
        case Kr83m:
        {
            double totQ = energy * 1000. / Wq_eV;
            if (energy == 9.4) {
                double m1 = 99678. - 21574. * log10(dfield);
                double m2a = 47.364 + (131.69 - 47.364) / pow(1. + pow(dfield / 71.368, 2.4130), 0.060318);
                double deltaT_ns = rand_exponential(deltaT_ns_halflife);
		printf("%.6f\t",deltaT_ns);
                Nph = (m1 * pow(2. * deltaT_ns + 10., -1.5) + m2a) / 0.21 / 1.05;
            }
            else{
                Nph = energy * (0.51987 + (1.0036 - 0.51987) / (1. + pow(dfield / 309.98, 1.0844)))*65.5/1.05;
            }
            Ne = totQ - Nph;
        } break;
  default: //beta, CH3T
    {
      double QyLvllowE = 1e3/Wq_eV+6.5*(1.-1./(1.+pow(dfield/47.408,1.9851)));
      double QyLvlmedE =32.988-32.988/(1.+pow(dfield/(0.026715*exp(density/0.33926)),0.6705));
      double DokeBirks = 1652.264+(1.415935e10-1652.264)/(1.+pow(dfield/0.02673144,1.564691));
      double totQuanta = energy * 1e3 / ( Wq_eV+(12.578-Wq_eV)/(1.+pow(energy/1.6,3.5)) );
      double LET_power = -2.;
      if ( density < 1. ) LET_power = 2.;
      double QyLvlhighE =28.;
      if ( density > 3. ) QyLvlhighE=49.;
      double Qy = QyLvlmedE+(QyLvllowE-QyLvlmedE)/pow(1.+1.304*pow(energy,2.1393),0.35535)+QyLvlhighE/(1.+DokeBirks*pow(energy,LET_power));
      double Ly = totQuanta/energy-Qy;
      Ne = Qy * energy;
      Nph= Ly * energy;
      double NexONi = alpha*erf(0.05*energy);
      double Ni = totQuanta/(1.+NexONi); double Nex = totQuanta - Ni;
      double elecFrac = Ne / totQuanta;
      double recombProb = 1.-(NexONi+1.)*elecFrac;
    } break;
    }

    assert(Ne!=-999 && Nph !=-999);
    //if (Nph> m2 * energy) Nph= m2 * energy;
    //if (Ne > m2 * energy) Ne = m2 * energy;
    if (Nph <0.) Nph =0.;
    if (Ne < 0.) Ne = 0.;

    YieldResult result;
    result.PhotonYield=Nph;
    result.ElectronYield=Ne;
    

    return result;

}

void NESTcalc::SetRandomSeed(unsigned long int s){
    rng.seed(s);
}

NESTcalc::NESTcalc() {
    rng.seed(0);
}


