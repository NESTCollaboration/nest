
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




YieldResult NESTcalc::GetYields(INTERACTION_TYPE species, double energy, double density, double dfield) {

  
    double massNum = 4.;//TODO: make this flexible
    const double m2 = 78.324, m3 = 2., m4 = 2., m6 = 0., m8 = 2.;
    const double deltaT_ns_halflife = 154.4;
    double Ne = -999; double Nph=-999;

    const double Wq_eV = 1.9896 + (20.8 - 1.9896) / (1. + pow(density / 4.0434, 1.4407));
    double alpha = 0.067366 + density * 0.039693;
    switch (species) {
        case NR:
        case WIMP:
        case B8:
        case DD:
        case AmBe:
        case Cf:
        {
            double epsilon = 11.5 * energy * pow(54., (-7. / 3.));
            double densCorr = 13.7 / Wq_eV;
            double peter = 0.;
            double power = 0.28884 * pow(dfield, -0.045639);
            double Qy = 8.494 / pow(1. + pow(energy / 5.1215, 1.671), power);
            double totQ = 12.256 * pow(energy, 1.0770) - peter / epsilon;
            Nph = densCorr * (totQ - energy * Qy);
            double Ly = Nph / energy;
            Qy = Qy + (totQ / energy - (Ly + Qy));
            Ne = Qy * energy;
            double a = 0.393821, b = 0.498445;//new coeff for bi-exciton quenching
            Nph *= 1. / (1. + a * pow(epsilon, b));
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
            double m5 = 23.156 + (10.737 - 23.156) / (1. + pow(dfield / 34.195, .87459));
            double m7 = 66.825 + (829.25 - 66.825) / (1. + pow(dfield / 43.608, .83344));
            double densCorr = 0.32856 + 0.23187 * density;
            double totQ = energy * 1000. / Wq_eV;
            double Qy = m1 + (m2 - m1) / (1. + pow(energy / m3, m4)) + m5 + (m6 - m5) / (1. + pow(energy / m7, m8));
            Qy /= 73. * Wq_eV * 1e-3;
            double Ly = (73. - Qy) * densCorr / (73. * Wq_eV * 0.001);
            Qy = Qy + (totQ / energy - (Ly + Qy));
            Ne = Qy * energy;
            Nph = Ly * energy;
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
            double m1 = 37.609 + (9.4398 - 37.609) / (1. + pow(dfield / 95.192, .65711));
            double m5 = 24.159;
            double m7 = 27.663 + (694.46 - 27.663) / (1. + pow(dfield / 28.595, .84217));
            double densCorr = 0.32856 + 0.23187 * density;
            double totQ = energy * 1000. / Wq_eV;
            double Qy = m1 + (m2 - m1) / (1. + pow(energy / m3, m4)) + m5 + (m6 - m5) / (1. + pow(energy / m7, m8));
            Qy /= 73. * Wq_eV * 1e-3;
            double Ly = (73. - Qy) * densCorr / (73. * Wq_eV * 0.001);
            Qy = Qy + (totQ / energy - (Ly + Qy));
            Ne = Qy * energy;
            Nph = Ly * energy;
        } break;
    }

    assert(Ne!=-999 && Nph !=-999);
    if (Nph> m2 * energy) Nph= m2 * energy;
    if (Ne > m2 * energy) Ne = m2 * energy;
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


