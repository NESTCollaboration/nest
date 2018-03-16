
#include "NEST.hh"
#include <detector.hh>
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
  
  xyTry[0]= xTest;
  xyTry[1]= yTest;
  
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
  result.quanta=GetQuanta(result.yields,density);
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

QuantaResult NESTcalc::GetQuanta ( YieldResult yields, double density ) {
  
  QuantaResult result;
  int Nq_actual, Ne, Nph, Ni, Nex;
  
  double NexONi = yields.ExcitonRatio; double alf = 1./(1.+NexONi);
  double Nq_mean = yields.PhotonYield + yields.ElectronYield;
  
  if ( yields.Lindhard == 1. ) {
    
    double Fano = 0.12707-0.029623*density- //Fano factor is  << 1
      0.0057042*pow(density,2.)+ //~0.1 for GXe w/ formula from Bolotnikov et al. 1995
      0.0015957*pow(density,3.); //to get it to be ~0.03 for LXe (E Dahl Ph.D. thesis)
    Nq_actual = int(floor(rand_gauss(Nq_mean,sqrt(Fano*Nq_mean))+0.5));
    if ( Nq_actual < 0 ) Nq_actual = 0;
    
    Ni = BinomFluct(Nq_actual,alf);
    Nex= Nq_actual - Ni;
    
  }
  
  else {
    
    Ni = poisson_draw(Nq_mean*alf);
    Nex= poisson_draw(Nq_mean*NexONi*alf);
    Nq_actual = Nex + Ni;
    
  }
  
  if ( Nex < 0 ) Nex = 0;
  if ( Ni < 0 ) Ni = 0;
  if ( Nex > Nq_actual ) Nex = Nq_actual;
  if ( Ni > Nq_actual ) Ni = Nq_actual;
  
  double elecFrac = yields.ElectronYield / Nq_mean;
  if ( elecFrac > 1. ) elecFrac = 1.;
  if ( elecFrac < 0. ) elecFrac = 0.;
  
  double recombProb = 1.-(NexONi+1.)*elecFrac;
  if ( recombProb < 0. ) recombProb = 0.;
  if ( recombProb > 1. ) recombProb = 1.;
  
  double cc = 0.070, bb = 0.530;
  double aa = cc/pow(1.-bb,2.);
  double omega = -aa*pow(recombProb-bb,2.)+cc; if(omega<0.)omega=0.;
  
  if ( yields.Lindhard < 1. ) omega = 0.03;
  double Variance = recombProb*(1.-recombProb)*Ni+omega*omega*Ni*Ni;
  Ne = int(floor(rand_gauss((1.-recombProb)*Ni,sqrt(Variance))+0.5));
  if ( Ne < 0 ) Ne = 0;
  if ( Ne > Ni) Ne =Ni;
  
  Nph = Nq_actual - Ne;
  if ( Nph > Nq_actual ) Nph = Nq_actual;
  if ( Nph < Nex ) Nph = Nex;
  
  if ( (Nph+Ne) != (Nex+Ni) ) cout << "\nERROR: Quanta not conserved. Tell Matthew Immediately!\n";
  
  result.photons =Nph;
  result.electrons=Ne;
  
  return result; //quanta returned with recomb fluctuations
  
}

YieldResult NESTcalc::GetYields ( INTERACTION_TYPE species, double energy, double density, double dfield ) {
  
  double massNum = 4.; //TODO: make this flexible
  const double m3 = 2., m4 = 2., m6 = 0.;
  double Ne = -999; double Nph = -999; double NexONi = -999; double m8 = 2., L = 1.;
  const double deltaT_ns_halflife = 154.4;
  
  double Wq_eV = 1.9896 + (20.8 - 1.9896) / (1. + pow(density / 4.0434, 1.4407));
  double alpha = 0.067366 + density * 0.039693, Ni, recombProb, Nq, Ly, Qy, ThomasImel;
  switch ( species ) {
  case NR:
  case WIMP:
  case B8:
  case DD:
  case AmBe:
  case Cf:
    {
      Nq = 12.6*pow(energy,1.05);
      ThomasImel = 0.0522*pow(dfield,-0.0694)*pow(density/2.9,0.3);
      Qy = 1. / (ThomasImel*sqrt(energy+9.75));
      Ly = Nq / energy - Qy;
      Ne = Qy * energy;
      Nph= Ly * energy;
      NexONi = 1.00*erf(0.01*energy);
      L = ( Nq / energy ) * Wq_eV * 1e-3;
    } break;
  case ion:
    {
      L = 0.96446 / (1. + pow(massNum * massNum / 19227., 0.99199));
      if (massNum == 4) L = 0.56136 * pow(energy, 0.056972);
      ThomasImel = 0.0067 / pow(1. + pow(dfield / 95.768, 8.5673), 0.060318) * pow(density / 2.857, 0.3);
      Nq = 1e3 * L * energy / Wq_eV;
      Ni = Nq / (1. + alpha);
      recombProb = 1. - log(1. + (ThomasImel / 4.) * Ni) / ((ThomasImel / 4.) * Ni);
      Nph= Nq * alpha / (1. + alpha) + recombProb*Ni;
      Ne = Nq - Nph;
    } break;
  case gammaRay:
    {
      double m1 = 33.951 + (3.3284 - 33.951) / (1. + pow(dfield / 165.34, .72665));
      double m2 = 1000 / Wq_eV;
      double m5 = 23.156 + (10.737 - 23.156) / (1. + pow(dfield / 34.195, .87459));
      double densCorr = 240720. / pow ( density, 8.2076 );
      double m7 = 66.825 + (829.25 - 66.825) / (1. + pow(dfield /densCorr,.83344));
      Nq = energy * 1000. / Wq_eV;
      if ( density < 1. ) m8 = -2.;
      Qy = m1 + (m2 - m1) / (1. + pow(energy / m3, m4)) + m5 + (m6 - m5) / (1. + pow(energy / m7, m8));
      Ly = Nq / energy - Qy;
      Ne = Qy * energy;
      Nph =Ly * energy;
      NexONi = alpha*erf(0.05*energy);
    } break;
  case Kr83m:
    {
      if ( energy == 9.4 ) {
	double deltaT_ns = rand_exponential ( deltaT_ns_halflife );
	Nq = energy * ( 1e3 / Wq_eV + 6.5 );
	double medTlevel = 47.8 + ( 69.201 - 47.8 ) / pow ( 1. + pow ( dfield / 250.13, 0.9 ), 1. );
	double highTrise = 1.15 + ( 1. - 1.15 ) / ( 1. + pow ( deltaT_ns / 1200., 18. ) );
	double lowTdrop = 14. * pow ( dfield, 0.19277 );
	printf ( "%.6f\t", deltaT_ns );
	Nph = energy*highTrise*(5.1e4*pow(2.*deltaT_ns+10.,-1.5)+medTlevel)/(1.+pow(deltaT_ns/lowTdrop,-3.));
	alpha = 0.;
      }
      else {
	Nq = energy * 1000. / Wq_eV;
	Nph= energy * ( 6. + ( 69.742 - 6. ) / pow ( 1. + pow ( dfield / 9.515, 1.9 ), 0.063032 ) );
      }
      Ne = Nq - Nph;
      NexONi = alpha*erf(0.05*energy);
    } break;
  default: //beta, CH3T
    {
      double QyLvllowE = 1e3/Wq_eV+6.5*(1.-1./(1.+pow(dfield/47.408,1.9851)));
      double QyLvlmedE =32.988-32.988/(1.+pow(dfield/(0.026715*exp(density/0.33926)),0.6705));
      double DokeBirks = 1652.264+(1.415935e10-1652.264)/(1.+pow(dfield/0.02673144,1.564691));
      Nq = energy * 1e3 / ( Wq_eV+(12.578-Wq_eV)/(1.+pow(energy/1.6,3.5)) );
      double LET_power = -2.;
      if ( density < 1. ) LET_power = 2.;
      double QyLvlhighE =28.;
      if ( density > 3. ) QyLvlhighE=49.;
      Qy = QyLvlmedE+(QyLvllowE-QyLvlmedE)/pow(1.+1.304*pow(energy,2.1393),0.35535)+QyLvlhighE/(1.+DokeBirks*pow(energy,LET_power));
      Ly = Nq / energy - Qy;
      Ne = Qy * energy;
      Nph= Ly * energy;
      NexONi = alpha*erf(0.05*energy);
    } break;
  }
  
  assert(Ne!=-999 && Nph!=-999
	 && NexONi!=-999);
  if ( Nph> energy / 7e-3 ) Nph= energy / 7e-3; //yields can never exceed 1 / [ W ~ 7 eV ]
  if ( Ne > energy / 7e-3 ) Ne = energy / 7e-3;
  if ( Nph < 0. ) Nph = 0.; if ( Ne < 0. ) Ne = 0.;
  if ( NexONi < 0. ) NexONi = 0.;
  if ( L < 0. ) L = 0.;
  if ( L > 1. ) L = 1.; //Lindhard Factor
  
  YieldResult result;
  result.PhotonYield = Nph;
  result.ElectronYield=Ne;
  result.ExcitonRatio =NexONi;
  result.Lindhard = L;
  return result; //everything needed to calculate fluctuations
  
}

void NESTcalc::SetRandomSeed ( unsigned long int s ) {
  rng.seed(s);
}

NESTcalc::NESTcalc ( ) {
  rng.seed(0);
}

vector<double> NESTcalc::GetS1 ( int Nph, double dz, double driftVelocity ) {
  
  vector<double> scintillation(8);  // return vector
  
  // Add some variability in g1 drawn from a polynomial spline fit
  double posDep = s1poly[0] + s1poly[1] * dz +
    s1poly[2] * pow(dz,2.)+
    s1poly[3] * pow(dz,3.)+
    s1poly[4] * pow(dz,4.);
  double dz_center = liquidBorder - driftVelocity * dtCntr; //go from t to z
  posDep /= s1poly[0]+
    s1poly[1] * pow(dz_center,1.)+
    s1poly[2] * pow(dz_center,2.)+
    s1poly[3] * pow(dz_center,3.)+
    s1poly[4] * pow(dz_center,4.); // Z is always in mm now never cm
  
  // generate a number of PMT hits drawn from a binomial distribution. Initialize number of photo-electrons
  int nHits=BinomFluct(Nph,g1*posDep), Nphe = 0;
  
  // Initialize the pulse area and spike count variables
  double pulseArea = 0., spike = 0., prob;
  
  // If single photo-electron efficiency is under 1 and the threshold is above 0 (some phe will be below threshold)
  if ( sPEthr > 0. && nHits < numPMTs ) {
    // Step through the pmt hits
    for ( int i = 0; i < nHits; i++ ) {
      // generate photo electron, integer count and area
      double phe1 = rand_gauss(1.,sPEres) + rand_gauss(noise[0],noise[1]); Nphe++;
      prob = rand_uniform();
      // zero the area if random draw determines it wouldn't have been observed.
      if ( prob > sPEeff ) { phe1 = 0.; } //add an else with Nphe++ if not doing mc truth
      // Generate a double photo electron if random draw allows it
      double phe2 = 0.;
      if ( rand_uniform() < P_dphe ) {
	// generate area and increment the photo-electron counter
	phe2 = rand_gauss(1.,sPEres) + rand_gauss(noise[0],noise[1]); Nphe++;
	// zero the area if phe wouldn't have been observed
	if ( rand_uniform() > sPEeff && prob > sPEeff ) { phe2 = 0.; } //add an else with Nphe++ if not doing mc truth
	// The dphe occurs simultaneously to the first one from the same source photon. If the first one is seen, so should be the second one
      }
      // Save the phe area and increment the spike count (very perfect spike count) if area is above threshold
      if ( (phe1+phe2) > sPEthr ) { spike++; pulseArea += phe1 + phe2; }
    }
  }
  else { // apply just an empirical efficiency by itself, without direct area threshold
    Nphe = nHits + BinomFluct(nHits,P_dphe); double eff=sPEeff; if(nHits>=numPMTs) eff=1.;
    pulseArea = rand_gauss(BinomFluct(Nphe,1.-(1.-eff)/(1.+P_dphe)),sPEres*sqrt(Nphe));
    spike = (double)nHits;
  }
  if ( pulseArea < 0. ) pulseArea = 0.;
  double pulseAreaC= pulseArea / posDep;
  double Nphd = pulseArea / (1.+P_dphe);
  double NphdC= pulseAreaC/ (1.+P_dphe);
  double spikeC = spike / posDep;
  
  scintillation[0] = nHits; scintillation[1] = Nphe;
  scintillation[2] = pulseArea; scintillation[3] = pulseAreaC;
  scintillation[4] = Nphd; scintillation[5] = NphdC;
  scintillation[6] = spike; scintillation[7] = spikeC;
  
  if ( spike > 10 )
    prob = 1.;
  else {
    if ( spike >= coinLevel ) { double numer, denom;
      for ( int i = spike; i > 0; i-- ) {
	denom += nCr ( numPMTs, i );
	if ( i >= coinLevel ) numer += nCr ( numPMTs, i );
      }
      prob = numer / denom;
    }
    else
      prob = 0.;
  }
  
  if ( rand_uniform() < prob ) // coincidence has to happen in different PMTs
    { ; }
  else { // some of these are set to -1 to flag them as having been below threshold
    //scintillation[0] *= -1.;
    //scintillation[1] *= -1.;
    scintillation[2] *= -1.;
    scintillation[3] *= -1.;
    //scintillation[4] *= -1.;
    //scintillation[5] *= -1.;
    scintillation[6] *= -1.;
    scintillation[7] *= -1.;
  }
  
  return scintillation;
  
}

vector<double> NESTcalc::GetS2 ( int Ne, double dt ) {
  
  vector<double> ionization(8);
  double alpha = 0.137, beta = 177., gamma = 45.7;
  double epsilon = 1.85 / 1.00126;
  
  double E_liq = E_gas / epsilon; //kV per cm
  double ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645; // arXiv:1710.11032
  if ( ExtEff > 1. ) ExtEff = 1.;
  if ( ExtEff < 0. ) ExtEff = 0.;
  int Nee = BinomFluct(Ne,ExtEff*exp(-dt/eLife_us));
  
  double elYield = Nee*
    (alpha*E_gas*1e3-beta*p_bar-gamma)*
    gasGap_mm*0.1; // arXiv:1207.2292
  int Nph = int(floor(rand_gauss(elYield,sqrt(s2Fano*elYield))+0.5));
  int nHits = BinomFluct(Nph,g1_gas);
  int Nphe = nHits + BinomFluct(nHits,P_dphe);
  double pulseArea=rand_gauss(Nphe,sPEres*sqrt(Nphe));
  double pulseAreaC= pulseArea/exp(-dt/eLife_us);
  double Nphd = pulseArea / (1.+P_dphe);
  double NphdC= pulseAreaC/ (1.+P_dphe);
  
  double S2b = rand_gauss(S2botTotRatio*pulseArea,sqrt(S2botTotRatio*pulseArea*(1.-S2botTotRatio)));
  double S2bc= S2b / exp(-dt/eLife_us); // for detectors using S2 bottom-only in their analyses
  
  ionization[0] = Nee; ionization[1] = Nph;
  ionization[2] = nHits; ionization[3] = Nphe;
  if ( s2_thr >= 0 ) {
    ionization[4] = pulseArea; ionization[5] = pulseAreaC;
    ionization[6] = Nphd; ionization[7] = NphdC;
  }
  else {
    ionization[4] = S2b; ionization[5] = S2bc;
    ionization[6] = S2b / (1.+P_dphe); ionization[7] = S2bc / (1.+P_dphe);
  }
  
  if ( Nphd < abs(s2_thr) ) ionization[0] *= -1.;
  return ionization;
  
}



long double NESTcalc::Factorial ( double x ) {
  
  return tgammal ( x + 1. );
  
}

double NESTcalc::nCr ( double n, double r ) {
  
  return Factorial(n) /
    ( Factorial(r) * Factorial(n-r) );
  
}
