
#include "NEST.hh"
//#include "/Volumes/USB20FD/z_git/NEST/LUX_Run03.hh"
#include "detector.hh"

using namespace NEST;
using namespace std;
using namespace DetectorEffects;

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

long NESTcalc::BinomFluct(long N0, double prob) {
  
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

NESTresult NESTcalc::FullCalculation(INTERACTION_TYPE species,double energy,double density,double dfield,
		  double A,double Z,vector<double> NuisParam){
  
  NESTresult result;
  result.yields = GetYields(species,energy,density,dfield,A,Z,NuisParam);
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
  
  double NexONi = yields.ExcitonRatio, Fano = 1.;
  double alf = 1./(1.+NexONi);
  double Nq_mean = yields.PhotonYield + yields.ElectronYield;
  
  if ( yields.Lindhard == 1. ) {
    
    Fano = 0.12707-0.029623*density- //Fano factor is  << 1
      0.0057042*pow(density,2.)+ //~0.1 for GXe w/ formula from Bolotnikov et al. 1995
      0.0015957*pow(density,3.); //to get it to be ~0.03 for LXe (E Dahl Ph.D. thesis)
    Nq_actual = int(floor(rand_gauss(Nq_mean,sqrt(Fano*Nq_mean))+0.5));
    if ( Nq_actual < 0 || Nq_mean == 0. ) Nq_actual = 0;
    
    Ni = BinomFluct(Nq_actual,alf);
    Nex= Nq_actual - Ni;
    
  }
  
  else {
    
    Ni = int(floor(rand_gauss(Nq_mean*alf,sqrt(Fano*Nq_mean*alf))+0.5)); if(Ni<0)Ni=0;
    Nex= int(floor(rand_gauss(Nq_mean*NexONi*alf,sqrt(Fano*Nq_mean*NexONi*alf))+0.5)); if(Nex<0)Nex=0;
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
  
  double ef = yields.ElectricField;
  double cc = 0.3+(2.419110e-2-0.3)/(1.+pow(ef/1.431556e4,0.5));
  double bb = 0.54;
  double aa = cc/pow(1.-bb,2.);
  double omega = -aa*pow(recombProb-bb,2.)+cc;
  if ( omega < 0.0 ) omega = 0.0;
  
  if ( yields.Lindhard < 1. )
    omega = 0.04*exp(-pow(elecFrac-0.5,2.)/0.17);
  double Variance = recombProb*(1.-recombProb)*Ni+omega*omega*Ni*Ni;
  Ne=int(floor(rand_gauss((1.-recombProb)*Ni,sqrt(Variance))+0.5));
  if ( Ne < 0 ) Ne = 0;
  if ( Ne > Ni) Ne =Ni;
  
  Nph = Nq_actual - Ne;
  if ( Nph > Nq_actual ) Nph = Nq_actual;
  if ( Nph < Nex ) Nph = Nex;
  
  if ( (Nph+Ne) != (Nex+Ni) )
    cerr << "\nERROR: Quanta not conserved. Tell Matthew Immediately!\n";
  
  result.photons =Nph;
  result.electrons=Ne;
  
  return result; //quanta returned with recomb fluctuations
  
}

YieldResult NESTcalc::GetYields ( INTERACTION_TYPE species, double energy, double density,
				  double dfield, double massNum, double atomNum, vector<double> NuisParam ) {
  
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
  case Cf: //this doesn't mean all NR is Cf, this is like a giant if statement. Same intrinsic yields, but different energy spectra (TestSpectra)
    {
      int massNumber; double ScaleFactor[2] = { 1., 1. };
      if ( massNum != 0. ) massNumber = int(massNum);
      else massNumber = SelectRanXeAtom ( rand_uniform() * 100.0 );
      ScaleFactor[0] = sqrt(MOLAR_MASS/(double)massNumber);
      ScaleFactor[1] = ScaleFactor[0];
      Nq = 12.6 * pow ( energy, 1.05 );
      ThomasImel = 0.0522*pow(dfield,-0.0694)*pow(density/2.9,0.3);
      Qy = 1. / (ThomasImel*sqrt(energy+9.75));
      Ly = Nq / energy - Qy;
      Ne = Qy * energy * ScaleFactor[1] * NuisParam[1];
      Nph= Ly * energy * ScaleFactor[0] * NuisParam[0];
      NexONi = 1.00*erf(0.01*energy); Nq = Nph + Ne;
      L = ( Nq / energy ) * Wq_eV * 1e-3;
    } break;
  case ion:
    {
      double A1 = massNum, A2 = SelectRanXeAtom(rand_uniform()*100.);
      double Z1 = atomNum, Z2 = ATOM_NUM;
      double Z_mean = pow(pow(Z1,(2./3.))+pow(Z2,(2./3.)),1.5);
      double E1c = pow(A1,3.)*pow(A1+A2,-2.)*pow(Z_mean,(4./3.))*pow(Z1,(-1./3.))*500.;
      double E2c = pow(A1+A2,2.)*pow(A1,-1.)*Z2*125.;
      double gamma = 4. * A1 * A2 / pow ( A1 + A2, 2. );
      double Ec_eV = gamma * E2c;
      double Constant = (2./3.)*(1./sqrt(E1c)+0.5*sqrt(gamma/Ec_eV));
      L = Constant * sqrt ( energy * 1e3 );
      double L_max = 0.96446 / ( 1.+pow ( massNum * massNum / 19227., 0.99199 ) );
      if ( atomNum == 2. && massNum == 4. ) L = 0.56136 * pow ( energy, 0.056972 );
      if ( L > L_max ) L = L_max;
      double densDep = pow(density/0.2679,-2.3245);
      double massDep = 0.02966094*exp(0.17687876*(massNum/4.-1.))+1.-0.02966094;
      double fieldDep= pow ( 1.+pow ( dfield/95., 8.7 ), 0.0592 );
      if ( density < 1. ) fieldDep = sqrt ( dfield );
      ThomasImel = 0.00625 * massDep / ( 1. + densDep ) / fieldDep;
      Wq_eV = 28.259+25.667*log10(density)
	-33.611*pow(log10(density),2.)
	-123.73*pow(log10(density),3.)
	-136.47*pow(log10(density),4.)
	-74.194*pow(log10(density),5.)
	-20.276*pow(log10(density),6.)
	-2.2352*pow(log10(density),7.);
      alpha = 0.64 / pow ( 1. + pow ( density / 10., 2. ), 449.61 );
      NexONi = alpha + 0.00178 * pow ( atomNum, 1.587 );
      Nq = 1e3 * L * energy / Wq_eV;
      Ni = Nq / ( 1. + NexONi );
      recombProb = 1. - log(1. + (ThomasImel / 4.) * Ni) / ((ThomasImel / 4.) * Ni);
      Nph = Nq * NexONi / ( 1. + NexONi ) + recombProb * Ni; Ne = Nq - Nph;
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
      Nq = energy * 1e3 / Wq_eV;//( Wq_eV+(12.578-Wq_eV)/(1.+pow(energy/1.6,3.5)) );
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
  
  assert(Ne!=-999 && Nph!=-999 && NexONi!=-999);
  if ( Nph> energy / 7e-3 ) Nph= energy / 7e-3; //yields can never exceed 1 / [ W ~ 7 eV ]
  if ( Ne > energy / 7e-3 ) Ne = energy / 7e-3;
  if ( Nph < 0. ) Nph = 0.; if ( Ne < 0. ) Ne = 0.;
  if ( NexONi < 0. ) NexONi = 0.;
  if ( L < 0. ) L = 0.;
  if ( L > 1. ) L = 1.; //Lindhard Factor
  if ( energy < 0.001*Wq_eV/L ) { Nph = 0.; Ne = 0.; }
  
  YieldResult result;
  result.PhotonYield = Nph;
  result.ElectronYield=Ne;
  result.ExcitonRatio =NexONi;
  result.Lindhard = L;
  result.ElectricField = dfield;
  return result; //everything needed to calculate fluctuations
  
}

void NESTcalc::SetRandomSeed ( unsigned long int s ) {
  
  rng.seed(s);
  
}

NESTcalc::NESTcalc ( ) {
  
  rng.seed(0);
  
}

vector<double> NESTcalc::GetS1 ( int Nph, double dx, double dy,
  double dz, double driftVelocity, double dV_mid, INTERACTION_TYPE type_num ) {
  
  vector<double> scintillation(9);  // return vector
  vector<double> newSpike(2); // for re-doing spike counting more precisely
  
  // Add some variability in g1 drawn from a polynomial spline fit
  double posDep = FitS1 ( dx, dy, dz );
  double dt = ( TopDrift - dz ) / driftVelocity;
  double dz_center = TopDrift - dV_mid * dtCntr; //go from t to z
  posDep /= FitS1 ( 0., 0., dz_center ); // XYZ always in mm now never cm
  
  // generate a number of PMT hits drawn from a binomial distribution. Initialize number of photo-electrons
  int nHits=BinomFluct(Nph,g1*posDep), Nphe = 0;
  
  // Initialize the pulse area and spike count variables
  double pulseArea = 0., spike = 0., prob;
  
  // If single photo-electron efficiency is under 1 and the threshold is above 0 (some phe will be below threshold)
  if ( sPEthr > 0. && nHits < numPMTs ) { // digital nHits eventually becomes spikes (spike++) based upon threshold
    // Step through the pmt hits
    for ( int i = 0; i < nHits; i++ ) {
      // generate photo electron, integer count and area
      double phe1 = rand_gauss(1.,sPEres) + rand_gauss(noise[0],noise[1]); Nphe++; if(phe1>DBL_MAX)phe1=1.;if(phe1<-DBL_MAX)phe1=0.;
      prob = rand_uniform();
      // zero the area if random draw determines it wouldn't have been observed.
      if ( prob > sPEeff ) { phe1 = 0.; } //add an else with Nphe++ if not doing mc truth
      // Generate a double photo electron if random draw allows it
      double phe2 = 0.;
      if ( rand_uniform() < P_dphe ) {
	// generate area and increment the photo-electron counter
	phe2 = rand_gauss(1.,sPEres) + rand_gauss(noise[0],noise[1]); Nphe++; if(phe2>DBL_MAX)phe2=1.;if(phe2<-DBL_MAX)phe2=0.;
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
  
  scintillation[0] = nHits; scintillation[1] = Nphe; //MC-true integer hits in same OR different PMTs, first without then with double phe's (Nphe > nHits)
  scintillation[2] = pulseArea; scintillation[3] = pulseAreaC; //floating real# smeared DAQ pulse areas in phe, uncorrected and XYZ corrected respectively
  scintillation[4] = Nphd; scintillation[5] = NphdC; //same as pulse areas except adjusted *downward* by constant for average 2-PE effect. (LUX phd units)
  scintillation[6] = spike; scintillation[7] = spikeC; //uncorrected and pos-corrected spike counts, both floats, but made more accurate later in GetSpike
  
  if ( spike < coinLevel ) //no chance of meeting coincidence requirement. Here, spike is still "perfect" (integer)
    prob = 0.;
  else {
    if ( spike > 10. ) prob = 1.;
    else {
      if ( coinLevel == 0 ) prob = 1.;
      else if ( coinLevel == 1 ) {
        if ( spike >= 1. ) prob = 1.;
        else prob = 0.;
      }
      else if ( coinLevel == 2 ) prob = 1.-pow((double)numPMTs,1.-spike);
      else {
        double numer = 0., denom = 0.;
        for ( int i = spike; i > 0; i-- ) {
          denom += nCr ( numPMTs, i );
          if ( i >= coinLevel ) numer += nCr ( numPMTs, i );
	}
        prob = numer / denom;
      } //end of case of coinLevel of 3 or higher
    } // end of case of spike is equal to 9 or lower
  } //the end of case of spike >= coinLevel
  
  newSpike = GetSpike ( Nph, dx, dy, dz, driftVelocity, dV_mid, scintillation ); // over-write spike with smeared version, ~real data
  scintillation[6] = newSpike[0]; // uncorr
  scintillation[7] = newSpike[1]; // 3-D corr spike RQ
  
  if ( rand_uniform() < prob ) // coincidence has to happen in different PMTs
    { ; }
  else { // some of these are set to -1 to flag them as having been below threshold
    scintillation[0] *= -1.; if ( scintillation[0] == 0. ) scintillation[0] = -DBL_MIN;
    scintillation[1] *= -1.; if ( scintillation[1] == 0. ) scintillation[1] = -DBL_MIN;
    scintillation[2] *= -1.; if ( scintillation[2] == 0. ) scintillation[2] = -DBL_MIN;
    scintillation[3] *= -1.; if ( scintillation[3] == 0. ) scintillation[3] = -DBL_MIN;
    scintillation[4] *= -1.; if ( scintillation[4] == 0. ) scintillation[4] = -DBL_MIN;
    scintillation[5] *= -1.; if ( scintillation[5] == 0. ) scintillation[5] = -DBL_MIN;
    scintillation[6] *= -1.; if ( scintillation[6] == 0. ) scintillation[6] = -DBL_MIN;
    scintillation[7] *= -1.; if ( scintillation[7] == 0. ) scintillation[7] = -DBL_MIN;
  }
  
  scintillation[8] =g1;
  
  return scintillation;
  
}

vector<double> NESTcalc::GetS2 ( int Ne, double dx, double dy, double dt,
				 double driftVelocity, bool IsInGasPhase ) {
  
  vector<double> ionization(9);
  double alpha = 0.137, beta = 177., gamma = 45.7;
  double epsilonRatio = 1.85 / 1.00126;
  if ( IsInGasPhase ) epsilonRatio = 1.;
  
  // Add some variability in g1_gas drawn from a polynomial spline fit
  double posDep = FitS2 ( dx, dy ); // XY is always in mm now, never cm
  posDep /= FitS2 ( 0., 0. );
  double dz = TopDrift - dt * driftVelocity;
  
  double E_liq = E_gas / epsilonRatio; //kV per cm
  double ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645; // arXiv:1710.11032
  if ( ExtEff > 1. || IsInGasPhase ) ExtEff = 1.;
  if ( ExtEff < 0. ) ExtEff = 0.;
  int Nee = BinomFluct(Ne,ExtEff*exp(-dt/eLife_us)); //MAKE this 1 for SINGLE e- DEBUGGING
  
  double elYield = (alpha*E_gas*1e3-beta*p_bar-gamma)* // arXiv:1207.2292
    ( anode - TopDrift ) * 0.1; //EL gap in mm -> cm, affecting S2 size linearly
  if ( (anode - TopDrift) <= 0. ) {
    cerr << "\tERR: The gas gap in the S2 calculation broke!!!!" << endl;
  }
  long Nph = long(floor(rand_gauss(elYield*double(Nee),
				   sqrt(s2Fano*elYield*double(Nee)))+0.5));
  long nHits = BinomFluct(Nph,g1_gas);
  long Nphe = nHits + BinomFluct(nHits,P_dphe);
  double pulseArea=rand_gauss(Nphe,sPEres*sqrt(Nphe));
  double pulseAreaC= pulseArea/exp(-dt/eLife_us);
  double Nphd = pulseArea / (1.+P_dphe);
  double NphdC= pulseAreaC/ (1.+P_dphe);
  
  double S2b = rand_gauss(S2botTotRatio*pulseArea,sqrt(S2botTotRatio*pulseArea*(1.-S2botTotRatio)));
  double S2bc= S2b / exp(-dt/eLife_us); // for detectors using S2 bottom-only in their analyses
  
  ionization[0] = Nee; ionization[1] = Nph; //integer number of electrons unabsorbed in liquid then getting extracted, followed by raw number of photons produced in the gas gap
  ionization[2] = nHits; ionization[3] = Nphe; //identical definitions to GetS1 follow, see above, except no spike, as S2 too big generally. S2 has more steps than S1 (e's 1st)
  if ( s2_thr >= 0 ) {
    ionization[4] = pulseArea; ionization[5] = pulseAreaC;
    ionization[6] = Nphd; ionization[7] = NphdC;
  }
  else { // the negative threshold is a polymorphic hidden feature: allows for header switching from total S2 to S2 bottom; doesn't mean literally negative, nor below threshold
    ionization[4] = S2b; ionization[5] = S2bc;
    ionization[6] = S2b / (1.+P_dphe); ionization[7] = S2bc / (1.+P_dphe);
  }
  
  if(pulseArea<fabs(s2_thr)) for(int i=0;i<8;i++) { ionization[i]*=-1.; if(ionization[i]==0.)ionization[i]=-DBL_MIN; }
  
  double SE = elYield* g1_gas;
  double g2 = ExtEff * SE;
  if ( s2_thr < 0 )
    g2 *= S2botTotRatio;
  ionization[8]=g2;
  
  if ( !dx && !dy && !dt && Ne == 1 ) {
    cout << endl << "g1 = " << g1 << " phd per photon\tg2 = " << g2 << " phd per electron (e-EE = ";
    cout << ExtEff*100. << "%, while SE_mean = " << SE << ")\t";
  }
  
  return ionization;
  
}

long double NESTcalc::Factorial ( double x ) {
  
  return tgammal ( x + 1. );
  
}

double NESTcalc::nCr ( double n, double r ) {
  
  return Factorial(n) /
    ( Factorial(r) * Factorial(n-r) );
  
}

DetectorParameters NESTcalc::GetDetector ( double xPos_mm, double yPos_mm, double zPos_mm,
					   bool IsInGasPhase ) {
  
  DetectorParameters detParam;
  vector<double> secondary(9);
  
  detParam.temperature = T_Kelvin;
  detParam.pressure = p_bar;
  detParam.GXeInterface = TopDrift;
  detParam.efFit = FitEF ( xPos_mm, yPos_mm, zPos_mm );
  detParam.rad = radius;
  detParam.dtExtrema[0] = dt_min;
  detParam.dtExtrema[1] = dt_max;
  
  if ( xPos_mm == 0. &&
       yPos_mm == 0. &&
       zPos_mm == detParam.GXeInterface / 2. ) {
    secondary = GetS2 ( 1, 0., 0., 0., 1., IsInGasPhase );
  }
  
  return detParam; //everything needed for testNEST to work
  
}

void NESTcalc::DriftRangeOverride ( double drift_low, double drift_high, DetectorParameters &detParam ) {
  
  // Grab previous drift time range in detector parameters
  double prev_dt_min = detParam.dtExtrema[0];
  double prev_dt_max = detParam.dtExtrema[1];
  
  // Reset drift time minimum and maximum
  detParam.dtExtrema[0] = drift_low;
  detParam.dtExtrema[1] = drift_high;
  
  // Ensure that we are not setting the new drift range outside the bounds of the previously set values.
  // This is the safest way to implement, so that we aren't working outside the bounds of the detector settings file.
  if (detParam.dtExtrema[0] < prev_dt_min || detParam.dtExtrema[1] > prev_dt_max || detParam.dtExtrema[0] > detParam.dtExtrema[1]) {
    cerr << "*** New drift time range completely outside of original fiducial! ***" << endl;
    exit(1);
  }
  
}

vector<double> NESTcalc::GetSpike ( int Nph, double dx, double dy, double dz,
  double driftSpeed, double dS_mid, vector<double> oldScint ) {
  
  vector<double> newSpike(2);
  
  if ( oldScint[7] > 70. ) {
    newSpike[0] = oldScint[6]; newSpike[1] = oldScint[7];
    return newSpike;
  }
  newSpike[0] = fabs(oldScint[6]);
  newSpike[0] = rand_gauss(newSpike[0],(sPEres/4.)*sqrt(newSpike[0]));
  if ( newSpike[0] < 0.0 ) newSpike[0] = 0.0;
  newSpike[1] = newSpike[0] / FitS1 ( dx, dy, dz ) * FitS1 ( 0., 0., TopDrift - dS_mid * dtCntr );
  
  return newSpike; // regular and position-corrected spike counts returned
  
}

int NESTcalc::SelectRanXeAtom ( double isotope ) {
  
  int A;
  if ( isotope > 0.000 && isotope <= 0.090 )
    A = 124;
  else if ( isotope > 0.090 && isotope <= 0.180 )
    A = 126;
  else if ( isotope > 0.180 && isotope <= 2.100 )
    A = 128;
  else if ( isotope > 2.100 && isotope <= 28.54 )
    A = 129;
  else if ( isotope > 28.54 && isotope <= 32.62 )
    A = 130;
  else if ( isotope > 32.62 && isotope <= 53.80 )
    A = 131;
  else if ( isotope > 53.80 && isotope <= 80.69 )
    A = 132;
  else if ( isotope > 80.69 && isotope <= 91.13 )
    A = 134;
  else
    A = 136;
  return A;

}

double NESTcalc::SetDensity ( double Kelvin, double bara ) { // currently only for fixed pressure (saturated vapor pressure); will add pressure dependence later
  
  if ( Kelvin < 161.40 ) { // solid Xenon
    cerr << "\nWARNING: SOLID PHASE. IS THAT WHAT YOU WANTED?\n";
    return 3.41; // from Yoo at 157K
    // other sources say 3.100 (Wikipedia, 'maximum') and 3.64g/mL at an unknown temperature
  }
  
  double VaporP_bar; //we will calculate using NIST
  if ( Kelvin < 289.7 ) VaporP_bar = pow(10.,4.0519-667.16/Kelvin);
  else VaporP_bar = DBL_MAX;
  if ( bara < VaporP_bar ) {
    double density = bara * 1e5 / ( Kelvin * 8.314 ); //ideal gas law approximation, mol/m^3
    density *= MOLAR_MASS * 1e-6;
    cerr << "\nWARNING: GAS PHASE. IS THAT WHAT YOU WANTED?\n"; return density; // in g/cm^3
  }
  
  return
    2.9970938084691329E+02 * exp ( -8.2598864714323525E-02 * Kelvin ) - 1.8801286589442915E+06 * exp ( - pow ( ( Kelvin - 4.0820251276172212E+02 ) / 2.7863170223154846E+01, 2. ) )
    - 5.4964506351743057E+03 * exp ( - pow ( ( Kelvin - 6.3688597345042672E+02 ) / 1.1225818853661815E+02, 2. ) )
    + 8.3450538370682614E+02 * exp ( - pow ( ( Kelvin + 4.8840568924597342E+01 ) / 7.3804147172071107E+03, 2. ) )
    - 8.3086310405942265E+02; // in grams per cubic centimeter based on zunzun fit to NIST data; will add gas later
  
}

double NESTcalc::SetDriftVelocity ( double Kelvin, double Density, double eField ) { //for liquid and solid only
  
  if ( Density < 1. ) return SetDriftVelocity_MagBoltz ( Density, eField );
  
  double speed = 0.0; // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
  int i, j; double vi, vf, slope, Ti, Tf, offset;
  
  double polyExp[11][7] = { { -3.1046, 27.037, -2.1668, 193.27, -4.8024, 646.04, 9.2471 }, //100K
                            { -2.7394, 22.760, -1.7775, 222.72, -5.0836, 724.98, 8.7189 }, //120
                            { -2.3646, 164.91, -1.6984, 21.473, -4.4752, 1202.2, 7.9744 }, //140
                            { -1.8097, 235.65, -1.7621, 36.855, -3.5925, 1356.2, 6.7865 }, //155
                            { -1.5000, 37.021, -1.1430, 6.4590, -4.0337, 855.43, 5.4238 }, //157, merging Miller with Yoo
                            { -1.4939, 47.879, 0.12608, 8.9095, -1.3480, 1310.9, 2.7598 }, //163, merging Miller with Yoo
                            { -1.5389, 26.602, -.44589, 196.08, -1.1516, 1810.8, 2.8912 }, //165
                            { -1.5000, 28.510, -.21948, 183.49, -1.4320, 1652.9, 2.884 }, //167
                            { -1.1781, 49.072, -1.3008, 3438.4, -.14817, 312.12, 2.8049 }, //184
                            {  1.2466, 85.975, -.88005, 918.57, -3.0085, 27.568, 2.3823 }, //200
                            { 334.60 , 37.556, 0.92211, 345.27, -338.00, 37.346, 1.9834 } }; //230
  
  double Temperatures[11] = { 100., 120., 140., 155., 157., 163., 165., 167., 184., 200., 230. };
  
  if ( Kelvin >= Temperatures[0] && Kelvin < Temperatures[1] ) i = 0;
  else if ( Kelvin >= Temperatures[1] && Kelvin < Temperatures[2] ) i = 1;
  else if ( Kelvin >= Temperatures[2] && Kelvin < Temperatures[3] ) i = 2;
  else if ( Kelvin >= Temperatures[3] && Kelvin < Temperatures[4] ) i = 3;
  else if ( Kelvin >= Temperatures[4] && Kelvin < Temperatures[5] ) i = 4;
  else if ( Kelvin >= Temperatures[5] && Kelvin < Temperatures[6] ) i = 5;
  else if ( Kelvin >= Temperatures[6] && Kelvin < Temperatures[7] ) i = 6;
  else if ( Kelvin >= Temperatures[7] && Kelvin < Temperatures[8] ) i = 7;
  else if ( Kelvin >= Temperatures[8] && Kelvin < Temperatures[9] ) i = 8;
  else if ( Kelvin >= Temperatures[9] && Kelvin <= Temperatures[10] ) i = 9;
  else {
    cerr << "\nERROR: TEMPERATURE OUT OF RANGE (100-230 K)\n";
  }
  
  j = i + 1;
  Ti = Temperatures[i];
  Tf = Temperatures[j];
  // functional form from http://zunzun.com
  vi = polyExp[i][0]*exp(-eField/polyExp[i][1])+polyExp[i][2]*exp(-eField/polyExp[i][3])+polyExp[i][4]*exp(-eField/polyExp[i][5])+polyExp[i][6];
  vf = polyExp[j][0]*exp(-eField/polyExp[j][1])+polyExp[j][2]*exp(-eField/polyExp[j][3])+polyExp[j][4]*exp(-eField/polyExp[j][5])+polyExp[j][6];
  if ( Kelvin == Ti ) return vi;
  if ( Kelvin == Tf ) return vf;
  if ( vf < vi ) {
    offset = (sqrt((Tf*(vf-vi)-Ti*(vf-vi)-4.)*(vf-vi))+sqrt(Tf-Ti)*(vf+vi))/(2.*sqrt(Tf-Ti));
    slope = -(sqrt(Tf-Ti)*sqrt((Tf*(vf-vi)-Ti*(vf-vi)-4.)*(vf-vi))-(Tf+Ti)*(vf-vi))/(2.*(vf-vi));
    speed = 1. / ( Kelvin - slope ) + offset;
  }
  else {
    slope = ( vf - vi ) / ( Tf - Ti );
    speed = slope * ( Kelvin - Ti ) + vi;
  }
  
  if ( speed <= 0. ) { cerr << "\nERROR: DRIFT SPEED NON-POSITIVE -- FIELD TOO LOW\n"; }
  return speed;
  
}

double NESTcalc::SetDriftVelocity_MagBoltz ( double density, double efieldinput ) //Nichole Barry UCD 2011
{
  density *= NEST_AVO / MOLAR_MASS;
  //Gas equation one coefficients (E/N of 1.2E-19 to 3.5E-19)
  double gas1a = 395.50266631436,
    gas1b = -357384143.004642, gas1c = 0.518110447340587;
  //Gas equation two coefficients (E/N of 3.5E-19 to 3.8E-17)
  double gas2a = -592981.611357632, gas2b = -90261.9643716643,
    gas2c = -4911.83213989609, gas2d = -115.157545835228,
    gas2f = -0.990440443390298, gas2g = 1008.30998933704, gas2h = 223.711221224885;
  double edrift = 0., gasdep = efieldinput / density, gas1fix = 0., gas2fix = 0.;
  
  if ( gasdep < 1.2e-19 && gasdep >= 0. ) edrift = 4e22 * gasdep;
  if ( gasdep < 3.5e-19 && gasdep >= 1.2e-19 ) {
    gas1fix = gas1b * pow ( gasdep, gas1c ); edrift = gas1a * pow ( gasdep, gas1fix );
  }
  if ( gasdep < 3.8e-17 && gasdep >= 3.5e-19 ) {
    gas2fix = log ( gas2g * gasdep );
    edrift = ( gas2a + gas2b * gas2fix + gas2c * pow ( gas2fix, 2. ) + gas2d * pow ( gas2fix, 3. )
               + gas2f * pow ( gas2fix, 4. ) ) * ( gas2h * exp ( gasdep ) );
  }
  if ( gasdep >= 3.8e-17 ) edrift = 6e21 * gasdep - 32279.;
  
  return edrift * 1e-5; // from cm/s into mm per microsecond
}

vector<double> NESTcalc::SetDriftVelocity_NonUniform ( double rho, bool IsInGasPhase,
  double z_step ) {
  
  vector<double> speedTable;
  DetectorParameters detParam;
  double driftTime, zz;
  
  for ( double pos_z = 0.0; pos_z < TopDrift; pos_z += z_step ) {
    
    driftTime = 0.0;
    for ( zz = pos_z; zz < TopDrift; zz += z_step ) {
      
      detParam = GetDetector ( 0., 0., zz, IsInGasPhase );
      if ( pos_z > gate ) {
	if ( !IsInGasPhase )
	  driftTime += z_step/SetDriftVelocity(T_Kelvin,rho,E_gas/(1.85/1.00126)*1e3);
	else // if gate == TopDrift properly set, shouldn't happen
	  driftTime += z_step/SetDriftVelocity(T_Kelvin,rho,E_gas*1e3);
      }
      else
	driftTime += z_step/SetDriftVelocity(T_Kelvin,rho,detParam.efFit);
      
    }
    
    speedTable.push_back ( ( zz - pos_z ) / driftTime ); //uses highest zz
    
  }
  
  return speedTable;
}

vector<double> NESTcalc::xyResolution ( double xPos_mm, double yPos_mm, double A_top ) {
  
  vector<double> xySmeared(2);
  A_top *= (1.-S2botTotRatio);
  
  double radius = sqrt(pow(xPos_mm,2.)+pow(yPos_mm,2.));
  double kappa = PosResBase + exp ( PosResExp * radius ); // arXiv:1710.02752
  double sigmaR = kappa / sqrt ( A_top ); // ibid.
  
  double phi = 2. * M_PI * rand_uniform();
  sigmaR = rand_gauss ( 0.0, sigmaR );
  double sigmaX = sigmaR * cos ( phi );
  double sigmaY = sigmaR * sin ( phi );
  
  xySmeared[0] = xPos_mm + sigmaX;
  xySmeared[1] = yPos_mm + sigmaY;
  
  return xySmeared; //new X and Y position in mm with empirical smearing. LUX Run03 example

}

void NESTcalc::SetDetector ( string paramName, double paramValue, double evtTime ) {
  
  if ( paramName == "g1" ) g1 = paramValue; //add time dependence as desired using evtTime, instead of fixed value
  else if ( paramName == "sPEres" ) sPEres = paramValue;
  else if ( paramName == "sPEthr" ) sPEthr = paramValue;
  else if ( paramName == "sPEeff" ) sPEeff = paramValue;
  else if ( paramName == "noise[0]" ) noise[0] = paramValue;
  else if ( paramName == "noise[1]" ) noise[1] = paramValue;
  else if ( paramName == "P_dphe" ) P_dphe = paramValue;
  else if ( paramName == "coinLevel" ) coinLevel = (int)paramValue;
  else if ( paramName == "numPMTs" ) numPMTs = (int)paramValue;
  else if ( paramName == "g1_gas" ) g1_gas = paramValue;
  else if ( paramName == "s2Fano" ) s2Fano = paramValue;
  else if ( paramName == "s2_thr" ) s2_thr = paramValue;
  else if ( paramName == "S2botTotRatio" ) S2botTotRatio = paramValue;
  else if ( paramName == "E_gas" ) E_gas = paramValue;
  else if ( paramName == "eLife_us" ) eLife_us = paramValue;
  else if ( paramName == "T_Kelvin" ) T_Kelvin = paramValue;
  else if ( paramName == "p_bar" ) p_bar = paramValue;
  else if ( paramName == "dtCntr" ) dtCntr = paramValue;
  else if ( paramName == "dt_min" ) dt_min = paramValue;
  else if ( paramName == "dt_max" ) dt_max = paramValue;
  else if ( paramName == "radius" ) radius = paramValue;
  else if ( paramName == "TopDrift" ) TopDrift = paramValue;
  else if ( paramName == "anode" ) anode = paramValue;
  else if ( paramName == "gate" ) gate = paramValue;
  else if ( paramName == "PosResExp" ) PosResExp = paramValue;
  else if ( paramName == "PosResBase" ) PosResBase = paramValue;
  else {
    cerr << "CAUTION: Invalid std::string parameter name. Using default detector effect parameters";
    cerr << "Check the detector.hh header file for valid variable names.\n";
  }

}
