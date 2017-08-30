
//#define DISCRIM
//#define FIT

#include "LUX.hh"

double rand_uniform ( ) {
  
  return ( double ) rand ( ) / ( double ) RAND_MAX;
  
}

double rand_gauss ( double mean, double sigma ) {
  
  double u = rand_uniform(), v = rand_uniform();
  return mean + sigma * sqrt(-2.*log(u))*cos(2.*M_PI*v);
  
}

double VonNeumann ( double xMin, double xMax, double yMin, double yMax );

int BinomFluct ( int N0, double prob ) {
  
  double mean = N0*prob;
  double sigma = sqrt(N0*prob*(1.-prob));
  int N1 = 0;
  
  if ( prob <= 0.00 ) return N1; if ( prob >= 1.00 ) return N0;
  
  if ( N0 < 10 ) {
    for(int i = 0; i < N0; i++) {
      if(rand_uniform() < prob) N1++;
    }
  }
  else {
    N1 =
      int(floor(rand_gauss(mean,sigma)+0.5));
  }
  
  if ( N1 > N0 ) N1 = N0; if ( N1 < 0 ) N1 = 0;
  
  return N1;
  
}

using namespace std;

typedef enum {
  
  //nuclear recoil
  NR = 0,
  WIMP = 1,
  B8 = 2,
  DD = 3,
  AmBe = 4,
  Cf = 5,
  ion = 6, //includes alphas, Pb-206
  
  //electron recoil
  gammaRay = 7,
  beta = 8,
  CH3T = 9,
  Kr83m = 10,
  
} INTERACTION_TYPES;

vector<double> GetQuanta ( int species, double energy, double density, double dfield ) {
  
  vector<double> yields(2);
  double massNum = 4., m2=78.324,m3=2.,m4=2.,m6=0.,m8=2., deltaT_ns = 400.;
  double Ly, Qy, totQ, Ne, Nph, Ni, recombProb, ThomasImel, densCorr;
  
  double Wq_eV = 1.9896+(20.8-1.9896)/(1.+pow(density/4.0434,1.4407));
  double alpha = 0.067366+density*0.039693;
  
  if ( species <= 5 ) {
    double epsilon = 11.5*energy*pow(54.,(-7./3.));
    densCorr = 13.7 / Wq_eV;
    double peter = 0.;
    double power = 0.28884*pow(dfield,-0.045639);
    Qy = 8.494/pow(1.+pow(energy/5.1215,1.671),power);
    totQ = 12.256*pow(energy,1.0770) - peter/epsilon;
    Nph = densCorr * (totQ - energy*Qy);
    Ly = Nph / energy;
    Qy = Qy + (totQ/energy - (Ly+Qy));
    Ne = Qy * energy;
    double a = 0.0, b = 1.0;
    Nph *= 1. / (1. + a*pow(epsilon,b));
  }
  
  else if ( species == 6 ) {
    double L = 0.96446/(1.+pow(massNum*massNum/19227.,0.99199));
    if ( massNum == 4 ) L = 0.56136*pow(energy,0.056972);
    double ThomasImel = 0.0067/pow(1.+pow(dfield/95.768,8.5673),0.060318)*pow(density/2.857,0.3);
    totQ = 1e3*L*energy/Wq_eV;
    Ni = totQ / (1.+alpha);
    recombProb = 1.-log(1.+(ThomasImel/4.)*Ni)/((ThomasImel/4.)*Ni);
    Nph = totQ * alpha / (1.+alpha) + recombProb*Ni; Ne = totQ - Nph;
  }
  
  else if ( species == 7 ) {
    double m1 = 35.028+(5.7254-35.028)/(1.+pow(dfield/78.904,.60422));
    double m5 = 21.416+(10.737-21.416)/(1.+pow(dfield/220.71,1.7433));
    double m7 = 66.825+(829.25-66.825)/(1.+pow(dfield/43.608,.83344));
    densCorr = 0.32856+0.23187*density;
    totQ = energy * 1000. / Wq_eV;
    Qy = m1 + (m2-m1) / (1.+pow(energy/m3,m4)) + m5 + (m6-m5) / (1.+pow(energy/m7,m8));
    Qy /= 73.*Wq_eV*1e-3;
    Ly = (73.-Qy)*densCorr/(73.*Wq_eV*0.001);
    Qy = Qy + (totQ/energy - (Ly+Qy));
    Ne = Qy * energy; Nph = Ly * energy;
  }
  
  else if ( species == 10 ) {
    totQ = energy * 1000. / Wq_eV;
    if ( energy == 9.4 ) {
      double m1 = 99678.-21574.*log10(dfield);
      m2 = 47.364+(131.69-47.364)/pow(1.+pow(dfield/71.368,2.4130),0.060318);
      Nph = (m1*pow(2.*deltaT_ns+10.,-1.5)+m2)/0.21;
    }
    else
      Nph = energy*(0.51987+(1.0036-0.51987)/(1.+pow(dfield/309.98,1.0844)))*65.5;
    Ne = totQ - Nph;
  }
  
  else {
    double m1 = 37.609+(9.4398-37.609)/(1.+pow(dfield/95.192,.65711));
    double m5 = 24.159;
    double m7 = 27.663+(694.46-27.663)/(1.+pow(dfield/28.595,.84217));
    densCorr = 0.32856+0.23187*density;
    totQ = energy * 1000. / Wq_eV;
    Qy = m1 + (m2-m1) / (1.+pow(energy/m3,m4)) + m5 + (m6-m5) / (1.+pow(energy/m7,m8));
    Qy /= 73.*Wq_eV*1e-3;
    Ly = (73.-Qy)*densCorr/(73.*Wq_eV*0.001);
    Qy = Qy + (totQ/energy - (Ly+Qy));
    Ne = Qy * energy; Nph = Ly * energy;
  }
  
  if ( Ne > m2*energy) Ne = m2*energy;
  if ( Nph< 0. ) Nph= 0.;
  if ( Ne < 0. ) Ne = 0.;
  
  yields.insert(yields.begin() + 0, Nph);
  yields.insert(yields.begin() + 1, Ne );
  
  return yields;
  
}

int main ( int argc, char** argv ) {

  if (argc != 5) {
    cout << "Enter particle type, recoil energy [keV], "
      << "Xe density [g/cm^3] E-field (V/cm)." << endl << "Particle type "
      << "can be NR, WIMP, B8, DD, AmBe, Cf, ion, gamma, beta, CH3T, Kr83m."
      << endl;
    return 2;
  }  
  vector<double> quanta(2);

  string type = argv[1]; int type_num;
  if ( type == "NR" ) type_num = NR;
  else if ( type == "WIMP")type_num = WIMP;
  else if ( type == "B8" ) type_num = B8;
  else if ( type == "DD" ) type_num = DD;
  else if ( type == "AmBe")type_num = AmBe;
  else if ( type == "Cf" ) type_num = Cf;
  else if ( type == "ion") type_num = ion;
  else if ( type =="gamma")type_num = gammaRay;
  else if ( type =="beta" )type_num = beta;
  else if ( type =="CH3T" )type_num = CH3T;
  else if (type == "Kr83m")type_num = Kr83m;
  else {cout << "Invalid particle type.\n"; return 3;};
  
  double keV = atof(argv[2]);
  double rho = atof(argv[3]);
  double fld = atof(argv[4]);
  
  quanta = GetQuanta(type_num,keV,rho,fld);
  cout << quanta[0] << "\t" << quanta[1] << endl;
  
  return 1;

}
