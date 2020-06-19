
#include <random>
#include <iostream>
#include <cfloat>
#include <vector>

unsigned int NUMEVT;
#define S2BORD 1000 //ns
#define RISEDEF 0.05 //5% rise time used as the t0 point
#define UNC 0.35 //uncertainty in definition of rise t

using namespace std;

int main ( int argc, char** argv ) { //compile with g++ -Ofast PulseShape.cpp -o PulseShape.exe
  
  if ( argc < 3 ) {
    cerr << "2 args needed: S1 (non-c) in phe (NOT phd). 1st value is min, 2nd max." << endl;
    return EXIT_FAILURE;
  }
  
  double MaxPossibleTime = DBL_MAX, S1LimPE[2]={atof(argv[1]),atof(argv[2])}; //avoid Seg Fault!
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.,20.); //ns for ran offset. 2 bin default
  vector<long> num, win; char line[60];
  vector<double> tns, bot, top;
  FILE* ifp = fopen ( "photon_times.txt", "r" );
  double a, b, c, d, e; long i = 0;
  
  while ( 1 ) {
    if ( i == 0 ) {
      while ( getc ( ifp ) != '\n' ) { ; }
      i = 1;
      continue;
    }
    if ( fgets ( line, 40, ifp ) == NULL ) {
      NUMEVT = unsigned(a) + 1; break;
    }
    else ;
    sscanf ( line, "%lf\t%lf\t%lf\t%lf\t%lf", &a, &b, &c, &d, &e );
    if ( feof ( ifp ) )
      break;
    num.push_back(long(a));
    tns.push_back(b);
    if ( b > S2BORD )
      e = 0;
    bot.push_back(c);
    top.push_back(d);
    win.push_back(long(e));
    //cerr << a << " " << b << " " << c << " " << d << " " << e << endl;
  }
  fclose ( ifp );
  
  vector<double> S1tot(NUMEVT);
  vector<double> S2tot(NUMEVT);
  vector<double> S1f30(NUMEVT);
  vector<double> T0X(NUMEVT);
  
  for ( i = 0; i < num.size(); i++ ) {
    T0X[num[i]] = MaxPossibleTime;
    if ( tns[i] < S2BORD ) S1tot[num[i]] += top[i] + bot[i];
    else S2tot[num[i]] += top[i] + bot[i];
  }
  double fraction, soFar = 0.; normal_distribution<double> distribution2(RISEDEF,RISEDEF*UNC);
  for ( i = 0; i < num.size(); i++ ) {
    if ( tns[i] < S2BORD && T0X[num[i]] >= MaxPossibleTime ) {
      soFar += top[i] + bot[i];
      fraction = distribution2(generator); if ( fraction < 0. ) fraction = 0.;
      if ( soFar > (fraction*S1tot[num[i]]) && T0X[num[i]] >= MaxPossibleTime )
	{ T0X[num[i]] = tns[i]; soFar = 0.0; continue; }
    }
  }
  
  double area[2] = {0.0,0.0}; double range[4] = {-14.,134.,-8.,32.}; double rRange[4]={range[0],range[1],range[2],range[3]};
  for ( i = 0; i < num.size(); i++ ) {
    if ( T0X[num[i]] < MaxPossibleTime ) tns[i] -= T0X[num[i]];
    if ( tns[i] < (S2BORD-T0X[num[i]]) && win[i] == 1 ) {
      if ( tns[i] > rRange[0] && tns[i] < rRange[1] ) {
	area[0] += top[i] + bot[i];
	if ( tns[i] > rRange[2] && tns[i] < rRange[3] )
	  area[1] += top[i] + bot[i];
      }
    }
    else {
      double offset = distribution(generator); //consider unique offset for each range[]
      rRange[0] = range[0]+offset; rRange[1] = range[1]+offset; rRange[2] = range[2]+offset; rRange[3] = range[3]+offset;
      if ( S1f30[num[i]] == 0. && area[0] != 0. && area[1] != 0. && area[0] != area[1] ) S1f30[num[i]] = area[1]/area[0];
      area[0] = 0.0; area[1] = 0.0;
    }
  }
  
  double mean = 0.0; int j = 0; vector<double> NonBlank;
  cout << "evt#\tS1, pe\tpromptFrac\n";
  for ( i = 0; i < NUMEVT; i++ ) {
    //cerr << T0X[i] << endl;
    if ( S1tot[i] > S1LimPE[0] && S2tot[i] > 0. && S1f30[i] > 0. && S1f30[i] < 1. && S1tot[i] < S1LimPE[1] )
      { mean += S1f30[i]; j++; cout << i << "\t" << S1tot[i] << "\t" << S1f30[i] << endl; NonBlank.push_back(S1f30[i]); }
    else
      cout << i << endl;
  }
  if ( j == 0 ) { cerr << "no events in range!\n"; return EXIT_FAILURE; }
  mean /= double ( j );
  std::sort ( NonBlank.begin(), NonBlank.end() );
  double median = NonBlank[j/2];
  
  double sigma = 0.0;
  for ( i = 0; i < NUMEVT; i++ )
    if ( S1f30[i] != 0. && S1tot[i] > S1LimPE[0] && S1tot[i] < S1LimPE[1] ) sigma += ( mean - S1f30[i] ) * ( mean - S1f30[i] );
  sigma /= (double(j)-1.); sigma = sqrt(sigma); double skew = 0.0;
  for ( i = 0; i < NUMEVT; i++ )
    if ( S1f30[i] != 0. && S1tot[i] > S1LimPE[0] && S1tot[i] < S1LimPE[1] ) skew += pow ( ( S1f30[i] - mean ) / sigma, 3. );
  skew /= double(j);
  cerr << "MEDIAN\t\tSIGMA\t\tSKEWNESS\n";
  cerr << median << "\t" << sigma << "\t" << skew << endl;
  
  return EXIT_SUCCESS;

}
