
#include <random>
#include <iostream>
#include <cfloat>
#include <vector>

#define NUMEVT 10000
#define S2BORD 1000 //ns

double S1tot[NUMEVT];
double S2tot[NUMEVT];
double S1f30[NUMEVT];
double T05[NUMEVT];
using namespace std;

int main ( int argc, char** argv ) { //compile with g++ -Ofast PulseShape.cpp -o PulseShape.exe
  
  double MaxPossibleTime = DBL_MAX, S1LimPE[2]={atof(argv[1]),atof(argv[2])}; //avoid Seg Fault!
  default_random_engine generator;
  uniform_real_distribution<double> distribution(0.,15.); //ns for random offsets
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
    if ( fgets ( line, 40, ifp ) == NULL ) break;
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
  
  for ( i = 0; i < num.size(); i++ ) {
    T05[num[i]] = MaxPossibleTime;
    if ( tns[i] < S2BORD ) S1tot[num[i]] += top[i] + bot[i];
    else S2tot[num[i]] += top[i] + bot[i];
  }
  double soFar = 0.;
  for ( i = 0; i < num.size(); i++ ) {
    if ( tns[i] < S2BORD && T05[num[i]] >= MaxPossibleTime ) {
      soFar += top[i] + bot[i];
      if ( soFar > (0.05*S1tot[num[i]]) && T05[num[i]] >= MaxPossibleTime )
	{ T05[num[i]] = tns[i]; soFar = 0.0; continue; }
    }
  }
  
  double area[2] = {0.0,0.0}; double range[4] = {-14.,134.,-8.,32.}; double rRange[4]={range[0],range[1],range[2],range[3]};
  for ( i = 0; i < num.size(); i++ ) {
    //if ( T05[num[i]] < MaxPossibleTime ) tns[i] -= T05[num[i]];
    if ( tns[i] < S2BORD && win[i] == 1 ) {
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
  
  double mean = 0.0; int j = 0;
  cout << "evt#\tS1, pe\tpromptFrac\n";
  for ( i = 0; i < NUMEVT; i++ ) {
    //cerr << T05[i] << endl;
    if ( S1tot[i] > S1LimPE[0] && S2tot[i] > 0. && S1f30[i] > 0. && S1f30[i] < 1. && S1tot[i] < S1LimPE[1] )
      { mean += S1f30[i]; j++; cout << i << "\t" << S1tot[i] << "\t" << S1f30[i] << endl; }
    else
      cout << i << endl;
  } mean /= double ( j );
  
  double sigma = 0.0;
  for ( i = 0; i < NUMEVT; i++ )
    if ( S1f30[i] != 0. && S1tot[i] > S1LimPE[0] && S1tot[i] < S1LimPE[1] ) sigma += ( mean - S1f30[i] ) * ( mean - S1f30[i] );
  sigma /= (double(j)-1.); sigma = sqrt(sigma); double skew = 0.0;
  for ( i = 0; i < NUMEVT; i++ )
    if ( S1f30[i] != 0. && S1tot[i] > S1LimPE[0] && S1tot[i] < S1LimPE[1] ) skew += pow ( ( S1f30[i] - mean ) / sigma, 3. );
  skew /= double(j); cerr << "MEAN\t\tSIGMA\t\tSKEW\n";
  cerr << mean << "\t" << sigma << "\t" << skew << endl;
  
  return 0;
}
