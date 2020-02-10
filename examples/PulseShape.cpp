
#include <iostream>
#include <vector>

#define NUMEVT 10000
#define S2BORD 1000 //ns

double S1tot[NUMEVT];
double S2tot[NUMEVT];
double S1f30[NUMEVT];
using namespace std;

int main ( int argc, char** argv ) {
  
  vector<long> num, win; char line[60];
  vector<double> tns, bot, top;
  FILE* ifp = fopen ( "photon_times.txt", "r" );
  double a,b,c,d,e; long i = 0;
  
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
    if ( tns[i] < S2BORD )
      S1tot[num[i]] += top[i] + bot[i];
    else
      S2tot[num[i]] += top[i] + bot[i];
  }
  
  double area[2] = {0.0,0.0}; double range[4] = {-14.,134.,-8.,32.};
  for ( i = 0; i < num.size(); i++ ) {
    if ( tns[i] < S2BORD && win[i] == 1 ) {
      if ( tns[i] > range[0] && tns[i] < range[1] ) {
	area[0] += top[i] + bot[i];
	if ( tns[i] > range[2] && tns[i] < range[3] )
	  area[1] += top[i] + bot[i];
      }
    }
    else {
      if ( S1f30[num[i]] == 0. &&
	   area[0] != 0. && area[1] != 0. && area[0] != area[1] ) S1f30[num[i]] = area[1]/area[0];
      area[0] = 0.0; area[1] = 0.0;
    }
  }
  
  for ( i = 0; i < NUMEVT; i++ ) {
    if ( S1tot[i] > 0. && S2tot[i] > 0. && S1f30[i] > 0. && S1f30[i] < 1. )
      cout << i << "\t" << S1tot[i] << "\t" << S1f30[i] << endl;
    else
      cout << i << endl;
  }
  
  return 0;

}
