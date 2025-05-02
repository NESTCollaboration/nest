
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <random>
#include <vector>
#include "analysis_G2.hh"

unsigned long int NUMEVT;  // first global var: #lines from file
#define S2BORD 1000   // ns - border for S1/S2 discrimination
#define RISEDEF \
  0.05  // 5% rise time used as the t0 point. 10% for seminal Kwong paper
#define UNC 0.35  // uncertainty in definition of rise time
// 35% for Run04, but all-natural 0% better for Run03 and Xed
#define LOSTAREA \
  1.015  // O(1%) adjustment for untraceable tiny bug in S2 areas. But not
         // universal?

using namespace std;

int main ( int argc, char** argv ) {
  
  double S1LimPE[2] = { minS1, maxS1 };
  
  if ( verbosity < 2 ) {
    if ( argc < 3 ) {
      cerr << "2 args needed: S1 (non-c) in phe (NOT phd). 1st value is min, 2nd "
	      "max."
	   << endl;
      return EXIT_FAILURE;
    }
    else {
      S1LimPE[0] = atof(argv[1]); S1LimPE[1] = atof(argv[2]);  // avoid Seg Fault!
    }
  }
  
  double MaxPossibleTime = DBL_MAX; default_random_engine generator;
  uniform_real_distribution<double> distribution(
      0., 20.);  // ns for ran offset. 2 bin default. For Case: -7,0
  // for LUX Run03: better results with -4 to 16. Still 2 samples-wide, but
  // different centering
  vector<uint64_t> num, win;
  char line[60];
  vector<double> tns, bot, top;
  FILE* ifp = fopen("photon_times.txt", "r");
  double a, b, c, d, e;
  uint64_t i = 0;

  while (1) {
    if (i == 0) {
      while (getc(ifp) != '\n') {
        ;
      }
      i = 1;
      continue;
    }
    if (fgets(line, 40, ifp) == NULL) {
      NUMEVT = unsigned(a) + 1;
      if ( NUMEVT % 10 )
        cerr << "Warning! Check to see if last event(s) missing b/c S1=0.\n";
      break;
    } else ;
    if ( verbosity < 2 )
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf", &a, &b, &c, &d, &e);
    else
      sscanf(line, "%lf\t%lf", &a, &b);
    if (feof(ifp)) break;
    num.push_back(uint64_t(a));
    tns.push_back(b);
    if ( verbosity < 2 ) {
      if (b > S2BORD) e = 0;
      bot.push_back(c);
      top.push_back(d);
      win.push_back(uint64_t(e));
    }
    // cerr << a << " " << b << " " << c << " " << d << " " << e << endl;
  }
  fclose(ifp);
  
  if ( verbosity >= 2 ) { //for pseudo-S1-only spike mode
    double range[4] = {-200., 200., -69.7, 196.2};
    vector<double> means(NUMEVT, 0.);
    vector<double> count(NUMEVT, 0.);
    vector<double> S1total(NUMEVT, 0.);
    vector<double> S1inner(NUMEVT, 0.);
    for ( i = 0; i < num.size(); ++i ) {
      if ( tns[i] < 1e4 ) {
	++count[num[i]];
	means[num[i]] += tns[i];
      }
    }
    for ( i = 0; i < NUMEVT; ++i ) means[i] /= count[i];
    for ( i = 0; i < num.size(); ++i ) {
      tns[i] -= means[num[i]]; //correction on means
      if ( tns[i] >= range[0] && tns[i] <= range[1] )
	++S1total[num[i]];
      if ( tns[i] >= range[2] && tns[i] <= range[3] )
	++S1inner[num[i]];
    }
    cout << "evt#\tpromptFrac" << endl;
    for ( i = 0; i < NUMEVT; ++i ) {
      double promptFrac = 0.0;
      if ( S1total[i] > 0. )
	promptFrac = S1inner[i] / S1total[i];
      cout << i << "\t" << promptFrac << endl;
    }
    return EXIT_SUCCESS;
  }
  
  vector<double> S1tot(NUMEVT, 0.);
  vector<double> S2tot(NUMEVT, 0.);
  vector<double> S2max(NUMEVT, 0.);
  vector<double> CI_left(NUMEVT, 0.);
  vector<double> CI_right(NUMEVT, 0.);
  vector<double> S2width(NUMEVT, 0.);
  vector<double> S1f30(NUMEVT, 0.);
  vector<double> T0X(NUMEVT, 0.);
  vector<double> driftT(NUMEVT, 0.);

  for (i = 0; i < num.size(); ++i) {
    T0X[num[i]] = MaxPossibleTime;
    if (tns[i] < S2BORD)
      S1tot[num[i]] += top[i] + bot[i];
    else {
      S2tot[num[i]] += (top[i] + bot[i]) * LOSTAREA;
      if ((top[i] + bot[i]) > S2max[num[i]]) {
        S2max[num[i]] = top[i] + bot[i];
        driftT[num[i]] =
            tns[i] -
            S2BORD;  // minus 1 microsecond fudge factor to get mean drift
      } else
        ;
    }
  }
  double fraction = 0.0225,
         soFar = 0.;  // default Confidence Interval is 95.5% (~2-sigma)
  for (i = 0; i < num.size(); ++i) {
    if (tns[i] >= S2BORD && S2width[num[i]] == 0.) {
      soFar += (top[i] + bot[i]) * LOSTAREA;
      if (soFar > fraction * S2tot[num[i]] && CI_left[num[i]] == 0. &&
          CI_right[num[i]] == 0.)
        CI_left[num[i]] = tns[i];
      if (soFar > 0.500000 * S2tot[num[i]] && CI_right[num[i]] == 0. &&
          CI_left[num[i]] != 0.)
        CI_right[num[i]] = tns[i];  // left half only because of asymmetry
      if (CI_left[num[i]] != 0. && CI_right[num[i]] != 0.) {
        S2width[num[i]] = (CI_right[num[i]] - CI_left[num[i]]) / 2.;
        soFar = 0.0;
        continue;
      }  // div by 2 as 1s--1s=2sigma
    }
  }
  soFar = 0.00;
  normal_distribution<double> distribution2(RISEDEF, RISEDEF * UNC);
  for (i = 0; i < num.size(); ++i) {
    if (tns[i] < S2BORD && T0X[num[i]] >= MaxPossibleTime) {
      soFar += top[i] + bot[i];
      fraction = distribution2(generator);
      if (fraction < 0.) fraction = 0.;
      if (soFar > (fraction * S1tot[num[i]]) &&
          T0X[num[i]] >= MaxPossibleTime) {
        T0X[num[i]] = tns[i]; //time of "actual" S1 peak start
        soFar = 0.0;
        continue;
      }
    }
  }

  double area[2] = {0.0, 0.0};
  double range[4] = {-14., 134., -8., 32.};
  double rRange[4] = {range[0], range[1], range[2],
                      range[3]};  // Dahl: -50,300,-50,10
  for (i = 0; i < num.size(); ++i) {
    if (T0X[num[i]] < MaxPossibleTime) tns[i] -= T0X[num[i]];
    if (tns[i] < (S2BORD - T0X[num[i]]) && win[i] == 1) {
      if (tns[i] > rRange[0] && tns[i] < rRange[1]) {
        area[0] += top[i] + bot[i];
        if (tns[i] > rRange[2] && tns[i] < rRange[3])
          area[1] += top[i] + bot[i];
      }
    } else {
      double offset0 =
          0.0 +
          distribution(generator);  // consider unique offset for each
                                    // range[]. In original try makes bumps
      // double offset1 = 10.+distribution(generator);
      // double offset2 = 10.+distribution(generator);
      // double offset3 = 10.+distribution(generator);
      // randomized offset was to create extra smearing to better match real data with unknowns in the DAQ
      rRange[0] = range[0] + offset0;
      rRange[1] = range[1] + offset0;
      rRange[2] = range[2] + offset0;
      rRange[3] = range[3] + offset0;
      if (S1f30[num[i]] == 0. && area[0] != 0. && area[1] != 0. &&
          area[0] != area[1])
        S1f30[num[i]] = area[1] / area[0];
      area[0] = 0.0;
      area[1] = 0.0;
    }
  }

  double mean = 0.0;
  int j = 0;
  vector<double> NonBlank;
  cout << "evt#\tS1, pe\tpromptFrac\tdrift [ms]\tS2, pe\tWidth [ns]\n";
  for (i = 0; i < NUMEVT; ++i) {
    // cerr << T0X[i] << endl;
    if (S1tot[i] > S1LimPE[0] && S2tot[i] > 0. && S1f30[i] > 0. &&
        S1f30[i] < 1. && S1tot[i] < S1LimPE[1]) {
      mean += S1f30[i];
      ++j;
      cout << i << "\t" << S1tot[i] << "\t" << S1f30[i] << "\t"
           << driftT[i] / 1e6 << "\t" << S2tot[i] << "\t" << S2width[i] << endl;
      NonBlank.push_back(S1f30[i]);
    } else
      cout << i << endl;
  }
  if (j == 0) {
    cerr << "no events in range!\n";
    return EXIT_FAILURE;
  }
  mean /= double(j);
  std::sort(NonBlank.begin(), NonBlank.end());
  double median = NonBlank[j / 2];

  double sigma = 0.0;
  for (i = 0; i < NUMEVT; ++i)
    if (S1f30[i] != 0. && S1tot[i] > S1LimPE[0] && S1tot[i] < S1LimPE[1])
      sigma += (mean - S1f30[i]) * (mean - S1f30[i]);
  sigma /= (double(j) - 1.);
  sigma = sqrt(sigma);
  double skew = 0.0;
  for (i = 0; i < NUMEVT; ++i)
    if (S1f30[i] != 0. && S1tot[i] > S1LimPE[0] && S1tot[i] < S1LimPE[1])
      skew += pow((S1f30[i] - mean) / sigma, 3.);
  skew /= double(j);
  cerr << "MEDIAN\t\tSIGMA\t\tSKEWNESS\n";
  cerr << median << "\t" << sigma << "\t" << skew
       << endl;  // can switch to "mean," as desired, but change the header name
                 // too :)
  
  return EXIT_SUCCESS;
  
}
