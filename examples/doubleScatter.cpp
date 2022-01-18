/*
 * File:   doubleScatter.cpp
 * Author: Greg Rischbeter
 *
 * Created on January 18, 2022
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "RandomGen.hh"
#include "analysis.hh"

using namespace std;

int main(int argc, char **argv) {

  if (argc < 3) {

    cerr << "This program aims to take execNEST outputs, and combine energy "
            "deposits into Double Scatters. "
         << endl;
    cerr
        << endl
        << "Usage: " << endl
        << "./doubleScatter <nEvents> <Pulse Type> <File 1> <optional: File 2> "
        << endl;

    cerr << "  -- <nEvents> is the number of double scatters to stitch "
            "together. "
         << endl;
    cerr << "  -- <Pulse Type> is either 0 (for uncorrected pulses) or 1 (for "
            "corrected pulses); corrected pulses sum S2s and take an "
            "S2-weighted average S1c. "
         << endl;
    cerr << "  -- <File 1> is execNEST output events used to stitch together "
            "double scatters."
         << endl;
    cerr << "  -- <File 2> is optional, and will be the second scatter in each "
            "event, if provided."
         << endl;

    return 1;
  }
  // Load the number of events
  int nEvents = atoi(argv[1]);
  // RandomGen::rndm()->SetSeed(time(nullptr));

  // Check if using corrected or uncorrected output pulse areas
  bool corrected = true;
  if (atoi(argv[2]) == 0)
    corrected = false;

  if (verbosity) {
    if (corrected)
      cout << "S1c [phd]\tS2c [phd]" << endl;
    else
      cout << "S1 [phd]\tS2 [phd]" << endl;
  }

  bool sameFile = true;
  if (argc > 4) {
    sameFile = false;
    // cerr << "Using two files to generate double scatters." << endl;
  }
  // open the first file and get the S1, S2, and their correction factors
  FILE *file1 = fopen(argv[3], "r");
  double a, b, c, d, e, f, g, h, i, j, k, l, m, n;
  int ch, nLines = 0, o;
  vector<double> S1a, S2a;
  vector<double> S1a_corFactor, S2a_corFactor;
  if (verbosity) {
    while (EOF != (ch = getc(file1))) {
      if ('\n' == ch && nLines)
        break;
      else
        nLines = 0;
      if (']' == ch && nLines == 0)
        ++nLines;
    }
  }

  while (1) {
    int scan = fscanf(
        file1,
        "%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
        &a, &b, &c, &d, &e, &f, &g, &h, &i, &j, &k, &l, &m, &n);
    if (feof(file1))
      break;
    // fprintf(stderr,"%.6f\t%.6f\t%.6f\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
    S1a.push_back(fabs(j));
    S2a.push_back(fabs(m));
    S1a_corFactor.push_back(fabs(k / j));
    S2a_corFactor.push_back(fabs(n / m));
  }
  fclose(file1);

  if (!sameFile) {
    cerr << "Using two files to generate double scatters..." << endl;
    // Load the 2nd file if provided, and stitch the two together to make double
    // scatters
    FILE *file2 = fopen(argv[4], "r");
    ch, nLines = 0, o;
    double a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2, l2, m2, n2;
    vector<double> S1b, S2b;
    vector<double> S1b_corFactor, S2b_corFactor;
    if (verbosity) {
      while (EOF != (ch = getc(file2))) {
        if ('\n' == ch && nLines)
          break;
        else
          nLines = 0;
        if (']' == ch && nLines == 0)
          ++nLines;
      }
    }

    while (1) {
      int scan2 = fscanf(
          file2,
          "%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
          &a2, &b2, &c2, &d2, &e2, &f2, &g2, &h2, &i2, &j2, &k2, &l2, &m2, &n2);
      if (feof(file2))
        break;
      // fprintf(stderr,"%.6f\t%.6f\t%.6f\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
      S1b.push_back(fabs(j2));
      S2b.push_back(fabs(m2));
      S1b_corFactor.push_back(fabs(k2 / j2));
      S2b_corFactor.push_back(fabs(n2 / m2));
    }
    fclose(file2);

    auto size1 = S1a.size();
    auto size2 = S1b.size();
    double totalS1, totalS2;
    double thisS1a, thisS1b, thisS2a, thisS2b;
    double thisS1a_corFactor, thisS2a_corFactor, thisS1b_corFactor,
        thisS2b_corFactor;
    for (int ii = 0; ii < nEvents; ii++) {
      int i1 = RandomGen::rndm()->integer_range(0, (int)size1);
      int i2 = RandomGen::rndm()->integer_range(0, (int)size2);
      // Draw deposits from two different sets of S1 and S2 vectors
      thisS1a = S1a[i1];
      thisS1b = S1b[i2];
      thisS2a = S2a[i1];
      thisS2b = S2b[i2];
      thisS1a_corFactor = S1a_corFactor[i1];
      thisS1b_corFactor = S1b_corFactor[i2];
      thisS2a_corFactor = S2a_corFactor[i1];
      thisS2b_corFactor = S2b_corFactor[i2];

      if (corrected) {
        double cor_S2a = thisS2a * thisS2a_corFactor;
        double cor_S2b = thisS2b * thisS2b_corFactor;
        totalS2 = cor_S2a + cor_S2b;    // total corrected S2 area
        totalS1 = (thisS1a + thisS1b) * // total uncorrected S1 area
                  (thisS1a_corFactor * cor_S2a + thisS1b_corFactor * cor_S2b) /
                  totalS2; // S2c-weighted average correction factor
      } else {
        totalS1 = (thisS1a + thisS1b); // total uncorrected S1 area
        totalS2 = (thisS2a + thisS2b); // total uncorrected S1 area
      }

      cout << totalS1 << "\t" << totalS2
           << endl; //"\t" << thisS1a << "  " << thisS1b << " " << thisS2a << "
                    //" << thisS2b << endl;
    }
    return 0;
  } else {
    cerr << "Generating double scatters from file..." << endl;
    // Drawing both scatters from the same file...
    int nGenerated = 0;

    auto size1 = S1a.size();
    double totalS1, totalS2;
    double thisS1a, thisS1b, thisS2a, thisS2b;
    double thisS1a_corFactor, thisS2a_corFactor, thisS1b_corFactor,
        thisS2b_corFactor;

    while (nGenerated < nEvents) {

      int i1 = RandomGen::rndm()->integer_range(0, (int)size1);
      int i2 = RandomGen::rndm()->integer_range(0, (int)size1);
      if (i1 == i2)
        continue;
      // Draw deposits from the same S1 and S2 vectors
      thisS1a = S1a[i1];
      thisS1b = S1a[i2];
      thisS2a = S2a[i1];
      thisS2b = S2a[i2];
      thisS1a_corFactor = S1a_corFactor[i1];
      thisS1b_corFactor = S1a_corFactor[i2];
      thisS2a_corFactor = S2a_corFactor[i1];
      thisS2b_corFactor = S2a_corFactor[i2];

      if (corrected) {
        double cor_S2a = thisS2a * thisS2a_corFactor;
        double cor_S2b = thisS2b * thisS2b_corFactor;
        totalS2 = cor_S2a + cor_S2b;    // total corrected S2 area
        totalS1 = (thisS1a + thisS1b) * // total uncorrected S1 area
                  (thisS1a_corFactor * cor_S2a + thisS1b_corFactor * cor_S2b) /
                  totalS2; // S2c-weighted average correction factor
      } else {
        totalS1 = (thisS1a + thisS1b); // total uncorrected S1 area
        totalS2 = (thisS2a + thisS2b); // total uncorrected S1 area
      }

      cout << totalS1 << "\t" << totalS2
           << endl; //"\t" << thisS1a << "  " << thisS1b << " " << thisS2a << "
                    //" << thisS2b << endl;
      nGenerated++;
    }
    return 0;
  }
}
