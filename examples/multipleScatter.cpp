/*
 * File:   multipleScatter.cpp
 * Author: Greg Rischbeter
 *
 * Created on January 28, 2022
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
  if (argc < 4) {
    cerr << "This program aims to take execNEST outputs, and combine energy "
            "deposits into Double Scatters. "
         << endl;
    cerr << endl
         << "Usage: " << endl
         << "./multipleScatter <nEvents> <nScatters> <Pulse Type> <File 1> "
            "<optional: File 2> "
         << endl;

    cerr << "  -- <nEvents> is the number of events to stitch together. "
         << endl;
    cerr << "  -- <nScatters> is the number of energy deposits to combine; use "
            "-1 for a random exponential draw. "
         << endl;
    cerr << "  -- <Pulse Type> is either 0 (for uncorrected pulses) or 1 (for "
            "corrected pulses); corrected pulses sum S2s and take an "
            "S2-weighted average S1c. "
         << endl;
    cerr << "  -- <File 1> is execNEST output events used to stitch together "
            "multiple scatters."
         << endl;
    cerr << "  -- <File 2> is optional, and will be the last scatter in each "
            "event, if provided."
         << endl;
    return 1;
  }
  // Load the number of events
  int nEvents = atoi(argv[1]);
  int nScatters = atoi(argv[2]);
  // RandomGen::rndm()->SetSeed(time(nullptr));
  bool randomScatters = false;
  if (nScatters < 2) {
    randomScatters = true;
  }
  // define parameters for random exponential draw of nScatters
  // if using nScatters = -1.
  // These are just 0th order fill-ins to produce a spread of values.
  double halfLife = 1.;
  double minScatters = 1.5;  // int-rounding will make the true minimum = 2
  double maxScatters = 50.;
  // Check if using corrected or uncorrected output pulse areas
  bool corrected = true;
  if (atoi(argv[3]) == 0) corrected = false;
  if (verbosity > 0) {
    if (corrected)
      cout << "S1c [phd]\tS2c [phd]\tnScatters" << endl;
    else
      cout << "S1 [phd]\tS2 [phd]\tnScatters" << endl;
  }
  bool sameFile = true;
  if (argc > 5) {
    sameFile = false;
    // cerr << "Using two files to generate scatters." << endl;
  }
  // open the first file and get the S1, S2, and their correction factors
  FILE *file1 = fopen(argv[4], "r");
  double a, b, c, d, e, f, g, h, i, j, k, l, m, n;
  int ch, nLines = 0;
  vector<double> S1a, S2a;
  vector<double> S1a_corFactor, S2a_corFactor;
  if (verbosity > 0) {
    while (EOF != (ch = getc(file1))) {
      if ('\n' == ch && nLines)
        break;
      else
        nLines = 0;
      if (']' == ch && nLines == 0) ++nLines;
    }
  }
  while (1) {
    int scan = fscanf(
        file1,
        "%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
        &a, &b, &c, &d, &e, &f, &g, &h, &i, &j, &k, &l, &m, &n);
    if (feof(file1)) break;
    // fprintf(stderr,"%.6f\t%.6f\t%.6f\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
    if ((fabs(n) > 1.) && (fabs(m) > 1.)) {
      S1a.push_back(fabs(j));
      S2a.push_back(fabs(m));
      S1a_corFactor.push_back(fabs(k / j));
      S2a_corFactor.push_back(fabs(n / m));
    }
  }  // Done with the first file....
  fclose(file1);
  if (!sameFile) {
    cerr << "Using two files to generate scatters..." << endl;
    // Load the 2nd file if provided, and stitch the two together to make
    // multiple scatter events scatters
    FILE *file2 = fopen(argv[5], "r");
    nLines = 0;
    double a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2, l2, m2, n2;
    vector<double> S1b, S2b;
    vector<double> S1b_corFactor, S2b_corFactor;
    if (verbosity > 0) {
      while (EOF != (ch = getc(file2))) {
        if ('\n' == ch && nLines)
          break;
        else
          nLines = 0;
        if (']' == ch && nLines == 0) ++nLines;
      }
    }

    while (1) {
      int scan2 = fscanf(
          file2,
          "%lf\t%lf\t%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
          &a2, &b2, &c2, &d2, &e2, &f2, &g2, &h2, &i2, &j2, &k2, &l2, &m2, &n2);
      if (feof(file2)) break;
      // fprintf(stderr,"%.6f\t%.6f\t%.6f\t%.6f,%.6f,%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",a,b,c,d,e,f,g,h,i,j,k,l,m,n);
      if ((fabs(n2) > 1.) &&
          (fabs(m2) > 1.)) {      // prevent loading events where S2 ~ PHE_MIN
        S1b.push_back(fabs(j2));  // was causing memory leaks
        S2b.push_back(fabs(m2));
        S1b_corFactor.push_back(fabs(k2 / j2));
        S2b_corFactor.push_back(fabs(n2 / m2));
      }
    }  // Done with the second file
    fclose(file2);

    auto size1 = S1a.size();
    auto size2 = S1b.size();
    double totalS1, totalS2, totalCorFactor;
    double thisS1a, thisS1b, thisS2a, thisS2b;
    double thisS1a_corFactor, thisS2a_corFactor, thisS1b_corFactor,
        thisS2b_corFactor;
    for (int ii = 0; ii < nEvents; ii++) {
      while (nScatters < 2) {
        nScatters =  // 0th order attempt at reproducing multiple scatter
                     // behavior
            (int)(RandomGen::rndm()->rand_exponential(halfLife, minScatters,
                                                      maxScatters) +
                  0.5);
      }
      totalS1 = 0.;
      totalCorFactor = 0.;
      totalS2 = 0.;
      for (int jj = 0; jj < nScatters - 1; jj++) {
        // randomly draw an energy deposit from the first file
        int index = RandomGen::rndm()->integer_range(0, (int)size1 - 1);

        thisS1a = S1a[index];
        thisS2a = S2a[index];
        thisS1a_corFactor = S1a_corFactor[index];
        thisS2a_corFactor = S2a_corFactor[index];
        if (corrected) {
          double cor_S2a = thisS2a * thisS2a_corFactor;
          totalS2 += cor_S2a;
          totalS1 += thisS1a;
          totalCorFactor += thisS1a_corFactor * cor_S2a;
        } else {
          totalS1 += thisS1a;
          totalS2 += thisS2a;
        }
      }  // now on to the last scatter from the 2nd file
      int index2 = RandomGen::rndm()->integer_range(0, (int)size2 - 1);
      thisS1b = S1b[index2];
      thisS2b = S2b[index2];
      thisS1b_corFactor = S1b_corFactor[index2];
      thisS2b_corFactor = S2b_corFactor[index2];

      if (corrected) {
        double cor_S2b = thisS2b * thisS2b_corFactor;
        totalS2 += cor_S2b;
        totalS1 += thisS1b;
        totalCorFactor += thisS1b_corFactor * cor_S2b;
        totalS1 *= totalCorFactor / totalS2;
      } else {
        totalS1 += thisS1b;
        totalS2 += thisS2b;
      }
      cout << totalS1 << "\t" << totalS2 << "\t" << nScatters
           << endl;  //"\t" << thisS1a << "  " << thisS1b << " " << thisS2a << "
                     //" << thisS2b << endl;
      if (randomScatters) nScatters = -1;
    }
    return 0;
  } else {
    cerr << "Generating multiple scatters from file..." << endl;
    // Drawing both scatters from the same file...
    int nGenerated = 0;

    auto size1 = S1a.size();
    double totalS1, totalS2, totalCorFactor;
    double thisS1a, thisS1b, thisS2a, thisS2b;
    double thisS1a_corFactor, thisS2a_corFactor, thisS1b_corFactor,
        thisS2b_corFactor;

    while (nGenerated < nEvents) {
      while (nScatters < 2) {
        nScatters =  // 0th order attempt at reproducing multiple scatter
                     // behavior
            (int)(RandomGen::rndm()->rand_exponential(halfLife, minScatters,
                                                      maxScatters) +
                  0.5);
      }
      totalS1 = 0.;
      totalCorFactor = 0.;
      totalS2 = 0.;
      for (int jj = 0; jj < nScatters; jj++) {
        // randomly draw an energy deposit from the first file
        int index = RandomGen::rndm()->integer_range(0, (int)size1 - 1);

        thisS1a = S1a[index];
        thisS2a = S2a[index];
        thisS1a_corFactor = S1a_corFactor[index];
        thisS2a_corFactor = S2a_corFactor[index];
        if (corrected) {
          double cor_S2a = thisS2a * thisS2a_corFactor;
          totalS2 += cor_S2a;
          totalS1 += thisS1a;
          totalCorFactor += thisS1a_corFactor * cor_S2a;
        } else {
          totalS1 += thisS1a;
          totalS2 += thisS2a;
        }
      }
      if (corrected)
        totalS1 *= totalCorFactor / totalS2;  // S2c-weighted average corrected
                                              // S1

      cout << totalS1 << "\t" << totalS2 << "\t" << nScatters
           << endl;  //"\t" << thisS1a << "  " << thisS1b << " " << thisS2a << "
                     //" << thisS2b << endl;
      nGenerated++;
      if (randomScatters) nScatters = -1;
    }
    return 0;
  }
}
