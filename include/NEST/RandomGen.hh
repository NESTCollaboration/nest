//
// RandomGen.hh
//
// Created by Jacob Cutter, May 10 2018

#ifndef RANDOMGEN_HH
#define RANDOMGEN_HH 1
#include "xoroshiro.hh"
#include <math.h>
#include <stdlib.h>
#include <random>
#include <vector>

using namespace std;

class RandomGen {
 public:
  static RandomGen* rndm();
  void SetSeed(unsigned long int s);
  double rand_uniform();
  double rand_gauss(double mean, double sigma);
  double rand_exponential(double half_life);
  double rand_skewGauss( double xi, double omega, double alpha);
  int poisson_draw(double mean);
  int integer_range(int min, int max);
  vector<double> VonNeumann(double xMin, double xMax, double yMin, double yMax,
                            double xTest, double yTest, double fValue);
  int SelectRanXeAtom();

 private:
  // Random number generator object for this class only
  //  std::ranlux24 rng;
  xoroshiro128plus64 rng;

  RandomGen(){};                // private so that it cannot be manually called
  RandomGen(RandomGen const&);  // copy constructor is private
  void operator=(RandomGen const&);  // assignment operator is private
  static RandomGen* m_pInstance;     // private class pointer instance
};

#endif
