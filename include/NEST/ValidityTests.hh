#ifndef VALIDITYTESTS_HH
#define VALIDITYTESTS_HH 1

#include <cfloat>
#include <math.h>
#include <iostream>
using namespace std;

class ValidityTests {
 public:
  static bool nearlyEqual(double a, double b, double epsilon = 1e-9);
};

#endif
