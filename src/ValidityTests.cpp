
#include "ValidityTests.hh"

using namespace std;

bool ValidityTests::nearlyEqual(double a, double b, double epsilon) {
  // Handles cases of equality statements for floats
  if (a == b) {  // shortcut, handles infinities
    return true;
  }
  double absA = fabs(a);
  double absB = fabs(b);
  double diff = fabs(a - b);

  if (a == 0 || b == 0 || (absA + absB < DBL_MIN)) {
    // a or b is zero or both are extremely close to it
    // relative error is less meaningful here
    return diff < (epsilon * DBL_MIN);
  }  // use relative error
  return diff / fmin((absA + absB), DBL_MAX) < epsilon;
}
