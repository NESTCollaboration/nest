
#include "RandomGen.hh"

using namespace std;

// Global static pointer used to ensure a single instance of the class.
RandomGen* RandomGen::m_pInstance = nullptr;

// Only allow one instance of class to be generated.
RandomGen* RandomGen::rndm() {
  if (!m_pInstance) m_pInstance = new RandomGen;
  return m_pInstance;
}

std::uint64_t splitmix64(std::uint64_t z) {
  z += 0x9e3779b97f4a7c15;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

void RandomGen::SetSeed(uint64_t s) {
  uint64_t s1 = splitmix64(s);
  rng = xoroshiro128plus64(s1, splitmix64(s1));
}

double RandomGen::rand_uniform() {
  return (static_cast<double>(rng()) - xoroshiro128plus64_min) /
         xoroshiro128plus64_minmax;
}

double RandomGen::rand_gauss(double mean, double sigma) {
  // std::normal_distribution<double> norm(mean, sigma);
  // return norm(rng);
  double u = rand_uniform(), v = rand_uniform();
  return mean + sigma * sqrt2 * sqrt(-log(u)) * cos(two_PI * v);
}

double RandomGen::rand_zero_trunc_gauss(double mean, double sigma) {
  double r = rand_gauss(mean, sigma);
  while (r <= 0) {
    r = rand_gauss(mean, sigma);
  }
  return r;
}

double RandomGen::rand_exponential(double half_life) {
  double r = rand_uniform();
  return -log(1 - r) * half_life / log2;
}

double RandomGen::rand_skewGauss(double xi, double omega, double alpha) {
  double delta = alpha / sqrt(1 + alpha * alpha);
  double gamma1 = four_minus_PI_div_2 *
                  (pow(delta * sqrt2_div_PI, 3.) /
                   pow(1 - 2. * delta * delta / M_PI, 1.5));  // skewness
  double muz = delta * sqrt2_div_PI;
  double sigz = sqrt(1. - muz * muz);
  double m_o;
  if (alpha > 0.) {
    m_o = muz - 0.5 * gamma1 * sigz - 0.5 * exp(-two_PI / alpha);
  }
  if (alpha < 0.) {
    m_o = muz - 0.5 * gamma1 * sigz + 0.5 * exp(two_PI / alpha);
  }
  double mode = xi + omega * m_o;
  // the height should be the value of the PDF at the mode
  double height = exp(-0.5 * (pow((mode - xi) / omega, 2.))) /
                  (sqrt2_PI * omega) *
                  erfc(-1. * alpha * (mode - xi) / omega / sqrt2);
  bool gotValue = false;
  double minX = xi - 6. * omega;
  double maxX =
      xi + 6. * omega;  // +/- 6sigma should be essentially +/- infinity
                        //  can increase these for even better accuracy, at the
                        //  cost of speed
  double testX, testY, testProb;
  while (gotValue == false) {
    testX = minX + (maxX - minX) * RandomGen::rndm()->rand_uniform();
    testY = height *
            RandomGen::rndm()->rand_uniform();  // between 0 and peak height
    // calculate the value of the skewGauss PDF at the test x-value
    testProb = exp(-0.5 * (pow((testX - xi) / omega, 2.))) /
               (sqrt2_PI * omega) *
               erfc(-1. * alpha * (testX - xi) / omega / sqrt2);
    if (testProb >= testY) {
      gotValue = true;
    }
  }
  return testX;
}

int RandomGen::poisson_draw(double mean) {
  return std::poisson_distribution<int>(mean)(rng);
}

int64_t RandomGen::binom_draw(int64_t N0, double prob) {
  return std::binomial_distribution<int64_t>(N0, prob)(rng);
}

int RandomGen::integer_range(int min, int max) {
  return std::uniform_int_distribution<int>(min, max)(rng);
}

vector<double> RandomGen::VonNeumann(double xMin, double xMax, double yMin,
                                     double yMax, double xTest, double yTest,
                                     double fValue) {
  vector<double> xyTry(3);

  xyTry[0] = xTest;
  xyTry[1] = yTest;

  if (xyTry[1] > fValue) {
    xyTry[0] = xMin + (xMax - xMin) * RandomGen::rndm()->rand_uniform();
    xyTry[1] = yMin + (yMax - yMin) * RandomGen::rndm()->rand_uniform();
    xyTry[2] = 1.;
  } else
    xyTry[2] = 0.;

  return xyTry;  // doing a vector means you can return 2 values at the same
                 // time
}

int RandomGen::SelectRanXeAtom() {  // to determine the isotope of Xe
  int A;
  double isotope = rand_uniform() * 100.;

  if (isotope > 0.000 && isotope <= 0.090)
    A = 124;
  else if (isotope > 0.090 && isotope <= 0.180)
    A = 126;
  else if (isotope > 0.180 && isotope <= 2.100)
    A = 128;
  else if (isotope > 2.100 && isotope <= 28.54)
    A = 129;
  else if (isotope > 28.54 && isotope <= 32.62)
    A = 130;
  else if (isotope > 32.62 && isotope <= 53.80)
    A = 131;
  else if (isotope > 53.80 && isotope <= 80.69)
    A = 132;
  else if (isotope > 80.69 && isotope <= 91.13)
    A = 134;
  else
    A = 136;
  return A;
}
