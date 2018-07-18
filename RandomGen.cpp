
#include "RandomGen.hh"

using namespace std;

// Global static pointer used to ensure a single instance of the class.
RandomGen* RandomGen::m_pInstance = NULL;

// Only allow one instance of class to be generated.
RandomGen* RandomGen::rndm() {
  if (!m_pInstance) m_pInstance = new RandomGen;
  return m_pInstance;
}

void RandomGen::SetSeed(unsigned long int s) { rng.seed(s); }

double RandomGen::rand_uniform() {
  return (double)(rng() - rng.min()) / (double)(rng.max() - rng.min());
}

double RandomGen::rand_gauss(double mean, double sigma) {
  double u = rand_uniform(), v = rand_uniform();
  return mean + sigma * sqrt(-2. * log(u)) * cos(2. * M_PI * v);
}

double RandomGen::rand_exponential(double half_life) {
  double r = rand_uniform();
  return log(1 - r) * -1 * half_life / log(2.);
}

int RandomGen::poisson_draw(double mean) {
  std::poisson_distribution<int> distribution(mean);
  return distribution(rng);
}

int RandomGen::integer_range(int min, int max) {
  std::uniform_int_distribution<int> distribution(min, max);
  return distribution(rng);
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

int RandomGen::SelectRanXeAtom() {
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
