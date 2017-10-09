
#include <NEST.hh>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace NEST;

/*
 * 
 */


int main(int argc, char** argv) {

    // User can modify
    std::string output_file_name ("rand_exp_ns_arr.txt");
    int nIterations = 200000;
    double half_life = 100.0;

    // Setup
    NEST::NESTcalc n;
    n.rng.seed(std::random_device()());
    double time_ns[nIterations];
    ofstream fout;
    fout.open(output_file_name);

    // Run many iterations to find photon timing
    for (int j = 0; j < nIterations; j++) {
        time_ns[j] = n.rand_exponential(half_life);
	fout << time_ns[j] << "\n";
    }

    fout.close();

    return 0;

}
