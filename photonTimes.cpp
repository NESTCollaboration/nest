
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
    INTERACTION_TYPE type_num = Kr83m;
    bool exciton = true;
    std::string output_file_name ("time_ns_arr.txt");
    int nIterations = 200000;

    // Setup
    NEST::NESTcalc n;
    double time_ns[nIterations];
    ofstream fout;
    fout.open(output_file_name);

    // Run many iterations to find photon timing
    for (int j = 0; j < nIterations; j++) {
        time_ns[j] = n.PhotonTime(type_num, exciton);
	fout << time_ns[j] << "\n";
    }

    fout.close();

    return 0;

}
