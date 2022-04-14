/**
 * @file legacyLArNEST.cpp
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 
 * @date 2022-04-14
 */
#include <fstream>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "LArDetector.hh"
#include "LArNEST.hh"

int main(int argc, char* argv[])
{
    // this sample program takes several arguments
    // including
    //  1) num_events - the number of events to generate
    //  2) energy - the energy of the depositions in MeV
    //  3) pdgcode - the pdg code of the particle creating the deposition
    //  4) density - the density of the LAr
    //  5) eField -  the value of the electric field in V/cm
    //  6) track_length - the size of the track length for the deposition
    // (optional) 7) seed - a seed for the random number generator
    if (argc < 7) 
    {
        std::cerr << "Error! This program expects six inputs!" << std::endl;
        exit(0);
    }
    size_t num_events = atoi(argv[1]);
    double energy = atof(argv[2]);
    std::vector<double> energy_vals;
    if (energy != 0.0) {
        energy_vals.resize(num_events, energy);
    }
    else {
        double start_val = .0001;
        double end_val = 1;
        double step_size = (end_val - start_val)/num_events;
        for (size_t i = 0; i < num_events; i++)
        {
            energy_vals.emplace_back(start_val + step_size * i);
        }
    }
    int pdgcode = atoi(argv[3]);
    double density = atof(argv[4]);
    double eField = atof(argv[5]);
    double track_length = atof(argv[6]);
    uint64_t seed = 0;
    if (argc > 7) {
        seed = atoi(argv[7]);
    }

    LArDetector* detector = new LArDetector();
    // Construct NEST class using detector object
    NEST::LArNEST larnest(detector);
    // Set random seed of RandomGen class
    RandomGen::rndm()->SetSeed(seed);

    // Construct NEST objects for storing calculation results
    NEST::QuantaResult quanta;
    std::ofstream output_file;
    output_file.open("legacy_larnest_output.csv");
    output_file << "energy,photons,electrons\n";
    // iterate over number of events
    for (size_t i = 0; i < num_events; i++)
    {
        quanta = larnest.LegacyCalculation(
            pdgcode=pdgcode,
            energy=energy_vals[i],
            density=density,
            eField=eField,
            track_length=track_length
        );
        output_file << energy_vals[i] << "," << quanta.photons << "," << quanta.electrons << "\n";
    }
    output_file.close();
    return 0;
}