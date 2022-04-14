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
#include <string>
#include <iostream>
#include <stdlib.h>

#include "LArDetector.hh"
#include "LArNEST.hh"

int main(int argc, char* argv[])
{
    // this sample program takes several arguments
    // including
    //  1) num_events - the number of events to generate
    //  2) pdgcode - the pdg code of the particle creating the deposition
    //  3) density - the density of the LAr
    //  4) track_length - the size of the track length for the deposition
    // (optional) 5) seed - a seed for the random number generator
    if (argc < 5) 
    {
        std::cerr << "Error! This program expects six inputs!" << std::endl;
        exit(0);
    }
    size_t num_events = atoi(argv[1]);

    // set up energy steps
    double energy = atof(argv[2]);
    int num_energy_steps = 10000;
    std::vector<double> energy_vals;
    double start_val = .0001;
    double end_val = 1;
    double step_size = (end_val - start_val)/num_energy_steps;
    for (size_t i = 0; i < num_energy_steps; i++)
    {
        energy_vals.emplace_back(start_val + step_size * i);
    }
    
    int pdgcode = atoi(argv[2]);
    double density = atof(argv[3]);
    std::vector<double> electric_field = {
        1, 50, 100, 200, 500, 1000, 1500, 2000
    };

    double track_length = atof(argv[4]);
    uint64_t seed = 0;
    if (argc > 4) {
        seed = atoi(argv[4]);
    }

    LArDetector* detector = new LArDetector();
    // Construct NEST class using detector object
    NEST::LArNEST larnest(detector);
    // Set random seed of RandomGen class
    RandomGen::rndm()->SetSeed(seed);

    // Construct NEST objects for storing calculation results
    NEST::QuantaResult quanta;
    std::ofstream output_file;
    output_file.open("legacy_larnest_output_" + std::to_string(pdgcode) + ".csv");
    output_file << "energy,efield,photons,electrons\n";
    // iterate over number of events
    for (size_t v = 0; v < electric_field.size(); v++)
    {
        for (size_t i = 0; i < num_energy_steps; i++)
        {
            for (size_t j = 0; j < num_events; j++)
            {
                quanta = larnest.LegacyCalculation(
                    pdgcode,
                    energy_vals[i],
                    density,
                    electric_field[v],
                    track_length
                );
                output_file << energy_vals[i] << ","; 
                output_file << electric_field[v] << ",";
                output_file << quanta.photons << "," << quanta.electrons << "\n";
            }
        }
    }
    output_file.close();
    return 0;
}