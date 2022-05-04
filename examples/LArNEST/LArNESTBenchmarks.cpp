/**
 * @file legacyLArNEST.cpp
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief Benchmarks for LAr ER model.  This program generates ...
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
    //  1) num_events - the number of events to generate for each field/energy value
    //  2) type - NR or ER for now 
    //  3) density - the density of the LAr to use
    // (optional) 4) seed - a seed for the random number generator
    if (argc < 4) 
    {
        std::cerr << "Error! This program expects three inputs!" << std::endl;
        exit(0);
    }
    size_t num_events = atoi(argv[1]);

    // set up energy steps
    int num_energy_steps = 50000;
    std::vector<double> energy_vals; 
    double start_val = .1;
    double end_val = 1000;
    double step_size = (end_val - start_val)/num_energy_steps;
    for (size_t i = 0; i < num_energy_steps; i++)
    {
        energy_vals.emplace_back(start_val + step_size * i);
    }
    
    std::string particle_type = argv[2];
    std::vector<double> electric_field;
    NEST::LArInteraction incident_particle;
    if (particle_type == "NR")
    {
        electric_field.insert(electric_field.end(), {
            1, 50, 100, 200, 500, 1000, 1500, 2000
        });
        incident_particle = NEST::LArInteraction::NR;
    }
    else if (particle_type == "ER")
    {
        electric_field.insert(electric_field.end(), {
            1, 100, 200, 600, 1000, 1500, 2500, 6000, 9000
        });
        incident_particle = NEST::LArInteraction::ER;
    }
    else
    {
        electric_field.insert(electric_field.end(), {
            1, 50, 100, 500, 1000, 5000, 10000, 20000
        });
        incident_particle = NEST::LArInteraction::Alpha;
    }
    
    double density = atof(argv[3]);
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
    NEST::LArNESTResult result;
    std::ofstream output_file;
    if (particle_type == "NR") {
        output_file.open("larnest_benchmarks_nr.csv");
    }
    else if (particle_type == "ER") {
        output_file.open("larnest_benchmarks_er.csv");
    }
    else {
        output_file.open("larnest_benchmarks_alpha.csv");
    }
    output_file << "energy,efield,TotalYield,QuantaYield,LightYield,Nph,Ne,Nex,Nion\n";
    // iterate over electric field values
    for (size_t v = 0; v < electric_field.size(); v++)
    {
        // iterate over energy values
        for (size_t i = 0; i < num_energy_steps; i++)
        {
            std::vector<double> TotalYield(num_events, 0.0);
            std::vector<double> QuantaYield(num_events, 0.0);
            std::vector<double> LightYield(num_events, 0.0);
            std::vector<double> Nph(num_events, 0.0);
            std::vector<double> Ne(num_events, 0.0);
            std::vector<double> Nex(num_events, 0.0);
            std::vector<double> Nion(num_events, 0.0);
            // collect statistics for each event
            for (size_t j = 0; j < num_events; j++)
            {
                result = larnest.FullCalculation(
                    incident_particle,
                    energy_vals[i],
                    electric_field[v],
                    density,
                    false
                );
                TotalYield[j] = result.yields.TotalYield;
                QuantaYield[j] = result.yields.QuantaYield;
                LightYield[j] = result.yields.LightYield;
                Nph[j] = result.yields.Nph;                
                Ne[j] = result.yields.Ne;
                Nex[j] = result.yields.Nex;
                Nion[j] = result.yields.Nion;
            }
            double TotalYield_mean = std::accumulate(TotalYield.begin(), TotalYield.end(), 0.0) / num_events;
            double QuantaYield_mean = std::accumulate(QuantaYield.begin(), QuantaYield.end(), 0.0) / num_events;
            double LightYield_mean = std::accumulate(LightYield.begin(), LightYield.end(), 0.0) / num_events;
            double Nph_mean = std::accumulate(Nph.begin(), Nph.end(), 0.0) / num_events;
            double Ne_mean = std::accumulate(Ne.begin(), Ne.end(), 0.0) / num_events;
            double Nex_mean = std::accumulate(Nex.begin(), Nex.end(), 0.0) / num_events;
            double Nion_mean = std::accumulate(Nion.begin(), Nion.end(), 0.0) / num_events;
            output_file << energy_vals[i] << ","; 
            output_file << electric_field[v] << ",";
            output_file << TotalYield_mean << "," << QuantaYield_mean << "," << LightYield_mean << ",";
            output_file << Nph_mean << "," << Ne_mean << "," << Nex_mean << "," << Nion_mean << "\n";
        }
    }
    output_file.close();
    return 0;
}