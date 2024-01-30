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

int main(int argc, char* argv[]) {
  // this sample program takes can take two arguments,
  // (optional) 1) density - the density of the LAr to use
  // (optional) 2) seed - a seed for the random number generator

  double density = 1.393;
  uint64_t seed = 0;
  uint64_t num_events = 100;
  if (argc > 1) {
    density = atof(argv[1]);
  }
  if (argc > 2) {
    seed = atoi(argv[2]);
  }
  if (argc > 3) {
    num_events = atoi(argv[3]);
  }

  // set up energy steps
  int num_energy_steps = 50000;
  std::vector<double> energy_vals;
  double start_val = .1;
  double end_val = 1000;
  double step_size = (end_val - start_val) / num_energy_steps;

  for (size_t i = 0; i < num_energy_steps; i++) {
    energy_vals.emplace_back(start_val + step_size * i);
  }
  std::vector<NEST::LArInteraction> particle_types = {
      NEST::LArInteraction::NR, NEST::LArInteraction::ER,
      NEST::LArInteraction::Alpha};
  std::vector<std::string> particle_type = {"NR", "ER", "Alpha"};

  LArDetector* detector = new LArDetector();
  // Construct NEST class using detector object
  NEST::LArNEST larnest(detector);
  // Set random seed of RandomGen class
  RandomGen::rndm()->SetSeed(seed);

  // Construct NEST objects for storing calculation results
  NEST::LArNESTResult result;
  std::ofstream output_file;

  output_file.open("fluctuation_benchmarks.csv");
  output_file << "type,energy,efield,TotalYield,QuantaYield,LightYield,Nph,Ne,"
                 "Nex,Nion,TotalYield_std,QuantaYield_std,LightYield_std,Nph_"
                 "std,Ne_std,Nex_std,Nion_std\n";

  // iterate over electric field values
  for (size_t k = 0; k < particle_types.size(); k++) {
    std::vector<double> electric_field;
    if (particle_types[k] == NEST::LArInteraction::NR) {
      electric_field.insert(electric_field.end(), {1, 50, 100, 200, 250, 500,
                                                   1000, 1500, 1750, 2000});
    } else if (particle_types[k] == NEST::LArInteraction::ER) {
      electric_field.insert(electric_field.end(), {1, 100, 200, 600, 1000, 1500,
                                                   2500, 6000, 9000, 9500});
    } else {
      electric_field.insert(electric_field.end(),
                            {1, 50, 100, 500, 1000, 5000, 10000, 20000});
    }
    // iterate over number of events
    for (size_t v = 0; v < electric_field.size(); v++) {
      for (size_t i = 0; i < num_energy_steps; i++) {
        std::vector<double> Nph(num_events);
        std::vector<double> Ne(num_events);
        std::vector<double> Nex(num_events);
        std::vector<double> Nion(num_events);
        // collect statistics for each event
        for (size_t j = 0; j < num_events; j++) {
          result = larnest.FullCalculation(particle_types[k], energy_vals[i], 0,
                                           electric_field[v], density, false);
          Nph[j] = double(result.fluctuations.NphFluctuation);
          Ne[j] = double(result.fluctuations.NeFluctuation);
          Nex[j] = double(result.fluctuations.NexFluctuation);
          Nion[j] = double(result.fluctuations.NionFluctuation);
        }
        double Nph_mean =
            std::accumulate(Nph.begin(), Nph.end(), 0.0) / double(num_events);
        double Ne_mean =
            std::accumulate(Ne.begin(), Ne.end(), 0.0) / double(num_events);
        double Nex_mean =
            std::accumulate(Nex.begin(), Nex.end(), 0.0) / double(num_events);
        double Nion_mean =
            std::accumulate(Nion.begin(), Nion.end(), 0.0) / double(num_events);

        double Nq_mean = Nph_mean + Ne_mean;
        double Nq_std = 0.0;
        double Nph_std = 0.0;
        double Ne_std = 0.0;
        double Nex_std = 0.0;
        double Nion_std = 0.0;

        for (size_t ii = 0; ii < Nph.size(); ii++) {
          Nq_std += (Nph[ii] + Ne[ii] - Nq_mean) * (Nph[ii] + Ne[ii] - Nq_mean);
        }
        Nq_std = std::sqrt(Nq_std / double(num_events));

        for (size_t ii = 0; ii < Nph.size(); ii++) {
          Nph_std += (Nph[ii] - Nph_mean) * (Nph[ii] - Nph_mean);
        }
        Nph_std = std::sqrt(Nph_std / double(num_events));

        for (size_t ii = 0; ii < Ne.size(); ii++) {
          Ne_std += (Ne[ii] - Ne_mean) * (Ne[ii] - Ne_mean);
        }
        Ne_std = std::sqrt(Ne_std / double(num_events));

        for (size_t ii = 0; ii < Nex.size(); ii++) {
          Nex_std += (Nex[ii] - Nex_mean) * (Nex[ii] - Nex_mean);
        }
        Nex_std = std::sqrt(Nex_std / double(num_events));

        for (size_t ii = 0; ii < Nion.size(); ii++) {
          Nion_std += (Nion[ii] - Nion_mean) * (Nion[ii] - Nion_mean);
        }
        Nion_std = std::sqrt(Nion_std / double(num_events));

        output_file << particle_type[k] << ",";
        output_file << energy_vals[i] << ",";
        output_file << electric_field[v] << ",";
        output_file << (Nq_mean / energy_vals[i]) << ",";
        output_file << (Ne_mean / energy_vals[i]) << ",";
        output_file << (Nph_mean / energy_vals[i]) << ",";
        output_file << Nph_mean << ",";
        output_file << Ne_mean << ",";
        output_file << Nex_mean << ",";
        output_file << Nion_mean << ",";
        output_file << (Nq_std / energy_vals[i]) << ",";
        output_file << (Ne_std / energy_vals[i]) << ",";
        output_file << (Nph_std / energy_vals[i]) << ",";
        output_file << Nph_std << ",";
        output_file << Ne_std << ",";
        output_file << Nex_std << ",";
        output_file << Nion_std << "\n";
      }
    }
  }
  output_file.close();
  return 0;
}