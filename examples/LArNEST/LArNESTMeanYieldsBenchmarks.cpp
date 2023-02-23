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
  if (argc > 1) {
    density = atof(argv[1]);
  }
  if (argc > 2) {
    seed = atoi(argv[2]);
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

  output_file.open("mean_yields_benchmarks.csv");
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
    for (size_t v = 0; v < electric_field.size(); v++) {
      // iterate over energy values
      for (size_t i = 0; i < num_energy_steps; i++) {
        result = larnest.FullCalculation(particle_types[k], energy_vals[i], 0,
                                         electric_field[v], density, false);
        output_file << particle_type[k] << ",";
        output_file << energy_vals[i] << ",";
        output_file << electric_field[v] << ",";
        output_file << result.yields.TotalYield << ",";
        output_file << result.yields.QuantaYield << ",";
        output_file << result.yields.LightYield << ",";
        output_file << result.yields.Nph << ",";
        output_file << result.yields.Ne << ",";
        output_file << result.yields.Nex << ",";
        output_file << result.yields.Nion << ",0,0,0,0,0,0,0\n";
      }
    }
  }
  output_file.close();
  return 0;
}