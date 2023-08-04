/**
 * @file LArNESTNeutronCapture.cpp
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief
 * @version 0.1
 * @date 2022-05-23
 */
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "LArDetector.hh"
#include "LArNEST.hh"

int main(int argc, char* argv[]) {
  double density = 1.393;
  uint64_t seed = 0;
  if (argc > 1) {
    density = atof(argv[1]);
  }
  if (argc > 2) {
    seed = atoi(argv[2]);
  }

  std::vector<double> electric_field = {1,    100,  200,  500,  600, 1000,
                                        1500, 2500, 6000, 9000, 9500};

  std::vector<double> energy_vals = {
      167.3,  516.0,  348.7,  867.3,  837.7,  1186.8, 1354.0, 1044.3,
      1881.5, 2229.5, 2566.1, 2432.5, 2781.8, 2842.5, 3111.3, 1972.6,
      2291.6, 2810.5, 3564.7, 3405.5, 2668.1, 3564.5, 2614.3, 3451.8,
      4102.5, 1828.8, 2130.7, 2668.1, 2771.8, 3089.4, 3150.2, 3365.5,
      3405.3, 3700.4, 4745.0, 5063.7, 5582.0, 6098.9};

  LArDetector* detector = new LArDetector();
  // Construct NEST class using detector object
  NEST::LArNEST larnest(detector);
  // Set random seed of RandomGen class
  RandomGen::rndm()->SetSeed(seed);

  // Construct NEST objects for storing calculation results
  NEST::LArNESTResult result;
  std::ofstream output_file;

  output_file.open("neutron_capture.csv");
  output_file << "type,energy,efield,TotalYield,QuantaYield,LightYield,Nph,Ne,"
                 "Nex,Nion\n";

  for (size_t v = 0; v < electric_field.size(); v++) {
    // iterate over energy values
    for (size_t i = 0; i < energy_vals.size(); i++) {
      result = larnest.FullCalculation(NEST::LArInteraction::ER, energy_vals[i],
                                       0, electric_field[v], density, false);
      output_file << "ER,";
      output_file << energy_vals[i] << ",";
      output_file << electric_field[v] << ",";
      output_file << result.yields.TotalYield << ",";
      output_file << result.yields.QuantaYield << ",";
      output_file << result.yields.LightYield << ",";
      output_file << result.yields.Nph << ",";
      output_file << result.yields.Ne << ",";
      output_file << result.yields.Nex << ",";
      output_file << result.yields.Nion << "\n";
    }
  }
  output_file.close();

  return 0;
}