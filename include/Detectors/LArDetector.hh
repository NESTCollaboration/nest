/**
 * @file LArDetector.hh
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief
 * @version
 * @date 2022-04-14
 */
#pragma once

#include "VDetector.hh"

/**
 * @brief An example LAr TPC detector.
 *
 */
class LArDetector : public VDetector {
 public:
  LArDetector();
  ~LArDetector() override = default;
  void Initialization() override;

 private:
};