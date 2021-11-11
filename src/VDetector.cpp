//
// VDetector.cpp
//
// Adapted from Quentin Riffard by Jacob Cutter, May 8, 2018

// *********************************************************************
// THIS DEFAULT VIRTUAL DETECTOR SHOULD ONLY BE MODIFIED BY DEVELOPERS.
// PLEASE DEFINE YOUR OWN DETECTOR (see DetectorExample_XENON10.hh).
// *********************************************************************

#include "VDetector.hh"

VDetector::VDetector() { Initialization(); }

VDetector::~VDetector() = default;

void VDetector::Initialization() {}
