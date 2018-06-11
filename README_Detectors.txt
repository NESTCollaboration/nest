
Creating Detectors:

The "Detectors" folder contains a number of header files:
	- VDetector.hh
	- DetectorExample_XENON10.hh

"VDetector.hh" is a virtual detector class that serves as the base class for building noble element detectors.
The implementation file, "VDetector.cpp" which lives in the main directory, contains some of the associated member functions.
These VDetector files should not be edited (except by developers), as they serve as a default inherited class.

"DetectorExample_XENON10.hh" is the default example detector file, which serves as a template for users to create their own detector. 
To create a detector, one can do (within the Detectors/ folder):
	
	cp DetectorExample_XENON10.hh MyDetector.hh

Then, one can edit the parameters/functions in MyDetector.hh as desired. SEE Detector.txt for more information on important values!
One would start by replacing every instance of "DetectorExample_XENON10" with "MyDetector":

	e.g. 
	
	// DetectorExample_XENON10.hh
	--> // MyDetector.hh

	#ifndef DetectorExample_XENON10_hh
	--> #ifndef MyDetector_hh

	etc.

Once MyDetector.hh is edited, one should edit "testNEST.cpp" in the following way:

	#include "Detectors/DetectorExample_XENON10.hh"
	--> #include "Detectors/MyDetector.hh"

	DetectorExample_XENON10* detector = new DetectorExample_XENON10();
	--> MyDetector* detector = new MyDetector();

Assuming of course that MyDetector.hh contains a class constructor for the class "MyDetector",
this loads the header file and creates a MyDetector object for manipulation.

Then, in testNEST.cpp (or in whatever implementation links the NEST class), one can access various member variables using the 
"get_" functions contained in VDetector.cpp (since the MyDetector class inherits all VDetector members):

	e.g. double g1 = detector->get_g1();
	e.g. double temperature = detector->get_T_Kelvin();

NOTE: Make sure to recompile after changing detector settings!


Modifying Detector Parameters:

Be careful not to get too overwhelmed with all of the detector parameters. The most important ones to change are these to give you a first pass:

- g1
- g1_gas
- E_gas
- geometry parameters, especially Z dependences, but most importantly of all the TopDrift and anode, since they set the size of the S2 gas gap

Reasons: the first sets S1 size of course, while the latter 3 control the S2 size. When running NEST: particle type, energy, field inputs most imp't

* Increases in either g1_gas and/or E_gas and/or the distance from "TopDrift" up to the anode increase the SE (single e-) size
* However, increasing the E_gas also increases the extraction efficiency
* Use case: if you are trying to increase g2, but you know your gas gap and your SE size, you must increase E_gas, but decrease g1_gas to compensate

How to calculate g2 from the detector parameters:
  ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645 = -0.03754*pow(E_gas/1.848,2.)+0.52660*(E_gas/1.848)-0.84645
  SE_size = g1_gas*gasGap_mm*0.1*(0.137*E_gas*1000-177*p_bar-45.7)
  g2 = ExtEff*SE_size
where E_gas, g1_gas, gasGap_mm, and p_bar are provided by the settings file

The next most important thing is the electric field of course, since that changes the recombination probability and thus the "slosh" between S1 & S2.
You can set it with either FitEF, or as one of the run-time input arguments to the testNEST default executable example.
If you follow all of these instructions, you should be getting at least the 0th-order answer. Temperature, pressure, etc. should all be 1st or 2nd-O

OptTrans and SinglePEWaveForm involve quick and dirty ray tracing and pulse shaping, to be used with useTiming = true in analysis.hh
FitS1, FitS2, and FitEF are for 3-D XY,Z position dependencies including field fringing. The first two can just return 1. (get corrected out anyhow)

NEST v2 is now predictive of things like SE size and extraction efficiency so that you don't have to put in manually.
But instead put in the more fundamental parameters listed above, and check the output to see if it is what you want (g2 is broken down: SE x e- ext eff)

For these functions and for all the #'s if you're confused or don't know something:

- Just stick with the XENON10 defaults, which are mixed with some published LUX values.
- Why? They are there because even though they're old detectors now, their #'s ~representative of any detector past/present/future, big/small even


