Installation Instructions:

1. Make an empty directory somewhere outside the NEST source code directory, and go there.
   For example: (from the source code directory) "mkdir ../nest-build; cd ../nest-build"

2. In your new directory, run: "cmake -DCMAKE_INSTALL_PREFIX=${PWD} ../NobleElementSimulationTechnique" (or whatever the path to your source directory is)

3. make; make install


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

Then, one can edit the parameters/functions in MyDetector.hh as desired.
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


Analysis Settings:

analysis.hh -- for specifying output style/analysis parameters
	- MCtruthE: whether true energy or "reconstructed"/statistically measured energy is expressed
	- MCtruthPos: whether true position or "reconstructed"/statistically fluctuated positions are expressed
	- useTiming: photon arrival times + pulse shapes used
	- usePE: whether pulse areas are given in pe, phd, or spike count.
	- useS2: whether bands are expressed in log(S2/S1) or log(S2)
	- min/maxS1, numBins: controls S1 analysis threshold cut and binning
	- min/maxS2: S2 analysis threshold cut
	- z_step: mm steps for integrating non-uniform field
	- E_step: keV steps for integrating WIMP spectrum

NOTE: Make sure to recompile after changing analysis settings!

How to calculate g2 from the detector parameters:
  ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645 = -0.03754*pow(E_gas/1.848,2.)+0.52660*(E_gas/1.848)-0.84645
  SE_size = g1_gas*gasGap_mm*0.1*(0.137*E_gas*1000-177*p_bar-45.7)
  g2 = ExtEff*SE_size
where E_gas, g1_gas, gasGap_mm, and p_bar are provided by the settings file


Run Instructions:

This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector:
	./testNEST numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

For 8B, numEvts (integer) is replaced with kg-days of exposure (floating-point #) with all other input parameters the same (unchanged). For WIMPs:
	./testNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

For cosmic-ray muons or other similar particles with elongated track lengths:
	./testNEST numEvts {MIP} LET[MeV*cm^2/gram] step_size[cm] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

NOTES:
	- If you want to use the default drift field from detector settings, put -1 as the argument for field_drift.
	- If you want to randomly distribute events in space, rather than specify a point source, put -1 as the argument for x,y,z.

	- If you want to use the ROOT tools provided (and have ROOT installed already on your machine of course!) then to compile
	g++ -g -Wno-deprecated-declarations -Ofast -o rootNEST `root-config --cflags --libs` rootNEST.cpp
	Wno-... flag is optional, gets rid of annoying warning that's ROOT's fault not you. Ofast optional too, for speed. Might be -O3 on your machine.
