# nest

Noble Element Simulation Technique (nest) is used to simulate noble-element energy deposition microphysics.

[![Build Status](https://travis-ci.com/NESTCollaboration/nest.svg?token=8i3psWNJAskpVjC6qe3w&branch=master)](https://travis-ci.com/NESTCollaboration/nest)
[![DOI](https://zenodo.org/badge/96344242.svg)](https://zenodo.org/badge/latestdoi/96344242)

(For Python bindings, see the related [nestpy project](https://github.com/NESTCollaboration/nestpy))

### Table of Contents

1. [ Getting the Repository ](#get)
2. [ Building NEST ](#build)
	* [ Prerequisites ](#prereq)
	* [ Running CMake ](#config)
	* [ Compiling NEST ](#compile)
	* [ Linking NEST in a Separate CMake Project ](#linknest)
3. [ Simulation Settings ](#settings)
	* [ Creating a Custom Detector ](#detector)
	* [ Modifying the testNEST Output (Analysis Settings) ](#analysis)
4. [ Running NEST ](#run)
	* [ Running in Different Modes ](#modes)
	* [ Input Arguments ](#inputs)
	* [ Other User-Modifiable Parameters ](#usermod)
5. [ Useful Tools ](#tools)
	* [ Using bareNEST ](#barenest)
	* [ Running rootNEST ](#rootnest)
6. [ GEANT4 Integration ](#geant)
7. [ Need Help with Detector Parameters? ](#params)
8. [ Versioning ](#versions)
9. [ Contact (Authors) ](#contact)
10. [ Citation ](#citation)
11. [ License ](#license)
12. [ Python ](#python)


<a name="get"></a>
## Getting the Repository

In the terminal, one can clone this repository by typing the command:

`git clone https://personal_username@github.com/NESTCollaboration/nest.git`

This uses the HTTPS protocol. For environments (e.g. computing clusters) where one has to use the SSH protocol:

`git clone git@github.com:NESTCollaboration/nest.git`

Anyone in the "NESTCollaboration" organization should be able to develop (push changes to the remote repository).

Please contact Matthew Szydagis about becoming involved in development before merging with the master branch. 


<a name="build"></a>
## Building NEST

<a name="prereq"></a>
### Prerequisites

In order to compile and run NEST, the following are required:

* CMake 2.8 (or higher)
* gcc (with C++11 support, 4.8.1 or higher)

The following are optional, depending on intended use of NEST:

* GEANT4
* ROOT (version 6 or higher, known CMake issues with ROOT 5).

<a name="config"></a>
### Running CMake

1. Make empty build and install directories somewhere outside the NEST source code directory, 
	and go to the build directory. 

	For example, from the source code directory: 

	```
	mkdir ../build; mkdir ../install; cd ../build
	```

2. From within the build directory, configure CMake. On systems where the C compiler/C++ compiler 
	paths are already specified, one can simply do:

	```
	cmake -DCMAKE_INSTALL_PREFIX=[path to install directory] ../relative/path/nest
	```

	To compile with GEANT4 integration, you will have to include an extra argument:

	```
	cmake -DG4=ON -DCMAKE_INSTALL_PREFIX=[path to install directory] ../relative/path/nest
	```

	To compile with ROOT compilation, you will have to include an extra argument:

	```
	cmake -DBUILD_ROOT=ON -DCMAKE_INSTALL_PREFIX=[path to install directory] ../relative/path/nest
	```
	
	If the compiler paths cannot be found (for example, on some cluster environments), one should do:

	```
	cmake -DCMAKE_C_COMPILER="/path/to/bin/gcc" -DCMAKE_CXX_COMPILER="/path/to/bin/g++" -DCMAKE_INSTALL_PREFIX=${PWD} ../relative/path/nest
	```

If CMake throws errors, it is likely that either your environment variables are not set correctly 
or you do not have the necessary prerequisites installed.

<a name="compile"></a>
### Compiling NEST

From within the build directory (where the CMake commands were executed), one can now compile NEST.

To generate the NEST library, testNEST executable, and various tools (bareNEST, rootNEST): 

```	
make; make install
```

The compiled executables and libraries, and all the headers, should now be present in the install 
directory ("../install" if next to the build directory). 

After changes to the source code, one should do a clean make for good measure:

```
make clean; make; make install
```

<a name="linknest"></a>
### Linking NEST in a Separate CMake Project

CMake configuration files are present in the install directory. 
You may use these to include NEST in cmake-configured projects.  

As an example, if you would like to link NEST in a project which generates an executable called "MyApp"
from the script "MyApp.cc", add the following lines to your CMakeLists.txt:

```
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 11)

find_package(NEST REQUIRED)
include_directories(${NEST_INCLUDE_DIRS})

add_executable(MyApp MyApp.cc)
target_link_libraries(MyApp NEST)
```

Then, include the following flag when running CMake on your project:

```
-DNEST_DIR=[path to NEST install directory]
```


<a name="settings"></a>
## Simulation Settings

<a name="detector"></a>
### Creating a Custom Detector

The "include/Detectors/" folder contains a number of header files:

* VDetector.hh
* DetectorExample_XENON10.hh

The VDetector files (.cpp, .hh) which serve as the base detector class **should not be edited** 
(except by NEST developers). VDetector is a virtual class inherited by your custom detector.

"DetectorExample_XENON10.hh" is an example of a custom detector file, and is currently set as the default in testNEST.
It serves as a template for users to create their own detector. Follow the steps below (where "MyDetector" is an example).

1. Within the "include/Detectors/" folder, create your own header using the template:

	```	
	cp DetectorExample_XENON10.hh MyDetector.hh
	```

2. Edit the parameters/functions in MyDetector.hh as desired. 
	At minimum, one has to replace every instance of "DetectorExample_XENON10" with "MyDetector": 

	```cpp	
	// DetectorExample_XENON10.hh 
	
	--> // MyDetector.hh
	```

	```cpp
	#ifndef DetectorExample_XENON10_hh 
	
	--> #ifndef MyDetector_hh
	```

	Note that you can add your own custom functions (within constraints) to this detector class, 
	which can be accessed through the detector object within testNEST.

3. Once your header file is done, one should edit "src/testNEST.cpp" in the following way:

	```cpp
	#include "DetectorExample_XENON10.hh" 
	
	--> #include "MyDetector.hh"
	```

	```cpp
	DetectorExample_XENON10* detector = new DetectorExample_XENON10(); 
	
	--> MyDetector* detector = new MyDetector();
	```

	You've now included your custom header file, and created a "MyDetector" object for manipulation.	

4. In testNEST.cpp (or in whatever implementation links the NEST class), one can access any detector parameter
	using the corresponding "get" functions and reset parameters using their "set" functions:

	```cpp
	// Grab detector parameter "g1"
	double temp_g1 = detector->get_g1();

	// Modify detector parameter "g1" internally
	detector->set_g1(0.2);
	```

WARNING: Whenever the detector headers or testNEST is modified, **make sure to do a clean recompile**!

<a name="analysis"></a>
### Modifying the testNEST Output (Analysis Settings)

The file "include/NEST/analysis.hh" is for specifying the output style and setting various analysis parameters.

An explanation of the various parameters:

* verbosity: if set to **false**, testNEST **only** outputs the columns of yield data
* MCtruthE: whether true energy or "reconstructed"/statistically measured energy is expressed
* MCtruthPos: whether true position or "reconstructed"/statistically fluctuated positions are expressed
* useTiming: photon arrival times + pulse shapes used
* usePD: whether pulse areas are given in pe, phd, or spike count.
* useS2: whether bands are expressed in log(S2/S1) or log(S2)
* min/maxS1, numBins: controls S1 analysis threshold cut and binning
* min/maxS2: S2 analysis threshold cut
* logMin/logMax: min and max log(S2/S1) or log(S2) admitted into analysis including limit
* z_step: mm steps for integrating non-uniform field
* E_step: keV steps for integrating WIMP spectrum

WARNING: Whenever you modify this header, **make sure to do a clean recompile**!


<a name="run"></a>
## Running testNEST

The testNEST executable is a comprehensive, powerful tool that will implement the NEST class in a 
pre-structured way. It loads pre-defined detector settings (including electric field, optics, etc. from the 
"VDetector" inheriting class) and source spectra (from the "TestSpectra" class) to do fast detector simulations.

Running testNEST without arguments will remind you of the necessary inputs:

```
./testNEST
```

<a name="modes"></a>
### Running in Different Modes

This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector:

```
./testNEST numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
```

To simulate 8B, numEvts (integer) is replaced with kg-days of exposure (floating-point #) with all other input parameters the same (unchanged). 

To simulate WIMPs:

```
./testNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
```

To simulate cosmic-ray muons or other similar particles with elongated track lengths:

```
./testNEST numEvts {MIP} LET[MeV*cm^2/gram] x,y-position[mm](Initial) field_drift[V/cm] x,y,z-position[mm](Final) {optional:seed}
```

<a name="inputs"></a>
### Input Arguments

* **field_drift**: if set to -1, the electric field map from the detector settings is used.
	If set to any other value, that drift field (in V/cm) will be used throughout the detector.
* **x,y,z-position**: if set to -1, each event's position is drawn randomly within the detector.
	Otherwise, the **exact** syntax "x,y,z" (with commas) produces events at this coordinate.
* **seed**: if set to -1, the random seed is internally set based on the timestamp.
	If another number is provided (optional), that number is the random seed (for debugging purposes).

<a name="usermod"></a>
### Other User-Modifiable Parameters

* **NuisParam**: an array currently located in line 43 of testNEST.cpp.
	These change the mean light and charge yields of nuclear recoils (separately) as E-independent multiplicative factors.


<a name="tools"></a>
## Useful Tools and Examples

The "examples/" folder contains two very useful codes: bareNEST and rootNEST.

<a name="barenest"></a>
### Using bareNEST

bareNEST is a minimal implementation of NEST, which does not contain the various complications
of testNEST (e.g. field non-uniformities, muon tracks, etc.). It provides the "bare minimum" function calls
needed to demonstrate the yield calculations step-by-step.

This tool provides the basics, showing you how to implement the NEST class in your own code in a no-frills 
approach.

NOTE: It is currently hard-coded to take no inputs, as this source code is a skeleton intended for 
adaptation/customization.

<a name="rootnest"></a>
### Running rootNEST

rootNEST is an extremely diverse and powerful tool, but requires compilation against ROOT libraries to work at all.

* In "default" mode (i.e., no #define pre-compiler directives commented in at the top) it executes Gaussian fits 
	to a histogrammed band in log(S2/S1) or log(S2) for you vs. S1. 
	- If you give it 2 arguments it'll give you the leakage of ER events below a smooth fit to the Gaussian means 
		of the NR band, so it shows your background discrimination for WIMPs. 
	- If the second argument is NR and first ER raw data (produced first with testNEST) instead of the other way 
		around then it will spit out NR "acceptance" below the centroid.

* In FIT mode (#define and #ifdef FIT) it takes 2 arguments and compares the NEST band to real data and 
	gives you the goodness of fit

* In either FIT or normal mode, if the number of bins in analysis.hh is set to 1, then it assumes you want 
	a Gaussian fit (S1, S2, and energy) for a mono-energy calibration peak

* In LIMIT mode (#define and #ifdef LIMIT) it takes 2 arguments:
	1. A file of NR efficiency vs. energy, which can be provided by real data or by running testNEST or 
		similar code based on it repeatedly and building a new text file. The code will perform a smooth fit 
		for you to the efficiency.
	2. 0, 1, or 2 for spin-independent or spin-dependent neutron or proton respectively. 
		It will ask you questions like #kg-days on screen

The following example tab-delimited plain ASCII text files are provided to go with rootNEST.

* Xe10_ERBand_Luiz.txt and Xe10_NRBand_Luiz.txt graciously produced by Prof. Luiz de Viveiros of Penn State 
	for XENON10 for use in FIT mode to compare to NEST MC output.
	
	This XENON10 data comes in the following format:
	
	```
	Bin Center  Bin Actual  Gaus Mean  Mean Error  Gaus Sigma  Sig Error
	```

* LUX_Run03.txt comes from Phys. Rev. Lett. 116, 161301 (2016) (or arXiv:1512.03506) Fig. 1 thick black band, 
	used without uncertainty here. Needed for LIMIT mode. Do NOT treat as "official" LUX numbers.

NOTE: While rootNEST can be built using NEST's built-in CMake (see setup instructions), you can compile 
yourself by doing:

```
g++ -g -Wno-c++11-extensions -Wno-deprecated-declarations -Ofast `root-config --cflags` -o rootNEST rootNEST.cpp `root-config --libs`
```

The argument "-Ofast" may be "-O3" on your machine (optimization flag). 


<a name="geant"></a>
## GEANT4 Integration

1. Make sure you have built NEST with the proper G4 option in the CMake stage (see "Running CMake").

2. Include the following lines in the CMakeLists.txt of your GEANT4-based simulation project:

	```
	find_package(NEST REQUIRED)
	include_directories(${NEST_INCLUDE_DIRS})
	```

3. In your simulation's build directory, use whatever cmake command you usually do to build the simulation, but 
	include an additional flag:

	```
	-DNEST_DIR=[path to NEST install directory]
	```

4. You must set NESTStackingAction to be your simulation's Stacking Action. 
	You probably already have a stacking action, which is set in a line that probably looks like:

	```cpp
	SetUserAction(new myStackingAction());
	```

	If you already have a stacking action, you must modify it so it inherits from NESTStackingAction:

	```cpp
	class myStackingAction : public NESTStackingAction
	```

	If you do not already have a stacking action, you still probably have class that inherits from 
	G4VUserActionInitialization. In this class's Build() method, do:

	```cpp
	SetUserAction(new NESTStackingAction());
	```

5. In your physics list, you'll want code that looks like this:

	```cpp
	NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic,[your electric field]);
	if (theNEST2ScintillationProcess->IsApplicable(*particle)) {
		pmanager->AddProcess(theScintillationProcess, ordDefault+1, ordInActive, ordDefault+1);
	}
	```

6. In your physics list's ConstructParticle() method, you must call: 
        ```cpp
        NEST::NESTThermalElectron::Definition();
        ```


<a name="params"></a>
## Need Help with Detector Parameters?

Be careful not to get too overwhelmed with all of the detector parameters. 
The most important ones to change are these to give you a first pass:

* g1
* g1_gas
* E_gas
* geometry parameters, especially Z dependences, but most importantly of all 
	the TopDrift and anode, since they set the size of the S2 gas gap

Reasons: the first sets S1 size of course, while the latter 3 control the S2 size. 
When running NEST: particle type, energy, field inputs are most important.

* Increases in either g1_gas and/or E_gas and/or the distance from "TopDrift" up to the anode increase the 
	SE (single e-) size
* However, increasing the E_gas also increases the extraction efficiency
* Use case: if you are trying to increase g2, but you know your gas gap and your SE size, you must increase 
	E_gas, but decrease g1_gas to compensate

How to calculate g2 from the detector parameters:

```cpp
ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645 = -0.03754*pow(E_gas/1.848,2.)+0.52660*(E_gas/1.848)-0.84645
SE_size = g1_gas*gasGap_mm*0.1*(0.137*E_gas*1000-177*p_bar-45.7)
g2 = ExtEff*SE_size
```

where "E_gas", "g1_gas", "gasGap_mm", and "p_bar" are in the detector settings file (e.g. "MyDetector.hh").

The next most important thing is the electric field of course, since that changes the recombination probability 
and thus the "slosh" between S1 & S2. You can set it with either FitEF, or as one of the run-time input arguments 
to the testNEST default executable example. If you follow all of these instructions, you should be getting at 
least the 0th-order answer. Temperature, pressure, etc. should all be 1st- or 2nd-order.

OptTrans and SinglePEWaveForm involve quick and dirty ray tracing and pulse shaping, to be used with 
"useTiming = true" in analysis.hh.

FitS1, FitS2, and FitEF are for 3-D X,Y,Z position dependencies including field fringing. The first two can 
just return 1 (get corrected out anyhow).

NEST v2 is now predictive of things like SE size and extraction efficiency so that you don't have to put in manually.
But instead put in the more fundamental parameters listed above, and check the output to see if it is what you 
want (g2 is broken down: SE x e- ext eff).

For these functions and for all the #'s if you're confused or don't know something:

Just stick with the XENON10 defaults, which are mixed with some published LUX values. Why? They are there because 
even though they're old detectors now, their numbers are ~ representative of detectors past/present/future, big/small even.


<a name="versions"></a>
## Versioning

For the versions available, see the [tags on this repository](https://github.com/NESTCollaboration/nest/tags). 


<a name="contact"></a>
## Contact (Authors)

If you have questions, please contact Prof. Matthew Szydagis, mszydagis@albany.edu.

See also the list of [contributors](https://github.com/orgs/NESTCollaboration/people) who participate in this project.


<a name="citation"></a>
## Citation

For now, please cite the version of the code that you are using by using the [Zenodo DOI](https://doi.org/10.5281/zenodo.2535713), where at the link you can get various bibliography styles including BibTeX.


<a name="license"></a>
## License

This code is licensed under the Open Source [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) License. It is available within the repository as the LICENSE file.

<a name="python"></a>
## Python

Python bindings are provided using pybind11 and distributed through the Python Package Index:

```
pip install nestpy
```

For more information, please see [the bindings project](https://github.com/NESTCollaboration/nestpy).
