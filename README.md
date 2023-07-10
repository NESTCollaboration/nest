# nest

Noble Element Simulation Technique (nest) is used to simulate noble-element energy deposition microphysics.

[![Build Status](https://www.travis-ci.com/NESTCollaboration/nest.svg?branch=master)](https://www.travis-ci.com/NESTCollaboration/nest)
[![DOI](https://zenodo.org/badge/96344242.svg)](https://zenodo.org/badge/latestdoi/96344242)

(For Python bindings, see the related [nestpy project](https://github.com/NESTCollaboration/nestpy))

(For Information Regarding a Specific NEST Release, please see the [NEST Release Notes](https://github.com/NESTCollaboration/nest/releases)

### Table of Contents

1. [ Getting the Repository ](#get)
2. [ Building NEST ](#build)
	* [ Prerequisites ](#prereq)
	* [ Running CMake ](#config)
	* [ Compiling NEST ](#compile)
	* [ Linking NEST in a Separate CMake Project ](#linknest)
3. [ Simulation Settings ](#settings)
	* [ Creating a Custom Detector ](#detector)
	* [ Modifying the execNEST Output (Analysis Settings) ](#analysis)
	* [Using Alternate NR and ER Yields Models](#yieldModels)
4. [ Running NEST ](#run)
	* [ Running in Different Modes ](#modes)
	* [ Input Arguments ](#inputs)
	* [ Other User-Modifiable Parameters ](#usermod)
	* [Running with Pulse Timing](#timing)
	* [ Running loopNEST ](#loop)
5. [ Useful Tools ](#tools)
	* [ Using bareNEST ](#barenest)
	* [ Running rootNEST ](#rootnest)
	* [Custom Energy Spectrum Example: 220RnCalib.sh](#220rn)
	* [Generating S1s and S2s for Multiple Scatters](#multiple)
6. [ GEANT4 Integration ](#geant)
7. [ Need Help with Detector Parameters? ](#params)
8. [ Versioning ](#versions)
9. [ Contact (Authors) ](#contact)
10. [ Citation ](#citation)
11. [ License ](#license)
12. [ Python ](#python)
13. [ Garfield++ Integration ](#garfield)


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

* CMake 3.11 (or higher: tested with 3.17.5 on Linux CentOS 7 and 3.21.3 on Mac OSX 10.14.6 and 3.25.1 on Linux Mint 19.3 Cinnamon)
* gcc (with C++17 support, with GCC 8 or higher; tested with 8.3.1-9.4.0 on Linux and Apple LLVM version 10.0.1, clang-1001.0.46.4)

The following are optional, depending on intended use of NEST:

* GEANT4 (tested as high as Geant4.10.7.p01)
* ROOT (version 6.22/04 or higher, as it breaks on 6.22/03, ever since C++17 switchover: tested on 6.24/04 most recently on Mac OSX, 6.24/06 on CentOS 7)

<a name="config"></a>
### Running CMake

1. Make empty build and install directories somewhere outside the NEST source code directory, 
	and go to the build directory. 

	For example, from the source code directory: 

	```
	mkdir ../build; mkdir ../install; cd ../build
	```

2. From within the build directory, configure CMake. On systems where the C compiler/C++ compiler 
	paths are already specified, one can simply do (sometimes not cmake but cmake3 e.g. on Kitware on CentOS 7):

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

To generate the NEST library, execNEST executable (***THE MOST IMPORTANT, SWISS-ARMY KNIFE TOOLBOX***), and various tools (bareNEST, rootNEST): 

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
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 11)

find_package(NEST REQUIRED)
include_directories(${NEST_INCLUDE_DIRS})

add_executable(MyApp MyApp.cc)
target_link_libraries(MyApp NEST::NESTCore)
```

Then, include the following flag when running CMake on your project:

```
-DCMAKE_PREFIX_PATH=[path to NEST install directory]
```


<a name="settings"></a>
## Simulation Settings

<a name="detector"></a>
### Creating a Custom Detector

The "include/Detectors/" folder contains a number of header files:

* VDetector.hh
* DetectorExample_XENON10.hh
* LUX_Run03.hh (defualt)

The VDetector files (.cpp, .hh) which serve as the base detector class **should not be edited** 
(except by NEST developers). VDetector is a virtual class inherited by your custom detector.

"DetectorExample_XENON10.hh" and "LUX_Run03.hh" are examples of a custom detector file, and the latter is currently 
set as the default in execNEST, the executable which basically does every analysis for you all in one place.
It serves as a template for users to create their own detector. Follow the steps below (where "MyDetector" is an example).

1. Within the "include/Detectors/" folder, create your own header using the template:

	```	
	cp LUX_Run03.hh MyDetector.hh
	```

2. Edit the parameters/functions in MyDetector.hh as desired. 
	At minimum, one has to replace every instance of "DetectorExample_RUN03" with "MyDetector": 

	```cpp	
	// DetectorExample_RUN03.hh 
	
	--> // MyDetector.hh
	```

	```cpp
	#ifndef DetectorExample_LUX_RUN03_hh 
	
	--> #ifndef MyDetector_hh
	```

	Note that you can add your own custom functions (within constraints) to this detector class, 
	which can be accessed through the detector object within execNEST.

3. Once your header file is done, one should edit "src/execNEST.cpp" in the following way:

	```cpp
	#include "LUX_Run03.hh" 
	
	--> #include "MyDetector.hh"
	```

	```cpp
	DetectorExample_LUX_RUN03 detector = new DetectorExample_LUX_RUN03(); 
	
	--> MyDetector* detector = new MyDetector();
	```

	You've now included your custom header file, and created a "MyDetector" object for manipulation.	

4. In execNEST.cpp (or in whatever implementation links the NEST class), one can access any detector parameter
	using the corresponding "get" functions and reset parameters using their "set" functions:

	```cpp
	// Grab detector parameter "g1"
	double temp_g1 = detector->get_g1();

	// Modify detector parameter "g1" internally
	detector->set_g1(0.2);
	```

WARNING: Whenever the detector headers or execNEST is modified, **make sure to do a clean recompile**!

<a name="analysis"></a>
### Modifying the execNEST Output (Analysis Settings)

The file "include/NEST/analysis.hh" is for specifying the output style and setting various analysis parameters.

An explanation of the various parameters:

* verbosity: if set to **-1**, execNEST **only** outputs the columns of yield data
* MCtruthE: whether true energy or "reconstructed"/statistically measured energy is expressed
* MCtruthPos: whether true position or "reconstructed"/statistically fluctuated positions are expressed
* s1CalculationMode: set the S1 calculation mode options are: 
  * `Full` [Default]: calculating the pulse area by looping over all the pmt hits.
  * `Parametric`: calculating the pulse area by using a parametric equation
  * `Hybrid`: Using Full and Parametric with a transition point: n_pmt_hits > n_pmts
  * `Waveform`: calculating the pulse area with the Full calculation mode and the waveform
* s2CalculationMode: S2 calculation mode options are:
  * `Full` [Default]: calculate only the pulse area
  * `Waveform`: calculate the pulse area and the waveform
  * `WaveformWithEtrain`: calculate the pulse area and the waveform with etrain
* usePD: whether pulse areas are given in pe, phd, or spike count.
* useS2: whether bands are expressed in log(S2/S1) or log(S2)
* min/maxS1, numBins: controls S1 analysis threshold cut and binning
* min/maxS2: S2 analysis threshold cut
* logMin/logMax: min and max log(S2/S1) or log(S2) admitted into analysis including limit
* z_step: mm steps for integrating non-uniform field
* E_step: keV steps for integrating WIMP spectrum

WARNING: Whenever you modify this header, **make sure to do a clean recompile**!

Only **verbosity**, **MCtruthE**, and **MCtruthPos** will change the event-by-event output for execNEST. 
The min/max S1,S2,log values will only effect post-simulation analyses, such as event selection for energy 
reconstruction and efficiency calculations in execNEST, and binning for creating bands and leakage calculations in rootNEST.

**Note**: There can be negative S1 and S2 pulse areas in the execNEST output. Negative pulse areas simply flag an event that is below thresholds set in the Detector settings file. 
For S1s, negative areas indicate that the PMT coincidence hasn't been met. For a PMT hit to count towards the coincidence requirement, the pulse area in that channel must be greater than the detector->sPEthr (single photon threshold) parameter. 
For S2s, the S2 area must be larger than the S2 threshold (generally O(100 phd)). 
Although these events are below threshold, they're still present in a detector, therefore they are not removed from the output.
For the actual pulse area of these events, take the magnitude of the S1 or S2. 

<a name="yieldModels"></a>
## Using Alternate NR and ER Yields Models

As of NESTv2.1.0, additional yields models have been provided for nuclear recoils and beta electronic recoils. 
Since NESTv2.0.0, the NR yields model has changed. To use the original model see lines 660-661 in src/NEST.cpp:

```
   return GetYieldNR(energy, density, dfield, massNum,NuisParam);
   //return GetYieldNROld ( energy, 1 );
```
By commenting out the first of those lines, and un-commenting the second, you will have the original NR model.

For the ER beta model, extensive work was done in arXiv:1910.04211 to create a LUX-specific beta model. 
To use this model, see lines 674-675 of src/NEST.cpp:

```
   return GetYieldBeta(energy,density,dfield);
   //return GetYieldBetaGR(energy,density,dfield,ERYieldsParam);
```
By commenting out the first of those lines, and un-commenting the second, you will have the LUX-specific yield model.


<a name="run"></a>
## Running execNEST

The execNEST executable is a comprehensive, powerful tool that will implement the NEST class in a 
pre-structured way. It loads pre-defined detector settings (including electric field, optics, etc. from the 
"VDetector" inheriting class) and source spectra (from the "TestSpectra" class) to do fast detector simulations.

Running execNEST without arguments will remind you of the necessary inputs:

```
./execNEST
```

<a name="modes"></a>
### Running in Different Modes

This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector:

```
./execNEST numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
```
To simulate time-dependent 83m-Kr decays -- 83m-Kr produces yields via a 32.1 keV $\gamma$ followed by a 9.4 keV gamma --  E_max[keV] is replaced with the time between the decays in ns; E_min[keV] is replaced with either 9.4, 32.1, or 41.5 [keV]. Example of 9.4 keV decay 250ns (NOTE: this is the MAX time separation along an exponential) after the inital 32.1 keV follows below. The minimum deltaT (minTimeSeparation) defaults to 100 (so, merged pulses are included) but you can change that at the top of execNEST.cpp (make =MAX to fix t)

```
./execNEST numEvts Kr83m 9.4 250 field_drift[V/cm] x,y,z-position[mm] {optional:seed}
```

To simulate 8B, numEvts (integer) is replaced with kg-days of exposure (floating-point #) with all other input parameters the same (unchanged). 

To simulate WIMPs:

```
./execNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
```

To simulate cosmic-ray muons or other similar particles with elongated track lengths:

```
./execNEST numEvts {MIP} LET[MeV*cm^2/gram] x,y-position[mm](Initial) field_drift[V/cm] x,y,z-position[mm](Final) {optional:seed}
```

<a name="inputs"></a>
### Input Arguments

* **field_drift**: if set to -1, the electric field map from the detector settings is used.
	If set to any other value, that drift field (in V/cm) will be used throughout the detector.
	The latter method makes the code much faster, so always use it if the field is uniform.
* **x,y,z-position**: if set to -1, each event's position is drawn randomly within the detector.
	Otherwise, the **exact** syntax "x,y,z" (with commas) produces events at that coordinate.
	"-999,-999,Z" will lead to Random X, Random Y, but fixed Z.
* **seed**: if set to -1, the random seed is internally set based on the timestamp.
  	If another number is provided (optional), that number is the random seed (for debugging purposes).
	If no seed is specified, then the default random seed is always 0 (for reproducibility purposes).

<a name="usermod"></a>
### Other User-Modifiable Parameters

* **NuisParam**: an array currently located in line 43 of execNEST.cpp.
	These change the mean light and charge yields of nuclear recoils (separately) as E-independent multiplicative factors.

<a name="timing"></a>
### Running with Pulse Timing

In include/NEST/analysis.hh, the `s1CalculationMode` and `s2CalculationMode` flags allow users to get S1 and S2 photon arrival times for each simulated event.
`s1CalculationMode=Full` and `s2CalculationMode=Full` (default) will not return timing info. If verbosity=1 and `s1CalculationMode=Waveform`/`s2CalculationMode=Waveform`, running execNEST will create a file called 
"photon_times.txt" in the user's build directory with S1 and S2 top/bottom PMT arrival times for each event. If `s2CalculationMode=WaveformWithEtrain`, 
approximated e-trains will be included for the output events.

The script, "pulseShape.cpp", in the "examples" directory will take the photon_times.txt output, and calculate pulse shape parameters
for each event in the file. This must be compiled separately from the other NEST executables. 

<a name="loop"></a>
### Running loopNEST

New as of v2.1.0, users have the ability to transform execNEST into _loopNEST_. This is declared via the **loopNEST** flag
in include/NEST/analysis.hh . **Using loopNEST will change the functionality of all execNEST inputs** and therefor should
only be used by experienced NEST users. The loopNEST event parameters are hard-coded into lines 70-140 of src/execNEST.cpp. 

loopNEST is a powerful tool that allows NEST to be iterated over many times, with different model parameters or detector parameters.
In essence, the user can vary select NEST parameters over a loop of values, and in tandem with FIT mode of rootNEST (see below), get a best fit model to their provided data set. See arXiv:1910.04211 for an example in the literature where this is useful. 

A bash script, "loopNEST.sh" is provided in the examples directory. It serves as an example of how one would use loopNEST to 
find best-fit model parameters for a data set.

<a name="tools"></a>
## Useful Tools and Examples

The "examples/" folder contains a few very useful codes, and most importantly, it is where the rootNEST.cpp code lives. 
Additionally, several other scripts and examples for performing advanced analyses are included. 

<a name="barenest"></a>
### Using bareNEST

**NOTE**: bareNEST is deprecated as of NEST v2.1.0. It will eventually be phased-out from the standard NEST compilation.

bareNEST is a minimal implementation of NEST, which does not contain the various complications
of execNEST (e.g. field non-uniformities, muon tracks, etc.). It provides the "bare minimum" function calls
needed to demonstrate the yield calculations step-by-step. 

It is currently hard-coded to take no inputs, as this source code is a skeleton intended for 
adaptation/customization.

<a name="rootnest"></a>
### Running rootNEST

rootNEST is an extremely diverse and powerful tool, but requires compilation against ROOT libraries to work at all.
**NEW TO v2.1.0**: the rootNEST mode is set in include/NEST/analysis.hh. The **mode** flag in analysis.hh has 3 options:


* **mode = 0**: "default" mode. It executes Gaussian fits to a histogrammed band in log(S2/S1) or log(S2) for you vs. S1. 
	- If you give it only one argument (being the file path to a execNEST output file) with **verbosity=1**, it will 
	        print out S1 vs. log(S2/S1) or S1 vs. log(S2) band information for the given S1,S2,log ranges in analysis.hh .
	- If you give it 2 arguments it'll give you the leakage of ER events below a smooth fit to the Gaussian means 
		of the NR band, so it shows your background discrimination for WIMPs. 
	- If the second argument is NR and first ER raw data (produced first with execNEST) instead of the other way 
		around then it will spit out NR "acceptance" below the centroid.

* **mode = 1**: "FIT" mode. It takes 2 arguments and compares the NEST band to real data and 
	gives you the goodness of fit. Note that the data file must be formatted to match rootNEST band outputs. 
	See "LUXRun03_CH3TBand_SpikeFull.txt" and "LUXRun03_DDBand_SpikeFull.txt" for example format.

* In either FIT or default mode, if the number of bins in analysis.hh is set to 1, then it assumes you want 
	a Gaussian fit (S1, S2, and energy) for a mono-energy calibration peak

* **mode = 2**: "LIMIT" mode. It takes 2 arguments:
	1. A file of NR efficiency vs. energy, which can be provided by real data or by running execNEST or 
		similar code based on it repeatedly and building a new text file. The code will perform a smooth fit 
		for you to the efficiency. See "LUXRun03_betaEff_Simulated.txt" or "LUXRun03_nrEff_Simulated.txt" for 
		formatting examples. 
	2. 0, 1, or 2 for spin-independent or spin-dependent neutron or proton respectively. 
		It will ask you questions like #kg-days on screen

There is a new **skewness** parameter in analysis.hh that allows rootNEST to fit skew-Gaussian to the log(S2/S1) 
or log(S2) distributions for bands and goodness-of-fit calculation. skewness=0 will perform Gaussian fits;
skewness=1 will perform skew-Gaussian fits; skewness=2 will perform more rigorous skew-Gaussian fits 
for additional accuracy. 

The following additional example tab-delimited plain ASCII text files are provided to go with rootNEST.

* Xe10_ERBand_Luiz.txt and Xe10_NRBand_Luiz.txt graciously produced by Prof. Luiz de Viveiros of Penn State 
	for XENON10 for use in FIT mode to compare to NEST MC output.
	
	This XENON10 data comes in the following format:
	
	```
	Bin Center  Bin Actual  Gaus Mean  Mean Error  Gaus Sigma  Sig Error
	```

NOTE: While rootNEST can be built using NEST's built-in CMake (see setup instructions), you can compile 
yourself by doing:

```
g++ -g -Wno-c++11-extensions -Wno-deprecated-declarations -Ofast `root-config --cflags` -o rootNEST rootNEST.cpp `root-config --libs`
```

The argument "-Ofast" may be "-O3" on your machine (optimization flag). 


<a name="220rn"></a>
## Custom Energy Spectrum Example: 220RnCalib.sh

While TestSpectra.cpp provides many useful energy spectra to be used directly in execNEST, it is often useful to use a custom
energy spectrum, or one from a physics signal that is not included in TestSpectra.cpp. The script, "220RnCalib.sh",
in the examples directory provides one method of creating a custom spectrum. 

This example uses the 212Pb ER spectrum ( from the 220Rn chain: a prevalent background in xenon TPCs ) and writes the execNEST
output to a single file: "Pb212.dat". This script uses the 212Pb energy spectrum, and calls execNEST for a given energy bin.
The number of events in each execNEST command is the probability for a 212Pb event to recoil in that energy bin, multiplied 
by the total desired number of output events (1M in this case). 

For any general spectrum, the user will want to run execNEST for every relevant bin of some energy spectrum PDF:

```
./execNEST numTotalEvts*PDF_value interactionType bin_min[keV] bin_max[keV] drift_field[V/cm] x,y,z_pos[mm] {optional:seed} >> outputfile 
```

Note the funnel command used here, ```>>``` appends the command's output to an existing file (if the file already exists). 
One must make sure they use proper ```rm outputfile``` commands when trying to re-run. See the "220RnCalib.sh" example.

<a name="multiple"></a>
## Generating S1s and S2s for 2+ Scatters

Included in the examples directory, the multipleScatter.cpp script produces the ```multipleScatter``` executable.
This tool uses one or two execNEST output files to stitch together energy depositions into multiple scatter events, and will return either corrected or uncorrected pulse areas.
Running the executable without arguments prints the usage:

    ./path/to/multipleScatter <nEvents> <nScatters> <Pulse Type> <File 1> <optional: File 2>

    -- <nEvents> is the number of scatters to stitch together per event.
    -- <nScatters> is the number of depositions to use per event; using -1 gives a random exponential draw for each event.
    -- <Pulse Type> is either 0 (for uncorrected pulses) or 1 (for corrected pulses); corrected pulses sum S2s and take an S2-weighted average S1c.
    -- <File 1> is execNEST output events used to stitch together energy deposits.
    -- <File 2> is optional, and will be the last scatter in each event, if provided.

The exponential draw used when running with ```nScatters = -1``` is a 0th order attempt at providing similar behavior to multiple scatters in LXe. 
The exponential uses a half-life of 1 scatter, between 2 and 50 scatters. These can be changed by the user in the script for max flexibility. 

Using two files can be useful to model certain types of interactions. For example, using a beta ER file with gamma ER events will mimic a Compton scatter followed by photoabsorption. Mixing ER and NR can be used for inelastic scattering of n's or WIMPs, or the Migdal effect.
The executable will print S1 and S2 areas to screen in addition to the number of scatters used in the event, which can be funneled into a text document.

<a name="geant"></a>
## GEANT4 Integration

1. Make sure you have built NEST with the proper G4 option in the CMake stage (see "Running CMake").

2. Include the following lines in the CMakeLists.txt of your GEANT4-based simulation project:

	```
	find_package(NEST REQUIRED)
	include_directories(${NEST_INCLUDE_DIRS})
    ...
    target_link_libraries(... NEST::NESTCore NEST::NESTG4)
	```

3. In your simulation's build directory, use whatever cmake command you usually do to build the simulation, but 
	include an additional flag:

	```
	-DCMAKE_PREFIX_PATH=[path to NEST install directory]
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
	Put this line after all other lines in Build(), else you may only ever run precisely one event.

5. In your physics list, you'll want code that looks like this in your G4VPhysicsConstructer::ConstructProcess() method:

	```cpp
	NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic,[your VDetector]);
	if (theNEST2ScintillationProcess->IsApplicable(*particle)) {
	   pmanager->AddProcess(theScintillationProcess, ordDefault + 1, ordInActive, ordDefault + 1);
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
to the execNEST default executable example. If you follow all of these instructions, you should be getting at 
least the 0th-order answer. Temperature, pressure, etc. should all be 1st- or 2nd-order.

OptTrans and SinglePEWaveForm involve quick and dirty ray tracing and pulse shaping, to be used with 
"s1CalculationMode/s2CalculationMode = Waveform" in analysis.hh.

FitS1, FitS2, and FitEF are for 3-D X,Y,Z position dependencies including field fringing. The first two can 
just return 1 (get corrected out anyhow).

NEST v2 is now predictive of things like SE size and extraction efficiency so that you don't have to put in manually.
But instead put in the more fundamental parameters listed above, and check the output to see if it is what you 
want (g2 is broken down: SE x e- ext eff).

For these functions and for all the #'s if you're confused or don't know something:

Just stick with the LUX or XENON10 defaults. Why? They are there because even though they're old detectors now,
their numbers are ~ representative of detectors past/present/future, big/small even.

There is an additional boolean detector parameter called **OldW13eV**. In order to incorporate the new measurement
of the work function in LXe by EXO-200 (arXiv:1908.04128), the OldW13eV flag (if false) will boost yields
to reconcile the discrepancy. This does not serve to solve the discrepancy, but provides NEST the flexibility to match
this new result.

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




<a name="garfield"></a>
## Garfield++ Integration

NEST can produce tables of transport parameters for liquid xenon for use with the Garfield++ package. These tables currently only contain the drift velocity and lateral/longitudinal diffusion information as a function of electric field. They can be used as .gas files in Garfield++ for use with the AvalancheMC process, and are intended to allow a user to make Garfield drift electrons through a "gas" that has the same transport properties as LXe. 

The script GenerateGarfieldGasTableForLiquidNoble.cpp (and its executable GenerateGasTable) will take inputs and produce a gas file for use with Garfield. To execute the script, build NEST, then execute the command

`./GarfieldppIntegration/GenerateGasTable Xe <nFieldPoints> <minField_VperCm> <maxField_VperCm> <Logarithmic:0or1> <Temperature_K>`

from the build directory. Generating tables for additional noble elements is not currently supported, but may be added in the future under popular demand.

Another note: for this to work, the temperature you use in your garfield .cpp file must be 293K, and the pressure must be 1350 torr. This gas file generation assumes these numbers and corrects for them so that the transport coefficients end up correct when the table is passed into Garfield++. If you use anything other than 1350 torr, your transport properties (diffusion, drift velocity, etc.) will be wrong!!!

For more information, contact Ryan Linehan (rlinehan@stanford.edu)
