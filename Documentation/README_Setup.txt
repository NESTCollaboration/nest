Installation Instructions:

1. Make an empty directory either inside or outside the NEST source code directory, and go there.
	For example, from the source code directory: 

	>> mkdir build; cd build

2. In the new build directory, configure CMake. On systems where the C compiler/C++ compiler 
paths are already specified, one can simply do:

	>> cmake -DCMAKE_INSTALL_PREFIX=${PWD} ../relative/path/NobleElementSimulationTechnique

If the compiler paths cannot be found (for example, on some cluster environments), one should do:

	>> cmake -DCMAKE_C_COMPILER="/path/to/bin/gcc" -DCMAKE_CXX_COMPILER="/path/to/bin/g++" -DCMAKE_INSTALL_PREFIX=${PWD} ../relative/path/NobleElementSimulationTechnique

3. To build the source code into a library and generate the testNEST executable:
	
	>> make; make install

4. After changes to the source code, one should do a clean make for good measure:

	>> make clean; make; make install


Run Instructions:

This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector:
	./testNEST numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

For 8B, numEvts (integer) is replaced with kg-days of exposure (floating-point #) with all other input parameters the same (unchanged). For WIMPs:
	./testNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

For cosmic-ray muons or other similar particles with elongated track lengths:
	./testNEST numEvts {MIP} LET[MeV*cm^2/gram] x,y-position[mm](Initial) field_drift[V/cm] x,y,z-position[mm](Final) {optional:seed}

NOTES:
	- If you want to use the default drift field from detector settings, put -1 as the argument for field_drift.
	- If you want to randomly distribute events in space, rather than specify a point source, put -1 as the argument for x,y,z.

	- If you want to use the ROOT tools provided (and have ROOT installed already on your machine of course!) then to compile
	g++ -g -Wno-deprecated-declarations -Ofast -o rootNEST `root-config --cflags --libs` rootNEST.cpp
	Wno-... flag is optional, gets rid of annoying warning that's ROOT's fault not you. Ofast optional too, for speed. Might be -O3 on your machine.
	You can now also use #define FIT mode with rootNEST to fit a band. Example XENON10 band files courtesy of Luiz de Viveiros are included in the 6-column format:
	Bin Center  Bin Actual  Gaus Mean  Mean Error  Gaus Sigma  Sig Error  (in files Xe10_ERBand_Luiz.txt and Xe10_NRBand_Luiz.txt).
	With #define LIMIT, you can start calculating an exclusion curve, by first providing an efficiency curve (arXiv:1512.03506 provided as an example <- LUX Run03)
	
	- NuisParam is available to change the mean light and charge yields of nuclear recoils (separately) as E-independent multiplicative factors
