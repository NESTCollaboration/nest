Installation instructions:

1. Make an empty directory somewhere outside the NEST source code directory, and go there.
   For example: (from the source code directory) "mkdir ../nest-build; cd ../nest-build"

2. In your new directory, run: "cmake -DCMAKE_INSTALL_PREFIX=${PWD} ../NobleElementSimulationTechnique" (or whatever the path to your source directory is)

3. make; make install

Run instructions, as of May 4, 2018:
 
./testNEST
This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector.

numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
for 8B or WIMPs, numEvts is kg-days of exposure

exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
for cosmic-ray muons or other similar particles with elongated track lengths...

numEvts {MIP} LET[MeV*cm^2/gram] step_size[cm] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
