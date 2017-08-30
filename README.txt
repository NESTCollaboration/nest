Installation instructions:

1. Make an empty directory somewhere outside the NEST source code directory, and go there. For example: (from the source code directory) "mkdir ../nest-build; cd ../nest-build"

2. In your new directory, run: "cmake -DCMAKE_INSTALL_PREFIX=${PWD} ../NobleElementSimulationTechnique" (or whatever the path to your source directory is)

3. make; make install

Run instructions:
 
./testNEST numEvts type_interaction E_min[keV] E_max[keV] density[g/cm^3] field_drift[V/cm]
