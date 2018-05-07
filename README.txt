Installation Instructions:

1. Make an empty directory somewhere outside the NEST source code directory, and go there.
   For example: (from the source code directory) "mkdir ../nest-build; cd ../nest-build"

2. In your new directory, run: "cmake -DCMAKE_INSTALL_PREFIX=${PWD} ../NobleElementSimulationTechnique" (or whatever the path to your source directory is)

3. make; make install


Settings Files:

Two header files are needed in the main directory to be included in the NEST code:
  1) detector.hh -- for specifying detector parameters
  2) analysis.hh -- for specifying output style/analysis parameters
		- MCtruthE: whether true energy or "reconstructed"/statistically measured energy is expressed
		- MCtruthPos: whether true position or "reconstructed"/statistically fluctuated positions are expressed
		- useTiming: photon arrival times + pulse shapes used
		- usePE: whether pulse areas are given in pe, phd, or spike count.
		- useS2: whether bands are expressed in log(S2/S1) or log(S2)
		- min/maxS1, numBins: controls S1 analysis threshold cut and binning
		- min/maxS2: S2 analysis threshold cut
		- z_step: mm steps for integrating non-uniform field
		- E_step: keV steps for integrating WIMP spectrum

NOTE: Make sure to recompile after changing analysis or detector settings!

How to calculate g2 from the detector parameters:
  ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645 = -0.03754*pow(E_gas/1.848,2.)+0.52660*(E_gas/1.848)-0.84645
  SE_size = g1_gas*gasGap_mm*0.1*(0.137*E_gas*1000-177*p_bar-45.7)
  g2 = ExtEff*SE_size
where E_gas, g1_gas, gasGap_mm, and p_bar are provided by the settings file


Run Instructions:

This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector:
	./testNEST numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

For 8B or WIMPs, numEvts is kg-days of exposure:
	./testNEST exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}
 
For cosmic-ray muons or other similar particles with elongated track lengths:
	./testNEST numEvts {MIP} LET[MeV*cm^2/gram] step_size[cm] field_drift[V/cm] x,y,z-position[mm] {optional:seed}

NOTES:
	- If you want to use the default drift field from detector settings, put -1 as the argument for field_drift.
	- If you want to randomly distribute events in space, rather than specify a point source, put -1 as the argument for x,y,z.


