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

