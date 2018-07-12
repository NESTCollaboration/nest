
 The following tools are provided in this sub-directory: bareNEST and rootNEST

 bareNEST is a simplified example based on testNEST, which strips all complications like muons and field fringing.
 It provides the basics showing you how to implement the NEST class in your own code in a no-frills approach, just to get the simple yields and quanta.

 rootNEST is an extremely diverse and powerful tool, but requires compilation against ROOT libraries to work at all.
 Instructions are in the setup readme for compiling. What it does:

 * In "default" mode (i.e., no #define pre-compiler directives commented in at the top) it executes Gaussian fits to a histogrammed band in log(S2/S1) or log(S2) for you vs. S1
  - If you give it 2 arguments it'll give you the leakage of ER events below a smooth fit to the Gaussian means of the NR band, so it shows your background discrimination for WIMPs
  - If the second argument is NR and first ER raw data (produced first with testNEST) instead of the other way around then it will spit out NR "acceptance" below the centroid

 * In FIT mode (#define and #ifdef FIT) it takes 2 arguments and compares the NEST band to real data and gives you the goodness of fit

 * In either FIT or normal mode, if the number of bins in analysis.hh is set to 1, then it assumes you want a Gaussian fit (S1, S2, and energy) for a mono-energy calibration peak

 * in LIMIT mode (#define and #ifdef LIMIT) it takes 2 arguments
  - one is a file of NR efficiency vs. energy, which can be provided by real data or by running testNEST or similar code based on it repeatedly and building a new text file
  - code will perform a smooth fit for you to the efficiency
  - the second argument is 0, 1, or 2 for spin-independent or spin-dependent neutron or proton respectively. It will ask you questions like #kg-days on screen

  The following example tab-delimited plain ASCII text files are provided to go with rootNEST above.

 + Xe10_ERBand_Luiz.txt and Xe10_NRBand_Luiz.txt graciously produced by Prof. Luiz de Viveiros of Penn State for XENON10 for use in FIT mode to compare to NEST MC output
 + LUX_Run03.txt comes from Phys. Rev. Lett. 116, 161301 (2016) Fig. 1 thick black band, used without uncertainty here. Needed for LIMIT mode. Do NOT treat as "official" LUX #'s

 If you have questions, please contact Prof. Matthew Szydagis, mszydagis@albany.edu
