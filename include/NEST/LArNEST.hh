/**
 * @file LArNEST.hh
 * @author NEST Collaboration
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Justin Mueller [Justin.Mueller@colostate.edu]
 * @author Michael Mooney [mrmooney@colostate.edu]
 * @brief 
 * @version 
 * @date 2022-04-13
 */
#pragma once

#include "NEST.hh"

// #define HIGH_E_NR 330.
// #define W_DEFAULT \
//   13.4  // default work func, in eV. arXiv:1611.10322. +/- 0.35. 19.5-19.6 eV
//         // for LAr
// #define W_SCINT \
//   8.5e-3  // the *max* possible energy of 1 scint phot, keV. Make this at least
//           // 10 eV for LAr
// #define NEST_AVO 6.0221409e+23
// #define ATOM_NUM \
//   18.  // period to make float. 18 for LAr. If changed here go to TestSpectra.hh
//        // too

// #define PHE_MIN 1e-6  // area
// #define FIELD_MIN 1.  // min elec field to make S2 (in V/cm)
// #define DENSITY 2.90  // g/cm^3, ref density for dependent effects. ~1.4 for LAr

// #define EPS_GAS \
//   1.00126  // poly-morphic: make negative to use Aprile/PandaX instead of
//            // LLNL/PIXeY's e- EE
// // for GAr it is 1.000574 at least at room T (doi.org/10.1103/PhysRev.34.615)
// #define EPS_LIQ \
//   1.85  // LXe dielectric constant explicitly NOT 1.96 (old). Update thx to Dan
//         // M. LAr 1.325

// #define SAMPLE_SIZE 10  // nano-seconds
// #define PULSE_WIDTH 10  // nano-seconds
// #define PULSEHEIGHT \
//   0.005                  // threshold height, in PE, for writing to photon_times
// #define SPIKES_MAXM 120  // above this switch to pulse area (70 phd in 1 array)
// #define PHE_MAX 180      // saturation threshold, in PE per bin i.e. sample

namespace larnest
{
    /**
     * @brief 
     * 
     */
    class LArNEST : public NEST::NESTcalc
    {
    public:
        explicit LArNEST(VDetector *detector);
    private:
        
    };
}