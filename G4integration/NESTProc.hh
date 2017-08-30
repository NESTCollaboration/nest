#ifndef NESTPROC_h
#define NESTPROC_h 1

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
//#include "G4ThermalElectron.hh"
#include "G4ProductionCuts.hh"
#include "G4MaterialCutsCouple.hh"
//#include "Analysis.hh"
#include "NEST.hh"
//#include "LUXSimManager.hh"

#define AVO 6.022e23 //Avogadro's number (#/mol)
#define EMASS 9.109e-31*kg
#define MillerDriftSpeed true

#define GASGAP 0.25*cm //S2 generation region
#define BORDER 0*cm //liquid-gas border z-coordinate

#define QE_EFF 1 //a base or maximum quantum efficiency
#define phe_per_e 1 //S2 gain for quick studies

// different field regions, for gamma-X studies
//#define WIN 0*mm //top Cu block (also, quartz window)
//#define TOP 0 //top grid wires
//#define ANE 0 //anode mesh
//#define SRF 0 //liquid-gas interface
//#define GAT 0 //gate grid
//#define CTH 0 //cathode grid
//#define BOT 0 //bottom PMT grid
//#define PMT 0 //bottom Cu block and PMTs

namespace NEST{
    struct Vertex {
        public:

            Vertex(double E, INTERACTION_TYPE species_in, std::array<double, 3>pos_in, double t_in, double density_in) : energy(E), species(species_in), time(t_in), pos(pos_in), size(0), density(density_in) { ; }
            Vertex(const Vertex& orig) : energy(orig.energy), species(orig.species), time(orig.time), pos(orig.pos), size(orig.size) { ; }

            static const Vertex merge(Vertex va, Vertex vb);

            std::array<double, 3> getPos() const {
                return pos;
            }

            INTERACTION_TYPE getSpecies() const {
                return species;
            }
            double gett() const{return time;}
            double getEnergy() const {
                return energy;
            }
            double getDensity() const {
                return density;
            }
        private:
            double energy; // keV
            INTERACTION_TYPE species;      // -1 for electron/positron/gamma, A for (neutral) nucleus
            double time;   // ns
            std::array<double, 3> pos; // mm
            double size;   // mm
            double density; // g/cc

        };
    std::vector<Vertex> cluster(std::vector<Vertex>);
    
    class NESTProc : public G4VRestDiscreteProcess {
    private:
        const G4int electron_cut_index;
    public: // constructor and destructor

        NESTProc(const G4String& processName = "S1",
                G4ProcessType type = fElectromagnetic);
        ~NESTProc();

    public: // methods, with descriptions
        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable', for any particle type except for an
        // 'opticalphoton' and for short-lived particles

        G4double GetMeanFreePath(const G4Track& aTrack,
                G4double,
                G4ForceCondition*);
        // Returns infinity; i. e. the process does not limit the step, but 
        // sets the 'StronglyForced' condition for the DoIt to be invoked at
        // every step.

        G4double GetMeanLifeTime(const G4Track& aTrack,
                G4ForceCondition*);
        // Returns infinity; i. e. the process does not limit the time, but
        // sets the 'StronglyForced' condition for the DoIt to be invoked at
        // every step.

        // For in-flight particles losing energy (or those stopped)
        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                const G4Step& aStep);
        G4VParticleChange* AtRestDoIt(const G4Track& aTrack,
                const G4Step& aStep);

        // These are the methods implementing the scintillation process.

        void SetTrackSecondariesFirst(const G4bool state);
        // If set, the primary particle tracking is interrupted and any
        // produced scintillation quanta are tracked next. When all have been
        // tracked, the tracking of the primary resumes.

        G4bool GetTrackSecondariesFirst() const;
        // Returns the boolean flag for tracking secondaries first.

        void SetScintillationYieldFactor(const G4double yieldfactor);
        // Called to set the scintillation quantum yield factor, useful for
        // shutting off scintillation entirely, or for producing a universal 
        // re-scaling to for example represent detector effects. Internally is
        // used for Lindhard yield factor for NR. Default should be user-set
        // to be 1 (for ER) in your simulation -- see NEST readme

        G4double GetScintillationYieldFactor() const;
        // Returns the quantum (photon/electron) yield factor. See above.



    protected:
        G4bool fTrackSecondariesFirst; // see above
        //bools for tracking some special particle cases

        std::unique_ptr<NEST::NESTcalc> fNESTcalc = NULL;
        std::vector<NEST::Vertex> unmerged_steps;



        G4double YieldFactor; // turns scint. on/off

    };

    ////////////////////
    // Inline methods
    ////////////////////

    inline
    G4bool NESTProc::IsApplicable(const G4ParticleDefinition& aParticleType) {
        if (aParticleType.GetParticleName() == "opticalphoton") return false;
        if (aParticleType.IsShortLived()) return false;
        if (aParticleType.GetParticleName() == "thermalelectron") return false;
        //if(abs(aParticleType.GetPDGEncoding())==2112 || //neutron (no E-dep.)
        return true;
    }

    inline
    void NESTProc::SetTrackSecondariesFirst(const G4bool state) {
        fTrackSecondariesFirst = state;
    }

    inline
    G4bool NESTProc::GetTrackSecondariesFirst() const {
        return fTrackSecondariesFirst;
    }

    inline
    void NESTProc::SetScintillationYieldFactor(const G4double yieldfactor) {
        YieldFactor = yieldfactor;
    }

    inline
    G4double NESTProc::GetScintillationYieldFactor() const {
        return YieldFactor;
    }
}





#endif /* NESTPROC_h */
