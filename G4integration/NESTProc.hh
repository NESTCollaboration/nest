#ifndef NESTPROC_h
#define NESTPROC_h 1

#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VRestDiscreteProcess.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "templates.hh"
//#include "G4ThermalElectron.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCuts.hh"
#include "NEST.hh"

#define AVO 6.022e23  // Avogadro's number (#/mol)
#define EMASS 9.109e-31 * kg
#define MillerDriftSpeed true

namespace NEST {

struct Hit {
 public:
  Hit(double E, double t, G4ThreeVector xyz) : E(E), t(t), xyz(xyz){};
  double E;
  double t;
  G4ThreeVector xyz;
};

struct Lineage {
 public:
  Lineage(INTERACTION_TYPE type) : type(type){};
  INTERACTION_TYPE type = NoneType;
  std::vector<Hit> hits;
  double density = -1;
  int A = -1;
  int Z = -1;
  NESTresult result;
  bool result_calculated = false;
};

class NESTProc : public G4VRestDiscreteProcess {
 public:  // constructor and destructor
  NESTProc(const G4String& processName = "S1",
           G4ProcessType type = fElectromagnetic, double efield = 0,
           std::unique_ptr<VDetector> detector =
               std::unique_ptr<VDetector>(new VDetector()));
  ~NESTProc();

 public:  // methods, with descriptions
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  // Returns true -> 'is applicable', for any particle type except for an
  // 'opticalphoton' and for short-lived particles

  G4double GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*);
  // Returns infinity; i. e. the process does not limit the step, but
  // sets the 'StronglyForced' condition for the DoIt to be invoked at
  // every step.

  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*);
  // Returns infinity; i. e. the process does not limit the time, but
  // sets the 'StronglyForced' condition for the DoIt to be invoked at
  // every step.

  // For in-flight particles losing energy (or those stopped)
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep);

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

  Lineage GetChildType(const G4Track* aTrack, const G4Track* sec) const;
  double efield = 0;
  G4Track* MakePhoton(G4ThreeVector xyz, double t);

 protected:
  G4bool fTrackSecondariesFirst;  // see above
  // bools for tracking some special particle cases

  std::unique_ptr<NEST::NESTcalc> fNESTcalc = NULL;
  std::vector<NEST::Lineage> lineages;
  std::vector<NEST::Lineage> lineages_prevEvent;
  std::vector<G4Track> photons_prevEvent;
  std::map<std::tuple<int, CLHEP::Hep3Vector, CLHEP::Hep3Vector>,
           long unsigned int>
      track_lins;

  G4double YieldFactor;  // turns scint. on/off
};

////////////////////
// Inline methods
////////////////////

inline G4bool NESTProc::IsApplicable(
    const G4ParticleDefinition& aParticleType) {
  if (aParticleType.GetParticleName() == "opticalphoton") return false;
  if (aParticleType.IsShortLived()) return false;
  if (aParticleType.GetParticleName() == "thermalelectron") return false;
  // if(abs(aParticleType.GetPDGEncoding())==2112 || //neutron (no E-dep.)
  return true;
}

inline void NESTProc::SetScintillationYieldFactor(const G4double yieldfactor) {
  YieldFactor = yieldfactor;
}

inline G4double NESTProc::GetScintillationYieldFactor() const {
  return YieldFactor;
}
}

#endif /* NESTPROC_h */
