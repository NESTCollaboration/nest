#ifndef NESTStackingAction_h
#define NESTStackingAction_h 1

#include "globals.hh"

#include "G4ClassificationOfNewTrack.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4UserStackingAction.hh"

/// Control which particles are tracked by G4

class NESTStackingAction : public G4UserStackingAction {
 public:
  NESTStackingAction();
  virtual ~NESTStackingAction();
  static G4StackManager* savedManager;
  static NESTStackingAction* theStackingAction;

  bool isUrgentEmpty();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
  virtual void NewStage();

 private:
};

#endif
