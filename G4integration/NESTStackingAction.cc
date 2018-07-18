#include "NESTStackingAction.hh"
#include "G4OpticalPhoton.hh"

NESTStackingAction::NESTStackingAction() {}
NESTStackingAction::~NESTStackingAction() {}

NESTStackingAction* NESTStackingAction::theStackingAction = 0;
G4StackManager* NESTStackingAction::savedManager = 0;

G4ClassificationOfNewTrack NESTStackingAction::ClassifyNewTrack(
    const G4Track* track) {
  savedManager = stackManager;
  if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    return fWaiting;
  }
  return fUrgent;
}

bool NESTStackingAction::isUrgentEmpty() {
  return savedManager->GetNUrgentTrack() == 0;
}
void NESTStackingAction::NewStage() {}