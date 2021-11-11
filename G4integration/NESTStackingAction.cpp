#include "NESTStackingAction.hh"
#include "G4OpticalPhoton.hh"
#include <iostream>

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
  if (!savedManager) {
    std::cerr
        << "savedManager not set by NESTStackingAction::ClassifyNewTrack(). "
           "Did you set up NESTStackingAction as your stacking action? Did you "
           "override ClassifyNewTrack and forget to set savedManager?"
        << std::endl;
  }
  return savedManager->GetNUrgentTrack() == 0;
}
void NESTStackingAction::NewStage() {}