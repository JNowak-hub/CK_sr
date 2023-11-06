// $Id: AccTrackingAction.cc based GEANT4, written by Adam Konefal (akonefal@us.edu.pl, may 2006

#include "AccTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

void AccTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
// Creating trajectory only for primaries or neutrons
//  if(aTrack->GetDefinition()->GetParticleName() == "neutron" ||
//     aTrack->GetParentID()==0)
//  { fpTrackingManager->SetStoreTrajectory(true); }
//  else
//  { fpTrackingManager->SetStoreTrajectory(false); }


  if(aTrack->GetParentID()==0)
  { fpTrackingManager->SetStoreTrajectory(true); }
  else
  { fpTrackingManager->SetStoreTrajectory(true); }
}


