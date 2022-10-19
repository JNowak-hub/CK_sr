// $Id: Accelerator.cc based on GEANT4, written by Adam Konefal, june 2003 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "AccRunAction.hh"
#include "AccEventAction.hh"
//#include "AccVisManager.hh"
#include "AccDetectorConstruction.hh"
#include "ML2PhysicsList.hh"
#include "AccPrimaryGeneratorAction.hh"
#include "AccSteppingAction.hh"
#include "AccTrackingAction.hh"
#include "MySession.hh"
//#include "UserHookForAbortState.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif   
 

int main()
{
//for RandGenerator
//G4long myseed = 1;
//CLHEP::HepRandom::setTheSeed(myseed);

//int numberOfEvent = 5;
// Constructing the default run manager
G4RunManager* runManager = new G4RunManager;
//CLHEP::HepRandom::restoreEngineStatus();  

// seting mandatory initialization classes
runManager->SetUserInitialization(new AccDetectorConstruction);
runManager->SetUserInitialization(new ML2PhysicsList);

// seting user's action classes
runManager->SetUserAction(new AccPrimaryGeneratorAction);
runManager->SetUserAction(new AccEventAction);
runManager->SetUserAction(new AccTrackingAction);
runManager->SetUserAction(new AccSteppingAction);
runManager->SetUserAction(new AccRunAction);
 

#ifdef G4VIS_USE
// Initializing visualization manager
G4VisManager* visManager = new G4VisExecutive;// AccVisManager;
visManager->Initialize();
#endif   
 
 
// Initializing G4 kernel
runManager->Initialize();

// geting the pointer to the UI manager and set verbosities
G4UImanager* UI = G4UImanager::GetUIpointer();
// Opening session for writing data to output file
MySession * LoggedSession = new MySession;
UI->SetCoutDestination(LoggedSession);

// commends to set output data
UI->ApplyCommand("/run/verbose 0");
UI->ApplyCommand("/event/verbose 0");
UI->ApplyCommand("/tracking/verbose 0");
   
// starting a run
// runManager->BeamOn(numberOfEvent);
UI->ApplyCommand("/control/execute simulation.mac");
delete LoggedSession;

// opening visualization session
  G4UIsession* session = new G4UIterminal;

  UI->ApplyCommand("/control/execute vis.mac");

 session->SessionStart();
// deleting writing session
 delete session;
 
// job termination
#ifdef G4VIS_USE
delete visManager;
#endif
// deleting session;
delete runManager;
return 0;
}

