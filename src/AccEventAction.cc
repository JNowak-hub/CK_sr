// $Id: AccEventAction.cc based on GEANT4, written by Adam Konefal (akonefal@us.edu.pl), may 2006

#include "AccEventAction.hh"
//#include "AccDetdoseHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include <fstream>
#include <iostream>

AccEventAction::AccEventAction():printModulo(1000000)       
{;}

AccEventAction::~AccEventAction()
{;}

void AccEventAction::BeginOfEventAction(const G4Event* evt)
{

  
 G4int numberOfTesting=evt->GetEventID()+1;
 if(numberOfTesting%printModulo == 0)
  {
///CLHEP::HepRandom::saveEngineStatus();
std::ofstream infoFile;
infoFile.open("info",std::ios::app);  
infoFile << numberOfTesting / 1000000 << G4endl; 
  }
  
 }


void AccEventAction::EndOfEventAction(const G4Event* evt)
{
  if(evt->GetEventID()  ==  9999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
   if(evt->GetEventID()  ==  19999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
    if(evt->GetEventID()  ==  29999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
    if(evt->GetEventID()  ==  39999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
   if(evt->GetEventID()  ==  49999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
    if(evt->GetEventID()  ==  59999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
    if(evt->GetEventID()  ==  69999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
   if(evt->GetEventID()  ==  79999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
   if(evt->GetEventID()  ==  89999999)
  {
  int n_evt = evt->GetEventID () + 1;
  G4cout << "Number of primary particles: " << n_evt << G4endl;
  }
  
   if(evt->GetEventID() == 99999999)
  {
  int n_evt  = evt->GetEventID () +1;
  G4cout << "Number of primary particles: " << n_evt << " - end of session "<< G4endl;
  }  
  // get number of stored trajectories
  //
//  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
//  G4int n_trajectories = 0;
//  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // extract the trajectories and draw them
  //
/*  if (G4VVisManager::GetConcreteInstance())
    {
     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                            ((*(evt->GetTrajectoryContainer()))[i]);
          trj->DrawTrajectory(50);
        }
    }*/
}