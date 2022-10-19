// $Id: AccSteppingAction.cc based on GEANT4, written by Adam Konefal (akonefal@us.edu.pl), may 2006

#include "AccSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include <fstream>
#include <iostream>
#include "G4SystemOfUnits.hh"

AccSteppingAction::AccSteppingAction()
{;}
AccSteppingAction::~AccSteppingAction()
{;}
void AccSteppingAction::UserSteppingAction(const G4Step* aStep)
{
 G4SteppingManager * SM = fpSteppingManager;

G4double EKin;
G4double xval;
G4double yval, zval, xdir, ydir, zdir;

if(aStep->GetTrack()->GetDefinition()->GetParticleName() == "gamma")// "neutron")
  {
  G4VPhysicalVolume * thePrePV = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume * thePostPV = aStep->GetPostStepPoint()->GetPhysicalVolume();
    if(thePostPV != 0)
    {
    G4String thePrePVName = thePrePV->GetName();
    G4String thePostPVName = thePostPV->GetName();
//----------------     
     if(thePrePVName !=  thePostPVName)
      {
        if(thePostPVName ==  "Detector")
         {
        EKin = aStep->GetTrack()->GetKineticEnergy();
	xval = aStep->GetTrack()->GetPosition().x();
        yval = aStep->GetTrack()->GetPosition().y();
        zval = aStep->GetTrack()->GetPosition().z();
        xdir = aStep->GetTrack()->GetMomentumDirection().x();
        ydir = aStep->GetTrack()->GetMomentumDirection().y();
        zdir = aStep->GetTrack()->GetMomentumDirection().z();
	 {
	 std::ofstream specFile;
	 specFile.open("rozklady_fotony",std::ios::app);
	 specFile << EKin / keV <<"     "<< xval / mm <<"     "<< yval / mm <<"     "<< zval / mm <<"     "<< xdir / mm <<"     "<< ydir / mm <<"     "<< zdir / mm << G4endl;
	//G4cout <<  EKin / keV  << G4endl;
	 }
	
    	}
      }
//-----------------      
    } 
  }  



G4double EKin1;
G4double xval1;
G4double yval1, zval1, xdir1, ydir1, zdir1;

if(aStep->GetTrack()->GetDefinition()->GetParticleName() == "e-")
  {
  G4VPhysicalVolume * thePrePV = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume * thePostPV = aStep->GetPostStepPoint()->GetPhysicalVolume();
    if(thePostPV != 0)
    {
    G4String thePrePVName = thePrePV->GetName();
    G4String thePostPVName = thePostPV->GetName();
//----------------     
     if(thePrePVName !=  thePostPVName)
      {
        if(thePostPVName ==  "Detector")
         {
        EKin = aStep->GetTrack()->GetKineticEnergy();
	xval = aStep->GetTrack()->GetPosition().x();
        yval = aStep->GetTrack()->GetPosition().y();
        zval = aStep->GetTrack()->GetPosition().z();
        xdir = aStep->GetTrack()->GetMomentumDirection().x();
        ydir = aStep->GetTrack()->GetMomentumDirection().y();
        zdir = aStep->GetTrack()->GetMomentumDirection().z();
	 {
	 std::ofstream specFile;
	 specFile.open("rozklady_electron",std::ios::app);
	 specFile << EKin1 / keV <<"     "<< xval1 / mm <<"     "<< yval1 / mm <<"     "<< zval1 / mm <<"     "<< xdir1 / mm <<"     "<< ydir1 / mm <<"     "<< zdir1 / mm << G4endl;
	//G4cout <<  EKin / keV  << G4endl;
	 }
	
    	}
      }
//-----------------      
    } 
  }  



}
//    EventNo = evt->GetEventID() + 1;
//    theNameOfProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
//    theStepNumber = aStep->GetTrack()->GetCurrentStepNumber();
//    xval = aStep->GetTrack()->GetPosition().x();
//    xval = aStep->GetPreStepPoint()->GetPosition().x();
//    yval = aStep->GetTrack()->GetPosition().y();
//    yval = aStep->GetPreStepPoint()->GetPosition().y();
//    zval = aStep->GetTrack()->GetPosition().z();
//    zval = aStep->GetPreStepPoint()->GetPosition().z();
//  EKin = aStep->GetTrack()->GetKineticEnergy();
//    EKin = aStep->GetPreStepPoint()->GetKineticEnergy();
//    xdir = aStep->GetTrack()->GetMomentumDirection().x();
//    ydir = aStep->GetTrack()->GetMomentumDirection().y();
//    zdir = aStep->GetTrack()->GetMomentumDirection().z();
//    xdir = aStep->GetPreStepPoint()->GetMomentumDirection().x();
//    ydir = aStep->GetPreStepPoint()->GetMomentumDirection().y();
//    zdir = aStep->GetPreStepPoint()->GetMomentumDirection().z();
