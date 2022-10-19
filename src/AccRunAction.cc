// $Id: AccRunAction.cc based on GEANT4, written by Adam Konefal (akonefal@us.edu.pl), may 2006

#include "AccRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VMarker.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4ios.hh"


AccRunAction::AccRunAction()
{
}

AccRunAction::~AccRunAction()
{}

void AccRunAction::BeginOfRunAction(const G4Run* aRun)
{
 
//  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
//G4RunManager::GetRunManager()->SetRandomNumberStore(true);

if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
     
// Markers for visualization

G4Circle circle( G4Point3D (0.,0.,0.));
circle.SetScreenDiameter (0.5);
circle.SetFillStyle (G4Circle::filled);
G4Colour colour(1.,0.,0.);
G4VisAttributes attribs(colour);
circle.SetVisAttributes(attribs);
G4VVisManager::GetConcreteInstance()->Draw(circle);  

G4Circle circle1( G4Point3D (0.,0.,1000.));
circle1.SetScreenDiameter (0.5);
circle1.SetFillStyle (G4Circle::filled);
G4Colour colour1(1.,0.,0.);
G4VisAttributes attribs1(colour1);
circle1.SetVisAttributes(attribs1);
G4VVisManager::GetConcreteInstance()->Draw(circle1);  

    }
}

void AccRunAction::EndOfRunAction(const G4Run* )
 {
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");          

 }
}


