// $Id: AccSteppingAction.hh based on GEANT4, by Adam Konefal (akonefal@us.edu.pl), may 2006

#ifndef AccSteppingAction_h
#define AccSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class AccSteppingAction : public G4UserSteppingAction
{
  public:
    AccSteppingAction();
    virtual ~AccSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};
#endif
