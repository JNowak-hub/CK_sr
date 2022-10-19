// $Id: AccRunAction.hh based on GEANT4, by Adam Konefal (akonefal@us.edu.pl), june 2003

#ifndef AccRunAction_h
#define AccRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class AccRunAction : public G4UserRunAction
{
  public:
    AccRunAction();
   ~AccRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

};

#endif

