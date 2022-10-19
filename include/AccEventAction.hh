// $Id: AccEventAction.hh based on GEANT4,by Adam Konefal (akonefal@us.edu.pl), may 2003

#ifndef AccEventAction_h
#define AccEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class AccEventAction : public G4UserEventAction
{
 public:
    AccEventAction();
    ~AccEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    void SetPrintModulo(G4int val){printModulo = val;};        
  public:
    G4double totE;
    G4int printModulo;
};


#endif

    
