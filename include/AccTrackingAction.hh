// $Id: AccTrackingAction.hh based on GEANT4, by Adam Konefal (akonefal@us.edu.pl), june 2003

#ifndef AccTrackingAction_h
#define AccTrackingAction_h 1

#include "G4UserTrackingAction.hh"


class AccTrackingAction : public G4UserTrackingAction 
{

  public:
    AccTrackingAction(){};
    virtual ~AccTrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);

};

#endif
