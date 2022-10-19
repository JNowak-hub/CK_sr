// $Id: AccPrimaryGeneratorAction.hh based on GEANT4, by Adam Konefal (akonefal@us.edu.pl), june 2003

#ifndef AccPrimaryGeneratorAction_h
#define AccPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class AccPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    AccPrimaryGeneratorAction();
    ~AccPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


