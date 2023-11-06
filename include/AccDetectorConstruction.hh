// $Id: AccDetectorConstruction.hh based on GEANT4, by Adam Konefal (akonefal@us.edu.pl), june 2003

#ifndef AccDetectorConstruction_H
#define AccDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class G4UserLimits;


class AccDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    AccDetectorConstruction();
    ~AccDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();
    G4bool IsUseUserLimits() {return fUseUserLimits;}

  private:
    void DefineMaterials();
    G4VPhysicalVolume* ConstructDetector();
    G4LogicalVolume* logicsphere;
    G4bool fUseUserLimits;
    G4UserLimits* theUserLimitsForSphere;
    G4double theMinEkineCutsInSphere;

};

#endif

