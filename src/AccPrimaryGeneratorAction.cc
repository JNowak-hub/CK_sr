// $Id: AccPrimaryGeneratorAction.cc, 2005
// -------------------------------------------------------------------
//    by Adam Konefal, based on GEANT 4
// -------------------------------------------------------------------

#include "AccPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <math.h>
#include <fstream>
#include <iostream>
#include "G4SystemOfUnits.hh"

using namespace std;

AccPrimaryGeneratorAction::AccPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
}

AccPrimaryGeneratorAction::~AccPrimaryGeneratorAction()
{
  delete particleGun;
}

void AccPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

}





