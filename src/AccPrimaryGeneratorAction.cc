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

//--------------losowanie energii 
//a: G4double E = (G4UniformRand() * (6700.0 - 40.0)) + 40.0;//20;
//   G4double u = G4UniformRand();
 // G4double g_E = 1.18158 * E * exp(-0.43 * E); 
 //G4double g_E = 3.303938 * E * exp(-1.32574 * E); 

//G4double a =	0.4283;//-3.37822734463395;
//G4double b =	-0.0707;//5.02293044376451;
//G4double c =	-0.736871824788955;

// dla widmo własne
/*
G4double a=-1138.895539026920;
G4double b=158.992446879125;
G4double c=5.637078192402;
G4double d=-0.026370988678;


G4double e =	6.765635028886;
G4double f =	-1.949779946847;
G4double g =	-1.497957354646;
G4double h =	-0.023419170223;
G4double i =	-1.432688801107;
G4double j =	-1.933924650206;
*/
/*
//------od 40 keV do 135 keV
G4double a=0.000000000423803;
G4double b=-0.000000159354536;
G4double c=0.000019719065437;
G4double d=-0.000836388341764;
G4double e=0.012649967292856;
G4double f=-0.043700657562441;
//------od 135 keV do 6.7 MeV
G4double k=6636.6747675079;
G4double h=-1.5201714661;
G4double i=-47.1786367373;
G4double j=-0.7254756429;
*/



/*
G4double e =	2152.136770450120;
G4double f =	-2.014448134633;
G4double g =	-10.560598571779;
G4double h =	-0.098772702788;
*/

/*
// dla widma literaturowego 2
G4double a = -670.30654453751600;
G4double b = 196.29039465063900;
G4double c = -10.98602550519320;
G4double d = 0.74307143664545;

G4double e = 88231.715484015600;
G4double f = -3.716374976830;
G4double g = -12.179986797937;
G4double h = -0.241856341018;
G4double i = 0.046614037000;
G4double j = -1.4211592185994000;
*/
/*
	
G4double d =	3458.884163;//690.481300650671;
G4double e = 	1.34965010101977;
G4double f =	-6.96478269889126;
G4double g =	0.325481188832102;//0.3854811888321020;//0.325481188832102;
G4double h =	-0.00000155813025;
G4double i =	4.95192269501133;
*/
/*
G4double a = 0.2909824701981;
G4double b = 0.0364591575867;
G4double c = -0.0113506073307;
	
	
G4double d = 690.4728016095960;
G4double e = 2.1776209437827;
G4double f = -9.6560734038034;
G4double g = 0.325413449;
*/
/* G4double a=1.2138817058;
 G4double c=-0.7380934141;
 G4double b=-0.2024339491;

 G4double e=-2.6403009412;
 G4double f=4.3103922901;
 G4double g=-0.4422687726;
*/
/*
 G4double Ep = 135.0; //0.547094; //0.469796;
G4double g_E;

if(E <= Ep)
{g_E = (a * pow(E,5)) + (b * pow(E,4)) + (c * pow(E,3)) + (d * pow(E,2)) +  (e * E) + f;} 
if(E > Ep)
{
g_E = k * pow(E,h) * exp(i * pow(E,j));
}// +  i * pow(E,j) ;}

G4double energy =  E;
if(u<g_E) 
{ 
particleGun->SetParticleEnergy(energy*keV);*/
//ofstream spectFileE1;
//spectFileE1.open("spect_E",ios::app);
//spectFileE1 << energy << G4endl;
//}

// else  goto a;


//  particleGun->SetParticleEnergy(1500*keV);

//--------------------- kierunki emisji neutron
  G4double pole = 6; // <-- tu podawac wielkosc pola napromieniania w cm 
//  G4double x1, y1; 
//b:// G4double x1 = (G4UniformRand() - 0.5) * 2 * (pole/2);
  // G4double y1 = (G4UniformRand() - 0.5) * 2 * (pole/2);

   G4double x =  CLHEP::RandGauss::shoot(0.0,2.55);//FWHM = 6mm
   G4double y =  CLHEP::RandGauss::shoot(0.0,2.55);

//  G4double r = (G4UniformRand() - 0.5) * 2  * (pole/2);

//  G4double r = CLHEP::RandGauss::shoot(0.0,2.5)*cm;
//if((r*r)<=((pole/2)*(pole/2)))
//{   
//    G4double fi = G4UniformRand()  * 360;
//G4double fi = 0;
//    x1 = r * sin(fi*M_PI/180);
//    y1 = r * cos(fi*M_PI/180);
//}
//else goto b;

if(((x*x)+(y*y))<=((pole/2)*(pole/2)))
{ 
 particleGun->SetParticlePosition(G4ThreeVector(x*mm,y*mm,-3.0*mm));
// particleGun->SetParticlePosition(G4ThreeVector(0*cm,0*cm,-35*cm));
/*
ofstream spectFileE1;
spectFileE1.open("ROZ",ios::app);
spectFileE1 << x1 / cm << "     " << y1 / cm << G4endl;
*/
}

//else goto b;
 
 particleGun->SetParticleMomentumDirection(G4ThreeVector(0*cm,0*cm,1000*mm));


//----------------------------------------------------------------
 particleGun->GeneratePrimaryVertex(anEvent);
}





