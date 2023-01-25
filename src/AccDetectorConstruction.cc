// $Id: AccDetectorConstruction.cc based on GEANT4, written by Adam Konefal (akonefal@us.edu.pl), may 2006

#include "AccDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

AccDetectorConstruction::AccDetectorConstruction()
{;}

AccDetectorConstruction::~AccDetectorConstruction()
{;}

G4VPhysicalVolume* AccDetectorConstruction::Construct()
{
 
 //------------------------------------------------------- rotations

  G4RotationMatrix* rot90alongY = new G4RotationMatrix();
  rot90alongY->rotateY(90.*deg);
  G4RotationMatrix* rot270alongZandY = new G4RotationMatrix();
  rot270alongZandY->rotateZ(270.*deg).rotateY(270.*deg);
  G4RotationMatrix* rot270alongZand90alongYand180alongX = new G4RotationMatrix();
  rot270alongZand90alongYand180alongX->rotateZ(270.*deg).rotateY(90.*deg).rotateZ(180.*deg);
  G4RotationMatrix* rot90alongYand90alongX = new G4RotationMatrix();
  rot90alongYand90alongX->rotateY(90.*deg).rotateX(90.*deg);
  G4RotationMatrix* rot90alongXand90alongY = new G4RotationMatrix();
  rot90alongXand90alongY->rotateX(90.*deg).rotateY(180.*deg);
  G4RotationMatrix* rot90alongZ = new G4RotationMatrix();
  rot90alongZ->rotateZ(90.*deg); 
  G4RotationMatrix* rot270alongZ = new G4RotationMatrix();
  rot270alongZ->rotateZ(270.*deg); 
  G4RotationMatrix* rot15alongY = new G4RotationMatrix();
  rot15alongY->rotateY(15.*deg); 
  G4RotationMatrix* rot330alongZ = new G4RotationMatrix();
  rot330alongZ->rotateZ(330.*deg); 
  G4RotationMatrix* rot30alongZ = new G4RotationMatrix();
  rot30alongZ->rotateZ(30.*deg); 
  G4RotationMatrix* rot345alongY = new G4RotationMatrix();
  rot15alongY->rotateY(345.*deg); 
  G4RotationMatrix* rot90alongX = new G4RotationMatrix();
  rot90alongX->rotateX(90.*deg); 
  G4RotationMatrix* rot270alongX = new G4RotationMatrix();
  rot270alongX->rotateX(270.*deg); 
  G4RotationMatrix* rot90alongXand180alongZ = new G4RotationMatrix();
  rot90alongXand180alongZ->rotateX(90.*deg).rotateZ(180.*deg); 

 //------------------------------------------------------ materials

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4int ncomponents, natoms, iz, n;
  G4double density, fractionmass;
  G4String name, symbol;
  G4double temperature, pressure, abundance;
  


// ------------------------------------------------------ define Elements

  
   G4Element* elN = new G4Element("Nitrogen","N",7., 14.01*g/mole);
 
  a = 16.00*g/mole;
   G4Element* elO = new G4Element(name="Oxygen",symbol="O", z=8., a);
   G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., 1.00794*g/mole);
   G4Element* elC = new G4Element(name="Carbon", symbol="C", z=1., 12.0107*g/mole);
/* 
  a = 1.01*g/mole;
   G4Element* elH = new G4Element(name="Hydrogen",symbol="H", z=1., a);

  a = 54.94*g/mole;
   G4Element* elMn = new G4Element(name="Manganesse",symbol="Mn", z=25., a);
   
  a = 28.09*g/mole;
   G4Element* elSi = new G4Element(name="Silicon",symbol="Mn", z=14., a);
   
  a = 52.00*g/mole;
   G4Element* elCr = new G4Element(name="Chromium",symbol="Cr", z=24., a);
   
  a = 58.7*g/mole;
   G4Element* elNi = new G4Element(name="Nickel",symbol="Ni", z=28., a);
   

  a = 55.85*g/mole;
   G4Element* elFe = new G4Element(name="Iron",symbol="Fe", z=26., a);

  a = 26.98*g/mole;
   G4Element* elAl = new G4Element(name="Aluminium",symbol="Al", z=13., a);
   
 
 a = 40.078*g/mole;
   G4Element* elCa = new G4Element(name="Calcium",symbol="Al", z=20., a);

  a = 137.372*g/mole;
   G4Element* elBa = new G4Element(name="Barium",symbol="Ba", z=56., a);
   
  a = 32.066*g/mole;
   G4Element* elS = new G4Element(name="Sulphur",symbol="S", z=16., a);

  a=22.99*g/mole;
   G4Element* elNa=new G4Element(name="Sodium",symbol="Na",z=11.,a);
 
  a=24.305*g/mole;
   G4Element* elMg=new G4Element(name="Magnesium",symbol="Mg",z=12.,a);

  a=30.974*g/mole;
   G4Element* elP=new G4Element(name="Phosphorus",symbol="P",z=15.,a);
 
  a=35.453*g/mole;
   G4Element* elCl=new G4Element(name="Chlorine",symbol="Cl",z=17.,a);

  a =39.098*g/mole;
   G4Element* elK=new G4Element(name="Potassium",symbol="K",z=19.,a);

  a=65.38*g/mole;
   G4Element* elZn=new G4Element(name="Zinc",symbol="Zn",z=30.,a);
*/

 // ------------ define materials from elements
 
  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="air", density, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  G4Material* plexiGlass = new G4Material(name="plexiglass", 1.19*g/cm3, ncomponents=3);
  plexiGlass->AddElement(elH, fractionmass=0.080538);
  plexiGlass->AddElement(elC, fractionmass=0.599848);
  plexiGlass->AddElement(elO, fractionmass=0.319614);
/*
//
  density = 0.9*g/cm3;
//  G4Material* paraffin = new G4Material(name="paraffin", density, ncomponents=2);
//  paraffin->AddElement(elC, natoms= 20);
//  paraffin->AddElement(elH, natoms= 42); 

  density = 1.0*g/cm3;
  G4Material* H2O = new G4Material(name="water", density, ncomponents=2);
  H2O->AddElement(elH, natoms= 2);
  H2O->AddElement(elO, natoms= 1); 

// definition of stainless steel

  density = 8.02*g/cm3;
  G4Material* sSteel = new G4Material("stainless_steel",density,5);
  sSteel->AddElement(elMn,0.02);
  sSteel->AddElement(elSi,0.01);
  sSteel->AddElement(elCr,0.19);
  sSteel->AddElement(elNi,0.1);
  sSteel->AddElement(elFe,0.68);

// definition of concrete

  density = 2.3*g/cm3;
  G4Material* concrete = new G4Material("normal_concrete",density,6);
  concrete->AddElement(elSi,0.227915);
  concrete->AddElement(elO,0.60541);
  concrete->AddElement(elH,0.09972);
  concrete->AddElement(elCa,0.04986);
  concrete->AddElement(elAl,0.014245);
  concrete->AddElement(elFe,0.00285);

// definition of BaSO4 concrete
 

  density = 3.5*g/cm3;
  G4Material* Baconcrete = new G4Material("Ba_concrete",density,8);
  Baconcrete->AddElement(elSi,0.1139575);
  Baconcrete->AddElement(elO,0.636105);
  Baconcrete->AddElement(elH,0.04986);
  Baconcrete->AddElement(elCa,0.02493);
  
Baconcrete->AddElement(elAl,0.0071225);
  Baconcrete->AddElement(elFe,0.001425);
  Baconcrete->AddElement(elBa,0.0833);
  Baconcrete->AddElement(elS,0.0833);
*/
// definition of natural W, by relative abundance

G4Isotope* W2 = new G4Isotope(name="W182", iz=74, n=182, a=181.9*g/mole);
G4Isotope* W3 = new G4Isotope(name="W183", iz=74, n=183, a=182.9*g/mole);
G4Isotope* W4 = new G4Isotope(name="W184", iz=74, n=184, a=183.9*g/mole);  
G4Isotope* W6 = new G4Isotope(name="W186", iz=74, n=186, a=185.9*g/mole); 

G4Element* natW = new G4Element(name="natural W", symbol="W", ncomponents=4);
natW->AddIsotope(W2, abundance= 26.3*perCent);
natW->AddIsotope(W3, abundance= 14.3*perCent);
natW->AddIsotope(W4, abundance= 30.8*perCent);
natW->AddIsotope(W6, abundance= 28.6*perCent);
 
// a = 183.85*g/mole;
 density = 19.3*g/cm3; 
 G4Material* W = new G4Material(name="Tungsten",  density, 1);
 W->AddElement(natW, fractionmass=1);
/*
// definition of natural Pb, by relative abundance

G4Isotope* Pb4 = new G4Isotope(name="Pb204", iz=82, n=204, a=203.97*g/mole);
G4Isotope* Pb6 = new G4Isotope(name="Pb206", iz=82, n=206, a=205.97*g/mole);
G4Isotope* Pb7 = new G4Isotope(name="Pb207", iz=82, n=207, a=206.98*g/mole);  
G4Isotope* Pb8 = new G4Isotope(name="Pb208", iz=82, n=208, a=207.98*g/mole); 

G4Element* natPb = new G4Element(name="natural Pb", symbol="Pb", ncomponents=4);
natPb->AddIsotope(Pb4, abundance= 1.5*perCent);
natPb->AddIsotope(Pb6, abundance= 24.1*perCent);
natPb->AddIsotope(Pb7, abundance= 22.1*perCent);
natPb->AddIsotope(Pb8, abundance= 52.3*perCent);
 
// a = 207.2*g/mole;
 density = 11.35*g/cm3; 
 G4Material* Pb = new G4Material(name="Lead",  density, 1);
 Pb->AddElement(natPb, fractionmass=1);


// definition of natural Fe, by relative abundance

G4Isotope* Fe4 = new G4Isotope(name="Fe54", iz=26, n=54, a=53.94*g/mole);
G4Isotope* Fe6 = new G4Isotope(name="Fe56", iz=26, n=56, a=55.935*g/mole);
G4Isotope* Fe7 = new G4Isotope(name="Fe57", iz=26, n=57, a=56.935*g/mole);  
G4Isotope* Fe8 = new G4Isotope(name="Fe58", iz=26, n=58, a=57.933*g/mole); 

G4Element* natFe = new G4Element(name="natural Fe", symbol="Fe", ncomponents=4);
natFe->AddIsotope(Fe4, abundance= 5.8*perCent);
natFe->AddIsotope(Fe6, abundance= 91.2*perCent);
natFe->AddIsotope(Fe7, abundance= 2.2*perCent);
natFe->AddIsotope(Fe8, abundance= .8*perCent);
 
 density = 7.874*g/cm3; 
 G4Material* Fe = new G4Material(name="Iron",  density, 1);
 Fe->AddElement(natFe, fractionmass=1);
*/

// definition natural Cu, by relative abundance

G4Isotope* Cu3 = new G4Isotope(name="Cu63", iz=29, n=63, a=62.93*g/mole);
G4Isotope* Cu5 = new G4Isotope(name="Cu65", iz=29, n=65, a=64.928*g/mole);

G4Element* natCu = new G4Element(name="natural Cu", symbol="Cu", ncomponents=2);
natCu->AddIsotope(Cu3, abundance= 69.17*perCent);
natCu->AddIsotope(Cu5, abundance= 30.83*perCent);

 density = 8.96*g/cm3; 
 G4Material* Cu = new G4Material(name="Copper", density, 1);
 Cu->AddElement(natCu, fractionmass=1);
/*
 a = 180.95*g/mole;
 density = 16.6*g/cm3; 
 G4Material* elTa = new G4Material(name="Tantalum", z=73., a, density);

 
  a = 12.01*g/mole;
  density = 2.267*g/cm3;
  G4Material* C = new G4Material(name="Carbon",z=6., a, density);

 a = 115*g/mole;
 density = 7.3*g/cm3;
 G4Material* In = new G4Material(name="Ind", z=49., a, density);
*/

  // definition of vacuum

  density = 1.2-20*g/cm3; //universe_mean_density; // from PhysicalConstants.h
  pressure = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* vacuum = new G4Material(name = "Galactic", z=1., a=1.01*g/mole,
  density, kStateGas, temperature, pressure);

  
  //------------------------------------------------------ volumes

  //--------------------------------------------- World  (world volume)
 
  G4double World_x = 1.2*m;
  G4double World_y = 1.2*m;
  G4double World_z = 1.2*m;
  
  G4Box* solidWorld = new G4Box("World", // name of the volume
                                World_x,
                                World_y,
                                World_z); // size of the volume
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, // name of the solid
                                                    Air, // material of the volume
				                                            "World"); // name of the volume
						    
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0, // no rotation
                                                  G4ThreeVector(), // position at (0,0,0)
                                                 "World",  // name of the volume
                                                  logicWorld, // name of the logical volume
					                                        0, // its mother volume
					                                        false, // no boolean operation
					                                        0); // copy number

 
 // ------------------------------------------------------phantom
 G4double phantomX = 67.5/2*cm;
 G4double phantomY = 64.5/2*cm;
 G4double phantomZ = 56./2*cm;
 
 G4Box* phantom = new G4Box("Phantom", phantomX, phantomY, phantomY);
 G4LogicalVolume* phantomVolume = new G4LogicalVolume(phantom, plexiGlass, "Phantom");
 G4VPhysicalVolume* phantomPhysVolume = new G4PVPlacement(0, G4ThreeVector(),phantomVolume, "Phantom",logicWorld,false,0);

 G4VisAttributes* phantomVisAtt = new G4VisAttributes(G4Colour(0,0,1));
 phantomVisAtt->SetForceSolid(false);
 phantomVolume->SetVisAttributes(phantomVisAtt);

 G4double phantomWaterX = 66./2*cm;
 G4double phantomWaterY = 63./2*cm;
 G4double phantomWaterZ = 41./2*cm;

 G4Box* phantomWater = new G4Box("Phantom", phantomX, phantomY, phantomY);
 G4LogicalVolume* phantomWaterVolume = new G4LogicalVolume(phantom, plexiGlass, "Phantom");
 G4VPhysicalVolume* phantomWaterPhysVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., -7.5*cm) ,phantomWaterVolume, "Phantom", phantomVolume, false, 0);

 G4VisAttributes* phantomWaterVisAtt = new G4VisAttributes(G4Colour(0,0,2));
 phantomWaterVisAtt->SetForceSolid(true);
 phantomWaterVolume->SetVisAttributes(phantomVisAtt);
 return physWorld;
 }




