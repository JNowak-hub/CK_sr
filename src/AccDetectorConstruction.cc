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
/*
//  density = 0.9*g/cm3;
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
                                World_x,World_y,World_z); // size of the volume
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

 
 // ------------------------------------------------------target
 
   G4double target_pRMin = 0.0*cm;   
   G4double target_pRMax = (0.76/2)*cm;   
   G4double target_pDz   = (0.4/2)*cm;   
   G4double target_pSPhi = 0.*deg;   
   G4double target_pDPhi = 360.*deg; 
   G4Tubs* solidtarget = new G4Tubs("target",target_pRMin,target_pRMax,target_pDz,     
                                     target_pSPhi,target_pDPhi);
   G4LogicalVolume* logictarget = new G4LogicalVolume(solidtarget,Cu,"target");
   G4double targetPos_x = 0.0*cm;
   G4double targetPos_y = 0.0*cm;
   G4double targetPos_z = 0.0*cm;
    
   G4VPhysicalVolume* phystarget
                  = new G4PVPlacement(0,
                                G4ThreeVector(targetPos_x,targetPos_y,targetPos_z),
				logictarget,"target",logicWorld,false,0);  

      
   G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(1,1,0));
   targetVisAtt->SetForceSolid(true);
   logictarget->SetVisAttributes(targetVisAtt);
   

   
   // ------------------------------------------------------collimator W1
 
   G4double colW1_pRMin = (1.86/2)*cm;   
   G4double colW1_pRMax = (18.14/2)*cm;   
   G4double colW1_pDz   = (5.28/2)*cm;   
   G4double colW1_pSPhi = 0.*deg;   
   G4double colW1_pDPhi = 360.*deg; 
   G4Tubs* solidcolW1 = new G4Tubs("colW1",colW1_pRMin,colW1_pRMax,colW1_pDz,     
                                     colW1_pSPhi,colW1_pDPhi);
   G4LogicalVolume* logiccolW1 = new G4LogicalVolume(solidcolW1,W,"colW1");
   G4double colW1Pos_x = 0.0*cm;
   G4double colW1Pos_y = 0.0*cm;
   G4double colW1Pos_z = (-5.28/2+0.2)*cm;
    
   G4VPhysicalVolume* physcolW1
                  = new G4PVPlacement(0,
                                G4ThreeVector(colW1Pos_x,colW1Pos_y,colW1Pos_z),
				logiccolW1,"colW1",logicWorld,false,0);  

      
   G4VisAttributes* colW1VisAtt = new G4VisAttributes(G4Colour(0,0,1));
   colW1VisAtt->SetForceSolid(true);
   logiccolW1->SetVisAttributes(colW1VisAtt);

 
   // ------------------------------------------------------collimator W2a
 
   G4double colW2a_pRMin = (0.57/2)*cm;   
   G4double colW2a_pRMax = (18.14/2)*cm;   
   G4double colW2a_pDz   = (1.43/2)*cm;   
   G4double colW2a_pSPhi = 0.*deg;   
   G4double colW2a_pDPhi = 360.*deg; 
   G4Tubs* solidcolW2a = new G4Tubs("colW2a",colW2a_pRMin,colW2a_pRMax,colW2a_pDz,     
                                     colW2a_pSPhi,colW2a_pDPhi);
   G4LogicalVolume* logiccolW2a = new G4LogicalVolume(solidcolW2a,W,"colW2a");
   G4double colW2aPos_x = 0.0*cm;
   G4double colW2aPos_y = 0.0*cm;
   G4double colW2aPos_z = (1.43/2+0.2)*cm;
    
   G4VPhysicalVolume* physcolW2a
                  = new G4PVPlacement(0,
                                G4ThreeVector(colW2aPos_x,colW2aPos_y,colW2aPos_z),
				logiccolW2a,"colW2a",logicWorld,false,0);  

      
   G4VisAttributes* colW2aVisAtt = new G4VisAttributes(G4Colour(1,0,0));
   colW2aVisAtt->SetForceSolid(true);
   logiccolW2a->SetVisAttributes(colW2aVisAtt);


   // ------------------------------------------------------collimator W2b
 
  G4double phiStart = 0;
  G4double phiTotal = 360*deg;
  G4int numZPlanes = 2;
  const G4double zPlane2b[] = { (1.43+0.2)*cm,(1.43+0.2+6.7)*cm};
  const G4double rInner2b[] = { (0.57/2)*cm, (6.43/2)*cm};
  const G4double rOuter2b[] = { (18.14/2)*cm, (18.14/2)*cm};

  G4Polycone* solidpcolW2b = new G4Polycone("W2b",phiStart,phiTotal,numZPlanes,
                                                           zPlane2b,rInner2b,rOuter2b);
							   
   G4LogicalVolume* logiccolW2b = new G4LogicalVolume(solidpcolW2b,W,"W2b");				     
   
   G4double pcolW2bPos_x = 0.0*cm;
   G4double pcolW2bPos_y = 0.0*cm;
   G4double pcolW2bPos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW2b
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW2bPos_x,pcolW2bPos_y,pcolW2bPos_z),
				logiccolW2b,"W2b",logicWorld,false,0);

      
   G4VisAttributes* colW2bVisAtt = new G4VisAttributes(G4Colour(0,1,0));
   colW2bVisAtt->SetForceSolid(true);
   logiccolW2b->SetVisAttributes(colW2bVisAtt);
   

   // ------------------------------------------------------collimator W3
 
 // G4double phiStart = 0;
//  G4double phiTotal = 360*deg;
//  G4int numZPlanes = 2;
  const G4double zPlaneW3[] = { (1.43+0.2+6.7)*cm,(1.43+0.2+6.7+3.57)*cm};
  const G4double rInnerW3[] = { (1.29/2)*cm, (1.86/2)*cm};
  const G4double rOuterW3[] = { (15.57/2)*cm, (11.43/2)*cm};

  G4Polycone* solidpcolW3 = new G4Polycone("W3",phiStart,phiTotal,numZPlanes,
                                                           zPlaneW3,rInnerW3,rOuterW3);
							   
   G4LogicalVolume* logiccolW3 = new G4LogicalVolume(solidpcolW3,W,"W3");				     
   
   G4double pcolW3Pos_x = 0.0*cm;
   G4double pcolW3Pos_y = 0.0*cm;
   G4double pcolW3Pos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW3
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW3Pos_x,pcolW3Pos_y,pcolW3Pos_z),
				logiccolW3,"W3",logicWorld,false,0);

      
   G4VisAttributes* colW3VisAtt = new G4VisAttributes(G4Colour(1,0,1));
   colW3VisAtt->SetForceSolid(true);
   logiccolW3->SetVisAttributes(colW3VisAtt);

 // ------------------------------------------------------collimator W4a
 
 // G4double phiStart = 0;
 // G4double phiTotal = 360*deg;
 // G4int numZPlanes = 2;
  const G4double zPlaneW4a[] = { (1.43+0.2+6.7+3.57+6.0)*cm, (1.43+0.2+6.7+3.57+6.0+1.0)*cm};
  const G4double rInnerW4a[] = { (1.29/2)*cm, (1.86/2)*cm};
  const G4double rOuterW4a[] = { (4.57/2)*cm, (4.57/2)*cm};

  G4Polycone* solidpcolW4a = new G4Polycone("W4a",phiStart,phiTotal,numZPlanes,
                                                           zPlaneW4a,rInnerW4a,rOuterW4a);
							   
   G4LogicalVolume* logiccolW4a = new G4LogicalVolume(solidpcolW4a,W,"W4a");				     
   
   G4double pcolW4aPos_x = 0.0*cm;
   G4double pcolW4aPos_y = 0.0*cm;
   G4double pcolW4aPos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW4a
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW4aPos_x,pcolW4aPos_y,pcolW4aPos_z),
				logiccolW4a,"W4a",logicWorld,false,0);

      
   G4VisAttributes* colW4aVisAtt = new G4VisAttributes(G4Colour(1,1,0));
   colW4aVisAtt->SetForceSolid(true);
   logiccolW4a->SetVisAttributes(colW4aVisAtt);

// ------------------------------------------------------collimator W4b
 
 // G4double phiStart = 0;
 // G4double phiTotal = 360*deg;
 // G4int numZPlanes = 2;
  const G4double zPlaneW4b[] = { (1.43+0.2+6.7+3.57+6.0+1.0)*cm, (1.43+0.2+6.7+3.57+6.0+1.0+2.29)*cm};
  const G4double rInnerW4b[] = { (1.71/2)*cm, (1.86/2)*cm};
  const G4double rOuterW4b[] = { (14.7/2)*cm, (14.7/2)*cm};

  G4Polycone* solidpcolW4b = new G4Polycone("W4b",phiStart,phiTotal,numZPlanes,
                                                           zPlaneW4b,rInnerW4b,rOuterW4b);
							   
   G4LogicalVolume* logiccolW4b = new G4LogicalVolume(solidpcolW4b,W,"W4b");				     
   
   G4double pcolW4bPos_x = 0.0*cm;
   G4double pcolW4bPos_y = 0.0*cm;
   G4double pcolW4bPos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW4b
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW4bPos_x,pcolW4bPos_y,pcolW4bPos_z),
				logiccolW4b,"W4b",logicWorld,false,0);

      
   G4VisAttributes* colW4bVisAtt = new G4VisAttributes(G4Colour(1,1,0));
   colW4bVisAtt->SetForceSolid(true);
   logiccolW4b->SetVisAttributes(colW4bVisAtt);

// ------------------------------------------------------collimator W4c
 
 // G4double phiStart = 0;
 // G4double phiTotal = 360*deg;
 // G4int numZPlanes = 2;
  const G4double zPlaneW4c[] = { (1.43+0.2+6.7+3.57+6.0+1.0+2.29)*cm, (1.43+0.2+6.7+3.57+6.0+1.0+2.29+0.5)*cm};
  const G4double rInnerW4c[] = { (1.86/2)*cm, (1.93/2)*cm};
  const G4double rOuterW4c[] = { (8.43/2)*cm, (8.42/2)*cm};

  G4Polycone* solidpcolW4c = new G4Polycone("W4c",phiStart,phiTotal,numZPlanes,
                                                           zPlaneW4c,rInnerW4c,rOuterW4c);
							   
   G4LogicalVolume* logiccolW4c = new G4LogicalVolume(solidpcolW4c,W,"W4c");				     
   
   G4double pcolW4cPos_x = 0.0*cm;
   G4double pcolW4cPos_y = 0.0*cm;
   G4double pcolW4cPos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW4c
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW4cPos_x,pcolW4cPos_y,pcolW4cPos_z),
				logiccolW4c,"W4c",logicWorld,false,0);

      
   G4VisAttributes* colW4cVisAtt = new G4VisAttributes(G4Colour(1,1,0));
   colW4cVisAtt->SetForceSolid(true);
   logiccolW4c->SetVisAttributes(colW4cVisAtt);

// ------------------------------------------------------collimator W4d
 
 // G4double phiStart = 0;
 // G4double phiTotal = 360*deg;
 // G4int numZPlanes = 2;
  const G4double zPlaneW4d[] = { (1.43+0.2+6.7+3.57+6.0+1.0+2.29+0.5)*cm, (1.43+0.2+6.7+3.57+6.0+1.0+2.29+0.5+4.71)*cm};
  const G4double rInnerW4d[] = { (1.93/2)*cm, (2.29/2)*cm};
  const G4double rOuterW4d[] = { (6.43/2)*cm, (6.43/2)*cm};

  G4Polycone* solidpcolW4d = new G4Polycone("W4d",phiStart,phiTotal,numZPlanes,
                                                           zPlaneW4d,rInnerW4d,rOuterW4d);
							   
   G4LogicalVolume* logiccolW4d = new G4LogicalVolume(solidpcolW4d,W,"W4d");				     
   
   G4double pcolW4dPos_x = 0.0*cm;
   G4double pcolW4dPos_y = 0.0*cm;
   G4double pcolW4dPos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW4d
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW4dPos_x,pcolW4dPos_y,pcolW4dPos_z),
				logiccolW4d,"W4d",logicWorld,false,0);

      
   G4VisAttributes* colW4dVisAtt = new G4VisAttributes(G4Colour(1,1,0));
   colW4dVisAtt->SetForceSolid(true);
   logiccolW4d->SetVisAttributes(colW4dVisAtt);

 // ------------------------------------------------------collimator W5
 
 // G4double phiStart = 0;
 // G4double phiTotal = 360*deg;
//  G4int numZPlanes = 2;
  const G4double zPlane[] = { (1.43+0.2+6.7+3.57+6.0+1.0+2.29+0.5+4.71+6.0)*cm,(1.43+0.2+6.7+3.57+6.0+1.0+2.29+0.5+4.71+6.0+9.22)*cm};
  const G4double rInner[] = { (2.45/2)*cm, (3.02/2)*cm};
  const G4double rOuter[] = { (4.0/2)*cm, (5.7/2)*cm};

  G4Polycone* solidpcolW5 = new G4Polycone("W5",phiStart,phiTotal,numZPlanes,
                                                           zPlane,rInner,rOuter);
							   
   G4LogicalVolume* logiccolW5 = new G4LogicalVolume(solidpcolW5,W,"W5");				     
   
   G4double pcolW5Pos_x = 0.0*cm;
   G4double pcolW5Pos_y = 0.0*cm;
   G4double pcolW5Pos_z = 0.0*cm;
   G4VPhysicalVolume* physpcolW5
                  = new G4PVPlacement(0,
                                G4ThreeVector(pcolW5Pos_x,pcolW5Pos_y,pcolW5Pos_z),
				logiccolW5,"W5",logicWorld,false,0);

      
   G4VisAttributes* colW5VisAtt = new G4VisAttributes(G4Colour(0,1,1));
   colW5VisAtt->SetForceSolid(true);
   logiccolW5->SetVisAttributes(colW5VisAtt);

  //   ----------------------------------------------detectors

 // Detector for spectrum measurements
 
   G4double Det_pRMin = 0.0*cm;   
   G4double Det_pRMax = 3.0*cm;   
   G4double Det_pDz   = .05*mm;   
   G4double Det_pSPhi = 0.*deg;    
   G4double Det_pDPhi = 360.*deg;   
   G4Tubs* solidDet = new G4Tubs("Detector",Det_pRMin,Det_pRMax,Det_pDz,
                                       Det_pSPhi,Det_pDPhi);
  
  G4LogicalVolume* 
                     logicDet = new G4LogicalVolume(solidDet,Air,"Detector");
  G4double DetPos_x = 0.0*cm;
  G4double DetPos_y = 0.0*cm;
  G4double DetPos_z = 59.9*cm;
  G4VPhysicalVolume* physDet
    = new G4PVPlacement(0,
                        G4ThreeVector(DetPos_x,DetPos_y,DetPos_z),
                        logicDet,"Detector",logicWorld,false,0); 


   G4VisAttributes* DetVisAtt = new G4VisAttributes(G4Colour(1,1,1));
   DetVisAtt->SetForceSolid(true);
   logicDet->SetVisAttributes(DetVisAtt);

 return physWorld;
 }




