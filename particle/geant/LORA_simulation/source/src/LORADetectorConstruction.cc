#include "LORADetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"

LORADetectorConstruction::LORADetectorConstruction()
:  world_log(0), cover_log(0),
cavity_log(0),scintillator_log_1(0), scintillator_log_2(0),scintillator_log_3(0),scintillator_log_4(0),
world_phys(0), cover_phys(0),
cavity_phys(0), scintillator_phys_1(0),scintillator_phys_2(0),scintillator_phys_3(0),scintillator_phys_4(0)
{;}

LORADetectorConstruction::~LORADetectorConstruction()
{
}

G4VPhysicalVolume* LORADetectorConstruction::Construct()
{
   
	//------------------------------------------------------ materials
	
	G4double a;  // atomic mass
	G4double z;  // atomic number
	G4double density, nel;
   
	
	//Elements
   G4Element* H = new G4Element("Hydrogen", "H", z=1., a=  1.01*g/mole);
   G4Element* C = new G4Element("Carbon",   "C", z=6., a= 12.01*g/mole);
	G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
	G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
	
	G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
	Air->AddElement(N, 70*perCent);
	Air->AddElement(O, 30*perCent);
	
	//Aluminium
	G4Material* Al = new G4Material("Aluminium", z=13., a= 27.0*g/mole, density= 2.70*g/cm3);
   
	//Scintillator polyvinyl-toluene: CH2CH(C6H4CH3)n [= 9 C + 10 H ]

   G4Material* Scinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
   Scinti->AddElement(C, 9);
   Scinti->AddElement(H, 10);
	
	// Print all the materials defined.
	//
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	
   //------------------------------------------------------ volumes
	
	//--------------world
	
	G4double world_length= 350*cm; 
	G4double world_height= 150*cm;
	G4Box* world_box = new G4Box("world_box",0.5*world_length, 0.5*world_length, 0.5*world_height);
	
	world_log = new G4LogicalVolume(world_box, Air, "world_log", 0, 0, 0);
	
	world_phys = new G4PVPlacement(0, G4ThreeVector(), world_log, "world",0, false, 0);
	
	//---------------cover
	
	G4double box_width = 98*cm;
	G4double box_height = 9*cm;
   G4double box_length = 125*cm;
	G4double box_thickness = 1*mm;
   
	G4Box* cover_box = new G4Box("cover_box", 0.5*box_length, 0.5*box_width, 0.5*box_height);
	
	cover_log = new G4LogicalVolume(cover_box, Al, "cover_log", 0, 0, 0);
	
	cover_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0.5*(box_height-world_height)), cover_log, "cover", world_log, false, 0);
             
   G4Box* cavity_box = new G4Box("cavity_box", 0.5*box_length-box_thickness, 0.5*box_width-box_thickness, 0.5*box_height-box_thickness);
   
   cavity_log = new G4LogicalVolume(cavity_box, Air, "chamber_log", 0, 0, 0);
   
   cavity_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0), cavity_log, "cavity", cover_log, false, 0);
                                  
	//---------------scintillators
	
	G4double scintillator_length=47.5*cm;
   G4double scintillator_thickness=3*cm;
	G4double sep=1*cm;
   
	G4Box* scintillator_box = new G4Box("scintillator_box", 0.5*scintillator_length, 0.5*scintillator_length, 0.5*scintillator_thickness);
	
	scintillator_log_1 = new G4LogicalVolume(scintillator_box, Scinti, "scintillator_log_1", 0, 0, 0);
   scintillator_log_2 = new G4LogicalVolume(scintillator_box, Scinti, "scintillator_log_2", 0, 0, 0);
   scintillator_log_3 = new G4LogicalVolume(scintillator_box, Scinti, "scintillator_log_3", 0, 0, 0);
   scintillator_log_4 = new G4LogicalVolume(scintillator_box, Scinti, "scintillator_log_4", 0, 0, 0);
   
   scintillator_phys_1 = new G4PVPlacement(0, G4ThreeVector(-0.5*(scintillator_length+sep),-0.5*(scintillator_length+sep),0), scintillator_log_1, "scintillator_1", cavity_log, false, 0);
	scintillator_phys_2 = new G4PVPlacement(0, G4ThreeVector(-0.5*(scintillator_length+sep),0.5*(scintillator_length+sep),0), scintillator_log_2, "scintillator_2", cavity_log, false, 0);
   scintillator_phys_3 = new G4PVPlacement(0, G4ThreeVector(0.5*(scintillator_length+sep),-0.5*(scintillator_length+sep),0), scintillator_log_3, "scintillator_3", cavity_log, false, 0);
   scintillator_phys_4 = new G4PVPlacement(0, G4ThreeVector(0.5*(scintillator_length+sep),0.5*(scintillator_length+sep),0), scintillator_log_4, "scintillator_4", cavity_log, false, 0);
   
	return world_phys;
}

