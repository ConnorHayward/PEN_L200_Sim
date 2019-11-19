#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SiliconPlateConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"

#include "G4Transform3D.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSEnergyDeposit.hh"
#include <G4VPrimitiveScorer.hh>

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VoxelLimits.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4GDMLParser.hh"

#include <G4VisAttributes.hh>
#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
Constructs DetectorConstruction, defines default values.
*/
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),fPBox(nullptr), fLBox(nullptr),
  fBox(nullptr)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fTargetMPT = new G4MaterialPropertiesTable();
  fExpHall_x = fExpHall_y = fExpHall_z = 2.0*m;
  fTargetName = "holder";
  fThickness = 1*mm;
  fTargetThickness = 3*mm;
  fDetectorType = 0;
  fABSL = 1;
  fRES=4.0;
  fLY=10500./MeV;
  fDetectorName = "6pmt_coverage_pe";
  fVolName = "World";
  DefineMaterials();
  SetTargetMaterial("PEN");
  SetWorldMaterial("Air");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
Sets thickness of target.
*/
void DetectorConstruction::SetSize(G4double value){
  fTargetThickness=value;
  if(fBox){
    fBox->SetZHalfLength(fTargetThickness/2);
  }
  UpdateGeometry();

  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetLY(G4double value){
  fLY=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetRes(G4double value){
  fRES=value;
  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets which detector geometry is used.
*/
void DetectorConstruction::SetDetectorType(G4int value){
  fDetectorType=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetABS(G4double value){
  fABSL=value;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetDetectorName(G4String name){
  fDetectorName=name;

  UpdateGeometry();
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of target.
*/
void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fTargetMaterial = pttoMaterial;
    fTargetName = fTargetMaterial->GetName();
    if ( fLBox ) { fLBox->SetMaterial(fTargetMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Sets material of world volume.
*/
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if ( fWLBox ) { fWLBox->SetMaterial(fWorldMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

/*
Defines materials used in simulation. Sets material properties for PEN and other optical components.
*/
void DetectorConstruction::DefineMaterials(){// ------------- Materials -------------
  G4double a, z, density;
  G4int nelements;

  G4NistManager* man = G4NistManager::Instance();

  fLAr = man->FindOrBuildMaterial("G4_lAr");
  fSi = man->FindOrBuildMaterial("G4_Si");
  fCu = man->FindOrBuildMaterial("G4_Cu");
  fGe = man->FindOrBuildMaterial("G4_Ge");
  fNylon = man->FindOrBuildMaterial("G4_NYLON-8062");

  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);

  // Scintillators

  fPEN = new G4Material("PEN", density= 1.3*g/cm3, nelements=3);
  G4int number_of_atoms;
  fPEN->AddElement(O, number_of_atoms=4);
  fPEN->AddElement(H, number_of_atoms=10);
  fPEN->AddElement(C, number_of_atoms=14);

  G4double wavelength;
  char filler;
  G4double varabsorlength;
  G4double ems;
  G4double rindex;

  G4double absEnergy[102]  = {0};
  G4double abs[102]={0};
  G4double emission[102]={0};
  G4double rIndex[102]={0};
  G4double rIndex_fAir[102]={0};
  G4double ems_abs[102]={0};

  G4int absEntries = 0;
  ifstream ReadAbs;

  G4String abs_file = "../input_files/Exp4.csv";
  G4double emission_fibre[102]={0};
  ReadAbs.open(abs_file);
  G4double var = GetABS();
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>varabsorlength>>filler>>ems>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      absEnergy[absEntries] = (1240/wavelength)*eV;
      fABSL = 1;
      abs[absEntries] = varabsorlength*mm;
      emission[absEntries] = ems;
      rIndex[absEntries] = 1.65;
      rIndex_fAir[absEntries]=1.0;
      ems_abs[absEntries]=ems;
      emission_fibre[absEntries]=1.0;
      absEntries++;
    }
  }

  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  const G4int nEntries1 = sizeof(absEnergy)/sizeof(G4double);
  assert(sizeof(rIndex) == sizeof(absEnergy));
  assert(sizeof(abs) == sizeof(absEnergy));
  assert(sizeof(emission) == sizeof(absEnergy));
  assert(sizeof(rIndex_fAir == sizeof(absEnergy)));

  G4double Fiber_energy[] = {2.00*eV,2.87*eV,3.2*eV, 8*eV, 9*eV};
  G4double Fiber_abslength[]={9.00*m, 9.00*m, 0.01*mm,0.01*mm, 0.01*mm};

  // fTargetMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  // fTargetMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true); // *
  // fTargetMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  // fTargetMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  fTargetMPT->AddProperty("WLSCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("WLSABSLENGTH",   absEnergy, emission, nEntries1)->SetSpline(true);

  // fTargetMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  // fTargetMPT->AddConstProperty("RESOLUTIONSCALE",4.0); // * 1, 4, 8
  // fTargetMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  // fTargetMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  // fTargetMPT->AddConstProperty("YIELDRATIO",0.05);
  // fTargetMPT->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);

  fPEN->SetMaterialPropertiesTable(fTargetMPT);

  G4MaterialPropertiesTable* lARMPT = new G4MaterialPropertiesTable();
  abs_file = "../input_files/lArScint.csv";

  ReadAbs.open(abs_file);
  absEntries = 0;
  G4double absEnergylArScint[43] = {};
  G4double emissionlArScint[43] = {};
  G4double abslArScint[43] = {};

  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>ems;//>>filler>>varabsorlength;
      if(ReadAbs.eof()){
        break;
      }
      absEnergylArScint[absEntries] = (1240/wavelength)*eV;
      emissionlArScint[absEntries] = ems;
      G4cout <<1240/wavelength << "  " << ems << "  " << varabsorlength <<G4endl;
      abslArScint[absEntries] = varabsorlength*m;

      absEntries++;
    }
  }
  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  G4double absEnergyLAr[]  = {1.0*eV, 1.9255*eV, 2.145*eV, 2.2704*eV, 2.4378*eV, 2.6085*eV,2.845*eV,3.0515*eV,3.397*eV,5.0*eV};
  G4double absLAr[]={10*m,10*m,10*m,10*m,10*m,10*m,10*m,10*m};
  G4double rIndexLAr[]={1.2295,1.2295, 1.2303, 1.2308,1.2316,1.2324,1.2336,1.2347,1.2367,1.2367};

  lARMPT->AddProperty("RINDEX",       absEnergyLAr, rIndexLAr, 10)->SetSpline(true);
  lARMPT->AddProperty("ABSLENGTH",    absEnergyLAr, absLAr, 10)->SetSpline(true);
  lARMPT->AddProperty("FASTCOMPONENT",absEnergylArScint, emissionlArScint, 43)->SetSpline(true);
  lARMPT->AddProperty("SLOWCOMPONENT",absEnergylArScint, emissionlArScint, 43)->SetSpline(true);
  lARMPT->AddConstProperty("SCINTILLATIONYIELD", 51./MeV);
  lARMPT->AddConstProperty("FASTTIMECONSTANT", 6.2*us);
  lARMPT->AddConstProperty("SLOWTIMECONSTANT",1300*us);
  lARMPT->AddConstProperty("YIELDRATIO",0.05);
  lARMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLAr->SetMaterialPropertiesTable(lARMPT);

  G4double absEnergyNylon[] = {1.9255*eV, 2.145*eV, 2.2704*eV, 2.4378*eV, 2.6085*eV,2.845*eV,3.0515*eV,3.397*eV};
  G4double rIndexNylon[]={1.5395, 1.5303, 1.5308,1.5316,1.5324,1.5336,1.5347,1.5367};
  G4double absNylon[]={1*m, 1*m, 1*m,1*m,1*m,1*m,1*m,1*m};
  G4MaterialPropertiesTable* nylonMPT = new G4MaterialPropertiesTable();
  nylonMPT->AddProperty("RINDEX", absEnergyNylon, rIndexNylon, 8)->SetSpline(true);
  nylonMPT->AddProperty("ABSLENGTH", absEnergyNylon, absNylon, 8)->SetSpline(true);
  fNylon->SetMaterialPropertiesTable(nylonMPT);

}

void DetectorConstruction::SetVolName(G4ThreeVector thePoint){
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* myVolume= theNavigator->LocateGlobalPointAndSetup(thePoint);
  fVolName =  myVolume->GetName();
}

void DetectorConstruction::SetPropertyTable(G4Material* material, G4MaterialPropertiesTable* table){
  material->SetMaterialPropertiesTable(table);
}

void DetectorConstruction::UpdateGeometry(){
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

/*
Clears stored geometry, then constructs all volumes that can be used in the simulation.

Builds and places volumes in world.

Defines detector sensitivities and properties.
*/
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

    G4GDMLParser parser;
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

// ------------- Volumes --------------

// The experimental Hall
  fWorldBox = new G4Tubs("World",0,0.26*m,0.5*m,0.*deg,360.*deg);
  fWLBox = new G4LogicalVolume(fWorldBox,fLAr,"World",0,0,0);
  fWPBox = new G4PVPlacement(0,G4ThreeVector(),fWLBox,"World",0,false,0);

  double rod_sides = 2.5*mm;
  double rod_length = 48*mm;

  G4Tubs* rod = new G4Tubs("target",0,1.5*mm,32.5*cm,0.*deg,360.*deg);
  G4LogicalVolume* rod_log = new G4LogicalVolume(rod,fCu,"target",0,0,0);

  G4double detector_height = 43*mm;
  G4Tubs* geDetector = new G4Tubs("Crystal", 0, 31.5*mm, detector_height, 0.*deg, 360.*deg);
  G4LogicalVolume* ge_log = new G4LogicalVolume(geDetector, fGe, "crystal_log",0,0,0);

  G4Tubs* hit_tracker = new G4Tubs("tracker",0.25*m,0.26*m,0.5*m,0.*deg,360.*deg);
  G4LogicalVolume* tracker_log = new G4LogicalVolume(hit_tracker,fNylon, "tracker",0,0,0);

  // Cables
  G4RotationMatrix* cableRotation = new G4RotationMatrix();
  cableRotation->rotateX(90.*deg);
  cableRotation->rotateY(90.*deg);
  cableRotation->rotateZ(0.*deg);

  G4RotationMatrix* cableRotationTwo = new G4RotationMatrix();
  cableRotationTwo->rotateX(90.*deg);
  cableRotationTwo->rotateY(90.*deg);
  cableRotationTwo->rotateY(180.*deg);

  G4VSolid* toroidSolid =  new G4Torus("ToroidSolid", 0*mm, 0.15*mm, 6*mm, 0.*deg, 140.*deg);
  G4LogicalVolume* toroidLog = new  G4LogicalVolume( toroidSolid, fPEN, "ToroidLog", 0, 0, 0);

  // Copper cables (vertical part)
  G4RotationMatrix* VcoppercableRotation = new G4RotationMatrix();
  VcoppercableRotation->rotateZ(90.*deg);
  VcoppercableRotation->rotateX(0.*deg);
  VcoppercableRotation->rotateY(90.*deg);


  G4VSolid* VcoppercableSolid = new G4Box("VCopperCableSolid", 44.95*mm, 5.0*mm, 0.05*mm);
  G4LogicalVolume* VcoppercableLog = new G4LogicalVolume(VcoppercableSolid, fCu, "VCopperCableLog", 0, 0, 0);

  // Copper cables (horizontal part)
  G4VSolid* H1coppercableSolid = new G4Box("H1CopperCableSolid", 8.5*mm, 20.8*mm, 0.05*mm);
  G4LogicalVolume* H1coppercableLog = new G4LogicalVolume(H1coppercableSolid, fCu, "H1CopperCableLog", 0, 0, 0);
  G4VSolid* H2coppercableSolid = new G4Box("H2CopperCableSolid", 8.5*mm, 8.8*mm, 0.05*mm);
  G4LogicalVolume* H2coppercableLog = new G4LogicalVolume(H2coppercableSolid, fCu, "H2CopperCableLog", 0, 0, 0);

  //G4VPhysicalVolume* toroidPhys = new  G4PVPlacement( cableRotation, G4ThreeVector(0, -1*mm, -0.3*mm), toroidLog, "aToroidPhys", fWLBox, 0, 0);
  // Mini Shroud

  G4Tubs* miniShroud = new G4Tubs("miniShroud", 55*mm, 55.5*mm, 0.4*m, 0.*deg, 360.*deg);
  G4LogicalVolume* shroud_log = new G4LogicalVolume(miniShroud, fNylon, "shroud",0,0,0);

  // Contact

  G4Tubs* contact = new G4Tubs("contact", 0, 2.5*mm, 0.15*mm, 0.*deg, 360.*deg);
  G4LogicalVolume* contact_log = new G4LogicalVolume(contact, fGe, "shroud",0,0,0);

  // Silicon Plates

  fSiliconPlate_h = 0.75*mm;

  SiliconPlateConstruction plate;
  G4VSolid* final_plate = plate.ConstructPlate();
  G4LogicalVolume* plate_log = new G4LogicalVolume(final_plate,fPEN,"plate");
  G4cout <<  G4BestUnit(plate_log->GetMass(true),"Mass") << G4endl;

  // --------------Detectors--------------
  char filler;
  G4double wavelength;

  G4double energyPerfect[7]  = {0.*eV, 1.*eV, 2.*eV, 3.*eV, 3.5*eV, 4.*eV, 10.*eV};
  G4double effPerfect[7]  = {1, 1, 1, 1, 1, 1, 1};
  G4double reflPerfect[7]  = {0, 0, 0, 0,0,0, 0};

  const G4int nPMT_EFF = sizeof(energyPerfect)/sizeof(G4double);

  G4OpticalSurface* perfect_optsurf = new G4OpticalSurface("perfect",glisur,polished, dielectric_metal);
  G4MaterialPropertiesTable* detector_MT = new G4MaterialPropertiesTable();
  detector_MT->AddProperty("EFFICIENCY", energyPerfect, effPerfect,nPMT_EFF)->SetSpline(true);
  detector_MT->AddProperty("REFLECTIVITY", energyPerfect, reflPerfect,nPMT_EFF)->SetSpline(true);
  perfect_optsurf->SetMaterialPropertiesTable(detector_MT);
  new G4LogicalSkinSurface("tracker_surf",tracker_log,perfect_optsurf);

  G4VisAttributes* detectorAttr = new G4VisAttributes(G4Colour::Grey());
  detectorAttr->SetVisibility(true);
  detectorAttr->SetForceSolid(true);
  ge_log->SetVisAttributes(detectorAttr);
  toroidLog->SetVisAttributes(detectorAttr);

  G4VisAttributes* rodAttr = new G4VisAttributes(G4Colour::Brown());
  rodAttr->SetVisibility(true);
  rodAttr->SetForceSolid(true);
  rod_log->SetVisAttributes(rodAttr);
  VcoppercableLog->SetVisAttributes(rodAttr);
  H1coppercableLog->SetVisAttributes(rodAttr);
  H2coppercableLog->SetVisAttributes(rodAttr);

  G4VisAttributes* plateAttr = new G4VisAttributes(G4Colour::Blue());
  plateAttr->SetVisibility(true);
  plateAttr->SetForceSolid(true);
  plate_log->SetVisAttributes(plateAttr);

  G4VisAttributes* contactAttr = new G4VisAttributes(G4Colour::Red());
  contactAttr->SetVisibility(true);
  contactAttr->SetForceSolid(true);
  contact_log->SetVisAttributes(contactAttr);

  G4VPhysicalVolume* crystal_placement;
  G4VPhysicalVolume* rod_placement;
  G4VPhysicalVolume* tracker_placement;
  G4VPhysicalVolume* shroud_placement;
  G4VPhysicalVolume* contact_placement;
  G4VPhysicalVolume* cable_placement;
  G4VPhysicalVolume* Vcoppercable_placement;
  G4VPhysicalVolume* H1coppercable_placement;
  G4VPhysicalVolume* H2coppercable_placement;
  G4double tmpR = 46.*mm + 1*mm + 1.5*mm;

  /*
  0 - PMT on base of tile, collimator included.
  */

  /* Middle String */
  for(int i = 0; i < 4; i++){
    fPBox             = new G4PVPlacement(0, G4ThreeVector(0,0,-fSiliconPlate_h-2.4*mm+i*10*cm),plate_log,"plate",fWLBox,false,i,false);
    crystal_placement = new G4PVPlacement(0, G4ThreeVector(0,0,detector_height+i*10*cm),ge_log,"crystal",fWLBox,false,i,false);
    contact_placement = new G4PVPlacement(0, G4ThreeVector(0,0.5*cm,i*10*cm-0.15*mm),contact_log, "contact",fWLBox, i,false);
    contact_placement = new G4PVPlacement(0, G4ThreeVector(0,1.8*cm,i*10*cm-0.15*mm),contact_log, "contact",fWLBox, i,false);
    fPBox             = new G4PVPlacement(0, G4ThreeVector(0,0,fSiliconPlate_h+2*detector_height+2.4*mm+i*10*cm),plate_log,"top-plate",fWLBox,false,i,false);
    cable_placement   = new  G4PVPlacement( cableRotation, G4ThreeVector(0, -1*mm, -0.3*mm+i*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
    cable_placement   = new  G4PVPlacement( cableRotationTwo, G4ThreeVector(0, 24*mm, -0.3*mm+i*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
    const G4double VertBarAngle = ((G4double) i) * 120.*deg;
    const G4ThreeVector VertBarTranslation (/* x */ tmpR * std::cos(VertBarAngle),  /* y */ tmpR * std::sin(VertBarAngle), /* z */ +11.25*cm);
    rod_placement     = new G4PVPlacement(0, VertBarTranslation, rod_log,"support rod",fWLBox,false,i,false);
    Vcoppercable_placement = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(0, 45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm + i*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
    Vcoppercable_placement = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(0, -45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm + i*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
    H1coppercable_placement = new G4PVPlacement(0,G4ThreeVector(0, -24.2, -2.4 - 2*fSiliconPlate_h  - 0.05*mm + i*10*cm), H1coppercableLog, "H1CopperCablePhys", fWLBox, 0, 0);
    H2coppercable_placement = new G4PVPlacement(0,G4ThreeVector(0, 36.2, -2.4 - 2*fSiliconPlate_h - 0.05*mm + i*10*cm), H2coppercableLog, "H2CopperCablePhys", fWLBox, 0, 0);

  }

  shroud_placement = new G4PVPlacement(0, G4ThreeVector(),shroud_log,"shroud",fWLBox,false,0,false);

  for(int i = 1; i < 3; i++){
    fPBox             = new G4PVPlacement(0, G4ThreeVector(0,0,-fSiliconPlate_h-2.4*mm-i*10*cm),plate_log,"plate",fWLBox,false,3+i,false);
    crystal_placement = new G4PVPlacement(0, G4ThreeVector(0,0,detector_height-i*10*cm),ge_log,"crystal",fWLBox,false,3+i,false);
    contact_placement = new G4PVPlacement(0, G4ThreeVector(0,0.5*cm,-i*10*cm-0.15*mm),contact_log, "contact",fWLBox, 3+i,false);
    contact_placement = new G4PVPlacement(0, G4ThreeVector(0,1.8*cm,-i*10*cm-0.15*mm),contact_log, "contact",fWLBox, 3+i,false);
    fPBox             = new G4PVPlacement(0, G4ThreeVector(0,0,fSiliconPlate_h+2*detector_height+2.4*mm-i*10*cm),plate_log,"top-plate",fWLBox,false,3+i,false);
    cable_placement   = new  G4PVPlacement( cableRotation, G4ThreeVector(0, -1*mm, -0.3*mm-i*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
    cable_placement   = new  G4PVPlacement( cableRotationTwo, G4ThreeVector(0, 24*mm, -0.3*mm-i*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
    Vcoppercable_placement = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(0, 45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm - i*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
    Vcoppercable_placement = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(0, -45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm - i*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
    H1coppercable_placement = new G4PVPlacement(0,G4ThreeVector(0, -24.2, -2.4 - 2*fSiliconPlate_h  - 0.05*mm - i*10*cm), H1coppercableLog, "H1CopperCablePhys", fWLBox, 0, 0);
    H2coppercable_placement = new G4PVPlacement(0,G4ThreeVector(0, 36.2, -2.4 - 2*fSiliconPlate_h - 0.05*mm - i*10*cm), H2coppercableLog, "H2CopperCablePhys", fWLBox, 0, 0);

  }

  G4double radius = 15*cm;
  G4double x;
  G4double y;
  int detcounter = 6;
  int rod_counter = 4;
  int shroud_counter = 1;
  double PI = 3.1415926535897;

  for(int i=0;i<6;i++){
    x = radius * cos(2 * PI * i/6);
    y = radius * sin(2 * PI * i/6);

    for(int j = 0; j < 4; j++){
      fPBox             = new G4PVPlacement(0, G4ThreeVector(x,y,-fSiliconPlate_h+j*10*cm-2.4*mm),plate_log,"plate",fWLBox,false,detcounter+i+j,false);
      crystal_placement = new G4PVPlacement(0, G4ThreeVector(x,y,detector_height+j*10*cm),ge_log,"crystal",fWLBox,false,detcounter+i+j,false);
      contact_placement = new G4PVPlacement(0, G4ThreeVector(x,y+0.5*cm,j*10*cm-0.15*mm),contact_log, "contact",fWLBox, detcounter+i+j,false);
      contact_placement = new G4PVPlacement(0, G4ThreeVector(x,y+1.8*cm,j*10*cm-0.15*mm),contact_log, "contact",fWLBox, 1+detcounter+i+j,false);
      fPBox             = new G4PVPlacement(0, G4ThreeVector(x,y,fSiliconPlate_h+2*detector_height+2.4*mm+j*10*cm),plate_log,"top-plate",fWLBox,false,detcounter+i+j,false);
      const G4double VertBarAngle = ((G4double) j) * 120.*deg;
      const G4ThreeVector VertBarTranslation (/* x */ x+tmpR * std::cos(VertBarAngle),  /* y */ y+tmpR * std::sin(VertBarAngle), /* z */ +11.25*cm);
      rod_placement     = new G4PVPlacement(0, VertBarTranslation, rod_log,"support rod",fWLBox,false,rod_counter+i+j,false);
      cable_placement   = new  G4PVPlacement( cableRotation, G4ThreeVector(x, y-1*mm, -0.3*mm+j*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
      cable_placement   = new  G4PVPlacement( cableRotationTwo, G4ThreeVector(x, y+24*mm, -0.3*mm+j*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
      Vcoppercable_placement  = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(x, y+45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm + j*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
      Vcoppercable_placement  = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(x, y-45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm + j*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
      H1coppercable_placement = new G4PVPlacement(0,G4ThreeVector(x,y -24.2, -2.4 - 2*fSiliconPlate_h  - 0.05*mm + j*10*cm), H1coppercableLog, "H1CopperCablePhys", fWLBox, 0, 0);
      H2coppercable_placement = new G4PVPlacement(0,G4ThreeVector(x,y+ 36.2, -2.4 - 2*fSiliconPlate_h - 0.05*mm + j*10*cm), H2coppercableLog, "H2CopperCablePhys", fWLBox, 0, 0);
    }

    for(int j = 1; j < 3; j++){
      fPBox             = new G4PVPlacement(0, G4ThreeVector(x,y,-fSiliconPlate_h-j*10*cm-2.4*mm),plate_log,"plate",fWLBox,false,3+detcounter+i+j,false);
      crystal_placement = new G4PVPlacement(0, G4ThreeVector(x,y,detector_height-j*10*cm),ge_log,"crystal",fWLBox,false,3+detcounter+i+j,false);
      contact_placement = new G4PVPlacement(0, G4ThreeVector(x,y+0.5*cm,-j*10*cm-0.15*mm),contact_log, "contact",fWLBox, detcounter+i+j,false);
      contact_placement = new G4PVPlacement(0, G4ThreeVector(x,y+1.8*cm,-j*10*cm-0.15*mm),contact_log, "contact",fWLBox, detcounter+i+j,false);
      fPBox             = new G4PVPlacement(0, G4ThreeVector(x,y,fSiliconPlate_h+2*detector_height-j*10*cm+2.4*mm),plate_log,"top-plate",fWLBox,false,3+detcounter+i+j,false);
      cable_placement   = new  G4PVPlacement( cableRotation, G4ThreeVector(x, y-1*mm, -0.3*mm-j*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
      cable_placement   = new  G4PVPlacement( cableRotationTwo, G4ThreeVector(x, y+24*mm, -0.3*mm-j*10*cm), toroidLog, "ToroidPhys", fWLBox, 0, 0);
      Vcoppercable_placement  = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(x, y+45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm - j*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
      Vcoppercable_placement  = new G4PVPlacement(VcoppercableRotation,G4ThreeVector(x, y-45.05*mm, 44.95*mm - 2.4*mm - 2*fSiliconPlate_h - 0.05*mm - j*10*cm), VcoppercableLog, "VCopperCablePhys", fWLBox, 0, 0);
      H1coppercable_placement = new G4PVPlacement(0,G4ThreeVector(x,y -24.2, -2.4 - 2*fSiliconPlate_h  - 0.05*mm - j*10*cm), H1coppercableLog, "H1CopperCablePhys", fWLBox, 0, 0);
      H2coppercable_placement = new G4PVPlacement(0,G4ThreeVector(x,y+ 36.2, -2.4 - 2*fSiliconPlate_h - 0.05*mm - j*10*cm), H2coppercableLog, "H2CopperCablePhys", fWLBox, 0, 0);
    }
    shroud_placement = new G4PVPlacement(0, G4ThreeVector(x,y,0),shroud_log,"shroud",fWLBox,false,shroud_counter+1,false);
    detcounter = detcounter+5;
    rod_counter = rod_counter+3;
  }

  tracker_placement = new G4PVPlacement(0, G4ThreeVector(0,0,0),tracker_log,"tracker",fWLBox,false,0,false);
  return fWPBox;
}
