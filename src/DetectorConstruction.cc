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

  // fAir
  //
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

  fAir = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  fAir->AddElement(N, 70.*perCent);
  fAir->AddElement(O, 30.*perCent);

  G4NistManager* man = G4NistManager::Instance();

  fGlass = man->FindOrBuildMaterial("G4_Pyrex_Glass");
  fTeflon = man->FindOr

		  			break;BuildMaterial("G4_TEFLON");
  fLAr = man->FindOrBuildMaterial("G4_lAr");
  fAl = man->FindOrBuildMaterial("G4_Al");
  fSi = man->FindOrBuildMaterial("G4_Si");
  fPb = man->FindOrBuildMaterial("G4_Pb");
  fCu = man->FindOrBuildMaterial("G4_Cu");
  fGe = man->FindOrBuildMaterial("G4_Ge");

  G4MaterialPropertiesTable* glassMPT = new G4MaterialPropertiesTable();
  G4double photonEnergy2[] = {2.479684*eV, 2.610194*eV, 2.755204*eV, 2.883353*eV, 2.917275*eV, 3.099605*eV};
  G4double refractiveIndex4[] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double absorption[] = {10*m,10*m,10*m,10*m,10*m,10*m};
  assert(sizeof(refractiveIndex4) == sizeof(photonEnergy2));
  glassMPT->AddProperty("RINDEX",photonEnergy2,refractiveIndex4,6)->SetSpline(true);
  glassMPT->AddProperty("ABSLENGTH",photonEnergy2,absorption,6);
  fGlass->SetMaterialPropertiesTable(glassMPT);

  // Water
  //
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);

  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
  water->AddElement(H, 2);
  water->AddElement(O, 1);

  G4Element* C = new G4Element("Carbon", "C", z=12, a=12*g/mole);
  G4Element* Pb = new G4Element("Lead", "Pb", z=87, a=207*g/mole);

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

  /* Scintillation properties */
  G4String abs_file = "../input_files/Exp4.csv";
  G4double emission_fibre[102]={0};
  ReadAbs.open(abs_file);
  G4double var = GetABS();
  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>ems>>filler>>varabsorlength>>filler>>rindex;
      if(ReadAbs.eof()){
        break;
      }
      absEnergy[absEntries] = (1240/wavelength)*eV;
      abs[absEntries] = 10*varabsorlength*cm;
      emission[absEntries] = ems;
      rIndex[absEntries] = 1.65; // 1.4906 for PMMA
      rIndex_fAir[absEntries]=1.0;
      ems_abs[absEntries]=0.02;
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

  fTargetMPT->AddProperty("RINDEX",       absEnergy, rIndex, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("ABSLENGTH",    absEnergy, abs, nEntries1)->SetSpline(true); // *
  fTargetMPT->AddProperty("FASTCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);
  fTargetMPT->AddProperty("SLOWCOMPONENT",absEnergy, emission, nEntries1)->SetSpline(true);

  fTargetMPT->AddConstProperty("SCINTILLATIONYIELD",10500./MeV); // * 2.5 * PEN = PS, 10*PEN=PS
  fTargetMPT->AddConstProperty("RESOLUTIONSCALE",4.0); // * 1, 4, 8
  fTargetMPT->AddConstProperty("FASTTIMECONSTANT", 5.198*ns);
  fTargetMPT->AddConstProperty("SLOWTIMECONSTANT",24.336*ns);
  fTargetMPT->AddConstProperty("YIELDRATIO",0.05);

  G4double wls_wavelength[3] = {2.5*eV,3.1*eV,9.6863*eV};
  G4double wls_length[3] = {10*m,5*cm, 0.01*mm};
  G4double wls_output[3] = {1,0,0};

    // G4cout << "---------------------------------------------------------------------------------" << G4endl;
    // G4cout << "---------------------------------------------------------------------------------" << G4endl;
    // G4cout << "---------------------------------------------------------------------------------" << G4endl;
    // G4cout << "---------------------------------------------------------------------------------" << G4endl;
    // G4cout << "---------------------------------------------------------------------------------" << G4endl;
    fTargetMPT->AddProperty("WLSABSLENGTH",wls_wavelength, wls_length, 3);//->SetSpline(true);
    fTargetMPT->AddProperty("WLSCOMPONENT",wls_wavelength, wls_output, 3);//->SetSpline(true);
    fTargetMPT->AddConstProperty("WLSTIMECONSTANT", 50*ns);


  fPEN->SetMaterialPropertiesTable(fTargetMPT);

  density = universe_mean_density;    //from PhysicalConstants.h
  fVacuum = new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                           kStateGas,2.73*kelvin,3.e-18*pascal);
  //
  // fAir
  G4MaterialPropertiesTable* worldMPT = new G4MaterialPropertiesTable();
  worldMPT->AddProperty("RINDEX", absEnergy, rIndex_fAir, nEntries1)->SetSpline(true);

  fAir->SetMaterialPropertiesTable(worldMPT);
  fVacuum->SetMaterialPropertiesTable(worldMPT);

  G4MaterialPropertiesTable* lARMPT = new G4MaterialPropertiesTable();
  abs_file = "../input_files/lArScint.csv";

  ReadAbs.open(abs_file);
  absEntries = 0;
  G4double absEnergylArScint[43] = {};
  G4double emissionlArScint[43] = {};

  if(ReadAbs.is_open())
  {
    while(!ReadAbs.eof())
    {
      ReadAbs>>wavelength>>filler>>ems;
      if(ReadAbs.eof()){
        break;
      }
      absEnergylArScint[absEntries] = (1240/wavelength)*eV;
      emissionlArScint[absEntries] = ems;
      absEntries++;
    }
  }
  else G4cout<<"Error opening file: " <<abs_file<<G4endl;
  ReadAbs.close();
  absEntries--;

  G4double absEnergyLAr[]  = {1.9255*eV, 2.145*eV, 2.2704*eV, 2.4378*eV, 2.6085*eV,2.845*eV,3.0515*eV,3.397*eV};
  G4double absLAr[]={10*m,10*m,10*m,10*m,10*m,10*m,10*m,10*m};
  G4double rIndexLAr[]={1.2295, 1.2303, 1.2308,1.2316,1.2324,1.2336,1.2347,1.2367};

  lARMPT->AddProperty("RINDEX",       absEnergyLAr, rIndexLAr, 8)->SetSpline(true);
  lARMPT->AddProperty("ABSLENGTH",    absEnergyLAr, absLAr, 8)->SetSpline(true);
  lARMPT->AddProperty("FASTCOMPONENT",absEnergylArScint, emissionlArScint, 43)->SetSpline(true);
  lARMPT->AddProperty("SLOWCOMPONENT",absEnergylArScint, emissionlArScint, nEntries1)->SetSpline(true);
  lARMPT->AddConstProperty("SCINTILLATIONYIELD", 51.*MeV);
  lARMPT->AddConstProperty("FASTTIMECONSTANT", 6.2*us);
  lARMPT->AddConstProperty("SLOWTIMECONSTANT",1300*us);
  lARMPT->AddConstProperty("YIELDRATIO",0.05);
  lARMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLAr->SetMaterialPropertiesTable(lARMPT)

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

  double thickness_foil = 2.5*mm;
  double holeHeight = fThickness;
  G4Box* absorberPb_box = new G4Box("absorber",15*mm,15*mm, 1.5*mm);
  G4Tubs* hole = new G4Tubs("hole",0*mm,1.5*mm,holeHeight, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator = new G4SubtractionSolid("collimator",absorberPb_box,hole);

  G4Box* absorber_foil = new G4Box("foil",15*mm,15*mm,thickness_foil);
  G4Tubs* hole_foil = new G4Tubs("hole_foil",0*mm,1.5*mm,thickness_foil, 0.*deg, 360.*deg);
  G4SubtractionSolid* collimator_foil = new G4SubtractionSolid("collimator_foil",absorber_foil,hole_foil);

  double target_width = 15*mm;
  fBox = new G4Box("target", target_width, target_width, fTargetThickness);
  fLBox = new G4LogicalVolume(fBox,fScintilator, "target",0,0,0);
  double position = 0;

  double rod_sides = 2.5*mm;
  double rod_length = 48*mm;

  G4Tubs* rod = new G4Tubs("target",0,1.5*mm,18*cm,0.*deg,360.*deg);
  G4LogicalVolume* rod_log = new G4LogicalVolume(rod,fCu,"target",0,0,0);

  G4double detector_height = 14*mm;
  G4Tubs* geDetector = new G4Tubs("Crystal", 0, 31.5*mm, detector_height, 0.*deg, 360.*deg);
  G4LogicalVolume* ge_log = new G4LogicalVolume(geDetector, fGe, "crystal_log",0,0,0);

  G4Tubs* hit_tracker = new G4Tubs("tracker",0.25*m,0.26*m,0.5*m,0.*deg,360.*deg);
  G4LogicalVolume* tracker_log = new G4LogicalVolume(hit_tracker,fLAr, "tracker",0,0,0);

  G4NistManager* man = G4NistManager::Instance();

  // Silicon Plates

  fSiliconPlate_h = 1.5*mm;
  fHolderWidth = 90.00*mm;

  G4double trg_b = (fHolderWidth/sqrt(2.))*mm;
  G4double trg_h = 45.*mm;
  G4double rect_x = 2.*(trg_h + 9.5*mm);
  G4double rect_y = fHolderWidth;
  G4double rect_bc_x = 28.*mm;
  G4double cut_x = 20.*mm;
  G4double cut_y = 31.*mm;

  SiliconPlateConstruction plate;
  G4VSolid* final_plate = plate.ConstructPlate();
  G4LogicalVolume* plate_log = new G4LogicalVolume(final_plate,fTargetMaterial,"plate");
  G4cout <<  G4BestUnit(plate_log->GetMass(true),"Mass") << G4endl;
  G4ThreeVector point = G4ThreeVector(0,0,5*cm);

  // Dummy plate for side detector

  G4VSolid* plate_side_det = new G4Box("side_det",rect_x/2,trg_h+1*mm,0.45*fSiliconPlate_h);
  G4VSolid* side_det = new G4SubtractionSolid("side_det",plate_side_det,final_plate);

  G4Navigator* pointNavigator = new G4Navigator();
  pointNavigator->SetWorldVolume(fWPBox);
  pointNavigator->LocateGlobalPointAndSetup(point);

  // --------------Detectors--------------
  char filler;
  G4double wavelength;

  G4double cath_eff;
  G4double photocath_energy[57];
  G4double photocath_EFF[57];
  G4double perfect_EFF[57];
  G4double perfect_REFL[57];
  G4double photocath_REFL[57]={0};
  G4String pmt_file = "../input_files/pmtQE.csv";

  ifstream ReadEff;
  G4int effCounter = 0;
  ReadEff.open(pmt_file);

  if(ReadEff.is_open())
  {
    while(!ReadEff.eof())
    {
      ReadEff>>wavelength>>filler>>cath_eff;
      if(ReadEff.eof()){
        break;
      }
      photocath_energy[57-effCounter] = (1240/wavelength)*eV;
      photocath_EFF[57-effCounter] = cath_eff;
      perfect_EFF[57-effCounter] = 1;
      perfect_REFL[57-effCounter] = 0;
      effCounter++;
    }
  }

  else G4cout<<"Error opening file: " <<pmt_file<<G4endl;
  ReadEff.close();
  effCounter--;

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

  G4RotationMatrix* rotationMatrix = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix1 = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* rotationMatrix2 = new G4RotationMatrix(0,0,0);

  G4VisAttributes* detectorAttr = new G4VisAttributes(G4Colour::Grey());
  detectorAttr->SetVisibility(true);
  detectorAttr->SetForceSolid(true);
  ge_log->SetVisAttributes(detectorAttr);

  G4VisAttributes* rodAttr = new G4VisAttributes(G4Colour::Brown());
  rodAttr->SetVisibility(true);
  rodAttr->SetForceSolid(true);
  rod_log->SetVisAttributes(rodAttr);

  G4VisAttributes* plateAttr = new G4VisAttributes(G4Colour::Blue());
  plateAttr->SetVisibility(true);
  plateAttr->SetForceSolid(true);
  plate_log->SetVisAttributes(plateAttr);

  G4VPhysicalVolume* crystal_placement;
  G4VPhysicalVolume* rod_placement;
  G4VPhysicalVolume* tracker_placement;
  G4double tmpR = 46.*mm + 1*mm + 1.5*mm;

  /*
  0 - PMT on base of tile, collimator included.
  */

  /* Middle String */
  for(int i = 0; i < 4; i++){
    fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,-fSiliconPlate_h+i*5*cm),plate_log,"plate",fWLBox,false,i,false);
    crystal_placement = new G4PVPlacement(0, G4ThreeVector(0,0,detector_height+i*5*cm),ge_log,"crystal",fWLBox,false,i,false);
    fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,fSiliconPlate_h+2*detector_height+i*5*cm),plate_log,"top-plate",fWLBox,false,i,false);
    const G4double VertBarAngle = ((G4double) i) * 120.*deg;
    const G4ThreeVector VertBarTranslation (/* x */ tmpR * std::cos(VertBarAngle),  /* y */ tmpR * std::sin(VertBarAngle), /* z */ 0.);
    rod_placement = new G4PVPlacement(0, VertBarTranslation, rod_log,"support rod",fWLBox,false,i,false);

  }

  for(int i = 1; i < 3; i++){
    fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,-fSiliconPlate_h-i*5*cm),plate_log,"plate",fWLBox,false,3+i,false);
    crystal_placement = new G4PVPlacement(0, G4ThreeVector(0,0,detector_height-i*5*cm),ge_log,"crystal",fWLBox,false,3+i,false);
    fPBox = new G4PVPlacement(0, G4ThreeVector(0,0,fSiliconPlate_h+2*detector_height-i*5*cm),plate_log,"top-plate",fWLBox,false,3+i,false);
  }


  /*
  radius = 15*cm
  x = radius * cos(n / num * 2 * PI);
  y = radius * sin(n / num * 2 * PI);
  */
  G4double radius = 15*cm;
  G4double x;
  G4double y;
  int detcounter = 6;
  int rod_counter = 4;
  double PI = 3.1415926535897;

  for(int i=0;i<6;i++){
    x = radius * cos(2 * PI * i/6);
    y = radius * sin(2 * PI * i/6);

    for(int j = 0; j < 4; j++){
      fPBox = new G4PVPlacement(0, G4ThreeVector(x,y,-fSiliconPlate_h+j*5*cm),plate_log,"plate",fWLBox,false,detcounter+i+j,false);
      crystal_placement = new G4PVPlacement(0, G4ThreeVector(x,y,detector_height+j*5*cm),ge_log,"crystal",fWLBox,false,detcounter+i+j,false);
      fPBox = new G4PVPlacement(0, G4ThreeVector(x,y,fSiliconPlate_h+2*detector_height+j*5*cm),plate_log,"top-plate",fWLBox,false,detcounter+i+j,false);
      const G4double VertBarAngle = ((G4double) j) * 120.*deg;
      const G4ThreeVector VertBarTranslation (/* x */ x+tmpR * std::cos(VertBarAngle),  /* y */ y+tmpR * std::sin(VertBarAngle), /* z */ 0.);
      rod_placement = new G4PVPlacement(0, VertBarTranslation, rod_log,"support rod",fWLBox,false,rod_counter+i+j,false);
    }

    for(int j = 1; j < 3; j++){
      fPBox = new G4PVPlacement(0, G4ThreeVector(x,y,-fSiliconPlate_h-j*5*cm),plate_log,"plate",fWLBox,false,3+detcounter+i+j,false);
      crystal_placement = new G4PVPlacement(0, G4ThreeVector(x,y,detector_height-j*5*cm),ge_log,"crystal",fWLBox,false,3+detcounter+i+j,false);
      fPBox = new G4PVPlacement(0, G4ThreeVector(x,y,fSiliconPlate_h+2*detector_height-j*5*cm),plate_log,"top-plate",fWLBox,false,3+detcounter+i+j,false);
    }
    detcounter = detcounter+5;
    rod_counter = rod_counter+3;
  }

  tracker_placement = new G4PVPlacement(0, G4ThreeVector(0,0,0),tracker_log,"tracker",fWLBox,false,0,false);
  return fWPBox;
}
