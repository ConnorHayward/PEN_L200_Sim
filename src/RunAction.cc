// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"

#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "EventAction.hh"

#include <string>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
 : G4UserRunAction(),fDetector(det),fPrimary(kin),
   fTimer(0)
{
  fTimer = new G4Timer;
  fMan = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fTimer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  fMan = G4AnalysisManager::Instance();

  G4String targetName = fDetector->GetTargetMaterialName();
  G4String targetThickness = G4BestUnit(fDetector->GetTargetSize(),"Length");
  G4int sourceName = fPrimary->GetSourceType();
  G4int directory = fDetector->GetDetectorType();
  G4String directorName;

  G4String position = G4BestUnit((fDetector->GetTargetSize()-fPrimary->GetPosition()),"Length");


  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
  fFileName = "../Data/"+targetName+"_test"+"_Plates.csv";
  fMan->SetVerboseLevel(0);
  fMan->OpenFile(fFileName);

 // fMan->CreateNtuple("EscapedDemo","Demo");
 // fMan->CreateNtupleDColumn("Hits");
 // fMan->CreateNtupleDColumn("Wavelength");
 // fMan->FinishNtuple();

 fMan->CreateNtuple("Ge_Energy","Test");
 fMan->CreateNtupleDColumn("Event No");
 fMan->CreateNtupleSColumn("Isotope");
 fMan->CreateNtupleSColumn("Material");
 fMan->CreateNtupleSColumn("Volume");
 fMan->CreateNtupleDColumn("Energy [MeV]");
 fMan->CreateNtupleDColumn("Veto");
 fMan->FinishNtuple();

 fMan->CreateNtuple("Ge_Energy_lAr","Test");
 fMan->CreateNtupleDColumn("Energy [MeV]");
 fMan->FinishNtuple();
}

void RunAction::SetFileName(G4String fileName)
{
  fFileName = fileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;
    G4cout << "End Run" << G4endl;
    fMan->Write();
    fMan->CloseFile();
    delete fMan;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
