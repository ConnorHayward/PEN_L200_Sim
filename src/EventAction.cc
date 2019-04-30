#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* myEvent)
{
	fDetectedPhotons = 0;
	fDepositedEnergy = 0;
	fProducedPhotons = 0;
	fTopPhoton =0;
	fBottomPhoton=0;
	fSidePhoton=0;
	fWavelength=0;
	fEscapedPhoton = 0;
	fAbsorbedPhoton=0;
	if (myEvent->GetEventID() % 1000 == 0)
		G4cout << "event no.: " << myEvent->GetEventID() << G4endl;
}

void EventAction::EndOfEventAction(const G4Event* myEvent)
{

	auto analysisManager = G4AnalysisManager::Instance();

	 if(fDetectedPhotons > 0){
 		analysisManager->FillNtupleDColumn(0,0,1);
	 	analysisManager->FillNtupleDColumn(0,1,fWavelength);
	// 	analysisManager->FillNtupleDColumn(0,2,fSidePhoton);
	// 	analysisManager->FillNtupleDColumn(0,3,fDetectedPhotons);
	// 	analysisManager->FillNtupleDColumn(0,4,fDepositedEnergy);
	 	analysisManager->AddNtupleRow(0);
	 }
}
