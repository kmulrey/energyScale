#include "LORASteppingAction.hh"
#include "LORARunAction.hh"
#include "LORAEventAction.hh"


#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "LORADetectorConstruction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"


LORASteppingAction::LORASteppingAction(LORARunAction* run, LORAEventAction* event)
:runAction(run), eventAction(event)
{ }


void LORASteppingAction::UserSteppingAction(const G4Step* aStep)
{
	const G4VPhysicalVolume* prevolume = aStep->GetPreStepPoint()->GetPhysicalVolume();   
   G4double EdepStep = aStep->GetTotalEnergyDeposit();
	//G4double EnonIondepStep = aStep->GetNonIonizingEnergyDeposit();
             
   if (prevolume->GetName()=="scintillator_1" || prevolume->GetName()=="scintillator_2" || prevolume->GetName()=="scintillator_3" || prevolume->GetName()=="scintillator_4") {
      runAction->AddEdep(EdepStep);
      //runAction->AddnonIonEdep(EnonIondepStep);
      eventAction->AddEdep(EdepStep);
      //eventAction->AddnonIonEdep(EnonIondepStep);
   }

}




