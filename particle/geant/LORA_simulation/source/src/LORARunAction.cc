#include "LORARunAction.hh"
#include "LORADetectorConstruction.hh"
#include "LORAPrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include <iomanip>


LORARunAction::LORARunAction(LORADetectorConstruction* det, LORAPrimaryGeneratorAction* kin)
:detector(det), primary(kin)
{ }


LORARunAction::~LORARunAction()
{ }


void LORARunAction::BeginOfRunAction(const G4Run* aRun)
{  

  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  edep = 0.;
  nonionedep = 0.; 
}

void LORARunAction::EndOfRunAction(const G4Run* aRun)
{ 

}

