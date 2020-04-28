#include "LORAEventAction.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

LORAEventAction::LORAEventAction()
:drawFlag("none"),printModulo(10000)
{

}


LORAEventAction::~LORAEventAction()
{
  //delete eventMessenger;
}


void LORAEventAction::BeginOfEventAction(const G4Event* evt)
{          
 TotalEnergyDeposit = 0.;
 TotalNonIonizingEnergyDeposit = 0.;
}


void LORAEventAction::EndOfEventAction(const G4Event* evt)
{

}




