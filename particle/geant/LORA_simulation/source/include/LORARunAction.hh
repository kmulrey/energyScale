#ifndef LORARunAction_h
#define LORARunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class LORADetectorConstruction;
class LORAPrimaryGeneratorAction;

class LORARunAction : public G4UserRunAction
{
  public:
    LORARunAction(LORADetectorConstruction*, LORAPrimaryGeneratorAction*);
   ~LORARunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
        
    void AddEdep(G4double val)     { edep += val;}
    void AddnonIonEdep(G4double val)     { nonionedep += val;}
                                  
  private:
    G4double        edep, nonionedep;
    LORADetectorConstruction*   detector;
    LORAPrimaryGeneratorAction* primary;
};

#endif

