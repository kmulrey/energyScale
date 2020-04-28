#ifndef LORAPrimaryGeneratorAction_h
#define LORAPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class LORAPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
   LORAPrimaryGeneratorAction();
   LORAPrimaryGeneratorAction(G4bool, G4double);
   LORAPrimaryGeneratorAction(G4String, G4double, G4double, G4double);
   LORAPrimaryGeneratorAction(G4String, G4double, G4double, G4double, G4double, G4double);
   ~LORAPrimaryGeneratorAction();
   
public:
   void SetNewParticle(G4String, G4double, G4double, G4double);
   void GeneratePrimaries(G4Event* anEvent);
	
	G4ParticleGun* GetParticleGun() {return particleGun;}
   
private:
   G4ParticleGun* particleGun;
   G4bool genRandomDirs;
};

#endif

