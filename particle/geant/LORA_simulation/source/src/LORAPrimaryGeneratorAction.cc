#include "LORAPrimaryGeneratorAction.hh"

#include "Randomize.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

using namespace CLHEP;

LORAPrimaryGeneratorAction::LORAPrimaryGeneratorAction()
{
   G4int n_particle = 1;
   particleGun = new G4ParticleGun(n_particle);
   
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4String particleName;
   particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="mu-")); //mu-
   particleGun->SetParticleEnergy(4.0*GeV);
   particleGun->SetParticlePosition(G4ThreeVector(-10.0*cm, -10.0*cm, 0.375*m));
   particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1)); //vertical
   genRandomDirs=false;
}

LORAPrimaryGeneratorAction::LORAPrimaryGeneratorAction(bool randomdirs, G4double energy=4.0*GeV)
{
   G4int n_particle = 1;
   particleGun = new G4ParticleGun(n_particle);
   
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4String particleName;
   particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="mu-")); //mu-
   particleGun->SetParticleEnergy(energy);
   particleGun->SetParticlePosition(G4ThreeVector(-10.0*cm, -10.0*cm, 0.375*m));
   particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1)); //vertical
   genRandomDirs=randomdirs;
}

LORAPrimaryGeneratorAction::LORAPrimaryGeneratorAction(G4String pname="mu-", G4double zenith=0, G4double azimuth=0, G4double energy=4.0*GeV)
{
   G4int n_particle = 1;
   particleGun = new G4ParticleGun(n_particle);
   
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   particleGun->SetParticleDefinition(particleTable->FindParticle(pname)); //mu-
   particleGun->SetParticleEnergy(energy);
   
   G4double xground=RandFlat::shoot(-75.0,75.0)*cm;
   G4double yground=RandFlat::shoot(-75.0,75.0)*cm;
   
   G4ThreeVector loc(xground+sin(zenith)*cos(azimuth)*100*cm,yground+sin(zenith)*sin(azimuth)*100*cm,(-75+cos(zenith)*100)*cm);
   G4ThreeVector v(-sin(zenith)*cos(azimuth),-sin(zenith)*sin(azimuth),-cos(zenith));

   particleGun->SetParticlePosition(loc);
   particleGun->SetParticleMomentumDirection(v);
   
   genRandomDirs=false;
   
}

LORAPrimaryGeneratorAction::LORAPrimaryGeneratorAction(G4String pname="mu-", G4double zenith=0, G4double azimuth=0, G4double energy=4.0*GeV, double xground=0, double yground=0)
{
   G4int n_particle = 1;
   particleGun = new G4ParticleGun(n_particle);
   
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   particleGun->SetParticleDefinition(particleTable->FindParticle(pname)); //mu-
   particleGun->SetParticleEnergy(energy);
   
   G4ThreeVector loc(xground+sin(zenith)*cos(azimuth)*100*cm,yground+sin(zenith)*sin(azimuth)*100*cm,(-75+cos(zenith)*100)*cm);
   G4ThreeVector v(-sin(zenith)*cos(azimuth),-sin(zenith)*sin(azimuth),-cos(zenith));
   
   particleGun->SetParticlePosition(loc);
   particleGun->SetParticleMomentumDirection(v);
   
   genRandomDirs=false;
   
}


LORAPrimaryGeneratorAction::~LORAPrimaryGeneratorAction()
{
   delete particleGun;
}

void LORAPrimaryGeneratorAction::SetNewParticle(G4String pname="mu-", G4double zenith=0, G4double azimuth=0, G4double energy=4.0*GeV)
{
   // std::cout << "SetNewParticle with 4 parameters " << pname << " " << energy << "\n";
    
   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   particleGun->SetParticleDefinition(particleTable->FindParticle(pname)); //mu-
   particleGun->SetParticleEnergy(energy);
   
   G4double xground=RandFlat::shoot(-75.0,75.0)*cm;
   G4double yground=RandFlat::shoot(-75.0,75.0)*cm;
   
   G4ThreeVector loc(xground+sin(zenith)*cos(azimuth)*100*cm,yground+sin(zenith)*sin(azimuth)*100*cm,(-75+cos(zenith)*100)*cm);
   G4ThreeVector v(-sin(zenith)*cos(azimuth),-sin(zenith)*sin(azimuth),-cos(zenith));
   
   particleGun->SetParticlePosition(loc);
   particleGun->SetParticleMomentumDirection(v);
   
   genRandomDirs=false;   
}

void LORAPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   //Integrated intensity of vertical muons above 1 GeV/c at sea level = 70 m-2 s-1 sr-1
   //I(theta) = I(0) * cos^2(theta)
   //exception: E > E_PION and theta < 70
   
   //dOmega = sin Theta dTheta dPhi
   
   if (genRandomDirs) {
      G4double azimuth=RandFlat::shoot(0.0,6.2831853);
      G4double zenith=RandFlat::shoot(0.0,1.5707963268);
      while (RandFlat::shoot(0.0,1.0)>2*cos(zenith)*cos(zenith)*cos(zenith)*sin(zenith))
      {
         zenith=RandFlat::shoot(0.0,1.5707963268);
         //keep taking random values; its pdf 
      }
   
      //zenith=0;
      G4double xground=RandFlat::shoot(-75.0,75.0)*cm;
      G4double yground=RandFlat::shoot(-75.0,75.0)*cm;

      //This includes the projection effect.
      //We can choose a random ground location on a fixed area
      G4ThreeVector loc(xground+sin(zenith)*cos(azimuth)*100*cm,yground+sin(zenith)*sin(azimuth)*100*cm,(-75+cos(zenith)*100)*cm);
      G4ThreeVector v(-sin(zenith)*cos(azimuth),-sin(zenith)*sin(azimuth),-cos(zenith));
      particleGun->SetParticlePosition(loc);
      particleGun->SetParticleMomentumDirection(v);      
   }
   particleGun->GeneratePrimaryVertex(anEvent);
}
