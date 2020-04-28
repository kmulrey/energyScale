#include "globals.hh"
#include "LORAPhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

LORAPhysicsList::LORAPhysicsList()
{
	defaultCutValue = 1.0*cm;
	SetVerboseLevel(1);
}

LORAPhysicsList::~LORAPhysicsList()
{;}

void LORAPhysicsList::ConstructParticle()
{
   
	ConstructBosons();
	ConstructLeptons();
	ConstructMesons();
	ConstructBaryons();
	
}

void LORAPhysicsList::ConstructBosons()
{
	// gamma
	G4Gamma::GammaDefinition();
}

void LORAPhysicsList::ConstructLeptons()
{
	// leptons
	//  e+/-
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	// mu+/-
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();
	// nu_e
	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	// nu_mu
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void LORAPhysicsList::ConstructMesons()
{
	//  mesons
	//    light mesons
	G4PionPlus::PionPlusDefinition();
	G4PionMinus::PionMinusDefinition();
	G4PionZero::PionZeroDefinition();
	G4Eta::EtaDefinition();
	G4EtaPrime::EtaPrimeDefinition();
	G4KaonPlus::KaonPlusDefinition();
	G4KaonMinus::KaonMinusDefinition();
	G4KaonZero::KaonZeroDefinition();
	G4AntiKaonZero::AntiKaonZeroDefinition();
	G4KaonZeroLong::KaonZeroLongDefinition();
	G4KaonZeroShort::KaonZeroShortDefinition();
}

void LORAPhysicsList::ConstructBaryons()
{
	//  baryons
	G4Proton::ProtonDefinition();
	G4AntiProton::AntiProtonDefinition();
	
	G4Neutron::NeutronDefinition();
	G4AntiNeutron::AntiNeutronDefinition();
}

void LORAPhysicsList::ConstructProcess()
{
   // Define transportation process
   
   AddTransportation();
   ConstructEM();
   ConstructGeneral();
   AddStepMax();	
}

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh" //needed?

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh" // needed?

#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh" // To get full atomic relaxtation model

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh" // needed?

#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh" //needed?

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh" // low energy ionization

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh" // needed?

#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

void LORAPhysicsList::ConstructEM()
{
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		
		if (particleName == "gamma") {
			// gamma         
			G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
			G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePhotoElectricModel();
			thePhotoElectricEffect->SetEmModel(theLivermorePhotoElectricModel);
			pmanager->AddDiscreteProcess(thePhotoElectricEffect);
			
			G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
			G4LivermoreComptonModel* theLivermoreComptonModel = new G4LivermoreComptonModel();
			theComptonScattering->SetEmModel(theLivermoreComptonModel);
			pmanager->AddDiscreteProcess(theComptonScattering);
			
			G4GammaConversion* theGammaConversion = new G4GammaConversion();
			G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermoreGammaConversionModel();
			theGammaConversion->SetEmModel(theLivermoreGammaConversionModel);
			pmanager->AddDiscreteProcess(theGammaConversion);
			
			G4RayleighScattering* theRayleigh = new G4RayleighScattering();
			G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
			theRayleigh->SetEmModel(theRayleighModel);
			pmanager->AddDiscreteProcess(theRayleigh);
			
			//pmanager->AddProcess(new G4StepLimiter(), -1, -1, 5); //What does this mean?
			
		} else if (particleName == "e-") {
			//electron
			G4eMultipleScattering* msc = new G4eMultipleScattering();
			msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,                   -1, 1, 1);
			
			// Ionisation
			G4eIonisation* eIoni = new G4eIonisation();
			G4LivermoreIonisationModel* theLivermoreIonisationModel = new G4LivermoreIonisationModel();
			eIoni->SetEmModel(theLivermoreIonisationModel);
			eIoni->SetStepFunction(0.2, 100*um); //     
			pmanager->AddProcess(eIoni,                 -1, 2, 2);
			
			// Bremsstrahlung
			G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
			eBrem->SetEmModel(new G4LivermoreBremsstrahlungModel());
			pmanager->AddProcess(eBrem,         -1,-3, 3);
			
			//pmanager->AddProcess(new G4StepLimiter(), -1, -1, 4);
			
		} else if (particleName == "e+") {
			//positron
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
			pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
			pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);
			
		} else if( particleName == "mu+" || 
                particleName == "mu-"    ) {
			//muon  
			pmanager->AddProcess(new G4MuMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4MuIonisation,         -1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung,     -1, 3, 3);
			pmanager->AddProcess(new G4MuPairProduction,     -1, 4, 4);       
			
		} else if( particleName == "proton" ||
                particleName == "pi-" ||
                particleName == "pi+"    ) {
			//proton  
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
			pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
			pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);       
			
		} else if( particleName == "alpha" || 
                particleName == "He3" )     {
			//alpha 
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
			
		} else if( particleName == "GenericIon" ) { 
			//Ions 
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
			
		} else if ((!particle->IsShortLived()) &&
                 (particle->GetPDGCharge() != 0.0) && 
                 (particle->GetParticleName() != "chargedgeantino")) {
			//all others charged particles except geantino
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
		}
	}
}

#include "G4Decay.hh"

void LORAPhysicsList::ConstructGeneral()
{
	// Add Decay Process
	G4Decay* theDecayProcess = new G4Decay();
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		if (theDecayProcess->IsApplicable(*particle)) { 
			pmanager ->AddProcess(theDecayProcess);
			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
			pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
		}
	}
}

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void LORAPhysicsList::AddStepMax()
{
	// Step limitation seen as a process
	G4StepLimiter* stepLimiter = new G4StepLimiter();
	////G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();
	
	theParticleIterator->reset();
	while ((*theParticleIterator)()){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		
		if (particle->GetPDGCharge() != 0.0)
      {
			pmanager ->AddDiscreteProcess(stepLimiter);
			////pmanager ->AddDiscreteProcess(userCuts);
      }
	}
}

void LORAPhysicsList::SetCuts()
{
	//G4VUserPhysicsList::SetCutsWithDefault method sets 
	//the default cut value for all particle types 
	//
	//SetCutsWithDefault();
	defaultCutValue = 1.0*mm;
	SetCutValue(0.01*mm, "e-");
	//G4Electron::GetProductionCutsTable()->SetEnergyRange(10*eV,100*MeV);
	if (verboseLevel>0) DumpCutValuesTable();
}

