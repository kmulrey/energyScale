
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "LORADetectorConstruction.hh"
#include "LORAPhysicsList.hh"
#include "LORAPrimaryGeneratorAction.hh"

#include "LORARunAction.hh"
#include "LORAEventAction.hh"
#include "LORASteppingAction.hh"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc,char** argv)
{
   // Construct the default run manager
   //
   G4RunManager* runManager = new G4RunManager;
   
   // set mandatory initialization classes
   //
   LORADetectorConstruction* detector = new LORADetectorConstruction;
   runManager->SetUserInitialization(detector);
   //
   G4VUserPhysicsList* physics = new LORAPhysicsList;
   runManager->SetUserInitialization(physics);
	
   // set mandatory user action class
   //
   cout << "1" << endl;
   //LORAPrimaryGeneratorAction* gen_action = new LORAPrimaryGeneratorAction(false, 4.0*GeV);
   LORAPrimaryGeneratorAction* gen_action = new LORAPrimaryGeneratorAction();
   //LORAPrimaryGeneratorAction* gen_action;

   cout << "2" << endl;

   runManager->SetUserAction(gen_action);
	
	// set user action classes
	LORARunAction*   run;
	LORAEventAction* event;
	
	runManager->SetUserAction(run = new LORARunAction(detector,gen_action)); 
	runManager->SetUserAction(event = new LORAEventAction());
	//runManager->SetUserAction(new VGCTrackingAction(gen_action,run,histo));
	runManager->SetUserAction(new LORASteppingAction(run, event));
   
   // Initialize G4 kernel
   //
   runManager->Initialize();
   
   // Get the pointer to the UI manager and set verbosities
   //
   G4UImanager* UI = G4UImanager::GetUIpointer();
   
   UI->ApplyCommand("/control/verbose 0");
	UI->ApplyCommand("/run/verbose 0"); 
	UI->ApplyCommand("/event/verbose 0"); 
	UI->ApplyCommand("/tracking/verbose 0");
   UI->ApplyCommand("/process/verbose 0");
   
	if (argc==3)   // batch mode  
   {
      char* inFile = argv[1];
      char* outFile = argv[2];
      
      G4String line;
      ifstream myfile(inFile);
      double penergy, zenith, azimuth;
      if (myfile.is_open())
      {
         getline (myfile,line);
         istringstream iss(line);
         int evtno, pid;
         
         iss >> evtno >> pid >> penergy >> zenith >> azimuth;
         
         double shower_axis_x = sin(zenith)*cos(azimuth);
         double shower_axis_y = sin(zenith)*sin(azimuth);
         
         const int nobins=200;
         
         double TotalArea[nobins];
         double DepositedEnergyTotal[nobins];
         double DepositedEnergyGamma[nobins];
         double DepositedEnergyElecPosi[nobins];
         double DepositedEnergyMuons[nobins];
         double DepositedEnergyHadrons[nobins];
         for (int i=0; i<nobins; i++) {
            TotalArea[i]=25*3.14159265*((i+1)*(i+1)-i*i); //m^2
            DepositedEnergyTotal[i]=0;
            DepositedEnergyGamma[i]=0;
            DepositedEnergyElecPosi[i]=0;
            DepositedEnergyMuons[i]=0;
            DepositedEnergyHadrons[i]=0;
         }
         while ( getline (myfile,line) )
         //for (int i=0; i<1000; i++)
         {
            //getline (myfile,line);
            istringstream iss(line);
            int cors_id, id;
            double px, py, pz, x, y, t, w;
            iss >> cors_id >> id >> px >> py >> pz >> x >> y >> t >> w;
            
            double R2 = x*x+y*y;
            double inner = x*shower_axis_x + y*shower_axis_y;
            double AxisDist = sqrt(R2-inner*inner);
            int DistBin = AxisDist/500; // bin 0 = 0-5 m, bin 1 = 5-10m, etc.
            double energy = sqrt(px*px+py*py+pz*pz)*GeV;
            
            G4String nextParticle="ignore";
            if (id==1) {nextParticle="gamma";}
            if (id==2) {nextParticle="e+";}
            if (id==3) {nextParticle="e-";}
            if (id==5) {nextParticle="mu+";}
            if (id==6) {nextParticle="mu-";}
            if (id==13) {nextParticle="neutron";}
            if (id==14) {nextParticle="proton";}
            
            if (nextParticle!="ignore" && DistBin<nobins) {
               //cout << nextParticle << " at " << AxisDist/100. << " m (bin " << DistBin << "), " << energy << " GeV" << endl;
               // GeV BUG !!! : gen_action->SetNewParticle(nextParticle,zenith,azimuth,energy*GeV);
               gen_action->SetNewParticle(nextParticle,zenith,azimuth,energy);
               runManager->BeamOn(1); //this gives the particle a random position in an area of 2.25 m^2
               double Deposit = w*(event->GetEnergyDeposit());    // what is w ? Weight, from thinning? (AC)
               //G4cout << "Deposit: " << Deposit/w << " (weight=" << w << ")" << G4endl;
               DepositedEnergyTotal[DistBin]+=Deposit;
                  if (id==1) { DepositedEnergyGamma[DistBin]+=Deposit; }
                  else if (id==2 || id==3) { DepositedEnergyElecPosi[DistBin]+=Deposit; }
                  else if (id==5 || id==6) { DepositedEnergyMuons[DistBin]+=Deposit; }
                  else if (id==13 || id==14) { DepositedEnergyHadrons[DistBin]+=Deposit; }
            }
         }
         myfile.close();
         //cout << "Bin   Total (Gev)  Gamma(%)   ElecPosi(%) Muons(%) Hadrons(%)" << endl;
         ofstream myfile(outFile);
         
         for (int i=0; i<nobins; i++) {
            //We put all particles in a distance ring into an area of 2.25 m^2
            //Correct particle density is achieved by applying a factor 2.25/A_ellipse
            //A "ring" is in fact an ellipse, which has Area A_ellipse=A_ring/cos(zenith)
            
            // (AC): What is actually written out here??
            // Deposited energy in MeV? Then the 2.25 cos theta / totalarea doesnt make sense?
            // or energy density in MeV / m^2 ? It seems MeV / m^2, I would say, as then
            // the area correction factor is right.
            // The comparison is with LORA MeV in detector, however...
             
            myfile << i << "     " << 2.25*cos(zenith)*DepositedEnergyTotal[i]/TotalArea[i] << "    " << int(100.*DepositedEnergyGamma[i]/DepositedEnergyTotal[i])
               << "   " << int(100.*DepositedEnergyElecPosi[i]/DepositedEnergyTotal[i]) 
               << "   " << int(100.*DepositedEnergyMuons[i]/DepositedEnergyTotal[i])
            << "   " << int(100.*DepositedEnergyHadrons[i]/DepositedEnergyTotal[i]) << endl;
         
         }
         myfile.close();
      }

   }

   delete runManager;
   
   return 0;
}

