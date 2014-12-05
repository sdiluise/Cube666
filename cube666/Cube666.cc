/**
 *
 *
 *
 *
 *
 * Created on: Dec 5, 2014
 *     Khoi Nguyen:  khoi.nguyen.nguyen@cern.ch
 *     Silvestro di Luise:   Silvestro.Di.Luise@cern.ch
 *
 *              
 *
 *
 *
 *
 */



//#include "preparation.hh"

//#include "CLHEP"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


#include "Cube666_DetectorConstruction.hh"
#include "Cube666_PhysicsList.hh"
#include "Cube666_PrimaryGeneratorAction.hh"
#include "Cube666_RunAction.hh"
#include "Cube666_EventAction.hh"
#include "Cube666_TrackingAction.hh"
#include "Cube666_SteppingAction.hh"
#include "Cube666_Digitizer.hh"

#include "TStopwatch.h"

using namespace std;

#include "Cube666_Analysis.hh"


#define DEBUG true
#if DEBUG
#define D(x) cout<<x<<endl;
#else
#define D(x)
#endif

int main(int argc, char** argv){
  //
  //   TStopwatch timer;
  //   timer.Start();
  //

  //test
  TStopwatch timer;
  timer.Start();
  
  if(argc > 3)
  Cube666_Analysis::fFilename = argv[3];

  G4RunManager* runManager = new G4RunManager; 
  runManager->SetVerboseLevel(0);


//   G4VUserPhysicsList* physics = new Cube666_PhysicsList(turnOnRayleigh,turnOnBoundaryProcess,
// 						     turnOnWLSProcess,turnOnScintProcess); 
  G4VUserPhysicsList* physics = new Cube666_PhysicsList();


  physics->SetVerboseLevel(0);

  D("--Init Physics--"); runManager->SetUserInitialization(physics);

  D("--Detector Constr--"); runManager->SetUserInitialization(new Cube666_DetectorConstruction); 
  D("--Primary Gen Action--"); runManager->SetUserAction(new Cube666_PrimaryGeneratorAction);
  D("--Run Action --"); runManager->SetUserAction(new Cube666_RunAction);
  D("--Event Action--"); runManager->SetUserAction(new Cube666_EventAction);
  D("--Tracking Action--"); runManager->SetUserAction(new Cube666_TrackingAction);
  D("--Stepping Action--"); runManager->SetUserAction(new Cube666_SteppingAction);

  D("----Settings end -----");



  /////////////////////////////////
  /////////////TEST////////////////
  /////////////////////////////////

  

  G4int visualization=0;
  if(argv[1]) visualization = atoi(argv[1]);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4VisManager* visManager = new G4VisExecutive;
  G4UIExecutive* ui = new G4UIExecutive(argc,argv);

  //cout<<"visManager->GetCurrentViewer() "<<visManager->GetCurrentViewer()<<endl; getchar();
  if(visualization){

    D(" ---- Visualization ----");
    visManager->Initialize();
    //UImanager->ApplyCommand("/vis/open HepRep 800x600-0+0");
    UImanager->ApplyCommand("/vis/open OGL 800x600-0+0");

    UImanager->ApplyCommand("/vis/viewer/set/style surface");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 89.9 -90.");

    //UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 20 40.");
    UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 89 -90");

    UImanager->ApplyCommand("/run/initialize");
    UImanager->ApplyCommand("/vis/drawVolume");
    UImanager->ApplyCommand("/vis/scene/add/trajectories");
    UImanager->ApplyCommand("/vis/scene/add/hits");
    UImanager->ApplyCommand("/vis/geometry/set/colour pmt ! green");
    UImanager->ApplyCommand("/vis/geometry/set/colour pmtCoat ! red");


    D(" ---- Visualization end ----");
  }




  //set random seed for random number generator
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  // G4long seed = time(0);
  // CLHEP::HepRandom::setTheSeed(seed);


  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed = time(0);
  G4Random::setTheSeed(seed);



  /////////////////////////////////
  ///////////END TEST//////////////
  /////////////////////////////////
    
  if(VERBOSE) runManager->SetVerboseLevel(4);

  D("--Initialize--"); runManager->Initialize(); D("--- Initialize End -- ");

  G4String cmd = "/control/execute ";
  G4String macro = argv[2];

  D("--Apply Command Start --"); 
  
  TStopwatch time_run; time_run.Start();

  UImanager->ApplyCommand(cmd + macro); 

  time_run.Stop();

  cout<<"\n \t BeamOn time: "
      <<"CPU time "<<time_run.CpuTime()<<"\t real time "<<time_run.RealTime()<<endl;

  D("--- Apply Command End -- ");

  //end test


  //runManager->BeamOn(nevents);
  if(visualization) ui->SessionStart();

  delete ui;
  delete visManager;
  delete runManager;

  timer.Stop();
  cout<<"finished running "
      <<"\t CPU time "<<timer.CpuTime()<<"\t real time "<<timer.RealTime()<<endl;

  return 0;
}

