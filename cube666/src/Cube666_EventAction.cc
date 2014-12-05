#include "Cube666_EventAction.hh"
#include "Cube666_Analysis.hh"

#include "Cube666_Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "Cube666_PrimaryGeneratorAction.hh"
#include "G4DigiManager.hh"
#include "Cube666_PMTHit.hh"


Cube666_EventAction::Cube666_EventAction()
  : G4UserEventAction(){;}


Cube666_EventAction::~Cube666_EventAction(){;}


void Cube666_EventAction::BeginOfEventAction(const G4Event* event){
 
 if(!(event->GetEventID()%1)){ G4cout<<" ====> Begin of event "<<event->GetEventID()<<G4endl; //getchar();
 }
  
  Cube666_Analysis::getInstance()->BeginOfEvent();
  return;
}


void Cube666_EventAction::EndOfEventAction(const G4Event* event){

  Cube666_Analysis::getInstance()->EndOfEvent();

  if(!(event->GetEventID()%1)){ G4cout<<" ====> End of event "<<event->GetEventID()<<G4endl; //getchar();
  }

  return;
}
