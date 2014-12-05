#include "Cube666_RunAction.hh"

Cube666_RunAction::Cube666_RunAction()
  :G4UserRunAction(),fAna(NULL){}

Cube666_RunAction::~Cube666_RunAction(){
  if(fAna) delete fAna;
}



void Cube666_RunAction::BeginOfRunAction(const G4Run* run){

  G4cout<<"-----> Starting Run: "<<run->GetRunID()<<G4endl;

  fAna = Cube666_Analysis::getInstance();
  fAna->BeginOfRun();

  return;
}


void Cube666_RunAction::EndOfRunAction(const G4Run* run){
  fAna->EndOfRun();
  G4cout<<"finished processing run "<<run->GetRunID()<<G4endl;
  G4cout<<"output stored in "<< Cube666_Analysis::fFilename <<G4endl;
  return;
}
