#ifndef _ARDM_RUNACTION_
#define _ARDM_RUNACTION_ 1


#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "Cube666_Analysis.hh"


class Cube666_RunAction : public G4UserRunAction{

private:
  Cube666_Analysis* fAna;

public:
  Cube666_RunAction();
  ~Cube666_RunAction();

  void BeginOfRunAction(const G4Run* run);
  void EndOfRunAction(const G4Run* run);

};


#endif //_ARDM_RUNACTION_
