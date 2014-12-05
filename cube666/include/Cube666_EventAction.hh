#ifndef _ARDM_EVENTACTION_
#define _ARDM_EVENTACTION_ 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"


class Cube666_EventAction : public G4UserEventAction{

public:
  Cube666_EventAction();
  ~Cube666_EventAction();

  void BeginOfEventAction(const G4Event* event);
  void EndOfEventAction(const G4Event* event);

};




#endif // _ARDM_EVENTACTION_ 
