#ifndef _ARDM_TRACKINGACTION_
#define _ARDM_TRACKINGACTION_ 1

#include "G4UserTrackingAction.hh"
#include "G4Track.hh"

class Cube666_TrackingAction : public G4UserTrackingAction{
  
  void verboseInfo(const G4Track* track);

public:
  Cube666_TrackingAction();
  ~Cube666_TrackingAction();
  
  void PreUserTrackingAction(const G4Track* track);
  void PostUserTrackingAction(const G4Track* track);

};



#endif // _ARDM_TRACKINGACTION_
