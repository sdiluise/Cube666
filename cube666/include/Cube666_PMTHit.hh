#ifndef _ARDM_PMTHIT_
#define _ARDM_PMTHIT_ 1

#include "G4VHit.hh"

#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"


class Cube666_PMTHit : public G4VHit{

private:
  G4ThreeVector fHitPos;
  G4double      fHitGlobalTime;
  G4double      fHitEnergy;
  G4int         fHitPMTid;  // <--> which PMT does the photon hit ?
  G4int         fIsDetected;
  
public:
  Cube666_PMTHit();
  ~Cube666_PMTHit();

  inline void setHitPos(G4ThreeVector hitpos) { fHitPos = hitpos;};
  inline G4ThreeVector getHitPos() { return fHitPos;};

  inline void setHitGlobalTime(G4double globalTime) { fHitGlobalTime = globalTime;};
  inline G4double getHitGlobalTime() { return fHitGlobalTime; }

  inline void setHitEnergy(G4double hitEnergy) { fHitEnergy = hitEnergy;};
  inline G4double getHitEnergy() { return fHitEnergy; }

  inline void setHitPMTid(G4int pmtid){ fHitPMTid = pmtid;};
  inline G4int getHitPMTid(){return fHitPMTid;};

  inline void setIsDetected(G4int isDetected){ fIsDetected = isDetected;};
  inline G4int getIsDetected(){ return fIsDetected;};
};


typedef G4THitsCollection<Cube666_PMTHit> Cube666_PMTHitsCollection;







#endif // _ARDM_PMTHIT_
