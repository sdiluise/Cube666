#ifndef _SENSITIVEPMT_
#define _SENSITIVEPMT_ 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Cube666_PMTHit.hh"
#include "G4VSensitiveDetector.hh"
#include "Cube666_Digitizer.hh"

#include "vector"
using namespace std;


class Cube666_SensitivePMT : public G4VSensitiveDetector{

protected:
  G4int    fnpmt;
  G4String fDetName;
  G4String fHitsColName;
  
  G4ThreeVector fNormal;
  vector<G4ThreeVector> fPMTvector;

  Cube666_PMTHitsCollection* fHitsCol;

  static const G4double fTimeSample;
  static const G4int    fNTimeSamples;

  G4double detectionProb(G4double phi);
  G4double detectionProb(G4ThreeVector hitpos,G4ThreeVector rPMT);
  G4double detectionProb(G4ThreeVector hitpos,G4int pmtID);

  G4int detected(G4ThreeVector hitpos,G4ThreeVector rPMT);
  G4int detected(G4ThreeVector hitpos,G4int pmtID);

  static Cube666_Digitizer* digitizer;

public:

  Cube666_SensitivePMT(G4int npmt,G4String detName,G4String hitsColName);
  ~Cube666_SensitivePMT();

  void Initialize(G4HCofThisEvent* hc);
  G4bool ProcessHits(G4Step* astep, G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* hc);
  
};




#endif //_SENSITIVEPMT_
