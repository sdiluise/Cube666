#ifndef _ARDM_TEST_
#define _ARDM_TEST_ 1

#include "preparation.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"

#include "vector"

using namespace std;

class Cube666_Test{

private:
  
  Cube666_Test();
  static Cube666_Test* fInstance;

public:
  ~Cube666_Test(){;};
  static Cube666_Test* getInstance();
  
  void reset();

  G4int fNDirectPhoton;
  G4int fTotNDirectPhotonPMT[NPMT];

  G4int fNPhotonHittingPMT;
  G4int fTotNPhotonHittingPMT[NPMT];
};


#endif
