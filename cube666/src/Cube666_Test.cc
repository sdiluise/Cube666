#include "Cube666_Test.hh"
#include "preparation.hh"
#include "Cube666_Analysis.hh"

Cube666_Test* Cube666_Test::fInstance = 0;

Cube666_Test::Cube666_Test(){
  reset();
}


Cube666_Test* Cube666_Test::getInstance(){

  if(!fInstance) fInstance = new Cube666_Test;
  return fInstance;
}


void Cube666_Test::reset(){
  fNDirectPhoton = 0;
  reset1DArray<G4int>(fTotNDirectPhotonPMT,NPMT,0);

  fNPhotonHittingPMT = 0;
  reset1DArray<G4int>(fTotNPhotonHittingPMT,NPMT,0);

  return;
}


