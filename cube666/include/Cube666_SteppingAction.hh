#ifndef _ARDM_STEPPING_ACTION_
#define _ARDM_STEPPING_ACTION_ 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

class Cube666_SteppingAction : public G4UserSteppingAction{

  void ioniProc_DoIt(const G4Step* step,G4String particleName = "e-");
  void scattProc_DoIt(const G4Step* step);
  void boundaryProc_DoIt(const G4Step* step);
  void WLSProc_DoIt(const G4Step* step);
  void scintillationProc_DoIt(const G4Step* step);
  void bremsstrahlungProc_DoIt(const G4Step* step,G4String particleName="e-");
  void multiScattProc_DoIt(const G4Step* step,G4String particleName="e-");
  void elasticScattProc_DoIt(const G4Step* step,G4String particleName="neutron");
  void inelasticScattProc_DoIt(const G4Step* step,G4String particleName="neutron");
  void neutronCaptureProc_DoIt(const G4Step* step,G4String particleName="neutron");


  void opPhoton_DoIt(const G4Step* step);

  void gamma_DoIt(const G4Step* step);

  void neutron_DoIt(const G4Step* step); //neutron interaction in makrolon

  void neutronShield_DoIt(const G4Step* step);  //test neutron shielding in LAr, makrolon, real ArDM geometry
  void neutronShield2_DoIt(const G4Step* step); //test neutron shield alone, only shield, independent of ArDM
  void neutronShield3_DoIt(const G4Step* step); //test neutron shield in LAr, only shield, independent of ArDM

  void neutron_Erecoil(const G4Step* step);
  void neutronInteraction_DoIt(const G4Step* step);
  void neutronBackground_DoIt(const G4Step* step);

  void photoElectricProc_DoIt(const G4Step* step);
  void comptonProc_DoIt(const G4Step* step);
  void pairProductionProc_DoIt(const G4Step* step);
 
  void stepVerboseInfo(const G4Step* step);
  void stepVerboseInfo_secondaries(const G4Step* step);
  void printTrackInfo(G4Track* track,const G4Step* step);
  void printParentInfo(const G4Step* step);

  G4int isParticle(const char* particleName,const G4Step* step);
  G4int isProcess(const char* processName,const G4Step* step);
  G4int isVolume(const char* volumeName,const G4Step* step);
  G4int isNextVolume(const char* nextVolumeName,const G4Step* step);

public:
  Cube666_SteppingAction();
  ~Cube666_SteppingAction();

  void UserSteppingAction(const G4Step* step);  

};

#endif // _ARDM_STEPPING_ACTION__
