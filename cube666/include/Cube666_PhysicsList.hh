#ifndef _ARDM_PHYSICSLIST_
#define _ARDM_PHYSICSLIST_ 1



#include "G4VUserPhysicsList.hh"
//#include "QGSP_BIC_EMY.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"

class Cube666_PhysicsList : public G4VUserPhysicsList{

private: 

  void electronPhysics(G4ProcessManager* pmanager);          //e-
  void positronPhysics(G4ProcessManager* pmanager);          //e+
  void muonPhysics(G4ProcessManager* pmanager);           //mu-
  void opPhotonPhysics(G4ProcessManager* pmanager);          //optical photon
  void gammaPhysics(G4ProcessManager* pmanager);             //high energy photon
  void protonPhysics(G4ProcessManager* pmanager);            //proton
  void neutronPhysics(G4ProcessManager* pmanager);           //neutron
  void ionPhysics(G4ProcessManager* pmanager);               //ions : alpha,deuteron,triton,He3,generic ion
  void hadronPhysics(G4ProcessManager* pmanager);            //proton, neutron, pion, kaon, ...
  void addScintillation();                                       //add scint. process to all particles except for opPhoton
  //void otherPhysics(G4ProcessManager* pmanager);             //other physics ??

public:
  Cube666_PhysicsList();
  ~Cube666_PhysicsList();

  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();


//   void ConstructOpticalProcess();
//   void ConstructEMProcess();
};






#endif // _ARDM_PHYSICSLIST_
