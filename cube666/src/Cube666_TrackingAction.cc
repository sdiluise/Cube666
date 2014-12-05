#include "Cube666_TrackingAction.hh"
#include "Cube666_Analysis.hh"
#include "preparation.hh"
#include "Cube666_Test.hh"


#include "G4TrackVector.hh"
#include "G4TrackingManager.hh"

#include "Cube666_Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"


#define DEBUG false
#if DEBUG
#define D(x) cout<<x<<endl;
#else
#define D(x)
#endif


Cube666_TrackingAction::Cube666_TrackingAction()
  : G4UserTrackingAction(){;}


Cube666_TrackingAction::~Cube666_TrackingAction(){;}


void Cube666_TrackingAction::PreUserTrackingAction(const G4Track* track){

  D(__FILE__<<"::"<<__FUNCTION__<<" track: "<<track);

  //consider primary particles only (ID==0)
  //if(track->GetParentID()) return;

  Cube666_Analysis* ana = Cube666_Analysis::getInstance();

  //#if VERBOSE
  verboseInfo(track);
  //#endif //VERBOSE



#if TEST_BRANCH0
#if TEST_BRANCH00

  //if(track->GetTrackID() == 1 && track->GetParentID() == 0){
  if(track->GetParentID() == 0){

    G4ThreeVector pos = track->GetPosition();

    ana->fPosx=pos.x()/mm;
    ana->fPosy=pos.y()/mm;
    ana->fPosz=pos.z()/mm;
  }

#endif //TEST_BRANCH00
#endif //TEST_BRANCH0

  
  return;
}

void Cube666_TrackingAction::PostUserTrackingAction(const G4Track* track){

#if VERBOSE
    cout<<"Cube666_TrackingAction::PostUserTrackingAction(..) ... position 1"<<endl;
#endif //VERBOSE

  return;
}




void Cube666_TrackingAction::verboseInfo(const G4Track* track){

  //G4Track *tr = (G4Track*)track->GetPosition();

  //G4ThreeVector pos = (G4Track*)tr->GetPosition();

  cout<<__FUNCTION__<<" : "<<track<<endl;

  cout<<"trackInfo: "
      <<" "<<track->GetParticleDefinition()->GetParticleName()
    //<<"\t particleID "<<track->GetParticleDefinition()->GetPDGEncoding()
      <<"\t charge "<<track->GetParticleDefinition()->GetPDGCharge()
      <<"\t parentID "<<track->GetParentID()
      <<"\t Ekin [MeV] "<<track->GetDynamicParticle()->GetKineticEnergy()/MeV
      <<endl;
  //if(track->GetDynamicParticle()->GetKineticEnergy()/eV < 1) track->SetTrackStatus(fKillTrackAndSecondaries);    

  //if(track->GetParentID()){
    cout<<"    position: "
	<<" "<<track->GetPosition().x()/mm
	<<" "<<track->GetPosition().y()/mm
	<<" "<<track->GetPosition().z()/mm
	<<" [mm]"
	<<endl;
  //}

  if(track->GetParticleDefinition()->GetParticleName() == "e-"){
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    cout<<" secondaries* "<<secondaries<<endl;
    if(secondaries) cout<<"number of secondary particles : "<<secondaries->size()<<endl;;
  }

  cout<<__FUNCTION__<<" -- end --"<<endl;

  return;
}




