#include "Cube666_SteppingAction.hh"
#include "Cube666_Analysis.hh"
#include "G4StepPoint.hh"
#include "Cube666_Test.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4DynamicParticle.hh"

#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "Cube666_Trajectory.hh"
#include "G4TrackVector.hh"

#define DEBUG false
#if DEBUG
#define D(x) cout<<x<<endl;
#else
#define D(x)
#endif


Cube666_SteppingAction::Cube666_SteppingAction()
  :G4UserSteppingAction(){;}

Cube666_SteppingAction::~Cube666_SteppingAction(){;}


void Cube666_SteppingAction::UserSteppingAction(const G4Step* step){

  D(" "<<__FUNCTION__);

  G4StepPoint* pPreStepPoint  = step->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = step->GetPostStepPoint();
  
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = step->GetDeltaPosition().unit();
  G4double      t0 = pPreStepPoint->GetGlobalTime();

  //cout<<"22222 "<<x0.x()/mm<<" "<<x0.y()/mm<<" "<<x0.z()/mm<<endl;

  const char* checkParticle = "alpha";
  const char* checkProcess  = "Ioni";

  Cube666_Analysis* ana = Cube666_Analysis::getInstance();

  scattProc_DoIt(step);
  boundaryProc_DoIt(step);
  WLSProc_DoIt(step);

  //G4String particleName  = "e-";
  G4String particleName  = "alpha";
  //G4String particleName  = "proton";


  opPhoton_DoIt(step);

  gamma_DoIt(step);

  neutron_DoIt(step);


#if VERBOSE
  stepVerboseInfo(step);
#endif //VERBOSE


  return;
}


void Cube666_SteppingAction::boundaryProc_DoIt(const G4Step* step){

  /*
    //checking boundary process

  static G4OpBoundaryProcess* boundaryProc=NULL;
  if(!boundaryProc){
    G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nproc = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    for(G4int i=0;i<nproc;i++){
      if((*pv)[i]->GetProcessName() == "OpBoundary"){
	boundaryProc = (G4OpBoundaryProcess*)(*pv)[i];
	break;
      }
    }
  }


  G4OpBoundaryProcessStatus boundaryProcStatus = Undefined;
  if(boundaryProc) boundaryProcStatus = boundaryProc->GetStatus();
  
  if(boundaryProcStatus != 4) return;
  */

  return;
}



void Cube666_SteppingAction::WLSProc_DoIt(const G4Step* step){

  //old code deleted
  //do something here
 
  return;
}



void Cube666_SteppingAction::opPhoton_DoIt(const G4Step* step){
 
  
  
  D("\t "<<__FUNCTION__);

  if(!isParticle("opticalphoton",step)) return;


  G4StepPoint* postStepPoint = step->GetPostStepPoint();

#if TEST_BRANCH0
#if TEST_BRANCH00

  //recording direct (UV-)light (ene > 6eV)


  if(isVolume("pmtCoat",step) || isVolume("pmtCathode",step) ){

    G4String volName = step->GetTrack()->GetVolume()->GetName();
    
    //if(isNextVolume("pmtCoat",step) || isNextVolume("pmtCathode",step) ){
    //G4String volName = step->GetTrack()->GetNextVolume()->GetName();
    
    G4String pmtIDstr = (volName.contains("pmtCoat")) ? volName(7,2) : volName(10,2);
    int pmtID = atoi(pmtIDstr.data());
    
    Cube666_Test* testobj = Cube666_Test::getInstance();
    
    if(step->GetTrack()->GetTotalEnergy() > 6*eV ){
      testobj->fNDirectPhoton++;
      testobj->fTotNDirectPhotonPMT[pmtID]++;
    }
    
    
    //this is to avoid double counting fTotNPhotonHittingPMT
    //only count the photon if it stops / is killed on the PMTs
    if(step->GetTrack()->GetTrackStatus() == fStopAndKill ||
       step->GetTrack()->GetTrackStatus() == fKillTrackAndSecondaries){
      
      testobj->fNPhotonHittingPMT++;
      testobj->fTotNPhotonHittingPMT[pmtID]++;
    }
    
    
  }
  
  
  
#endif //TEST_BRANCH00
#endif //TEST_BRANCH0
  



  return;
}






void Cube666_SteppingAction::stepVerboseInfo(const G4Step* step){
    cout//<<"in Cube666_SteppingAction::UserSteppingAction(..)  .... "
      //<<"SteppingAction(..): "
      //<<"StepInfo : "
      //<<" "<<step->GetTrack()->GetTrackStatus()
      //<<"\t "
      <<step->GetTrack()->GetParticleDefinition()->GetParticleName()
      //<<"\t particleID "<<step->GetTrack()->GetParticleDefinition()->GetPDGEncoding()
      //<<"\t charge "<<step->GetTrack()->GetParticleDefinition()->GetPDGCharge()
      //<<"\t parentID "<<step->GetTrack()->GetParentID()
      //<<"\t primPart "<<step->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetPDGcode()
      //<<"\t primPart* "<<step->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()
      //<<"\t daughter "<<step->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetDaughter()->GetPDGcode()
      //<<"\t daughter* "<<step->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetDaughter()
      //<<"\t nextPart* "<<step->GetTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetNext()
      <<"\t "<<step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
      <<"\t steplength [cm] "<<step->GetStepLength()/cm
      //<<"\t Etot [keV] "<<step->GetTrack()->GetDynamicParticle()->GetTotalEnergy()/keV
      //<<"\t Ekin [keV] "<<step->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/keV
      //<<"\t preStepEkin [keV] "<<step->GetPreStepPoint()->GetKineticEnergy()/keV
      //<<"\t postStepEkin [keV] "<<step->GetPostStepPoint()->GetKineticEnergy()/keV
      //<<"\t DelE [keV] "<<(step->GetPreStepPoint()->GetKineticEnergy()-step->GetPostStepPoint()->GetKineticEnergy())/keV
      //<<"\t Edep [keV] "<<step->GetTotalEnergyDeposit()/keV
      //<<"\t DelE - Edep "<<(step->GetPreStepPoint()->GetKineticEnergy()-step->GetPostStepPoint()->GetKineticEnergy())/keV - step->GetTotalEnergyDeposit()/keV
      //<<"\t EdepNonIoni [keV] "<<step->GetNonIonizingEnergyDeposit()/keV
      //<<"\t EdepRad [keV] "<<-(step->GetDeltaEnergy()/keV) - (step->GetTotalEnergyDeposit()/keV)
      //<<"\t EdepIoni [keV] "<<step->GetTotalEnergyDeposit()/keV - step->GetNonIonizingEnergyDeposit()/keV
      //<<"\t test2 "<<step->GetTotalEnergyDeposit() + step->GetDeltaEnergy()
      <<"\t "<<step->GetTrack()->GetVolume()->GetName();

    if(step->GetTrack()->GetNextVolume())
      cout<<"\t "<<step->GetTrack()->GetNextVolume()->GetName();

    //cout<<"\t "<<step->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();

    //G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    //cout<<"\t "<<pos.x()<<" "<<pos.y()<<" "<<pos.z();

    //cout<<"\t trackStatus "<<step->GetTrack()->GetTrackStatus();


    cout<<endl;

  return;
}





void Cube666_SteppingAction::scattProc_DoIt(const G4Step* step){

  D("\t "<<__FUNCTION__);

  if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() !="OpRayleigh") return;

//   G4StepPoint* preStepPoint  = step->GetPreStepPoint();
//   G4StepPoint* postStepPoint = step->GetPostStepPoint();
//   G4ThreeVector prevMomDir   = preStepPoint->GetMomentumDirection();
//   G4ThreeVector currentMomDir = postStepPoint->GetMomentumDirection();
  
//   G4double scattAngle = prevMomDir.angle(currentMomDir);
//   //if(scattAngle>1e-4){
//   if(scattAngle>0) Cube666_Test::getInstance()->fNscatt++;
  
//   Cube666_Analysis* ana = Cube666_Analysis::getInstance();
//   ana->fScattAngle.push_back(scattAngle);
  
//   G4ThreeVector prevPol    = preStepPoint->GetPolarization();
//   G4ThreeVector currentPol = postStepPoint->GetPolarization();
//   ana->fPolAngle.push_back(prevPol.angle(currentPol));
//   ana->fNewPol_oldMom_Angle.push_back(currentPol.angle(prevMomDir));
//   ana->fNewMom_oldMom_Angle.push_back(currentMomDir.angle(prevMomDir));
//   ana->fNewMom_oldPol_Angle.push_back(currentMomDir.angle(prevPol));
//   ana->fNewPol_oldPol_Angle.push_back(currentPol.angle(prevPol));
  
  return;
}





void Cube666_SteppingAction::gamma_DoIt(const G4Step* step){

  D("\t "<<__FUNCTION__);

  if(!isParticle("gamma",step)) return;

  //old code deleted
  //do something with gamma

  return;
} 









void Cube666_SteppingAction::stepVerboseInfo_secondaries(const G4Step* step){

  const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
  if(!secondary) return;

  G4int nsec = (*secondary).size();
  if(!nsec) return;

  cout<<step->GetTrack()->GetDefinition()->GetParticleName()
      <<"\t "<<step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
      <<"\t preStepEkin [keV] "<<step->GetPreStepPoint()->GetKineticEnergy()/keV
      <<"\t postStepEkin [keV] "<<step->GetPostStepPoint()->GetKineticEnergy()/keV
      <<"\t Edep [keV] "<<step->GetTotalEnergyDeposit()/keV
      <<"\t "<<nsec<<" secondary particle(s) "
      <<endl;

  for(G4int i=0;i<nsec;i++) printTrackInfo((*secondary)[i],step);  
  //for(G4int i=0;i<nsec;i++){ printTrackInfo((*secondary)[i],step);  if(!(i%10))getchar();}

  cout<<endl<<endl;

  return;
}



void Cube666_SteppingAction::printTrackInfo(G4Track* track,const G4Step* step){
  
  cout<<track->GetParticleDefinition()->GetParticleName()
      <<"\t Ekin [keV] "<<track->GetKineticEnergy()/keV

      <<endl;

}




void Cube666_SteppingAction::printParentInfo(const G4Step* step){

  G4int parentID = step->GetTrack()->GetParentID();
  
  G4TrajectoryContainer* container = G4RunManager::GetRunManager()->GetCurrentEvent()->GetTrajectoryContainer();
  if(!container) return;
  
  
  G4int ntrj = container->size();
  for(G4int i=0;i<ntrj;i++){
    Cube666_Trajectory* trj = (Cube666_Trajectory*)((*container)[i]);
    if(trj && trj->GetTrackID() == parentID){
      cout<<" currentPart "<<step->GetTrack()->GetParticleDefinition()->GetParticleName()
	//<<"\t Ekin [keV] "<<step->GetTrack()->GetKineticEnergy()/keV
	  <<"\t preStepEkin [keV] "<<step->GetPreStepPoint()->GetKineticEnergy()/keV
	  <<"\t postStepEkin [keV] "<<step->GetPostStepPoint()->GetKineticEnergy()/keV
	  <<"\t parentName "<<trj->GetParticleName()
	  <<"\t parentEkin [keV] "<<trj->GetEkin()/keV
	  <<endl;
      return;
    }
  }

  
  return;

}






void Cube666_SteppingAction::neutron_DoIt(const G4Step* step){

  if(!isParticle("neutron",step) ) return;

  //old code removed
  //do something with neutron


  return;
}




G4int Cube666_SteppingAction::isParticle(const char* particleName,const G4Step* step){

  if(step->GetTrack()->GetParticleDefinition()->GetParticleName().contains(particleName)) return 1;
  return 0;
}



G4int Cube666_SteppingAction::isProcess(const char* processName,const G4Step* step){

  if(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().contains(processName)) return 1;
  return 0;
}




G4int Cube666_SteppingAction::isVolume(const char* volumeName,const G4Step* step){

  if(!step->GetTrack()->GetVolume()) return 0;
  if(!step->GetTrack()->GetVolume()->GetName().contains(volumeName)) return 0;
  return 1;
}




G4int Cube666_SteppingAction::isNextVolume(const char* nextVolumeName,const G4Step* step){

  if(!step->GetTrack()->GetNextVolume()) return 0;
  if(!step->GetTrack()->GetNextVolume()->GetName().contains(nextVolumeName)) return 0;
  return 1;
}




