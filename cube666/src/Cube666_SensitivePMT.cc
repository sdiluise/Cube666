#include "Cube666_SensitivePMT.hh"
#include "preparation.hh"
#include "Cube666_Analysis.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"


const G4double Cube666_SensitivePMT :: fTimeSample   = TIMESAMPLE;
const G4int    Cube666_SensitivePMT :: fNTimeSamples = NTIMESAMPLE;
Cube666_Digitizer* Cube666_SensitivePMT:: digitizer     = NULL;

Cube666_SensitivePMT::Cube666_SensitivePMT(G4int npmt,G4String detName,G4String hitsColName)
  :fnpmt(npmt),
   fDetName(detName),
   fHitsColName(hitsColName),
   G4VSensitiveDetector(detName){
  
  fNormal.set(0.,0.,1.);
  
  collectionName.insert(hitsColName);
  fPMTvector=setPMTVector();

  //fPMTvector is not used to place the PMTarray in any physical volume,
  //it is rather used to calculate the angle, at which the photon arrives on the PMT surface,
  //so it should points to the center of the sensitive spherical part of the PMT,
  //and not to the center of the pmt as described in the approx_pmt_... part in preparation.hh,
  //or in preration.cc :: constructApproxPMT(...).
  //
  //the xy coordinate are the same for both center of sensitive part and center of pmt itself,
  //but z-coord. is different.
  //
  //--> so now, reset the z-coord.
  //

  digitizer = new Cube666_Digitizer("PMTDigitizer");
  
}  

  
  
Cube666_SensitivePMT::~Cube666_SensitivePMT(){
  fPMTvector.clear();
}



G4double Cube666_SensitivePMT::detectionProb(G4double phi){
  
  phi *= 180/PI; //convert from radian to degree.
  //cout<<"in Cube666_SensitivePMT::detectionProb(G4double phi) ... phi = "<<phi<<endl;

  if(phi > 46.5) return 0;
  return (2587 + 1134/(1 + (phi/27)*(phi/27)))/3722;  
}


G4double Cube666_SensitivePMT::detectionProb(G4ThreeVector hitpos, G4ThreeVector rPMT){
  return detectionProb(fNormal.angle(hitpos-rPMT));
}


G4double Cube666_SensitivePMT::detectionProb(G4ThreeVector hitpos,G4int pmtID){
  return detectionProb(hitpos,fPMTvector[pmtID]);
}


G4int Cube666_SensitivePMT::detected(G4ThreeVector hitpos,G4ThreeVector rPMT){
  
  if(drand48() < detectionProb(hitpos,rPMT)) return 1;
  else return 0;
}


G4int Cube666_SensitivePMT::detected(G4ThreeVector hitpos,G4int pmtID){

#if VERBOSE
  cout<<"pmtID "<<pmtID<<"\t fPMTvector.size() "<<fPMTvector.size()<<endl;
#endif //VERBOSE

  return detected(hitpos,fPMTvector[pmtID]);
}





void Cube666_SensitivePMT::Initialize(G4HCofThisEvent* hc){
  fHitsCol = new Cube666_PMTHitsCollection(fDetName,fHitsColName);//"PMTarray","PMTHits");
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsColName);
  hc->AddHitsCollection(hcID,fHitsCol);
  verboseLevel = 0;

#if VERBOSE
  cout<<"hitsColName = "<<fHitsColName <<"\t id "<<hcID<<"\t nHitsColInThisEvt "<<G4SDManager::GetSDMpointer()->GetHCtable()->entries() <<endl;
  //getchar();
#endif //VERBOSE

  return;
}




G4bool Cube666_SensitivePMT::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist){

  if(!aStep->GetTrack()->GetParticleDefinition()->GetParticleName().contains("opticalphoton")) return true;

  if(aStep->GetTrack()->GetDynamicParticle()->GetTotalEnergy() > PHOTON_ENERGY_THRESHOLD){
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);    
    return true;
  }

  G4int photonDetected = 0;


#if TEST_BRANCH0
    
  Cube666_PMTHit* newhit = new Cube666_PMTHit;
  G4Track* track = aStep->GetTrack();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();

  //G4int globalTime = (int)(track->GetGlobalTime()/ns/fTimeSample);
  G4int globalTime = track->GetGlobalTime()/ns;
  G4String volName = track->GetVolume()->GetName();

  //10 = length of string "pmtCathode", take the last 2 character
  //e.g. : volName = "pmtCathode3"  --> pmtIDstr = volName.substr(10,2) = "3"  --> pmtID = 3
  //       volName = "pmtCathode"   --> pmtIDstr = volName.substr(10,2) = "17" --> pmtID = 17


  G4String pmtIDstr = volName(10,2); //<--> volName.substr(10,2)

  G4int pmtID = atoi(pmtIDstr.data());

  //check if the photon is detected
  //store only detected photons !
  photonDetected = detected(hitpos,pmtID);

  newhit->setHitGlobalTime(globalTime);
  newhit->setHitPMTid(pmtID);
  newhit->setIsDetected(photonDetected);

  fHitsCol->insert(newhit);


#if VERBOSE
  cout<<"in Cube666_SensitivePMT::ProcessHits(..) ... volName "<<volName<<"\t pmtID "<<pmtID<<"\t isDetected "<<photonDetected <<endl;
#endif //VERBOSE


#endif //TEST_BRANCH0

  //photon falls onto PMT's surface --> kill it !
  aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

  return true;
}



void Cube666_SensitivePMT::EndOfEvent(G4HCofThisEvent* hc){

  digitizer->Digitize_simple(fHitsColName);

  return;
}




