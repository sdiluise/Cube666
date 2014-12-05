#include "Cube666_Analysis.hh"
#include "preparation.hh"
#include "Cube666_Test.hh"

#include "algorithm"

#include "TH1F.h"
#include "TF1.h"



Cube666_Analysis* Cube666_Analysis::fInstance = 0;

G4String Cube666_Analysis::fFilename="";

Cube666_Analysis::Cube666_Analysis(){


#if TEST_BRANCH0
  fNtimeSample = NTIMESAMPLE;
  fTimeSample  = TIMESAMPLE;
  fNPMT        = NPMT;

#endif //TEST_BRANCH0
  
}


Cube666_Analysis::~Cube666_Analysis(){
  if(fOutputFile) delete fOutputFile;
  

}



Cube666_Analysis* Cube666_Analysis::getInstance(){
  if(!fInstance) fInstance = new Cube666_Analysis;
  return fInstance;
}


void Cube666_Analysis::BeginOfRun(){
  fAnaTree = NULL;
  fOutputFile = NULL;
  bookTree();

  return; 
}


void Cube666_Analysis::EndOfRun(){

  fAnaTree->Write();
  fOutputFile->Close();

  if(VERBOSE) G4cout<<"in Cube666_Analysis::EndOfRun(), tree written."<<G4endl;
  return;
}


void Cube666_Analysis::BeginOfEvent(){

  Reset();
}


void Cube666_Analysis::Reset(){

#if TEST_BRANCH0
  
  fNphotonDetected = 0;
  fNphotonEmitted=0;
  

#if TEST_BRANCH00
//   //for safety reason, not really necessary !
//   reset1DArray<int>(fTotNphotonPMT,NPMT,0);
//   reset2DArray<int>(fNphotonPMT,NPMT,NTIMESAMPLE,0);

  fNDirectPhoton = 0;
  fNPhotonHittingPMT = 0;

  reset1DArray<G4int>(fTotNphotonPMT,NPMT,0);
  for(int i=0;i < NPMT;i++) reset1DArray<G4int> (fNphotonPMT[i],NTIMESAMPLE,0);
  for(int i=0;i < NPMT;i++) reset1DArray<float> (fRawData[i],NTIMESAMPLE,0);

  reset1DArray<G4int>(fTotNDirectPhotonPMT,NPMT,0);
  reset1DArray<G4int>(fTotNPhotonHittingPMT,NPMT,0);
 
  fPosx = DEFAULTVALUE; 
  fPosy = DEFAULTVALUE;
  fPosz = DEFAULTVALUE;

#endif

#endif //TEST_BRANCH0


  Cube666_Test::getInstance()->reset();
  return;

}


void Cube666_Analysis::EndOfEvent(){
  Cube666_Test* testObj = Cube666_Test::getInstance();



#if TEST_BRANCH0
#if TEST_BRANCH00

  copy1DArray(testObj->fTotNDirectPhotonPMT,fTotNDirectPhotonPMT,NPMT);
  fNDirectPhoton = testObj->fNDirectPhoton;

  copy1DArray(testObj->fTotNPhotonHittingPMT,fTotNPhotonHittingPMT,NPMT);
  fNPhotonHittingPMT = testObj->fNPhotonHittingPMT;


#endif //TEST_BRANCH00
#endif //TEST_BRANCH0


#if TEST_BRANCH0

  fAnaTree->Fill();

#endif //TEST_BRANCH0



#if VERBOSE
  G4cout<<"in Cube666_Analysis::EndOfEvent(), tree filled."<<G4endl;

#endif // VERBOSE


  return;
}



void Cube666_Analysis::bookTree(){

  cout<<"booking tree .... "<<endl;

  if(!fOutputFile) fOutputFile = NULL;

  if(fFilename==""){
    ostringstream str;
    str<<ARDM_OUTPUT_TREE_DIR_DEFAULT<<"Cube666_Analysis_tree.root";
    fFilename=str.str();
  }

  cout<<"creating tree  "<<fFilename<<endl;

  fOutputFile = new TFile(fFilename.c_str(),"RECREATE");
  //fAnaTree    = new TTree("tree",outputFilename.c_str());
  fAnaTree    = new TTree("Data",fFilename.c_str());
  
  
#if TEST_BRANCH0
  fAnaTree->Branch("fNPMT",            &fNPMT);
  fAnaTree->Branch("fNtimeSample",     &fNtimeSample);
  fAnaTree->Branch("fNphotonDetected", &fNphotonDetected);
  fAnaTree->Branch("fTimeSample",      &fTimeSample);
  fAnaTree->Branch("fNphotonEmitted",        &fNphotonEmitted);


#if TEST_BRANCH00  

  fAnaTree->Branch("fTotNphotonPMT", &fTotNphotonPMT, Form("fTotNphotonPMT[%i]/I",fNPMT));

  //recording direct light
  fAnaTree->Branch("fNDirectPhoton",       &fNDirectPhoton);
  fAnaTree->Branch("fTotNDirectPhotonPMT", &fTotNDirectPhotonPMT, Form("fTotNDirectPhotonPMT[%i]/I",fNPMT));
  fAnaTree->Branch("fNPhotonHittingPMT",   &fNPhotonHittingPMT);
  fAnaTree->Branch("fTotNPhotonHittingPMT",&fTotNPhotonHittingPMT,Form("fTotNPhotonHittingPMT[%i]/I",fNPMT));

  ostringstream branchname,leafname;

  for(int i=0;i<fNPMT;i++){
    branchname.str("");
    leafname.str("");
    branchname<<"fNphotonPMT"<<i;
    leafname<<"fNphotonPMT["<< fNtimeSample <<"]/I";
    fAnaTree->Branch(branchname.str().c_str(),&fNphotonPMT[i],leafname.str().c_str());
  }

  
  
  for(int i=0;i<fNPMT;i++){
    branchname.str("");
    leafname.str("");
    branchname<<"eRawData"<<i;
    leafname<<"eRawData["<< fNtimeSample <<"]/F";
    fAnaTree->Branch(branchname.str().c_str(),&fRawData[i],leafname.str().c_str());
  }


  fAnaTree->Branch("fPosx",    &fPosx);
  fAnaTree->Branch("fPosy",    &fPosy);
  fAnaTree->Branch("fPosz",    &fPosz);


#endif //TEST_BRANCH00
  
#endif //TEST_BRANCH0
  

  return;
}



