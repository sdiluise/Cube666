
#include "preparation.hh"
#include "Cube666_Digitizer.hh"

#include "G4DigiManager.hh"
#include "Cube666_PMTHit.hh"
#include "Cube666_Analysis.hh"

#include "TRandom.h"

#include "assert.h"


#define ARDM_PMT_RISETIME (12*ns)
#define ARDM_TRANSIT_TIME (68*ns) 
#define ARDM_SIGMA_TRANSIT_TIME (2.8*ns) 
#define ARDM_MEAN_PEAK_WIDTH (28.72*ns)
#define ARDM_MEAN_FWHM (23.2 *ns)
#define ARDM_PMT_TAU (7.51*ns)




int Cube666_Digitizer::fNTimeSamples = NTIMESAMPLE;
double* Cube666_Digitizer::fPMTCalib = NULL;
double* Cube666_Digitizer::fPMTCalibSigma = NULL;
double* Cube666_Digitizer::fPMTDarkCountRate = NULL;
double* Cube666_Digitizer::fPMTWhiteNoiseSigma = NULL;


Cube666_Digitizer::Cube666_Digitizer(G4String name,G4String pmtCalibFilename)
  :G4VDigitizerModule(name){


  fPMTCalib = new double[NPMT];
  fPMTCalibSigma = new double[NPMT];
  fPMTDarkCountRate = new double[NPMT];
  fPMTWhiteNoiseSigma = new double[NPMT];

  if(pmtCalibFilename != ""){
    
    //calibration file for the PMTs should have the following format :
    //pmtid,pmtCalibConstant,error_on_pmtCalibConstant,darkCountRate,whiteNoiseSigma
    //unwanted lines can be commented out with #
    //delimiter can be ',' , ' ', '\t' 
    //the 2nd argument in readTextfile_float is the delimiter

    vector<vector<double> > calib;// = readTextfile_float(pmtCalibFilename,',');
    for(int i=0; i < NPMT;i++){
      if(i < calib.size()){
	fPMTCalib[i] = calib[i][1];
	fPMTCalibSigma[i] = calib[i][2];
	fPMTDarkCountRate[i] = calib[i][3];
	fPMTWhiteNoiseSigma[i] = calib[i][4];

      }else{
	fPMTCalib[i] = 1;
	fPMTCalibSigma[i] = 0;
	fPMTDarkCountRate[i] = 0;
	fPMTWhiteNoiseSigma[i] = 0;      
      }
    
    }
  
  }else{
  
    for(int i=0; i < NPMT;i++){
      fPMTCalib[i] = 1;
      fPMTCalibSigma[i] = 0;
      fPMTDarkCountRate[i] = 0;
      fPMTWhiteNoiseSigma[i] = 0;
    }
  
  }

}



Cube666_Digitizer::~Cube666_Digitizer(){;}






void Cube666_Digitizer::reset(float* array, int nentries){

  for(int i=0;i<nentries;i++) array[i] = 0;
  return;
}





void Cube666_Digitizer::Digitize(G4HCofThisEvent* hc){
  
  Digitize_simple(hc);
  //Digitize_pmtRes(hc);
  
  return; 
}





void Cube666_Digitizer::Digitize(){;}




void Cube666_Digitizer::Digitize_simple(Cube666_PMTHit* hit){

  if(!hit->getIsDetected()) return;


  int pmtid = hit->getHitPMTid();
  
  Cube666_Analysis* ana = Cube666_Analysis::getInstance();

#if TEST_BRANCH0
  
  ana->fNphotonDetected++;
  
  ana->fTotNphotonPMT[pmtid]++;

#if TEST_BRANCH00
  
  int globalTime = (int)(hit->getHitGlobalTime()/TIMESAMPLE);
  if(globalTime < fNTimeSamples) ana->fNphotonPMT[pmtid][globalTime]++;
    
#endif //TEST_BRANCH00
    
  

#endif //TEST_BRANCH0


  return;
}









void Cube666_Digitizer::Digitize_simple(Cube666_PMTHitsCollection* hitsCol){

  int nhits = hitsCol->entries();
  for(int i=0; i < nhits; i ++ )   Digitize_simple((*hitsCol)[i]);
  
  return;
}







void Cube666_Digitizer::Digitize_simple(G4String hitsColName){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();

#if VERBOSE
  cout<<"hitsColName "<<hitsColName<<endl; //getchar();
#endif //VERBOSE

  G4int hitColID = digiMan->GetHitsCollectionID(hitsColName);

  Cube666_PMTHitsCollection* hitCol = (Cube666_PMTHitsCollection*) (digiMan->GetHitsCollection(hitColID));
  if(hitCol) Digitize_simple(hitCol);

  return; 
}





void Cube666_Digitizer::Digitize_simple(G4HCofThisEvent* hc){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
  
  for(int hci=0; hc->GetNumberOfCollections();hci++){
    Cube666_PMTHitsCollection* hitCol = (Cube666_PMTHitsCollection*) (digiMan->GetHitsCollection(hci));
    if(hitCol) Digitize_simple(hitCol);
  }

  return; 
}








void Cube666_Digitizer::Digitize_pmtRes_try1(Cube666_PMTHit* hit,float* rawDataArray){

  //simulting PMT (voltage) response to one single photon hit
  //right flank : exponentially decaying funtion,
  //left  flank : gaussian


  if(!hit->getIsDetected()) return;

  double peakStartTime = hit->getHitGlobalTime();
  double peakTime      = peakStartTime + ARDM_PMT_RISETIME;
  double peakStopTime  = peakStartTime + ARDM_MEAN_PEAK_WIDTH;

  int pmtid = hit->getHitPMTid();
  double pmtRes = gRandom->Gaus(fPMTCalib[pmtid],fPMTCalibSigma[pmtid]);
  while(pmtRes <= 0  ) pmtRes = gRandom->Gaus(fPMTCalib[pmtid],fPMTCalibSigma[pmtid]);


  double tau = ARDM_PMT_TAU;
  double Vo = pmtRes * TIMESAMPLE / (tau + sqrt(PI/2) * ARDM_SIGMA_TRANSIT_TIME ) ;

  int peakStartTimeInt = (int) (peakStartTime / TIMESAMPLE);
  int peakTimeInt      = (int) (peakTime      / TIMESAMPLE);
  int peakStopTimeInt  = (int) (peakStopTime  / TIMESAMPLE);

#if VERBOSE
  cout<<"peakStartTimeInt "<<peakStartTimeInt
      <<"\t peakTimeInt "<<peakTimeInt
      <<"\t peakStopTimeInt "<<peakStopTimeInt
      <<"\t Vo "<<Vo
      <<"\t pmtRes "<<pmtRes
      <<endl;
  
#endif //VERBOSE

  for(int i=peakStartTimeInt; i < peakTimeInt; i++) rawDataArray[i] += Vo * TMath::Gaus((i+1./2)*TIMESAMPLE, peakTime,ARDM_SIGMA_TRANSIT_TIME,1);
  for(int i=peakTimeInt; i <= peakStopTimeInt; i++) rawDataArray[i] += Vo * exp( - ((i+1./2)*TIMESAMPLE - peakTime)/tau);


  //for(int i=0;i<fNTimeSamples;i++){ if(rawDataArrayy[i]) cout<< i <<"\t " << rawDataArray[i]<<endl;  }

  return ;

}






void Cube666_Digitizer::Digitize_pmtRes_try2(Cube666_PMTHit* hit,float* rawDataArray){

  //simulting PMT (voltage) response to one single photon hit
  //pmtResponse = log-normal distribution
  //V = Vo * exp( - log(t/tau) * log(t/tau) / 2 / sigma /sigma  )

  //or more precisely :
  //
  //V = Vo * exp( - log( (t-peakStartTime)/tau) * log( (t-peakStartTime)/tau) / 2 / sigma / sigma )
  //
  //can try : tau ~ FWHM ~ 4*5.8 ns , sigma ~ .25 (from real data fitting)
  //
  //Vo varies from pmt to pmt
  //tau, sigma can be assumed to be the same for all PMTs


  if(!hit->getIsDetected()) return;

  double peakStartTime = hit->getHitGlobalTime();
  double peakStopTime  = peakStartTime + ARDM_MEAN_PEAK_WIDTH;

  int pmtid = hit->getHitPMTid();
  double pmtRes = gRandom->Gaus(fPMTCalib[pmtid],fPMTCalibSigma[pmtid]);
  while(pmtRes <= 0  ) pmtRes = gRandom->Gaus(fPMTCalib[pmtid],fPMTCalibSigma[pmtid]);

  //double Vo = pmtRes * TIMESAMPLE / (tau + sqrt(PI/2) * ARDM_SIGMA_TRANSIT_TIME ) ;
  double tau = ARDM_MEAN_FWHM;
  double sigma = .25;

  double Vo = pmtRes * TIMESAMPLE / (tau * sigma * sqrt(2*PI) * exp(sigma*sigma/2));

  int peakStartTimeInt = (int) (peakStartTime / TIMESAMPLE);
  int peakStopTimeInt  = (int) (peakStopTime  / TIMESAMPLE);

#if VERBOSE
  cout<<"peakStartTimeInt "<<peakStartTimeInt
      <<"\t peakStopTimeInt "<<peakStopTimeInt
      <<"\t Vo "<<Vo
      <<"\t pmtRes "<<pmtRes
      <<endl;
  
#endif //VERBOSE

  for(int i=peakStartTimeInt; i <= peakStopTimeInt; i++) 
    rawDataArray[i] += Vo * TMath::Gaus(log( ( (i+1./2)*TIMESAMPLE - peakStartTime) / tau),0,sigma);


  return ;
}





void Cube666_Digitizer::Digitize_pmtRes(Cube666_PMTHit* hit,float* rawDataArray){

  //Digitize_pmtRes_try1(hit,rawDataArray);
  Digitize_pmtRes_try2(hit,rawDataArray);
  return;
}






void Cube666_Digitizer::Digitize_pmtRes(Cube666_PMTHit* hit){

  //the electron cloud from dynodes needs some transit time ts to reach the anode,
  //the transit time ts follows a gaussian distribution with mean mean_ts and width sigma_ts 
  //given by the technical sheet of the PMT.
  //
  //--> how to simulate the voltage signal of a single photon ? 
  //
  //i.   photon hits PMT surface at time t0
  //ii.  electron cloud hits the anode at time t1 = t0 + gauss(mean=transitTime,sigma=sigma_ts)
  //       --> peakTime = t1
  //iii. peakStartTime = peakTime - riseTime
  //iv.  stopTime = t2 = peakTime + mean(stopTime - peakTime) <-- mean(stopTime - peakTime) from real data
  //v.   get pmtRes = gauss(mean_pmtRes,sigma_pmtRes) <-- mean, sigma from real data for each PMT
  //vi.  now we have to throw random numbers to generate the amplitude of the voltage signal for each time sample,
  //     so that the peakIntegral we get from that is equal to the pmtRes in (v.)
  //
  //since the transit time of the electron cloud is small compared to the processing time of the amplifier after that,
  //in order to simplify things, instead of proceed like described above, we can do the following :
  //
  //i.   photon hits surface at time t0 =: startTime
  //ii.  t1 = t0 + rise time  =: peakTime
  //iii. t2 = t1 + mean(stopTime - peakTime) =: stopTime
  //iv.  pmtRes = gauss(mean_pmtRes,sigma_pmtRes)
  //v.   now we have to simulate the amplitude of the voltage signal :
  //      1. left  flank (startTime < t < peakTime) = (half-)gaussian with mean = peakTime, sigma = transit time spread
  //      2. right flank (peakTime  < t < stopTime) = exponential decay with tau = 1/RC <-- R,C from the circuit of the PMT, emperical value.
  //      
  //vi.  left blank  = Vo * exp( - (t-peakTime) * (t-peakTime) / 2/sigma_ts/sigma_s)
  //     right blank = Vo * exp( - (t-peakTime)/tau)
  //
  //vii. integral(left blank) = int(left blank)_peakStartTime ^peakTime ~= int(left blank)_minusInf ^peakTime = Vo * sqrt(2pi) * sigma_ts
  //     integral(left blank) = int(left blank)_peakTime ^peakStopTime ~= int(left blank)_peakTime ^inf = Vo * tau
  // 
  //     --> integral = pmtRes = Vo * (tau + sqrt(2pi) * sigma_ts)
  //     --> Vo = pmtRes / (tau + sqrt(2pi) * sigma_ts)
  //
  //



#if TEST_BRANCH00

  int pmtid = hit->getHitPMTid();
  Digitize_pmtRes(hit,Cube666_Analysis::getInstance()->fRawData[pmtid]);
  
#endif //TEST_BRANCH00

  return;


  // if(!hit.getIsDetected()) return;

  // Cube666_Analysis* ana = Cube666_Analysis::getInstance();

  // double peakStartTime = hit.getHitGlobalTime();
  // double peakTime      = peakStartTime + ARDM_PMT_RISETIME;
  // double peakStopTime  = peakStartTime + ARDM_MEAN_PEAK_WIDTH;

  // int pmtid = hit.getHitPMTid();
  // double pmtRes = gRandom->Gaus(fPMTCalib[pmtid],fPMTCalibSigma[pmtid]);

  // double tau = ARDM_PMT_TAU;
  // double Vo = pmtRes / (tau + sqrt(2*PI) * ARDM_SIGMA_TRANSIT_TIME) ;

  // int peakStartTimeInt = (int) (peakStartTime / TIMESAMPLE);
  // int peakTimeInt      = (int) (peakTime      / TIMESAMPLE);
  // int peakStopTimeInt  = (int) (peakStopTime  / TIMESAMPLE);

  // for(int i=peakStartTimeInt; i < peakTimeInt; i++) ana->fRawData[pmtid][i] += Vo * TMath::Gaus((i+1./2)*TIMESAMPLE, peakTime,ARDM_SIGMA_TRANSIT_TIME,1);
  // for(int i=peakTimeInt; i <= peakStopTimeInt; i++) ana->fRawData[pmtid][i] += Vo * exp( - ((i+1./2)*TIMESAMPLE - peakTime)/tau);

  // return ;
}









void Cube666_Digitizer::Digitize_pmtRes(Cube666_PMTHitsCollection* hitcol){

  int nhits = hitcol->entries();
  for(int i=0;i<nhits;i++) Digitize_pmtRes((*hitcol)[i]);

  return;
}




void Cube666_Digitizer::Digitize_pmtRes(G4String hitsColName){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();

  G4int hitColID = digiMan->GetHitsCollectionID(hitsColName);
  Cube666_PMTHitsCollection* hitCol = (Cube666_PMTHitsCollection*) (digiMan->GetHitsCollection(hitColID));
  Digitize_pmtRes(hitCol);
  return;
}





void Cube666_Digitizer::Digitize_pmtRes(G4HCofThisEvent* hc){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
  
  for(int hci=0; hci < hc->GetNumberOfCollections(); hci++){
    Cube666_PMTHitsCollection* hitCol = (Cube666_PMTHitsCollection*) (digiMan->GetHitsCollection(hci));
    Digitize_pmtRes(hitCol);
  }

  return;
}









void Cube666_Digitizer::Digitize_darkCurrent(int pmtid,float* rawDataArray){

  assert(pmtid >=0 && pmtid < NPMT);

  int nelectrons = gRandom->Poisson(fPMTDarkCountRate[pmtid]);

  //cout<<"nelectrons "<<nelectrons<<endl;

  if(nelectrons <= 0) return;

  double timespan = TIMESAMPLE * fNTimeSamples;
  for(int i=0; i < nelectrons;i++){
    double time = gRandom->Uniform(0,timespan);
    
    Cube666_PMTHit* hit;
    hit->setHitGlobalTime(time);
    hit->setHitPMTid(pmtid);
    hit->setIsDetected(1);

    Digitize_pmtRes(hit,rawDataArray);  
  }

  return;
}








void Cube666_Digitizer::Digitize_darkCurrent(int pmtid){

  //probability that an electron gets out of the PMT cathode by itself is given by dark count rate
  //the array PMT_DARK_COUNT_RATE above gives the mean number of "dark electrons"
  //in a time span, which corresponds to the timespan we would record for a normal event ( = 2048 samples * 4ns/sample)
  //--> from real data
  //
  //--> throw a poisson distribution around that mean value to get the number of dark electrons to be emitted.
  //
  //pulse shape for "dark electron" is the same for the one induced by an optical photon.
  //the probability that a dark electron appears is uniform over the whole time range of an optical-photon-event.


  assert(pmtid >=0 && pmtid < NPMT);

  int nelectrons = gRandom->Poisson(fPMTDarkCountRate[pmtid]);

  if(nelectrons <= 0) return;

  double timespan = TIMESAMPLE * fNTimeSamples;
  for(int i=0; i < nelectrons;i++){
    double time = gRandom->Uniform(0,timespan);
    
    Cube666_PMTHit* hit;
    hit->setHitGlobalTime(time);
    hit->setHitPMTid(pmtid);
    hit->setIsDetected(1);

    Digitize_pmtRes(hit);  
  }

  return;
}









void Cube666_Digitizer::Digitize_whiteNoise(int pmtid,float* rawDataArray){

#if TEST_BRANCH00
  assert(pmtid >=0 && pmtid < NPMT);

  for(int i=0;i < fNTimeSamples;i++){
    rawDataArray[i] += gRandom->Gaus(0,fPMTWhiteNoiseSigma[pmtid]);  
  }
#endif //TEST_BRANCH00

  return;
}








void Cube666_Digitizer::Digitize_whiteNoise(int pmtid){

#if TEST_BRANCH00
  assert(pmtid >=0 && pmtid < NPMT);

  Digitize_whiteNoise(pmtid,Cube666_Analysis::getInstance()->fRawData[pmtid]);
#endif //TEST_BRANCH00

  return;
}










void Cube666_Digitizer::Digitize_pmtRes(int* array,int pmtid,int nentries,float* rawDataArray){

  //array = int array[],
  //= number of photons detected in each time sample by PMT pmtid

  assert(pmtid >=0 && pmtid < NPMT);

  reset(rawDataArray,nentries);

  for(int i=0;i<nentries;i++){
    int nhits = array[i];
    //cout<<"pmt "<<pmtid<<"\t sample "<<i<<"\t nhist "<<nhits<<endl;
    for(int ii=0;ii<nhits;ii++){
      Cube666_PMTHit* hit;
      hit->setHitGlobalTime(i*TIMESAMPLE);
      hit->setHitPMTid(pmtid);
      hit->setIsDetected(1);
      Digitize_pmtRes(hit,rawDataArray);
    }  
  }


#if VERBOSE
  for(int i=0;i<nentries;i++){
    cout<<"pmtid "<<pmtid<<"\t sample "<< i << "\t rawData "<<rawDataArray[i]<<endl;
  }

#endif //VERBOSE


  return;
}









void Cube666_Digitizer::Digitize_pmtRes(int* array,int pmtid,int nentries){



#if TEST_BRANCH00
  //array = int array[],
  //= number of photons detected in each time sample by PMT pmtid

  Digitize_pmtRes(array,pmtid,nentries,Cube666_Analysis::getInstance()->fRawData[pmtid]);

#endif //TEST_BRANCH00

  return;
}











void Cube666_Digitizer::Digitize_pmtRes_fromTree(string filepath){

  TFile* inputFile = new TFile(filepath.c_str());
  assert(inputFile->IsOpen());

  TTree* inputTree = (TTree*)inputFile->Get("tree");
  if(!inputTree) inputTree = (TTree*)inputFile->Get("Data");
  assert(inputTree);

  inputTree->SetBranchStatus("*RawData*",0);

  string outputFilename = filepath.substr(0,filepath.rfind(".root")) + "_digitized.root";

  TFile* outputFile = new TFile(outputFilename.c_str(),"recreate");
  TTree* outputTree = inputTree->CloneTree(0);


  float fRawData[NPMT][NTIMESAMPLE];

  ostringstream branchname,leafname;
  for(int i=0;i<NPMT;i++){
    branchname.str("");
    leafname.str("");

    //branchname = "eRawData" (instead of "fRawData") just to be consistent with the naming in realData.
    branchname<<"eRawData"<<i;
    leafname<<"eRawData["<< fNTimeSamples <<"]/F";
    outputTree->Branch(branchname.str().c_str(),&fRawData[i],leafname.str().c_str());
  }
  

  int fNphotonPMT[NPMT][NTIMESAMPLE];
  for(int i=0; i < NPMT; i++){
    branchname.str("");
    branchname << "fNphotonPMT" << i;
    inputTree->SetBranchAddress(branchname.str().c_str(),&fNphotonPMT[i]);
  }


  int nevents = inputTree->GetEntries();

  for(int evi=0;evi < nevents; evi++){
    if(!(evi%1000)) cout<<"event "<<evi<<endl;
    //if(evi > 0) break;
    inputTree->GetEntry(evi);


    for(int pmti=0;pmti < NPMT;pmti++){
      //reset(fRawData[pmti],fNTimeSamples);
      Digitize_pmtRes(fNphotonPMT[pmti],pmti,fNTimeSamples,fRawData[pmti]);

#if DARK_CURRENT_DIGITIZATION

      Digitize_darkCurrent(pmti,fRawData[pmti]);

#endif //DARK_CURRENT_DIGITIZATION


#if WHITE_NOSIE_DIGITIZATION

      Digitize_whiteNoise(pmti,fRawData[pmti]);

#endif //WHITE_NOISE_DIGITIZATION

    }


    outputTree->Fill();
  }
  
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  delete inputFile;
  delete outputFile;

  return;
}
