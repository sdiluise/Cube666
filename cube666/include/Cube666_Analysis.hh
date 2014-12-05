#ifndef _ARDM_ANALYSIS_
#define _ARDM_ANALYSIS_ 1

#include "preparation.hh"

#include "globals.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"


#include "TTree.h"
#include "TFile.h"
#include "TArrayF.h"
#include "TClonesArray.h"
#include "TGraph2D.h"
#include "TRandom.h"



#include <vector>
using namespace std;



class Cube666_Analysis {

private:
  Cube666_Analysis();

  static Cube666_Analysis* fInstance;

  TFile* fOutputFile;
  TTree* fAnaTree;


  //************
  //lightmap = vector<vector<TGraph2D> >
  //1st index = pmtid
  //2nd index = z-position id
  //
  //
  //************
  //for each PMT pmtid, and each z-position zposid,
  //the TGraph2D stores the fraction of light detected by the PMT pmtid as a function of the xy-position of the interaction point.
  //this fraction of detected light can be calculated with respect to 
  //   1. either the total number of photons emitted (nPhotonDetected_pmtid / sum_nPhotonEmitted)
  //   2. or     the total number of photons detected by both PTMArrays together (nPhotonDetected_pmtid / (sum_nPhotonDetected_by_topArray + sum_nPhotonDetected_by_btmArray ) )
  //
  //************* the code implemented here assumes that option (1) be used. ********************
  //
  //
  //*************
  //how to get TGraph2D :
  //for each z-position, we move the source along an xy-grid on that z-plane,
  //at each xy-point, we simulate 1k events, 10k photons / event,
  //and count then the number of photons detected on each PMT in each event.
  //
  //
  //*****
  //for option (1) mentioned above, we proceed as follows :
  //
  //we plot then the distribution  nPhotonDetected_pmtid / 10k 
  //--> it's a gaussian distribution 
  //(note : this is a gaussian since the denominator is a constant ! if it were not a constant, it would be way more complicated !
  // the ratio Z = X / Y between two iid (independently identically distributed) gaussian variable X and Y is a cauchy distribution, or even more complicated
  // depending on the means of X,Y-distributions and their correlation !!!)
  //
  //--> take the mean to feed to TGraph2D, i.e. tgraph2d->SetPoint(n,xpos,ypos,gaussian_mean)
  //--> what we're after is the following :
  //
  //in later simulation, instead of tracing all the photons emitted in scintillation process, 
  //which can be computationally very expensive if the number of photons emitted is large,
  //we will use the lightmap to calculate the number of photons detected by each PMT in an event.
  //
  //in that case, on event-by-event basis, we would generate a random number according to the distribution of nPhotonDetected_pmtid / 10k
  //to get the fraction of photons (out of all emitted photons) that PMT pmtid would detect.
  //i.e. we would need the mean and the spread (sigma) of the nPhotonDetected_pmtid / 10k distribution.
  //
  //each xy-point has its own sigma.
  //
  //we plot all those sigma and fit the distribution with again a gaussian, take the mean of that gaussian, and store it in lightmap_err.
  //later whenever we need this sigma, we will make the assumption that all the xy-points have the same sigma if they lie on the same z-plane.
  //
  //--> the mean  is fed to lightmap (TGraph2D)
  //--> and sigma is fed to lightmap_err
  //
  //
  //
  //*****
  //for option (2) :
  //
  //we plot the distribution nPhotonDetected_pmtid / (sum_nPhotonDetected_by_topArray + sum_nPhotonDetected_by_btmArray )
  //
  //both the numerator and denominator are gaussian distributed. their ratio follows theoretically a cauchy (or even way more complicated) distribution !
  //
  //but in reality, the ratio-distribution can be very well approximated by another gaussian distribution.
  //so to save us all the theoretical fuss, we will just use a gaussian to fit this distribution.
  //
  //following the same philosophy explained above for option (1), 
  //we will store the mean of the distributions in lightmap, and the widths (sigma) in lightmap_err  
  //
  //
  //
  //
  //*******
  //together with vector<vector<TGraph2D> > lightmap and vector<vector<double> > lightmap_err.
  //we use vector<vector<double> > lightmap_zpos  
  //the index convention is the same for lightmap.
  //
  //*******
  //lightmap_zpos stores all the zpositions of the lightmap (grid) for individual pmt.
  //in normal case, all the 24 vector<double> stored in lightmap should be the same.
  //but still we employ this "representation" to cope with the contingency, that one of the zposition might be missing.
  //
  //
  //
  //



public:
  ~Cube666_Analysis();
  static G4String fFilename;

  //methods
  
  //get singleton
  static Cube666_Analysis* getInstance();

  void simulate_timeStructure(double Edep,G4ThreeVector interactionPos,double interaction_globalTime,double wval=W0_LAR,
			      double componentRatio=LAR_YIELDRATIO,double tau1=LAR_FASTTIMECONSTANT, double tau3=LAR_SLOWTIMECONSTANT);

  
  void BeginOfRun();
  void EndOfRun();
  
  void BeginOfEvent();
  void EndOfEvent();
  void Reset();
  void bookTree();

  
#if TEST_BRANCH0
  //data members
  G4double fTimeSample;
  G4int fNPMT;
  G4int fNevents;
  G4int fNtimeSample;
  G4int fNphotonDetected;   //total number of detected photons by the PMT array in 1 event
  G4int fNphotonEmitted;

  //truth position of the source where particles are emitted                                                                                                          
  G4double fPosx;
  G4double fPosy;
  G4double fPosz;


  G4int fTotNphotonPMT[NPMT]; //total number of detected photons by 1 PMT in 1 event


#if TEST_BRANCH00


  //restructuring anaTree

  //test
  //can't create a branch from an object of a user-defined class 
  //<-- this works in stand-alone ROOT-macro, but not in GEANT4 with embedded ROOT --> go figure !!
  //fail to build dictionaries for new classes
  //--> use the following extremely ugly code instead !

  //creating branches storing informations recorded by each PMT

  G4int fNphotonPMT[NPMT][NTIMESAMPLE]; //Nphotons detected by 1 PMT in 1 timesample in 1 event                                   
  float fRawData[NPMT][NTIMESAMPLE];    //this is similar to fNphotonPMT, but only makes sense when we digitize a photo-electron peak

  //recording direct light                                                                                                                                            
  //i.e. number of UV-photons falling onto PMTs' surfaces                                                                                                              

  G4int fNDirectPhoton;                    //total number of UV-photons falling onto PMT in one array in 1 event   
  G4int fTotNDirectPhotonPMT[NPMT];        //total number UV-photons falling onto 1 PMT in 1 event                                                              
  G4int fNPhotonHittingPMT;                //total number of photons falling onto PMT in one array in 1 event (visible + UV light !)                                
  G4int fTotNPhotonHittingPMT[NPMT];       //total number photons falling onto 1 PMT in 1 event                                                                

#endif //TEST_BRANCH00

#endif //TEST_BRANCH0


  //end test

};

#endif //_ARDM_ANALYSIS_

