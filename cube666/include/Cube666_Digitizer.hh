#ifndef _ARDM_DIGITIZER_
#define _ARDM_DIGITIZER_



#include "G4VDigi.hh"
#include "G4VDigitizerModule.hh"
#include "preparation.hh"
#include "Cube666_PMTHit.hh"
#include "G4HCofThisEvent.hh"

using namespace std;

class Cube666_Digitizer : public G4VDigitizerModule {

  static int fNTimeSamples;

  static double* fPMTCalib;
  static double* fPMTCalibSigma;
  static double* fPMTDarkCountRate;
  static double* fPMTWhiteNoiseSigma;

  void Digitize_simple(Cube666_PMTHit* hit); //simple digitization, just record global time of the hit
  void Digitize_simple(Cube666_PMTHitsCollection* hitcol); //simple digitization, just record global time of the hit

  void Digitize_pmtRes(Cube666_PMTHit* hit); //take pmtResponse into account. simulate analog signal for a hit.
  void Digitize_pmtRes_try1(Cube666_PMTHit* hit,float* array); //take pmtResponse into account. simulate analog signal for a hit.
  void Digitize_pmtRes_try2(Cube666_PMTHit* hit,float* array); //take pmtResponse into account. simulate analog signal for a hit.
  void Digitize_pmtRes(Cube666_PMTHit* hit,float* array); //take pmtResponse into account. simulate analog signal for a hit.

  void Digitize_pmtRes(Cube666_PMTHitsCollection* hitscol); //take pmtResponse into account. simulate analog signal for all hits.

  void Digitize_darkCurrent(int pmtid);
  void Digitize_darkCurrent(int pmtid,float* rawDataArray);

  void Digitize_whiteNoise(int pmtid);
  void Digitize_whiteNoise(int pmtid,float* rawDataArray);

  void Digitize_pmtRes(int* array, int pmtid, int nentries = NTIMESAMPLE); //arry = int array[] = number of photons detected in each timesample by 1 PMT

  void Digitize_pmtRes(int* array, int pmtid, int nentries, float* rawDataArray); //arry = int array[] = number of photons detected in each timesample by 1 PMT

public :

  Cube666_Digitizer(G4String name="PMTDigitizer",G4String pmtCalibFilename="");
  ~Cube666_Digitizer();

  void reset(float* array, int nentries);
  //void reset();

  void Digitize_simple(G4String hitsColName);
  void Digitize_pmtRes(G4String hitsColName);

  void Digitize_simple(G4HCofThisEvent* hc);
  void Digitize_pmtRes(G4HCofThisEvent* hc);
  
  void Digitize_pmtRes_fromTree(string filepath);
  void Digitize(G4HCofThisEvent* hc);

  void Digitize();

  static void setNTimeSamples(int ntimeSamples){ fNTimeSamples = ntimeSamples;};
};





#endif // _ARDM_DIGITIZER_
