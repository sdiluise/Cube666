#ifndef _ARDM_DETECTORCONSTRUCTION_
#define _ARDM_DETECTORCONSTRUCTION_ 1


#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VSensitiveDetector.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"

using namespace std;

class Cube666_DetectorConstruction : public G4VUserDetectorConstruction{

protected:
  G4NistManager* fNist;
  G4Material*    fVacuum;
  G4Material*    fWorldMat;
  G4Material*    fTankMat;
  G4Material*    fLAr;
  G4Material*    fGAr;
  G4Material*    fPMTMat;
  G4Material*    fWLSMat;
  G4Material*    fCathodeGridMat;
  G4Material*    fPMTCathodeMat;
  G4Material*    fPMTBaseMat;


  G4VPhysicalVolume* fWorldPhys;
  G4VPhysicalVolume* fTankPhys; //put 6 walls together to make a tank
  G4VPhysicalVolume* fLArColPhys;
  G4VPhysicalVolume* fGArColPhys;
  vector<G4VPhysicalVolume*> fCathodeGridPhys;
  vector<G4VPhysicalVolume*> fPMTArrayPhys;
  vector<G4VPhysicalVolume*> fPMTCoatArrayPhys;
  vector<G4VPhysicalVolume*> fPMTCathodeArrayPhys;



  G4LogicalBorderSurface* fLAr_GAr_surf;
  G4LogicalBorderSurface* fGAr_LAr_surf;
  G4LogicalBorderSurface* fLAr_tank_surf;
  G4LogicalBorderSurface* fGAr_tank_surf;
  vector<G4LogicalBorderSurface*> fPMT_PMTCoat_surf;
  vector<G4LogicalBorderSurface*> fPMTCoat_PMT_surf;
  vector<G4LogicalBorderSurface*> fLAr_PMTCoat_surf;
  vector<G4LogicalBorderSurface*> fPMTCoat_LAr_surf;
  vector<G4LogicalBorderSurface*> fLAr_PMT_surf;
  vector<G4LogicalBorderSurface*> fPMT_LAr_surf;


  
  G4VPhysicalVolume* constructWorld(G4Material* fWorldMat,G4double world_half_size);
  void addTank();
  void addLArColumn();
  void addGArColumn();
  void addPMT();
  void addPMTCoat();
  void addPMTCathode();
  void addCathodeGrid();
  
  void build_surfaces();

  G4LogicalBorderSurface* build_LAr_tank_surf();
  G4LogicalBorderSurface* build_GAr_tank_surf();
  G4LogicalBorderSurface* build_LAr_GAr_surf(G4String order="LAr_GAr"); 
  vector<G4LogicalBorderSurface*> build_LAr_PMTCoat_surf(G4String order="LAr_PMTCoat");
  vector<G4LogicalBorderSurface*> build_PMT_PMTCoat_surf(G4String order="PMT_PMTCoat");
  vector<G4LogicalBorderSurface*> build_LAr_PMT_surf(G4String order="LAr_PMT");


  
  vector<G4ThreeVector> setCathodeWireVector(G4String xyAxis);
  G4LogicalVolume* constructPMT(G4Material* fPMTMat);
  vector<G4VPhysicalVolume*> placePMT(G4LogicalVolume* fPMTLog,G4VPhysicalVolume* fMotherPhys,vector<G4ThreeVector> rPMT);
  


  //void setMatPropTab_detMat(); //LAr
  G4Material* getTankMat(); //stainless steel

  void setMatPropTab_WLS(); //TPB
  void setMatPropTab_LAr(); //LAr
  void setMatPropTab_GAr(); //GAr
  //set the scintillation properties of argon (liquid and gas)
  void setArScintProperty(G4MaterialPropertiesTable* propTab,G4String medium="LAr");
  void setMatPropTab_PMTMat(); //PMTMat
  G4Material* getFR4();

  G4double GArRefIndex(G4double wavelength); //wavelength in micrometer !!
  G4double LArRefIndex(G4double wavelength); //wavelength in micrometer !!
  G4double GArEpsilon(G4double wavelength);  //wavelength in micrometer !!
  G4double LArEpsilon(G4double wavelength);  //wavelength in micrometer !!
  G4double RayleighAttenuationLength_LAr(G4double wavelength); //wavelength in micrometer !!
  G4double RayleighAttenuationLength_GAr(G4double wavelength); //wavelength in micrometer !!

  G4double getPhotonWavelength(G4double photonEnergy); //the returned wavelength is measured in micrometer !!

  G4VSensitiveDetector* getSD();


public:
  Cube666_DetectorConstruction();
  ~Cube666_DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  
};







#endif // _ARDM_DETECTORCONSTRUCTION_
