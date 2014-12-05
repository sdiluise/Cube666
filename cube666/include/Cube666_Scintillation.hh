//inspired by the code for class DMArScint on the SVN

#ifndef _ARDM_SCINTILLATION_
#define _ARDM_SCINTILLATION_ 1


#include "globals.hh"
#include "G4Scintillation.hh"
#include "G4ProcessType.hh"
#include "G4VParticleChange.hh"

#include "vector"
using namespace std;

class Cube666_Scintillation : public G4Scintillation{

private:
  G4double NDriftedElectrons_neutronLike(G4double Edeposit);
  G4double NDriftedElectrons_electronLike(G4double Edeposit);
  G4double NPhotons_neutronLike(G4double Edeposit);
  G4double NPhotons_electronLike(G4double Edeposit);
  vector<G4double> quenchingFactor_electronLike();

  G4double ArDMScint_sample_time(G4double tau1,G4double tau2);
  G4double ArDMScint_single_exp(G4double t,G4double tau2); 
  G4double ArDMScint_bi_exp(G4double t,G4double tau1,G4double tau2);  
  G4double getScintillationYield(G4DynamicParticle* aparticle,
				 G4MaterialPropertiesTable* aMaterialPropertiesTable,
				 G4double TotalEnergyDeposit);


  //G4double NPhotons_GAr(G4double Edeposit,G4double excitationRatio);
  G4double NPhotons_GAr(G4double Edeposit);

public:
  Cube666_Scintillation(G4String procName="Cube666_Scint",G4ProcessType type=fElectromagnetic);
  ~Cube666_Scintillation();
  
  //redefine the virtual method PostStepDoIt(..) of class G4Scintillation
  //to take into account quenching effect

  G4VParticleChange* PostStepDoIt_tracing_opPhoton(const G4Track& aTrack, const G4Step& aStep);  
  G4VParticleChange* PostStepDoIt_using_lightmap(const G4Track& aTrack, const G4Step& aStep);  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);  

};








#endif //_ARDM_SCINTILLATION_

