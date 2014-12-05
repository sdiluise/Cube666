//code copied from the old code in SVN !!

#include "preparation.hh"
#include "Cube666_Scintillation.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4PhysicsVector.hh"

#include "Cube666_Analysis.hh"
#include "Cube666_Test.hh"
#include "Cube666_Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"
#include "G4TrackVector.hh"

#define DEBUG false
#if DEBUG
#define D(x) cout<<x<<endl;
#else
#define D(x)
#endif


Cube666_Scintillation::Cube666_Scintillation(G4String procName,G4ProcessType type):
  G4Scintillation(procName,type)
{

    
  SetFiniteRiseTime(false);
  //SetTrackSecondariesFirst(true); 
  //<-- commenting this out reduces the run-time significantly
  

}



Cube666_Scintillation::~Cube666_Scintillation(){;}




G4VParticleChange* Cube666_Scintillation::PostStepDoIt_tracing_opPhoton(const G4Track& aTrack, const G4Step& aStep){


  //code copied from method G4Scintillation::AtRestDoIt(..),
  //only the part for calculating the number of emitted photon is changed
  //to take quenching effect into account.

  D("\t -->"<<__FILE__<<"::"<<__FUNCTION__);

  aParticleChange.Initialize(aTrack);

  D("\t\t part: "<<aTrack.GetDefinition()->GetParticleName());

  //scintillation process happens only in liquid or gaseous argon
  
  if( !(aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() == "LAr"
       || aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() == "GAr"
       )){


    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  


  D("\t DumpPhysicsTable ");

  //DumpPhysicsTable();
  
  D("\t DumpPhysicsTable -- end --");
  
  

  //test
  //retrieving Cube666_Analysis singleton
  //Cube666_Analysis* ana     = Cube666_Analysis::getInstance();
  //Cube666_Test*     testObj = Cube666_Test::getInstance();
  //Cube666_Analysis singleton retrieved
  //end test


  G4double TotalEnergyDeposit;

//   if(aTrack.GetParticleDefinition()->GetParticleName() == "neutron" &&
//      aStep.GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().contains("Elastic")){
    
//     TotalEnergyDeposit = aStep.GetPreStepPoint()->GetKineticEnergy() - aStep.GetPostStepPoint()->GetKineticEnergy();

//   }else if(aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName().contains("[0.0]"))
//     //skip all the ions which are products of elastic neutron-scatt. off nuclei (Ar40[0.0], Fe56[0.0], Ni ...)
//     return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
//   else
//     TotalEnergyDeposit = aStep.GetTotalEnergyDeposit(); 


    TotalEnergyDeposit = aStep.GetTotalEnergyDeposit(); 

    D(" total energy deposit "<<TotalEnergyDeposit);
  //getchar();
  //test
  //Edep < 2eV --> there shouldn't be any photon generated !
  //this should speed up the run-time a tiny bit
    if(TotalEnergyDeposit < 200.*eV){
      D(" exit");
      return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

  //end test

  
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material*        aMaterial = aTrack.GetMaterial();
  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double      t0 = pPreStepPoint->GetGlobalTime();

  G4ThreeVector mom = pPreStepPoint->GetMomentum();

  //cout<<"Track-Step: "<<aTrack.GetDefinition()->GetParticleName()<<" -pos- "<<x0.x()<<" "<<x0.y()<<" "<<x0.z()
  //   <<" -mom- "<<mom.x()/MeV<<" "<<mom.y()/MeV<<" "<<mom.z()/MeV<<" Edep: "<<TotalEnergyDeposit/MeV<<endl;  

  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();


  if(VERBOSE)
    cout<<"Cube666_Scint.::PostStepDoIt(..) .... aMaterialPropertiesTable "<<aMaterialPropertiesTable<<endl;

  if (!aMaterialPropertiesTable)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  
  const G4MaterialPropertyVector* Fast_Intensity = 
    aMaterialPropertiesTable->GetProperty("FASTCOMPONENT"); 
  const G4MaterialPropertyVector* Slow_Intensity =
    aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

  if (!Fast_Intensity && !Slow_Intensity )
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  
  G4int nscnt = 1;
  if (Fast_Intensity && Slow_Intensity) nscnt = 2;
  



  //start of modified part
  G4int modification=1;
  G4double MeanNumberOfPhotons = DEFAULTVALUE; 
  G4double YieldRatio = DEFAULTVALUE;
  G4double ResolutionScale = DEFAULTVALUE;
  G4MaterialPropertyVector *Scint_Yield_Vector = NULL;

  //start of modified code
  //calculate the MeanNumberOfPhotons in a different way.
  //the following code is copied from the SVN
  
  ResolutionScale    = aMaterialPropertiesTable->GetConstProperty("RESOLUTIONSCALE");
  
  
  //to-do : implement code for scintillation in LAr and GAr separately !
  
  //if(1 || aTrack.GetVolume()->GetName().contains("ArCol") /* || aTrack.GetVolume()->GetName() == "GArCol"*/){
  //     cout<<"volume "<<aTrack.GetVolume()->GetName()
  // 	<<"\t material "<<aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()
  // 	<<endl; getchar();
  
  //if(aTrack.GetVolume()->GetName() == "LArCol"){
  if(aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() == "LAr"){
    //scintillation happens only in liquid argon volume or gaseous argon volume
    //for the time being, implement it only for LArCol.
    
    if(aTrack.GetDefinition()->GetParticleName().contains("[0.0]")){//<-- from SVN, what is Ar40[0.0] exactly ??
      //neutron elastic interaction
      //neutron interacts with matter 
      //-->generates these generic ions : (??? <-- correct ??)
      //Ar40[0.0], Fe56[0.0], Ni , Cr ... (<-- from tank material !)
      //so Ar40[0.0] is the product of neutron interaction with argon atoms
      
      //number of photons emitted  (Quenching)
      MeanNumberOfPhotons = NPhotons_neutronLike(TotalEnergyDeposit);
      YieldRatio = SCINTILLATION_YIELD_RATIO; //WARP <-- taken from SVN-code

    }else if(aTrack.GetDefinition()->GetParticleName() == "e-" ||
	     aTrack.GetDefinition()->GetParticleName() == "gamma"){ //electron/gamma interaction
      MeanNumberOfPhotons = NPhotons_electronLike(TotalEnergyDeposit);		
      YieldRatio = aMaterialPropertiesTable->GetConstProperty("ELECTRONEXCITATIONRATIO");

    }
    //       else if(aTrack.GetDefinition()->GetParticleName() == "neutron"){ //neutron interaction
    // 	MeanNumberOfPhotons = NPhotons_neutronLike(TotalEnergyDeposit);		
    // 	YieldRatio = aMaterialPropertiesTable->GetConstProperty("NEUTRONEXCITATIONRATIO");
    //       }
    else{

      MeanNumberOfPhotons = TotalEnergyDeposit/(W_GAMMA);
      YieldRatio = .35; //<-- arbitrary number ! for the time being, this is just a "place-holder" !
    }
    
    //}else if(aTrack.GetVolume()->GetName() == "GArCol"){
  }else if(aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() == "GAr"){
    
    //for the time being, for all particle types, meanNumberOfPhotons = Edep / W-val in GAr
    if(aTrack.GetDefinition()->GetParticleName() == "alpha"  
       || aTrack.GetDefinition()->GetParticleName() == "e-"     
       || aTrack.GetDefinition()->GetParticleName() == "gamma" 
       || aTrack.GetDefinition()->GetParticleName().contains("[0.0]") //ion
       ){
      MeanNumberOfPhotons = NPhotons_GAr(TotalEnergyDeposit);		
      YieldRatio          = GAR_YIELDRATIO_ALPHA; //<-- make this a property of the material !
      //cout<<"MeanNumberOfPhotons "<<MeanNumberOfPhotons<<endl; getchar();
    }
    
  }
  
  //MeanNumberOfPhotons = ceil(MeanNumberOfPhotons);
  //cout<<aStep.GetTrack()->GetParticleDefinition()->GetParticleName()<<"\t MeanNumberOfPhotons "<<MeanNumberOfPhotons<<endl;


  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////

  //MeanNumberOfPhotons = 1;//TotalEnergyDeposit/(W_GAMMA);
  //YieldRatio = .35;
  

  //cout<<" MeanNumberOfPhotons "<<MeanNumberOfPhotons<<endl;
  

  if(MeanNumberOfPhotons < 0) return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  
  
  G4int NumPhotons;
  
  if (MeanNumberOfPhotons > 10.){
    G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
    NumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);

  }else{
    NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
  }
  

  if (NumPhotons <= 0){
    
    aParticleChange.SetNumberOfSecondaries(0);
    
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  aParticleChange.SetNumberOfSecondaries(NumPhotons);
  
  //
  //if (fTrackSecondariesFirst) {
  if( GetTrackSecondariesFirst() ){
    if (aTrack.GetTrackStatus() == fAlive ) aParticleChange.ProposeTrackStatus(fSuspend);
  }
  
  
  D("NPhotons: "<<NumPhotons);

#if TEST_BRANCH0
  
  //ana->fNphotonEmitted +=NumPhotons;

#endif //TEST_BRANCH0

  G4int materialIndex = aMaterial->GetIndex();
  
  D(" Material "<<materialIndex<<" "<<aMaterial->GetName());

  //return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  // Retrieve the Scintillation Integral for this material  
  // new G4PhysicsOrderedFreeVector allocated to hold CII's
  
  G4int Num = NumPhotons;
  

  for (G4int scnt = 1; scnt <= nscnt; scnt++) {
    
    G4double ScintillationTime = 0.*ns;
    G4double ScintillationRiseTime = 0.*ns;
    //
    G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;
    
   
    if (scnt == 1) {
      if (nscnt == 1) {
	if(Fast_Intensity){

	  ScintillationTime   = aMaterialPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
	  //
	  //if (fFiniteRiseTime) ScintillationRiseTime = aMaterialPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");
	  if (GetFiniteRiseTime()) ScintillationRiseTime = aMaterialPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");
	  //
	  //ScintillationIntegral = 0;//(G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
	  ScintillationIntegral = (G4PhysicsOrderedFreeVector*)((*fFastIntegralTable)(materialIndex));
	}

	if(Slow_Intensity){
	  ScintillationTime   = aMaterialPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");
	  //
	  //if (fFiniteRiseTime) ScintillationRiseTime = aMaterialPropertiesTable->GetConstProperty("SLOWSCINTILLATIONRISETIME");
	  if (GetFiniteRiseTime()) ScintillationRiseTime = aMaterialPropertiesTable->GetConstProperty("SLOWSCINTILLATIONRISETIME");
	  //ScintillationIntegral = 0;//(G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
	  ScintillationIntegral = (G4PhysicsOrderedFreeVector*)((*fSlowIntegralTable)(materialIndex));
	}

      }else {
	//G4double YieldRatio = aMaterialPropertiesTable->GetConstProperty("YIELDRATIO");

	//in the constructor of G4Scintillation, ExcitationRatio is set to 1.0
	// n gamma emitted in the fast comp
	if ( GetScintillationExcitationRatio()  == 1.0 ) Num = G4int (std::min(YieldRatio,1.0) * NumPhotons);
	else 	                      Num = G4int (std::min(GetScintillationExcitationRatio(),1.0) * NumPhotons);

	ScintillationTime   = aMaterialPropertiesTable->GetConstProperty("FASTTIMECONSTANT");
	if (GetFiniteRiseTime()) ScintillationRiseTime = aMaterialPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");
	//ScintillationIntegral =0;// (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
	ScintillationIntegral = (G4PhysicsOrderedFreeVector*)((*fFastIntegralTable)(materialIndex));
      }

    }else {
      Num = NumPhotons - Num;
      ScintillationTime   =   aMaterialPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");
      if (GetFiniteRiseTime()) ScintillationRiseTime = aMaterialPropertiesTable->GetConstProperty("SLOWSCINTILLATIONRISETIME");
      //ScintillationIntegral = 0;//(G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
      ScintillationIntegral = (G4PhysicsOrderedFreeVector*)((*fSlowIntegralTable)(materialIndex));
    }
    

    //cout<<"Scintillation Integral: "<<ScintillationIntegral<<endl;

    if (!ScintillationIntegral) continue;
    

    // Max Scintillation Integral
    //
    G4double CIImax = ScintillationIntegral->GetMaxValue();


    //cout<<"Num = "<<Num<<endl;//getchar();


    //test
    for (G4int i = 0; i < Num; i++) {

      G4double sampledEnergy = 9.68;

      
    }
    

    //
    for (G4int i = 0; i < Num; i++) {
      
      //
      // Determine photon energy
      //
      G4double CIIvalue = G4UniformRand()*CIImax;
      G4double sampledEnergy = ScintillationIntegral->GetEnergy(CIIvalue);
      
      // Generate random photon direction
      
      G4double cost = 1. - 2.*G4UniformRand();
      G4double sint = std::sqrt((1.-cost)*(1.+cost));
      
      G4double phi = twopi*G4UniformRand();
      G4double sinp = std::sin(phi);
      G4double cosp = std::cos(phi);
      
      G4double px = sint*cosp;
      G4double py = sint*sinp;
      G4double pz = cost;
      
      // Create photon momentum direction vector 
      
      G4ParticleMomentum photonMomentum(px, py, pz);

      // Determine polarization of new photon 
      
      G4double sx = cost*cosp;
      G4double sy = cost*sinp; 
      G4double sz = -sint;
      
      G4ThreeVector photonPolarization(sx, sy, sz);
      
      G4ThreeVector perp = photonMomentum.cross(photonPolarization);
      
      phi = twopi*G4UniformRand();
      sinp = std::sin(phi);
      cosp = std::cos(phi);
      
      photonPolarization = cosp * photonPolarization + sinp * perp;
      
      photonPolarization = photonPolarization.unit();
      
      // Generate a new photon:
      
      G4DynamicParticle* aScintillationPhoton =	new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),photonMomentum);
      aScintillationPhoton->SetPolarization(photonPolarization.x(),photonPolarization.y(),photonPolarization.z());
      aScintillationPhoton->SetKineticEnergy(sampledEnergy);
      
      // Generate new G4Track object:
      
      G4double rand;
      
      if (aParticle->GetDefinition()->GetPDGCharge() != 0) rand = G4UniformRand();
      else                                                 rand = 1.0;
      
      G4double delta = rand * aStep.GetStepLength();

      G4double deltaTime = delta / ((pPreStepPoint->GetVelocity()+ pPostStepPoint->GetVelocity())/2.);
      
      // emission time distribution
      if (ScintillationRiseTime==0.0) deltaTime = deltaTime - ScintillationTime * std::log( G4UniformRand() );

      //deltaTime = deltaTime + sample_time(ScintillationRiseTime, ScintillationTime);
      else deltaTime = deltaTime + ArDMScint_sample_time(ScintillationRiseTime, ScintillationTime);

      G4double aSecondaryTime = t0 + deltaTime;

      G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();
      
      G4Track* aSecondaryTrack = new G4Track(aScintillationPhoton,aSecondaryTime,aSecondaryPosition);
      

      //cout<<"SecondaryTrack: "<<i<<" "<<aSecondaryTrack->GetDefinition()->GetParticleName()<<" -pos- "<<x0.x()<<" "<<x0.y()<<" "<<x0.z()<<" Kin:"<<sampledEnergy/eV<<endl;  

      aSecondaryTrack->SetTouchableHandle(aStep.GetPreStepPoint()->GetTouchableHandle());
      // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);
      aSecondaryTrack->SetParentID(aTrack.GetTrackID());

      aParticleChange.AddSecondary(aSecondaryTrack);

      //cout<<"NSEC: "<<aParticleChange.GetNumberOfSecondaries()<<endl;
    }

  }//ncnt
  
  D("\t "<<__FUNCTION__<<" end return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep)");
 
  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}






G4VParticleChange* Cube666_Scintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){

  D(__FUNCTION__);

  return PostStepDoIt_tracing_opPhoton(aTrack,aStep);

}









G4double Cube666_Scintillation::NDriftedElectrons_neutronLike(G4double Edeposit){//Edeposit given in keV !!

  if(ELECTRIC_FIELD_STRENGTH == 0){
    if(VERBOSE) 
      cout<<"Cube666_Scintillation::NDriftedElectrons_neutronLike .... "
	  <<"ELECTRIC_FIELD_STRENGTH == 0 ... return 0."
	  <<endl;

    return 0.;
  }

  //ionization due to slow ions
  //values obtained after fit to data following Lindhard theory
  //WARP proposal
  G4double a_lindhard     = A_LINDHARD;
  G4double alpha_lindhard = ALPHA_LINDHARD;
  G4double nuclearQuenching = a_lindhard*pow(Edeposit/keV,alpha_lindhard);
  
  G4double w0 = W0;
  G4double NElec = nuclearQuenching*Edeposit/w0; 
  
  //recombination : box model of Thomas and Imel assumed
  G4double electric_field_strength = ELECTRIC_FIELD_STRENGTH;
  G4double Xi = CX/electric_field_strength;
  NElec *= log(1+Xi)/Xi;

  if(VERBOSE) 
    cout<<"in Cube666_Scintillation::NDriftedElectrons_neutronLike : "
	<<"NElec "<<NElec
	<<"\t quenchingFactor "<<nuclearQuenching
	<<"\t Edep/W0 "<<Edeposit/w0
	<<"\t Efield/(V/cm) "<<electric_field_strength/(volt/cm)
	<<"\t Xi "<<Xi
	<<endl;
  return NElec;
}


G4double Cube666_Scintillation::NPhotons_neutronLike(G4double Edeposit){

  //approximation:
  //actual number of photons emitted nphotons = number of excited Ar - number of e- drifting away from e-ion pairs


  //q = quenching factor
  //number of excited Ar : Nex = q*Edeposit/W_GAMMA

  //number of e-ion pairs: Npairs = q*Edeposit/W0
  //number of electrons pulled out of e-ion pairs due to E-field : Ne = Npairs * E/CX * ln(1 + CX/E)
  //CX = 1856 (cm/keV) <-- WARP proposal.


  //--> nphotons = Nex - Ne


  G4double ndriftedElectrons = NDriftedElectrons_neutronLike(Edeposit); 

  //in SVN-code, parameters needed for calculating nuclearQuenching 
  //are different in 2 methods: DMArScint::NumberOfElectrons(..) and DMArScint::NumberOfPhotons(..)
  //<-- why ??

  G4double w_gamma = W_GAMMA;
  G4double a_lindhard = A_LINDHARD;
  G4double alpha_lindhard = ALPHA_LINDHARD;
  G4double nuclearQuenching   = a_lindhard*pow(Edeposit/keV,alpha_lindhard);
  G4double nphotons = nuclearQuenching*Edeposit/w_gamma - ndriftedElectrons;

  if(VERBOSE)
  cout<<"in Cube666_Scintillation::NPhotons_neutronLike(..) ... "
      <<"Edep [keV] "<<Edeposit/keV
      <<"\t ndriftedElectrons "<<ndriftedElectrons
      <<"\t nphotons "<<nphotons
      <<"\t nuclearQuenching "<<nuclearQuenching
      <<"\t w_gamma [keV] "<<w_gamma/keV
      <<endl;

  return nphotons;



  //another question apart from the discrepency in calculating the quenching factor mentioned above:
  //according to section 6.2 in Lilian Kaufmann's dissertation:
  //
  //nphotons = nexcitations - ndriftedElectrons
  //
  //nexcitations = quenchingFactor*Edeposit/W_GAMMA <-- implicit assumption:
  //100% of the deposited energy (after quenching) is used to excite Ar atoms
  //
  //ndriftedElectrons is calculated as in Cube666_Scintillation::NDriftedElectrons_neutronLike(..)
  //with the, again implicit, assumption that
  //100% of the deposited energy (after quenching) is used to create e-ion pairs !
  //
  //<-- does it make sense ???
  //am i misunderstanding the calculation in the SVN-code / in Kaufmann's thesis !??
  //


  //<--- answer: i'm misunderstanding the meaning of W_GAMMA and W0 !!
  //W_GAMMA and W0 are not the excitation or ionisation energy of Ar, respectively !!
  //
  //the deposited energy is used to i. excite Ar-atoms and ii. ionize Ar-atoms !
  //the 2 W-values are determined statistically (when the excitation and ionization happen at the same time) 
  //and therefore does take into account
  //the "interference" of the 2 effects !!
  //W0 is therefore larger then the ionisation energy of a single Ar-atom !
}



vector<G4double> Cube666_Scintillation::quenchingFactor_electronLike(){
  vector<G4double> q;

  //dEdx in MeV/cm
  G4double dEdx[] = {15.,9.,6.7,5.4,4.6,4.1,3.7,3.4,3.1,2.9}; //10 entries ! <-- from SVN, why these numbers ???
  G4int nentries = sizeof(dEdx)/sizeof(G4double);

//   if(ELECTRIC_FIELD_STRENGTH == 0){
//     for(int i=0;i<nentries;i++) q.push_back(1.);
//     return q;
//   }

  //constants taken from Lilian Kaufmann's thesis
  G4double A = .8;
  G4double k = .0486 *kilovolt/MeV; 
  G4double electric_field_strength = ELECTRIC_FIELD_STRENGTH;

  //dEdx in MeV/cm
  for(int i=0;i<nentries;i++){
    //G4double qi = A/(1 + k/(kilovolt/MeV)*dEdx[i]/(electric_field_strength/(kilovolt/cm)));
    //G4double qi = A/(1 + k/(kilovolt/MeV)*dEdx[i]/(electric_field_strength/(volt/cm)));
    G4double qi = A/(1 + k*dEdx[i]*MeV/cm/electric_field_strength);
    q.push_back(qi);
  }

  return q;
}



G4double Cube666_Scintillation::NDriftedElectrons_electronLike(G4double Edeposit){//Edeposit given in keV !!
  //copied from SVN-code

  //lindhard:
  //quenching factor q = A/(1 + k/E*dEdx)
  //A = .8
  //k = .00486 kV/MeV
  //dEdx energy loss
  

  //[k] = kV/MeV
  //[E] = kV/cm
  //[dEdx] = MeV/cm

  //dEdx changes with E --> calculate q, and then number of drifted electrons, over energy steps !!

  G4double Erest = Edeposit;
  G4double Estep = 10. *keV;
  G4double w0 = W0;

  vector<G4double> q = quenchingFactor_electronLike();

  G4int i=0;
  G4double quenching=0;

  if(Erest<= 100. *keV){
    while(Erest >= Estep){//this should happen 9 times !
      quenching += q[i];
      i++;
      Erest -= Estep;
    }
    return (Estep*quenching + Erest*q[i])/w0;

  }else return (Edeposit*q[q.size()-1]/w0);
  
}



G4double Cube666_Scintillation::NPhotons_electronLike(G4double Edeposit){//Edeposit given in keV !!
//   //i don't understand the SVN code for this method !!

//   //G4double NPhot = Er/WG;
//   //G4double nele0= Er/W0; //Maximum number of electrons, no quenching!
//   // G4double kie=0.; //Parameter not used so far...
//   //  NPhot = NPhot - (NElectrons - (kie*nele0));//reduction due to drift field
//   G4double WG_2 = 25.0e-6;  //<-- what is this energy ?? 
//   G4double NPhot = Er/WG_2; //<-- why ?? what's to do with the number of electrons drifting away ??.
//   return NPhot;

  G4double efield = ELECTRIC_FIELD_STRENGTH;
  G4double w_gamma = W_GAMMA;
  if(!efield) return Edeposit/w_gamma;

  G4double w0 = W0;

  vector<G4double> q = quenchingFactor_electronLike();
  G4int i=0;
  G4double quenching=0;

//   G4double Erest = Edeposit;
//   G4double Estep = 10. *keV;

//   G4double nphotons = 0;
// //   if(Erest<= 100. *keV){
// //     while(Erest >= Estep){//this should happen 9 times !
// //       quenching += q[i];
// //       i++;
// //       Erest -= Estep;
// //     }
// //     nphotons = (Estep*quenching + Erest*q[i])*(1./W_GAMMA - 1./W0);

// //   }else nphotons = (Edeposit*q[q.size()-1]*(1./W_GAMMA - 1./W0));

//   if(Erest<= 100. *keV){
//     while(Erest >= Estep){//this should happen 9 times !
//       nphotons += q[i]*(Estep/eV)*(1./(W_GAMMA/eV) - 1./(W0/eV));
//       i++;
//       Erest -= Estep;
//     }
//     nphotons += Erest/eV*q[i]*(1./(W_GAMMA/eV) - 1./(W0/eV));

//   }else nphotons = (Edeposit/eV*q[q.size()-1]*(1./(W_GAMMA/eV) - 1./(W0/eV)));



  //test
  G4double Erest = Edeposit;
  G4double Estep = 10. *keV;

  G4double nphotons = 0;

  if(Erest<= 100. *keV){
    while(Erest >= Estep){//this should happen 9 times !
      nphotons += q[i]*Estep*(1./w_gamma - 1./w0);
      i++;
      Erest -= Estep;
    }
    nphotons += Erest*q[i]*(1./w_gamma - 1./w0);

  }else nphotons = Edeposit*q[q.size()-1]*(1./w_gamma - 1./w0);



  //end test


  return nphotons;
}



G4double Cube666_Scintillation::ArDMScint_sample_time(G4double tau1,G4double tau2){
  //the code for this function is exactly the same as the one of the private (!!) method G4Scintillation::sample_time(G4double,G4double)
  //it is used in the virtual method G4Scintillation::PostStepDoIt(..)

  //i modify this method in the new class Cube666_Scintillation (derived from G4Scintillation)
  //in PostStepDoIt(..), i want to use G4Scintillation::sample_time(..)
  //but it's not possible since the method was declared as private in the mother class G4Scintillation
  //and therefore not accessible from the daughter class Cube666_Scintillation (<-- is it correct ??)
  //that's why i create another method which is identical with G4Scintillation::sample_time(G4double,G4double)

  //tau1: rise time
  //tau2: decay time

  cout<<"*******************"<<endl;

  while(1){
    //2 random numbers
    G4double ran1 = G4UniformRand();
    G4double ran2 = G4UniformRand();

    //
    //exponential distribution as envelope function: very efficient
    //
    G4double d = (tau1+tau2)/tau2;
    //make sure the envelope function is
    //always larger than the bi-exponential
    G4double t = -1.0*tau2*std::log(1-ran1);
    G4double g = d*ArDMScint_single_exp(t,tau2);
    if(ran2 <= ArDMScint_bi_exp(t,tau1,tau2)/g) return t;
  }
  return -1.0;
}


G4double Cube666_Scintillation::ArDMScint_single_exp(G4double t,G4double tau2){
  //used in ArDMScint_sample_time
  //the same reason as explained in ArDMScint_sample_time
  return std::exp(-1.0*t/tau2)/tau2;
}



G4double Cube666_Scintillation::ArDMScint_bi_exp(G4double t,G4double tau1,G4double tau2){
  //used in ArDMScint_sample_time
  //the same reason as explained in ArDMScint_sample_time
  return std::exp(-1.0*t/tau2)*(1-std::exp(-1.0*t/tau1))/tau2/tau2*(tau1+tau2);
}




G4double Cube666_Scintillation::getScintillationYield(G4DynamicParticle* aParticle,
						   G4MaterialPropertiesTable* aMaterialPropertiesTable,
						   G4double TotalEnergyDeposit){ 
  
  
  G4double ScintillationYield = 0.;
  if (GetScintillationByParticleType()) {
    // The scintillation response is a function of the energy
    // deposited by particle types.
    
    // Get the definition of the current particle
    G4MaterialPropertyVector *Scint_Yield_Vector = NULL;
    G4ParticleDefinition *pDef = aParticle->GetDefinition();
    
    // Obtain the G4MaterialPropertyVectory containing the
    // scintillation light yield as a function of the deposited
    // energy for the current particle type
      
    if(pDef==G4Proton::ProtonDefinition())  
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("PROTONSCINTILLATIONYIELD");
    
    else if(pDef==G4Deuteron::DeuteronDefinition())
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("DEUTERONSCINTILLATIONYIELD");
    
    else if(pDef==G4Triton::TritonDefinition())
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("TRITONSCINTILLATIONYIELD");
    
    else if(pDef==G4Alpha::AlphaDefinition())
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("ALPHASCINTILLATIONYIELD");
    
    // Ions (particles derived from G4VIon and G4Ions)
    // and recoil ions below tracking cut from neutrons after hElastic
    else if(pDef->GetParticleType()== "nucleus" || 
	    pDef==G4Neutron::NeutronDefinition()) 
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("IONSCINTILLATIONYIELD");
    
    // Electrons (must also account for shell-binding energy
    // attributed to gamma from standard PhotoElectricEffect)
    else if(pDef==G4Electron::ElectronDefinition() ||
	    pDef==G4Gamma::GammaDefinition())
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("ELECTRONSCINTILLATIONYIELD");
    
    // Default for particles not enumerated/listed above
    else
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("ELECTRONSCINTILLATIONYIELD");
    
    // If the user has not specified yields for (p,d,t,a,carbon)
    // then these unspecified particles will default to the 
    // electron's scintillation yield
    if(!Scint_Yield_Vector){
      Scint_Yield_Vector = aMaterialPropertiesTable->GetProperty("ELECTRONSCINTILLATIONYIELD");
    }
    
    // Throw an exception if no scintillation yield is found
    if (!Scint_Yield_Vector) {
      G4cerr << "\nG4Scintillation::PostStepDoIt(): "
	     << "Request for scintillation yield for energy deposit"
	     <<" and particle type without correct entry in MaterialPropertiesTable\n"
	     << "ScintillationByParticleType requires at minimum that ELECTRONSCINTILLATIONYIELD is set by the user\n"
	     << G4endl;
      G4Exception("G4Scintillation::PostStepDoIt","No correct entry in MaterialPropertiesTable",
		  FatalException,"Missing MaterialPropertiesTable entry.");
    }
    
    if (verboseLevel>1) {
      G4cout << "\n"
	     << "Particle = " << pDef->GetParticleName() << "\n"
	     << "Energy Dep. = " << TotalEnergyDeposit/MeV << "\n"
	     << "Yield = " 
	     << Scint_Yield_Vector->Value(TotalEnergyDeposit) 
	//<< Scint_Yield_Vector->GetProperty(TotalEnergyDeposit) 
	     << "\n" << G4endl;
    }
    
    // Obtain the scintillation yield using the total energy
    // deposited by the particle in this step.
    
    // Units: [# scintillation photons]
    ScintillationYield = Scint_Yield_Vector->Value(TotalEnergyDeposit);
  } else {
    // The default linear scintillation process
    ScintillationYield = aMaterialPropertiesTable->GetConstProperty("SCINTILLATIONYIELD");
    
    // Units: [# scintillation photons / MeV]
    //MY added 
    ScintillationYield *= GetScintillationYieldFactor();
  }
  
  return ScintillationYield;
}




G4double Cube666_Scintillation::NPhotons_GAr(G4double Edeposit){

  G4double w0     = W_GAMMA_GAR;

  //test
  return Edeposit/w0;
  //end test


  G4double Efield = ELECTRIC_FIELD_STRENGTH;
  G4double nphoton = 0;



  if(Efield == 0.){
    //if Efield == 0
    //nphotons = Edposit/w0;
    nphoton = Edeposit/w0;
  
  }else{
    //Efield = 0 --> no quenching by external field
    //Edeposit is expended for :
    //Ni   ionisations --> with following recombination and finally emission of UV photon
    //Nex  excitations --> with following de-excitation and finally emission of (infrared + UV <-- ??) photon 
    //kinetic energy of electrons kicked out of Ar-atoms in the above Ni ionisations
    //
    //--> Edeposit = Ni*Ei + Nex*Eex + Ni*Ekin_e
    //
    //we're interested in the number of UV photons emitted
    //
    //according to "total ionization in gases by high-energy particles : an appraisal of our understanding"
    //http://www.sciencedirect.com/science/article/pii/0020708X61901089 -- R.L. Platzman
    //(section 3. magnitudes of W in pure gases -- (A) noble gases)
    //
    //about 18% of the absorbed energy goes into Ekin_e (--> heat)
    //
    //--> we can now calculate Ni if we know the ratio Nex/Ni = excitationRatio
    //
    //
  }


  return nphoton;
}













