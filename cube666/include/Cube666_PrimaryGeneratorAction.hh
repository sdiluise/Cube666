#ifndef _PRIMARY_GENERATORACTION_
#define _PRIMARY_GENERATORACTION_ 1


#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "preparation.hh"

class Cube666_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{

private:
  G4ParticleGun* fParticleGun;
  G4ThreeVector getRandomPolarization(G4ThreeVector particleMomentumDirection);

  static Cube666_PrimaryGeneratorAction* fInstance;

public:
  Cube666_PrimaryGeneratorAction();
  ~Cube666_PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event* event);

  static Cube666_PrimaryGeneratorAction* getInstance();
  G4GeneralParticleSource* fParticleSource;

};




#endif //_PRIMARY_GENERATORACTION_ 
