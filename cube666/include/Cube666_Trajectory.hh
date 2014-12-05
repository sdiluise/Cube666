#ifndef _ARDM_TRAJECTORY_
#define _ARDM_TRAJECTORY_ 1

#include "globals.hh"

#include "G4VTrajectory.hh"
#include "G4TrajectoryPoint.hh"

#include "G4Track.hh"
#include "vector"
using namespace std;

typedef vector<G4TrajectoryPoint*> Cube666_TrajectoryPointContainer;


class Cube666_Trajectory : public G4VTrajectory{

  G4int         fTrackID;
  G4int         fParentID;
  G4String      fParticleName;
  G4double      fParticleCharge;
  G4int         fPDGEncoding;
  G4ThreeVector fInitialMomentum;
  G4ThreeVector fInitialMomentumDirection;
  G4String      fVolumeName;   
  
  Cube666_TrajectoryPointContainer* fPositionRecord;
  
  G4double      fEkin;
  G4double      fSteplength;

public:
  Cube666_Trajectory(); 
  Cube666_Trajectory(const G4Track* track);
  Cube666_Trajectory(Cube666_Trajectory&  trajectory); 
  ~Cube666_Trajectory();

  G4int         GetTrackID()                   const { return fTrackID;                  };
  G4int         GetParentID()                  const { return fParentID;                 }; 
  G4String      GetParticleName()              const { return fParticleName;             };
  G4double      GetCharge()                    const { return fParticleCharge;           };  
  G4int         GetPDGEncoding()               const { return fPDGEncoding;              };
  G4ThreeVector GetInitialMomentum()           const { return fInitialMomentum;          };
  G4ThreeVector GetInitialMomentumDirection()  const { return fInitialMomentumDirection; };
  G4int         GetPointEntries()              const { return fPositionRecord->size();   };
  G4VTrajectoryPoint* GetPoint(G4int i)        const { return (*fPositionRecord)[i];     };
  G4String      GetVolumeName()                const { return fVolumeName;               };

  Cube666_TrajectoryPointContainer* GetPositionRecord(){ return fPositionRecord;      };

  G4double      GetEkin()                      const { return fEkin;                     };
  G4double      GetSteplength()                const { return fSteplength;               };


  void          AppendStep(const G4Step* step);
  void          MergeTrajectory(G4VTrajectory* secondTrajectory);

  void*  operator new(size_t);
  void   operator delete(void*);
  G4bool operator ==(const G4VTrajectory& right) const{
    return (this==&right);
  };


};

#endif // _ARDM_TRAJECTORY_
