1. you may want to change some paths in 

cube666/CMakeLists.txt  

and in 


cube666/include/preparation.hh  (search for _HOME_ in there.)

before running. (you can search for "sdiluise" to see which paths need to be changed.)




2. there are 2 convenient scripts for rebuilding / relinking (e.g. after adding a new class) :

rebuild.sh
relink.sh

you just have to source (one of ) them . After changing the paths in the cmake file, you can source the rebuild.sh script and it will rebuild and compile for you.




3. to compile :

cd cube666_build
gmake -j4




4. to run :

cd cube666_build
./Cube666 0 run.mac


the 1st argument governs the visualization, 1 = visualization on, 0 = vis. off.
you can change the run.mac macro to change the particle you emit or change the position, volume, surface, …. where they are emitted.



5. You will also have to set the environment variables pointing to the data files, otherwise it may complain. following is a snippet of G4 env. variables from my home on the brutus cluster just to get an idea where the data files may reside :

G4LEVELGAMMADATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/PhotonEvaporation2.3
G4NEUTRONXSDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4NEUTRONXS1.2
G4LEDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4EMLOW6.32
G4NEUTRONHPDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4NDL4.2
G4RADIOACTIVEDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/RadioactiveDecay3.6
G4PIIDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4PII1.3
G4SAIDXSDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/G4SAIDDATA1.1
G4REALSURFACEDATA=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/share/Geant4-9.6.2/data/RealSurface1.0


Setup the environment according to you directory PATHs

This is an example script you can find in cube666_build/ConfigData.sh:

export G4LEDATA=/Users/sdiluise/geant4.10.00.p02-build/data/G4EMLOW6.35
export G4LEVELGAMMADATA=/Users/sdiluise/geant4.10.00.p02-build/data/PhotonEvaporation3.0/
export G4NEUTRONXSDATA=/Users/sdiluise/geant4.10.00.p02-build/data/G4NEUTRONXS1.4/
export G4SAIDXSDATA=/Users/sdiluise/geant4.10.00.p02-build/data/G4SAIDDATA1.1/
export G4NEUTRONHPDATA=/Users/sdiluise/geant4.10.00.p02-build/data/G4NDL4.4/
