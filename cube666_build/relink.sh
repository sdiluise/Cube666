
G4LIB=/cluster/home04/phys/khnguyen/work/geneve/software/geant4/geant4.9.6.p02-install/lib64/Geant4-9.6.2/
CUBE666_SRC=/cluster/home04/phys/khnguyen/work/geneve/analysis/mc/ArDM/code/cube666/

cmake -DGEANT4_USING_OPENGL_X11=ON -DGeant4_DIR=$G4LIB $CUBE666_SRC

