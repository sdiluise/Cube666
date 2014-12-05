G4BUILD=/Users/sdiluise/cube666/cube666_build/
G4INSTALL=/Users/sdiluise/geant4.10.00.p02-install/
CUBE666_SRC=/Users/sdiluise/cube666/cube666/

cd $G4BUILD

#cmake -DCMAKE_INSTALL_PREFIX=$G4INSTALL $CUBE666_SRC
cmake -DCMAKE_INSTALL_PREFIX=$G4INSTALL -DGEANT4_USE_QT=ON $CUBE666_SRC

make -j4 VERBOSE=1

make install

