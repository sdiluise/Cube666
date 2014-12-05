#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Cube666)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
#include_directories(${PROJECT_SOURCE_DIR}/include )

##tell cmake where to pick root-classes up
##/eth/sw/ROOT/5.30.00/i386-linux/root/include
#set(ROOT_header /Users/knguyen/work/geneve/software/root/root/include)

#ORIG
#set(ROOT_header $ENV{HOME}/work/geneve/software/root/root/include/)
#MY
set(ROOT_header /Applications/ROOT534_19/root/include)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_header})

#############
####test#####
#############

#include_directories(`root-config --cflags`)
#include_directories(`root-config --libs`)

#############
##end test###
#############

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

##pick root-headerfiles in /eth/sw/ROOT/5.30.00/i386-linux/root/include/
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${ROOT_header}/*.h)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#add_executable(exampletest exampletest.cc ${sources} ${headers})
add_executable(Cube666 Cube666.cc ${sources} ${headers})
#target_link_libraries(exampletest ${Geant4_LIBRARIES})

##link to root-library.
##/eth/sw/ROOT/5.30.00/i386-linux/root/lib
#target_link_libraries(exampletest ${Geant4_LIBRARIES} -L/Users/knguyen/work/geneve/software/root/root/lib)
##ROOT_lib = `root-config --libs`

#ORIG
#set(ROOT_lib -L$ENV{HOME}/work/geneve/software/root/root/lib/ -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -pthread -lm -ldl -rdynamic)

#MY
set(ROOT_lib -L/Applications/ROOT534_19/root/lib -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -pthread -lm -ldl -rdynamic)



target_link_libraries(Cube666 ${Geant4_LIBRARIES} ${ROOT_lib})


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build test. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Cube666_SCRIPTS
#   Cube666.in
#   Cube666.out
#  exampletest.in
#  exampletest.out
#  init.mac
#  init_vis.mac
#  run1.mac
#  run2.mac
#  vis.mac
  )

foreach(_script ${Cube666_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(test DEPENDS exampletest)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS Cube666 DESTINATION bin)


