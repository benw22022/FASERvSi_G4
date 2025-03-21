#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(FASERvSi_baseline)

# If we have a conda env activated - look in there for the dependancies
if(DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Using Conda environment at: $ENV{CONDA_PREFIX}")
  include_directories("$ENV{CONDA_PREFIX}/include")
  link_directories("$ENV{CONDA_PREFIX}/lib")
  set(CMAKE_BUILD_RPATH "$ENV{CONDA_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH "$ENV{CONDA_PREFIX}/lib")
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")

  # Find and setup Expat
  find_package(EXPAT REQUIRED)
  set(EXPAT_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
  set(EXPAT_LIBRARY "$ENV{CONDA_PREFIX}/lib/libexpat.so")
  include_directories(${EXPAT_INCLUDE_DIR})

  # Find and setup Motif
  set(MOTIF_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
  set(MOTIF_LIBRARIES "$ENV{CONDA_PREFIX}/lib/libXm.so")
  include_directories(${MOTIF_INCLUDE_DIR})
  link_directories(${MOTIF_LIBRARIES})
endif()


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
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find HepMC3 (required package)
#
find_package(HepMC3 REQUIRED)

##----------------------------------------------------------------------------
## Find Pythia6 (optional package)
##
#find_package(Pythia6 QUIET)
#if(PYTHIA6_FOUND)
#  message(STATUS "G4 Examples: Pythia6 found. --> FASERvSi_baseline example with Pythia6 enabled.") 
#  add_definitions(-DG4LIB_USE_PYTHIA) 
#else()
#  set(PYTHIA6_LIBRARIES "")  
#endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${HEPMC3_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(FASERvSi_baseline FASERvSi_baseline.cc ${sources} ${headers})

if(DEFINED ENV{CONDA_PREFIX})
  target_link_libraries(FASERvSi_baseline ${Geant4_LIBRARIES}
                        ${HEPMC3_LIBRARIES} ${HEPMC3_FIO_LIBRARIES} ${HEPMC3_LIB}
                        ${PYTHIA6_LIBRARIES} ${EXPAT_LIBRARY} "$ENV{CONDA_PREFIX}/lib/libstdc++.so")
else()
  target_link_libraries(FASERvSi_baseline ${Geant4_LIBRARIES}
                        ${HEPMC3_LIBRARIES} ${HEPMC3_FIO_LIBRARIES} ${HEPMC3_LIB}
                        ${PYTHIA6_LIBRARIES} ${EXPAT_LIBRARY})
endif()


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build FASERvSi_baseline. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(FASERvSi_baseline
    init_vis.mac vis.mac
  )

foreach(_script ${FASERvSi_baseline})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS FASERvSi_baseline DESTINATION bin)

