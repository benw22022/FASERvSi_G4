//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file eventgenerator/HepMC/HepMCEx01/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "DetectorParameters.hh"
#include "G4PVDivision.hh"
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  G4int GetNlayers() const { return fNLayers; }
  G4int GetModulesPerLayer() const { return fModulesPerLayer; }

private:
  void DefineMaterials();

  G4Material* fAir;
  G4Material* fAr;
  G4Material* fSilicon;
  G4Material* fScinti;
  G4Material* fTungsten;
  G4Material* fVacuum;
  
  G4LogicalVolume* fTarget_log;
  G4LogicalVolume* fSCT_strip_log;
  G4LogicalVolume* fTruthTracker_log;

  std::vector<G4VPhysicalVolume*> fSD_phys;
  std::vector<G4VPhysicalVolume*> fTarget_phys;

  G4int fModulesPerLayer{8};
  G4int fNLayers{0};

  bool detectorCreated{false};

  void ConstructSDandField() override;
  
  DetectorConstructionMessenger* messenger;

};

#endif
