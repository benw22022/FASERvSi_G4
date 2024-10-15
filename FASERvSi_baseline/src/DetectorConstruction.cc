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
/// \file eventgenerator/HepMC/HepMCEx01/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "DetectorConstruction.hh"
#include "DetectorParameters.hh"
#include "TrackerParametrisation.hh"
#include "TrackerSD.hh"
#include "Detector.hh"
#include <string>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
    messenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
    delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterials()
{
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------
  G4NistManager* nistManager = G4NistManager::Instance();
  fAir = nistManager->FindOrBuildMaterial("G4_AIR");
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");

  G4double a, z, density;
  G4int nel;

  // Scintillator
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  fScinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  fScinti-> AddElement(elC, 9);
  fScinti-> AddElement(elH, 10);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
    
  //-------------------------------------------------------------------------
  // Materials
  //-------------------------------------------------------------------------
  DefineMaterials();

  //-------------------------------------------------------------------------
  // Detector geometry
  //-------------------------------------------------------------------------

  //------------------------------ experimental hall
  G4Box* experimentalHall_box = new G4Box("expHall_b", DetectorParameters::Get()->fexpHall_x, DetectorParameters::Get()->fexpHall_y, DetectorParameters::Get()->fexpHall_z);
  G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(experimentalHall_box, fAir,"expHall_L", 0,0,0);
  
  G4VPhysicalVolume * experimentalHall_phys = new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log, "expHall_P", 0, false,0);
  G4VisAttributes* experimentalHallVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));

  experimentalHallVisAtt-> SetForceWireframe(true);
  experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);

  // Loop over number of SCT layers assume; assume Tungsten and SCT modules are directly adjacent
  for (unsigned int i{0}; i < DetectorParameters::Get()->fnumSCTLayers; i++)
  {   
    
    //------------------------------ Targets (tungsten sheets)
    std::string targetLogVolName = "Target" + std::to_string(i+1) + "_log";
    std::string targetPhysVolName = "Target" + std::to_string(i+1) + "_phys";

    G4Box* Target_box = new G4Box("Target_box", DetectorParameters::Get()->fdetWidth/2, DetectorParameters::Get()->fdetHeight/2, DetectorParameters::Get()->ftungstenThickness/2);
    G4LogicalVolume* Target_log = new G4LogicalVolume(Target_box, fTungsten, "Target_log");
    Target_log->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
    
    G4double targetOffset = DetectorParameters::Get()->ftargetStartPosZ + i * (DetectorParameters::Get()->fSCTThickness + DetectorParameters::Get()->ftungstenThickness);

    G4VPhysicalVolume * Target_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, targetOffset), Target_log, targetPhysVolName, experimentalHall_log, false, 0);
    
    fTarget_log.push_back(Target_log);
    fTarget_phys.push_back(Target_phys);
    

    //------------------------------ Sensitive detectors
    std::string SDLogVolName = "SD" + std::to_string(i+1) + "_log";
    std::string SDPhysVolName = "SD" + std::to_string(i+1) + "_phys";

    G4Box* SD_box = new G4Box("SD_box", DetectorParameters::Get()->fdetWidth/2, DetectorParameters::Get()->fdetHeight/2, DetectorParameters::Get()->fSCTThickness/2);
    G4LogicalVolume* SD_log = new G4LogicalVolume(SD_box, fAir, "SD1_log");
    SD_log->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
    
    G4double SDOffset = DetectorParameters::Get()->ftargetStartPosZ + DetectorParameters::Get()->ftungstenThickness + i * (DetectorParameters::Get()->fSCTThickness + DetectorParameters::Get()->ftungstenThickness);
    G4VPhysicalVolume * SD_phys = new G4PVPlacement(0, G4ThreeVector(0,  0, SDOffset), SD_log, SDPhysVolName, experimentalHall_log, false, 0);

    fSD_log.push_back(SD_log);
    fSD_phys.push_back(SD_phys);
  }
    
  return experimentalHall_phys;
}


void DetectorConstruction::ConstructSDandField(){
  
  G4SDManager *sdman = G4SDManager::GetSDMpointer();

  for (unsigned int i{0}; i<fSD_phys.size(); i++)
  {
    std::string detName = "sensDet" + std::to_string(i);

    Detector* sensDet = new Detector(detName);
    fSD_log[i]->SetSensitiveDetector(sensDet);
    sdman->AddNewDetector(sensDet);
  }

}