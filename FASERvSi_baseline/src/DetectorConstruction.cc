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
#include "G4GDMLParser.hh"

#include "SCTModuleGeometry.hh"

#include "DetectorConstruction.hh"
#include "DetectorParameters.hh"
#include "SCTModuleDetector.hh"


#include <string>
#include <fstream>

void checkOverlaps(G4VPhysicalVolume* physvol)
{
  if (physvol->CheckOverlaps())
  {
    std::string message = "Overlap detected in " + physvol->GetName();
    G4Exception("G4PVPlacement::CheckOverlaps()", physvol->GetName(), RunMustBeAborted, message.c_str());
  }
}



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

  // Vacuum
  G4double universe_mean_density = 1.e-25*g/cm3;
  G4Element* elN  = new G4Element("Nitrogen","N",  z=7.,  a= 14.00674*g/mole);
  G4Element* elO  = new G4Element("Oxygen",  "O",  z=8.,  a= 15.9994*g/mole);
  fVacuum = new G4Material("universe_mean_density", 1.e-25*g/cm3, nel=2);
  fVacuum-> AddElement(elN, .7);
  fVacuum-> AddElement(elO, .3);

  // Scintillator
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  fScinti = new G4Material("Scintillator", density= 1.032*g/cm3, nel=2);
  fScinti-> AddElement(elC, 9);
  fScinti-> AddElement(elH, 10);
}

/*
Construct tracking layer with modules oriented with long edge in y-direction
____________________________
|     ||     ||     ||     |
|  1  ||  2  ||  3  ||  4  |
|     ||     ||     ||     |
|_____||_____||_____||_____|
|     ||     ||     ||     |
|   5 ||  6  ||  7  ||  8  |
|     ||     ||     ||     |
|_____||_____||_____||_____|
*/

G4LogicalVolume* constructVertTrackingLayerLogical(SCTModuleGeometry& sctModule)
{
  G4Box* sct_module_box = sctModule.GetModuleBox();
  G4Box* tracking_layer_box = new G4Box("tracking_layer_b",
     sct_module_box->GetXHalfLength()*4 + 0.1*mm, 
     sct_module_box->GetYHalfLength()*2 + 0.1*mm, 
     sct_module_box->GetZHalfLength()+0.1*mm);

  G4LogicalVolume* tracking_layer_log = new G4LogicalVolume(tracking_layer_box, 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "tracking_layer_log", 0,0,0);
  
  G4ThreeVector m1_transl = G4ThreeVector(-3*sct_module_box->GetXHalfLength(), sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m2_transl = G4ThreeVector(-sct_module_box->GetXHalfLength(), sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m3_transl = G4ThreeVector(sct_module_box->GetXHalfLength(), sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m4_transl = G4ThreeVector(3*sct_module_box->GetXHalfLength(),sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m5_transl = G4ThreeVector(-3*sct_module_box->GetXHalfLength(), -sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m6_transl = G4ThreeVector(-sct_module_box->GetXHalfLength(), -sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m7_transl = G4ThreeVector(sct_module_box->GetXHalfLength(), -sct_module_box->GetYHalfLength(), 0);
  G4ThreeVector m8_transl = G4ThreeVector(3*sct_module_box->GetXHalfLength(), -sct_module_box->GetYHalfLength(), 0);

  sctModule.PlaceModule(tracking_layer_log, m1_transl, 0, 1);
  sctModule.PlaceModule(tracking_layer_log, m2_transl, 0, 2);
  sctModule.PlaceModule(tracking_layer_log, m3_transl, 0, 3);
  sctModule.PlaceModule(tracking_layer_log, m4_transl, 0, 4);
  sctModule.PlaceModule(tracking_layer_log, m5_transl, 0, 5);
  sctModule.PlaceModule(tracking_layer_log, m6_transl, 0, 6);
  sctModule.PlaceModule(tracking_layer_log, m7_transl, 0, 7);
  sctModule.PlaceModule(tracking_layer_log, m8_transl, 0, 8);

  G4VisAttributes* tracking_layerVisAtt = new G4VisAttributes(G4Colour::Blue());
  tracking_layerVisAtt->SetForceWireframe(true);
  tracking_layerVisAtt->SetVisibility(false);
  tracking_layer_log->SetVisAttributes(tracking_layerVisAtt);

  return tracking_layer_log;
}

/*
Construct tracking layer with modules oriented with long edge in x-direction
----------------------------
|            |             |
|      1     |      2      |
----------------------------
|            |             |
|      3     |      4      |
----------------------------
|            |             |
|      5     |     6       |
----------------------------
|            |             |
|      7     |     8       |
----------------------------
*/
G4LogicalVolume* constructHozTrackingLayerLogical(SCTModuleGeometry& sctModule)
{
  G4Box* sct_module_box = sctModule.GetModuleBox();
  G4Box* tracking_layer_box = new G4Box("tracking_layer_b",
     sct_module_box->GetYHalfLength()*2 + 0.1*mm, 
     sct_module_box->GetXHalfLength()*4 + 0.1*mm, 
     sct_module_box->GetZHalfLength()+0.1*mm);

  G4LogicalVolume* tracking_layer_log = new G4LogicalVolume(tracking_layer_box, 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "tracking_layer_log", 0,0,0);
  
  G4RotationMatrix* rot90deg = new G4RotationMatrix();
  rot90deg->rotateZ(90 * deg);
  G4ThreeVector m1_transl = G4ThreeVector(-sct_module_box->GetYHalfLength(), 3*sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m2_transl = G4ThreeVector(sct_module_box->GetYHalfLength(), 3*sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m3_transl = G4ThreeVector(-sct_module_box->GetYHalfLength(), sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m4_transl = G4ThreeVector(sct_module_box->GetYHalfLength(), sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m5_transl = G4ThreeVector(-sct_module_box->GetYHalfLength(), -sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m6_transl = G4ThreeVector(sct_module_box->GetYHalfLength(), -sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m7_transl = G4ThreeVector(-sct_module_box->GetYHalfLength(), -3*sct_module_box->GetXHalfLength(), 0);
  G4ThreeVector m8_transl = G4ThreeVector(sct_module_box->GetYHalfLength(), -3*sct_module_box->GetXHalfLength(), 0);

  sctModule.PlaceModule(tracking_layer_log, m1_transl, rot90deg, 1);
  sctModule.PlaceModule(tracking_layer_log, m2_transl, rot90deg, 2);
  sctModule.PlaceModule(tracking_layer_log, m3_transl, rot90deg, 3);
  sctModule.PlaceModule(tracking_layer_log, m4_transl, rot90deg, 4);
  sctModule.PlaceModule(tracking_layer_log, m5_transl, rot90deg, 5);
  sctModule.PlaceModule(tracking_layer_log, m6_transl, rot90deg, 6);
  sctModule.PlaceModule(tracking_layer_log, m7_transl, rot90deg, 7);
  sctModule.PlaceModule(tracking_layer_log, m8_transl, rot90deg, 8);

  G4VisAttributes* tracking_layerVisAtt = new G4VisAttributes(G4Colour::Blue());
  tracking_layerVisAtt->SetForceWireframe(true);
  tracking_layerVisAtt->SetVisibility(false);
  tracking_layer_log->SetVisAttributes(tracking_layerVisAtt);

  return tracking_layer_log;
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

  //* experimental hall
  G4Box* experimentalHall_box = new G4Box("expHall_b", DetectorParameters::Get()->fexpHall_x/2, DetectorParameters::Get()->fexpHall_y/2, DetectorParameters::Get()->fexpHall_z/2);
  G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(experimentalHall_box, fAir,"expHall_L", 0,0,0);
  G4VPhysicalVolume* experimentalHall_phys = new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log, "expHall_P", 0, false,0);
  G4VisAttributes* experimentalHallVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
  experimentalHallVisAtt->SetForceWireframe(true);
  experimentalHall_log->SetVisAttributes(experimentalHallVisAtt);

  //* SCT module and tracking layers
  SCTModuleGeometry sctModule = SCTModuleGeometry();
  // fSCT_strip_log = sctModule.GetStripLogical();
  fSCT_strip_log = sctModule.GetTruthTrackerPlaneLogical();
  G4LogicalVolume* tracking_hoz_layer_log = constructHozTrackingLayerLogical(sctModule);
  G4LogicalVolume* tracking_vert_layer_log = constructVertTrackingLayerLogical(sctModule);
  G4Box* tracking_layer_box = dynamic_cast<G4Box*>(tracking_vert_layer_log->GetSolid());
  G4double tracking_layer_thickness = 2*tracking_layer_box->GetZHalfLength();
  G4double tungsten_thickness = 2*DetectorParameters::Get()->ftungstenThickness;  
  
  //* Tungsten target
  G4Box* Target_box = new G4Box("Target_box", DetectorParameters::Get()->fdetWidth/2, DetectorParameters::Get()->fdetHeight/2, DetectorParameters::Get()->ftungstenThickness/2);
  G4LogicalVolume* Target_log = new G4LogicalVolume(Target_box, fTungsten, "Target_log");
  G4VisAttributes* TargetVisAtt =  new G4VisAttributes(G4Colour::Red());
  TargetVisAtt->SetForceWireframe(true);
  Target_log->SetVisAttributes(TargetVisAtt);

  //* Place layers and targets 
  G4double pos = DetectorParameters::Get()->ftargetStartPosZ;
  G4double target_mass = 0*g;
  for (unsigned int i{0}; i < DetectorParameters::Get()->fnumSCTLayers; i++)
  {
    G4VPhysicalVolume* Target_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, pos), Target_log, "Target_phys", experimentalHall_log, false, i);
    auto solid = Target_log->GetSolid();
    auto material = Target_log->GetMaterial();
    G4double volume = solid->GetCubicVolume();
    G4double density = material->GetDensity();  // g/cm³
    G4double mass = density/(g/cm3) * volume/cm3;  // grams
    target_mass = target_mass + mass;

    checkOverlaps(Target_phys);

    G4VPhysicalVolume* SD_phys;
    pos += tungsten_thickness/2 + tracking_layer_thickness/2;
    if (i%2 == 0)
    {
      SD_phys = new G4PVPlacement(0, G4ThreeVector(0,  0, pos), tracking_hoz_layer_log, "HozLayer_phys", experimentalHall_log, false, i);
    }
    else
    {
      SD_phys = new G4PVPlacement(0, G4ThreeVector(0,  0, pos), tracking_vert_layer_log, "VertLayer_phys", experimentalHall_log, false, i);
    }
    pos += tungsten_thickness/2 + tracking_layer_thickness/2;
    checkOverlaps(SD_phys);
  }

  // Print mass and length of the detector
  G4cout << "Detector length = " << pos - DetectorParameters::Get()->ftargetStartPosZ << " mm" << G4endl;
  G4cout << "Tungsten target mass = " << target_mass << " g" << G4endl;

  //* GDML dump
  G4GDMLParser* gdmlParser = new G4GDMLParser();  
  std::remove("FASERvSi.gdml"); // delete file
  gdmlParser->Write("FASERvSi.gdml", experimentalHall_phys, true);
  delete gdmlParser;

  return experimentalHall_phys;
}


void DetectorConstruction::ConstructSDandField(){
  
  G4SDManager *sdman = G4SDManager::GetSDMpointer();
  std::string detName = "strip_detector";
  SCTModuleDetector* sensDet = new SCTModuleDetector(detName);

  G4cout << "Attaching sensitive detector to logical volume " << fSCT_strip_log->GetName() << " at address " << fSCT_strip_log << G4endl;
  fSCT_strip_log->SetSensitiveDetector(sensDet);
  sdman->AddNewDetector(sensDet);

}