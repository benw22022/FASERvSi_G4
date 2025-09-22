#include "SCTModuleGeometry.hh"
#include <cmath>


SCTModuleGeometry::SCTModuleGeometry()
{   
    // Calculate bounding box dimensions
    G4double x_bounds = fStripLength * sin(fStereoAngle) + fPlaneWidth * cos(fStereoAngle);
    G4double y_bounds = fStripLength * cos(fStereoAngle) + fPlaneWidth * sin(fStereoAngle);
    G4double z_bounds = (fPlaneThickness * 2) + fPlaneSeparation;

    fModuleBoundingBox = new G4Box("module_bounding_box", x_bounds/2, y_bounds/2, z_bounds/2);  // Oriented with long edge in y-direction
    fModule_log = new G4LogicalVolume(fModuleBoundingBox, fAir, "module_log", 0,0,0);

    G4Box* fStrip_plane = new G4Box("strip_plane", fPlaneWidth/2, fStripLength/2, fPlaneThickness/2);  // Oriented with long edge in y-direction
    G4Box* fStrip_indiv = new G4Box("strip_indiv", fPlaneWidth/2, (fPlaneWidth/fNstrips)/2, fPlaneThickness/2);  // Oriented with long edge in y-direction
    fStrip_plane_log = new G4LogicalVolume(fStrip_plane, fSilicon, "strip_plane_log", 0,0,0);
    fStrip_indiv_log = new G4LogicalVolume(fStrip_indiv, fSilicon, "strip_inidv_log", 0,0,0);
    fStrip_div = new G4PVDivision("strip_div", fStrip_indiv_log, fStrip_plane_log, kXAxis, fNstrips, 0 );
    
    // Create a tracking plane that we can set to be senstive to record truth hits
    G4Box* fTruthTrackerPlane = new G4Box("truth_tracker_plane", x_bounds/2, y_bounds/2, 0.1*mm);  // Very thin box
    fTruthTrackerPlane_log = new G4LogicalVolume(fTruthTrackerPlane, fAir, "truth_tracker_plane_log", 0,0,0);
    G4VPhysicalVolume* truth_tracker_plane_phys = new G4PVPlacement(
        0, 
        G4ThreeVector(0, 0, 0), 
        fTruthTrackerPlane_log, 
        "truth_tracker_plane_phys", 
        fModule_log, 
        false, 
        0);
        
    // Place the two strip planes inside the module bounding box
    G4VPhysicalVolume* strip_plane_side1_phys = new G4PVPlacement(
        0, 
        G4ThreeVector(0, 0, -fPlaneThickness/2 - fPlaneSeparation/2), 
        fStrip_plane_log, 
        "SCT_front_phys", 
        fModule_log, 
        true, 
        0);

    // Apply stereo angle rotation around Z-axis for the back strip plane
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateZ(fStereoAngle);  
    G4VPhysicalVolume* strip_plane_side2_phys = new G4PVPlacement(
        rot, 
        G4ThreeVector(0, 0, fPlaneThickness/2 + fPlaneSeparation/2), 
        fStrip_plane_log, 
        "SCT_back_phys", 
        fModule_log,
        false, 
        1);

    // Set visibility attributes for bounding box
    G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour::Brown());
    boxVisAtt->SetForceWireframe(true);
    boxVisAtt->SetVisibility(false);
    fModule_log->SetVisAttributes(boxVisAtt);
    fTruthTrackerPlane_log->SetVisAttributes(boxVisAtt);

    // Set visibility attributes for silicon planes
    G4VisAttributes* strip_planeVisAtt = new G4VisAttributes(G4Colour::Green());
    // strip_planeVisAtt->SetForceWireframe(true);
    strip_planeVisAtt->SetForceSolid(true);
    fStrip_plane_log->SetVisAttributes(strip_planeVisAtt);

    // Make individual strips invisible in visualization (too many to display usefully)
    G4VisAttributes* strip_invis = new G4VisAttributes();
    strip_invis->SetVisibility(false);
    fStrip_indiv_log->SetVisAttributes(strip_invis);
}

G4VPhysicalVolume* SCTModuleGeometry::PlaceModule(G4LogicalVolume* mother_log, const G4ThreeVector& position, G4RotationMatrix* rotation, G4int copyNo) const
{
    G4VPhysicalVolume* module_phys = new G4PVPlacement(
        rotation,
        position,
        fModule_log,
        "SCT_module_phys",
        mother_log,
        false,
        copyNo);
    return module_phys;
}