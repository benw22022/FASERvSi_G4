#ifndef SCTModuleGeometry_hh
#define SCTModuleGeometry_hh

#include <vector>
#include <algorithm>

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4PVDivision.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"

class SCTModuleGeometry {
  
    public:
        SCTModuleGeometry();
        ~SCTModuleGeometry(){};

        G4Box* GetModuleBox() const { return fModuleBoundingBox; };
        G4LogicalVolume* GetModuleLogical() const { return fModule_log; };
        G4LogicalVolume* GetStripLogical() const { return fStrip_indiv_log; };
        G4LogicalVolume* GetStripPlaneLogical() const { return fStrip_plane_log; };
        G4LogicalVolume* GetTruthTrackerPlaneLogical() const { return fTruthTrackerPlane_log; };
        G4VPhysicalVolume* PlaceModule(G4LogicalVolume* mother_log, const G4ThreeVector& position, G4RotationMatrix* rotation, G4int copyNo) const;

    private:
        G4int fNstrips = 768;
        G4double fStereoAngle = 40 * mrad;
        G4double fStripLength = 128.05*mm;
        G4double fPlaneWidth = 63.56*mm;
        G4double fPlaneThickness = 285*1e-3*mm;
        G4double fPlaneSeparation = 6.51*mm;

        G4Box* fModuleBoundingBox;
        G4LogicalVolume* fTruthTrackerPlane_log;
        G4LogicalVolume* fModule_log;
        G4LogicalVolume* fStrip_plane_log;
        G4LogicalVolume* fStrip_indiv_log;
        G4PVDivision* fStrip_div;

        //TODO Make materials centrally defined
        G4Material* fAir = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
        G4Material* fSilicon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
};


#endif