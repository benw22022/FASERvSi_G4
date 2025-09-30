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

        // Geometry parameters: static so that they can be accessed without an instance of the class
        static constexpr G4double strips() { return fNstrips; }  // Number of strips per plane
        static constexpr G4double stereoAngle() { return fStereoAngle; }  // Stereo angle between the two strip planes in a module
        static constexpr G4double stripLength() { return fStripLength; }
        static constexpr G4double planeWidth() { return fPlaneWidth; }
        static constexpr G4double planeThickness() { return fPlaneThickness; } // 285 microns
        static constexpr G4double planeSeparation() { return fPlaneSeparation; } // Separation between the two strip planes in a module
        static constexpr G4double stripWidth() { return fPlaneWidth / fNstrips; } // Width of an individual strip
    
    private:
        //https://www.sciencedirect.com/science/article/pii/S016890020601388X?ref=pdf_download&fr=RR-2&rr=98740f01e8e5a935
        static constexpr G4int fNstrips = 770; // 768 active strips + first/last for the bias voltage
        static constexpr G4double fStereoAngle = 40 * mrad;
        static constexpr G4double fStripLength = 128.05*mm;
        static constexpr G4double fPlaneWidth = 63.56*mm;
        static constexpr G4double fPlaneThickness = 285*1e-3*mm;
        static constexpr G4double fPlaneSeparation = fPlaneThickness*2;  //From p52 of ID TDR https://cds.cern.ch/record/331063/files/ATLAS-TDR-4-Volume-I.pdf - "rectangles are then sandwiched together, 0.1 mm apart"
        // Paki says that this is impossible - should be > 0.3 mm

        static constexpr G4double fModuleHeight = 7.08*mm;

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