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
        //https://pdf.sciencedirectassets.com/271580/1-s2.0-S0168900206X0709X/1-s2.0-S016890020601388X/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEAYaCXVzLWVhc3QtMSJIMEYCIQD4NL5X0Ae2zKd98BywUh%2FqUFnIGGaNhmEaf%2FHA3KqUIQIhAMPyK27OLraGi%2BL4H6ejEGlp%2F15Wvb1fRDhm05FDArFwKrsFCI%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQBRoMMDU5MDAzNTQ2ODY1IgzJ2XQZiZT%2BNBFPPJAqjwWQjk3dbRcjYWxdDNzHhXjbqLhy9linLjKNA6H46sftf8odO9uIH0tdrTB1uldjf10c9g%2FQiQGdCUSzzSJAQwOm4n3%2FYk4G9cWTGjyRMDIfsz4iBiyCTcQRuY65KR9JgDLD1mwsuKlFLLXzsO9OBmRp%2BwTXK1AVnDUCnm0vUIiDfccktp%2B4wMFB63ffduvVPZYIdxFE2ymQRHkw7NXvqJwrRWSuf1ox1nfF4eY1fA7pSCZjVETO9CGBHw2oHzZi%2FF4gluTBMKLS%2FcLYpEbu2bQyTOf7twZEHnBYxYwl5DXQqi%2FeaKDnbppCnS%2FQuUzcQvyClc9HvRy3pU2G8TRf3udqaqDxJ2sr1dVs14NUUtm5IL%2FjNlfyM7VhDCpamoQTl1e1dHIQuqpHu4ZKV2cpcuuTrZF3n4msVgfTci4c8C8tnUcdwb7AhuKsXofZYUuhjNfCl%2FkkezxeXNPH5V4QGC5eFiytX3MTFGljB0RK9mMVjLnTjOurAaBfwQaC9Bch5XMOStSUTwduNvLy725JN5zNGjOsa0phaMEr2wFPntCdjdgBvrLJSZ0s3OZDPlPejtMmatUwlRwn3e34d6a5IYb1XiE76iaq3s%2Fr%2BWxHbbDbAL0QcF5mbDTB2qkJ%2BveugJbzzozSvLsYy3t4svqZneFFOWIY73h1Xg0irHX2hxjMWopQWanhciKhl%2FHZFpmK9hgqed0AIMsKl1p63pdUWJ5tBuA02wihlz3FTw3MXg1lvCAlTeQGwnJpIFfpMuqRT%2Bv5%2Fgwy53sYweaHbWQrEnSKC2ehlPHiSRGiwymo0W8VtsR04iWEfAIHwXXIDZhaUVAkwER8I%2BVzPJSd5p90TysDetHl7ez0ZgB58e5DDdveMOG52sYGOrABzPVu1dbTBZ6I5XJ3KKjFO7wAujuF1GRhxCyvSEkHenT9U7U%2BYNPxNUWHnWtpO50r0O3EEQhtr9TptZ2VHsQ54gaKHN%2BuE7XvmXLq%2FrUcZ18Ld8InLf%2BEuhLZmOfGho0dWULklNrq7FBdECa4XIaeZeQBqfWBFK2qy1NSyvs7B5nzNQFtZNbF%2FpizRhON1Bz2IZskXEW69708Oax7hGmoqLPY8FQJlkvCgCSR63PKjYc%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250926T150221Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYURJHOZKX%2F20250926%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=5d9f8657664c87f50622aa814b881250f03381a17479183fa4019d3818620e20&hash=e6fb7217b34f3923002a11cc646038f5258f0b980251cad566496424688f10f2&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S016890020601388X&tid=spdf-ca231d76-bb2f-47d0-8058-6e21ddcf2bc8&sid=14c994ae6eb7f64f9b095f400535f8c7436agxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=1d005d5950500e5a575e&rr=9853a451a9cc362e&cc=gb
        static constexpr G4int fNstrips = 770; // 768 active strips + first/last for the bias voltage
        static constexpr G4double fStereoAngle = 40 * mrad;
        static constexpr G4double fStripLength = 128.05*mm;
        static constexpr G4double fPlaneWidth = 63.56*mm;
        static constexpr G4double fPlaneThickness = 285*1e-3*mm;
        static constexpr G4double fPlaneSeparation = 0.1*mm;  //From p52 of ID TDR https://cds.cern.ch/record/331063/files/ATLAS-TDR-4-Volume-I.pdf - "rectangles are then sandwiched together, 0.1 mm apart"
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