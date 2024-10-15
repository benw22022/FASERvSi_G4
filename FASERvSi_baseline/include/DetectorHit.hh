#ifndef EXN04DETECTORHIT_HH
#define EXN04DETECTORHIT_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4Types.hh"

#include <vector>


class DetectorHit : public G4VHit {
public:
  DetectorHit();
  ~DetectorHit(){};
  /// Draw pixels
  void Draw();
  ///// Get hit ID calculated as 1000 * sensorID + cellID
  //G4int ID() { return 1000 * fCopyNumSensor + fCopyNumCell; }
  inline void SetPosition(G4double x, G4double y, G4double z) {
    fPosX = x;
    fPosY = y;
    fPosZ = z;
  }

  inline void SetPDGID(G4int pdgid){
    fPDGID = pdgid;
  }

  inline void SetEnergy(G4double E){
    fEnergy = E;
  }

  /// Get hit X position
  inline G4double GetX() const { return fPosX; }
  /// Get hit Y position
  inline G4double GetY() const { return fPosY; }
  /// Get hit Z position
  inline G4double GetZ() const { return fPosZ; }
  /// Get hit pdgID
  inline G4double GetPDGID() const { return fPDGID; }
  /// Get hit Energy
  inline G4double GetEnergy() const { return fEnergy; }


private:
  /// Position along x axis
  G4double fPosX = -1;
  /// Position along y axis
  G4double fPosY = -1;
  /// Position along z axis
  G4double fPosZ = -1;
  /// PDGID
  G4int fPDGID = -999;
  /// Energy
  G4double fEnergy = -999.;

};

typedef G4THitsCollection<DetectorHit> DetectorHitCollection;

#endif /* EXN04DETECTORHIT_HH */
