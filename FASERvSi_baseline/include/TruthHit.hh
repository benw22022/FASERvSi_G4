#ifndef TruthHit_HH
#define TruthHit_HH

#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include  "G4LorentzVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Colour.hh"

#include <vector>
#include <iostream>


class TruthHit : public G4VHit {
public:
  TruthHit(){};
  ~TruthHit(){};

  inline void* operator new(size_t);
  inline void operator delete(void*);
  friend std::ostream& operator<<(std::ostream& os, const TruthHit& hit);
  bool operator<(const TruthHit& other) const;
  void Draw() override;
  
  
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

  inline void SetCharge(G4double charge){
    fCharge = charge;
  }

  inline void SetPx(G4double px){
    fPx = px;
  }
  
  inline void SetPy(G4double py){
    fPy = py;
  }

  inline void SetPz(G4double pz){
    fPz = pz;
  }
  
  inline void SetT(G4double t){
    fT = t;
  }
  inline void SetMass(G4double mass){
    fMass = mass;
  }
  inline void SetTrackID(G4int trackID){
    fTrackID = trackID;
  }
  inline void SetParentID(G4int parentID){
    fParentID = parentID;
  }
  
  inline void SetModuleNumber(G4int moduleNumber){
    fModuleNumber = moduleNumber;
  }
  inline void SetLayerNumber(G4int layerNumber){
    fLayerNumber = layerNumber;
  }

  /// Get hit X position
  inline G4double GetX() const { return fPosX; }
  /// Get hit Y position
  inline G4double GetY() const { return fPosY; }
  /// Get hit Z position
  inline G4double GetZ() const { return fPosZ; }
  /// Get hit time
  inline G4double GetT() const { return fT; }
  /// Get hit pdgID
  inline G4double GetPDGID() const { return fPDGID; }
  /// Get hit Energy
  inline G4double GetEnergy() const { return fEnergy; }
  /// Get hit Charge
  inline G4double GetCharge() const { return fCharge; }
  /// Get hit px
  inline G4double GetPx() const { return fPx; }
  /// Get hit py
  inline G4double GetPy() const { return fPy; }
  /// Get hit pz
  inline G4double GetPz() const { return fPz; }
  /// Get hit mass
  inline G4double GetMass() const { return fMass; }
  /// Get hit track ID
  inline G4int GetTrackID() const { return fTrackID; }
  /// Get hit parent ID
  inline G4int GetParentID() const { return fParentID; }
  /// Get module number within a layer (0 to 8) 
  inline G4int GetModuleNumber() const { return fModuleNumber; }
  /// Get layer number (0 to NTrackingLayers)
  inline G4int GetLayerNumber() const { return fLayerNumber; }

  inline G4ThreeVector GetTrackVertex() const {return fTrackVertex;};
  inline G4LorentzVector GetTrackP4() const {return fTrackP4;};
  inline G4int GetIsPrimaryTrack() const { return fIsPrimaryTrack; }
  inline G4int GetIsSecondaryTrack() const { return fIsSecondaryTrack; }

  inline void SetTrackVertex(const G4ThreeVector& vertex) { fTrackVertex = vertex; }
  inline void SetTrackP4(const G4LorentzVector& p4) { fTrackP4 = p4; }
  inline void SetIsPrimaryTrack(G4int isPrimary) { fIsPrimaryTrack = isPrimary; }
  inline void SetIsSecondaryTrack(G4int isSecondary) { fIsSecondaryTrack = isSecondary; }
  inline void SetColour(const G4Colour& colour){ fColour = colour; }
  inline G4Colour GetColour() const {return fColour; }

  inline void SetTruthHitID(G4long id) { fTruthHitID = id; }
  inline G4long GetTruthHitID() const { return fTruthHitID; }

private:
  /// Position along x axis
  G4double fPosX = -1;
  /// Position along y axis
  G4double fPosY = -1;
  /// Position along z axis
  G4double fPosZ = -1;
  /// Get time of hit
  G4double fT = -1;
  /// PDGID
  G4int fPDGID = -999;
  /// Energy
  G4double fEnergy = -999.;
  //Charge
  G4double fCharge = -999.;
  // Theta
  G4double fTheta = -999.;
  // Momentum
  G4double fPx = -999.;
  G4double fPy = -999.;
  G4double fPz = -999.;
  // Mass
  G4double fMass = -999.;
  // Track ID
  G4int fTrackID = -999;
  // Parent ID
  G4int fParentID = -999;
  
  G4int fModuleNumber = -1;    // Module number within a layer (0 to 8)
  G4int fLayerNumber = -1;     // Number of the tracking layer (0 to NTrackingLayers)

  G4ThreeVector fTrackVertex{-999., -999., -999.};
  G4LorentzVector fTrackP4{0,0,0,0};
  G4int fIsPrimaryTrack = 0; // 1 if primary track, 0 otherwise
  G4int fIsSecondaryTrack = 0; // 1 if secondary track, 0 otherwise
  G4long fTruthHitID = -1;

  G4Colour fColour = G4Colour::Brown();

};

using TruthHitCollection = G4THitsCollection<TruthHit>;
extern G4ThreadLocal G4Allocator<TruthHit> *truthHitAllocator;

inline void* TruthHit::operator new(size_t) {
  if (!truthHitAllocator) {
    truthHitAllocator = new G4Allocator<TruthHit>;
  }
  return truthHitAllocator->MallocSingle();
}

inline void TruthHit::operator delete(void* aHit) {
  if (!truthHitAllocator) {
    truthHitAllocator = new G4Allocator<TruthHit>;
  }
  truthHitAllocator->FreeSingle((TruthHit*) aHit);
}

#endif