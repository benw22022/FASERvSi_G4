#include "TruthHit.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

G4ThreadLocal G4Allocator<TruthHit>* truthHitAllocator = nullptr;


std::ostream& operator<<(std::ostream& os, const TruthHit& hit) {
    os << "Truth Hit ["
       << "pos : " << hit.GetX() << ", " << hit.GetY() << ", " << hit.GetZ()
       << ", Module: " << hit.GetModuleNumber()
       << ", Layer: " << hit.GetLayerNumber()
       << ", PDG: " << hit.GetPDGID()
       << "]";
    return os;
}

/// @brief  Comparison operator to allow sorting of hits. Means that hits can be stored in ordered containers like std::set
/// @param other 
/// @return 
bool TruthHit::operator<(const TruthHit& other) const {
    if (fPosX != other.fPosX)
        return fPosX < other.fPosX;
    if (fPosY != other.fPosY)
        return fPosY < other.fPosY;
    if (fPosZ != other.fPosZ)   
        return fPosZ < other.fPosZ;
    else return false;
}


void TruthHit::Draw() {
    G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
  
    G4ThreeVector pos(fPosX, fPosY, fPosZ);
    G4Circle circle(pos);
    circle.SetScreenSize(3); // pixels
    circle.SetFillStyle(G4Circle::filled);

    // G4cout << "Drawing truth hit at " << fPosX << ", " << fPosY << ", " << fPosZ << std::endl;
    G4VisAttributes attribs(fColour);
    attribs.SetVisibility(true);
    circle.SetVisAttributes(attribs);
    // visManager->Draw(circle);
}