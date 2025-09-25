#include "SCTModuleHit.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

G4ThreadLocal G4Allocator<SCTModuleHit>* hitAllocator = nullptr;


std::ostream& operator<<(std::ostream& os, const SCTModuleHit& hit) {
    os << "Hit ["
       << "Strip: " << hit.GetStripNumber()
       << ", Side: " << hit.GetStripSide()
       << ", Module: " << hit.GetModuleNumber()
       << ", Layer: " << hit.GetLayerNumber()
       << ", PDG: " << hit.GetPDGID()
       << "]";
    return os;
}

/// @brief  Comparison operator to allow sorting of hits. Means that hits can be stored in ordered containers like std::set
/// @param other 
/// @return 
bool SCTModuleHit::operator<(const SCTModuleHit& other) const {
    if (fLayerNumber != other.fLayerNumber)
        return fLayerNumber < other.fLayerNumber;
    if (fModuleNumber != other.fModuleNumber)
        return fModuleNumber < other.fModuleNumber;
    if (fStripSide != other.fStripSide)
        return fStripSide < other.fStripSide;
    return fStripNumber < other.fStripNumber;
}


void SCTModuleHit::Draw() {
    G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
    if (!visManager) return;
  
    G4ThreeVector pos(fPosX, fPosY, fPosZ);
    G4Circle circle(pos);
    circle.SetScreenSize(5); // pixels
    circle.SetFillStyle(G4Circle::filled);
  
    // Color by truth match status
    G4Colour colour =  fColour;
    G4cout << "Drawing hit at " << fPosX << ", " << fPosY << ", " << fPosZ  << " isReco = " << fIsReco << std::endl;
    G4VisAttributes attribs(colour);
    attribs.SetVisibility(true);
    circle.SetVisAttributes(attribs);
    visManager->Draw(circle);
    // if (fIsReco){
    //     G4cout << "Drawing reco hit at " << fPosX << ", " << fPosY << ", " << fPosZ << std::endl;
    //     visManager->Draw(circle);
    // }
}