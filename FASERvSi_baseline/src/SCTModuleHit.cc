#include "SCTModuleHit.hh"

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