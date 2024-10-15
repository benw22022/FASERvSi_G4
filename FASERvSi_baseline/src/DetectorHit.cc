#include "DetectorHit.hh"

#include "G4Transform3D.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"

#include <cstdlib>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorHit::DetectorHit(){}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorHit::Draw() {
  G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    if (!pVVisManager->FilterHit(*this))
      return;
//    G4Transform3D trans(G4RotationMatrix(),
//                        G4ThreeVector(fPosX * CLHEP::cm, fPosY * CLHEP::cm, fPosZ * CLHEP::cm));
//    G4VisAttributes attribs;
//    auto solid = HexagonSolid("dummy", 0.3 * CLHEP::mm, 0.6496345 * CLHEP::cm);
//    G4Colour colour(1, 1, 1);
//    attribs.SetColour(colour);
//    attribs.SetForceSolid(true);
//    pVVisManager->Draw(*solid, attribs, trans);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
