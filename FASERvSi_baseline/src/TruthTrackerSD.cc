#include "TruthTrackerSD.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "TruthHit.hh"


TruthTrackerSD::TruthTrackerSD(G4String name) :
  G4VSensitiveDetector(name) {
  G4cout << "creating a sensitive detector with name: " << name << G4endl;
  collectionName.insert(name);
}

TruthTrackerSD::~TruthTrackerSD(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TruthTrackerSD::Initialize(G4HCofThisEvent *HCE) {
  G4cout << "Initializing TruthTrackerSD" << G4endl;
  fHitCollection = new TruthHitsCollection(GetName(), collectionName[0]);

  auto *runManager = G4RunManager::GetRunManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TruthTrackerSD::EndOfEvent(G4HCofThisEvent *HCE) {

  if (fHCID < 0) 
  { 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }

  HCE->AddHitsCollection(fHCID, fHitCollection);
  std::cout << "TruthSD::EndOfEvent Number of hits in this event: " << fNHits << std::endl;
  fNHits = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool TruthTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist){
  // G4cout << "Processing hit in TruthTrackerSD" << G4endl;
  G4Track* track = aStep->GetTrack();
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
  G4double charge = track->GetDynamicParticle()->GetCharge();

  if (charge == 0) return false; // skip neutral particles, they don't hit

  G4ThreeVector posHit = preStepPoint->GetPosition();
  G4int pdgid = track->GetParticleDefinition()->GetPDGEncoding();
  G4double energy = track->GetDynamicParticle()->Get4Momentum().e();
  G4double px = track->GetDynamicParticle()->Get4Momentum().px();
  G4double py = track->GetDynamicParticle()->Get4Momentum().py();
  G4double pz = track->GetDynamicParticle()->Get4Momentum().pz();
  G4double m = track->GetDynamicParticle()->Get4Momentum().m();
  G4double time = track->GetDynamicParticle()->Get4Momentum().t();
  
  G4VPhysicalVolume* physVol = preStepPoint->GetPhysicalVolume();
  G4int module_number = preStepPoint->GetTouchableHandle()->GetCopyNumber(0);
  G4int layer_number = preStepPoint->GetTouchableHandle()->GetCopyNumber(1);

  G4TouchableHistory* touchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4ThreeVector sensorCenterGlobal = touchable->GetTranslation();
  G4double sensorCentreZ = sensorCenterGlobal.z();

  TruthHit* tmpHit = new TruthHit();

  tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, sensorCentreZ/mm); // in mm
  tmpHit->SetPDGID(pdgid);
  tmpHit->SetEnergy(energy/GeV);
  tmpHit->SetCharge(charge);
  tmpHit->SetPx(px/GeV);
  tmpHit->SetPy(py/GeV);
  tmpHit->SetPz(pz/GeV);
  tmpHit->SetMass(m/GeV);
  tmpHit->SetTrackID(track->GetTrackID());
  tmpHit->SetParentID(track->GetParentID());
  tmpHit->SetModuleNumber(module_number);
  tmpHit->SetLayerNumber(layer_number);
  tmpHit->SetT(time/ns);
  tmpHit->SetTrackVertex(track->GetVertexPosition()/mm);
  tmpHit->SetTrackP4(track->GetDynamicParticle()->Get4Momentum()/GeV);
  if (track->GetParentID() == 0) {
    tmpHit->SetIsPrimaryTrack(1);
    tmpHit->SetIsSecondaryTrack(0);
  }
  else {
    tmpHit->SetIsPrimaryTrack(0);
    tmpHit->SetIsSecondaryTrack(1);
  }

  G4cout << "Truth hit: " << *tmpHit << G4endl;

  fHitCollection->insert(tmpHit);
  fNHits++;
  
  return 0;
}