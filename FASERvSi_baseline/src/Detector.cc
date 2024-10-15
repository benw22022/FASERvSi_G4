#include "Detector.hh"
#include "G4SystemOfUnits.hh"

Detector::Detector(G4String name) :
  G4VSensitiveDetector(name) {
  G4cout << "creating a sensitive detector with name: " << name << G4endl;
  collectionName.insert(name);
}

Detector::~Detector(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Detector::Initialize(G4HCofThisEvent *HCE) {
  fHitCollection = new DetectorHitCollection(GetName(), collectionName[0]);

  if (fHCID < 0)
    fHCID = GetCollectionID(0);
  HCE->AddHitsCollection(fHCID, fHitCollection);

  fTmpHits.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Detector::EndOfEvent(G4HCofThisEvent *) {
  for (auto it = fTmpHits.begin(); it != fTmpHits.end(); ++it)
    fHitCollection->insert(*it);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool Detector::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist){

  G4Track* track = aStep->GetTrack();
  //track->SetTrackStatus(fStopAndKill);
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  G4ThreeVector posHit = preStepPoint->GetPosition();
  G4int pdgid = track->GetParticleDefinition()->GetPDGEncoding();
  G4double energy = track->GetDynamicParticle()->Get4Momentum().e();
  //G4cout << "track pos: " << posHit << ", pdgid: " << pdgid << G4endl;

  DetectorHit* tmpHit = new DetectorHit();

  tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, posHit[2]/mm); // in mm
  tmpHit->SetPDGID(pdgid);
  tmpHit->SetEnergy(energy/GeV);
  
  fTmpHits.push_back(tmpHit);

  return 0;
}
