#include "SCTModuleDetector.hh"
#include "G4SystemOfUnits.hh"
#include "SCTModuleHit.hh"

SCTModuleDetector::SCTModuleDetector(G4String name) :
  G4VSensitiveDetector(name) {
  G4cout << "creating a sensitive detector with name: " << name << G4endl;
  collectionName.insert(name);
}

SCTModuleDetector::~SCTModuleDetector(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SCTModuleDetector::Initialize(G4HCofThisEvent *HCE) {
  G4cout << "Initializing SCTModuleDetector" << G4endl;
  fHitCollection = new SCTModuleHitsCollection(GetName(), collectionName[0]);

  if (fHCID < 0) { fHCID = GetCollectionID(0); }
  HCE->AddHitsCollection(fHCID, fHitCollection);
  fTrackIDRecord.clear(); // Clear the track ID record for each event
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void SCTModuleDetector::EndOfEvent(G4HCofThisEvent *) {
// }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool SCTModuleDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist){
  // G4cout << "Processing hit in SCTModuleDetector" << G4endl;
  G4Track* track = aStep->GetTrack();
  //track->SetTrackStatus(fStopAndKill);
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  G4String volName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  // G4cout << "Hit volume: " << volName << G4endl;

  G4ThreeVector posHit = preStepPoint->GetPosition();
  G4int pdgid = track->GetParticleDefinition()->GetPDGEncoding();
  G4double energy = track->GetDynamicParticle()->Get4Momentum().e();
  G4double px = track->GetDynamicParticle()->Get4Momentum().px();
  G4double py = track->GetDynamicParticle()->Get4Momentum().py();
  G4double pz = track->GetDynamicParticle()->Get4Momentum().pz();
  G4double m = track->GetDynamicParticle()->Get4Momentum().m();
  G4double charge = track->GetDynamicParticle()->GetCharge();
  G4double time = track->GetDynamicParticle()->Get4Momentum().t();
  G4ThreeVector delta_momentum = aStep->GetDeltaMomentum();
  G4double delta_energy = aStep->GetDeltaEnergy();
  G4int sensor_id = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

  G4int strip_number = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
  G4int strip_side = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
  G4int module_number = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);
  G4int layer_number = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);

  
  G4TouchableHandle touchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector sensorCenterGlobal = touchable->GetTranslation();
  G4double sensorCentreZ = sensorCenterGlobal.z();

  // for (int i = 0; i <= touchable->GetHistoryDepth(); ++i) {
  //   G4String volName = touchable->GetVolume(i)->GetName();
  //   G4int copyNum = touchable->GetCopyNumber(i);
  //   G4cout << "Level " << i << ": " << volName << " (copy " << copyNum << ")" << G4endl;
  // }

  SCTModuleHit* tmpHit = new SCTModuleHit();

  // tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, posHit[2]/mm); // in mm
  // fix the hit z-position to be the centre of the sensor - this way every hit on the same sensor has the same z-pos

  tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, sensorCentreZ/mm); // in mm
  tmpHit->SetPDGID(pdgid);
  tmpHit->SetEnergy(energy/GeV);
  tmpHit->SetCharge(charge);
  tmpHit->SetPx(px/GeV);
  tmpHit->SetPy(py/GeV);
  tmpHit->SetPz(pz/GeV);
  tmpHit->SetDeltaPx(delta_momentum.x()/GeV);
  tmpHit->SetDeltaPy(delta_momentum.y()/GeV);
  tmpHit->SetDeltaPz(delta_momentum.z()/GeV);
  tmpHit->SetDeltaE(delta_energy/GeV);
  tmpHit->SetMass(m/GeV);
  tmpHit->SetTrackID(track->GetTrackID());
  tmpHit->SetParentID(track->GetParentID());
  tmpHit->SetCopyNumSensor(sensor_id);
  tmpHit->SetStripNumber(strip_number);
  tmpHit->SetStripSide(strip_side);
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

  // Check if this sensor has already been hit by this track
  if (fTrackIDRecord.find(sensor_id) != fTrackIDRecord.end()) 
  {
    std::vector<G4int> tracks_that_hit_sensor = fTrackIDRecord[sensor_id];
    if (std::find(tracks_that_hit_sensor.begin(), tracks_that_hit_sensor.end(), track->GetTrackID()) != tracks_that_hit_sensor.end()) 
    {
      return 0; // This track has already hit this sensor, so we skip this hit
    }
    fTrackIDRecord[sensor_id].push_back(track->GetTrackID());
  }
  else
  {
    fTrackIDRecord[sensor_id] = {track->GetTrackID()};
  }

  fHitCollection->insert(tmpHit);
  
  return 0;
}