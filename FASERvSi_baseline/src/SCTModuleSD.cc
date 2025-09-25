#include "SCTModuleDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "SCTModuleGeometry.hh"
#include "SCTModuleHit.hh"
#include "reco/Channel.hh"



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

  auto *runManager = G4RunManager::GetRunManager();
  fDetector = (DetectorConstruction*) (runManager->GetUserDetectorConstruction());


  // if (fHCID < 0) { fHCID = GetCollectionID(0); }
  // HCE->AddHitsCollection(fHCID, fHitCollection);
  // fTrackIDRecord.clear(); // Clear the track ID record for each event
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCTModuleDetector::EndOfEvent(G4HCofThisEvent *HCE) {

  if (fHCID < 0) 
  { 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }

  HCE->AddHitsCollection(fHCID, fHitCollection);
  fTrackIDRecord.clear(); // Clear the track ID record for each event
  std::cout << "SCTModuleSD::EndOfEvent Number of hits in this event: " << fNHits << std::endl;
  fNHits = 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool SCTModuleDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist){
  // G4cout << "Processing hit in SCTModuleDetector" << G4endl;
  G4Track* track = aStep->GetTrack();
  //track->SetTrackStatus(fStopAndKill);
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
  G4double charge = track->GetDynamicParticle()->GetCharge();

  if (charge == 0) return false; // skip neutral particles, they don't hit

  G4String volName = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
  // G4cout << "Hit volume: " << volName << G4endl;

  G4ThreeVector posHit = preStepPoint->GetPosition();
  G4int pdgid = track->GetParticleDefinition()->GetPDGEncoding();
  G4double energy = track->GetDynamicParticle()->Get4Momentum().e();
  G4double px = track->GetDynamicParticle()->Get4Momentum().px();
  G4double py = track->GetDynamicParticle()->Get4Momentum().py();
  G4double pz = track->GetDynamicParticle()->Get4Momentum().pz();
  G4double m = track->GetDynamicParticle()->Get4Momentum().m();
  // G4double charge = track->GetDynamicParticle()->GetCharge();
  G4double time = track->GetDynamicParticle()->Get4Momentum().t();
  G4ThreeVector delta_momentum = aStep->GetDeltaMomentum();
  G4double delta_energy = aStep->GetDeltaEnergy();
  G4VPhysicalVolume* physVol = preStepPoint->GetPhysicalVolume();


  G4int strip_number = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
  G4int strip_side = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
  G4int module_number = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);
  G4int layer_number = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);

  Channel channel;
  channel.setStrip(strip_number);
  channel.setSide(strip_side);
  channel.setModule(module_number);
  channel.setLayer(layer_number);

  // G4TouchableHandle touchable = preStepPoint->GetTouchableHandle();
  G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4ThreeVector sensorCenterGlobal = touchable->GetTranslation();
  G4double sensorCentreZ = sensorCenterGlobal.z();

  // G4cout << G4endl;
  // for (int i = 0; i <= touchable->GetHistoryDepth(); ++i) {
  //   G4String volName = touchable->GetVolume(i)->GetName();
  //   G4int copyNum = touchable->GetCopyNumber(i);
  //   G4cout << "Level " << i << ": " << volName << " (copy " << copyNum << ")" << G4endl;
  // }

  SCTModuleHit* tmpHit = new SCTModuleHit();

  // tmpHit->SetPosition(posHit[0]/mm, posHit[1]/mm, posHit[2]/mm); // in mm
  // fix the hit z-position to be the centre of the sensor - this way every hit on the same sensor has the same z-pos
  G4ThreeVector strip_global_pos = aStep->GetPreStepPoint()->GetTouchable()->GetTranslation();
  const G4RotationMatrix* strip_global_rotation = aStep->GetPreStepPoint()->GetTouchable()->GetRotation();


  for (G4int i = 0; i < touchable->GetHistory()->GetDepth(); ++i) {
    G4VPhysicalVolume* pv = touchable->GetVolume(i);
    G4ThreeVector pos = pv->GetTranslation();
    G4RotationMatrix* rot = pv->GetRotation();
    G4cout << "Level " << i << ": " << pv->GetName() << " at " << pos << G4endl;
}

  //* Get the geometry of the strip that was hit
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4TouchableHistory* touchable1 = (G4TouchableHistory*)(preStep->GetTouchable());
  G4VPhysicalVolume* stripPhys = touchable1->GetVolume();

  G4ThreeVector stripTranslation = stripPhys->GetTranslation();
  const G4RotationMatrix* stripRotation = stripPhys->GetRotation();

  //* Transform the local coordinates of the strip ends to the global coordinate system
  G4ThreeVector localStart(0, -SCTModuleGeometry::stripLength()/2, 0);
  G4ThreeVector localEnd(0, SCTModuleGeometry::stripLength()/2, 0);
  G4ThreeVector globalStart = touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(localStart);
  G4ThreeVector globalEnd   = touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(localEnd);
  
  tmpHit->SetStripEnds(std::make_pair(globalStart, globalEnd));
  tmpHit->SetStripCentre(strip_global_pos);
  tmpHit->SetStripRotation(*strip_global_rotation);
  
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
  tmpHit->SetStripNumber(strip_number);
  tmpHit->SetStripSide(strip_side);
  tmpHit->SetModuleNumber(module_number);
  tmpHit->SetLayerNumber(layer_number);
  tmpHit->SetT(time/ns);
  tmpHit->SetTrackVertex(track->GetVertexPosition()/mm);
  tmpHit->SetTrackP4(track->GetDynamicParticle()->Get4Momentum()/GeV);
  tmpHit->SetPhysVol(physVol);
  if (track->GetParentID() == 0) {
    tmpHit->SetIsPrimaryTrack(1);
    tmpHit->SetIsSecondaryTrack(0);
  }
  else {
    tmpHit->SetIsPrimaryTrack(0);
    tmpHit->SetIsSecondaryTrack(1);
  }

  // // Check if this sensor has already been hit by this track
  // G4int sensor_id = channel.value(); // Unique identifier for the sensor based on strip, side, module, layer
  // if (fTrackIDRecord.find(sensor_id) != fTrackIDRecord.end()) 
  // {
  //   std::vector<G4int> tracks_that_hit_sensor = fTrackIDRecord[sensor_id];
  //   if (std::find(tracks_that_hit_sensor.begin(), tracks_that_hit_sensor.end(), track->GetTrackID()) != tracks_that_hit_sensor.end()) 
  //   {
  //     return 0; // This track has already hit this sensor, so we skip this hit
  //   }
  //   fTrackIDRecord[sensor_id].push_back(track->GetTrackID());
  // }
  // else
  // {
  //   fTrackIDRecord[sensor_id] = {track->GetTrackID()};
  // }

  fHitCollection->insert(tmpHit);
  fNHits++;
  
  return 0;
}