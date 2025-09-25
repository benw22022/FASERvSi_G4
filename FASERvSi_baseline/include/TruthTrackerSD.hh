#ifndef TruthTrackerSD_HH
#define TruthTrackerSD_HH

#include "TruthHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include  "G4LorentzVector.hh"
#include "DetectorConstruction.hh"

typedef G4THitsCollection<TruthHit> TruthHitsCollection;

class TruthTrackerSD : public G4VSensitiveDetector
{
public:
    TruthTrackerSD(G4String);
    ~TruthTrackerSD();

  void Initialize(G4HCofThisEvent *HCE);

  void EndOfEvent(G4HCofThisEvent *HCE);

  G4bool ProcessHits(G4Step*, G4TouchableHistory*);

private:
  /// Hit collection stored in the event, filled in at the end of event based
  /// on temporary hits
  TruthHitsCollection *fHitCollection{nullptr};
  /// ID of hit collection
  G4int fHCID = -1;
  G4int fNHits = 0;
};

#endif