//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file eventgenerator/HepMC/HepMCEx01/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
//

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VHitsCollection.hh"
#include "EventAction.hh"
#include "DetectorHit.hh"
#include "MuonHit.hh"
#include "TrackerHit.hh"

#include "DetectorParameters.hh"
#include <string>
#include "Randomize.hh"
#include "EventInformation.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::EventAction()
 : G4UserEventAction(),
   ftrackerCollID(-1),
   fcalorimeterCollID(-1),
   fmuonCollID(-1)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*)
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  

  G4int nPerTree=25;
  
  auto hce = evt->GetHCofThisEvent();
  auto sdManager = G4SDManager::GetSDMpointer();
  G4int collId;
  G4int offset=0;

  // Choose a random number to be the event number
  G4int eventNumber = static_cast<G4int>(G4UniformRand() * 90000000) + 10000000;

  // Fill truth tree
  EventInformation* eventInfo = static_cast<EventInformation*>(evt->GetUserInformation());

  man->FillNtupleIColumn(0, 0,  eventNumber);
  man->FillNtupleDColumn(0, 1,  eventInfo->GetVertexX());
  man->FillNtupleDColumn(0, 2,  eventInfo->GetVertexY());
  man->FillNtupleDColumn(0, 3,  eventInfo->GetVertexZ());
  man->FillNtupleDColumn(0, 4,  eventInfo->GetNeutrinoE());
  man->FillNtupleDColumn(0, 5,  eventInfo->GetNeutrinoPx());
  man->FillNtupleDColumn(0, 6,  eventInfo->GetNeutrinoPy());
  man->FillNtupleDColumn(0, 7,  eventInfo->GetNeutrinoPz());
  man->FillNtupleIColumn(0, 8,  eventInfo->GetNeutrinoPDG());
  man->FillNtupleIColumn(0, 9,  eventInfo->GetTargetPDG());
  man->FillNtupleIColumn(0, 10, eventInfo->isCCInteraction());
  man->FillNtupleIColumn(0, 11, eventInfo->GetCCLeptonPDG());
  man->FillNtupleDColumn(0, 12, eventInfo->GetCCLeptonE());
  man->FillNtupleDColumn(0, 13, eventInfo->GetCCLeptonPx());
  man->FillNtupleDColumn(0, 14, eventInfo->GetCCLeptonPy());
  man->FillNtupleDColumn(0, 15, eventInfo->GetCCLeptonPz());
  man->AddNtupleRow(0);

  G4cout << "EventInformation: " << \
  eventNumber << ", " <<\
  "vx = " << eventInfo->GetVertexX() << " mm, " << \
  "vy = " << eventInfo->GetVertexY() << " mm, " << \
  "vz = " << eventInfo->GetVertexZ() << " mm, " << \
  "nuE = " << eventInfo->GetNeutrinoE() << " GeV, " << \
  "nuPx = " << eventInfo->GetNeutrinoPx() << " GeV, " << \
  "nuPy = " << eventInfo->GetNeutrinoPy() << " GeV, " <<\
  "nuPz = " << eventInfo->GetNeutrinoPz() << " GeV, " <<\
  "nu pdgc" << eventInfo->GetNeutrinoPDG() << ", " <<\
  "target pdgc" << eventInfo->GetTargetPDG() << G4endl;
  
  // Fill Hits trees
  for (unsigned int det_idx{0}; det_idx < DetectorParameters::Get()->fnumSCTLayers; det_idx++)
  {   
    G4int tree_idx = det_idx + 1; // Need to offset by one since tree `0` is the `truth` tree
    std::string detName = "sensDet" + std::to_string(det_idx);  
    G4cout << ">>>     Filling ntuple for " << detName << " >>> Event " << evt->GetEventID() << "  EventID = " << eventNumber << G4endl;    

    
    collId = sdManager->GetCollectionID(detName);
    auto hc = hce->GetHC(collId);
    if (!hc){
      G4cout << ">>>     ERROR didn't find " << detName << G4endl;
      return;
    }

    if (hc->GetSize() == 0)
    {
      G4cout << ">>>     WARNING: Empty Hits Collection in " << detName << G4endl;
      man->FillNtupleIColumn(tree_idx,0, eventNumber);
      man->FillNtupleDColumn(tree_idx,1,-999);
      man->FillNtupleDColumn(tree_idx,2,-999);
      man->FillNtupleDColumn(tree_idx,3,-999);
      man->FillNtupleDColumn(tree_idx,4,-999);
      man->FillNtupleDColumn(tree_idx,5,-999);
      man->FillNtupleDColumn(tree_idx,6, -999);
      man->FillNtupleDColumn(tree_idx,7, -999);
      man->FillNtupleDColumn(tree_idx,8, -999);
      man->FillNtupleDColumn(tree_idx,9, -999);
      man->FillNtupleDColumn(tree_idx,10,-999);
      man->FillNtupleDColumn(tree_idx,11,-999);
      man->AddNtupleRow(tree_idx);
    }

    for (unsigned int i = 0; i < hc->GetSize(); ++i) 
    {
      auto hit = static_cast<DetectorHit *>(hc->GetHit(i));

      double pre_x  = hit->GetPreStepPosition().x();
      double pre_y  = hit->GetPreStepPosition().y();
      double pre_z  = hit->GetPreStepPosition().z();
      double post_x = hit->GetPostStepPosition().x();
      double post_y = hit->GetPostStepPosition().y();
      double post_z = hit->GetPostStepPosition().z();

      if (hit->GetEnergy() < 1)continue;
      man->FillNtupleIColumn(tree_idx,0, eventNumber);
      man->FillNtupleDColumn(tree_idx,1,hit->GetX());
      man->FillNtupleDColumn(tree_idx,2,hit->GetY());
      man->FillNtupleDColumn(tree_idx,3,hit->GetZ());
      man->FillNtupleDColumn(tree_idx,4,hit->GetEnergy());
      man->FillNtupleDColumn(tree_idx,5,hit->GetPDGID());
      man->FillNtupleDColumn(tree_idx,6, hit->GetPreStepPosition().x());
      man->FillNtupleDColumn(tree_idx,7, hit->GetPreStepPosition().y());
      man->FillNtupleDColumn(tree_idx,8, hit->GetPreStepPosition().z());
      man->FillNtupleDColumn(tree_idx,9, hit->GetPostStepPosition().x());
      man->FillNtupleDColumn(tree_idx,10, hit->GetPostStepPosition().y());
      man->FillNtupleDColumn(tree_idx,11, hit->GetPostStepPosition().z());
      man->AddNtupleRow(tree_idx);
    }
  }
}
