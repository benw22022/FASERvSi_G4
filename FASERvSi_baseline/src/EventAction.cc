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
  
  for (unsigned int tree_idx{0}; tree_idx < DetectorParameters::Get()->fnumSCTLayers; tree_idx++)
  {   

    std::string detName = "sensDet" + std::to_string(tree_idx);  
    // man->FillNtupleIColumn(tree_idx, 0, evt->GetEventID());
    
    G4cout << ">>>     Filling ntuple for " << detName << ">>> Event " << evt->GetEventID() << G4endl;
    
    
    collId = sdManager->GetCollectionID(detName);
    auto hc = hce->GetHC(collId);
    if (!hc){
      G4cout << ">>>     ERROR didn't find " << detName << G4endl;
      return;
    }

    man->FillNtupleIColumn(tree_idx,0, evt->GetEventID());
    
    if (hc->GetSize() == 0)
    {
      G4cout << ">>>     WARNING: Empty Hits Collection in " << detName << G4endl;
      // man->FillNtupleIColumn(tree_idx, 0, evt->GetEventID());
      man->FillNtupleDColumn(tree_idx,1,-999);
      man->FillNtupleDColumn(tree_idx,2,-999);
      man->FillNtupleDColumn(tree_idx,3,-999);
      man->FillNtupleDColumn(tree_idx,4,-999);
      man->FillNtupleDColumn(tree_idx,5,-999);
      man->AddNtupleRow(tree_idx);
    }

    for (unsigned int i = 0; i < hc->GetSize(); ++i) 
    {
      auto hit = static_cast<DetectorHit *>(hc->GetHit(i));
      if (hit->GetEnergy() < 1)continue;
      man->FillNtupleDColumn(tree_idx,1,hit->GetX());
      man->FillNtupleDColumn(tree_idx,2,hit->GetY());
      man->FillNtupleDColumn(tree_idx,3,hit->GetZ());
      man->FillNtupleDColumn(tree_idx,4,hit->GetEnergy());
      man->FillNtupleDColumn(tree_idx,5,hit->GetPDGID());
      man->AddNtupleRow(tree_idx);
    }
  }
}
