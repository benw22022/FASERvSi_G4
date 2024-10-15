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
  man = G4AnalysisManager::Instance();
  man->OpenFile("output.root");

  //man->SetFirstNtupleColumnId(1);
  
  for (unsigned int i{0}; i < DetectorParameters::Get()->fnumSCTLayers; i++)
  {   

    std::string treeName = "Hits" + std::to_string(i+1);

    man->CreateNtuple(treeName, treeName);
    man->CreateNtupleIColumn("fEvent");
    man->CreateNtupleDColumn("ep_x");
    man->CreateNtupleDColumn("ep_y");
    man->CreateNtupleDColumn("ep_z");
    man->CreateNtupleDColumn("ep_E");
    man->CreateNtupleDColumn("em_x");
    man->CreateNtupleDColumn("em_y");
    man->CreateNtupleDColumn("em_z");
    man->CreateNtupleDColumn("em_E");

    man->CreateNtupleDColumn("mp_x");
    man->CreateNtupleDColumn("mp_y");
    man->CreateNtupleDColumn("mp_z");
    man->CreateNtupleDColumn("mp_E");
    man->CreateNtupleDColumn("mm_x");
    man->CreateNtupleDColumn("mm_y");
    man->CreateNtupleDColumn("mm_z");
    man->CreateNtupleDColumn("mm_E");

    man->CreateNtupleDColumn("hp_x");
    man->CreateNtupleDColumn("hp_y");
    man->CreateNtupleDColumn("hp_z");
    man->CreateNtupleDColumn("hp_E");
    man->CreateNtupleDColumn("hm_x");
    man->CreateNtupleDColumn("hm_y");
    man->CreateNtupleDColumn("hm_z");
    man->CreateNtupleDColumn("hm_E");

    //man->FinishNtuple(0);
    man->FinishNtuple();
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction()
{
  man->Write();
  man->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*)
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* evt)
{
  
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;

  G4int nPerTree=25;
  
    
  auto hce = evt->GetHCofThisEvent();
  auto sdManager = G4SDManager::GetSDMpointer();
  G4int collId;
  G4int offset=0;
  
  for (unsigned int i{0}; i < DetectorParameters::Get()->fnumSCTLayers; i++)
  {   

    std::string detName = "sensDet" + std::to_string(i);  
    
    G4cout << ">>>     Filling ntuple for " << detName << " " << i << G4endl;
    
    
    collId = sdManager->GetCollectionID(detName);
    auto hc = hce->GetHC(collId);
    if (!hc){
      G4cout << ">>>     ERROR didn't find " << detName << G4endl;
      return;
    }
    
    G4int count_charge=0;
    G4int count_e=0;
    G4int count_m=0;
    G4int count_h=0;
    G4int charge_e=0;
    G4int charge_m=0;
    G4int charge_h=0;
    
    
    man->FillNtupleIColumn(i,0,evt->GetEventID());
    
    for (unsigned int i = 0; i < hc->GetSize(); ++i) {
      auto hit = static_cast<DetectorHit *>(hc->GetHit(i));
      //G4cout << "hitX: " << hit->GetX() << G4endl;
      //G4cout << "hit x: " << hit->GetX() << " y: " << hit->GetY() << " z: " << hit->GetZ() << " E: " << hit->GetEnergy() << " ID: " << hit->GetPDGID() << G4endl;     
      if (hit->GetEnergy() < 1)continue;
      
      if (hit->GetPDGID() == -11){
	man->FillNtupleDColumn(i,1,hit->GetX());
	man->FillNtupleDColumn(i,2,hit->GetY());
	man->FillNtupleDColumn(i,3,hit->GetZ());
	man->FillNtupleDColumn(i,4,hit->GetEnergy());
	count_charge-=1;
	count_e+=1;
	charge_e-=1;
      }else if(hit->GetPDGID() == 11){
	man->FillNtupleDColumn(i,5,hit->GetX());
	man->FillNtupleDColumn(i,6,hit->GetY());
	man->FillNtupleDColumn(i,7,hit->GetZ());
	man->FillNtupleDColumn(i,8,hit->GetEnergy());
	count_charge+=1;
	count_e+=1;
	charge_e+=1;
      }else if (hit->GetPDGID() == -13){
	man->FillNtupleDColumn(i,9,hit->GetX());
	man->FillNtupleDColumn(i,10,hit->GetY());
	man->FillNtupleDColumn(i,11,hit->GetZ());
	man->FillNtupleDColumn(i,12,hit->GetEnergy());
	count_charge-=1;
	count_m+=1;
	charge_m-=1;
      }else if(hit->GetPDGID() == 13){
	man->FillNtupleDColumn(i,13,hit->GetX());
	man->FillNtupleDColumn(i,14,hit->GetY());
	man->FillNtupleDColumn(i,15,hit->GetZ());
	man->FillNtupleDColumn(i,16,hit->GetEnergy());
	count_charge+=1;
	count_m+=1;
	charge_m+=1;
      }else if ((hit->GetPDGID() < 0 && hit->GetPDGID() > -6) || hit->GetPDGID() < -25){
	man->FillNtupleDColumn(i,17,hit->GetX());
	man->FillNtupleDColumn(i,18,hit->GetY());
	man->FillNtupleDColumn(i,19,hit->GetZ());
	man->FillNtupleDColumn(i,20,hit->GetEnergy());
	count_charge-=1;
	count_h+=1;
	charge_h-=1;
      }else if ((hit->GetPDGID() > 0 && hit->GetPDGID() < 6) || hit->GetPDGID() > 25){
	man->FillNtupleDColumn(i,21,hit->GetX());
	man->FillNtupleDColumn(i,22,hit->GetY());
	man->FillNtupleDColumn(i,23,hit->GetZ());
	man->FillNtupleDColumn(i,24,hit->GetEnergy());
	count_charge+=1;
	count_h+=1;
	charge_h+=1;
      }
    }
    
    if ( (count_e==2 && charge_e==0) ||
	 (count_m==2 && charge_m==0) ||
	 (count_h==2 && charge_h==0) ){
       man->AddNtupleRow(i);
    }else{
      //G4cout << "WARNING: Did not find pair of opposite charge particles:"  << G4endl;
      //G4cout << "   count_e = " << count_e << " charge_e = " << charge_e << G4endl;
      //G4cout << "   count_m = " << count_m << " charge_m = " << charge_m << G4endl;
      //G4cout << "   count_h = " << count_h << " charge_h = " << charge_h << G4endl; 
    }
    
  }


  
  /*
  if( ftrackerCollID<0 || fcalorimeterCollID<0 || fmuonCollID<0) return;

  G4HCofThisEvent* HCE = evt-> GetHCofThisEvent();
  TrackerHitsCollection* THC = NULL;
  CalorimeterHitsCollection* CHC = NULL;
  MuonHitsCollection* MHC = NULL;

  if( HCE ) {
    THC = (TrackerHitsCollection*)(HCE->GetHC(ftrackerCollID));
    CHC = (CalorimeterHitsCollection*)(HCE->GetHC(fcalorimeterCollID));
    MHC = (MuonHitsCollection*)(HCE->GetHC(fmuonCollID));
  }

  if( THC ) {
    G4int n_hit = THC-> entries();
    G4cout << "     " << n_hit
         << " hits are stored in TrackerHitsCollection." << G4endl;
  }

  if( CHC ) {
    G4int n_hit = CHC-> entries();
    G4cout << "     " << n_hit
         << " hits are stored in CalorimeterHitsCollection." << G4endl;
    G4double totE = 0;
    for( int i = 0; i < n_hit; i++ ) {
      totE += (*CHC)[i]-> GetEdep();
    }
    G4cout << "     Total energy deposition in calorimeter : "
         << totE / GeV << " (GeV)" << G4endl;
  }

  if( MHC ) {
    G4int n_hit = MHC-> entries();
    G4cout << "     " << n_hit
         << " hits are stored in MuonHitsCollection." << G4endl;
  }
  */
}

