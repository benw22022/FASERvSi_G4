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
/// \file eventgenerator/HepMC/HepMCEx01/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "DetectorParameters.hh"
#include <unordered_set>
#include <numeric>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction()
 : G4UserRunAction()
{
  man = G4AnalysisManager::Instance();
  messenger = new RunActionMessenger(this);
  foutputFileName = "output.root";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun-> GetRunID() << " start." << G4endl;
  
  G4RunManager::GetRunManager()-> SetRandomNumberStore(true);

  man->OpenFile(foutputFileName);

  man->CreateNtuple("Hits", "Hits");
  man->CreateNtupleIColumn("fEvent");           // 0
  man->CreateNtupleDColumn("vertex_x");         // 1
  man->CreateNtupleDColumn("vertex_y");         // 2
  man->CreateNtupleDColumn("vertex_z");         // 3
  man->CreateNtupleDColumn("nu_E");             // 4
  man->CreateNtupleDColumn("nu_px");            // 5
  man->CreateNtupleDColumn("nu_py");            // 6
  man->CreateNtupleDColumn("nu_pz");            // 7
  man->CreateNtupleIColumn("nu_pdgc");          // 8
  man->CreateNtupleIColumn("target_pdgc");      // 9
  man->CreateNtupleIColumn("isCC");             // 10
  man->CreateNtupleIColumn("cclepton_pdgc");    // 11
  man->CreateNtupleDColumn("cclepton_E");       // 12
  man->CreateNtupleDColumn("cclepton_px");      // 13
  man->CreateNtupleDColumn("cclepton_py");      // 14
  man->CreateNtupleDColumn("cclepton_pz");      // 15
  man->CreateNtupleDColumn("x",    m_hits_x);
  man->CreateNtupleDColumn("y",    m_hits_y);
  man->CreateNtupleDColumn("z",    m_hits_z);
  man->CreateNtupleDColumn("E",    m_hits_E);
  man->CreateNtupleIColumn("pdgc", m_hits_pdgc);
  man->CreateNtupleDColumn("charge", m_hits_charge);
  man->CreateNtupleIColumn("layer", m_hits_layernum);
  man->FinishNtuple();  
}

void RunAction::EndOfRunAction(const G4Run*) 
{
  man->Write();
  man->CloseFile();
}


G4String RunAction::GetOutputFileName() const
{
  return foutputFileName;
}


void RunAction::SetOutputFileName(G4String fname)
{
  foutputFileName = fname;
}


void RunAction::FillHitsRow(G4double x, G4double y, G4double z, G4double E, G4int pdgc, G4double charge, G4int layernum)
{
  m_hits_x.push_back(x);
  m_hits_y.push_back(y);
  m_hits_z.push_back(z);
  m_hits_E.push_back(E);
  m_hits_pdgc.push_back(pdgc);
  m_hits_charge.push_back(charge);
  m_hits_layernum.push_back(layernum);
}


void RunAction::ClearHits()
{
  m_hits_x.clear();
  m_hits_y.clear();
  m_hits_z.clear();
  m_hits_E.clear();
  m_hits_pdgc.clear();
  m_hits_charge.clear();
  m_hits_layernum.clear();
}


std::set<G4int> RunAction::FindHitsToMerge(G4double xtol, G4double ytol) const
{
  std::set<G4int> hits_to_be_merged = std::set<G4int>();

  if (m_hits_x.size() < 2) return hits_to_be_merged; 

  for (unsigned int i{0}; i < m_hits_x.size(); i++)
  {
    bool can_break{false};
    
    if (m_hits_charge.at(i) == 0){ continue; } // Hits must be charged for merging 

    for (unsigned int j{i+1}; j < m_hits_x.size(); j++)
    {
  
      if (m_hits_charge.at(j) == 0){ continue; } // Hits must be charged for merging 
      
      G4double x_dist = abs(m_hits_x.at(i) - m_hits_x.at(j));
      G4double y_dist = abs(m_hits_y.at(i) - m_hits_y.at(j));

      if (x_dist < xtol and y_dist < ytol and m_hits_z.at(i) == m_hits_z.at(j))
      {
        hits_to_be_merged.emplace(i);
        hits_to_be_merged.emplace(j);
        can_break = true;
      }
    }
    if (can_break) return hits_to_be_merged; 
  }
  return hits_to_be_merged;
}

void RemoveElementsFromVec(const std::set<G4int>& hits_to_be_merged, std::vector<G4double>& a) {
    // Copy elements to a vector and sort them in reverse order
    std::vector<G4int> hits_to_be_merged_vec(hits_to_be_merged.begin(), hits_to_be_merged.end());
    std::sort(hits_to_be_merged_vec.rbegin(), hits_to_be_merged_vec.rend());

    // Iterate and remove elements
    for (const auto& index : hits_to_be_merged_vec) {
        if (index < a.size()) {
            a.erase(a.begin() + index);
        }
    }
}


void RemoveElementsFromVec(const std::set<G4int>& hits_to_be_merged, std::vector<G4int>& a) {
    // Copy elements to a vector and sort them in reverse order
    std::vector<G4int> hits_to_be_merged_vec(hits_to_be_merged.begin(), hits_to_be_merged.end());
    std::sort(hits_to_be_merged_vec.rbegin(), hits_to_be_merged_vec.rend());

    // Iterate and remove elements
    for (const auto& index : hits_to_be_merged_vec) {
        if (index < a.size()) {
            a.erase(a.begin() + index);
        }
    }
}


G4int RunAction::MergeHits(G4double xtol, G4double ytol)
{   
  
  G4int number_of_merged_hits = 0;

  while (true)
  {
    std::set<G4int> hits_to_be_merged = this->FindHitsToMerge(xtol, ytol);
      
    if (hits_to_be_merged.size() == 0) { break; }
    
    std::vector<G4double> xpos;
    std::vector<G4double> ypos;
    std::vector<G4double> zpos;
    std::vector<G4double> charge;
    std::vector<G4double> energy;
    std::vector<G4int> pdgc;
    std::vector<G4int> layer;
    
    for (const int& elem : hits_to_be_merged) {
      xpos.push_back(m_hits_x.at(elem));
      ypos.push_back(m_hits_y.at(elem));
      zpos.push_back(m_hits_z.at(elem));
      charge.push_back(m_hits_charge.at(elem));
      pdgc.push_back(m_hits_pdgc.at(elem));
      energy.push_back(m_hits_E.at(elem));
      layer.push_back(m_hits_layernum.at(elem));
      number_of_merged_hits++;
    }

    G4double merged_hit_x      = std::reduce(xpos.begin(), xpos.end()) / static_cast<float>(xpos.size());
    G4double merged_hit_y      = std::reduce(ypos.begin(), ypos.end()) / static_cast<float>(ypos.size());
    G4double merged_hit_z      = std::reduce(zpos.begin(), zpos.end()) / static_cast<float>(zpos.size());
    G4double merged_hit_charge = std::reduce(charge.begin(), charge.end());
    G4double merged_hit_energy = std::reduce(energy.begin(), energy.end());
    G4double merged_hit_layer  = layer.at(0);

    // Construct new PDGC for merged hits
    G4int merged_hit_pdgc{22222222};
    for (const int& p : pdgc) 
    {
      if (abs(p) == 11 || abs(p) == 13 || abs(p) == 15)
      { 
        merged_hit_pdgc = 11111111;
        break;
      }
    }

    RemoveElementsFromVec(hits_to_be_merged, m_hits_x);
    RemoveElementsFromVec(hits_to_be_merged, m_hits_y);
    RemoveElementsFromVec(hits_to_be_merged, m_hits_z);
    RemoveElementsFromVec(hits_to_be_merged, m_hits_E);
    RemoveElementsFromVec(hits_to_be_merged, m_hits_pdgc);
    RemoveElementsFromVec(hits_to_be_merged, m_hits_charge);
    RemoveElementsFromVec(hits_to_be_merged, m_hits_layernum);
    
    m_hits_x.push_back(merged_hit_x);
    m_hits_y.push_back(merged_hit_y);
    m_hits_z.push_back(merged_hit_z);
    m_hits_E.push_back(merged_hit_energy);
    m_hits_pdgc.push_back(merged_hit_pdgc);
    m_hits_charge.push_back(merged_hit_charge);
    m_hits_layernum.push_back(merged_hit_layer);
  }

  return number_of_merged_hits;
}