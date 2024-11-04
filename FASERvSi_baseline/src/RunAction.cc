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

  man->CreateNtuple("truth", "truth");
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
  man->FinishNtuple();

  for (unsigned int i{0}; i < DetectorParameters::Get()->fnumSCTLayers; i++)
  {   
    std::string treeName = "Hits" + std::to_string(i+1);
    std::cout << "Making tree " << treeName << std::endl;

    man->CreateNtuple(treeName, treeName);
    man->CreateNtupleIColumn("fEvent");         // 0           
    man->CreateNtupleDColumn("x");              // 1      
    man->CreateNtupleDColumn("y");              // 2      
    man->CreateNtupleDColumn("z");              // 3      
    man->CreateNtupleDColumn("E");              // 4      
    man->CreateNtupleDColumn("pdgc");           // 5         
    man->CreateNtupleDColumn("prestep_x");      // 6        
    man->CreateNtupleDColumn("prestep_y");      // 7         
    man->CreateNtupleDColumn("prestep_z");      // 8     
    man->CreateNtupleDColumn("poststep_x");     // 9         
    man->CreateNtupleDColumn("poststep_y");     // 10        
    man->CreateNtupleDColumn("poststep_z");     // 11            
    man->FinishNtuple();
  }  
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