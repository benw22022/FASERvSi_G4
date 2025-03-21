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
// * institutes,nor the evtActcies providing financial support for this *
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
/// \file Rungenerator/HepMC/HepMCEx01/src/RunActionMessenger.cc
/// \brief Implementation of the RunActionMessenger class
//
//
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "RunActionMessenger.hh"
#include "RunAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMessenger::RunActionMessenger
                             (RunAction* evtAct)
  : act(evtAct)
{
  dir= new G4UIdirectory("/ntuple/");
  dir-> SetGuidance("Ntuple writing");

//   verbose=
//     new G4UIcmdWithAnInteger("/ntuple/verbose", this);
//   verbose-> SetGuidance("Set verbose level");
//   verbose-> SetParameterName("verboseLevel", false, false);
//   verbose-> SetRange("verboseLevel>=0 && verboseLevel<=1");

  open= new G4UIcmdWithAString("/ntuple/output", this);
  open-> SetGuidance("path to output NTuple");
  open-> SetParameterName("Ntuple output name", true, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunActionMessenger::~RunActionMessenger()
{
  // delete verbose;
  delete open;

  delete dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
//   if (command==verbose) {
//     int level= verbose-> GetNewIntValue(newValues);
//     gen-> SetVerboseLevel(level);
//   } else 
  if (command==open) {
    act-> SetOutputFileName(newValues);
    G4cout << "Writing to NTuple: " << act-> GetOutputFileName() << G4endl;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String RunActionMessenger::GetCurrentValue(G4UIcommand* command)
{
    G4String cv;

//   if (command == verbose) {
//     cv = verbose-> ConvertToString(gen-> GetVerboseLevel());
//   } else  
    if (command == open) 
    {
        cv = act-> GetOutputFileName();
    }
    return cv;
}
