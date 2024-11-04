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
/// \file eventgenerator/HepMC/HepMCEx01/include/EventInformation.hh
/// \brief Definition of the EventInformation class
//
//

#ifndef EventInformation_h
#define EventInformation_h 1

#include "G4VUserEventInformation.hh"
#include "globals.hh"
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "HepMCG4Interface.hh"

#include "G4RunManager.hh"
#include "G4LorentzVector.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "EventInformation.hh"



class EventInformation : public G4VUserEventInformation {

  public:
    EventInformation();
    EventInformation(G4int nu_pdgid, G4int target_pdgid, G4LorentzVector nu_p4, G4LorentzVector vpos);
    ~EventInformation() override;
    
    G4int GetNeutrinoPDG() const;
    G4int GetTargetPDG() const;
    G4LorentzVector GetNeutrinoP4() const;
    G4LorentzVector GetVertexPos() const;
    
    void SetNeutrinoPDG(G4int pdgc);
    void SetTargetPDG(G4int pdgc);
    void SetNeutrinoP4(G4LorentzVector p4);
    void SetVertexPos(G4LorentzVector  vpos);
    
    static bool isPDGNeutrino(G4int pdgc);

    G4double GetNeutrinoPx() const;
    G4double GetNeutrinoPy() const;
    G4double GetNeutrinoPz() const;
    G4double GetNeutrinoE() const;

    G4double GetVertexX() const;
    G4double GetVertexY() const;
    G4double GetVertexZ() const;

    bool isCCInteraction() const;
    void SetIsCCInteraction(bool isCC);

    void SetCCLeptonPDG(G4int pdgc);
    G4int GetCCLeptonPDG();
    void SetCCLeptonP4(G4LorentzVector p4);
    G4LorentzVector GetCCLeptonP4();

    G4double GetCCLeptonPx() const;
    G4double GetCCLeptonPy() const;
    G4double GetCCLeptonPz() const;
    G4double GetCCLeptonE() const;

    void Print() const override;

  private:
    G4int fnu_pdgc;
    G4int ftarget_pdgc;
    G4LorentzVector fvertex_pos;
    G4LorentzVector fnu_p4;
    G4LorentzVector fCCLepton_p4{0,0,0,0};
    G4int fCCLepton_pdgc{-999};
    bool fisCCInteraction{0};

};

#endif
