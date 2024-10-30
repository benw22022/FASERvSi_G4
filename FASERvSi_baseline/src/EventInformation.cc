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

#include <stdexcept>
#include "EventInformation.hh"
#include "G4SystemOfUnits.hh"


EventInformation::EventInformation()
{
    
}

EventInformation::EventInformation(G4int nu_pdgc, G4int target_pdgc, G4LorentzVector nu_p4, G4LorentzVector vpos)
: ftarget_pdgc{target_pdgc}, fnu_p4{nu_p4}, fvertex_pos{vpos}
{
    SetNeutrinoPDG(nu_pdgc);
} 

EventInformation::~EventInformation()
{
}


G4int EventInformation::GetNeutrinoPDG() const
{
    return fnu_pdgc;
}


G4int EventInformation::GetTargetPDG() const
{
    return ftarget_pdgc;
}

G4LorentzVector EventInformation::GetNeutrinoP4() const
{
    return fnu_p4;
}

G4LorentzVector EventInformation::GetVertexPos() const
{
    return fvertex_pos;
}

void EventInformation::SetNeutrinoPDG(G4int pdgc)
{   
    if (isPDGNeutrino(pdgc))
    {
        fnu_pdgc = pdgc;
    }
    else
    {
        throw std::invalid_argument("Tried to set EventInformation with an invalid neutrino PDG code. Code must be (-) 12, 14 or 16. Code supplied was " + std::to_string(pdgc));
    }

}

void EventInformation::SetTargetPDG(G4int pdgc)
{
    ftarget_pdgc = pdgc;
}

void EventInformation::SetNeutrinoP4(G4LorentzVector p4)
{
    fnu_p4 = p4;
}


void EventInformation::SetVertexPos(G4LorentzVector vpos)
{
    fvertex_pos = vpos;
}


bool EventInformation::isPDGNeutrino(G4int pdgc)
{
    if (abs(pdgc) == 12 || abs(pdgc) == 14 || abs(pdgc) == 16)
    {
        return true;
    }
    return false;
}


G4double EventInformation::GetNeutrinoPx() const
{
    return fnu_p4.px()*GeV;
}

G4double EventInformation::GetNeutrinoPy() const
{
    return fnu_p4.py()*GeV;
}

G4double EventInformation::GetNeutrinoPz() const
{
    return fnu_p4.pz()*GeV;
}

G4double EventInformation::GetNeutrinoE() const
{
    return fnu_p4.e()*GeV;  
}

G4double EventInformation::GetVertexX() const
{
    return fvertex_pos.x()*mm;
}

G4double EventInformation::GetVertexY() const
{
    return fvertex_pos.y()*mm;
}

G4double EventInformation::GetVertexZ() const
{
    return fvertex_pos.z()*mm;
}


void EventInformation::Print() const
{
    G4cout << "EventInformation:\n"\
    << "Target: " << ftarget_pdgc << "\n"\
    << "Neutrino: " << fnu_pdgc << " " << fnu_p4 << "\n"\
    << "Vertex: " << fvertex_pos << G4endl;
}