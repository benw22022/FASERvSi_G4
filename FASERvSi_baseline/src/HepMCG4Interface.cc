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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMCG4Interface.cc
/// \brief Implementation of the HepMCG4Interface class
//
//

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Interface::HepMCG4Interface()
  : hepmcEvent(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HepMCG4Interface::~HepMCG4Interface()
{
  // delete *hepmcEvent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool HepMCG4Interface::CheckVertexInsideWorld
                         (const G4ThreeVector& pos) const
{
  G4Navigator* navigator= G4TransportationManager::GetTransportationManager()
                                                 -> GetNavigatorForTracking();

  G4VPhysicalVolume* world= navigator-> GetWorldVolume();
  G4VSolid* solid= world-> GetLogicalVolume()-> GetSolid();
  EInside qinside= solid-> Inside(pos);

  if( qinside != kInside) return false;
  else return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Interface::HepMC2G4(const std::shared_ptr<HepMC3::GenEvent> hepmcevt,
                                G4Event* g4event)
{
  // std::cout << "In HepMC2G4 " << std::endl;
  int vtx_counter{0};
  for (const auto& vertex : hepmcevt->vertices()) {
    
    // std::cout << "Checking out vertex " << vtx_counter << std::endl;
    vtx_counter++;

    // real vertex?
    G4bool qvtx=false;

    EventInformation* eventInfo = new EventInformation();

    int par_counter{0};
    for (const auto& particle : vertex->particles_in())  {

      if (particle->end_vertex() && particle->status()==4) {
        qvtx=true;
        if (EventInformation::isPDGNeutrino(particle->pdg_id()))
        {
          eventInfo->SetNeutrinoPDG(particle->pdg_id());
          HepMC3::FourVector particle_p4 = particle->momentum();
          G4LorentzVector p4(particle_p4.px(), particle_p4.py(), particle_p4.pz(), particle_p4.e());
          eventInfo->SetNeutrinoP4(p4);
        }
        else //! Assuming only 2-body collision
        {
          eventInfo->SetTargetPDG(particle->pdg_id());
        }
        // break;
      }
    }
    if (!qvtx) continue;

    // std::cout << "Found a vertex! " << vtx_counter << std::endl;

    // check world boundary
    HepMC3::FourVector pos = vertex->position();
    G4LorentzVector xvtx(pos.x(), pos.y(), pos.z(), pos.t());
    if (! CheckVertexInsideWorld(xvtx.vect()*mm))
    { 
      std::cout << "WARNING: tried to generate vertex outside of world volume!" << std::endl;
      std::cout << "WARNING: position was (" << pos.x() << ", "<<  pos.y() << ", " << pos.z() << ", " << pos.t() << ")" << std::endl;
      continue;
    }

    eventInfo->SetVertexPos(xvtx);

    // std::cout << "Vertex inside world " << vtx_counter << std::endl;

    // create G4PrimaryVertex and associated G4PrimaryParticles
    G4PrimaryVertex* g4vtx =
      new G4PrimaryVertex(xvtx.x()*mm, xvtx.y()*mm, xvtx.z()*mm,
                          xvtx.t()*mm/c_light);
      
    G4cout << "Placing vertex at (" <<  xvtx.x()*mm << ", " << xvtx.y()*mm << ", " << xvtx.z()*mm << ", " << xvtx.t()*mm/c_light << ")" << G4endl;

    for (const auto& particle : vertex->particles_out())  {
        if( particle->status() == 1 ||  particle->status() == 5 )
        {
          G4int pdgcode = particle->pdg_id();
          pos = particle->momentum();
          G4LorentzVector p(pos.px(), pos.py(), pos.pz(), pos.e());
          G4PrimaryParticle* g4prim = new G4PrimaryParticle(pdgcode, p.x()*GeV, p.y()*GeV, p.z()*GeV);

          // std::cout << "Setting primary particle: " << "pdgc = " << pdgcode << std::endl;
          g4vtx->SetPrimary(g4prim);

          if (particle->status() == 5) //! My special HEPMC code to indicate the leptons from nuCC interactions
          {
            eventInfo->SetIsCCInteraction(true);
            eventInfo->SetCCLeptonPDG(particle->pdg_id());
            eventInfo->SetCCLeptonP4(p);
          }
        }
    }
    
    // std::cout << "Setting primary vertex" << std::endl;
    
    g4event->SetUserInformation(eventInfo);
    g4event->AddPrimaryVertex(g4vtx);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::shared_ptr<HepMC3::GenEvent> HepMCG4Interface::GenerateHepMCEvent()
{
  std::shared_ptr<HepMC3::GenEvent> anevent = std::make_shared<HepMC3::GenEvent>();
  return anevent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HepMCG4Interface::GeneratePrimaryVertex(G4Event* anEvent)
{
  // delete previous event object
  // delete *hepmcEvent;

  // generate next event
  hepmcEvent = GenerateHepMCEvent();
  if(! hepmcEvent) {
    G4cout << "HepMCInterface: no generated particles. run terminated..."
           << G4endl;
    G4RunManager::GetRunManager()-> AbortRun();
    return;
  }
  HepMC2G4(hepmcEvent, anEvent);
}
