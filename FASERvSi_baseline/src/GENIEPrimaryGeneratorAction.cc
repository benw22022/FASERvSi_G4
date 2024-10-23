#include "GENIEPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction.hh"

#include <G4GeneralParticleSource.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include <TMath.h>
#include <TTree.h>
#include <TFile.h>

#include <vector>

GENIEPrimaryGeneratorAction::GENIEPrimaryGeneratorAction(G4GeneralParticleSource* gps)
  : fGPS(gps)
{
}

GENIEPrimaryGeneratorAction::~GENIEPrimaryGeneratorAction()
{
}

void GENIEPrimaryGeneratorAction::ShootParticle(G4Event* anEvent, G4int pdg, TLorentzVector p4)
{
  // unknown pgd codes in GENIE --> skip it!
  // ref: https://internal.dunescience.org/doxygen/ConvertMCTruthToG4_8cxx_source.html
  // This has been a known issue with GENIE
  const int genieLo = 2000000001;
  const int genieHi = 2000000202;
  if ( pdg >= genieLo && pdg <= genieHi) {
    std::cout<<"This unknown PDG code ["<<pdg<<"] was present in the GENIE input, "
             <<"but will not be processed by Geant4."
             <<std::endl;
    return;
  }
  
  G4ParticleDefinition* particleDefinition;
  if ( pdg == 0) {
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
  } else {
    particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
  }
  
  if ( pdg > 1000000000) { // If the particle is a nucleus
    // If the particle table doesn't have a definition yet, ask the ion
    // table for one. This will create a new ion definition as needed.
    if (!particleDefinition) {
      int Z = (pdg % 10000000) / 10000; // atomic number
      int A = (pdg % 10000) / 10;       // mass number
      particleDefinition = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z, A, 0.);
    }
  }
 
  /*G4cout << "Particle PDG " << pdg << " mass " << particleDefinition->GetPDGMass()*MeV << G4endl;
  G4cout << "p4  " << p4.X() << " " << p4.Y() << " " << p4.Z() << " " << p4.E() << G4endl;
  G4cout << "kinE " << (p4.E()*GeV - particleDefinition->GetPDGMass()*MeV)  << G4endl;*/

  // load the gun...
  fGPS->SetParticleDefinition(particleDefinition);
  fGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(( p4.E()*GeV - particleDefinition->GetPDGMass()*MeV));  // kinetic energy
  fGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(p4.X(), p4.Y(), p4.Z()));
  fGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
  fGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords( G4ThreeVector(fneuX4.X(), fneuX4.Y(), fneuX4.Z()) );
  fGPS->GeneratePrimaryVertex(anEvent); // ...and shoot!

}

void GENIEPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent,
    G4String filename, G4int startIdx, G4int nuVtxOpt)
{
  fFileName    = filename;
  // fEvtStartIdx = startIdx;

  // int fevtID = anEvent->GetEventID();
  std::cout << "oooOOOooo Event # " << fevtID << " oooOOOooo" <<std::endl;
  G4cout << "GeneratePrimaries from file " << fFileName
    << ", fevtID starts from "<< fEvtStartIdx
    << ", now at " << fEvtStartIdx + fevtID<<G4endl;

  fgFaserFile = new TFile(fFileName, "read");
  if (!fgFaserFile->IsOpen()) {
    std::cout<<"Cannot open gFaser file : "<<fFileName<<std::endl;
    exit(1);
  }
  fgFaserTree = (TTree*)fgFaserFile->Get("gFaser");
  if (!fgFaserTree) {
    std::cout<<"No gFaser event tree in input file : "<<fFileName<<std::endl;
    exit(1);
  }
  int nev = fgFaserTree->GetEntries();
  std::cout<<"Input gFaser tree has "<<nev<<((nev==1)? " entry." : " entries.")<<std::endl;
  
  // define the branches we are interested in from the  tree
  const int kNPmax = 250;

  double vx, vy, vz;                                  // vertex position
  std::vector<int>    *pdgc = nullptr;                             // PDG code
  std::vector<double> *E = nullptr;                             // Particle 4-momenta
  std::vector<double> *px = nullptr; 
  std::vector<double> *py = nullptr; 
  std::vector<double> *pz = nullptr;                            
  std::vector<int>    *status = nullptr;                           // Genie status code; 0 == initial state, 1 == stable final state
  int nFinalStateParticles;                           // Number of stable final state particles
  // std::vector<int> firstDaughter, lastDaughter;       // TODO: Keep track of truth record 
  // std::vector<int> firstMother, lastMother;

  fgFaserTree->SetBranchAddress("vx",&vx);
  fgFaserTree->SetBranchAddress("vy",&vy);
  fgFaserTree->SetBranchAddress("vz",&vz);
  fgFaserTree->SetBranchAddress("E",&E);
  fgFaserTree->SetBranchAddress("px",&px);
  fgFaserTree->SetBranchAddress("py",&py);
  fgFaserTree->SetBranchAddress("pz",&pz);
  fgFaserTree->SetBranchAddress("pdgc",&pdgc);
  fgFaserTree->SetBranchAddress("status",&status);
  // fgFaserTree->SetBranchAddress("firstMother",&firstMother);
  // fgFaserTree->SetBranchAddress("lastMother",&lastMother);

  if (fEvtStartIdx + fevtID >= nev) {
    std::cout<<"** event index beyond range !! **"<<std::endl;
  }
  
  // fetch a single entry from GENIE input file
  std::cout<<"** Getting entry " << fEvtStartIdx + fevtID<<  "  **"<<std::endl;
  fgFaserTree->GetEntry(fEvtStartIdx + fevtID); 
  std::cout<<"** Got entry " << fEvtStartIdx + fevtID <<  "  **"<<std::endl;

  // compute/repackage what is not directly available from the tree
  // position is randomly extracted in the detector fiducial volume
  // or set set to the center according to config parameter
 
  fneuIdx = fEvtStartIdx + fevtID;
  // fneuP4.SetPxPyPzE(pxv,pyv,pzv,Ev);
  // ffslP4.SetPxPyPzE(pxl,pyl,pzl,El);
  
  // fint_type = DecodeInteractionType(cc,nc,em);
  // fscattering_type = DecodeScatteringType(qel,res,dis,coh,dfr,imd,imdanh,mec,nuel,singlek,amnugamma);

  fneuX4.SetX(vx);
  fneuX4.SetY(vy);
  fneuX4.SetZ(vz);
  fneuX4.SetT(0.);

  // for (int ipar=0; ipar<nFinalStateParticles; ++ipar) {

  //   // Only fire the stable final state particles from the gun
  //   if (status.at(ipar) == 1){
  //     TLorentzVector p( px.at(ipar), py.at(ipar), pz.at(ipar), E.at(ipar) );
  //     std::cout << "Shooting a particle" << std::endl;
  //     ShootParticle( anEvent, pdgc.at(ipar), p);
  //   }
  // }
  for (int ipar=0; ipar<nFinalStateParticles; ++ipar) {

    // Only fire the stable final state particles from the gun
    if (status->at(ipar) == 1){
      TLorentzVector p( px->at(ipar), py->at(ipar), pz->at(ipar), E->at(ipar) );
      std::cout << "Shooting a particle" << std::endl;
      ShootParticle( anEvent, pdgc->at(ipar), p);
    }
  }

  fevtID++;
  fgFaserFile->Close();


  // TODO: Maybe add this option back in?
  // G4Random::setTheSeed(fEvtStartIdx+fevtID+1);
  // switch(nuVtxOpt) {
  //   case 0:
  //     fneuX4.SetX(0.*m);
  //     fneuX4.SetY(0.*m);
  //     fneuX4.SetZ(GeometricalParameters::Get()->GetFLArEPosition().z() -
  //               GeometricalParameters::Get()->GetFLArEFidVolSize().z()/2);
  //     fneuX4.SetT(0.);
  //     break;
  //   case 1:
  //     fneuX4.SetX(GeometricalParameters::Get()->GetFLArEPosition().x() +
  //                (G4UniformRand()-0.5) * GeometricalParameters::Get()->GetFLArEFidVolSize().x());
  //     fneuX4.SetY(GeometricalParameters::Get()->GetFLArEPosition().y() +
  //                (G4UniformRand()-0.5) * GeometricalParameters::Get()->GetFLArEFidVolSize().y());
  //     fneuX4.SetZ(GeometricalParameters::Get()->GetFLArEPosition().z() +
  //                (G4UniformRand()-0.5) * GeometricalParameters::Get()->GetFLArEFidVolSize().z());
  //     fneuX4.SetT(0.);
  //     break;
  // }

  // now we shoot all the final state particles into the Geant4 event
  // - final state lepton (if NC, it's the neutrino!)
  // - all final state hadrons
  
}

// G4int GENIEPrimaryGeneratorAction::DecodeInteractionType(bool cc, bool nc, bool em)
// {
//   // using ref: https://internal.dunescience.org/doxygen/Generator_2src_2Framework_2Interaction_2InteractionType_8h_source.html
//   // returning the value of the corresponding enum member from the flags that are in output.
//   // there are no flags in the  tree for some of the interaction types...
 
//   enum InteractionType {
//     kIntNull   = 0,
//     kIntEM,          //
//     kIntWeakCC,      //
//     kIntWeakNC,      //
//     kIntWeakMix,     // CC + NC + interference
//     kIntDarkMatter,  //
//     kIntNDecay,      //
//     kIntNOsc,        //
//     kIntNHL,         //
//     kIntDarkNC       // 
//   };
	
//   if (cc) return kIntWeakCC;
//   else if (nc) return kIntWeakNC;
//   else if (em) return kIntEM;

//   return kIntNull;

// }

// G4int GENIEPrimaryGeneratorAction::DecodeScatteringType(bool qel, bool res, bool dis, bool coh, bool dfr, 
// 						       bool imd, bool imdanh, bool mec, bool nuel,
// 	       					       bool singlek, bool amnugamma)
// {
//   // using ref: https://internal.dunescience.org/doxygen/ScatteringType_8h_source.html
//   // returning the value of the corresponding enum member from the flags that are in output.
//   // there are no flags in the  tree for some of the scattering types...
  
//   enum ScatteringType {
//     kScUnknown = -100,
//     kScNull = 0,
//     kScQuasiElastic,
//     kScSingleKaon,
//     kScDeepInelastic,
//     kScResonant,
//     kScCoherentProduction,
//     kScDiffractive,
//     kScNuElectronElastic,
//     kScInverseMuDecay,
//     kScAMNuGamma,
//     kScMEC,
//     kScCoherentElastic,
//     kScInverseBetaDecay,
//     kScGlashowResonance,
//     kScIMDAnnihilation,
//     kScPhotonCoherent,
//     kScPhotonResonance,
//     kScSinglePion,
//     kScDarkMatterElastic = 101,
//     kScDarkMatterDeepInelastic,
//     kScDarkMatterElectron,
//     kScNorm
//   };

//   if(qel) return kScQuasiElastic;
//   else if(res) return kScResonant;
//   else if(dis) return kScDeepInelastic;
//   else if(coh) return kScCoherentProduction;
//   else if(dfr) return kScDiffractive;
//   else if(imd) return kScInverseMuDecay;
//   else if(imdanh) return kScIMDAnnihilation;
//   else if(mec) return kScMEC;
//   else if(nuel) return kScNuElectronElastic;
//   else if(singlek) return kScSingleKaon;
//   else if(amnugamma) return kScAMNuGamma;

//   return kScUnknown;
// }
