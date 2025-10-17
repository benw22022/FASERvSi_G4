#include <vector>
#include <functional>
#include <iostream>
#include <string>
#include <map>
#include <iomanip>
#include <random>

#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <G4Poisson.hh>
#include <G4Trajectory.hh>
#include <G4LorentzVector.hh>
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"


#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TString.h>
#include <Math/ProbFunc.h>

#include "EventInformation.hh"
#include "AnalysisManager.hh"
#include "reco/Barcode.hh"
// #include "reco/SpacePoint.hh"
#include "FPFParticle.hh"
#include "SCTModuleHit.hh"
#include "TruthHit.hh"

//---------------------------------------------------------------------
//---------------------------------------------------------------------
// AnalysisManager "singleton" instance
// once initialized, can be used to point to AnalysisManager
// from anywhere else in the codebase
AnalysisManager *AnalysisManager::fInstance = 0;

AnalysisManager *AnalysisManager::GetInstance()
{
  if (!fInstance)
  {
    G4cout << "AnalysisManager: Re-initialization" << G4endl;
    fInstance = new AnalysisManager();
  }
  return fInstance;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
AnalysisManager::AnalysisManager()
{
  fFile = nullptr;
  fFilename = "test.root";

  fMessenger = new AnalysisManagerMessenger(this);

  fEvt = nullptr;
  fTrk = nullptr;
  fPrim = nullptr;
  fRecoHitsTree = nullptr;
  fTruthHitsTree = nullptr;
  fActsParticlesTree = nullptr;
  
  fSaveTrack = false;
}

AnalysisManager::~AnalysisManager() {}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookEvtTree()
{
  fEvt = new TTree("event", "event info");
  fEvt->Branch("evtID", &evtID, "evtID/I");
  fEvt->Branch("vtxID", &vertexID, "vtxID/I");
  fEvt->Branch("weight", &weight, "weight/D");
  fEvt->Branch("genType", &genType);
  fEvt->Branch("processName", &processName);
  fEvt->Branch("initPDG", &initPDG, "initPDG/I");
  fEvt->Branch("initX", &initX, "initX/D");
  fEvt->Branch("initY", &initY, "initY/D");
  fEvt->Branch("initZ", &initZ, "initZ/D");
  fEvt->Branch("initT", &initT, "initT/D");
  fEvt->Branch("initPx", &initPx, "initPx/D");
  fEvt->Branch("initPy", &initPy, "initPy/D");
  fEvt->Branch("initPz", &initPz, "initPz/D"); 
  fEvt->Branch("initE", &initE, "initE/D");
  fEvt->Branch("initM", &initM, "initM/D");
  fEvt->Branch("initQ", &initQ, "initQ/D");
  fEvt->Branch("intType", &intType, "intType/I");
  fEvt->Branch("scatteringType", &scatteringType, "scatteringType/I");
  fEvt->Branch("fslPDG", &fslPDG, "fslPDG/I");
  fEvt->Branch("tgtPDG", &tgtPDG, "tgtPDG/I");
  fEvt->Branch("tgtA", &tgtA, "tgtA/I");
  fEvt->Branch("tgtZ", &tgtZ, "tgtZ/I");
  fEvt->Branch("hitnucPDG", &hitnucPDG, "hitnucPDG/I");
  fEvt->Branch("xs", &xs, "xs/D");
  fEvt->Branch("Q2", &Q2, "Q2/D");
  fEvt->Branch("xBj", &xBj, "xBj/D");
  fEvt->Branch("y", &y, "y/D");
  fEvt->Branch("W", &W, "W/D");
}

void AnalysisManager::bookPrimTree()
{
  fPrim = new TTree("primaries", "primaries info");
  fPrim->Branch("evtID", &evtID, "evtID/I");
  fPrim->Branch("vtxID", &primVtxID, "vtxID/I");
  fPrim->Branch("PDG", &primPDG, "PDG/I");
  fPrim->Branch("trackID", &primTrackID, "trackID/I");
  fPrim->Branch("barcode", &primParticleID, "bardcode/I");
  fPrim->Branch("mass", &primM, "mass/F");
  fPrim->Branch("charge", &primQ, "charge/F");
  fPrim->Branch("Vx", &primVx, "Vx/F"); // position
  fPrim->Branch("Vy", &primVy, "Vy/F");
  fPrim->Branch("Vz", &primVz, "Vz/F");
  fPrim->Branch("Vt", &primVt, "Vt/F");
  fPrim->Branch("Px", &primPx, "Px/F"); // momentum
  fPrim->Branch("Py", &primPy, "Py/F");
  fPrim->Branch("Pz", &primPz, "Pz/F");
  fPrim->Branch("E", &primE, "E/F");    // initial total energy
  fPrim->Branch("KE", &primKE, "KE/F"); // initial kinetic energy
  fPrim->Branch("Eta", &primEta, "Eta/F");
  fPrim->Branch("Phi", &primPhi, "Phi/F");
  fPrim->Branch("Pt", &primPt, "Pt/F");
  fPrim->Branch("P", &primP, "P/F");
  fPrim->Branch("Theta", &primTheta, "Theta/F");
}

void AnalysisManager::bookTrkTree()
{
  fTrk = new TTree("trajectories", "trajectories info");
  fTrk->Branch("evtID", &evtID, "evtID/I");
  fTrk->Branch("trackTID", &trackTID, "trackTID/I");
  fTrk->Branch("trackPID", &trackPID, "trackPID/I");
  fTrk->Branch("trackPDG", &trackPDG, "trackPDG/I");
  fTrk->Branch("trackKinE", &trackKinE, "trackKinE/D");
  fTrk->Branch("trackNPoints", &trackNPoints, "trackNPoints/I");
  fTrk->Branch("trackPointX", &trackPointX);
  fTrk->Branch("trackPointY", &trackPointY);
  fTrk->Branch("trackPointZ", &trackPointZ);
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::bookHitsTrees()
{
  // create subdirectory in file
  fHits = fFile->mkdir("Hits","Hits output",kTRUE);
  fFile->cd(fHits->GetName());

  //* Reco Hits Tree [i == unsigned int; F == float; l == Long unsigned 64 int]
  fRecoHitsTree = new TTree("recoHits", "recoHitsTree");
  fRecoHitsTree->Branch("event_id", &recoHitsEventID, "event_id/i");
  fRecoHitsTree->Branch("hit_x", &recoHitsX);
  fRecoHitsTree->Branch("hit_y", &recoHitsY);
  fRecoHitsTree->Branch("hit_z", &recoHitsZ);
  fRecoHitsTree->Branch("hit_pdgc", &recoHitsPDGC);
  fRecoHitsTree->Branch("hit_isTruthMatched", &recoHitsIsTruthMatched);
  fRecoHitsTree->Branch("hit_moduleID", &recoHitsModuleID);
  fRecoHitsTree->Branch("hit_layerID", &recoHitsLayerID);
  fRecoHitsTree->Branch("hit_trackID", &recoHitsTrackID);
  fRecoHitsTree->Branch("hit_parentID", &recoHitsParentID);
  fRecoHitsTree->Branch("hit_truthHitID", &recoHitsTruthHitID);
  fRecoHitsTree->Branch("hit_px", &recoHitsPx);
  fRecoHitsTree->Branch("hit_py", &recoHitsPy);
  fRecoHitsTree->Branch("hit_pz", &recoHitsPz);
  fRecoHitsTree->Branch("hit_energy", &recoHitsE);
  fRecoHitsTree->Branch("hit_mass", &recoHitsMass);
  fRecoHitsTree->Branch("hit_charge", &recoHitsCharge);

  //* Truth Hits Tree [i == unsigned int; F == float; l == Long unsigned 64 int]
  fTruthHitsTree = new TTree("truthHits", "truthHitsTree");
  fTruthHitsTree->Branch("event_id", &truthHitsEventID, "event_id/i");
  fTruthHitsTree->Branch("hit_x", &truthHitsX);
  fTruthHitsTree->Branch("hit_y", &truthHitsY);
  fTruthHitsTree->Branch("hit_z", &truthHitsZ);
  fTruthHitsTree->Branch("hit_pdgc", &truthHitsPDGC);
  fTruthHitsTree->Branch("hit_moduleID", &truthHitsModuleID);
  fTruthHitsTree->Branch("hit_layerID", &truthHitsLayerID);
  fTruthHitsTree->Branch("hit_trackID", &truthHitsTrackID);
  fTruthHitsTree->Branch("hit_parentID", &truthHitsParentID);
  fTruthHitsTree->Branch("hit_px", &truthHitsPx);
  fTruthHitsTree->Branch("hit_py", &truthHitsPy);
  fTruthHitsTree->Branch("hit_pz", &truthHitsPz);
  fTruthHitsTree->Branch("hit_energy", &truthHitsE);
  fTruthHitsTree->Branch("hit_mass", &truthHitsMass);
  fTruthHitsTree->Branch("hit_charge", &truthHitsCharge);
  fTruthHitsTree->Branch("hit_theta", &truthHitsTheta);
  fTruthHitsTree->Branch("hit_truthHitID", &truthHitsID);


  //* Acts truth particle tree
  fActsParticlesTree = new TTree("particles", "ActsParticlesTree");
  fActsParticlesTree->Branch("event_id", &recoHitsEventID, "event_id/i");
  fActsParticlesTree->Branch("particle_id", &ActsParticlesParticleId);
  fActsParticlesTree->Branch("particle_type", &ActsParticlesParticleType);
  fActsParticlesTree->Branch("process", &ActsParticlesProcess);
  fActsParticlesTree->Branch("vx", &ActsParticlesVx);
  fActsParticlesTree->Branch("vy", &ActsParticlesVy);
  fActsParticlesTree->Branch("vz", &ActsParticlesVz);
  fActsParticlesTree->Branch("vt", &ActsParticlesVt);
  fActsParticlesTree->Branch("px", &ActsParticlesPx);
  fActsParticlesTree->Branch("py", &ActsParticlesPy);
  fActsParticlesTree->Branch("pz", &ActsParticlesPz);
  fActsParticlesTree->Branch("m", &ActsParticlesM);
  fActsParticlesTree->Branch("q", &ActsParticlesQ);
  fActsParticlesTree->Branch("eta", &ActsParticlesEta);
  fActsParticlesTree->Branch("phi", &ActsParticlesPhi);
  fActsParticlesTree->Branch("pt", &ActsParticlesPt);
  fActsParticlesTree->Branch("p", &ActsParticlesP);
  fActsParticlesTree->Branch("vertex_primary", &ActsParticlesVertexPrimary);
  fActsParticlesTree->Branch("vertex_secondary", &ActsParticlesVertexSecondary);
  fActsParticlesTree->Branch("particle", &ActsParticlesParticle);
  fActsParticlesTree->Branch("generation", &ActsParticlesGeneration);
  fActsParticlesTree->Branch("sub_particle", &ActsParticlesSubParticle);
  fActsParticlesTree->Branch("e_loss", &ActsParticlesELoss);
  fActsParticlesTree->Branch("total_x0", &ActsParticlesPathInX0);
  fActsParticlesTree->Branch("total_l0", &ActsParticlesPathInL0);
  fActsParticlesTree->Branch("number_of_hits", &ActsParticlesNumberOfHits);
  fActsParticlesTree->Branch("outcome", &ActsParticlesOutcome);

  fFile->cd();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::BeginOfRun()
{
  G4cout << "Run has been started, preparing output" << G4endl;

  if (fFile)
    delete fFile;

  // Preparing output file
  fFile = new TFile(fFilename.c_str(), "RECREATE");
  
  // Booking common output trees
  bookEvtTree();
  bookPrimTree();
  if (fSaveTrack) bookTrkTree();

  bookHitsTrees();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::EndOfRun()
{
  G4cout << "Run has ended, closing output" << G4endl;
  // save common trees at the top of the output file
  fFile->cd();
  fEvt->Write();
  fPrim->Write();
  if (fSaveTrack) fTrk->Write();

  fFile->cd(fHits->GetName());
  fRecoHitsTree->Write();
  fTruthHitsTree->Write();
  fActsParticlesTree->Write();
  fFile->cd(); // go back to top

  fFile->Close();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::BeginOfEvent()
{
  G4cout << "Starting new event, resetting variables" << G4endl;
  // reset vectors that need to be cleared for a new event
  // only reset arrays or vectors, tipically no need for other defaults

  primaries.clear();
  primaryIDs.clear();

  // track ID to primary ancestor association
  trackToPrimaryAncestor.clear();

  trackPointX.clear();
  trackPointY.clear();
  trackPointZ.clear();


  recoHitsX.clear();
  recoHitsY.clear();
  recoHitsZ.clear();
  recoHitsPDGC.clear();
  recoHitsIsTruthMatched.clear();
  recoHitsModuleID.clear();
  recoHitsLayerID.clear();
  recoHitsTrackID.clear();
  recoHitsParentID.clear();
  recoHitsTruthHitID.clear();
  recoHitsPx.clear();
  recoHitsPy.clear();
  recoHitsPz.clear();
  recoHitsE.clear();
  recoHitsCharge.clear();
  recoHitsMass.clear();

  truthHitsX.clear();
  truthHitsY.clear();
  truthHitsZ.clear();
  truthHitsPDGC.clear();
  truthHitsModuleID.clear();
  truthHitsLayerID.clear();
  truthHitsTrackID.clear();
  truthHitsParentID.clear();
  truthHitsPx.clear();
  truthHitsPy.clear();
  truthHitsPz.clear();
  truthHitsE.clear();
  truthHitsMass.clear();
  truthHitsCharge.clear();
  truthHitsID.clear();
  truthHitsTheta.clear();


  ActsParticlesParticleId.clear();
  ActsParticlesParticleType.clear();
  ActsParticlesProcess.clear();
  ActsParticlesVx.clear();
  ActsParticlesVy.clear();
  ActsParticlesVz.clear();
  ActsParticlesVt.clear();
  ActsParticlesPx.clear();
  ActsParticlesPy.clear();
  ActsParticlesPz.clear();
  ActsParticlesM.clear();
  ActsParticlesQ.clear();
  ActsParticlesEta.clear();
  ActsParticlesPhi.clear();
  ActsParticlesPt.clear();
  ActsParticlesP.clear();
  ActsParticlesVertexPrimary.clear();
  ActsParticlesVertexSecondary.clear();
  ActsParticlesParticle.clear();
  ActsParticlesGeneration.clear();
  ActsParticlesSubParticle.clear();
  ActsParticlesELoss.clear();
  ActsParticlesPathInX0.clear();
  ActsParticlesPathInL0.clear();
  ActsParticlesNumberOfHits.clear();
  ActsParticlesOutcome.clear();
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::EndOfEvent(const G4Event *event)
{
  G4cout << "Ending event, filling output trees" << G4endl;
  /// evtID
  evtID = event->GetEventID();

  // FILL EVENT TREE
  FillEventTree(event);

  //-----------------------------------------------------------

  // FILL PRIMARIES/TRAJECTORIES TREE
  FillPrimariesTree(event);
  if(fSaveTrack) FillTrajectoriesTree(event);

  //-----------------------------------------------------------

  // Get the hit collections
  // If there is no hit collection, there is nothing to be done
  fHCofEvent = event->GetHCofThisEvent();
  if (!fHCofEvent)
  {
    G4cout << "No hits recorded in any sensitive volume --> nothing to save!" << G4endl;
    return;
  }

  FillHitsOutput();

}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillEventTree(const G4Event *event)
{
  G4cout << "Filling event tree" << G4endl;
  EventInformation* eventInfo = static_cast<EventInformation*>(event->GetUserInformation());
  eventInfo->Print();
  auto metadata = eventInfo->GetEventMetadata();
  for(int i=0; i<metadata.size(); i++)
  {
    vertexID = i;
    weight = metadata[i].weight;
    genType = metadata[i].generatorType;
    processName = metadata[i].processName;
    initPDG = metadata[i].pdg;
    initX = metadata[i].x4.x();
    initY = metadata[i].x4.y();
    initZ = metadata[i].x4.z();
    initT = metadata[i].x4.t();
    initPx = metadata[i].p4.x();
    initPy = metadata[i].p4.y();
    initPz = metadata[i].p4.z();
    initE = metadata[i].p4.e();
    initM = metadata[i].mass;
    initQ = metadata[i].charge;
    intType = metadata[i].intType;     
    scatteringType = metadata[i].scatteringType;   
    fslPDG = metadata[i].fsl_pdg;           
    tgtPDG = metadata[i].tgt_pdg;  
    tgtZ = metadata[i].tgt_Z;     
    tgtA = metadata[i].tgt_A;     
    hitnucPDG = metadata[i].hitnuc_pdg;  
    xs = metadata[i].xs;
    Q2 = metadata[i].Q2;  
    xBj = metadata[i].xBj;
    y = metadata[i].y; 
    W = metadata[i].W; 

    fEvt->Fill();
  }
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillPrimariesTree(const G4Event *event)
{
  G4cout << "Filling primaries tree" << G4endl;
  nPrimaryVertex = event->GetNumberOfPrimaryVertex();
  G4cout << "\nNumber of primary vertices  : " << nPrimaryVertex << G4endl;
  
  /// loop over the vertices, and then over primary particles,
  /// neutrino truth info from event generator.
  for (G4int ivtx = 0; ivtx < event->GetNumberOfPrimaryVertex(); ++ivtx)
  {
    G4cout << "=== Vertex " << ivtx+1 << " of " << nPrimaryVertex << " -> " 
           << event->GetPrimaryVertex(ivtx)->GetNumberOfParticle() << " primaries ===" << G4endl;
    for (G4int ipp = 0; ipp < event->GetPrimaryVertex(ivtx)->GetNumberOfParticle(); ++ipp)
    {
      G4PrimaryParticle *primary_particle = event->GetPrimaryVertex(ivtx)->GetPrimary(ipp);
      if (primary_particle)
      {
 
        primVtxID = ivtx;
        primTrackID = ipp + 1; // confirm matches track id?

        auto particleId = ActsFatras::Barcode();
        particleId.setVertexPrimary(ivtx);
        particleId.setGeneration(0);
        particleId.setSubParticle(0);
        particleId.setParticle(primTrackID - 1);

        primParticleID = particleId.value();
        primPDG = primary_particle->GetPDGcode();
        primVx = event->GetPrimaryVertex(ivtx)->GetPosition().x();
        primVy = event->GetPrimaryVertex(ivtx)->GetPosition().y();
        primVz = event->GetPrimaryVertex(ivtx)->GetPosition().z();
        primVt = event->GetPrimaryVertex(ivtx)->GetT0();
        primPx = primary_particle->GetMomentum().x();
        primPy = primary_particle->GetMomentum().y();
        primPz = primary_particle->GetMomentum().z();
        primM = primary_particle->GetMass()/MeV;
        primQ = primary_particle->GetCharge();

        G4double energy = GetTotalEnergy(primPx, primPy, primPz, primM);
        G4LorentzVector p4(primPx,primPy,primPz,energy);
        primEta = p4.eta();
        primPhi = p4.phi();
        primPt = p4.perp();
        primP = p4.vect().mag();
        primE = energy;
        primKE = energy - primM;
        primTheta = p4.theta();

        // store a copy as a FPFParticle for further processing
        primaryIDs.push_back(primTrackID); //store to avoid duplicates
        primaries.push_back(FPFParticle(primPDG, 0, 
		                        primTrackID, primaryIDs.size()-1, 1,
		                        primM,
                            primVx, primVy, primVz, primVt,
                            primPx, primPy, primPz,energy));

        G4cout << G4endl;
        G4cout << "PrimaryParticleInfo: PDG code " << primPDG << G4endl
          << "Particle unique ID : " << primTrackID << G4endl
          << "Momentum : (" << primPx << ", " << primPy << ", " << primPz << ") MeV" << G4endl
          << "Vertex : (" << primVx << ", " << primVy << ", " << primVz << ") mm" << G4endl;

        fPrim->Fill();
      }
    }
  }

  G4cout << "\nNumber of primaries  : " << primaryIDs.size() << G4endl;
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillTrajectoriesTree(const G4Event* event)
{
  G4cout << "Filling trajectories tree" << G4endl;
  int count_tracks = 0;

  G4cout << "==== Saving track information to tree ====" << G4endl; 
  auto trajectoryContainer = event->GetTrajectoryContainer(); 
  if (!trajectoryContainer)
  {
    G4cout << "No tracks found: did you enable their storage with '/tracking/storeTrajectory 1'?" << G4endl;
    return;
  }

  for (size_t i = 0; i < trajectoryContainer->entries(); ++i) 
  { 
    auto trajectory = static_cast<G4Trajectory*>((*trajectoryContainer)[i]); 
    trackTID = trajectory->GetTrackID();
    trackPID = trajectory->GetParentID();
    trackPDG = trajectory->GetPDGEncoding(); 
    trackKinE = trajectory->GetInitialKineticEnergy(); 
    trackNPoints = trajectory->GetPointEntries(); 
    count_tracks++; 
    for (size_t j = 0; j < trackNPoints; ++j) 
    { 
      G4ThreeVector pos = trajectory->GetPoint(j)->GetPosition(); 
      trackPointX.push_back( pos.x() );
      trackPointY.push_back( pos.y() );
      trackPointZ.push_back( pos.z() );
    }
    fTrk->Fill();
    trackPointX.clear(); 
    trackPointY.clear();
    trackPointZ.clear();
  }
  G4cout << "Total number of recorded track: " << count_tracks << std::endl;
}


//---------------------------------------------------------------------
//---------------------------------------------------------------------

void AnalysisManager::FillHitsOutput()
{
  G4cout << "==== Filling Hits output trees ====" << G4endl;
  int nHits = 0;
  auto* hc = fHCofEvent->GetHC(0);
  // G4int nHC = fHCofEvent->GetNumberOfCollections();
  // for (G4int i = 0; i < nHC; ++i) {
  //     auto* hc = fHCofEvent->GetHC(i);
  //     auto* recoHitCollection = dynamic_cast<SCTModuleHitCollection*>(hc);
  // //     if (recoHitCollection && recoHitCollection->GetName() == "RecoSpacePoints") {
        
  //       G4cout << "Found hit collection: " << recoHitCollection->GetName() << G4endl;
  //       std::map<G4int, G4int> sub_part_map{};
  //       for (auto hit : *recoHitCollection->GetVector())
  //       {
  //         // std::cout << "Processing hit from track ID " << hit->GetTrackID() << " with PDG " << hit->GetPDGID() << " and charge " << hit->GetCharge() << std::endl;
  //         if (hit->GetCharge() == 0) { continue; } // skip neutral particles, they don't hit

  //         nHits++;
  //         recoHitsEventID = evtID;
  //         recoHitsX.push_back(hit->GetX());
  //         recoHitsY.push_back(hit->GetY());
  //         recoHitsZ.push_back(hit->GetZ());
  //         recoHitsPx.push_back(hit->GetPx());
  //         recoHitsPy.push_back(hit->GetPy());
  //         recoHitsPz.push_back(hit->GetPz());
  //         recoHitsE.push_back(hit->GetEnergy());
  //         recoHitsMass.push_back(hit->GetMass());
  //         recoHitsCharge.push_back(hit->GetCharge());
  //         recoHitsPDGC.push_back(hit->GetPDGID());
  //         recoHitsIsTruthMatched.push_back(hit->GetIsTruthMatched());
  //         recoHitsModuleID.push_back(hit->GetModuleNumber());
  //         recoHitsLayerID.push_back(hit->GetLayerNumber());
  //         recoHitsTrackID.push_back(hit->GetTrackID());
  //         recoHitsParentID.push_back(hit->GetParentID());
  //         recoHitsTruthHitID.push_back(hit->GetTruthHitID());
          

  //         //* Here we compute a unique particle ID following the Acts convention
  //         //* (I'm reusing some code I wrote for a different project (FPFSim) here)
  //         //* This will allow us to get to write out the primary/secondary particles of the event 
  //         auto particleId = ActsFatras::Barcode();
  //         particleId.setVertexPrimary(1);
  //         particleId.setVertexSecondary(0);
  //         particleId.setParticle(hit->GetTrackID() - 1); // The track ID is the primary particle index plus one
  //         particleId.setGeneration(hit->GetParentID());

  //         sub_part_map.try_emplace(hit->GetTrackID() - 1, sub_part_map.size());

  //         // This is a fudge - assumes that that the secondary particles are always sub-particles of the primary particle
  //         particleId.setSubParticle(hit->GetParentID() == 0 ? 0 : sub_part_map[hit->GetTrackID() - 1]);
        

  //         // Now fill the Acts particles tree
  //         bool isDuplicate = false;
  //         for (const auto &id : ActsParticlesParticleId)
  //         {
  //           if (id == particleId.value())
  //           {
  //             isDuplicate = true;
  //           }
  //         }
  //         if (isDuplicate) continue; // Skip this particle if it's already been added

  //         ActsParticlesParticleId.push_back(particleId.value());
  //         ActsParticlesParticleType.push_back(hit->GetPDGID());
  //         ActsParticlesProcess.push_back(0);
  //         ActsParticlesVx.push_back(hit->GetTrackVertex().x());
  //         ActsParticlesVy.push_back(hit->GetTrackVertex().y());
  //         ActsParticlesVz.push_back(hit->GetTrackVertex().z());
  //         ActsParticlesVt.push_back(0);
  //         ActsParticlesPx.push_back(hit->GetTrackP4().px());
  //         ActsParticlesPy.push_back(hit->GetTrackP4().py());
  //         ActsParticlesPz.push_back(hit->GetTrackP4().pz());
  //         ActsParticlesM.push_back(hit->GetTrackP4().m());
  //         ActsParticlesQ.push_back(hit->GetCharge());

  //         ActsParticlesEta.push_back(hit->GetTrackP4().eta());
  //         ActsParticlesPhi.push_back(hit->GetTrackP4().phi());
  //         ActsParticlesPt.push_back(pow(pow(hit->GetTrackP4().px(), 2) + pow(hit->GetTrackP4().py(), 2), 0.5));
  //         ActsParticlesP.push_back(pow(pow(hit->GetTrackP4().px(), 2) + pow(hit->GetTrackP4().py(), 2) + pow(hit->GetTrackP4().pz(), 2), 0.5));
  //         ActsParticlesVertexPrimary.push_back(hit->GetIsPrimaryTrack());     //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesVertexSecondary.push_back(hit->GetIsSecondaryTrack()); //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesParticle.push_back(1);                                 //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesGeneration.push_back(0);                               //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesSubParticle.push_back(0);                              //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesELoss.push_back(0);                                    //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesPathInX0.push_back(0);                                 //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesPathInL0.push_back(0);                                 //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesNumberOfHits.push_back(0);                             //? These variables need to be filled, but are unused by Acts
  //         ActsParticlesOutcome.push_back(0);                                  //? These variables need to be filled, but are unused by Acts
  //       } // end of loop over hits
  //       fActsParticlesTree->Fill();
  //       fRecoHitsTree->Fill();

  //       G4cout << "Total number of reconstructed hits: " << nHits << G4endl;
  //   } // Close Reco Hits block
  
    //* Fill the truth hits tree
    // else
    // {
      auto* truthHitCollection = dynamic_cast<TruthHitCollection*>(hc);
      if (truthHitCollection && truthHitCollection->GetName() == "truth_tracker") {
    
        G4cout << "Found hit collection: " << truthHitCollection->GetName() << G4endl;
        nHits = 0;
        for (auto hit : *truthHitCollection->GetVector())
        {
          nHits++;
          truthHitsEventID = evtID;
          truthHitsX.push_back(hit->GetX());
          truthHitsY.push_back(hit->GetY());
          truthHitsZ.push_back(hit->GetZ());
          truthHitsPx.push_back(hit->GetPx());
          truthHitsPy.push_back(hit->GetPy());
          truthHitsPz.push_back(hit->GetPz());
          truthHitsE.push_back(hit->GetEnergy());
          truthHitsMass.push_back(hit->GetMass());
          truthHitsCharge.push_back(hit->GetCharge());
          truthHitsPDGC.push_back(hit->GetPDGID());
          truthHitsModuleID.push_back(hit->GetModuleNumber());
          truthHitsLayerID.push_back(hit->GetLayerNumber());
          truthHitsTrackID.push_back(hit->GetTrackID());
          truthHitsParentID.push_back(hit->GetParentID());
          truthHitsID.push_back(hit->GetTruthHitID());
          truthHitsTheta.push_back(hit->GetTrackP4().theta());
        } // end of loop over hits
        fTruthHitsTree->Fill();

        G4cout << "Total number of truth hits: " << nHits << G4endl;
    
    
      }
    } // Close Truth Hits block

  // } // Close loop over hit collections
// }

float_t AnalysisManager::GetTotalEnergy(float_t px, float_t py, float_t pz, float_t m)
{
  return TMath::Sqrt(px * px + py * py + pz * pz + m * m);
}
