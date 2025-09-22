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
#include "FPFParticle.hh"
#include "SCTModuleHit.hh"

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
  fActsHitsTree = nullptr;
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

  //* Acts Hits Tree [i == unsigned int; F == float; l == Long unsigned 64 int]
  fActsHitsTree = new TTree("hits", "ActsHitsTree");
  fActsHitsTree->Branch("event_id", &ActsHitsEventID, "event_id/i");
  fActsHitsTree->Branch("geometry_id", &ActsHitsGeometryID, "geometryid/l");
  fActsHitsTree->Branch("particle_id", &ActsHitsParticleID, "particle_id/l");
  fActsHitsTree->Branch("tx", &ActsHitsX, "tx/F");
  fActsHitsTree->Branch("ty", &ActsHitsY, "ty/F");
  fActsHitsTree->Branch("tz", &ActsHitsZ, "tz/F");
  fActsHitsTree->Branch("tt", &ActsHitsT, "tt/F");
  fActsHitsTree->Branch("tpx", &ActsHitsPx, "tpx/F");
  fActsHitsTree->Branch("tpy", &ActsHitsPy, "tpy/F");
  fActsHitsTree->Branch("tpz", &ActsHitsPz, "tpz/F");
  fActsHitsTree->Branch("te", &ActsHitsE, "tpe/F");
  fActsHitsTree->Branch("deltapx", &ActsHitsDeltaPx, "deltapx/F");
  fActsHitsTree->Branch("deltapy", &ActsHitsDeltaPy, "deltapy/F");
  fActsHitsTree->Branch("deltapz", &ActsHitsDeltaPz, "deltapz/F");
  fActsHitsTree->Branch("deltae", &ActsHitsDeltaE, "deltae/F");
  fActsHitsTree->Branch("index", &ActsHitsIndex, "index/I");
  fActsHitsTree->Branch("volume_id", &ActsHitsVolumeID, "volume_id/i");
  fActsHitsTree->Branch("boundary_id", &ActsHitsBoundaryID, "boundary_id/i");
  fActsHitsTree->Branch("layer_id", &ActsHitsLayerID, "layer_id/i");
  fActsHitsTree->Branch("approach_id", &ActsHitsApproachID, "approach_id/i");
  fActsHitsTree->Branch("sensitive_id", &ActsHitsSensitiveID, "sensitive_id/i");

  //* Acts truth particle tree
  fActsParticlesTree = new TTree("particles", "ActsParticlesTree");
  fActsParticlesTree->Branch("event_id", &ActsHitsEventID, "event_id/i");
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
  fActsHitsTree->Write();
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
  FillHitsOutput();

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

  //-----------------------------------------------------------

  // FillHitsOutput();

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

  // loop over the detected Hits sensitive volumes
  int nHits = 0;
  auto sdManager = G4SDManager::GetSDMpointer();
  G4cout << "getting sdManager pointer" << G4endl;
  G4int sdId = sdManager->GetCollectionID("strip_detector");
  G4cout << "Looking for hit collection with ID " << sdId << G4endl;
  
  G4cout  << "fHCofEvent->GetHC(" << sdId<<") = " << fHCofEvent->GetHC(sdId) << G4endl;
  auto hitCollection = dynamic_cast<G4THitsCollection<SCTModuleHit>*>(fHCofEvent->GetHC(sdId));
  G4cout << "Done dynamic cast of hist collection" << G4endl;

  G4cout << "Found hit collection with ID " << sdId << G4endl;


  if (!hitCollection)
  {
    G4cout << "No hits recorded by " << "strip_detector" << G4endl;
    return;
  }
  
  std::map<G4int, G4int> sub_part_map{};
  for (auto hit : *hitCollection->GetVector())
  {
    std::cout << "Processing hit from track ID " << hit->GetTrackID() << " with PDG " << hit->GetPDGID() << " and charge " << hit->GetCharge() << std::endl;
    if (hit->GetCharge() == 0)
      continue; // skip neutral particles, they don't hit

    /*
    * A note on the ActsHitsGeometryID variable
      This variable in Acts keeps track of an Acts::GeometryIdentifier. This is essentially a long unsigned int, the bits of which are used to
      look up the the volume/layer/boundary/sensitive indices of a piece of geometry. In principle it should be possible to assign this variable
      here in GEANT4 but I don't understand the Acts code well enough to do it without adding Acts as a dependancy to this codebase.
      As a result I set `geometry_id` to zero and give the the resposibility of assigning this variable to the user during the reading of the `hits` tree.
    */

    nHits++;
    ActsHitsEventID = evtID;
    ActsHitsGeometryID = 0;

    int hitID = hit->GetTrackID();
    int nPrimaries = ActsParticlesParticleId.size();

    auto particleId = ActsFatras::Barcode();
    particleId.setVertexPrimary(1);
    particleId.setVertexSecondary(0);
    particleId.setParticle(hit->GetTrackID() - 1); // The track ID is the primary particle index plus one
    particleId.setGeneration(hit->GetParentID());

    sub_part_map.try_emplace(hit->GetTrackID() - 1, sub_part_map.size());

    // This is a fudge - assumes that that the secondary particles are always sub-particles of the primary particle
    particleId.setSubParticle(hit->GetParentID() == 0 ? 0 : sub_part_map[hit->GetTrackID() - 1]);
    ActsHitsParticleID = particleId.value();

    ActsHitsX = hit->GetX();
    ActsHitsY = hit->GetY();
    ActsHitsZ = hit->GetZ();
    ActsHitsT = hit->GetT();
    ActsHitsPx = hit->GetPx();
    ActsHitsPy = hit->GetPy();
    ActsHitsPz = hit->GetPz();
    ActsHitsE = hit->GetEnergy();
    ActsHitsDeltaPx = hit->GetDeltaPx();
    ActsHitsDeltaPy = hit->GetDeltaPy();
    ActsHitsDeltaPz = hit->GetDeltaPz();
    ActsHitsDeltaE = hit->GetDeltaE();
    ActsHitsIndex = hit->GetCopyNumSensor(); // index of layer: 0, 1, 2, ...

    // These variables I'm not 100% sure about. I reverse engineered them by matching them to how they're set when writing the hits from the particle gun in Acts
    // In principle with the right headers from Acts we could construct the geometry ID value here
    ActsHitsVolumeID = 1;
    ActsHitsBoundaryID = 0;
    ActsHitsLayerID = (hit->GetCopyNumSensor() + 1) * 2; // Acts specfic layer ID, goes 2, 4, 6, ...
    ActsHitsApproachID = 0;
    ActsHitsSensitiveID = 1;
    fActsHitsTree->Fill();

    // Now fill the Acts particles tree
    bool isDuplicate = false;
    for (const auto &id : ActsParticlesParticleId)
    {
      if (id == particleId.value())
      {
        isDuplicate = true;
      }
    }
    if (isDuplicate) continue; // Skip this particle if it's already been added

    ActsParticlesParticleId.push_back(particleId.value());
    ActsParticlesParticleType.push_back(hit->GetPDGID());
    ActsParticlesProcess.push_back(0);
    ActsParticlesVx.push_back(hit->GetTrackVertex().x());
    ActsParticlesVy.push_back(hit->GetTrackVertex().y());
    ActsParticlesVz.push_back(hit->GetTrackVertex().z());
    ActsParticlesVt.push_back(0);
    ActsParticlesPx.push_back(hit->GetTrackP4().px());
    ActsParticlesPy.push_back(hit->GetTrackP4().py());
    ActsParticlesPz.push_back(hit->GetTrackP4().pz());
    ActsParticlesM.push_back(hit->GetTrackP4().m());
    ActsParticlesQ.push_back(hit->GetCharge());

    ActsParticlesEta.push_back(hit->GetTrackP4().eta());
    ActsParticlesPhi.push_back(hit->GetTrackP4().phi());
    ActsParticlesPt.push_back(pow(pow(hit->GetTrackP4().px(), 2) + pow(hit->GetTrackP4().py(), 2), 0.5));
    ActsParticlesP.push_back(pow(pow(hit->GetTrackP4().px(), 2) + pow(hit->GetTrackP4().py(), 2) + pow(hit->GetTrackP4().pz(), 2), 0.5));
    ActsParticlesVertexPrimary.push_back(hit->GetIsPrimaryTrack());     //? These variables need to be filled, but are unused by Acts
    ActsParticlesVertexSecondary.push_back(hit->GetIsSecondaryTrack()); //? These variables need to be filled, but are unused by Acts
    ActsParticlesParticle.push_back(1);                                 //? These variables need to be filled, but are unused by Acts
    ActsParticlesGeneration.push_back(0);                               //? These variables need to be filled, but are unused by Acts
    ActsParticlesSubParticle.push_back(0);                              //? These variables need to be filled, but are unused by Acts
    ActsParticlesELoss.push_back(0);                                    //? These variables need to be filled, but are unused by Acts
    ActsParticlesPathInX0.push_back(0);                                 //? These variables need to be filled, but are unused by Acts
    ActsParticlesPathInL0.push_back(0);                                 //? These variables need to be filled, but are unused by Acts
    ActsParticlesNumberOfHits.push_back(0);                             //? These variables need to be filled, but are unused by Acts
    ActsParticlesOutcome.push_back(0);                                  //? These variables need to be filled, but are unused by Acts
  } // end of loop over hits
  fActsParticlesTree->Fill();

  G4cout << "Total Hits recorded hits: " << nHits << G4endl;
}

float_t AnalysisManager::GetTotalEnergy(float_t px, float_t py, float_t pz, float_t m)
{
  return TMath::Sqrt(px * px + py * py + pz * pz + m * m);
}
