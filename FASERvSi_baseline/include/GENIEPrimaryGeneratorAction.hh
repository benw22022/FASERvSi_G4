#ifndef GENIEPrimaryGeneratorAction_HH
#define GENIEPrimaryGeneratorAction_HH

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

class GENIEPrimaryGeneratorAction{
  public:
    GENIEPrimaryGeneratorAction(G4GeneralParticleSource*);
    ~GENIEPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent, G4String filename, G4int startIdx, G4int nuVtxOpt);
    void ShootParticle(G4Event* anEvent, G4int pdg, TLorentzVector p4);
    
    // G4int DecodeInteractionType(bool cc, bool nc, bool em);
    // G4int DecodeScatteringType(bool qel, bool res, bool dis, bool coh, bool dfr, 
		// 	       bool imd, bool imdanh, bool mec, bool nuel,
	  //                      bool singlek, bool amnugamma);

    G4int             GetNeuIdx() { return fneuIdx; };
    G4int             GetNeuPDG() { return fneuPDG; };
    TLorentzVector    GetNeuP4()  { return fneuP4; };
    TLorentzVector    GetNeuX4()  { return fneuX4; };
    G4int GetInteractionTypeId() { return fint_type; };
    G4int GetScatteringTypeId()  { return fscattering_type; };
    G4int             GetFSLPDG() { return ffslPDG; };
    TLorentzVector    GetFSLP4()  { return ffslP4; };
    G4double          GetGetW()   { return fW; };
    
  private:
    G4GeneralParticleSource* fGPS;

    G4String fFileName;
    G4int  fEvtStartIdx{0};
    G4int  fevtID{0};
    TFile* fgFaserFile;
    TTree* fgFaserTree;

    G4int fneuIdx;
    G4int fneuPDG;
    TLorentzVector fneuP4;
    TLorentzVector fneuX4;
    G4int fint_type;
    G4int fscattering_type;
    G4int ffslPDG;
    TLorentzVector ffslP4;
    G4double fW;
};

#endif
