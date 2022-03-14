///  /data/ojalvo/2021_Run3_Trigger/test/CMSSW_11_2_0/src/L1Trigger/Run3Ntuplizer/plugins/BoostedTauStudies.cc

// Get Gen Particles, miniTaus, boostedMiniTaus, boosted jets, regions, trigger towers
//Make gen particle tree - find higgs boson candidate
//    - store pt, eta, phi, mass, of Higgs boson candidate
//    - find matched taus, tau1 and tau2
//    - get tau1/tau2 pt, eta, phi, gen particle type (electron, muon, hadron)
//    - get delta eta, delta phi, deltaR between taus
//    - store gen-level plots of all variables , also 2D plots of tau Delta R vs higgs Pt
//Make reco particle tree - match gen taus, electrons, muons to reco candidates
//    - match gen tau to reco tau, gen electron to reco electron, etc.
//    - match to L1Boosted objects
//    - store 2D array of towers
//    - store 2D array of regions


// system include files
#include <memory>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "L1Trigger/Run3Ntuplizer/plugins/helpers.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace l1extra;
using namespace std;

bool compareByPtGenCand (genP4Cand i, genP4Cand j) {return(i.p4.pt() > j.p4.pt()); };
//
// class declaration
//

class BoostedTauStudies : public edm::EDAnalyzer {
public:
  explicit BoostedTauStudies(const edm::ParameterSet&);
  ~BoostedTauStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void zeroOutAllVariables();

private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es);      
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;


  edm::Service<TFileService> tfs_;  
  TTree* genTree;
  double genPt_1, genEta_1, genPhi_1, genM_1, genLegDR;
  double genPt_2, genEta_2, genPhi_2, genM_2;
  double genTauPt_1, genTauEta_1, genTauPhi_1, genTauM_1, genTauLegDR;
  double genTauPt_2, genTauEta_2, genTauPhi_2, genTauM_2;
  int genPdg_1, genPdg_2;
  double mom_Mass, mom_Pt, mom_Eta, mom_Phi;
  double momVis_Mass, momVis_Pt, momVis_Eta, momVis_Phi;
  int genDM_1, genDM_2;
  double genLegDEta, genLegDPhi;
  int mom_Pdg;
  int run, lumi, event;
  int nTauE, nTauMu, nTauH, nCandsFound;

  double recoPt_1, recoEta_1, recoPhi_1, recoM_1;
  double recoPt_2, recoEta_2, recoPhi_2, recoM_2;
  double recoDiCandPt, recoDiCandM;
  double recoDR, recoDEta, recoDPhi;
  double l1Pt, l1Eta, l1Phi;
  double l1BoostedPt, l1BoostedEta, l1BoostedPhi;

  edm::EDGetTokenT<reco::GenParticleCollection> genSource_;
  edm::EDGetTokenT<vector<pat::Tau> > boostedTauSource_;
  edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<vector<l1extra::L1JetParticle>> l1BoostedSource_;
  bool debug;

};
//vector<pat::Tau>                      "slimmedTausBoosted"        ""                "PAT"
BoostedTauStudies::BoostedTauStudies(const edm::ParameterSet& iConfig) :
  genSource_(        consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>( "genParticles"))),
  boostedTauSource_( consumes<vector<pat::Tau> >              ( edm::InputTag( "slimmedTausBoosted"  ,"","PAT"))),
  stage2JetToken_(   consumes<BXVector<l1t::Jet>>             ( edm::InputTag("caloStage2Digis","Jet","RECO"))),
  l1BoostedSource_(  consumes<vector<l1extra::L1JetParticle> >( edm::InputTag( "uct2016EmulatorDigis","Boosted","")))
{

  genTree = tfs_->make<TTree>("genParticleTree", "Gen Particle Tree ");
  genTree->Branch("run",        &run,     "run/I");
  genTree->Branch("lumi",       &lumi,    "lumi/I");
  genTree->Branch("event",      &event,   "event/I");
  genTree->Branch("mom_Pdg",    &mom_Pdg, "mom_Pdg/I");
  genTree->Branch("mom_Mass",   &mom_Mass,"mom_Mass/D");
  genTree->Branch("mom_Pt",     &mom_Pt,  "mom_Pt/D");
  genTree->Branch("mom_Eta",    &mom_Eta, "mom_Eta/D");
  genTree->Branch("mom_Phi",    &mom_Phi, "mom_Phi/D");
  genTree->Branch("momVis_Mass", &momVis_Mass,"momVis_Mass/D");
  genTree->Branch("momVis_Pt",   &momVis_Pt,  "momVis_Pt/D");
  genTree->Branch("momVis_Eta",  &momVis_Eta, "momVis_Eta/D");
  genTree->Branch("momVis_Phi",  &momVis_Phi, "momVis_Phi/D");

  genTree->Branch("genTauPt_1",     &genTauPt_1,    "genTauPt_1/D");
  genTree->Branch("genTauEta_1",    &genTauEta_1,   "genTauEta_1/D");
  genTree->Branch("genTauPhi_1",    &genTauPhi_1,   "genTauPhi_1/D");
  genTree->Branch("genTauM_1",      &genTauM_1,     "genTauM_1/D");

  genTree->Branch("genTauPt_2",     &genTauPt_2,    "genTauPt_2/D");
  genTree->Branch("genTauEta_2",    &genTauEta_2,   "genTauEta_2/D");
  genTree->Branch("genTauPhi_2",    &genTauPhi_2,   "genTauPhi_2/D");
  genTree->Branch("genTauM_2",      &genTauM_2,     "genTauM_2/D");

  genTree->Branch("genTauLegDR",    &genTauLegDR,   "genTauLegDR/D");

  genTree->Branch("genPt_1",     &genPt_1,    "genPt_1/D");
  genTree->Branch("genEta_1",    &genEta_1,   "genEta_1/D");
  genTree->Branch("genPhi_1",    &genPhi_1,   "genPhi_1/D");
  genTree->Branch("genM_1",      &genM_1,     "genM_1/D");
  genTree->Branch("genPdg_1",    &genPdg_1,   "genPdg_1/D");

  genTree->Branch("genPt_2",     &genPt_2,    "genPt_2/D");
  genTree->Branch("genEta_2",    &genEta_2,   "genEta_2/D");
  genTree->Branch("genPhi_2",    &genPhi_2,   "genPhi_2/D");
  genTree->Branch("genM_2",      &genM_2,     "genM_2/D");
  genTree->Branch("genPdg_2",    &genPdg_2,   "genPdg_2/D");

  genTree->Branch("genLegDR",    &genLegDR,   "genLegDR/D");
  genTree->Branch("genLegDEta",  &genLegDEta, "genLegDEta/D");
  genTree->Branch("genLegDPhi",  &genLegDPhi, "genLegDPhi/D");

  genTree->Branch("nTauE",   &nTauE,   "nTauE/I");
  genTree->Branch("nTauMu",  &nTauMu,  "nTauMu/I");
  genTree->Branch("nTauH",   &nTauH,   "nTauH/I");
  genTree->Branch("nGenCandsFound",   &nCandsFound,   "nCandsFound/I");

  genTree->Branch("recoPt_1",     &recoPt_1,    "recoPt_1/D");
  genTree->Branch("recoEta_1",    &recoEta_1,   "recoEta_1/D");
  genTree->Branch("recoPhi_1",    &recoPhi_1,   "recoPhi_1/D");
  genTree->Branch("recoM_1",      &recoM_1,     "recoM_1/D");

  genTree->Branch("recoPt_2",     &recoPt_2,    "recoPt_2/D");
  genTree->Branch("recoEta_2",    &recoEta_2,   "recoEta_2/D");
  genTree->Branch("recoPhi_2",    &recoPhi_2,   "recoPhi_2/D");
  genTree->Branch("recoM_2",      &recoM_2,     "recoM_2/D");

  genTree->Branch("recoDiCandPt", &recoDiCandPt, "recoDiCandPt/D");
  genTree->Branch("recoDiCandM",  &recoDiCandM,  "recoDiCandM/D");
  genTree->Branch("recoDR",    &recoDR,   "recoDR/D");
  genTree->Branch("recoDEta",  &recoDEta, "recoDEta/D");
  genTree->Branch("recoDPhi",  &recoDPhi, "recoDPhi/D");

  genTree->Branch("l1Pt",     &l1Pt,    "l1Pt/D");
  genTree->Branch("l1Eta",    &l1Eta,   "l1Eta/D");
  genTree->Branch("l1Phi",    &l1Phi,   "l1Phi/D");

  genTree->Branch("l1BoostedPt",     &l1BoostedPt,    "l1BoostedPt/D");
  genTree->Branch("l1BoostedEta",    &l1BoostedEta,   "l1BoostedEta/D");
  genTree->Branch("l1BoostedPhi",    &l1BoostedPhi,   "l1BoostedPhi/D");

}

void BoostedTauStudies::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  //make this into an input
  debug = true;
 //Get Gen Particles, miniTaus, boostedMiniTaus, boosted jets, regions, trigger tower
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

//-------------
//Get Gen Particles, miniTaus, boostedMiniTaus, boosted jets, regions, trigger towers
//Make gen particle tree - find higgs boson candidate
//    - store pt, eta, phi, mass, of Higgs boson 
//    - store pt, eta, phi, mass, of visible Higgs boson 
//    - find matched taus, tau1 and tau2
//    - get tau1/tau2 pt, eta, phi, gen particle type (electron, muon, hadron)
//    - get delta eta, delta phi, deltaR between taus
//    - store gen-level plots of all variables , also 2D plots of tau Delta R vs higgs Pt

  edm::Handle<vector<pat::Tau> > slimmedBoostedTaus;
  edm::Handle<reco::GenParticleCollection> genParticles;
  
  //need a p4vis for the boson
  vector<reco::GenParticle> genTaus;
  vector<reco::GenParticle> genElectrons;
  vector<reco::GenParticle> genMuons;
  vector<genP4Cand> genCandLegs;

  reco::Candidate::LorentzVector momVisP4;
  vector<genVisTau> genVisTaus;
  genVisTaus.clear(); genTaus.clear(); genElectrons.clear(); genMuons.clear(); genCandLegs.clear();
  zeroOutAllVariables();

  // Accessing L1boosted collection reco::Candidate::LorentzVector
  std::vector<reco::Candidate::LorentzVector> l1BoostedJets;
  l1BoostedJets.clear();
  edm::Handle<vector<l1extra::L1JetParticle>> l1Boosted;
  if(!evt.getByToken(l1BoostedSource_, l1Boosted)) cout<<"ERROR GETTING THE L1BOOSTED JETS"<<std::endl;
  evt.getByToken(l1BoostedSource_, l1Boosted);
  const vector<l1extra::L1JetParticle> &l1B = *l1Boosted;
  for(auto obj : l1B) {
    reco::Candidate::LorentzVector temp(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    l1BoostedJets.push_back(temp);
  }

  // Accessing existing L1 seed stored in MINIAOD
  std::vector<reco::Candidate::LorentzVector> l1Jets;
  l1Jets.clear();
  edm::Handle<BXVector<l1t::Jet>> stage2Jets;
  if(!evt.getByToken(stage2JetToken_, stage2Jets)) cout<<"ERROR GETTING THE STAGE 2 JETS"<<std::endl;
  evt.getByToken(stage2JetToken_, stage2Jets);
  const BXVector<l1t::Jet> &s2j = *stage2Jets;
  for(auto obj : s2j) {
    reco::Candidate::LorentzVector temp(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    l1Jets.push_back(temp);
  }

  if(evt.getByToken(genSource_, genParticles)){//Begin Getting Gen Particles
    for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      //Identify the Higgs (or other hard process)
      int genPdgId = genparticle->pdgId();
      if(genPdgId == 25){
	mom_Pdg = 25;
	mom_Mass = genparticle->mass();
	mom_Pt = genparticle->pt();
	mom_Eta = genparticle->eta();
	mom_Phi = genparticle->phi();
	if(debug) std::cout<<"Found Higgs, PdgID = "<< mom_Pdg <<
		    ", M = "<< mom_Mass <<
		    ", Pt = "<< mom_Pt <<
		    ", Eta = "<< mom_Eta <<
		    ", Phi = "<< mom_Phi << std::endl;

      }
      if(abs(genPdgId) == 15&&genparticle->status()==1){
	if(genparticle->pt()>genTauPt_1&&genparticle->pt()>genTauPt_2){
	  genTauM_1   = genparticle->mass();
	  genTauPt_1  = genparticle->pt();
	  genTauEta_1 = genparticle->eta();
	  genTauPhi_1 = genparticle->phi();
	}
	else if(genparticle->pt()<genTauPt_1&&genparticle->pt()>genTauPt_2){
	  genTauM_2   = genparticle->mass();
	  genTauPt_2  = genparticle->pt();
	  genTauEta_2 = genparticle->eta();
	  genTauPhi_2 = genparticle->phi();
	}
      }

      if(genTauPt_1>0&&genTauPt_2>0){
	genTauLegDR = abs(sqrt( ((genTauEta_1-genTauEta_2)*(genTauEta_1-genTauEta_2))+((genTauPhi_1-genTauPhi_2)*(genTauPhi_1-genTauPhi_2))));
      }

      //Check if this is needed
      //if ( genparticle->status() > 21 && genparticle->status() < 41){
        //genMother = genparticle->motherRef(0)->pdgId();
	//check if the object is a decay product of a tau
	if(genparticle->isDirectPromptTauDecayProductFinalState()){
	  //if it is an electron then add it to the electron collection
	  if(abs(genPdgId)==11){ // electron
	    genElectrons.push_back(*genparticle);
	    genP4Cand Temp;
	    Temp.isTauH=false;
	    Temp.pdgId=genPdgId;
	    Temp.decayMode=-10;
	    Temp.p4 = genparticle->p4();
	    genCandLegs.push_back(Temp);
	  }
	  if(abs(genPdgId)==13){ // muon
	    genMuons.push_back(*genparticle);
	    genP4Cand Temp;
	    Temp.isTauH=false;
	    Temp.pdgId=genPdgId;
	    Temp.decayMode=-10;
	    Temp.p4 = genparticle->p4();
	    genCandLegs.push_back(Temp);
	  }
	}

	if(abs(genparticle->pdgId())==15){
	  std::vector<const reco::GenParticle*> stableDaughters;
	  findDaughters(&*genparticle, stableDaughters,1);
	  bool isPromptEorMu = false;
	  for(auto daughter : stableDaughters)
	    if(daughter->isDirectPromptTauDecayProductFinalState()&&(daughter->pdgId()==11||daughter->pdgId()==13)){
	      if(debug)std::cout<<"Found a tau decaying to a muon or electron"<<std::endl;
	      isPromptEorMu = true;
	    }

	  // Find hadronic taus later
	  if(isPromptEorMu==false && genparticle->motherRef(0)->pdgId() ==25)
	    genTaus.push_back(*genparticle);	  
	}
	
	// }
      // now create the boson candidate
      if ( ( genparticle->fromHardProcessFinalState() && (genPdgId==11 || genPdgId==13) ) || ( genparticle->statusFlags().isDirectHardProcessTauDecayProduct() && !(genPdgId==12||genPdgId==14 ||genPdgId==16) )){ 

	if(debug)std::cout<<"building the visible boson"<<std::endl;
	momVisP4 += genparticle->p4();
      }
    }

    momVis_Mass = momVisP4.mass();
    momVis_Pt =   momVisP4.pt();
    momVis_Eta =  momVisP4.eta();
    momVis_Phi =  momVisP4.phi();
    
    // Process the taus and get the visible decay products
    for(auto genTau: genTaus){
      reco::Candidate::LorentzVector visGenTau= getVisMomentum(&genTau, &*genParticles);
      genVisTau Temp;
      int decayMode = GetDecayMode(&genTau);
      Temp.p4 = visGenTau;
      Temp.decayMode = decayMode;
      genVisTaus.push_back(Temp);

      genP4Cand Temp2;
      Temp2.isTauH=true;
      Temp2.pdgId=15;
      Temp2.decayMode=decayMode;
      Temp2.p4 = visGenTau;
      genCandLegs.push_back(Temp2);
    //cout<<"Tau Decay Mode "<<decayMode<<"tau vis pt: "<<genPt<<" genEta: "<<genEta<<" genPhi: "<<genPhi<<std::endl;
    }

    nTauE = genElectrons.size();
    nTauMu = genMuons.size();
    nTauH = genVisTaus.size();
    nCandsFound = genCandLegs.size();

    if(debug){
      std::cout<<"nTauE: "<<nTauE<<", nTauMu: "<<nTauMu<<", nTauH: "<<nTauH<<std::endl;
    }

    std::sort(genCandLegs.begin(),genCandLegs.end(),compareByPtGenCand);
    //Now Fill the tree with leg 1 and leg 2 also include number of legs for debugging... 

    //Leg1
    if(genCandLegs.size()>0){
      genPt_1 = genCandLegs.at(0).p4.pt();
      genEta_1 = genCandLegs.at(0).p4.eta();
      genPhi_1 = genCandLegs.at(0).p4.phi();
      genM_1 = genCandLegs.at(0).p4.mass();
      genPdg_1 = genCandLegs.at(0).pdgId;
      genDM_1 = genCandLegs.at(0).decayMode;
    }

    //Leg2
    if(genCandLegs.size()>1){
      genPt_2 = genCandLegs.at(1).p4.pt();
      genEta_2 = genCandLegs.at(1).p4.eta();
      genPhi_2 = genCandLegs.at(1).p4.phi();
      genM_2 = genCandLegs.at(1).p4.mass();
      genPdg_2 = genCandLegs.at(1).pdgId;
      genDM_2 = genCandLegs.at(1).decayMode;
      
      genLegDR = reco::deltaR(genCandLegs.at(0).p4, genCandLegs.at(1).p4);
      genLegDEta = abs(genCandLegs.at(0).p4.eta() - genCandLegs.at(1).p4.eta());
      genLegDPhi = abs(genCandLegs.at(0).p4.phi() - genCandLegs.at(1).p4.phi());

      // We have a funny case where three cands are somehow being filled. This needs further debugging as there should only be two.
      // To circumnavigate, let's see if the genLegDR == 0 and nCandsFound == 3
      if(genCandLegs.size()>2 && genLegDR < 0.001 && nCandsFound == 3){
	genPt_2 = genCandLegs.at(2).p4.pt();
	genEta_2 = genCandLegs.at(2).p4.eta();
	genPhi_2 = genCandLegs.at(2).p4.phi();
	genM_2 = genCandLegs.at(2).p4.mass();
	genPdg_2 = genCandLegs.at(2).pdgId;
	genDM_2 = genCandLegs.at(2).decayMode;
	
	genLegDR = reco::deltaR(genCandLegs.at(0).p4, genCandLegs.at(2).p4);
	genLegDEta = abs(genCandLegs.at(0).p4.eta() - genCandLegs.at(2).p4.eta());
	genLegDPhi = abs(genCandLegs.at(0).p4.phi() - genCandLegs.at(2).p4.phi());

      }
    }

  //get the boosted L1 objects
  //get boosted Reco Taus
  //match efficiency with the di-tau objects


  //get TPGs, build 4x4 grid with tpgs
  //get regions, buid 3x3 grid with regions
  //get Reco Taus
  //check to see if gen tau is boosted deltaR<1.0 
  //check to see if there are two hadronic tau candidates
  //match boosted tau pair to boosted vis candidate
  //find efficiency of boosted tau pair. 

  // match closest deltaR to gen leg cand 1
  // match closest deltaR to gen leg cand 2
  // count number of matching taus for leg 1 and leg 2
  //int nMatchLeg1=0;
  //int nMatchLeg2=0;
  int nRecoCandsFound = 0;
  //find  best di-tau candidate
  reco::Candidate::LorentzVector candTau1;
  reco::Candidate::LorentzVector candTau2;
  float minimumRecoGenDRMatch = 2;
  float minimumL1GenDRMatch = 2;

  if(evt.getByToken(boostedTauSource_, slimmedBoostedTaus)){//Begin Getting reco taus
    if(genLegDR<2){
      for(unsigned int i = 0; i < slimmedBoostedTaus->size(); i++){
	candTau1 = slimmedBoostedTaus->at(i).p4();
	for(unsigned int j = i+1; j < slimmedBoostedTaus->size(); j++){
	  candTau2 = slimmedBoostedTaus->at(j).p4();
	  reco::Candidate::LorentzVector recoDiTauCand = candTau1 + candTau2;

	  //requiring that the composite cand is closest to the momVisP4 and that the 
	  if(reco::deltaR(recoDiTauCand,momVisP4)<minimumRecoGenDRMatch && reco::deltaR(candTau1,candTau2)<2 ){
	    minimumRecoGenDRMatch = reco::deltaR(recoDiTauCand,momVisP4);	    
	    recoPt_1  = candTau1.pt();
	    recoEta_1 = candTau1.eta();
	    recoPhi_1 = candTau1.phi();
	    recoM_1   = candTau1.M();

	    recoPt_2  = candTau2.pt();
	    recoEta_2 = candTau2.eta();
	    recoPhi_2 = candTau2.phi();
	    recoM_2   = candTau2.M();

	    recoDR = reco::deltaR(candTau1, candTau2);
	    recoDEta = abs(candTau1.eta() - candTau2.eta());
	    recoDPhi = abs(candTau1.phi() - candTau2.phi());
	    recoDiCandPt = recoDiTauCand.pt();
	    recoDiCandM = recoDiTauCand.M();
	    nRecoCandsFound++;
	  }
	}
      }
      for(auto l1Jet : l1Jets){
	if(reco::deltaR(l1Jet,momVisP4)<minimumL1GenDRMatch&&l1Jet.pt()>l1Pt){
	  minimumL1GenDRMatch = reco::deltaR(l1Jet,momVisP4);
	  l1Pt = l1Jet.pt();
	  l1Eta = l1Jet.eta();
	  l1Phi = l1Jet.phi();
	}
      }
      minimumL1GenDRMatch = 2;
      for(auto l1Jet : l1BoostedJets){
	if(reco::deltaR(l1Jet,momVisP4)<minimumL1GenDRMatch&&l1Jet.pt()>l1BoostedPt){
	  if(debug)std::cout<<"Boosted Jet Pt: "<<l1Jet.pt()<<", Eta: "<<l1Jet.eta()<<", Phi: "<<l1Jet.phi()<<std::endl;
	  minimumL1GenDRMatch = reco::deltaR(l1Jet,momVisP4);
	  l1BoostedPt  = l1Jet.pt();
	  l1BoostedEta = l1Jet.eta();
	  l1BoostedPhi = l1Jet.phi();
	}
      }
    }
  }
  genTree->Fill();
  }



}




// ------------ method called once each job just before starting event loop  ------------
void 
BoostedTauStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BoostedTauStudies::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void
BoostedTauStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}
 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  BoostedTauStudies::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  BoostedTauStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  BoostedTauStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BoostedTauStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void BoostedTauStudies::zeroOutAllVariables(){
  genPt_1=-99, genEta_1=-99, genPhi_1=-99, genM_1=-99, genLegDR=-99;
  genPt_2=-99, genEta_2=-99, genPhi_2=-99, genM_2=-99;
  genTauPt_1=-99, genTauEta_1=-99, genTauPhi_1=-99, genTauM_1=-99, genTauLegDR=99;
  genTauPt_2=-99, genTauEta_2=-99, genTauPhi_2=-99, genTauM_2=-99;
  genPdg_1=-99, genPdg_2=-99;
  mom_Mass=-99, mom_Pt=-99, mom_Eta=-99, mom_Phi=-99;
  momVis_Mass=-99, momVis_Pt=-99, momVis_Eta=-99, momVis_Phi=-99;
  mom_Pdg = -99;
  genDM_1 = -99, genDM_2 = -99;
  genLegDR = 99, genLegDEta = 99, genLegDPhi = 99;
  recoPt_1 = 0, recoEta_1 = 0, recoPhi_1 = 0, recoM_1 = 0;
  recoPt_2 = 0, recoEta_2 = 0, recoPhi_2 = 0, recoM_2 = 0;
  recoDiCandPt = 0, recoDiCandM =0;
  recoDR = 99, recoDEta = 99, recoDPhi = 99;
  l1Pt = 0; l1Eta = -99; l1Phi = -99;
  l1BoostedPt = 0; l1BoostedEta = -99; l1BoostedPhi = -99;
}

BoostedTauStudies::~BoostedTauStudies(){

}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedTauStudies);

