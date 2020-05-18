/**
   
*/

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

#include "L1Trigger/L1TCaloLayer1/src/UCTParameters.hh"

#include "L1Trigger/L1TCaloLayer1/src/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTTower.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"

#include "L1Trigger/L1TCaloSummary/src/UCTObject.hh"
#include "L1Trigger/L1TCaloSummary/src/UCTSummaryCard.hh"
#include "L1Trigger/L1TCaloSummary/src/UCTGeometryExtended.hh"

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
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "L1Trigger/L1TCaloLayer1/src/L1UCTCollections.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Muon.h"


#include "L1Trigger/L1TCaloLayer1/src/L1TCaloLayer1FetchLUTs.hh"

#include <iostream>
#include <string>
#include <bitset>
#include <fstream>

using namespace l1tcalo;
using namespace l1extra;
using namespace std;


//
// class declaration
//

class BoostedJetStudies : public edm::EDAnalyzer {
public:
  explicit BoostedJetStudies(const edm::ParameterSet&);
  ~BoostedJetStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es);      
  virtual void beginJob() override;
  virtual void endJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  void print();
  bitset<12> mergeBits(bitset<12> bitsin);
  // ----------member data ---------------------------
  ofstream myfile;
  edm::InputTag genSrc_;
  edm::EDGetTokenT<EcalTrigPrimDigiCollection> ecalTPSource;
  std::string ecalTPSourceLabel;
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hcalTPSource;
  std::string hcalTPSourceLabel;

  edm::EDGetTokenT<vector<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcAK8_;

  std::vector< std::vector< std::vector < uint32_t > > > ecalLUT;
  std::vector< std::vector< std::vector < uint32_t > > > hcalLUT;
  std::vector< std::vector< uint32_t > > hfLUT;

  uint32_t nPumBins;

  std::vector< std::vector< std::vector < uint32_t > > > pumLUT;

  std::vector< UCTTower* > twrList;

  bool useLSB;
  bool useCalib;
  bool useECALLUT;
  bool useHCALLUT;
  bool useHFLUT;

  double caloScaleFactor;

  uint32_t jetSeed;
  uint32_t tauSeed;
  float tauIsolationFactor;
  uint32_t eGammaSeed;
  double eGammaIsolationFactor;

  bool verbose;

  UCTParameters uctParameters;
  UCTLayer1 *layer1;
  UCTSummaryCard *summaryCard;

  TH1F* nEvents;
  TH1F* recoJet_pt;
  TH1F* recoJet_eta;
  TH1F* recoJet_phi;

  TH1F* recoJetAK8_pt;
  TH1F* recoJetAK8_eta;
  TH1F* recoJetAK8_phi;

  TTree* l1Tree;
  int run, lumi, event;

  double genPt, genEta, genPhi;
  double recoPt, recoEta, recoPhi;
  double l1Pt, l1Eta, l1Phi;

  double genPt_1, genEta_1, genPhi_1;
  double recoPt_1, recoEta_1, recoPhi_1;
  double l1Pt_1, l1Eta_1, l1Phi_1;
  
  double genPt_2, genEta_2, genPhi_2;
  double recoPt_2, recoEta_2, recoPhi_2;
  double l1Pt_2, l1Eta_2, l1Phi_2;

  double genDeltaEta, genDeltaPhi, genDeltaR, genMass;
  double recoDeltaEta, recoDeltaPhi, recoDeltaR, recoMass;
  double l1DeltaEta, l1DeltaPhi, l1DeltaR, l1Mass;

  int l1NthJet_1, l1NthJet_2;
  int recoNthJet_1, recoNthJet_2;

  double vbfBDT;
  double recoPt_;

  int nGenJets, nRecoJets, nL1Jets;
  int l1Matched_1, l1Matched_2;
  void createBranches(TTree *tree);
  TTree* efficiencyTree;
  edm::Service<TFileService> tfs_;  

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BoostedJetStudies::BoostedJetStudies(const edm::ParameterSet& iConfig) :
  ecalTPSource(consumes<EcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalToken"))),
  ecalTPSourceLabel(iConfig.getParameter<edm::InputTag>("ecalToken").label()),
  hcalTPSource(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalToken"))),
  hcalTPSourceLabel(iConfig.getParameter<edm::InputTag>("hcalToken").label()),
  ecalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hcalLUT(28, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(256))),
  hfLUT(12, std::vector < uint32_t >(256)),
  nPumBins(iConfig.getParameter<unsigned int>("nPumBins")),
  pumLUT(nPumBins, std::vector< std::vector<uint32_t> >(2, std::vector<uint32_t>(13))),
  useLSB(iConfig.getParameter<bool>("useLSB")),
  useCalib(iConfig.getParameter<bool>("useCalib")),
  useECALLUT(iConfig.getParameter<bool>("useECALLUT")),
  useHCALLUT(iConfig.getParameter<bool>("useHCALLUT")),
  useHFLUT(iConfig.getParameter<bool>("useHFLUT")),
  caloScaleFactor(iConfig.getParameter<double>("caloScaleFactor")),
  jetSeed(iConfig.getParameter<unsigned int>("jetSeed")),
  tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
  tauIsolationFactor(iConfig.getParameter<double>("tauIsolationFactor")),
  eGammaSeed(iConfig.getParameter<unsigned int>("eGammaSeed")),
  eGammaIsolationFactor(iConfig.getParameter<double>("eGammaIsolationFactor")),
  verbose(iConfig.getParameter<bool>("verbose")),
  uctParameters(iConfig.getParameter<double>("activityFraction"), 
		iConfig.getParameter<double>("ecalActivityFraction"), 
		iConfig.getParameter<double>("miscActivityFraction")),
  jetSrc_(    consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsAK8"))),
  genSrc_((        iConfig.getParameter<edm::InputTag>( "genParticles")))
{
  myfile.open ("goodJetEvents.txt");

  std::vector<double> pumLUTData;
  char pumLUTString[10];
  for(uint32_t pumBin = 0; pumBin < nPumBins; pumBin++) {
    for(uint32_t side = 0; side < 2; side++) {
      if(side == 0) sprintf(pumLUTString, "pumLUT%2.2dp", pumBin);
      else sprintf(pumLUTString, "pumLUT%2.2dn", pumBin);
      pumLUTData = iConfig.getParameter<std::vector < double > >(pumLUTString);
      for(uint32_t iEta = 0; iEta < std::max((uint32_t) pumLUTData.size(), MaxUCTRegionsEta); iEta++) {
	pumLUT[pumBin][side][iEta] = (uint32_t) round(pumLUTData[iEta] / caloScaleFactor);
      }
      if(pumLUTData.size() != (MaxUCTRegionsEta))
	std::cerr << "PUM LUT Data size integrity check failed; Expected size = " << MaxUCTRegionsEta
		  << "; Provided size = " << pumLUTData.size()
		  << "; Will use what is provided :(" << std::endl;
    }
  }
  /*
  produces< L1CaloRegionCollection >();
  produces< L1EmParticleCollection >( "Isolated" ) ;
  produces< L1EmParticleCollection >( "NonIsolated" ) ;
  produces< L1JetParticleCollection >( "Central" ) ;
  produces< L1JetParticleCollection >( "Forward" ) ;
  produces< L1JetParticleCollection >( "Tau" ) ;
  produces< L1JetParticleCollection >( "IsoTau" ) ;
  produces< L1EtMissParticleCollection >( "MET" ) ;
  produces< L1EtMissParticleCollection >( "MHT" ) ;
  */
  layer1 = new UCTLayer1(&uctParameters);
  summaryCard = new UCTSummaryCard(layer1, &pumLUT, jetSeed, tauSeed, tauIsolationFactor, eGammaSeed, eGammaIsolationFactor);
  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	vector<UCTTower*> towers = regions[rgn]->getTowers();
	for(uint32_t twr = 0; twr < towers.size(); twr++) {
	  twrList.push_back(towers[twr]);
	}
      }
    }
  }
  // Initialize the Trees

  recoPt_      = iConfig.getParameter<double>("recoPtCut");
  nEvents      = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  recoJet_pt   = tfs_->make<TH1F>( "recoJet_pt" , "p_{t}", 300,  0., 300. );
  recoJet_eta  = tfs_->make<TH1F>( "recoJet_eta"  , "eta", 100,  -3, 3. );
  recoJet_phi  = tfs_->make<TH1F>( "recoJet_phi"  , "phi", 100,  -4, 4. );
  
  recoJetAK8_pt   = tfs_->make<TH1F>( "recoJetAK8_pt" , "p_{t}", 300,  0., 300. );
  recoJetAK8_eta  = tfs_->make<TH1F>( "recoJetAK8_eta"  , "eta", 100,  -3, 3. );
  recoJetAK8_phi  = tfs_->make<TH1F>( "recoJetAK8_phi"  , "phi", 100,  -4, 4. );

  efficiencyTree = tfs_->make<TTree>("efficiencyTree", "Gen Matched Jet Tree ");
  createBranches(efficiencyTree);
  
}

BoostedJetStudies::~BoostedJetStudies() {
  myfile.close();
  if(layer1 != 0) delete layer1;
  if(summaryCard != 0) delete summaryCard;
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void BoostedJetStudies::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  using namespace edm;

   nEvents->Fill(1);
   run = evt.id().run();
   lumi = evt.id().luminosityBlock();
   event = evt.id().event();
   Handle<L1CaloRegionCollection> regions;
   
   std::vector<pat::Jet> goodJets;
   std::vector<pat::Jet> goodJetsAK8;


  // Start Running Layer 1
  edm::Handle<EcalTrigPrimDigiCollection> ecalTPs;
  evt.getByToken(ecalTPSource, ecalTPs);
  edm::Handle<HcalTrigPrimDigiCollection> hcalTPs;
  evt.getByToken(hcalTPSource, hcalTPs);

  std::unique_ptr<L1CaloRegionCollection> rgnCollection (new L1CaloRegionCollection);
  std::unique_ptr<L1EmParticleCollection> iEGCands(new L1EmParticleCollection);
  std::unique_ptr<L1EmParticleCollection> nEGCands(new L1EmParticleCollection);
  std::unique_ptr<L1JetParticleCollection> iTauCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> nTauCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> cJetCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> fJetCands(new L1JetParticleCollection);
  std::unique_ptr<L1JetParticleCollection> bJetCands(new L1JetParticleCollection);
  std::unique_ptr<L1EtMissParticleCollection> metCands(new L1EtMissParticleCollection);
  std::unique_ptr<L1EtMissParticleCollection> mhtCands(new L1EtMissParticleCollection);

  uint32_t expectedTotalET = 0;

  if(!layer1->clearEvent()) {
    std::cerr << "UCT: Failed to clear event" << std::endl;
    exit(1);
  }

  for ( const auto& ecalTp : *ecalTPs ) {
    int caloEta = ecalTp.id().ieta();
    int caloPhi = ecalTp.id().iphi();
    int et = ecalTp.compressedEt();
    bool fgVeto = ecalTp.fineGrain();
    if(et != 0) {
      UCTTowerIndex t = UCTTowerIndex(caloEta, caloPhi);
      if(!layer1->setECALData(t,fgVeto,et)) {
	std::cerr << "UCT: Failed loading an ECAL tower" << std::endl;
	return;
      }
      expectedTotalET += et;
    }
  }

  for ( const auto& hcalTp : *hcalTPs ) {
    int caloEta = hcalTp.id().ieta();
    uint32_t absCaloEta = abs(caloEta);
    // Tower 29 is not used by Layer-1
    if(absCaloEta == 29) {
      continue;
    }
    // Prevent usage of HF TPs with Layer-1 emulator if HCAL TPs are old style
    else if(hcalTp.id().version() == 0 && absCaloEta > 29) {
      continue;
    }
    else if(absCaloEta <= 41) {
      int caloPhi = hcalTp.id().iphi();
      if(caloPhi <= 72) {
	int et = hcalTp.SOI_compressedEt();
	bool fg = hcalTp.SOI_fineGrain();
	if(et != 0) {
	  UCTTowerIndex t = UCTTowerIndex(caloEta, caloPhi);
	  uint32_t featureBits = 0;
	  if(fg) featureBits = 0x1F; // Set all five feature bits for the moment - they are not defined in HW / FW yet!
	  if(!layer1->setHCALData(t, featureBits, et)) {
	    std::cerr << "caloEta = " << caloEta << "; caloPhi =" << caloPhi << std::endl;
	    std::cerr << "UCT: Failed loading an HCAL tower" << std::endl;
	    return;
	    
	  }
	  expectedTotalET += et;
	}
      }
      else {
	std::cerr << "Illegal Tower: caloEta = " << caloEta << "; caloPhi =" << caloPhi << std::endl;	
      }
    }
    else {
      std::cerr << "Illegal Tower: caloEta = " << caloEta << std::endl;
    }
  }

  if(!layer1->process()) {
    std::cerr << "UCT: Failed to process layer 1" << std::endl;
    exit(1);
  }

  // Crude check if total ET is approximately OK!
  // We can't expect exact match as there is region level saturation to 10-bits
  // 1% is good enough
  int diff = abs((int)layer1->et() - (int)expectedTotalET);
  if(verbose && diff > 0.01 * expectedTotalET ) {
    print();
    std::cout << "Expected " 
	      << std::showbase << std::internal << std::setfill('0') << std::setw(10) << std::hex
	      << expectedTotalET << std::dec << std::endl;
  }
 
  UCTGeometry g;

  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	uint32_t rawData = regions[rgn]->rawData();
	uint32_t regionData = rawData & 0x0000FFFF;
	// uint32_t regionLoc = rawData >> LocationShift;
	// uint32_t regionET = rawData & 0x3FF;
	uint32_t crate = regions[rgn]->getCrate();
	uint32_t card = regions[rgn]->getCard();
	uint32_t region = regions[rgn]->getRegion();
	bool negativeEta = regions[rgn]->isNegativeEta();
	uint32_t rPhi = g.getUCTRegionPhiIndex(crate, card);
	// We want to reuse L1CaloRegion and L1CaloRegionDetID
	// We do not want to change those classes too much
	// We want comparison to legacy for Barrel and Endcap to work transparently
	// Noting that rEta is packed in 5 bits of L1CaloRegionDetID, we have a scheme!
	// We store the Barrel and Endcap regions in the same location as done for RCT
	// HF has changed in the upgrade, 6x2 HF regions instead of 4x2 in case of RCT
	// Note that for upgrade region numbers range 0-6 for Barrel/Endcap and 7-12 for HF
	// So, the scheme used for rEta for upgrade is:
	// rEta= 0- 3 for -HF regions 7-10
	// rEta= 4-10 for -B/E regions 0-6
	// rEta=11-17 for +B/E regions 0-6
	// rEta=18-23 for +HF regions 7-12
	// rEta=30 for -HF region 11
	// rEta=31 for -HF region 12
	uint32_t rEta = 10 - region;
	if(negativeEta && region == 11) rEta = 30;
	if(negativeEta && region == 12) rEta = 31;
	if(!negativeEta) rEta = 11 + region; // Positive eta portion is offset by 11
	rgnCollection->push_back(L1CaloRegion((uint16_t) regionData, (unsigned) rEta, (unsigned) rPhi, (int16_t) 0));
      }
    }
  }  
  //evt.put(std::move(rgnCollection), "");

  if(!summaryCard->process()) {
    std::cerr << "UCT: Failed to process summary card" << std::endl;
    exit(1);      
  }

  double pt = 0;
  double eta = -999.;
  double phi = -999.;
  double mass = 0;
  double caloScaleFactor = 0.5;
  /* Do not bother with all of these things...  
  std::list<UCTObject*> emObjs = summaryCard->getEMObjs();
  for(std::list<UCTObject*>::const_iterator i = emObjs.begin(); i != emObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    nEGCands->push_back(L1EmParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EmParticle::kNonIsolated));
  }
  std::list<UCTObject*> isoEMObjs = summaryCard->getIsoEMObjs();
  for(std::list<UCTObject*>::const_iterator i = isoEMObjs.begin(); i != isoEMObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    iEGCands->push_back(L1EmParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EmParticle::kIsolated));
  }
  std::list<UCTObject*> tauObjs = summaryCard->getTauObjs();
  for(std::list<UCTObject*>::const_iterator i = tauObjs.begin(); i != tauObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    nTauCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kTau));
  }
  std::list<UCTObject*> isoTauObjs = summaryCard->getIsoTauObjs();
  for(std::list<UCTObject*>::const_iterator i = isoTauObjs.begin(); i != isoTauObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    iTauCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kTau));
  }
  std::list<UCTObject*> centralJetObjs = summaryCard->getCentralJetObjs();
  for(std::list<UCTObject*>::const_iterator i = centralJetObjs.begin(); i != centralJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    cJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kCentral));
    if(pt > 150.) {
      std::cout << "Jet: pt = " << pt << " eta = " << eta << " phi = " << phi << std::endl;
    }
  }
  std::list<UCTObject*> forwardJetObjs = summaryCard->getForwardJetObjs();
  for(std::list<UCTObject*>::const_iterator i = forwardJetObjs.begin(); i != forwardJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());
    fJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kForward));
    }
  */


  Handle<vector<pat::Jet> > jetsAK8;

  if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting Reco Jets

    for (const pat::Jet &jet : *jetsAK8) {
      //recoJetAK8_pt->Fill( jetAK8.pt() );
      //recoJetAK8_eta->Fill( jetAK8.eta() );
      //recoJetAK8_phi->Fill( jetAK8.phi() );
      //get rid of the cruft for analysis to save disk space
      if(jet.pt() < recoPt_ ) 
	continue;

      if(jet.jetFlavourInfo().getbHadrons().size()<2)
	continue;

      if(jet.subjets("SoftDropPuppi").size()!=2)
	continue;
      
      std::cout<<"Found a jet with two subjets and two bhadrons"<<std::endl;

      myfile << run << ":" <<lumi <<":"<<event<<std::endl;

      goodJetsAK8.push_back(jet);
      
    }
  }
  else
    cout<<"Error getting AK8 jets"<<std::endl;


    //std::cout<<"AK8 jets size: "<<jetsAK8->size()<<std::endl;

  bitset<12> goodPattern5(string("001110111110"));

  std::list<UCTObject*> boostedJetObjs = summaryCard->getBoostedJetObjs();
  for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    const UCTObject* object = *i;
    pt = ((double) object->et()) * caloScaleFactor;
    eta = g.getUCTTowerEta(object->iEta());
    phi = g.getUCTTowerPhi(object->iPhi());

    std::cout<<"printing the tower ET:"<<std::endl;
    for(uint32_t iPhi = 0; iPhi < 12; iPhi++){
      for(uint32_t iEta = 0; iEta < 12; iEta++) {
	std::cout<< object->boostedJetTowers()[iEta+iPhi*12]<<setw(20)<<" ";
      }
      std::cout<<std::endl;
    }

    for(auto recoJet : goodJetsAK8){
      if(reco::deltaR(eta,phi, recoJet.eta(),recoJet.phi())<0.4){

	std::string etastring = object->activeTowerEta().to_string<char,std::string::traits_type,std::string::allocator_type>();
	std::string phistring = object->activeTowerPhi().to_string<char,std::string::traits_type,std::string::allocator_type>();
	
	std::cout<<"activeTowerEta: "<< etastring <<std::endl;

	bitset<12> bits_eta = object->activeTowerEta();
	bitset<12> bits_phi = object->activeTowerPhi();
	bitset<12> merged_bits_eta = mergeBits(bits_eta);

	std::cout<<"merged bits_eta: "<<merged_bits_eta.to_string<char,std::string::traits_type,std::string::allocator_type>()<<std::endl;
	std::cout<<std::endl;
	/*
	std::cout<<"activeTowerPhi: "<< phistring<<std::endl;

	for(unsigned int i = 0; i < 11; i++){
	  bitset<12> mask;
	  mask[i] = 1;
	  if(bits_phi[i+1]&bits_phi[i]==1)
	    merged_bits_phi[i] = bits_phi>>1;
	}
	std::cout<<"merged bits_phi: "<<bits_phi.to_string<char,std::string::traits_type,std::string::allocator_type>()<<std::endl;
	*/
	if(object->activeTowerEta()==goodPattern5||object->activeTowerPhi()==goodPattern5 ){
	  std::cout<<"Found subjets!!"<<std::endl;
	}

      }
    }

    bJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kCentral));// using kCentral for now, need a new type
    
  }

  /*
  const UCTObject* et = summaryCard->getET();
  pt = ((double) et->et()) * caloScaleFactor;
  double totET = pt;
  const UCTObject* met = summaryCard->getMET();
  pt = ((double) met->et()) * caloScaleFactor;
  eta = g.getUCTTowerEta(met->iEta());
  phi = g.getUCTTowerPhi(met->iPhi());
  metCands->push_back(L1EtMissParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EtMissParticle::kMET, totET));
  */
  
  /*
  const UCTObject* ht = summaryCard->getHT();
  pt = ((double) ht->et()) * caloScaleFactor;
  double totHT = pt;
  const UCTObject* mht = summaryCard->getMHT();
  pt = ((double) mht->et()) * caloScaleFactor;
  eta = g.getUCTTowerEta(mht->iEta());
  phi = g.getUCTTowerPhi(mht->iPhi());
  mhtCands->push_back(L1EtMissParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1EtMissParticle::kMHT, totHT));
  */

  /*
  evt.put(std::move(iEGCands), "Isolated");
  evt.put(std::move(nEGCands), "NonIsolated");
  evt.put(std::move(iTauCands), "IsoTau");
  evt.put(std::move(nTauCands), "Tau");
  evt.put(std::move(cJetCands), "Central");
  evt.put(std::move(fJetCands), "Forward");
  evt.put(std::move(bJetCands), "Boosted");
  evt.put(std::move(metCands), "MET");
  evt.put(std::move(mhtCands), "MHT");
  */

  // Finish Running Layer 1

  // Start Runing Analysis
  Handle<vector<pat::Jet> > jets;
  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Reco Taus
    for (const pat::Jet &jet : *jets) {
      recoJet_pt->Fill( jet.pt() );
      recoJet_eta->Fill( jet.eta() );
      recoJet_phi->Fill( jet.phi() );
      //get rid of the low pt stuff for analysis to save disk space
      if(jet.pt() > recoPt_ ) {
	goodJets.push_back(jet);

      }
    }
  }
  else
    cout<<"Error getting reco jets"<<std::endl;


  std::list<UCTObject*> l1JetsSorted;

  /*
  for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    const UCTObject* object = *i;

    //for(int j = 0; j < l1JetsSorted.size(); j++){
    for(std::vector<UCTObject>::const_iterator j =  l1JetsSorted.begin(); j != l1JetsSorted.end(); j++){
      //if(object->et() < j->et() ){
      //l1JetsSorted.insert(j,*object);
      //break;
      //}
      //l1JetsSorted.push_back(*object);
      
    }

  }
  */
  int k = 0;
  for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    const UCTObject* object = *i;

    std::cout<< k<<" boosted jet pt: "<< object->et()<<std::endl;
    k++;
    /*
    for(std::list<UCTObject*>::iterator j =  l1JetsSorted.begin(); j != l1JetsSorted.end(); j++){
    const UCTObject* object2 = *j;
      if(object->et() < object2->et() ){
	//l1JetsSorted.insert(j,*object);
	break;
      }
      // if the jet doesn't belong in the middle of the vector then put it at the end
      if(j!=l1JetsSorted.end() && (next(j))==l1JetsSorted.end()){
	l1JetsSorted.push_back(object);
      }
    }


    }
    */
  }
  //std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPtJets);

  //boostedJetObjs sort these guys by pt



  for(auto jet:goodJetsAK8){
    std::cout<<"N BHadrons: "<< jet.jetFlavourInfo().getbHadrons().size()<<std::endl;

    if(jet.jetFlavourInfo().getbHadrons().size()<2)
      continue;

    if(jet.subjets("SoftDropPuppi").size()!=2)
      continue;

    std::cout<<"Found a jet with two subjets and two bhadrons"<<std::endl;

    int i = 0;
    l1extra::L1JetParticle l1Jet;
    //for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
    //const UCTObject test = *i;
    /*
    if(l1JetsSorted.size() > 0){
      for(auto l1jet : l1JetsSorted){
        TLorentzVector temp;
        temp.SetPtEtaPhiE(l1jet.pt(),l1jet.eta(),l1jet.phi(),l1jet.et());
        l1Jets->push_back(temp);
        if(reco::deltaR(jet, recoJet_1)<0.4 && foundL1Jet_1 == 0 ){
          l1Pt_1  = l1jet.pt();
          l1Eta_1 = l1jet.eta();
          l1Phi_1 = l1jet.phi();
          l1NthJet_1 = i;
          l1NTau_1 = nTausInfo[i];
          foundL1Jet_1 = 1;
	  std::cout<<"bitset pattern Eta: "<< jet.activeTowerEta()<<std::endl;
	  std::cout<<"bitset pattern Phi: "<< jet.activeTowerPhi()<<std::endl;
	  
        }
        i++;
      }
      }*/
  }


    //take more variables from here: https://github.com/gouskos/HiggsToBBNtupleProducerTool/blob/opendata_80X/NtupleAK8/src/FatJetInfoFiller.cc#L215-L217
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools


  //Match to boosted jets and see if we can match subjettiness functions...
  //for(std::list<UCTObject*>::const_iterator i = boostedJetObjs.begin(); i != boostedJetObjs.end(); i++) {
  //}  
}

void BoostedJetStudies::print() {
  vector<UCTCrate*> crates = layer1->getCrates();
  for(uint32_t crt = 0; crt < crates.size(); crt++) {
    vector<UCTCard*> cards = crates[crt]->getCards();
    for(uint32_t crd = 0; crd < cards.size(); crd++) {
      vector<UCTRegion*> regions = cards[crd]->getRegions();
      for(uint32_t rgn = 0; rgn < regions.size(); rgn++) {
	if(regions[rgn]->et() > 10) {
	  int hitEta = regions[rgn]->hitCaloEta();
	  int hitPhi = regions[rgn]->hitCaloPhi();
	  vector<UCTTower*> towers = regions[rgn]->getTowers();
	  for(uint32_t twr = 0; twr < towers.size(); twr++) {
	    if(towers[twr]->caloPhi() == hitPhi && towers[twr]->caloEta() == hitEta) {
	      std::cout << "*";
	    }
	    if(towers[twr]->et() > 10) std::cout << *towers[twr];
	  }
	  std::cout << *regions[rgn];
	}
      }
      std::cout << *cards[crd];
    }
    std::cout << *crates[crt];
  }
  std::cout << *layer1;
}


/* The goal of this function is to take a bitwise pattern and merge it such that 
 * there are no patterns that look like 11 or 00. 
 * Example: 010010  -> 001010
 *          111111  -> 000001
 *          111110  -> 000010
 *          000001  -> 000001
 *          001100  -> 000010
 *          110011  -> 000101
 *          0011001100 -> 0000001010
 */
bitset<12> BoostedJetStudies::mergeBits(bitset<12> bitsin){
  bitset<12> finalbits = 0x0;
  // initialize bitmasks
  bitset<12> upper_mask = 0x1FFE;
  bitset<12> lower_mask = 0x0; //default all 0
  int i = 0;
  bitset<12> one = (0x1);
  // compare position j to j+1 so only loop at most 11 times.
  for(int j = 0; j<11; j++){
    std::cout<<"j: "<<j<<std::endl;
    if(bitsin[i+1]==bitsin[i]){
      bitset<12> x = bitsin & upper_mask;
      std::cout<<" bitsin & upper_mask = "<<bitsin.to_string<char,std::string::traits_type,std::string::allocator_type>()<<"&"<<
 	                                upper_mask.to_string<char,std::string::traits_type,std::string::allocator_type>()<<" = "<<
	                                         x.to_string<char,std::string::traits_type,std::string::allocator_type>()<< std::endl;
      x = x>>1;

      std::cout<<"x = x>>1: "<<x.to_string<char,std::string::traits_type,std::string::allocator_type>()<< std::endl;

      std::cout<<"finalbits = (finalbits & lower_mask)|x = ("<<finalbits.to_string<char,std::string::traits_type,std::string::allocator_type>();
      finalbits = (finalbits & lower_mask)|x;


      std::cout<<" & "<< lower_mask.to_string<char,std::string::traits_type,std::string::allocator_type>()
	       <<")|"<<           x.to_string<char,std::string::traits_type,std::string::allocator_type>()
	       <<" = "<<  finalbits.to_string<char,std::string::traits_type,std::string::allocator_type>() <<std::endl;

      //iterate upper mask
      upper_mask = upper_mask << 1;      
      std::cout<<"upper_mask << 1 = "<<upper_mask.to_string<char,std::string::traits_type,std::string::allocator_type>()<<std::endl;
      //iterate lower_mask (add a 1)
      lower_mask = (lower_mask << 1) | one;      
      std::cout<<"upper_mask << 1 = "<<lower_mask.to_string<char,std::string::traits_type,std::string::allocator_type>()<<std::endl;
    }
    else{
      //move on to next position
      i++;
    }
    
  }
  
  return finalbits;
}

void BoostedJetStudies::createBranches(TTree *tree){
    tree->Branch("run",     &run,     "run/I");
    tree->Branch("lumi",    &lumi,     "lumi/I");
    tree->Branch("event",   &event,    "event/I");

    tree->Branch("recoPt_1",      &recoPt_1,     "recoPt_1/D");
    tree->Branch("recoEta_1",     &recoEta_1,    "recoEta_1/D");
    tree->Branch("recoPhi_1",     &recoPhi_1,    "recoPhi_1/D");
    tree->Branch("recoNthJet_1",  &recoNthJet_1, "recoNthJet_1/I");

    tree->Branch("recoPt_2",      &recoPt_2,      "recoPt_2/D");
    tree->Branch("recoEta_2",     &recoEta_2,     "recoEta_2/D");
    tree->Branch("recoPhi_2",     &recoPhi_2,     "recoPhi_2/D");
    tree->Branch("recoNthJet_2",  &recoNthJet_2,  "recoNthJet_2/I");

    tree->Branch("recoDeltaEta",  &recoDeltaEta, "recoDeltaEta/D");
    tree->Branch("recoDeltaPhi",  &recoDeltaPhi, "recoDeltaPhi/D");
    tree->Branch("recoDeltaR",    &recoDeltaR,   "recoDeltaR/D");
    tree->Branch("recoMass",      &recoMass,     "recoMass/D");
      
    tree->Branch("l1Pt_1",        &l1Pt_1,       "l1Pt_1/D"); 
    tree->Branch("l1Eta_1",       &l1Eta_1,      "l1Eta_1/D");
    tree->Branch("l1Phi_1",       &l1Phi_1,      "l1Phi_1/D");
    tree->Branch("l1NthJet_1",    &l1NthJet_1,   "l1NthJet_1/I");

    tree->Branch("l1Pt_2",        &l1Pt_2,       "l1Pt_2/D"); 
    tree->Branch("l1Eta_2",       &l1Eta_2,      "l1Eta_2/D");
    tree->Branch("l1Phi_2",       &l1Phi_2,      "l1Phi_2/D");
    tree->Branch("l1NthJet_2",    &l1NthJet_2,   "l1NthJet_2/I");

    tree->Branch("l1DeltaEta",    &l1DeltaEta,   "l1DeltaEta/D");
    tree->Branch("l1DeltaPhi",    &l1DeltaPhi,   "l1DeltaPhi/D");
    tree->Branch("l1DeltaR",      &l1DeltaR,     "l1DeltaR/D");
    tree->Branch("l1Mass",        &l1Mass,       "l1Mass/D");

    tree->Branch("l1Matched_1",   &l1Matched_1, "l1Matched_1/I");
    tree->Branch("l1Matched_2",   &l1Matched_2, "l1Matched_2/I");
    tree->Branch("nRecoJets",     &nRecoJets,    "nRecoJets/I");
    tree->Branch("nL1Jets",       &nL1Jets,      "nL1Jets/I");
    tree->Branch("vbfBDT",        &vbfBDT,       "vbfBDT/D");
  }


// ------------ method called once each job just before starting event loop  ------------
void 
BoostedJetStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BoostedJetStudies::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void
BoostedJetStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  if(!L1TCaloLayer1FetchLUTs(iSetup, ecalLUT, hcalLUT, hfLUT, useLSB, useCalib, useECALLUT, useHCALLUT, useHFLUT)) {
    std::cerr << "L1TCaloLayer1::beginRun: failed to fetch LUTS - using unity" << std::endl;
  }
  for(uint32_t twr = 0; twr < twrList.size(); twr++) {
    twrList[twr]->setECALLUT(&ecalLUT);
    twrList[twr]->setHCALLUT(&hcalLUT);
    twrList[twr]->setHFLUT(&hfLUT);
  }
}
 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  BoostedJetStudies::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  BoostedJetStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  BoostedJetStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BoostedJetStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedJetStudies);
