#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

class AK4JetNtupleProducer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit AK4JetNtupleProducer(const edm::ParameterSet&);
  ~AK4JetNtupleProducer();
		
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void clearVars();  
  bool isNeutralPdg(int);

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfcand_token_;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genjet_token_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particletable_token_;
  edm::Handle<std::vector<Run3ScoutingParticle>> pfcands_;
  edm::Handle<reco::GenParticleCollection> gencands_;
  edm::Handle<reco::JetFlavourInfoMatchingCollection> genjets_;

  TTree* tree;

  std::vector<Float16_t> pfcand_pt_log_nopuppi;
  std::vector<Float16_t> pfcand_e_log_nopuppi;
  std::vector<Float16_t> pfcand_etarel;
  std::vector<Float16_t> pfcand_phirel;
  std::vector<Float16_t> pfcand_abseta;
  std::vector<Float16_t> pfcand_charge;
  std::vector<Float16_t> pfcand_isEl;
  std::vector<Float16_t> pfcand_isMu;
  std::vector<Float16_t> pfcand_isGamma;
  std::vector<Float16_t> pfcand_isChargedHad;
  std::vector<Float16_t> pfcand_isNeutralHad;
  std::vector<Float16_t> pfcand_lostInnerHits;
  std::vector<Float16_t> pfcand_normchi2;
  std::vector<Float16_t> pfcand_quality;
  std::vector<Float16_t> pfcand_dz;
  std::vector<Float16_t> pfcand_dzsig;
  std::vector<Float16_t> pfcand_dxy;
  std::vector<Float16_t> pfcand_dxysig;
  std::vector<Float16_t> pfcand_btagEtaRel;
  std::vector<Float16_t> pfcand_btagPtRatio;
  std::vector<Float16_t> pfcand_btagPParRatio;

  float j_pt;
  float j_eta;
  float j_phi;
  float j_mass;
  int j_no;
  int j_npfcands;

  int j_nCHadrons;
  int j_nBHadrons;
  int j_partonFlavour;
  int j_hadronFlavour;
  int sample_isQCD;

  int event_no;

  bool isQCD_ = false;

  float dR_= 0.4;
};

AK4JetNtupleProducer::AK4JetNtupleProducer(const edm::ParameterSet& iConfig):
  pfcand_token_(consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("pf_candidates"))),
  genjet_token_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  particletable_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>())
{

  isQCD_ = iConfig.getUntrackedParameter<bool>("isQCD", false);

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("Events", "Events");

  tree->Branch("pfcand_pt_log_nopuppi", &pfcand_pt_log_nopuppi);
  tree->Branch("pfcand_e_log_nopuppi", &pfcand_e_log_nopuppi);
  tree->Branch("pfcand_etarel", &pfcand_etarel);
  tree->Branch("pfcand_phirel", &pfcand_phirel);
  tree->Branch("pfcand_abseta", &pfcand_abseta);
  tree->Branch("pfcand_charge", &pfcand_charge);
  tree->Branch("pfcand_isEl", &pfcand_isEl);
  tree->Branch("pfcand_isMu", &pfcand_isMu);
  tree->Branch("pfcand_isGamma", &pfcand_isGamma);
  tree->Branch("pfcand_isChargedHad", &pfcand_isChargedHad);
  tree->Branch("pfcand_isNeutralHad", &pfcand_isNeutralHad);
  tree->Branch("pfcand_lostInnerHits", &pfcand_lostInnerHits);
  tree->Branch("pfcand_normchi2", &pfcand_normchi2);
  tree->Branch("pfcand_quality", &pfcand_quality);
  tree->Branch("pfcand_dz", &pfcand_dz);
  tree->Branch("pfcand_dzsig", &pfcand_dzsig);
  tree->Branch("pfcand_dxy", &pfcand_dxy);
  tree->Branch("pfcand_dxysig", &pfcand_dxysig);
  tree->Branch("pfcand_btagEtaRel", &pfcand_btagEtaRel);
  tree->Branch("pfcand_btagPtRatio", &pfcand_btagPtRatio);
  tree->Branch("pfcand_btagPParRatio", &pfcand_btagPParRatio);

  tree->Branch("j_pt", &j_pt);
  tree->Branch("j_eta", &j_eta);
  tree->Branch("j_phi", &j_phi);
  tree->Branch("j_mass", &j_mass);
  tree->Branch("j_no", &j_no);
  tree->Branch("j_npfcands", &j_npfcands);

  tree->Branch("j_nCHadrons", &j_nCHadrons);
  tree->Branch("j_nBHadrons", &j_nBHadrons);
  tree->Branch("j_partonFlavour", &j_partonFlavour);
  tree->Branch("j_hadronFlavour", &j_hadronFlavour);
  tree->Branch("sample_isQCD", &sample_isQCD);

  tree->Branch("event_no", &event_no);
}

AK4JetNtupleProducer::~AK4JetNtupleProducer() {
}

bool AK4JetNtupleProducer::isNeutralPdg(int pdgId) {
   const int neutralPdgs_array[] = {9, 21, 22, 23, 25, 12, 14, 16, 111, 130, 310, 311, 421, 511, 2112}; // gluon, gluon, gamma, Z0, higgs, electron neutrino, muon neutrino, tau neutrino, pi0, K0_L, K0_S; K0, neutron
   const std::vector<int> neutralPdgs(neutralPdgs_array, neutralPdgs_array + sizeof(neutralPdgs_array) / sizeof(int));
   if (std::find(neutralPdgs.begin(), neutralPdgs.end(), std::abs(pdgId)) == neutralPdgs.end())
     return false;
 
   return true;
}

void AK4JetNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // PF candidates
  iEvent.getByToken(pfcand_token_, pfcands_);
  // GEN jets
  iEvent.getByToken(genjet_token_, genjets_);
  // Particle table
  auto pdt = iSetup.getHandle(particletable_token_);
  const HepPDT::ParticleDataTable* particletable_ = pdt.product();

  // Create jets
  std::vector<fastjet::PseudoJet> j_part;
  j_part.reserve(pfcands_->size());
  int pfcand_i = 0;
  for (auto pfcands_iter = pfcands_->begin(); pfcands_iter != pfcands_->end(); ++pfcands_iter) {

    auto pfm = particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId())) != nullptr
                        ? particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId()))->mass()
                        : -99.f;
    if (pfm < -90) continue;

    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfm);
    j_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    j_part.back().set_user_index(pfcand_i);
    pfcand_i++;
  }

  fastjet::JetDefinition ak4_def = fastjet::JetDefinition(fastjet::antikt_algorithm, dR_);
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  fastjet::ClusterSequenceArea ak4_cs(j_part, ak4_def, area_def);
  std::vector<fastjet::PseudoJet> ak4_jets = fastjet::sorted_by_pt(ak4_cs.inclusive_jets(170.0));

  // Match jet-gen jet
  std::map<int, reco::JetFlavourInfo> genmatch_resultmap;
  std::vector<int> genmatch_unmatched;
  std::vector<std::tuple<int, int, float> > pairlist;

  int ak4_jet_idx = 0;

  for(unsigned int i=0; i<ak4_jets.size(); i++) {
    bool found_match = false;
    for(unsigned int j=0; j<genjets_->size(); j++) {
      float dR = reco::deltaR(ak4_jets[i].eta(), ak4_jets[i].phi(), (*genjets_)[j].first.get()->eta(), (*genjets_)[j].first.get()->phi());
      if(dR < dR_) {
        pairlist.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    ak4_jet_idx++;
    if(!found_match) {
       genmatch_unmatched.push_back(i);
    }
  }

  std::sort(pairlist.begin(), pairlist.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairlist.size() > 0) {
    reco::JetFlavourInfo genjet_assn = (*genjets_)[std::get<1>(pairlist[0])].second;
    genmatch_resultmap[std::get<0>(pairlist[0])] = genjet_assn;
    for(unsigned int k=1; k<pairlist.size(); k++) {
      if(std::get<0>(pairlist[k]) == std::get<0>(pairlist[0]) ||
         std::get<1>(pairlist[k]) == std::get<1>(pairlist[0])) {
        pairlist.erase(pairlist.begin() + k);
      }
    }
    pairlist.erase(pairlist.begin());
  }
  
  ak4_jet_idx = 0;
  for(auto &j: ak4_jets) {

    float etasign = j.eta() > 0 ? 1 : -1;
    float jet_px = j.pt() * cos(j.phi());
    float jet_py = j.pt() * sin(j.phi());
    float jet_pz = j.pt() * sinh(j.eta());
    math::XYZVector jet_dir_temp(jet_px, jet_py, jet_pz);
    math::XYZVector jet_dir = jet_dir_temp.Unit();
    TVector3 jet_dir3(jet_px, jet_py, jet_pz);

    const std::vector<fastjet::PseudoJet> constituents = j.constituents();
    for (auto &cand : constituents) {
      auto *reco_cand = dynamic_cast<const Run3ScoutingParticle*> (&pfcands_->at(cand.user_index()));
      float trk_px = reco_cand->trk_pt() * cos(reco_cand->trk_phi());
      float trk_py = reco_cand->trk_pt() * sin(reco_cand->trk_phi());
      float trk_pz = reco_cand->trk_pt() * sinh(reco_cand->trk_eta());
      math::XYZVector track_mom(trk_px, trk_py, trk_pz);
      TVector3 track_mom3(trk_px, trk_py, trk_pz);
      double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);

      float reco_cand_p = reco_cand->pt() * cosh(reco_cand->eta());
      auto rcm = particletable_->particle(HepPDT::ParticleID(reco_cand->pdgId())) != nullptr
                        ? particletable_->particle(HepPDT::ParticleID(reco_cand->pdgId()))->mass()
                        : -99.f;
      if (rcm < -90) continue;

      pfcand_e_log_nopuppi.push_back(log(sqrt(reco_cand_p*reco_cand_p + rcm*rcm)));
      pfcand_pt_log_nopuppi.push_back(log(reco_cand->pt()));
      pfcand_etarel.push_back(etasign * (reco_cand->eta() - j.eta()));
      pfcand_phirel.push_back(deltaPhi(reco_cand->phi(), j.phi()));
      pfcand_abseta.push_back(abs(reco_cand->eta()));
      if (isNeutralPdg(reco_cand->pdgId())) {
         pfcand_charge.push_back(0);
      } else {
         pfcand_charge.push_back(abs(reco_cand->pdgId())/reco_cand->pdgId());
      }
      pfcand_isEl.push_back(abs(reco_cand->pdgId()) == 11);
      pfcand_isMu.push_back(abs(reco_cand->pdgId()) == 13);
      pfcand_isGamma.push_back(abs(reco_cand->pdgId()) == 22);
      pfcand_isChargedHad.push_back(abs(reco_cand->pdgId()) == 211);
      pfcand_isNeutralHad.push_back(abs(reco_cand->pdgId()) == 130);
      pfcand_lostInnerHits.push_back(reco_cand->lostInnerHits());
      pfcand_normchi2.push_back(reco_cand->normchi2());
      pfcand_quality.push_back(reco_cand->quality());
      pfcand_dz.push_back(reco_cand->dz());
      pfcand_dzsig.push_back(reco_cand->dzsig());
      pfcand_dxy.push_back(reco_cand->dxy());
      pfcand_dxysig.push_back(reco_cand->dxysig());
      pfcand_btagEtaRel.push_back(reco::btau::etaRel(jet_dir, track_mom));
      pfcand_btagPtRatio.push_back(track_mom3.Perp(jet_dir3) / track_mag);
      pfcand_btagPParRatio.push_back(jet_dir.Dot(track_mom) / track_mag);
    }

    j_pt = j.pt();
    j_eta = j.eta();
    j_phi = j.phi();
    j_mass = j.m();

    j_no = ak4_jets.size();
    j_npfcands = constituents.size();

    sample_isQCD = isQCD_;

    if (std::find(genmatch_unmatched.begin(), genmatch_unmatched.end(), ak4_jet_idx) != genmatch_unmatched.end()) {
      j_nCHadrons = -99;
      j_nBHadrons = -99;
      j_partonFlavour = -99;
      j_hadronFlavour = -99;
    } else {
      reco::JetFlavourInfo flavJet = genmatch_resultmap[ak4_jet_idx];
      j_nCHadrons = flavJet.getcHadrons().size();
      j_nBHadrons = flavJet.getbHadrons().size();
      j_partonFlavour = flavJet.getPartonFlavour();
      j_hadronFlavour = flavJet.getHadronFlavour();
    }

    event_no = iEvent.id().event();

    ak4_jet_idx++;

    tree->Fill();	
    clearVars();
  }
}

void AK4JetNtupleProducer::clearVars(){
  pfcand_pt_log_nopuppi.clear();
  pfcand_e_log_nopuppi.clear();
  pfcand_etarel.clear();
  pfcand_phirel.clear();
  pfcand_abseta.clear();
  pfcand_charge.clear();
  pfcand_isEl.clear();
  pfcand_isMu.clear();
  pfcand_isGamma.clear();
  pfcand_isChargedHad.clear();
  pfcand_isNeutralHad.clear();
  pfcand_lostInnerHits.clear();
  pfcand_normchi2.clear();
  pfcand_quality.clear();
  pfcand_dz.clear();
  pfcand_dzsig.clear();
  pfcand_dxy.clear();
  pfcand_dxysig.clear();
  pfcand_btagEtaRel.clear();
  pfcand_btagPtRatio.clear();
  pfcand_btagPParRatio.clear();
}

void AK4JetNtupleProducer::beginJob() {
}

void AK4JetNtupleProducer::endJob() {
}

void AK4JetNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(AK4JetNtupleProducer);
