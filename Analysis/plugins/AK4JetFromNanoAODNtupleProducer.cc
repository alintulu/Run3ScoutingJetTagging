#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

class AK4JetFromNanoAODNtupleProducer :  public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit AK4JetFromNanoAODNtupleProducer(const edm::ParameterSet &);
  ~AK4JetFromNanoAODNtupleProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  std::map<int, reco::JetFlavourInfo> genmatch_resultmap;
  std::vector<int> genmatch_unmatched;

private:
  typedef edm::View<reco::Candidate> CandidateView;

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void clearVars();
  virtual void genJetMatching(edm::View<reco::Jet>, const reco::JetFlavourInfoMatchingCollection&);
 
  edm::Service<TFileService> fs;

  edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> genjets_token_;
  edm::Handle<CandidateView> pfcands_;
  edm::Handle<reco::JetFlavourInfoMatchingCollection> genjets_;
  
  edm::EDGetTokenT<edm::ValueMap<float>> normchi2_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dz_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dxy_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dzsig_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dxysig_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> lostInnerHits_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> quality_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkPt_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkEta_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkPhi_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> msd_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<float>> n2b1_value_map_token_;

  edm::Handle<edm::ValueMap<float>> normchi2_value_map_;
  edm::Handle<edm::ValueMap<float>> dz_value_map_;
  edm::Handle<edm::ValueMap<float>> dxy_value_map_;
  edm::Handle<edm::ValueMap<float>> dzsig_value_map_;
  edm::Handle<edm::ValueMap<float>> dxysig_value_map_;
  edm::Handle<edm::ValueMap<int>> lostInnerHits_value_map_;
  edm::Handle<edm::ValueMap<int>> quality_value_map_;
  edm::Handle<edm::ValueMap<float>> trkPt_value_map_;
  edm::Handle<edm::ValueMap<float>> trkEta_value_map_;
  edm::Handle<edm::ValueMap<float>> trkPhi_value_map_;
  edm::Handle<edm::ValueMap<float>> msd_value_map_;
  edm::Handle<edm::ValueMap<float>> n2b1_value_map_;
  
  bool isQCD_ = false;

private:
  TTree* tree;

  std::vector<Float16_t> pfcand_pt_log_nopuppi;
  std::vector<Float16_t> pfcand_pt_log;
  std::vector<Float16_t> pfcand_e_log_nopuppi;
  std::vector<Float16_t> pfcand_etarel;
  std::vector<Float16_t> pfcand_phirel;
  std::vector<Float16_t> pfcand_erel;
  std::vector<Float16_t> pfcand_erel_log;
  std::vector<Float16_t> pfcand_ptrel;
  std::vector<Float16_t> pfcand_ptrel_log;
  std::vector<Float16_t> pfcand_abseta;
  std::vector<Float16_t> pfcand_deltaR;
  std::vector<Float16_t> pfcand_mask;
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
};

AK4JetFromNanoAODNtupleProducer::AK4JetFromNanoAODNtupleProducer(const edm::ParameterSet &iConfig)
    : jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      genjets_token_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("gen_jets"))) {

  isQCD_ = iConfig.getUntrackedParameter<bool>("isQCD", false);

  const auto &normchi2_value_map_tag = iConfig.getParameter<edm::InputTag>("normchi2_value_map");
  if (!normchi2_value_map_tag.label().empty()) {
    normchi2_value_map_token_ = consumes<edm::ValueMap<float>>(normchi2_value_map_tag);
  }

  const auto &dz_value_map_tag = iConfig.getParameter<edm::InputTag>("dz_value_map");
  if (!dz_value_map_tag.label().empty()) {
    dz_value_map_token_ = consumes<edm::ValueMap<float>>(dz_value_map_tag);
  }

  const auto &dxy_value_map_tag = iConfig.getParameter<edm::InputTag>("dxy_value_map");
  if (!dxy_value_map_tag.label().empty()) {
    dxy_value_map_token_ = consumes<edm::ValueMap<float>>(dxy_value_map_tag);
  }

  const auto &dzsig_value_map_tag = iConfig.getParameter<edm::InputTag>("dzsig_value_map");
  if (!dzsig_value_map_tag.label().empty()) {
    dzsig_value_map_token_ = consumes<edm::ValueMap<float>>(dzsig_value_map_tag);
  }

  const auto &dxysig_value_map_tag = iConfig.getParameter<edm::InputTag>("dxysig_value_map");
  if (!dxysig_value_map_tag.label().empty()) {
    dxysig_value_map_token_ = consumes<edm::ValueMap<float>>(dxysig_value_map_tag);
  }

  const auto &lostInnerHits_value_map_tag = iConfig.getParameter<edm::InputTag>("lostInnerHits_value_map");
  if (!lostInnerHits_value_map_tag.label().empty()) {
    lostInnerHits_value_map_token_ = consumes<edm::ValueMap<int>>(lostInnerHits_value_map_tag);
  }

  const auto &quality_value_map_tag = iConfig.getParameter<edm::InputTag>("quality_value_map");
  if (!quality_value_map_tag.label().empty()) {
    quality_value_map_token_ = consumes<edm::ValueMap<int>>(quality_value_map_tag);
  }

  const auto &trkPt_value_map_tag = iConfig.getParameter<edm::InputTag>("trkPt_value_map");
  if (!trkPt_value_map_tag.label().empty()) {
    trkPt_value_map_token_ = consumes<edm::ValueMap<float>>(trkPt_value_map_tag);
  }

  const auto &trkEta_value_map_tag = iConfig.getParameter<edm::InputTag>("trkEta_value_map");
  if (!trkEta_value_map_tag.label().empty()) {
    trkEta_value_map_token_ = consumes<edm::ValueMap<float>>(trkEta_value_map_tag);
  }

  const auto &trkPhi_value_map_tag = iConfig.getParameter<edm::InputTag>("trkPhi_value_map");
  if (!trkPhi_value_map_tag.label().empty()) {
    trkPhi_value_map_token_ = consumes<edm::ValueMap<float>>(trkPhi_value_map_tag);
  }

  const auto &msd_value_map_tag = iConfig.getParameter<edm::InputTag>("msd_value_map");
  if (!msd_value_map_tag.label().empty()) {
    msd_value_map_token_ = consumes<edm::ValueMap<float>>(msd_value_map_tag);
  }

  const auto &n2b1_value_map_tag = iConfig.getParameter<edm::InputTag>("n2b1_value_map");
  if (!n2b1_value_map_tag.label().empty()) {
    n2b1_value_map_token_ = consumes<edm::ValueMap<float>>(n2b1_value_map_tag);
  }

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("Events", "Events");

  tree->Branch("pfcand_pt_log_nopuppi", &pfcand_pt_log_nopuppi);
  tree->Branch("pfcand_pt_log", &pfcand_pt_log);
  tree->Branch("pfcand_e_log_nopuppi", &pfcand_e_log_nopuppi);
  tree->Branch("pfcand_etarel", &pfcand_etarel);
  tree->Branch("pfcand_phirel", &pfcand_phirel);
  tree->Branch("pfcand_erel", &pfcand_erel);
  tree->Branch("pfcand_erel_log", &pfcand_erel_log);
  tree->Branch("pfcand_ptrel", &pfcand_ptrel);
  tree->Branch("pfcand_ptrel_log", &pfcand_ptrel_log);
  tree->Branch("pfcand_abseta", &pfcand_abseta);
  tree->Branch("pfcand_deltaR", &pfcand_deltaR);
  tree->Branch("pfcand_mask", &pfcand_mask);
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

AK4JetFromNanoAODNtupleProducer::~AK4JetFromNanoAODNtupleProducer() {}

void AK4JetFromNanoAODNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void AK4JetFromNanoAODNtupleProducer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // Input jets
  auto jets = iEvent.getHandle(jet_token_);
  // PF candidates
  iEvent.getByToken(pfcand_token_, pfcands_);
  // GEN jets with jet flavour
  iEvent.getByToken(genjets_token_, genjets_);
  
  iEvent.getByToken(normchi2_value_map_token_, normchi2_value_map_);
  iEvent.getByToken(dz_value_map_token_, dz_value_map_);
  iEvent.getByToken(dxy_value_map_token_, dxy_value_map_);
  iEvent.getByToken(dzsig_value_map_token_, dzsig_value_map_);
  iEvent.getByToken(dxysig_value_map_token_, dxysig_value_map_);
  iEvent.getByToken(lostInnerHits_value_map_token_, lostInnerHits_value_map_);
  iEvent.getByToken(quality_value_map_token_, quality_value_map_);
  iEvent.getByToken(trkPt_value_map_token_, trkPt_value_map_);
  iEvent.getByToken(trkEta_value_map_token_, trkEta_value_map_);
  iEvent.getByToken(trkPhi_value_map_token_, trkPhi_value_map_);
  iEvent.getByToken(msd_value_map_token_, msd_value_map_);
  iEvent.getByToken(n2b1_value_map_token_, n2b1_value_map_);

  // Match jet-gen jet
  genJetMatching(*jets, *genjets_);

  // Loop over jet
  for (std::size_t jet_n = 0; jet_n < jets->size(); jet_n++) {
    const auto &jet = (*jets)[jet_n];
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);
    sample_isQCD = isQCD_;

    if (std::find(genmatch_unmatched.begin(), genmatch_unmatched.end(), jet_n) != genmatch_unmatched.end()) {
      j_nCHadrons = -99;
      j_nBHadrons = -99;
      j_partonFlavour = -99;
      j_hadronFlavour = -99;
    } else {
      reco::JetFlavourInfo flavJet = genmatch_resultmap[jet_n];
      j_nCHadrons = flavJet.getcHadrons().size();
      j_nBHadrons = flavJet.getbHadrons().size();
      j_partonFlavour = flavJet.getPartonFlavour();
      j_hadronFlavour = flavJet.getHadronFlavour();
    }

    event_no = iEvent.id().event();
    j_pt = jet.pt();
    j_eta = jet.eta();
    j_phi = jet.phi();
    j_mass = jet.mass();
    j_no = jets->size();
    j_npfcands = jet.numberOfDaughters();

    const float etasign = jet.eta() > 0 ? 1 : -1;
    math::XYZVector jet_dir = jet.momentum().Unit();
    TVector3 jet_direction(jet.momentum().Unit().x(), jet.momentum().Unit().y(), jet.momentum().Unit().z());
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());

    std::vector<reco::CandidatePtr> daughters;
    for (const auto &dau : jet.daughterPtrVector()) { 
       auto cand = pfcands_->ptrAt(dau.key());
       daughters.push_back(cand);
    }
    std::sort(daughters.begin(), daughters.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
    
    for (const auto &cand : daughters) {
       const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

       auto candP4 = cand->p4();

       pfcand_charge.push_back(reco_cand->charge());
       pfcand_isEl.push_back(std::abs(reco_cand->pdgId()) == 11);
       pfcand_isMu.push_back(std::abs(reco_cand->pdgId()) == 13);
       pfcand_isChargedHad.push_back(std::abs(reco_cand->pdgId()) == 211);
       pfcand_isGamma.push_back(std::abs(reco_cand->pdgId()) == 22);
       pfcand_isNeutralHad.push_back(std::abs(reco_cand->pdgId()) == 130);
       pfcand_phirel.push_back(reco::deltaPhi(candP4, jet));
       pfcand_etarel.push_back(etasign * (candP4.eta() - jet.eta()));
       pfcand_deltaR.push_back(reco::deltaR(candP4, jet));
       pfcand_abseta.push_back(std::abs(candP4.eta()));
       pfcand_ptrel_log.push_back(std::log(candP4.pt() / jet.pt()));
       pfcand_ptrel.push_back(candP4.pt() / jet.pt());
       pfcand_erel_log.push_back(std::log(candP4.energy() / jet.energy()));
       pfcand_erel.push_back(candP4.energy() / jet.energy());
       pfcand_pt_log.push_back(std::log(candP4.pt()));
       pfcand_mask.push_back(1);
       pfcand_pt_log_nopuppi.push_back(std::log(cand->pt()));
       pfcand_e_log_nopuppi.push_back(std::log(cand->energy()));

       if ((*normchi2_value_map_)[cand] > 900) {
          pfcand_normchi2.push_back(0);
          pfcand_lostInnerHits.push_back(0);
          pfcand_quality.push_back(0);
          pfcand_dz.push_back(0);
          pfcand_dzsig.push_back(0);
          pfcand_dxy.push_back(0);
          pfcand_dxysig.push_back(0);
          pfcand_btagEtaRel.push_back(0);
          pfcand_btagPtRatio.push_back(0);
          pfcand_btagPParRatio.push_back(0);
      } else {
          pfcand_normchi2.push_back((*normchi2_value_map_)[cand]);
          pfcand_lostInnerHits.push_back((*lostInnerHits_value_map_)[cand]);
          pfcand_quality.push_back((*quality_value_map_)[cand]);
          pfcand_dz.push_back((*dz_value_map_)[cand]);
          pfcand_dzsig.push_back((*dzsig_value_map_)[cand]);
          pfcand_dxy.push_back((*dxy_value_map_)[cand]);
          pfcand_dxysig.push_back((*dxysig_value_map_)[cand]);
          float trk_px = (*trkPt_value_map_)[cand] * std::cos((*trkPhi_value_map_)[cand]);
          float trk_py = (*trkPt_value_map_)[cand] * std::sin((*trkPhi_value_map_)[cand]);
          float trk_pz = (*trkPt_value_map_)[cand] * std::sinh((*trkEta_value_map_)[cand]);
          math::XYZVector track_mom(trk_px, trk_py, trk_pz);
          TVector3 track_direction(trk_px, trk_py, trk_pz);
          double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);
          pfcand_btagEtaRel.push_back(reco::btau::etaRel(jet_dir, track_mom));
          pfcand_btagPtRatio.push_back(track_direction.Perp(jet_direction) / track_mag);
          pfcand_btagPParRatio.push_back(jet_dir.Dot(track_mom) / track_mag);
      }
    }
    tree->Fill();
    clearVars();
  }
}

void AK4JetFromNanoAODNtupleProducer::clearVars(){
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

void AK4JetFromNanoAODNtupleProducer::genJetMatching(edm::View<reco::Jet> jets, const reco::JetFlavourInfoMatchingCollection& genjets) {
  std::vector<std::tuple<int, int, float> > pairlist;

  for (unsigned int i=0; i<jets.size(); i++) {
    bool found_match = false;
    for (unsigned int j=0; j<genjets.size(); j++) {
      float dR = reco::deltaR(jets[i].eta(), jets[i].phi(), genjets[j].first.get()->eta(), genjets[j].first.get()->phi());
      if(dR < 0.4) {
        pairlist.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
  if(!found_match) {
    genmatch_unmatched.push_back(i);
    }
  }

  std::sort(pairlist.begin(), pairlist.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(pairlist.size() > 0) {
    reco::JetFlavourInfo genjet_assn = genjets[std::get<1>(pairlist[0])].second;
    genmatch_resultmap[std::get<0>(pairlist[0])] = genjet_assn;
    for(unsigned int k=1; k<pairlist.size(); k++) {
      if(std::get<0>(pairlist[k]) == std::get<0>(pairlist[0]) ||
         std::get<1>(pairlist[k]) == std::get<1>(pairlist[0])) {
        pairlist.erase(pairlist.begin() + k);
      }
    }
    pairlist.erase(pairlist.begin());
  }
}

void AK4JetFromNanoAODNtupleProducer::beginJob() {
}

void AK4JetFromNanoAODNtupleProducer::endJob() {
}

// define this as a plug-in
DEFINE_FWK_MODULE(AK4JetFromNanoAODNtupleProducer);
