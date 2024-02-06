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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "Run3ScoutingJetTagging/Analysis/interface/FatJetMatching.h"

using namespace deepntuples;

FatJetMatching<reco::Jet> ak8_reco_match;

class AK8JetFromNanoAODNtupleProducer :  public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit AK8JetFromNanoAODNtupleProducer(const edm::ParameterSet &);
  ~AK8JetFromNanoAODNtupleProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  std::map<int, reco::GenJet> genmatch_resultmap;
  std::vector<int> genmatch_unmatched;

private:
  typedef edm::View<reco::Candidate> CandidateView;

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void clearVars();
  virtual void genJetMatching(edm::View<reco::Jet>, const reco::GenJetCollection&);
 
  edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;
  edm::EDGetTokenT<reco::GenParticleCollection> gencand_token_;
  edm::EDGetTokenT<reco::GenJetCollection> genjet_token_;
  edm::Handle<CandidateView> pfcands_;
  edm::Handle<reco::GenParticleCollection> gencands_;
  edm::Handle<reco::GenJetCollection> genjets_;
  
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
  int mantissa_prescission_ = 10;

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

  float fj_pt;
  float fj_eta;
  float fj_phi;
  float fj_mass;
  float fj_msd;
  float fj_n2b1;
  int fj_no;
  int fj_npfcands;

  float fj_gen_mass;
  float fj_genjet_sdmass;

  int label_Top_bcq;
  int label_Top_bqq;
  int label_Top_bc;
  int label_Top_bq;
  int label_W_cq;
  int label_W_qq;
  int label_Z_bb;
  int label_Z_cc;
  int label_Z_qq;
  int label_H_bb;
  int label_H_cc;
  int label_H_qqqq;
  int label_H_tautau;
  int label_H_qq;
  int label_QCD_all;
  int sample_isQCD;

  int event_no;
};

AK8JetFromNanoAODNtupleProducer::AK8JetFromNanoAODNtupleProducer(const edm::ParameterSet &iConfig)
    : jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      gencand_token_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("gen_candidates"))),
      genjet_token_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("gen_jets"))) {

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
  fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
  fs->file().SetCompressionLevel(9);

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

  tree->Branch("fj_pt", &fj_pt);
  tree->Branch("fj_eta", &fj_eta);
  tree->Branch("fj_phi", &fj_phi);
  tree->Branch("fj_mass", &fj_mass);
  tree->Branch("fj_msd", &fj_msd);
  tree->Branch("fj_n2b1", &fj_n2b1);
  tree->Branch("fj_no", &fj_no);
  tree->Branch("fj_npfcands", &fj_npfcands);

  tree->Branch("fj_gen_mass", &fj_gen_mass);
  tree->Branch("fj_genjet_sdmass", &fj_genjet_sdmass);

  tree->Branch("label_Top_bcq", &label_Top_bcq);
  tree->Branch("label_Top_bqq", &label_Top_bqq);
  tree->Branch("label_Top_bc", &label_Top_bc);
  tree->Branch("label_Top_bq", &label_Top_bq);
  tree->Branch("label_W_cq", &label_W_cq);
  tree->Branch("label_W_qq", &label_W_qq);
  tree->Branch("label_Z_bb", &label_Z_bb);
  tree->Branch("label_Z_cc", &label_Z_cc);
  tree->Branch("label_Z_qq", &label_Z_qq);
  tree->Branch("label_H_bb", &label_H_bb);
  tree->Branch("label_H_cc", &label_H_cc);
  tree->Branch("label_H_qqqq", &label_H_qqqq);
  tree->Branch("label_H_tautau", &label_H_tautau);
  tree->Branch("label_H_qq", &label_H_qq);
  tree->Branch("label_QCD_all", &label_QCD_all);
  tree->Branch("sample_isQCD", &sample_isQCD);

  tree->Branch("event_no", &event_no);
}

AK8JetFromNanoAODNtupleProducer::~AK8JetFromNanoAODNtupleProducer() {}

void AK8JetFromNanoAODNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void AK8JetFromNanoAODNtupleProducer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // Input jets
  auto jets = iEvent.getHandle(jet_token_);
  // PF candidates
  iEvent.getByToken(pfcand_token_, pfcands_);
  // GEN candidates
  iEvent.getByToken(gencand_token_, gencands_);
  // GEN jets
  iEvent.getByToken(genjet_token_, genjets_);
  
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

    // Flavour matching
    auto ak8_label = ak8_reco_match.flavorLabel(jet, *gencands_, 0.8);

    label_Top_bcq = (ak8_label.first == FatJetMatching<reco::Jet>::Top_bcq);
    label_Top_bqq = (ak8_label.first == FatJetMatching<reco::Jet>::Top_bqq);
    label_Top_bc = (ak8_label.first == FatJetMatching<reco::Jet>::Top_bc);
    label_Top_bq = (ak8_label.first == FatJetMatching<reco::Jet>::Top_bq);
    label_W_cq = (ak8_label.first == FatJetMatching<reco::Jet>::W_cq);
    label_W_qq = (ak8_label.first == FatJetMatching<reco::Jet>::W_qq);
    label_Z_bb = (ak8_label.first == FatJetMatching<reco::Jet>::Z_bb);
    label_Z_cc = (ak8_label.first == FatJetMatching<reco::Jet>::Z_cc);
    label_Z_qq = (ak8_label.first == FatJetMatching<reco::Jet>::Z_qq);
    label_H_bb = (ak8_label.first == FatJetMatching<reco::Jet>::H_bb);
    label_H_cc = (ak8_label.first == FatJetMatching<reco::Jet>::H_cc);
    label_H_qqqq = (ak8_label.first == FatJetMatching<reco::Jet>::H_qqqq);
    label_H_tautau = (ak8_label.first == FatJetMatching<reco::Jet>::H_tautau);
    label_H_qq = (ak8_label.first == FatJetMatching<reco::Jet>::H_qq);
    label_QCD_all = (ak8_label.first == FatJetMatching<reco::Jet>::QCD_all);

    sample_isQCD = isQCD_;

    // Mass matching
    fj_gen_mass = (ak8_label.first < FatJetMatching<reco::Jet>::QCD_all && ak8_label.second) ? ak8_label.second->mass() : 0;

    if (std::find(genmatch_unmatched.begin(), genmatch_unmatched.end(), jet_n) != genmatch_unmatched.end()) {
      fj_genjet_sdmass = -99;
    } else {
      fj_genjet_sdmass = genmatch_resultmap[jet_n].mass();
    }

    event_no = iEvent.id().event();
    fj_pt = MiniFloatConverter::reduceMantissaToNbitsRounding(jet.pt(), mantissa_prescission_);
    fj_eta = MiniFloatConverter::reduceMantissaToNbitsRounding(jet.eta(), mantissa_prescission_);
    fj_phi = MiniFloatConverter::reduceMantissaToNbitsRounding(jet.phi(), mantissa_prescission_);
    fj_mass = MiniFloatConverter::reduceMantissaToNbitsRounding(jet.mass(), mantissa_prescission_);
    fj_msd = MiniFloatConverter::reduceMantissaToNbitsRounding((*msd_value_map_)[jets->ptrAt(jet_n)], mantissa_prescission_);
    fj_n2b1 = MiniFloatConverter::reduceMantissaToNbitsRounding((*n2b1_value_map_)[jets->ptrAt(jet_n)], mantissa_prescission_);
    fj_no = jets->size();
    fj_npfcands = jet.numberOfDaughters();

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
       pfcand_phirel.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(reco::deltaPhi(candP4, jet), mantissa_prescission_));
       pfcand_etarel.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(etasign * (candP4.eta() - jet.eta()), mantissa_prescission_));
       pfcand_deltaR.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(reco::deltaR(candP4, jet), mantissa_prescission_));
       pfcand_abseta.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(std::abs(candP4.eta()), mantissa_prescission_));
       pfcand_ptrel_log.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(std::log(candP4.pt() / jet.pt()), mantissa_prescission_));
       pfcand_ptrel.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(candP4.pt() / jet.pt(), mantissa_prescission_));
       pfcand_erel_log.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(std::log(candP4.energy() / jet.energy()), mantissa_prescission_));
       pfcand_erel.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(candP4.energy() / jet.energy(), mantissa_prescission_));
       pfcand_pt_log.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(std::log(candP4.pt()), mantissa_prescission_));
       pfcand_pt_log_nopuppi.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(std::log(cand->pt()), mantissa_prescission_));
       pfcand_e_log_nopuppi.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(std::log(cand->energy()), mantissa_prescission_));

       if ((*normchi2_value_map_)[cand] > 900) {
          pfcand_normchi2.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_lostInnerHits.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_quality.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_dz.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_dzsig.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_dxy.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_dxysig.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_btagEtaRel.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_btagPtRatio.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
          pfcand_btagPParRatio.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(0, mantissa_prescission_));
      } else {
          pfcand_normchi2.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*normchi2_value_map_)[cand], mantissa_prescission_));
          pfcand_lostInnerHits.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*lostInnerHits_value_map_)[cand], mantissa_prescission_));
          pfcand_quality.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*quality_value_map_)[cand], mantissa_prescission_));
          pfcand_dz.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*dz_value_map_)[cand], mantissa_prescission_));
          pfcand_dzsig.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*dzsig_value_map_)[cand], mantissa_prescission_));
          pfcand_dxy.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*dxy_value_map_)[cand], mantissa_prescission_));
          pfcand_dxysig.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding((*dxysig_value_map_)[cand], mantissa_prescission_));
          float trk_px = (*trkPt_value_map_)[cand] * std::cos((*trkPhi_value_map_)[cand]);
          float trk_py = (*trkPt_value_map_)[cand] * std::sin((*trkPhi_value_map_)[cand]);
          float trk_pz = (*trkPt_value_map_)[cand] * std::sinh((*trkEta_value_map_)[cand]);
          math::XYZVector track_mom(trk_px, trk_py, trk_pz);
          TVector3 track_direction(trk_px, trk_py, trk_pz);
          double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);
          pfcand_btagEtaRel.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(reco::btau::etaRel(jet_dir, track_mom), mantissa_prescission_));
          pfcand_btagPtRatio.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(track_direction.Perp(jet_direction) / track_mag, mantissa_prescission_));
          pfcand_btagPParRatio.push_back(MiniFloatConverter::reduceMantissaToNbitsRounding(jet_dir.Dot(track_mom) / track_mag, mantissa_prescission_));
      }
    }
    tree->Fill();
    clearVars();
  }
}

void AK8JetFromNanoAODNtupleProducer::clearVars(){
  pfcand_pt_log_nopuppi.clear();
  pfcand_pt_log.clear();
  pfcand_ptrel.clear();
  pfcand_ptrel_log.clear();
  pfcand_e_log_nopuppi.clear();
  pfcand_etarel.clear();
  pfcand_erel.clear();
  pfcand_erel_log.clear();
  pfcand_phirel.clear();
  pfcand_abseta.clear();
  pfcand_deltaR.clear();
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

void AK8JetFromNanoAODNtupleProducer::genJetMatching(edm::View<reco::Jet> jets, const reco::GenJetCollection& genjets) {
  std::vector<std::tuple<int, int, float> > pairlist;

  for (unsigned int i=0; i<jets.size(); i++) {
    bool found_match = false;
    for (unsigned int j=0; j<genjets.size(); j++) {
      float dR = reco::deltaR(jets[i].eta(), jets[i].phi(), genjets[j].eta(), genjets[j].phi());
      if(dR < 0.8) {
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
    reco::GenJet genjet_assn = genjets[std::get<1>(pairlist[0])];
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

void AK8JetFromNanoAODNtupleProducer::beginJob() {
}

void AK8JetFromNanoAODNtupleProducer::endJob() {
}

// define this as a plug-in
DEFINE_FWK_MODULE(AK8JetFromNanoAODNtupleProducer);
