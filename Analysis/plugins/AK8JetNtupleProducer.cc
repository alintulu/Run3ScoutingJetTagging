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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
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

#include "Run3ScoutingAnalysisTools/Analysis/interface/FatJetMatching.h"

using namespace deepntuples;

FatJetMatching<fastjet::PseudoJet> ak8_pdseudojet_match;

class AK8JetNtupleProducer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit AK8JetNtupleProducer(const edm::ParameterSet&);
  ~AK8JetNtupleProducer();
		
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();
  bool isNeutralPdg(int);
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfcand_token_;
  const edm::EDGetTokenT<reco::GenParticleCollection> gencand_token_;
  edm::Handle<std::vector<Run3ScoutingParticle>> pfcands_;
  edm::Handle<reco::GenParticleCollection> gencands_;
  edm::Handle<reco::GenJetCollection> genjets_;

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

  bool isQCD_ = false;
  const edm::EDGetTokenT<reco::GenJetCollection> genjet_token_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particletable_token_;

  float dR_= 0.8;
};

AK8JetNtupleProducer::AK8JetNtupleProducer(const edm::ParameterSet& iConfig):
  pfcand_token_(consumes<std::vector<Run3ScoutingParticle> > (iConfig.getParameter<edm::InputTag>("pf_candidates"))),
  gencand_token_(consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("gen_candidates"))),
  genjet_token_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("gen_jets"))),
  particletable_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>())
{

  isQCD_ = iConfig.getUntrackedParameter<bool>("isQCD", false);

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
  fs->file().SetCompressionLevel(9);

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

AK8JetNtupleProducer::~AK8JetNtupleProducer() {
}

bool AK8JetNtupleProducer::isNeutralPdg(int pdgId) {
   const int neutralPdgs_array[] = {9, 21, 22, 23, 25, 12, 14, 16, 111, 130, 310, 311, 421, 511, 2112}; // gluon, gluon, gamma, Z0, higgs, electron neutrino, muon neutrino, tau neutrino, pi0, K0_L, K0_S; K0, neutron
   const std::vector<int> neutralPdgs(neutralPdgs_array, neutralPdgs_array + sizeof(neutralPdgs_array) / sizeof(int));
   if (std::find(neutralPdgs.begin(), neutralPdgs.end(), std::abs(pdgId)) == neutralPdgs.end())
     return false;
 
   return true;
}

void AK8JetNtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // PF candidates
  iEvent.getByToken(pfcand_token_, pfcands_);
  // GEN candidates
  iEvent.getByToken(gencand_token_, gencands_);
  // GEN jets
  iEvent.getByToken(genjet_token_, genjets_);
  // Particle table
  auto pdt = iSetup.getHandle(particletable_token_);
  const HepPDT::ParticleDataTable* particletable_ = pdt.product();

  // Create AK8 Jet
  std::vector<fastjet::PseudoJet> fj_part;
  fj_part.reserve(pfcands_->size());
  int pfcand_i = 0;
  for (auto pfcands_iter = pfcands_->begin(); pfcands_iter != pfcands_->end(); ++pfcands_iter) {

    auto pfm = particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId())) != nullptr
                        ? particletable_->particle(HepPDT::ParticleID(pfcands_iter->pdgId()))->mass()
                        : -99.f;
    if (pfm < -90) continue;

    math::PtEtaPhiMLorentzVector p4(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfm);
    fj_part.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
    fj_part.back().set_user_index(pfcand_i);
    pfcand_i++;
  }

  fastjet::JetDefinition ak8_def = fastjet::JetDefinition(fastjet::antikt_algorithm, dR_);
  fastjet::GhostedAreaSpec area_spec(5.0,1,0.01);
  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  // Substructure
  double sd_z_cut = 0.10;
  double sd_beta = 0;
  fastjet::contrib::SoftDrop sd_groomer = fastjet::contrib::SoftDrop(sd_beta, sd_z_cut, dR_);
  fastjet::contrib::EnergyCorrelatorN2 N2 = fastjet::contrib::EnergyCorrelatorN2(1.0);

  fastjet::ClusterSequenceArea ak8_cs(fj_part, ak8_def, area_def);
  std::vector<fastjet::PseudoJet> ak8_jets = fastjet::sorted_by_pt(ak8_cs.inclusive_jets(170.0));

  // Match jet-gen jet
  std::map<int, reco::GenJet> genmatch_resultmap;
  std::vector<int> genmatch_unmatched;
  std::vector<std::tuple<int, int, float> > parlist;

  int ak8_jet_idx = 0;

  for(unsigned int i=0; i<ak8_jets.size(); i++) {
    bool found_match = false;
    for(unsigned int j=0; j<genjets_->size(); j++) {
      float dR = reco::deltaR(ak8_jets[i].eta(), ak8_jets[i].phi(), (*genjets_)[j].eta(), (*genjets_)[j].phi());
      if(dR < dR_) {
        parlist.push_back(std::make_tuple(i, j, dR));
        found_match = true;
      }
    }
    ak8_jet_idx++;
    if(!found_match) {
       genmatch_unmatched.push_back(i);
    }
  }

  std::sort(parlist.begin(), parlist.end(), [](std::tuple<int, int, float> t1, std::tuple<int, int, float> t2){ return std::get<2>(t1) < std::get<2>(t2); });

  while(parlist.size() > 0) {
    reco::GenJet genjet_assn = (*genjets_)[std::get<1>(parlist[0])];
    genmatch_resultmap[std::get<0>(parlist[0])] = genjet_assn;
    for(unsigned int k=1; k<parlist.size(); k++) {
      if(std::get<0>(parlist[k]) == std::get<0>(parlist[0]) ||
         std::get<1>(parlist[k]) == std::get<1>(parlist[0])) {
        parlist.erase(parlist.begin() + k);
      }
    }
    parlist.erase(parlist.begin());
  }
  
  ak8_jet_idx = 0;
  for(auto &j: ak8_jets) {

    // Flavour matching
    auto ak8_label = ak8_pdseudojet_match.flavorLabel(j, *gencands_, dR_);

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

    fj_pt = j.pt();
    fj_eta = j.eta();
    fj_phi = j.phi();
    fj_mass = j.m();

    fastjet::PseudoJet sd_ak8 = sd_groomer(j);
    fj_msd = sd_ak8.m();
    fj_n2b1 = N2(sd_ak8);

    fj_no = ak8_jets.size();
    fj_npfcands = constituents.size();

    label_Top_bcq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bcq);
    label_Top_bqq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bqq);
    label_Top_bc = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bc);
    label_Top_bq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Top_bq);
    label_W_cq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::W_cq);
    label_W_qq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::W_qq);
    label_Z_bb = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Z_bb);
    label_Z_cc = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Z_cc);
    label_Z_qq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::Z_qq);
    label_H_bb = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_bb);
    label_H_cc = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_cc);
    label_H_qqqq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_qqqq);
    label_H_tautau = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_tautau);
    label_H_qq = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::H_qq);
    label_QCD_all = (ak8_label.first == FatJetMatching<fastjet::PseudoJet>::QCD_all);
    sample_isQCD = isQCD_;

    fj_gen_mass = (ak8_label.first < FatJetMatching<fastjet::PseudoJet>::QCD_all && ak8_label.second) ? ak8_label.second->mass() : 0;

    if (std::find(genmatch_unmatched.begin(), genmatch_unmatched.end(), ak8_jet_idx) != genmatch_unmatched.end()) {
      fj_genjet_sdmass = -99;
    } else {
      fj_genjet_sdmass = genmatch_resultmap[ak8_jet_idx].mass();
    }

    event_no = iEvent.id().event();

    ak8_jet_idx++;

    tree->Fill();	
    clearVars();
  }
}

void AK8JetNtupleProducer::clearVars(){
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

void AK8JetNtupleProducer::beginJob() {
  
}

void AK8JetNtupleProducer::endJob() {
}

void AK8JetNtupleProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void AK8JetNtupleProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void AK8JetNtupleProducer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void AK8JetNtupleProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void AK8JetNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(AK8JetNtupleProducer);
