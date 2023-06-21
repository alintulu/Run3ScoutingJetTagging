/*
 * FatJetMatching.hh
 *
 *  Created on: Feb 1, 2017
 *      Author: hqu
 */

#ifndef FATJETHELPERS_INTERFACE_FATJETMATCHING_H_
#define FATJETHELPERS_INTERFACE_FATJETMATCHING_H_

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <unordered_set>
#include <utility>

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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Math/interface/deltaR.h"

namespace deepntuples {


namespace ParticleID{
enum PdgId { p_unknown, p_d, p_u, p_s, p_c, p_b, p_t, p_bprime, p_tprime,
  p_eminus = 11, p_nu_e, p_muminus, p_nu_mu, p_tauminus, p_nu_tau,
  p_tauprimeminus, p_nu_tauprime, p_g = 21, p_gamma, p_Z0,
  p_Wplus, p_h0, p_Zprime0 = 32, p_Zpprime0, p_Wprimeplus, p_H0,
  p_A0, p_Hplus, p_G = 39, p_R0 = 41, p_H30 = 45, p_A20 = 46,
  p_LQ, p_cluster = 91, p_string,
  p_pi0 = 111, p_rho0 = 113, p_klong = 130, p_piplus = 211, p_rhoplus = 213, p_eta = 221, p_omega = 223,
  p_kshort = 310, p_k0, p_kstar0 = 313, p_kplus = 321, p_kstarplus = 323, p_phi = 333,
  p_dplus = 411, p_d0 = 421, p_dsplus = 431, p_b0 =511, p_bplus = 521,
  p_bs0 = 531, p_bcplus = 541,
  p_neutron = 2112, p_proton = 2212,
  p_sigmaminus = 3112, p_lambda0 = 3122,
  p_sigma0 = 3212, p_sigmaplus = 3222, p_ximinus = 3312, p_xi0 = 3322, p_omegaminus = 3334,
  p_sigmac0 = 4112, p_lambdacplus = 4122, p_xic0 = 4132,
  p_sigmacplus = 4212, p_sigmacpp = 4222, p_xicplus = 4232, p_omegac0 = 4332,
  p_sigmabminus = 5112, p_lambdab0 = 5122, p_xibminus = 5132, p_sigmab0 = 5212, p_sigmabplus = 5222,
  p_xib0 = 5232, p_omegabminus = 5332,
};
}

template <typename T>
class FatJetMatching {
public:
  enum FatJetFlavor {
    Default = 0,
    Top = 1,
    W = 2,
    Z = 3,
    H = 4,
  };

  enum FatJetLabel {
    Invalid=0,
    Top_all=10, Top_bcq, Top_bqq, Top_bc, Top_bq, Top_bele, Top_bmu, Top_btau,
    W_all=20, W_cq, W_qq,
    Z_all=30, Z_bb, Z_cc, Z_qq,
    H_all=40, H_bb, H_cc, H_qq, H_qqqq, H_tautau,
    QCD_all=50, QCD_bb, QCD_cc, QCD_b, QCD_c, QCD_others
  };

public:
  FatJetMatching() {}
  virtual ~FatJetMatching() {}

  std::pair<int, const reco::GenParticle*> flavorLabel(const T jet, const reco::GenParticleCollection& genParticles, double distR);

private:
  std::pair<int, const reco::GenParticle*> top_label(const T jet, const reco::GenParticle *parton, double distR);
  std::pair<int, const reco::GenParticle*> w_label(const T jet, const reco::GenParticle *parton, double distR);
  std::pair<int, const reco::GenParticle*> z_label(const T jet, const reco::GenParticle *parton, double distR);
  std::pair<int, const reco::GenParticle*> higgs_label(const T jet, const reco::GenParticle *parton, double distR);
  std::pair<int, const reco::GenParticle*> qcd_label(const T jet, const reco::GenParticleCollection& genParticles, double distR);


private:
  void printGenInfoHeader() const;
  void printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) const;
  const reco::GenParticle* getFinal(const reco::GenParticle* particle);
  bool isHadronic(const reco::GenParticle* particle) const;
  std::vector<const reco::GenParticle*> getDaughterQuarks(const reco::GenParticle* particle);

private:
  double jetR_ = 0.8;
  bool   requiresQuarksContained_ = true;

  bool debug_ = false;
  std::unordered_set<const reco::GenParticle*> processed_;


};

}

#endif /* FATJETHELPERS_INTERFACE_FATJETMATCHING_H_ */
