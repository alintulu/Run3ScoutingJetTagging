_This repository creates jet-level Run 3 scouting ntuples for training mass regression with ParticleNet. It was originally created for the "Jet tagging using the Run3 Scouting data" hackathon on November 8-11 2021 and then updated._

# Download, compile and run code

The code was develop in `CMSSW_12_3_0`.

1. Log into lxplus
2. Prepare the CMSSW release

```
cmsrel CMSSW_12_3_0
cd CMSSW_12_3_0/src
cmsenv
```

3. Clone this repository and compile

```
git clone git@github.com:alintulu/Run3ScoutingJetTagging.git
scram b -j 8
```

# Produce ntuples

## AK4 producer

The AK4 producer creates ntuples which can be used for

- b tagging
- c tagging
- bb tagging
- cc tagging
- uds tagging
- g tagging

### Example

To create an example ntuple from a QCD datasets, make sure you have a valid grid proxy and

1. If you have [this PR merged](https://github.com/cms-sw/cmssw/pull/40438) in your CMSSW environment you can run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK4FromNanoAOD.py inputDataset="/QCD_"
```

2. Else run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK4.py inputDataset="/QCD_"
```

### Derivation of labels

For training, the labels (jet flavour) are assigned by the usage of 4 jet characteristics:
- `j_nCHadrons`: the number of c hadrons in the jet
- `j_nBHadrons`: the number of b hadrons in the jet
- `j_partonFlavour`: the parton flavour of the jet
- `j_hadronFlavour`: the hadron flavour of the jet

The assignment goes as follows:
- b tagged: `j_nBHadrons == 1`
- bb tagged: `j_nBHadrons > 1`
- c tagged: `j_nCHadrons == 1`
- cc tagged: `j_nCHadrons > 1`
- uds tagged: `(j_hadronFlavour==0) & (np.abs(j_partonFlavour)>0) & (np.abs(j_partonFlavour)<4)`
- g tagged: `(j_hadronFlavour==0) & (j_partonFlavour==21)`
- undefined: `(j_hadronFlavour==0) & (j_partonFlavour==0)`

The jet characteristics are assigned to each scouting jet as follows:
1. Find the jet flavour of GEN jets through the [ghost association method](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#Hadron_parton_based_jet_flavour)
2. Match scouting jets to GEN jets by requiring dR < 0.4
3. The four characteristics assosicated to the jet flavour are assigned to the scouting jet by the matched GEN jet

## AK8 producer

The AK8 producer creates ntuples which can be used for

- Higgs to bb vs QCD tagging
- Higgs to cc vs QCD tagging
- Higgs to udsuds vs QCD tagging
- Higgs to tautau vs QCD tagging
- Mass regression

### Examples

To create an example ntuple from a QCD datasets, make sure you have a valid grid proxy and

1. If you have [this PR merged](https://github.com/cms-sw/cmssw/pull/40438) in your CMSSW environment you can run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK8FromNanoAOD.py inputDataset="/QCD_"
```

2. Else run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK8.py inputDataset="/QCD_"
```

### Derivation of labels for flavour tagging

For training, the labels (jet flavour) are defined as:
- `label_H_bb`: the jet was created from a Higgs boson decaying to two b quarks
- `label_H_cc`: the jet was created from a Higgs boson decaying to two c quarks
- `label_H_qq`: the jet was created from a Higgs boson decaying to two light (uds) quarks
- `label_H_tautau`: the jet was created from a Higgs boson decaying to two tau leptons
- `label_QCD_all`: the code was unable to defined the origin of the jet, and hence labeled it as QCD

The jet labels are assigned to each scouting jet as follows:
1. For each jet, match a final state Higgs boson by requiring dR < 0.8 (if not, assign jet as QCD)
2. Find the final state daughters of the Higgs boson, require dR(jet, daughter) < 0.8 and tag the jet as `label_H_bb` if they are both b quarks etc
3. Tag the jet as QCD if final state daughters don't match to the jet

### Derivation of labels for mass regression

For training, the labels (jet mass) are defined as:
- `fj_gen_mass`: the generated mass of the GEN particle which the scouting jet originates from
- `fj_genjet_mass`: the mass of the GEN jet matched to the scouting jet

The former is used when the jet is a signal jet from a signal sample and the latter if the jet is a QCD jet from a QCD sample.

The jet labels are assigned to each scouting jet as follows:
1. For each signal jet, match a final state Higgs boson by requiring dR < 0.8 and assign `fj_gen_mass` the value of the generated Higgs mass
2. For signal jets with no matched Higgs boson and for QCD jets, `fj_gen_mass` is set to 0
3. For each jet, match a GEN jet by requiring dR < 0.8
4. Assign `fj_genjet_mass` to the value of the matched GEN jet
5. If there is no match, assign `fj_genjet_mass` to -99
