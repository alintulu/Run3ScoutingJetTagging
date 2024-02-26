# Contents of this repository

1. [Ntuples: How to produce training samples](#produce-ntuples)
   - [Getting started](#getting-started)
   - [AK4 producer](#ak4-producer)
        - [Which are the input samples?](#ak4-input-samples)
        - [How does one produce a ntuple?](#ak4-example)
        - [How are the labels derived?](#ak4-derivation-of-labels)
   - [AK8 producer](#ak8-producer)
        - [Which are the input samples?](#ak8-input-samples)
        - [How does one produce a ntuple?](#ak8-example)
        - [How are the labels for flavour tagging derived?](#ak8-derivation-of-labels-for-flavour-tagging)
        - [How are the labels for mass regression derived?](#ak8-derivation-of-labels-for-mass-regression)
2. [Training: How to run the training on GPU](#training)
   - [Problems with the training](#problems)
3. [Evaluation: How to evaluate the training](#evaluation)
4. [Calibration: How to calibrate the result](#calibration)

# Produce ntuples

## Getting started

The code was develop in `CMSSW_13_1_0_pre1`. While newer releases should also work, the have not been tested.

1. Log into lxplus
2. Prepare the CMSSW release

```
cmsrel CMSSW_13_1_0_pre1
cd CMSSW_13_1_0_pre1/src
cmsenv
```

3. Clone this repository and compile

```
git clone git@github.com:alintulu/Run3ScoutingJetTagging.git
scram b -j 8
```

## AK4 producer

The AK4 producer creates ntuples which can be used for

- b tagging
- c tagging
- bb tagging
- cc tagging
- uds tagging
- g tagging

### AK4 input samples 

**Please note**, before the Run3Winter24 MC campaign the Run 3 scouting collections only existed as AODSIM. However, the ntupiliser requires GEN information which is only available as MINIAODSIM. As a result, when using samples preceeding the Run3Winter24 MC campaigns  one needs to provide both the AODSIM and MINIAODSIM samples (see [here](https://github.com/alintulu/Run3ScoutingJetTagging/blob/main/Analysis/test/AK4.py#L24-L30) for an example). In contrast, for samples produces with the Run3Winter24 MC campaign or a later campaign only requires the MINIADOSIM sample.

AK4 flavour tagging is usually trained on three types of input samples:

1. BulkGravitonToHH

An example dataset is:

```
/BulkGravitonToHH_MX960_MH121_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v4/AODSIM
```

Here the generated mass of the BulkGraviton is 960 GeV, while the generated mass for the Higgs boson is 121 GeV. In order to avoid mass sculpting (which happens when the model learns to distinguish signal from background by looking at the mass of the particle), we use several different datasets with different generated masses for both particles.

2. QCD
3. TTbar

### AK4 example

To create an example ntuple from a QCD datasets, make sure you have a valid grid proxy and

1. If you have [this PR merged](https://github.com/cms-sw/cmssw/pull/40438) in your CMSSW environment you can run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK4FromNanoAOD.py inputDataset="/QCD_"
```

2. Else run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK4.py inputDataset="/QCD_"
```

### AK4 derivation of labels

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

### AK8 input samples

**Please note**, before the Run3Winter24 MC campaign the Run 3 scouting collections only existed as AODSIM. However, the ntupiliser requires GEN information which is only available as MINIAODSIM. As a result, when using samples preceeding the Run3Winter24 MC campaigns  one needs to provide both the AODSIM and MINIAODSIM samples (see [here](https://github.com/alintulu/Run3ScoutingJetTagging/blob/main/Analysis/test/AK4.py#L24-L30) for an example). In contrast, for samples produces with the Run3Winter24 MC campaign or a later campaign only requires the MINIADOSIM sample.

AK8 flavour tagging and mass regression is usually trained on two types of input samples:

1. BulkGravitonToHH

An example dataset is:

```
/BulkGravitonToHH_MX960_MH121_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v4/AODSIM
```

Here the generated mass of the BulkGraviton is 960 GeV, while the generated mass for the Higgs boson is 121 GeV. In order to avoid mass sculpting (which happens when the model learns to distinguish signal from background by looking at the mass of the particle), we use several different datasets with different generated masses for both particles.

2. QCD

### AK8 Example

To create an example ntuple from a QCD datasets, make sure you have a valid grid proxy and

1. If you have [this PR merged](https://github.com/cms-sw/cmssw/pull/40438) in your CMSSW environment you can run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK8FromNanoAOD.py inputDataset="/QCD_"
```

2. Else run
```
cmsRun Run3ScoutingJetTagging/Analysis/test/AK8.py inputDataset="/QCD_"
```

### AK8 derivation of labels for flavour tagging

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

### AK8 derivation of labels for mass regression

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

# Training

The training is performed using [Weaver](https://github.com/hqucms/weaver-core/). Follow [these instructions](https://github.com/hqucms/weaver-core/#set-up-your-environment) in order to set up the correct environment. The linked repository contains basic information about the configuration files needed to perform the training. More indepth information about ParticleNet can be found over at [CMS Machine Learning Documetation](https://cms-ml.github.io/documentation/inference/particlenet.html).

Two types of configuration files are needed by Weaver:

1. A data configuration YAML file, which describe how to process the input data
3. A network configuration python file, which defines the network architecture we are using for our model

The [Training](Training) folder contains the configuration files for

1. [AK4 flavour tagging](Training/AK4/flavour)
2. [AK8 flavour tagging](Training/AK8/flavour)
3. [AK8 mass regression](Training/AK8/massreg)

Each training type contains two folders called `data` and `network`. These hosts the data and network configuration files respectively. The network for AK4 and AK8 flavour tagging are the same, however the mass regression network differs in a few aspects. For flavour tagging, [the last unit of the network is the soft-max function](https://github.com/alintulu/Run3ScoutingJetTagging/blob/main/Training/AK8/flavour/network/flavour.py#L30C9-L30C37), which normalises the output and obtains a joint probability distribution for the output classes

```
'output_names': ['softmax'],
```

While for mass regression, [the soft-max unit is removed](https://github.com/alintulu/Run3ScoutingJetTagging/blob/main/Training/AK8/massreg/network/massreg.py#L31) which allows the output to be a single real number.

The training is performed by running the `train.sh` file. Set the desired batch size, number of epochs etc before starting the training.

## Problems

If you're training fails unexpectedly without an obvious cause, it may be because your input samples have too few events. This can be fixed by e.g. adjusting the batchsize, however, in our case simply merging the input samples helped.

# Evaluation

The output folder stores the trained PyTorch models after every epoch, as well as the log file that records the loss and accuracy in the runtime. The predict step also produces a predicted root file in the output folder, including the truth label, the predicted store, and several observer variables we provided in the data configuration file. With the predicted root file, you can evaluate the performance of your training. 

The [Evaluation](Evaluation) folder contains some notebooks that may be helpful when creating ROC curves for AK4 and AK8 flavour tagging as well plots to determine the performance of the training. **Please note**, the scripts are not complete but are added to this repository as they may provide a good starting point.

# Calibration

It is required to calibrate both the flavour tagging and mass regerssion algorithms. At the moment no complete example are stored in this repository. However, as starting point based on Section 4.2.2 of [AN-2021/004](https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2021/004) of how to calibrate the mass regression can be found in the [Calibration/AK8/massreg](Calibration/AK8/massreg) folder.
