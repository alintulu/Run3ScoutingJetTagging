import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *

from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register('inputDataset',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input dataset"
)

process = cms.Process("LL")
params.parseArguments()

process.load("PhysicsTools.NanoAOD.run3scouting_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring("/store/mc/Run3Summer22EEMiniAODv3/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/40000/04dd812c-f268-4ce5-82f1-079fbe668795.root"),
        secondaryFileNames = cms.untracked.vstring(
         	"/store/mc/Run3Summer22EEDRPremix/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/40000/6b3243ae-152c-4c60-bdf0-4894442c362e.root",
         	"/store/mc/Run3Summer22EEDRPremix/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/40000/b57139ff-9119-493f-9c74-8120938ff94e.root",
         	"/store/mc/Run3Summer22EEDRPremix/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/40000/d328cfe6-2695-4d94-8f6c-f5ca0ae67073.root",
         	"/store/mc/Run3Summer22EEDRPremix/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/40000/d4eca600-c846-494b-b922-207eb0512f20.root"
        )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("test.root")
)


# Create SoftDrop pruned GEN jets
process.load('PhysicsTools.NanoAOD.jetMC_cff')
process.genJetSequence = cms.Sequence(
   process.patJetPartonsNano+
   process.genJetFlavourAssociation
)

# Create ParticleNet ntuple
process.tree = cms.EDAnalyzer("AK4JetNtupleProducer",
      isQCD = cms.untracked.bool( '/QCD_' in params.inputDataset ),
      gen_jets = cms.InputTag( "genJetFlavourAssociation" ),
      pf_candidates = cms.InputTag( "hltScoutingPFPacker" ),
)

process.p = cms.Path(process.genJetSequence*process.tree)
