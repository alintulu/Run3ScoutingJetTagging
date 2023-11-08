from coffea.nanoevents import NanoAODSchema

class ScoutingNanoAODSchema(NanoAODSchema):
    """ScoutingNano schema builder
    """

    mixins = {
        "ScoutingParticle" : "PFCand",
        "ScoutingJet" : "Jet",
        "ScoutingFatJet" : "FatJet",
        "ScoutingPhoton" : "PtEtaPhiMCollection",
        "ScoutingElectron" : "PtEtaPhiMCollection",
        "ScoutingMuon" : "PtEtaPhiMCollection",
        "ScoutingTrack" : "PtEtaPhiMCollection",
        "ScoutingPrimaryVertex" : "PtEtaPhiMCollection",
        "ScoutingDisplacedVertex" : "PtEtaPhiMCollection",
        "ScoutingMET" : "MissingET",
        "ScoutingRho" : "Rho",
        "GenJet": "PtEtaPhiMCollection",
        "GenJetAK8": "PtEtaPhiMCollection",
        "GenDressedLepton": "PtEtaPhiMCollection",
        "GenIsolatedPhoton": "PtEtaPhiMCollection",
        "GenMET": "MissingET",
        "GenPart": "GenParticle",
        "SubGenJetAK8": "PtEtaPhiMCollection",
        "GenVisTau": "GenVisTau",
    }
    all_cross_references = {
        "GenPart_genPartIdxMother": "GenPart",
        "GenVisTau_genPartIdxMother": "GenPart",
        "ScoutingJet_genJetIdx" : "GenJet",
        "ScoutingFatJet_genJetAK8Idx" : "GenJetAK8",
    }

class JMENanoAODSchema(NanoAODSchema):
    """JMENano schema builder
    JMENano is an extended NanoAOD format that includes various jet collections down to low pt for JME studies
    More info at https://twiki.cern.ch/twiki/bin/viewauth/CMS/JMECustomNanoAOD
    Customization at https://github.com/nurfikri89/cmssw/blob/master/PhysicsTools/NanoAOD/python/custom_jme_cff.py
    """

    mixins = {
        **NanoAODSchema.mixins,
        "JetCHS": "Jet",
        "JetCalo": "Jet",
        #"JetPuppi": "Jet", # PUPPI is default "Jet" in JMENano, included here for completion
        "FatJetForJEC": "FatJet",
        "FatJetCHS": "FatJet",
        "TrigObjJMEAK4": "PtEtaPhiMCollection",
        "TrigObjJMEAK8": "PtEtaPhiMCollection"
    }

    all_cross_references = {
        **NanoAODSchema.all_cross_references,
        "JetCHS_genJetIdx": "GenJet",
        "JetCalo_genJetIdx": "GenJet",
        #"JetPuppi_genJetIdx": "GenJet", # PUPPI is default "Jet" in JMENano, included here for completion
        "FatJetForJEC_genJetIdx": "GenJetAK8ForJEC",
        "FatJetCHS_genJetIdx": "GenJetAK8ForJEC",
    }

class ScoutingJMENanoAODSchema(JMENanoAODSchema):

    mixins = {
        **JMENanoAODSchema.mixins,
        **ScoutingNanoAODSchema.mixins,
    }
    all_cross_references = {
        **JMENanoAODSchema.all_cross_references,
        **ScoutingNanoAODSchema.all_cross_references,
    }
