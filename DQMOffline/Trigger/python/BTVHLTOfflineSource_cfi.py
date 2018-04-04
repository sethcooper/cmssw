import FWCore.ParameterSet.Config as cms

BTVHLTOfflineSource = cms.EDAnalyzer(
    "BTVHLTOfflineSource",
    #
    dirname                 = cms.untracked.string("HLT/BTV"),
    processname             = cms.string("HLT"),
    verbose                 = cms.untracked.bool(False),
    #
    triggerSummaryLabel     = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerResultsLabel     = cms.InputTag("TriggerResults","","HLT"),
    onlineDiscrLabelPF      = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsPF", "probb"),
    onlineDiscrLabelCalo    = cms.InputTag("hltDeepCombinedSecondaryVertexBJetTagsCalo", "probb"),
    offlineDiscrLabelPF     = cms.InputTag("pfDeepCSVJetTags", "probb"),
    offlineDiscrLabelCalo   = cms.InputTag("pfDeepCSVJetTags", "probb"),
    hltFastPVLabel          = cms.InputTag("hltFastPrimaryVertex"),
    hltPFPVLabel            = cms.InputTag("hltVerticesPFSelector"),
    hltCaloPVLabel          = cms.InputTag("hltVerticesL3"),
    offlinePVLabel          = cms.InputTag("offlinePrimaryVertices"),
    turnon_threshold_low    = cms.float32(0.2),
    turnon_threshold_medium = cms.float32(0.5),
    turnon_threshold_high   = cms.float32(0.8),
    #
    pathPairs = cms.VPSet(
        cms.PSet(
            pathname = cms.string("HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_"),
            pathtype = cms.string("PF"),
        ),
        cms.PSet(
            pathname = cms.string("HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_"),
            pathtype = cms.string("Calo"),
        )
    )
)

btvHLTDQMSourceExtra = cms.Sequence(
)
