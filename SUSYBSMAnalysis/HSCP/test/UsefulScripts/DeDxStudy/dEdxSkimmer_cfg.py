import FWCore.ParameterSet.Config as cms

process = cms.Process("DEDXUNCSKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')


process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')


process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
      SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/CADE98C7-2D0F-E511-85EC-02163E0145C8.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/CC9920A9-200F-E511-A61E-02163E013614.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/D8FE05F1-220F-E511-8831-02163E01220A.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/E670676D-380F-E511-A3B4-02163E0136B6.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/EA257BC4-230F-E511-B463-02163E0145D2.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/EE5A5457-5B0F-E511-82C3-02163E0143CD.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/398/00000/F26AFA89-2D0F-E511-99C8-02163E013502.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/408/00000/4AE2D067-540F-E511-A489-02163E01180A.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/408/00000/EC359307-520F-E511-A5FE-02163E01460C.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/487/00000/D0156FCE-AE0F-E511-B9B2-02163E011B79.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/489/00000/E0678161-B10F-E511-B0FF-02163E0135D4.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/02C92A21-C410-E511-83B8-02163E014349.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/10AEC5C0-F60F-E511-999B-02163E01424D.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/3227E81F-EE0F-E511-BD4A-02163E014216.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/3E5D79DE-F20F-E511-AE3D-02163E013491.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/4287EE6D-F20F-E511-BA78-02163E01440B.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/626D3017-F20F-E511-A10D-02163E011B95.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/923EAEDD-F20F-E511-814F-02163E0136C3.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/A02C30F3-F50F-E511-987B-02163E0140E7.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/AECB4CB6-FA0F-E511-B27D-02163E0142DF.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/C89AEA36-F30F-E511-A0E4-02163E011A9B.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/EEA2C288-F00F-E511-97B4-02163E0140FA.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/491/00000/FE743D85-F10F-E511-A932-02163E0138C6.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/545/00000/A8CF8124-6310-E511-A98C-02163E0118BB.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/545/00000/E2EC53D4-4310-E511-80A6-02163E011B58.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/547/00000/525C5E42-4110-E511-A528-02163E011C9E.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/551/00000/3AE7296D-4610-E511-89C7-02163E0145C5.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/554/00000/5CCEF846-4E10-E511-9A78-02163E0140F7.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/571/00000/9232E9D0-5C10-E511-BE37-02163E011B58.root",
   "/store/data/Run2015A/SingleMu/RECO/PromptReco-v1/000/247/571/00000/ECE21045-6110-E511-991E-02163E011B58.root",
   ),
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)

#process.GlobalTag.globaltag = GTAG
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "GR_P_V56", '')


from RecoTracker.TrackProducer.TrackRefitter_cfi import *
process.RefitterForDeDx = TrackRefitter.clone(
      src = cms.InputTag("generalTracks"),
      NavigationSchool = cms.string("")                                   
)

from RecoTracker.DeDx.dedxEstimators_cff import *
process.dedxHitInfo = dedxHitInfo.clone()
process.dedxHitInfo.tracks=cms.InputTag("RefitterForDeDx")
process.dedxHitInfo.trajectoryTrackAssociation = cms.InputTag("RefitterForDeDx")
process.dedxHitInfo.minTrackPt = cms.double(0.0)

#make the pool output
process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring(
         "drop *",
         "keep EventAux_*_*_*",
         "keep LumiSummary_*_*_*",
         "keep *_RefitterForDeDx_*_DEDXUNCSKIM",
         "keep *_siPixelClusters_*_*",
         "keep *_siStripClusters_*_*",
         "keep *_dedxHitInfo_*_DEDXUNCSKIM",
    ),
    fileName = cms.untracked.string("dEdxSkim.root"),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('*')
    ),
)

#schedule the sequence
process.p = cms.Path(process.RefitterForDeDx * process.dedxHitInfo)
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.p, process.endPath1)
