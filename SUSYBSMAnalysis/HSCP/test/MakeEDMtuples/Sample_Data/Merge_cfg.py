import FWCore.ParameterSet.Config as cms
#process = cms.Process("MergeHLT")
process = cms.Process("MergeHSCP")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("FWCore.MessageService.MessageLogger_cfi")
from SUSYBSMAnalysis.HSCP.HSCPVersion_cff import *

process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
	   ['file:HSCP_1.root',
	    'file:HSCP_2.root',
	    'file:HSCP_3.root',
	    'file:HSCP_4.root',
	    'file:HSCP_5.root',
	    'file:HSCP_6.root',
	    'file:HSCP_7.root',
	    'file:HSCP_8.root',
	    'file:HSCP_9.root',
	    'file:HSCP_10.root',
	    'file:HSCP_11.root',
	    'file:HSCP_12.root',
	    'file:HSCP_13.root',
	    'file:HSCP_14.root',
	    'file:HSCP_15.root',
	    'file:HSCP_16.root',
	    'file:HSCP_17.root',
	    'file:HSCP_18.root',
	    'file:HSCP_19.root',
	    'file:HSCP_20.root',
	    'file:HSCP_21.root',
	    'file:HSCP_22.root',
	    'file:HSCP_23.root',
	    'file:HSCP_24.root',
	    'file:HSCP_25.root',
	    'file:HSCP_26.root',
	    'file:HSCP_27.root',
	    'file:HSCP_28.root',
	    'file:HSCP_29.root',
	    'file:HSCP_30.root',
	    'file:HSCP_31.root',
	    'file:HSCP_32.root']
   )
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.HSCPHLTDuplicate = cms.EDFilter("HSCPHLTFilter",
   RemoveDuplicates = cms.bool(True),
   TriggerProcess   = cms.string("HLT"),
   MuonTrigger1Mask    = cms.int32(0),  #Activated
   PFMetTriggerMask    = cms.int32(0),  #Activated
   L2MuMETTriggerMask  = cms.int32(0),  #Activated
)
process.DuplicateFilter = cms.Path(process.HSCPHLTDuplicate   )

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.HSCPHLTTriggerMuDeDx = process.hltHighLevel.clone()
process.HSCPHLTTriggerMuDeDx.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
process.HSCPHLTTriggerMuDeDx.andOr = cms.bool( True ) #OR
process.HSCPHLTTriggerMuDeDx.throw = cms.bool( False )
process.HSCPHLTTriggerMuDeDx.HLTPaths = ["HLT_Mu*_dEdx*"]
process.HSCPHLTTriggerMuDeDxFilter = cms.Path(process.HSCPHLTTriggerMuDeDx   )

process.HSCPHLTTriggerMetDeDx = process.HSCPHLTTriggerMuDeDx.clone() 
process.HSCPHLTTriggerMetDeDx.HLTPaths = ["HLT_MET*_dEdx*"]
process.HSCPHLTTriggerMetDeDxFilter = cms.Path(process.HSCPHLTTriggerMetDeDx   )

process.HSCPHLTTriggerHtDeDx = process.HSCPHLTTriggerMuDeDx.clone()
process.HSCPHLTTriggerHtDeDx.HLTPaths = ["HLT_HT*_dEdx*"]
process.HSCPHLTTriggerHtDeDxFilter = cms.Path(process.HSCPHLTTriggerHtDeDx   )

#process.HSCPHLTTriggerMu
if CMSSW4_2 or CMSSW4_4:
  #This needs to be done differently for 2011 data because HLT_Mu40 did not exist in trigger menu at beginning of run
  process.HSCPHLTTriggerMu = cms.EDFilter("HSCPHLTFilter",
     RemoveDuplicates = cms.bool(False),
     TriggerProcess  = cms.string("HLT"),
     MuonTrigger1Mask    = cms.int32(1),  #Activated
     PFMetTriggerMask    = cms.int32(0),  #Activated
     L2MuMETTriggerMask  = cms.int32(0),  #Activated
   )
else:
   process.HSCPHLTTriggerMu = process.HSCPHLTTriggerMuDeDx.clone()
   process.HSCPHLTTriggerMu.HLTPaths = ["HLT_Mu40_*"]

process.HSCPHLTTriggerMuFilter = cms.Path(process.HSCPHLTTriggerMu   )

process.HSCPHLTTriggerMet = process.HSCPHLTTriggerMuDeDx.clone()
process.HSCPHLTTriggerMet.HLTPaths = ["HLT_MET80_*"]
process.HSCPHLTTriggerMetFilter = cms.Path(process.HSCPHLTTriggerMet   )

process.HSCPHLTTriggerPFMet = process.HSCPHLTTriggerMuDeDx.clone()
#Name change from 2011 to 2012
if CMSSW4_2 or CMSSW4_4:
    process.HSCPHLTTriggerPFMet.HLTPaths = ["HLT_PFMHT150_*"]
else:
    process.HSCPHLTTriggerPFMet.HLTPaths = ["HLT_PFMET150_*"]
process.HSCPHLTTriggerPFMetFilter = cms.Path(process.HSCPHLTTriggerPFMet   )

process.HSCPHLTTriggerHt = process.HSCPHLTTriggerMuDeDx.clone()
process.HSCPHLTTriggerHt.HLTPaths = ["HLT_HT650_*"]
process.HSCPHLTTriggerHtFilter = cms.Path(process.HSCPHLTTriggerHt   )

if CMSSW4_2 or CMSSW4_4:
   #Needs to be done separately as had lower threshold prescaled trigger in menu in 2011
   process.HSCPHLTTriggerL2Mu = cms.EDFilter("HSCPHLTFilter",
     RemoveDuplicates = cms.bool(False),
     TriggerProcess   = cms.string("HLT"),
     MuonTrigger1Mask    = cms.int32(0),  #Activated
     PFMetTriggerMask    = cms.int32(0),  #Activated
     L2MuMETTriggerMask  = cms.int32(1),  #Activated
   )
else:
   process.HSCPHLTTriggerL2Mu = process.HSCPHLTTriggerMuDeDx.clone()
   process.HSCPHLTTriggerL2Mu.HLTPaths = ["HLT_L2Mu*MET*"]

process.HSCPHLTTriggerL2MuFilter = cms.Path(process.HSCPHLTTriggerL2Mu   )


process.HSCPHLTTriggerCosmic = process.HSCPHLTTriggerMuDeDx.clone()
process.HSCPHLTTriggerCosmic.HLTPaths = ["HLT_L2Mu*NoBPTX*"]
process.HSCPHLTTriggerCosmicFilter = cms.Path(process.HSCPHLTTriggerCosmic   )

process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring("keep *"),
    #fileName = cms.untracked.string('/uscmst1b_scratch/lpc1/3DayLifetime/farrell/HSCPEDMUpdateData2012_30Nov2012/XXX_OUTPUT_XXX.root'),
    fileName = cms.untracked.string('output.root'),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('DuplicateFilter')
    ),
)

process.endPath = cms.EndPath(process.Out)

#process.schedule = cms.Schedule(process.DuplicateFilter, process.HSCPHLTTriggerMuDeDxFilter, process.HSCPHLTTriggerMetDeDxFilter, process.HSCPHLTTriggerHtDeDxFilter, process.HSCPHLTTriggerMuFilter, process.HSCPHLTTriggerMetFilter, process.HSCPHLTTriggerPFMetFilter, process.HSCPHLTTriggerHtFilter, process.HSCPHLTTriggerL2MuFilter, process.HSCPHLTTriggerCosmicFilter, process.endPath)
process.schedule = cms.Schedule(process.DuplicateFilter, process.endPath)
