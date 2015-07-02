import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = False
isBckg = False
isData = True
isSkimmedSample = False
GTAG = 'GR_P_V56'
OUTPUTFILE = 'HSCP.root'

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
   '/store/data/Run2015A/SingleMu/AOD/PromptReco-v1/000/247/551/00000/0E25C85D-4610-E511-AF22-02163E0141D2.root',
   '/store/data/Run2015A/SingleMu/AOD/PromptReco-v1/000/247/554/00000/289A2FDD-4D10-E511-8ECA-02163E014113.root',
   '/store/data/Run2015A/SingleMu/AOD/PromptReco-v1/000/247/571/00000/7E30BF43-6110-E511-A437-02163E011C8D.root',
)

#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )
