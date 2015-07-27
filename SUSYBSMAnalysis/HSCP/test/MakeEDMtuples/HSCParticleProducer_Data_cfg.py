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
   'root://cms-xrd-global.cern.ch//store/data/Run2015B/ZeroBias/RECO/PromptReco-v1/000/251/252/00000/2E803A6D-8C27-E511-B186-02163E012073.root',
)

#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )
