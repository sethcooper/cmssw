import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = False
isBckg = True
isData = False
isSkimmedSample = False
GTAG = 'MCRUN2_74_V7::All'

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
'/store/relval/CMSSW_7_4_0/RelValZMM_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_GENSIM_7_1_15-v1/00000/0C579384-8EDD-E411-B509-0025905A611C.root',
)

#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/BuildHSCParticles/HSCParticleProducer_cfg.py') )
