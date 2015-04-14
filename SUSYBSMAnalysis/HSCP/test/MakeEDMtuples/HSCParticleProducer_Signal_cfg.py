import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = True
isBckg = False
isData = False
isSkimmedSample = False
GTAG = 'START72_V1::All'
InputFileList = cms.untracked.vstring()

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
   'root://cmseos:1094//eos/uscms/store/user/aackert/HSCP/AODgen/HSCPgluinoM1000_AOD-SIM.root',
#   'root://cmseos:1094//eos/uscms/store/user/aackert/HSCP/AODgen/HSCPgluinoM1500_AOD-SIM.root', 
#   'root://cmseos:1094//eos/uscms/store/user/aackert/HSCP/AODgen/HSCPgmstauM308_AOD-SIM.root', 
#   'root://cmseos:1094//eos/uscms/store/user/aackert/HSCP/AODgen/HSCPgmstauM494_AOD-SIM.root', 
#   'root://cmseos:1094//eos/uscms/store/user/aackert/HSCP/AODgen/HSCPstopM500_AOD-SIM.root', 
#   'root://cmseos:1094//eos/uscms/store/user/aackert/HSCP/AODgen/HSCPstopM900_AOD-SIM.root',
)


#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/BuildHSCParticles/HSCParticleProducer_cfg.py') )
