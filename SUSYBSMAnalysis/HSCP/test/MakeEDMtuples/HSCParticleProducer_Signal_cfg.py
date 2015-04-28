import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = True
isBckg = False
isData = False
isSkimmedSample = False
GTAG = 'MCRUN2_74_V7::All'
InputFileList = cms.untracked.vstring()

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgluinoM1000_AOD-SIM.root',
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgluinoM1500_AOD-SIM.root',
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPstopM500_AOD-SIM.root', 
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPstopM900_AOD-SIM.root',
#    'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgmstauM308_AOD-SIM.root',
   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgmstauM494_AOD-SIM.root',
)


#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )
