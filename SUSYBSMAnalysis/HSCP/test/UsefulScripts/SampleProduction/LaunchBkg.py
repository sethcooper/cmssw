#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob

"""
check that a file exist and is not corrupted
"""
def checkInputFile(url):
    if(url.startswith('/store')==True):
       url= 'root://eoscms//eos/cms'+url
    command_out = commands.getstatusoutput("root -l -b -q " + url)
    if(command_out[1].find("Error")>=0 or command_out[1].find("probably not closed")>=0 or command_out[1].find("Corrupted")>=0):return False
    return True


samples = [
#SampleName, generator config file, NJobs, NEvents/Job
  ['MC_13TeV_DYToMuMu' , 'ZMM_13TeV_TuneCUETP8M1_cfi', 4000, 2500],
  ['MC_13TeV_WToLNu' , 'WToLNu_13TeV_pythia8_cff', 4000, 2500],
]


if sys.argv[1]=='1':  #GEN-SIM-DIGI-RECO-AOD in one shot
   #this is common to all samples, so we need to do it only once
   os.system('cmsDriver.py --filein file:step1.root --fileout step2.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions auto:run2_mc --step HLT:GRun,RAW2DIGI,L1Reco,RECO --python_filename RECO_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n -1')

   for S in samples:
      if('HSCP' in S[1]):
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions auto:run2_mc --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 10')
      else:
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions auto:run2_mc --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 10')

      with open("GEN_SIM_Template_cfg.py", "a") as f:
         f.write('process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
         f.write('process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
         f.write('process.maxEvents.input = cms.untracked.int32(XXX_NEVENTS_XXX)\n')
#         f.write('process.RAWSIMoutput.fileName = cms.untracked.string("file:XXX_OUTPUT_XXX_XXX_I_XXX.root")\n')
         f.write('process.source.firstLuminosityBlock =  cms.untracked.uint32(1XXX_I_XXX)\n')
         f.close()


      with open("HSCPEDM_Template_cfg.py", "w") as f:
         f.write("import sys, os\n")
         f.write("import FWCore.ParameterSet.Config as cms\n")
         f.write("\n")
         if('HSCP' in S[1]):
            f.write("isSignal = True\n")
            f.write("isBckg = False\n")
         else:
            f.write("isSignal = False\n")
            f.write("isBckg = True\n")
         f.write("isData = False\n")
         f.write("isSkimmedSample = False\n")
         f.write("GTAG = 'MCRUN2_74_V8'\n")
         f.write("OUTPUTFILE = 'XXX_OUTPUT_XXX_XXX_I_XXX.root'\n")
         f.write("InputFileList = cms.untracked.vstring()\n")
         f.write("\n")
         f.write("InputFileList.extend(['file:"+"step2.root"+"'])\n")

         f.write("\n")
         f.write("#main EDM tuple cfg that depends on the above parameters\n")
         f.write("execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )\n")
         f.close()


      JobName = S[0]+"_SIMEDM"
      FarmDirectory = "FARM_"+JobName
      LaunchOnCondor.Jobs_NEvent = S[3]
      LaunchOnCondor.Jobs_Skip = 0
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
      for i in range(0,S[2]):
           LaunchOnCondor.Jobs_Count = i
           LaunchOnCondor.Jobs_Skip+=LaunchOnCondor.Jobs_NEvent
           LaunchOnCondor.SendCluster_Push  (["CMSSW", ["GEN_SIM_Template_cfg.py", "RECO_Template_cfg.py", "HSCPEDM_Template_cfg.py"] ])
           LaunchOnCondor.Jobs_FinalCmds = ['rm step1.root; rm step2.root']
      LaunchOnCondor.SendCluster_Submit()


elif sys.argv[1]=='2':  #MergeAll EDM files into one
   for S in samples:
      JobName = S[0]
      FarmDirectory = "FARM_"+JobName
      InputFiles = LaunchOnCondor.GetListOfFiles('"file:',os.getcwd()+"/FARM_"+S[0]+"_SIMEDM/outputs/*.root",'"')
      LaunchOnCondor.SendCMSMergeJob(FarmDirectory, JobName, InputFiles,  '"'+S[0]+'.root"', '"keep *"')



elif sys.argv[1]=='3':  #Transfert final EDM files from your place to CERN CMST3
   for S in samples:
      os.system('lcg-cp --verbose -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN=/storage/data/cms/store/user/quertenmont/15_05_20_HSCP_AODSIM/'+S[0]+'.root" "srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms//store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/'+S[0]+'.root" &> LOG_'+S[0]+'.log &')
