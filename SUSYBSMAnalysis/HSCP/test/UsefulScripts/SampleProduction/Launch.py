#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob

samples = [
#SampleName, generator config file, NJobs, NEvents/Job
  ['ZMM' , 'ZMM_13TeV_TuneCUETP8M1_cfi', 200, 500],
  ['Gluino400' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_400_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 50],
  ['Gluino1000', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_1000_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 50],
  ['Gluino1600', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_1600_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 50],
  ['Stop400' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_400_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 50],
  ['Stop1000' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_1000_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 50],
  ['Stop1600' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_1600_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 50],
  ['GMStau308', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_494_TuneZ2star_13TeV_pythia6_cff.py', 100, 50],
  ['GMStau494', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_494_TuneZ2star_13TeV_pythia6_cff.py', 100, 50],
  ['GMStau871', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_871_TuneZ2star_13TeV_pythia6_cff.py', 100, 50],
]


if sys.argv[1]=='1':
   #this is common to all samples, so we need to do it only once
   os.system('cmsDriver.py --filein file:step1.root --fileout XXX_OUTPUT_XXX_XXX_I_XXX.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions MCRUN2_74_V7 --step HLT:GRun,RAW2DIGI,L1Reco,RECO --python_filename RECO_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n -1')

   for S in samples:
      if('HSCP' in S[1]):
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions MCRUN2_74_V7 --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 10')
      else:
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions MCRUN2_74_V7 --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 10')

      with open("GEN_SIM_Template_cfg.py", "a") as f:
         f.write('process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
         f.write('process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
         f.write('process.maxEvents.input = cms.untracked.int32(XXX_NEVENTS_XXX)\n')
#         f.write('process.RAWSIMoutput.fileName = cms.untracked.string("file:XXX_OUTPUT_XXX_XXX_I_XXX.root")\n')
         f.write('process.source.firstLuminosityBlock =  cms.untracked.uint32(1XXX_I_XXX)\n')
         f.close()

      JobName = S[0]+"_SIMAOD"
      FarmDirectory = "FARM_"+JobName
      LaunchOnCondor.Jobs_NEvent = S[3]
      LaunchOnCondor.Jobs_Skip = 0
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
      for i in range(0,S[2]):
           LaunchOnCondor.Jobs_Count = i
           LaunchOnCondor.Jobs_Skip+=LaunchOnCondor.Jobs_NEvent
           LaunchOnCondor.SendCluster_Push  (["CMSSW", ["GEN_SIM_Template_cfg.py", "RECO_Template_cfg.py"] ])
           LaunchOnCondor.Jobs_FinalCmds = ['rm step1.root']
      LaunchOnCondor.SendCluster_Submit()




elif sys.argv[1]=='2':
   for S in samples:
      JobName = S[0]+"_EDM"
      FarmDirectory = "FARM_EDM"
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

      f = open("HSCPEDM_cfg.py", "w")
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
      f.write("GTAG = 'MCRUN2_74_V7::All'\n")
      f.write("OUTPUTFILE = 'XXX_OUTPUT_XXX.root'\n")
      f.write("InputFileList = cms.untracked.vstring()\n")
      f.write("\n")
      for i in range(0,S[2]):
         f.write("InputFileList.extend(['file:"+os.getcwd()+"/FARM_"+S[0]+"_SIMAOD/outputs/"+S[0]+"_SIMAOD_%04i.root'])\n" % i)
      f.write("\n")
      f.write("#main EDM tuple cfg that depends on the above parameters\n")
      f.write("execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )\n")
      f.close()

      LaunchOnCondor.SendCluster_Push  (["CMSSW", "HSCPEDM_cfg.py" ])
      LaunchOnCondor.SendCluster_Submit()
