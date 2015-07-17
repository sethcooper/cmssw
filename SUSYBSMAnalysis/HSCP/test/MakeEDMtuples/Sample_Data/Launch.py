#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor 
import glob


#the vector below contains the "TypeMode" of the analyses that should be run
#SamplesToMerge = [1,2,3]#,3,5]


#if len(sys.argv)==1:       
#	print "Please pass in argument a number between 1 and 5"
#        print "  1  - Merge SingleMu DoubleMu and MET"
#        print "  2  - Merge only SingleMu "
#        print "  3  - Merge only DoubleMu "
#	sys.exit()


CMSSW_VERSION = os.getenv('CMSSW_VERSION','CMSSW_VERSION')
if CMSSW_VERSION == 'CMSSW_VERSION':
  print 'please setup your CMSSW environement'
  sys.exit(0)

JobName = "MergeFiles"
FarmDirectory = "FARM"
LaunchOnCondor.Jobs_RunHere = 1
LaunchOnCondor.SendCluster_Create (FarmDirectory, JobName)
LaunchOnCondor.SendCluster_Push (["CMSSW", "Merge_cfg.py"])
LaunchOnCondor.SendCLuster_Submit ()

