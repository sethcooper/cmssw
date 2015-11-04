#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor
import glob

print("compile the Stability code")
os.system("sh " + os.getcwd() + "/StabilityCheck.sh ") #just compile

print 'STABILITY'
FarmDirectory = "FARM"
JobName = "HSCPStability"
LaunchOnCondor.Jobs_RunHere = 1
LaunchOnCondor.Jobs_Queue = "8nh"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
#LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StabilityCheck.C", '"ANALYSE"'])
#LaunchOnCondor.SendCluster_Push(["BASH", "sh " + os.getcwd()+"/StabilityCheck.sh " + os.getcwd()+"/pictures"])

NJobs = 100
for Job in range(0,NJobs) :
      LaunchOnCondor.SendCluster_Push(["BASH", "sh " + os.getcwd()+"/StabilityCheck.sh " + os.getcwd()+"/pictures " + str(Job) +" " + str(NJobs)])
LaunchOnCondor.SendCluster_Submit()
