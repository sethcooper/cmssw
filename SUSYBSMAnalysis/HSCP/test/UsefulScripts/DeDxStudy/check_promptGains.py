#!/bin/env/python

import urllib
import string
import os,sys,time
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob
import commands
import json
import collections # kind of map

prompt_gains      = "conddb_payloads.txt"
gains             = "/afs/cern.ch/cms/tracker/sistrvalidation/WWW/CalibrationValidation/ParticleGain"
sql3_command      = "sqlite3 -line"
gains_association = "gains_association.txt"
gainsDirs         = os.listdir(gains)
prompt_payloads   = []
prompt_runs       = []


print "reading the coddb, please wait ..."
os.system("conddb list SiStripApvGain_FromParticles_GR10_v1_express | grep SiStripApvGain > "+prompt_gains)
f = open(prompt_gains, 'r')
while True:
   line = f.readline().split()
   if len(line) == 0:
      break
   else:
      prompt_payloads = prompt_payloads + [line[3]]
      prompt_runs     = prompt_runs     + [line[0]]
f.close()

payload = ""
gainsDirs.sort()
os.system("rm -f "+gains_association)
fout = open(gains_association, "w")
listOfAppliedPromptRuns = []
listOfRuns = []
completeSplitting = []
listOfPromptGainDirs = []
listOfGainDirs = []
print "File ",gains_association," has been created!"
for i in range(0,len(prompt_payloads)):
   for gainsDir in gainsDirs:
      line = os.popen("sqlite3 -line "+gains+"/"+gainsDir+"/sqlite/Gains_Sqlite.db 'select * from IOV;' | grep -i payload").read().split()
      if len(line) == 0:
         continue
      elif prompt_payloads[i] == line[2]:
         if len(listOfAppliedPromptRuns) != 0:
            fout.write("_to_"+prompt_runs[i]+"\n")
         print "Match found!"
#         fout.write("Prompt payload of "+gainsDir+" is "+prompt_payloads[i]+" belonging to run "+prompt_runs[i]+'\n')
         fout.write(gainsDir+" used by Run_"+prompt_runs[i])
         listOfAppliedPromptRuns = listOfAppliedPromptRuns + [prompt_runs[i]]
         listOfPromptGainDirs = listOfPromptGainDirs + [gainsDir]
#fout.close()
listOfAppliedPromptRuns = listOfAppliedPromptRuns + ['999999']
fout.write("_to_999999\n")
fout.write("\n-----------------------------\n")
fout.write("\n")

#passedFirstOne=False
#for gainsDir in gainsDirs:
#   runs = gainsDir.split("_")
#   if len(runs) >= 4 and len(runs) < 8 and int(runs[1]) >= int(listOfAppliedPromptRuns[0]):
#      passedFirstOne=True
#   if passedFirstOne:
#      listOfRuns = listOfRuns + [runs[1]] + [runs[3]]
#      listOfGainDirs = listOfGainDirs + [gainsDir]
#   else: continue

print "PromptRuns are: ",listOfAppliedPromptRuns
print "They belong to: ",listOfPromptGainDirs
#print "All runs are: ",listOfRuns
#completeSplitting = listOfAppliedPromptRuns+listOfRuns
#completeSplitting.sort()
#print "All intervals are: ", completeSplitting
#print "This means we have ",len(completeSplitting)," intervals to consider"

#for i in range(1,len(completeSplitting)):
#   fout.write("Interval Run_"+completeSplitting[i-1]+"_to_"+completeSplitting[i]+"\n")
#fout.write("Interval Run_"+completeSplitting[len(completeSplitting)-1]+"_to_999999\n")
#fout.close()

fout.write("\n")
fout.write("\n-----------------------------\n")
fout.write("\n")

for gainsDir in gainsDirs:
   j = 0
   runs = gainsDir.split("_")
   if len(runs) >= 4 and len(runs) < 8:
      while j < len(listOfAppliedPromptRuns)-1:
         print "j = ",j
         CurrentPromptRuns = [listOfAppliedPromptRuns[j], listOfAppliedPromptRuns[j+1]]
         if int(runs[1]) >= int(CurrentPromptRuns[0]) and int(runs[3]) <= int(CurrentPromptRuns[1]):
            fout.write(gainsDir+' => Run_'+runs[1]+"_to_"+runs[3]+'='+' "'+gainsDir+'" / "'+listOfPromptGainDirs[j]+'"\n')
            j = j+1
         elif int(runs[1]) >= int(CurrentPromptRuns[0]) and int(runs[3]) > int(CurrentPromptRuns[1]):
            NextPromptRuns = [listOfAppliedPromptRuns[j+1], listOfAppliedPromptRuns[j+2]]
            fout.write(gainsDir+' => Run_'+runs[1]+"_to_"+CurrentPromptRuns[1]+'='+' "'+gainsDir+'" / "'+listOfPromptGainDirs[j]+'"\n')
            fout.write(gainsDir+' => Run_'+NextPromptRuns[0]+"_to_"+runs[3]+'='+' "'+gainsDir+'" / "'+listOfPromptGainDirs[j+1]+'"\n')
            j = j+1
         else:
            j = j+1
   else: continue

fout.close()
