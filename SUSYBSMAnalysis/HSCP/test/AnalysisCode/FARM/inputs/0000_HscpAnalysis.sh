#! /bin/sh
####################################
#        LaunchOnFarm Script       #
#     Loic.quertenmont@cern.ch     #
#            April 2010            #
####################################

export SCRAM_ARCH=slc6_amd64_gcc491
export BUILD_ARCH=slc5_amd64_gcc462
export VO_CMS_SW_DIR=/nfs/soft/cms
cd /afs/cern.ch/work/q/querten/public/15_03_12_HSCP_Run2Preparation/CMSSW_7_4_0/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis
eval `scramv1 runtime -sh`
root -l -b << EOF
   TString makeshared(gSystem->GetMakeSharedLib());
   makeshared.ReplaceAll("-W ", "-Wno-deprecated-declarations -Wno-deprecated -Wno-unused-local-typedefs -Wno-attributes ");
   makeshared.ReplaceAll("-Woverloaded-virtual ", " ");
   makeshared.ReplaceAll("-Wshadow ", " -std=c++0x -D__USE_XOPEN2K8 ");
   cout << "Compilling with the following arguments: " << makeshared << endl;
   gSystem->SetMakeSharedLib(makeshared);
   gSystem->SetIncludePath("-I$ROOFITSYS/include");
   gSystem->Load("libFWCoreFWLite");
   AutoLibraryLoader::enable();
   gSystem->Load("libDataFormatsFWLite.so");
   gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
   gSystem->Load("libDataFormatsVertexReco.so");
   gSystem->Load("libDataFormatsHepMCCandidate.so");
   gSystem->Load("libPhysicsToolsUtilities.so");
   gSystem->Load("libdcap.so");
   .x /afs/cern.ch/work/q/querten/public/15_03_12_HSCP_Run2Preparation/CMSSW_7_4_0/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/Analysis_Step3.C+("ANALYSE_0_to_0", 0)
   .q
EOF

mv HscpAnalysis* /afs/cern.ch/work/q/querten/public/15_03_12_HSCP_Run2Preparation/CMSSW_7_4_0/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/FARM/outputs/
