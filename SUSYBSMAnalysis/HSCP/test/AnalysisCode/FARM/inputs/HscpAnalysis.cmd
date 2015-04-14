####################################
#        LaunchOnFarm Script       #
#     Loic.quertenmont@cern.ch     #
#            April 2010            #
####################################

bsub -q 8nh -J HscpAnalysis0000_ '/afs/cern.ch/work/q/querten/public/15_03_12_HSCP_Run2Preparation/CMSSW_7_4_0/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/FARM/inputs/0000_HscpAnalysis.sh 0 ele'
