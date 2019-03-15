# EgammaAnalysis-TnPTreeProducer
TnP package for EGM
With modifications for HEEP7

Full recipe:
cmsrel CMSSW_10_2_5
cd CMSSW_10_2_5/src
cmsenv

git clone -b 10_2_5_tagAndProbeHEEP7  https://github.com/sethcooper/cmssw.git EgammaAnalysis/TnPTreeProducer

scram b -j8
