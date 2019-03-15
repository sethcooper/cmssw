from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys

config = config()

config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/user/s/scooper/work/private/cmssw/1025/LQElectronTagAndProbe/src/EgammaAnalysis/TnPTreeProducer/python/TnPTreeProducer_cfg.py'
#config.JobType.sendExternalFolder = True
config.JobType.pyCfgParams = ['isMC=False','doEleID=True','doPhoID=False','doTrigger=True','GT=102X_dataRun2_v8']

config.Data.outLFNDirBase = '/store/user/scooper/LQ/TagAndProbe/'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
#config.Data.unitsPerJob = 10000
#config.Data.splitting = 'EventAwareLumiBased'

config.General.requestName = 'Run2017B_31Mar2018'
config.Data.inputDataset = '/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD'
config.Data.publication = False

config.Site.storageSite = 'T2_CH_CERN'
