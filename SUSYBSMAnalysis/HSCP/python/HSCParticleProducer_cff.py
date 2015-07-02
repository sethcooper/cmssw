import FWCore.ParameterSet.Config as cms


####################################################################################
#   Save Isolation info around tracks
####################################################################################

from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *
from TrackingTools.TrackAssociator.default_cfi import *

TrackAssociatorParametersForHSCPIsol = TrackAssociatorParameterBlock.TrackAssociatorParameters.clone()
TrackAssociatorParametersForHSCPIsol.useHO = cms.bool(False)
TrackAssociatorParametersForHSCPIsol.CSCSegmentCollectionLabel     = cms.InputTag("cscSegments")
TrackAssociatorParametersForHSCPIsol.DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments")
TrackAssociatorParametersForHSCPIsol.EERecHitCollectionLabel       = cms.InputTag("reducedEcalRecHitsEE")
TrackAssociatorParametersForHSCPIsol.EBRecHitCollectionLabel       = cms.InputTag("reducedEcalRecHitsEB")
TrackAssociatorParametersForHSCPIsol.HBHERecHitCollectionLabel     = cms.InputTag("reducedHcalRecHits", "hbhereco")
TrackAssociatorParametersForHSCPIsol.HORecHitCollectionLabel       = cms.InputTag("reducedHcalRecHits", "horeco")

HSCPIsolation01 = cms.EDProducer("ProduceIsolationMap",
      inputCollection  = cms.InputTag("generalTracks"),
      IsolationConeDR  = cms.double(0.1),
      TkIsolationPtCut = cms.double(10),
      TKLabel          = cms.InputTag("generalTracks"),
      TrackAssociatorParameters=TrackAssociatorParametersForHSCPIsol,
)

HSCPIsolation03 = HSCPIsolation01.clone()
HSCPIsolation03.IsolationConeDR  = cms.double(0.3)

HSCPIsolation05 = HSCPIsolation01.clone()
HSCPIsolation05.IsolationConeDR  = cms.double(0.5)


####################################################################################
#   Save muon segments in a compressed format
####################################################################################

MuonSegmentProducer = cms.EDProducer("MuonSegmentProducer",
   CSCSegments        = cms.InputTag("cscSegments"),
   DTSegments         = cms.InputTag("dt4DSegments"),
)


####################################################################################
#   HSCParticle Producer
####################################################################################

#ALL THIS IS NEEDED BY ECAL BETA CALCULATOR (TrackAssociator)
from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *
from TrackingTools.TrackAssociator.default_cfi import * 
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi import *
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *

from SUSYBSMAnalysis.HSCP.HSCPSelections_cff import *
HSCParticleProducer = cms.EDFilter("HSCParticleProducer",
   TrackAssociatorParameterBlock, #Needed for ECAL/Track Matching

   #DOES THE PRODUCER ACT AS AN EDFILTER?
   filter = cms.bool(True),

   #WHAT (BETA) INFORMATION TO COMPUTE
   useBetaFromTk      = cms.bool(True),  #does nothing because we have all the info saved in EDM format
   useBetaFromMuon    = cms.bool(True),  #does nothing because we have all the info saved in EDM format
   useBetaFromRpc     = cms.bool(False), #must be updated for AOD
   useBetaFromEcal    = cms.bool(False), #must be updated for AOD and 74X

   #TAG OF THE REQUIRED INPUT COLLECTION (ONLY ACTIVATED CALCULATOR)
   tracks             = cms.InputTag("generalTracks"),
   tracksIsolation    = cms.InputTag("generalTracks"),
   muons              = cms.InputTag("muons"),
   MTmuons            = cms.InputTag("muons"),
   EBRecHitCollection = cms.InputTag("ecalRecHit:EcalRecHitsEB"),
   EERecHitCollection = cms.InputTag("ecalRecHit:EcalRecHitsEE"),
   rpcRecHits         = cms.InputTag("rpcRecHits"),

   #TRACK SELECTION FOR THE HSCP SEED
   minMuP             = cms.double(20),
   minTkP             = cms.double(20),
   maxTkChi2          = cms.double(20),
   minTkHits          = cms.uint32(3),
   minSAMuPt          = cms.double(70),
   minMTMuPt          = cms.double(70),

   #MUON/TRACK MATCHING THRESHOLDS (ONLY IF NO MUON INNER TRACK)
   minDR              = cms.double(0.1),
   maxInvPtDiff       = cms.double(0.005),
   minMTDR              = cms.double(0.3),

   #SELECTION ON THE PRODUCED HSCP CANDIDATES (WILL STORE ONLY INTERESTING CANDIDATES)
   SelectionParameters = cms.VPSet(
      HSCPSelectionDefault,
      HSCPSelectionMTMuonOnly,
      HSCPSelectionSAMuonOnly,
   ),
)

####################################################################################
#   HSCParticle Selector  (Just an Example of what we can do)
####################################################################################

HSCParticleSelector = cms.EDFilter("HSCParticleSelector",
   source = cms.InputTag("HSCParticleProducer"),
   filter = cms.bool(True),

   SelectionParameters = cms.VPSet(
      HSCPSelectionHighdEdx, #THE OR OF THE TWO SELECTION WILL BE APPLIED
      HSCPSelectionHighTOF,
   ),
)

####################################################################################
#   HSCP Candidate Sequence
####################################################################################

HSCParticleProducerSeq = cms.Sequence(HSCPIsolation01 * HSCPIsolation03 * HSCPIsolation05 * MuonSegmentProducer * HSCParticleProducer)

