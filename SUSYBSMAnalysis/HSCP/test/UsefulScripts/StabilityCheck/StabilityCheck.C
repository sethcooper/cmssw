
#include <exception>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TPaveText.h"


namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;


#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"

#endif


bool PassPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj, const reco::DeDxData& dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev);
bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold=false);


bool PassPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj, const reco::DeDxData& dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev)
{
   if(TypeMode==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;
   if(TypeMode==2 && hscp.type() != HSCParticleType::globalMuon)return false;
   reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;

   if(fabs(track->eta())>GlobalMaxEta) return false;
   if(track->found()<GlobalMinNOH)return false;
   if(track->hitPattern().numberOfValidPixelHits()<2)return false; 
   if(dedxSObj.numberOfMeasurements()<GlobalMinNOM)return false;
//   if(tof && tof->nDof()<GlobalMinNDOF && (dttof->nDof()<GlobalMinNDOFDT || csctof->nDof()<GlobalMinNDOFCSC) )return false;

   if(track->qualityMask()<GlobalMinQual )return false;
   if(track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(track->pt()<GlobalMinPt)return false;
   if(dedxSObj.dEdx()<GlobalMinIs)return false;
   if(dedxMObj.dEdx()<GlobalMinIm)return false;
//   if(tof && tof->inverseBeta()<GlobalMinTOF)return false;
//   if(tof && tof->inverseBetaErr()>GlobalMaxTOFErr)return false;

   fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
   vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
   if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return false;}
   const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return false;}

   double dz  = track->dz (vertexColl[0].position());
   double dxy = track->dxy(vertexColl[0].position());
   for(unsigned int i=1;i<vertexColl.size();i++){
      if(fabs(track->dz (vertexColl[i].position())) < fabs(dz) ){
         dz  = track->dz (vertexColl[i].position());
         dxy = track->dxy(vertexColl[i].position());
      }
   }
   double v3d = sqrt(dz*dz+dxy*dxy);
   if(v3d>GlobalMaxV3D )return false;

   fwlite::Handle<HSCPIsolationValueMap> IsolationH;
   IsolationH.getByLabel(ev, "HSCPIsolation", "R03"); //New format used for data since 17-07-2015
   if(!IsolationH.isValid()){
     IsolationH.getByLabel(ev, "HSCPIsolation03");//Old format used for first 2015B data, Signal and MC Backgrounds
     if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
   }
  const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();


   HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
    if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;

   double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
   if(EoP>GlobalMaxEIsol)return false;

   if((track->ptError()/track->pt())>GlobalMaxPterr)return false;
   if(std::max(0.0,track->pt() - track->ptError())<GlobalMinPt)return false;
   return true;
}


bool PassingTrigger(const fwlite::ChainEvent& ev, const std::string& TriggerName){
   edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
   if(!tr.isValid())         tr = ev.triggerResultsByName("MergeHLT");
   if(!tr.isValid())return false;

   if(TriggerName=="Any"){
      if(passTriggerPatterns(tr, "HLT_PFMET170_NoiseCleaned_v*"))return true;
      if(passTriggerPatterns(tr, "HLT_Mu45_eta2p1_v*"))return true;
      if(passTriggerPatterns(tr, "HLT_Mu50_v*"))return true;
   }else{
      if(passTriggerPatterns(tr, (TriggerName + "_v*").c_str()))return true;
   }
   return false;
}




void StabilityCheck(string DIRNAME="COMPILE", string OUTDIRNAME="pictures")
{
  printf("DIRNAME = %s\n", DIRNAME.c_str());
  if(DIRNAME=="COMPILE") return;
  OUTDIRNAME+="/";

   Event_Weight = 1;
   MaxEntry = -1;

   system((string("mkdir -p ") + OUTDIRNAME).c_str());

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.07);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505,"X");
   TH1::AddDirectory(kTRUE);

   std::map<unsigned int, unsigned int> RunBinIndex;
   unsigned int NextIndex=0;

   InitBaseDirectory();
   GetSampleDefinition(samples , DIRNAME+"/../../AnalysisCode/Analysis_Samples.txt");
   keepOnlyValidSamples(samples);
   vector<string> DataFileName;
   for(unsigned int s=0;s<samples.size();s++){
      if(samples[s].Type!=0)continue; //only data sample is considered
      GetInputFiles(samples[s], BaseDirectory, DataFileName, 0);
   }
//   DataFileName.clear();
//   DataFileName.push_back("root://eoscms//eos/cms/store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/Data_Run2015B_251253_251883.root");

   for(unsigned int f=0;f<DataFileName.size();f++){printf("file %i : %s\n", f, DataFileName[f].c_str());}

   std::vector<string> triggers;
   triggers.push_back("Any");
   triggers.push_back("HLT_Mu45_eta2p1");
   triggers.push_back("HLT_Mu50");
   triggers.push_back("HLT_PFMET170_NoiseCleaned");

   TProfile** NVertProf = new TProfile*[triggers.size()];
   TProfile** dEdxProf = new TProfile*[triggers.size()];
   TProfile** dEdxMProf = new TProfile*[triggers.size()];
   TProfile** dEdxMSProf = new TProfile*[triggers.size()];
   TProfile** dEdxMPProf = new TProfile*[triggers.size()];
   TProfile** dEdxMSCProf = new TProfile*[triggers.size()];
   TProfile** dEdxMPCProf = new TProfile*[triggers.size()];
   TProfile** dEdxMSFProf = new TProfile*[triggers.size()];
   TProfile** dEdxMPFProf = new TProfile*[triggers.size()];
   TProfile** PtProf   = new TProfile*[triggers.size()];
   TProfile** TOFAODProf   = new TProfile*[triggers.size()];
   TProfile** TOFAODDTProf   = new TProfile*[triggers.size()];
   TProfile** TOFAODCSCProf   = new TProfile*[triggers.size()];
   TProfile** TOFAODOverMinProf   = new TProfile*[triggers.size()];
   TProfile** TOFAODDTOverMinProf   = new TProfile*[triggers.size()];
   TProfile** TOFAODCSCOverMinProf   = new TProfile*[triggers.size()];
   TProfile** VertexAODProf   = new TProfile*[triggers.size()];
   TProfile** VertexAODDTProf   = new TProfile*[triggers.size()];
   TProfile** VertexAODCSCProf   = new TProfile*[triggers.size()];
   TProfile** TOFProf   = new TProfile*[triggers.size()];
   TProfile** TOFDTProf   = new TProfile*[triggers.size()];
   TProfile** TOFCSCProf   = new TProfile*[triggers.size()];
   TProfile** TOFOverMinProf   = new TProfile*[triggers.size()];
   TProfile** TOFDTOverMinProf   = new TProfile*[triggers.size()];
   TProfile** TOFCSCOverMinProf   = new TProfile*[triggers.size()];
   TProfile** VertexProf   = new TProfile*[triggers.size()];
   TProfile** VertexDTProf   = new TProfile*[triggers.size()];
   TProfile** VertexCSCProf   = new TProfile*[triggers.size()];
   TH1D**     Count    = new TH1D*    [triggers.size()];
   TH1D**     CountMuAOD  = new TH1D*    [triggers.size()];
   TH1D**     CountMu  = new TH1D*    [triggers.size()];
   TH1D**     HdEdx    = new TH1D*    [triggers.size()];
   TH1D**     HPt      = new TH1D*    [triggers.size()];
   TH1D**     HTOFAOD  = new TH1D*    [triggers.size()];
   TH1D**     HTOF      = new TH1D*    [triggers.size()];


   TFile* OutputHisto = new TFile((OUTDIRNAME + "/Histos.root").c_str(),"RECREATE");
   for(unsigned int i=0;i<triggers.size();i++){
      NVertProf[i] = new TProfile((triggers[i] + "NVertProf").c_str(), "NVertProf", 10000 ,0, 10000);
      dEdxProf[i] = new TProfile((triggers[i] + "dEdxProf").c_str(), "dEdxProf", 10000 ,0, 10000);
      dEdxMProf[i] = new TProfile((triggers[i] + "dEdxMProf").c_str(), "dEdxMProf", 10000 ,0, 10000);
      dEdxMSProf[i] = new TProfile((triggers[i] + "dEdxMSProf").c_str(), "dEdxMSProf", 10000 ,0, 10000);
      dEdxMPProf[i] = new TProfile((triggers[i] + "dEdxMPProf").c_str(), "dEdxMPProf", 10000 ,0, 10000);
      dEdxMSCProf[i] = new TProfile((triggers[i] + "dEdxMSCProf").c_str(), "dEdxMSCProf", 10000 ,0, 10000);
      dEdxMPCProf[i] = new TProfile((triggers[i] + "dEdxMPCProf").c_str(), "dEdxMPCProf", 10000 ,0, 10000);
      dEdxMSFProf[i] = new TProfile((triggers[i] + "dEdxMSFProf").c_str(), "dEdxMSFProf", 10000 ,0, 10000);
      dEdxMPFProf[i] = new TProfile((triggers[i] + "dEdxMPFProf").c_str(), "dEdxMPFProf", 10000 ,0, 10000);

      PtProf  [i] = new TProfile((triggers[i] + "PtProf"  ).c_str(), "PtProf"  , 10000 ,0, 10000);
      TOFAODProf  [i] = new TProfile((triggers[i] + "TOFAODProf"  ).c_str(), "TOFAODProf"  , 10000 ,0, 10000);
      TOFAODDTProf  [i] = new TProfile((triggers[i] + "TOFAODDTProf"  ).c_str(), "TOFAODDTProf"  , 10000 ,0, 10000);
      TOFAODCSCProf  [i] = new TProfile((triggers[i] + "TOFAODCSCProf"  ).c_str(), "TOFAODCSCProf"  , 10000 ,0, 10000);

      TOFAODOverMinProf  [i] = new TProfile((triggers[i] + "TOFAODOverMinProf"  ).c_str(), "TOFAODOverMinProf"  , 10000 ,0, 10000);
      TOFAODDTOverMinProf  [i] = new TProfile((triggers[i] + "TOFAODDTOverMinProf"  ).c_str(), "TOFAODDTOverMinProf"  , 10000 ,0, 10000);
      TOFAODCSCOverMinProf  [i] = new TProfile((triggers[i] + "TOFAODCSCOverMinProf"  ).c_str(), "TOFAODCSCOverMinProf"  , 10000 ,0, 10000);

      VertexAODProf  [i] = new TProfile((triggers[i] + "VertexAODProf"  ).c_str(), "VertexAODProf"  , 10000 ,0, 10000);
      VertexAODDTProf  [i] = new TProfile((triggers[i] + "VertexAODDTProf"  ).c_str(), "VertexAODDTProf"  , 10000 ,0, 10000);
      VertexAODCSCProf  [i] = new TProfile((triggers[i] + "VertexAODCSCProf"  ).c_str(), "VertexAODCSCProf"  , 10000 ,0, 10000);


      TOFProf  [i] = new TProfile((triggers[i] + "TOFProf"  ).c_str(), "TOFProf"  , 10000 ,0, 10000);
      TOFDTProf  [i] = new TProfile((triggers[i] + "TOFDTProf"  ).c_str(), "TOFDTProf"  , 10000 ,0, 10000);
      TOFCSCProf  [i] = new TProfile((triggers[i] + "TOFCSCProf"  ).c_str(), "TOFCSCProf"  , 10000 ,0, 10000);

      TOFOverMinProf  [i] = new TProfile((triggers[i] + "TOFOverMinProf"  ).c_str(), "TOFOverMinProf"  , 10000 ,0, 10000);
      TOFDTOverMinProf  [i] = new TProfile((triggers[i] + "TOFDTOverMinProf"  ).c_str(), "TOFDTOverMinProf"  , 10000 ,0, 10000);
      TOFCSCOverMinProf  [i] = new TProfile((triggers[i] + "TOFCSCOverMinProf"  ).c_str(), "TOFCSCOverMinProf"  , 10000 ,0, 10000);

      VertexProf  [i] = new TProfile((triggers[i] + "VertexProf"  ).c_str(), "VertexProf"  , 10000 ,0, 10000);
      VertexDTProf  [i] = new TProfile((triggers[i] + "VertexDTProf"  ).c_str(), "VertexDTProf"  , 10000 ,0, 10000);
      VertexCSCProf  [i] = new TProfile((triggers[i] + "VertexCSCProf"  ).c_str(), "VertexCSCProf"  , 10000 ,0, 10000);

      Count   [i] = new TH1D(    (triggers[i] + "Count"   ).c_str(), "Count"   , 10000 ,0, 10000);  Count  [i]->Sumw2();
      CountMuAOD[i] = new TH1D(    (triggers[i] + "CountMuAOD" ).c_str(), "CountMuAOD" , 10000 ,0, 10000);  CountMuAOD[i]->Sumw2();
      CountMu [i] = new TH1D(    (triggers[i] + "CountMu" ).c_str(), "CountMu" , 10000 ,0, 10000);  CountMu[i]->Sumw2();
      HdEdx   [i] = new TH1D(    (triggers[i] + "HdEdx"   ).c_str(), "HdEdx"   , 10000 ,0, 10000);  HdEdx  [i]->Sumw2();
      HPt     [i] = new TH1D(    (triggers[i] + "HPt"     ).c_str(), "HPt"     , 10000 ,0, 10000);  HPt    [i]->Sumw2();
      HTOFAOD  [i] = new TH1D(    (triggers[i] + "HTOFAOD"     ).c_str(), "HTOFAOD"     , 10000 ,0, 10000);  HTOFAOD [i]->Sumw2();
      HTOF     [i] = new TH1D(    (triggers[i] + "HTOF"     ).c_str(), "HTOF"     , 10000 ,0, 10000);  HTOF    [i]->Sumw2();
   }

   TypeMode      = 0;

   std::unordered_map<unsigned int,double> TrackerGains;
   double dEdxSF [2];
   dEdxSF [0] = 1.00000;
   dEdxSF [1] = 1.21836;
   TH3F* dEdxTemplates = loadDeDxTemplate("../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true);


   moduleGeom::loadGeometry("../../../data/CMS_GeomTree.root");
   muonTimingCalculator tofCalculator;
   tofCalculator.loadTimeOffset("../../../data/MuonTimeOffset.txt");
   unsigned int CurrentRun = 0;



   fwlite::ChainEvent ev(DataFileName);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Looping on Tree              :");
   int TreeStep = ev.size()/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t e=0;e<ev.size();e++){
      ev.to(e); 
      if(e%TreeStep==0){printf(".");fflush(stdout);}

      //if run change, update conditions
      if(CurrentRun != ev.eventAuxiliary().run()){
         CurrentRun = ev.eventAuxiliary().run();
         tofCalculator.setRun(CurrentRun);
      }


      if(RunBinIndex.find(ev.eventAuxiliary().run()) == RunBinIndex.end()){
         RunBinIndex[ev.eventAuxiliary().run()] = NextIndex;
         for(unsigned int i=0;i<triggers.size();i++){
            int Bin = HdEdx[i]->GetXaxis()->FindBin(NextIndex);
            char Label[2048]; sprintf(Label,"%6i",ev.eventAuxiliary().run());
            HdEdx[i]->GetXaxis()->SetBinLabel(Bin, Label);
            HPt[i]->GetXaxis()->SetBinLabel(Bin, Label);
            HTOFAOD[i]->GetXaxis()->SetBinLabel(Bin, Label);
            HTOF[i]->GetXaxis()->SetBinLabel(Bin, Label);
            Count[i]->GetXaxis()->SetBinLabel(Bin, Label);
            NVertProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxProf[i]->GetXaxis()->SetBinLabel(Bin, Label);      
            dEdxMProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxMSProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxMPProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxMSCProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxMPCProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxMSFProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            dEdxMPFProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            PtProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFAODProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFAODDTProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFAODCSCProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFAODOverMinProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFAODDTOverMinProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFAODCSCOverMinProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            VertexAODProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            VertexAODDTProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            VertexAODCSCProf[i]->GetXaxis()->SetBinLabel(Bin, Label);

            TOFProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFDTProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFCSCProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFOverMinProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFDTOverMinProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            TOFCSCOverMinProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            VertexProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            VertexDTProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
            VertexCSCProf[i]->GetXaxis()->SetBinLabel(Bin, Label);
         }
         NextIndex++;
      }
      unsigned int CurrentRunIndex = RunBinIndex[ev.eventAuxiliary().run()];

      if(!PassingTrigger(ev,"Any")){continue;} //need to pass at least one of the trigger, otherwise save time

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(ev,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      fwlite::Handle<DeDxHitInfoAss> dedxCollH;
      dedxCollH.getByLabel(ev, "dedxHitInfo");
      if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(ev, "muons",TOF_Label.c_str());
      if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}


      fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
      TOFDTCollH.getByLabel(ev, "muons",TOFdt_Label.c_str());
      if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
      TOFCSCCollH.getByLabel(ev, "muons",TOFcsc_Label.c_str());
      if(!TOFCSCCollH.isValid()){printf("Invalid CSCTOF collection\n");continue;}

      fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
      vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
      if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
      const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;

      fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
      fwlite::Handle<DTRecSegment4DCollection> DTSegmentCollHandle;            
      if(true){ //do not reocmpute TOF on MC background
         CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
         if(!CSCSegmentCollHandle.isValid()){printf("CSC Segment Collection not found!\n"); continue;}

         DTSegmentCollHandle.getByLabel(ev, "dt4DSegments");
         if(!DTSegmentCollHandle.isValid()){printf("DT Segment Collection not found!\n"); continue;}
      }



      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;
         reco::MuonRef muon = hscp.muonRef();         
         if(muon.isNull())continue; 
         if(!hscp.muonRef()->isStandAloneMuon())continue;

         const DeDxHitInfo* dedxHits = NULL;
         if(TypeMode!=3 && !track.isNull()) {
            DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());		 
            if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
         }

         bool useClusterCleaning = true;
         DeDxData dedxSObj = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode==5, false, TrackerGains.size()>0?&TrackerGains:NULL, true, true, 99, false, 1);
         DeDxData dedxMObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, TrackerGains.size()>0?&TrackerGains:NULL, true, true, 99, false, 1);
         DeDxData dedxMSObj = computedEdx(dedxHits, dEdxSF, NULL,          false,useClusterCleaning, false      , false, TrackerGains.size()>0?&TrackerGains:NULL, true, true, 99, false, 1);
         DeDxData dedxMPObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, TrackerGains.size()>0?&TrackerGains:NULL, false, true, 99, false, 1);

         const reco::MuonTimeExtra* tofaod = NULL;
         const reco::MuonTimeExtra* dttofaod = NULL;
         const reco::MuonTimeExtra* csctofaod = NULL;
         if(!hscp.muonRef().isNull()){ tofaod  = &TOFCollH->get(hscp.muonRef().key()); dttofaod  = &TOFDTCollH->get(hscp.muonRef().key()); csctofaod  = &TOFCSCCollH->get(hscp.muonRef().key());}

         const reco::MuonTimeExtra* tof = NULL;
         const reco::MuonTimeExtra* dttof = NULL;
         const reco::MuonTimeExtra* csctof = NULL;
         if(!hscp.muonRef().isNull() && hscp.muonRef()->isStandAloneMuon() ){
            const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
            const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
            tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, 1 ); //apply T0 correction on data but not on signal MC
            tof  = &tofCalculator.combinedTOF; dttof = &tofCalculator.dtTOF;  csctof = &tofCalculator.cscTOF;
         }          

         if(!PassPreselection(hscp, dedxSObj, dedxMObj, tof, dttof, csctof, ev)){continue;}

         for(unsigned int i=0;i<triggers.size();i++){
            if(!PassingTrigger(ev,triggers[i])){continue;}

            NVertProf[i]->Fill(CurrentRunIndex, vertexColl.size()); 


            if(tofaod && tofaod->nDof()>=GlobalMinNDOF && (dttofaod->nDof()>=GlobalMinNDOFDT || csctofaod->nDof()>=GlobalMinNDOFCSC) && tofaod->inverseBetaErr()<=GlobalMaxTOFErr){
               if(tofaod->inverseBeta()>=GlobalMinTOF)CountMuAOD[i]->Fill(CurrentRunIndex);
               if(tofaod->inverseBeta()>=GlobalMinTOF)TOFAODOverMinProf[i]->Fill(CurrentRunIndex, tofaod->inverseBeta());
               if(dttofaod->inverseBeta()>=GlobalMinTOF)TOFAODDTOverMinProf[i]->Fill(CurrentRunIndex, dttofaod->inverseBeta());
               if(csctofaod->inverseBeta()>=GlobalMinTOF)TOFAODCSCOverMinProf[i]->Fill(CurrentRunIndex, csctofaod->inverseBeta());
               TOFAODProf[i]->Fill(CurrentRunIndex, tofaod->inverseBeta());
               if(dttofaod->nDof()>=GlobalMinNDOFDT) TOFAODDTProf[i]->Fill(CurrentRunIndex, dttofaod->inverseBeta());
               if(csctofaod->nDof()>=GlobalMinNDOFCSC) TOFAODCSCProf[i]->Fill(CurrentRunIndex, csctofaod->inverseBeta());
               if(tofaod->inverseBeta() > 1.1 ) HTOFAOD[i]->Fill(CurrentRunIndex);            
               VertexAODProf[i]->Fill(CurrentRunIndex, tofaod->timeAtIpInOut());
               if(dttof->nDof()>=GlobalMinNDOFDT) VertexAODDTProf[i]->Fill(CurrentRunIndex, dttofaod->timeAtIpInOut());
               if(csctof->nDof()>=GlobalMinNDOFCSC) VertexAODCSCProf[i]->Fill(CurrentRunIndex, csctofaod->timeAtIpInOut());
            }

            if(tof && tof->nDof()>=GlobalMinNDOF && (dttof->nDof()>=GlobalMinNDOFDT || csctof->nDof()>=GlobalMinNDOFCSC) && tof->inverseBetaErr()<=GlobalMaxTOFErr){
               if(tof->inverseBeta()>=GlobalMinTOF)CountMu[i]->Fill(CurrentRunIndex);
               if(tof->inverseBeta()>=GlobalMinTOF)TOFOverMinProf[i]->Fill(CurrentRunIndex, tof->inverseBeta());
               if(dttof->inverseBeta()>=GlobalMinTOF)TOFDTOverMinProf[i]->Fill(CurrentRunIndex, dttof->inverseBeta());
               if(csctof->inverseBeta()>=GlobalMinTOF)TOFCSCOverMinProf[i]->Fill(CurrentRunIndex, csctof->inverseBeta());
               TOFProf[i]->Fill(CurrentRunIndex, tof->inverseBeta());
               if(dttof->nDof()>=GlobalMinNDOFDT) TOFDTProf[i]->Fill(CurrentRunIndex, dttof->inverseBeta());
               if(csctof->nDof()>=GlobalMinNDOFCSC) TOFCSCProf[i]->Fill(CurrentRunIndex, csctof->inverseBeta());
               if(tof->inverseBeta() > 1.1 ) HTOF[i]->Fill(CurrentRunIndex);            
               VertexProf[i]->Fill(CurrentRunIndex, tof->timeAtIpInOut());
               if(dttof->nDof()>=GlobalMinNDOFDT) VertexDTProf[i]->Fill(CurrentRunIndex, dttof->timeAtIpInOut());
               if(csctof->nDof()>=GlobalMinNDOFCSC) VertexCSCProf[i]->Fill(CurrentRunIndex, csctof->timeAtIpInOut());
            }

            if(hscp.trackRef()->pt() >60 ) HPt[i]->Fill(CurrentRunIndex);
            if(dedxSObj.dEdx() > 0.15 ) HdEdx[i]->Fill(CurrentRunIndex);
            Count[i]->Fill(CurrentRunIndex);

            dEdxProf[i]->Fill(CurrentRunIndex, dedxSObj.dEdx());
            dEdxMProf[i]->Fill(CurrentRunIndex, dedxMObj.dEdx());
            dEdxMSProf[i]->Fill(CurrentRunIndex, dedxMSObj.dEdx());
            dEdxMPProf[i]->Fill(CurrentRunIndex, dedxMPObj.dEdx());
            if(fabs(track->eta())<0.5){
            dEdxMSCProf[i]->Fill(CurrentRunIndex, dedxMSObj.dEdx());
            dEdxMPCProf[i]->Fill(CurrentRunIndex, dedxMPObj.dEdx());
            }
            if(fabs(track->eta())>1.5){
            dEdxMSFProf[i]->Fill(CurrentRunIndex, dedxMSObj.dEdx());
            dEdxMPFProf[i]->Fill(CurrentRunIndex, dedxMPObj.dEdx());
            }
            PtProf[i]->Fill(CurrentRunIndex, hscp.trackRef()->pt());
         }

      }
   }printf("\n");

   TCanvas* c1;
   TLegend* leg;

   for(unsigned int i=0;i<triggers.size();i++){
   c1 = new TCanvas("c1","c1",600,600);
   HdEdx[i]->Divide(Count[i]);
   HdEdx[i]->LabelsDeflate("X");
   HdEdx[i]->LabelsOption("av","X");
   HdEdx[i]->GetXaxis()->SetNdivisions(505);
   HdEdx[i]->SetTitle("");
   HdEdx[i]->SetStats(kFALSE);
   HdEdx[i]->GetXaxis()->SetTitle("");
   HdEdx[i]->GetYaxis()->SetTitle("Ratio over Threshold");
   HdEdx[i]->GetYaxis()->SetTitleOffset(0.9);
   HdEdx[i]->GetXaxis()->SetLabelSize(0.04);
   HdEdx[i]->SetLineColor(Color[0]);
   HdEdx[i]->SetFillColor(Color[0]);
   HdEdx[i]->SetMarkerSize(0.4);
   HdEdx[i]->SetMarkerStyle(Marker[0]);
   HdEdx[i]->SetMarkerColor(Color[0]);
   HdEdx[i]->Draw("E1");

   leg = new TLegend(0.55,0.86,0.79,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->AddEntry(HdEdx[i],"I_{as} > 0.15","P");
   leg->Draw();

   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"ROT_Is");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   HPt[i]->Divide(Count[i]);
   HPt[i]->LabelsDeflate("X");
   HPt[i]->LabelsOption("av","X");
   HPt[i]->GetXaxis()->SetNdivisions(505);
   HPt[i]->SetTitle("");
   HPt[i]->SetStats(kFALSE);
   HPt[i]->GetXaxis()->SetTitle("");
   HPt[i]->GetYaxis()->SetTitle("Ratio over Threshold");
   HPt[i]->GetYaxis()->SetTitleOffset(0.9);
   HPt[i]->GetXaxis()->SetLabelSize(0.04);
   HPt[i]->SetLineColor(Color[0]);
   HPt[i]->SetFillColor(Color[0]);
   HPt[i]->SetMarkerSize(0.4);
   HPt[i]->SetMarkerStyle(Marker[0]);
   HPt[i]->SetMarkerColor(Color[0]);
   HPt[i]->Draw("E1");

   leg = new TLegend(0.55,0.86,0.79,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->AddEntry(HPt[i],"p_{T} > 60 GeV/c","P");
   leg->Draw();
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"ROT_Pt");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   HTOFAOD[i]->Divide(CountMuAOD[i]);
   HTOFAOD[i]->LabelsDeflate("X");
   HTOFAOD[i]->LabelsOption("av","X");
   HTOFAOD[i]->GetXaxis()->SetNdivisions(505);
   HTOFAOD[i]->SetTitle("");
   HTOFAOD[i]->SetStats(kFALSE);
   HTOFAOD[i]->GetXaxis()->SetTitle("");
   HTOFAOD[i]->GetYaxis()->SetTitle("Ratio over Threshold");
   HTOFAOD[i]->GetYaxis()->SetTitleOffset(0.9);
   HTOFAOD[i]->GetXaxis()->SetLabelSize(0.04);
   HTOFAOD[i]->SetLineColor(Color[0]);
   HTOFAOD[i]->SetFillColor(Color[0]);
   HTOFAOD[i]->SetMarkerSize(0.4);
   HTOFAOD[i]->SetMarkerStyle(Marker[0]);
   HTOFAOD[i]->SetMarkerColor(Color[0]);
   HTOFAOD[i]->Draw("E1");

   leg = new TLegend(0.55,0.86,0.79,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->AddEntry(HTOF[i],"1/#beta > 1.1 (AOD)","P");
   leg->Draw();
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"ROT_TOF");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   HTOF[i]->Divide(CountMu[i]);
   HTOF[i]->LabelsDeflate("X");
   HTOF[i]->LabelsOption("av","X");
   HTOF[i]->GetXaxis()->SetNdivisions(505);
   HTOF[i]->SetTitle("");
   HTOF[i]->SetStats(kFALSE);
   HTOF[i]->GetXaxis()->SetTitle("");
   HTOF[i]->GetYaxis()->SetTitle("Ratio over Threshold");
   HTOF[i]->GetYaxis()->SetTitleOffset(0.9);
   HTOF[i]->GetXaxis()->SetLabelSize(0.04);
   HTOF[i]->SetLineColor(Color[0]);
   HTOF[i]->SetFillColor(Color[0]);
   HTOF[i]->SetMarkerSize(0.4);
   HTOF[i]->SetMarkerStyle(Marker[0]);
   HTOF[i]->SetMarkerColor(Color[0]);
   HTOF[i]->Draw("E1");

   leg = new TLegend(0.55,0.86,0.79,0.93,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->AddEntry(HTOF[i],"1/#beta > 1.1","P");
   leg->Draw();
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"ROT_TOF");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   c1->SetLogy(true);
   Count[i]->LabelsDeflate("X");
   Count[i]->LabelsOption("av","X");
   Count[i]->GetXaxis()->SetNdivisions(505);
   Count[i]->SetTitle("");
   Count[i]->SetStats(kFALSE);
   Count[i]->GetXaxis()->SetTitle("");
   Count[i]->GetYaxis()->SetTitle("#Tracks");
   Count[i]->GetYaxis()->SetTitleOffset(0.9);
   Count[i]->GetXaxis()->SetLabelSize(0.04);
   Count[i]->SetLineColor(Color[0]);
   Count[i]->SetFillColor(Color[0]);
   Count[i]->SetMarkerSize(0.4);
   Count[i]->SetMarkerStyle(Marker[0]);
   Count[i]->SetMarkerColor(Color[0]);
   Count[i]->Draw("E1");

   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Count");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   NVertProf[i]->LabelsDeflate("X");
   NVertProf[i]->LabelsOption("av","X");
   NVertProf[i]->GetXaxis()->SetNdivisions(505);
   NVertProf[i]->SetTitle("");
   NVertProf[i]->SetStats(kFALSE);
   NVertProf[i]->GetXaxis()->SetTitle("");
   NVertProf[i]->GetYaxis()->SetTitle("#RecoVertex");
   NVertProf[i]->GetYaxis()->SetTitleOffset(0.9);
   NVertProf[i]->GetXaxis()->SetLabelSize(0.04);
   NVertProf[i]->SetLineColor(Color[0]);
   NVertProf[i]->SetFillColor(Color[0]);
   NVertProf[i]->SetMarkerSize(0.4);
   NVertProf[i]->SetMarkerStyle(Marker[0]);
   NVertProf[i]->SetMarkerColor(Color[0]);
   NVertProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_NVert");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   dEdxProf[i]->LabelsDeflate("X");
   dEdxProf[i]->LabelsOption("av","X");
   dEdxProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxProf[i]->SetTitle("");
   dEdxProf[i]->SetStats(kFALSE);
   dEdxProf[i]->GetXaxis()->SetTitle("");
   dEdxProf[i]->GetYaxis()->SetTitle("dE/dx discriminator");
   dEdxProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxProf[i]->SetLineColor(Color[0]);
   dEdxProf[i]->SetFillColor(Color[0]);
   dEdxProf[i]->SetMarkerSize(0.4);
   dEdxProf[i]->SetMarkerStyle(Marker[0]);
   dEdxProf[i]->SetMarkerColor(Color[0]);
   dEdxProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_Is");
   delete c1;



   c1 = new TCanvas("c1","c1",600,600);
   dEdxMProf[i]->LabelsDeflate("X");
   dEdxMProf[i]->LabelsOption("av","X");
   dEdxMProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMProf[i]->SetTitle("");
   dEdxMProf[i]->SetStats(kFALSE);
   dEdxMProf[i]->GetXaxis()->SetTitle("");
   dEdxMProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMProf[i]->SetLineColor(Color[0]);
   dEdxMProf[i]->SetFillColor(Color[0]);
   dEdxMProf[i]->SetMarkerSize(0.4);
   dEdxMProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMProf[i]->SetMarkerColor(Color[0]);
   dEdxMProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_Im");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   dEdxMSProf[i]->LabelsDeflate("X");
   dEdxMSProf[i]->LabelsOption("av","X");
   dEdxMSProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMSProf[i]->SetTitle("");
   dEdxMSProf[i]->SetStats(kFALSE);
   dEdxMSProf[i]->GetXaxis()->SetTitle("");
   dEdxMSProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMSProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMSProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMSProf[i]->SetLineColor(Color[0]);
   dEdxMSProf[i]->SetFillColor(Color[0]);
   dEdxMSProf[i]->SetMarkerSize(0.4);
   dEdxMSProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMSProf[i]->SetMarkerColor(Color[0]);
   dEdxMSProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_ImS");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   dEdxMPProf[i]->LabelsDeflate("X");
   dEdxMPProf[i]->LabelsOption("av","X");
   dEdxMPProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMPProf[i]->SetTitle("");
   dEdxMPProf[i]->SetStats(kFALSE);
   dEdxMPProf[i]->GetXaxis()->SetTitle("");
   dEdxMPProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMPProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMPProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMPProf[i]->SetLineColor(Color[0]);
   dEdxMPProf[i]->SetFillColor(Color[0]);
   dEdxMPProf[i]->SetMarkerSize(0.4);
   dEdxMPProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMPProf[i]->SetMarkerColor(Color[0]);
   dEdxMPProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_ImP");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   dEdxMSCProf[i]->LabelsDeflate("X");
   dEdxMSCProf[i]->LabelsOption("av","X");
   dEdxMSCProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMSCProf[i]->SetTitle("");
   dEdxMSCProf[i]->SetStats(kFALSE);
   dEdxMSCProf[i]->GetXaxis()->SetTitle("");
   dEdxMSCProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMSCProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMSCProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMSCProf[i]->SetLineColor(Color[0]);
   dEdxMSCProf[i]->SetFillColor(Color[0]);
   dEdxMSCProf[i]->SetMarkerSize(0.4);
   dEdxMSCProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMSCProf[i]->SetMarkerColor(Color[0]);
   dEdxMSCProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_ImSC");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   dEdxMPCProf[i]->LabelsDeflate("X");
   dEdxMPCProf[i]->LabelsOption("av","X");
   dEdxMPCProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMPCProf[i]->SetTitle("");
   dEdxMPCProf[i]->SetStats(kFALSE);
   dEdxMPCProf[i]->GetXaxis()->SetTitle("");
   dEdxMPCProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMPCProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMPCProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMPCProf[i]->SetLineColor(Color[0]);
   dEdxMPCProf[i]->SetFillColor(Color[0]);
   dEdxMPCProf[i]->SetMarkerSize(0.4);
   dEdxMPCProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMPCProf[i]->SetMarkerColor(Color[0]);
   dEdxMPCProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_ImPC");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   dEdxMSFProf[i]->LabelsDeflate("X");
   dEdxMSFProf[i]->LabelsOption("av","X");
   dEdxMSFProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMSFProf[i]->SetTitle("");
   dEdxMSFProf[i]->SetStats(kFALSE);
   dEdxMSFProf[i]->GetXaxis()->SetTitle("");
   dEdxMSFProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMSFProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMSFProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMSFProf[i]->SetLineColor(Color[0]);
   dEdxMSFProf[i]->SetFillColor(Color[0]);
   dEdxMSFProf[i]->SetMarkerSize(0.4);
   dEdxMSFProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMSFProf[i]->SetMarkerColor(Color[0]);
   dEdxMSFProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_ImSF");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   dEdxMPFProf[i]->LabelsDeflate("X");
   dEdxMPFProf[i]->LabelsOption("av","X");
   dEdxMPFProf[i]->GetXaxis()->SetNdivisions(505);
   dEdxMPFProf[i]->SetTitle("");
   dEdxMPFProf[i]->SetStats(kFALSE);
   dEdxMPFProf[i]->GetXaxis()->SetTitle("");
   dEdxMPFProf[i]->GetYaxis()->SetTitle("dE/dx estimator");
   dEdxMPFProf[i]->GetYaxis()->SetTitleOffset(0.9);
   dEdxMPFProf[i]->GetXaxis()->SetLabelSize(0.04);
   dEdxMPFProf[i]->SetLineColor(Color[0]);
   dEdxMPFProf[i]->SetFillColor(Color[0]);
   dEdxMPFProf[i]->SetMarkerSize(0.4);
   dEdxMPFProf[i]->SetMarkerStyle(Marker[0]);
   dEdxMPFProf[i]->SetMarkerColor(Color[0]);
   dEdxMPFProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_ImPF");
   delete c1;



   c1 = new TCanvas("c1","c1",600,600);
   PtProf[i]->LabelsDeflate("X");
   PtProf[i]->LabelsOption("av","X");
   PtProf[i]->GetXaxis()->SetNdivisions(505);
   PtProf[i]->SetTitle("");
   PtProf[i]->SetStats(kFALSE);
   PtProf[i]->GetXaxis()->SetTitle("");
   PtProf[i]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   PtProf[i]->GetYaxis()->SetTitleOffset(0.9);
   PtProf[i]->GetXaxis()->SetLabelSize(0.04);
   PtProf[i]->SetLineColor(Color[0]);
   PtProf[i]->SetFillColor(Color[0]);
   PtProf[i]->SetMarkerSize(0.4);
   PtProf[i]->SetMarkerStyle(Marker[0]);
   PtProf[i]->SetMarkerColor(Color[0]);
   PtProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_Pt");
   delete c1;





   c1 = new TCanvas("c1","c1",600,600);
   TOFAODProf[i]->LabelsDeflate("X");
   TOFAODProf[i]->LabelsOption("av","X");
   TOFAODProf[i]->GetXaxis()->SetNdivisions(505);
   TOFAODProf[i]->SetTitle("");
   TOFAODProf[i]->SetStats(kFALSE);
   TOFAODProf[i]->GetXaxis()->SetTitle("");
   TOFAODProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFAODProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFAODProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFAODProf[i]->SetLineColor(Color[0]);
   TOFAODProf[i]->SetFillColor(Color[0]);
   TOFAODProf[i]->SetMarkerSize(0.4);
   TOFAODProf[i]->SetMarkerStyle(Marker[0]);
   TOFAODProf[i]->SetMarkerColor(Color[0]);
   TOFAODProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFAOD");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   TOFAODDTProf[i]->LabelsDeflate("X");
   TOFAODDTProf[i]->LabelsOption("av","X");
   TOFAODDTProf[i]->GetXaxis()->SetNdivisions(505);
   TOFAODDTProf[i]->SetTitle("");
   TOFAODDTProf[i]->SetStats(kFALSE);
   TOFAODDTProf[i]->GetXaxis()->SetTitle("");
   TOFAODDTProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFAODDTProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFAODDTProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFAODDTProf[i]->SetLineColor(Color[0]);
   TOFAODDTProf[i]->SetFillColor(Color[0]);
   TOFAODDTProf[i]->SetMarkerSize(0.4);
   TOFAODDTProf[i]->SetMarkerStyle(Marker[0]);
   TOFAODDTProf[i]->SetMarkerColor(Color[0]);
   TOFAODDTProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFAODDT");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   TOFAODCSCProf[i]->LabelsDeflate("X");
   TOFAODCSCProf[i]->LabelsOption("av","X");
   TOFAODCSCProf[i]->GetXaxis()->SetNdivisions(505);
   TOFAODCSCProf[i]->SetTitle("");
   TOFAODCSCProf[i]->SetStats(kFALSE);
   TOFAODCSCProf[i]->GetXaxis()->SetTitle("");
   TOFAODCSCProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFAODCSCProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFAODCSCProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFAODCSCProf[i]->SetLineColor(Color[0]);
   TOFAODCSCProf[i]->SetFillColor(Color[0]);
   TOFAODCSCProf[i]->SetMarkerSize(0.4);
   TOFAODCSCProf[i]->SetMarkerStyle(Marker[0]);
   TOFAODCSCProf[i]->SetMarkerColor(Color[0]);
   TOFAODCSCProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFAODCSC");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   TOFAODOverMinProf[i]->LabelsDeflate("X");
   TOFAODOverMinProf[i]->LabelsOption("av","X");
   TOFAODOverMinProf[i]->GetXaxis()->SetNdivisions(505);
   TOFAODOverMinProf[i]->SetTitle("");
   TOFAODOverMinProf[i]->SetStats(kFALSE);
   TOFAODOverMinProf[i]->GetXaxis()->SetTitle("");
   TOFAODOverMinProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFAODOverMinProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFAODOverMinProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFAODOverMinProf[i]->SetLineColor(Color[0]);
   TOFAODOverMinProf[i]->SetFillColor(Color[0]);
   TOFAODOverMinProf[i]->SetMarkerSize(0.4);
   TOFAODOverMinProf[i]->SetMarkerStyle(Marker[0]);
   TOFAODOverMinProf[i]->SetMarkerColor(Color[0]);
   TOFAODOverMinProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFAODOverMin");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   TOFAODDTOverMinProf[i]->LabelsDeflate("X");
   TOFAODDTOverMinProf[i]->LabelsOption("av","X");
   TOFAODDTOverMinProf[i]->GetXaxis()->SetNdivisions(505);
   TOFAODDTOverMinProf[i]->SetTitle("");
   TOFAODDTOverMinProf[i]->SetStats(kFALSE);
   TOFAODDTOverMinProf[i]->GetXaxis()->SetTitle("");
   TOFAODDTOverMinProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFAODDTOverMinProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFAODDTOverMinProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFAODDTOverMinProf[i]->SetLineColor(Color[0]);
   TOFAODDTOverMinProf[i]->SetFillColor(Color[0]);
   TOFAODDTOverMinProf[i]->SetMarkerSize(0.4);
   TOFAODDTOverMinProf[i]->SetMarkerStyle(Marker[0]);
   TOFAODDTOverMinProf[i]->SetMarkerColor(Color[0]);
   TOFAODDTOverMinProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFAODDTOverMin");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   TOFAODCSCOverMinProf[i]->LabelsDeflate("X");
   TOFAODCSCOverMinProf[i]->LabelsOption("av","X");
   TOFAODCSCOverMinProf[i]->GetXaxis()->SetNdivisions(505);
   TOFAODCSCOverMinProf[i]->SetTitle("");
   TOFAODCSCOverMinProf[i]->SetStats(kFALSE);
   TOFAODCSCOverMinProf[i]->GetXaxis()->SetTitle("");
   TOFAODCSCOverMinProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFAODCSCOverMinProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFAODCSCOverMinProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFAODCSCOverMinProf[i]->SetLineColor(Color[0]);
   TOFAODCSCOverMinProf[i]->SetFillColor(Color[0]);
   TOFAODCSCOverMinProf[i]->SetMarkerSize(0.4);
   TOFAODCSCOverMinProf[i]->SetMarkerStyle(Marker[0]);
   TOFAODCSCOverMinProf[i]->SetMarkerColor(Color[0]);
   TOFAODCSCOverMinProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFAODCSCOverMin");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   VertexAODProf[i]->LabelsDeflate("X");
   VertexAODProf[i]->LabelsOption("av","X");
   VertexAODProf[i]->GetXaxis()->SetNdivisions(505);
   VertexAODProf[i]->SetTitle("");
   VertexAODProf[i]->SetStats(kFALSE);
   VertexAODProf[i]->GetXaxis()->SetTitle("");
   VertexAODProf[i]->GetYaxis()->SetTitle("1/#beta");
   VertexAODProf[i]->GetYaxis()->SetTitleOffset(0.9);
   VertexAODProf[i]->GetXaxis()->SetLabelSize(0.04);
   VertexAODProf[i]->SetLineColor(Color[0]);
   VertexAODProf[i]->SetFillColor(Color[0]);
   VertexAODProf[i]->SetMarkerSize(0.4);
   VertexAODProf[i]->SetMarkerStyle(Marker[0]);
   VertexAODProf[i]->SetMarkerColor(Color[0]);
   VertexAODProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_VertexAOD");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   VertexAODDTProf[i]->LabelsDeflate("X");
   VertexAODDTProf[i]->LabelsOption("av","X");
   VertexAODDTProf[i]->GetXaxis()->SetNdivisions(505);
   VertexAODDTProf[i]->SetTitle("");
   VertexAODDTProf[i]->SetStats(kFALSE);
   VertexAODDTProf[i]->GetXaxis()->SetTitle("");
   VertexAODDTProf[i]->GetYaxis()->SetTitle("1/#beta");
   VertexAODDTProf[i]->GetYaxis()->SetTitleOffset(0.9);
   VertexAODDTProf[i]->GetXaxis()->SetLabelSize(0.04);
   VertexAODDTProf[i]->SetLineColor(Color[0]);
   VertexAODDTProf[i]->SetFillColor(Color[0]);
   VertexAODDTProf[i]->SetMarkerSize(0.4);
   VertexAODDTProf[i]->SetMarkerStyle(Marker[0]);
   VertexAODDTProf[i]->SetMarkerColor(Color[0]);
   VertexAODDTProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_VertexAODDT");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   VertexAODCSCProf[i]->LabelsDeflate("X");
   VertexAODCSCProf[i]->LabelsOption("av","X");
   VertexAODCSCProf[i]->GetXaxis()->SetNdivisions(505);
   VertexAODCSCProf[i]->SetTitle("");
   VertexAODCSCProf[i]->SetStats(kFALSE);
   VertexAODCSCProf[i]->GetXaxis()->SetTitle("");
   VertexAODCSCProf[i]->GetYaxis()->SetTitle("1/#beta");
   VertexAODCSCProf[i]->GetYaxis()->SetTitleOffset(0.9);
   VertexAODCSCProf[i]->GetXaxis()->SetLabelSize(0.04);
   VertexAODCSCProf[i]->SetLineColor(Color[0]);
   VertexAODCSCProf[i]->SetFillColor(Color[0]);
   VertexAODCSCProf[i]->SetMarkerSize(0.4);
   VertexAODCSCProf[i]->SetMarkerStyle(Marker[0]);
   VertexAODCSCProf[i]->SetMarkerColor(Color[0]);
   VertexAODCSCProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_VertexAODCSC");
   delete c1;








   c1 = new TCanvas("c1","c1",600,600);
   TOFProf[i]->LabelsDeflate("X");
   TOFProf[i]->LabelsOption("av","X");
   TOFProf[i]->GetXaxis()->SetNdivisions(505);
   TOFProf[i]->SetTitle("");
   TOFProf[i]->SetStats(kFALSE);
   TOFProf[i]->GetXaxis()->SetTitle("");
   TOFProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFProf[i]->SetLineColor(Color[0]);
   TOFProf[i]->SetFillColor(Color[0]);
   TOFProf[i]->SetMarkerSize(0.4);
   TOFProf[i]->SetMarkerStyle(Marker[0]);
   TOFProf[i]->SetMarkerColor(Color[0]);
   TOFProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOF");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   TOFDTProf[i]->LabelsDeflate("X");
   TOFDTProf[i]->LabelsOption("av","X");
   TOFDTProf[i]->GetXaxis()->SetNdivisions(505);
   TOFDTProf[i]->SetTitle("");
   TOFDTProf[i]->SetStats(kFALSE);
   TOFDTProf[i]->GetXaxis()->SetTitle("");
   TOFDTProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFDTProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFDTProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFDTProf[i]->SetLineColor(Color[0]);
   TOFDTProf[i]->SetFillColor(Color[0]);
   TOFDTProf[i]->SetMarkerSize(0.4);
   TOFDTProf[i]->SetMarkerStyle(Marker[0]);
   TOFDTProf[i]->SetMarkerColor(Color[0]);
   TOFDTProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFDT");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   TOFCSCProf[i]->LabelsDeflate("X");
   TOFCSCProf[i]->LabelsOption("av","X");
   TOFCSCProf[i]->GetXaxis()->SetNdivisions(505);
   TOFCSCProf[i]->SetTitle("");
   TOFCSCProf[i]->SetStats(kFALSE);
   TOFCSCProf[i]->GetXaxis()->SetTitle("");
   TOFCSCProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFCSCProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFCSCProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFCSCProf[i]->SetLineColor(Color[0]);
   TOFCSCProf[i]->SetFillColor(Color[0]);
   TOFCSCProf[i]->SetMarkerSize(0.4);
   TOFCSCProf[i]->SetMarkerStyle(Marker[0]);
   TOFCSCProf[i]->SetMarkerColor(Color[0]);
   TOFCSCProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFCSC");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   TOFOverMinProf[i]->LabelsDeflate("X");
   TOFOverMinProf[i]->LabelsOption("av","X");
   TOFOverMinProf[i]->GetXaxis()->SetNdivisions(505);
   TOFOverMinProf[i]->SetTitle("");
   TOFOverMinProf[i]->SetStats(kFALSE);
   TOFOverMinProf[i]->GetXaxis()->SetTitle("");
   TOFOverMinProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFOverMinProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFOverMinProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFOverMinProf[i]->SetLineColor(Color[0]);
   TOFOverMinProf[i]->SetFillColor(Color[0]);
   TOFOverMinProf[i]->SetMarkerSize(0.4);
   TOFOverMinProf[i]->SetMarkerStyle(Marker[0]);
   TOFOverMinProf[i]->SetMarkerColor(Color[0]);
   TOFOverMinProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFOverMin");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   TOFDTOverMinProf[i]->LabelsDeflate("X");
   TOFDTOverMinProf[i]->LabelsOption("av","X");
   TOFDTOverMinProf[i]->GetXaxis()->SetNdivisions(505);
   TOFDTOverMinProf[i]->SetTitle("");
   TOFDTOverMinProf[i]->SetStats(kFALSE);
   TOFDTOverMinProf[i]->GetXaxis()->SetTitle("");
   TOFDTOverMinProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFDTOverMinProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFDTOverMinProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFDTOverMinProf[i]->SetLineColor(Color[0]);
   TOFDTOverMinProf[i]->SetFillColor(Color[0]);
   TOFDTOverMinProf[i]->SetMarkerSize(0.4);
   TOFDTOverMinProf[i]->SetMarkerStyle(Marker[0]);
   TOFDTOverMinProf[i]->SetMarkerColor(Color[0]);
   TOFDTOverMinProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFDTOverMin");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   TOFCSCOverMinProf[i]->LabelsDeflate("X");
   TOFCSCOverMinProf[i]->LabelsOption("av","X");
   TOFCSCOverMinProf[i]->GetXaxis()->SetNdivisions(505);
   TOFCSCOverMinProf[i]->SetTitle("");
   TOFCSCOverMinProf[i]->SetStats(kFALSE);
   TOFCSCOverMinProf[i]->GetXaxis()->SetTitle("");
   TOFCSCOverMinProf[i]->GetYaxis()->SetTitle("1/#beta");
   TOFCSCOverMinProf[i]->GetYaxis()->SetTitleOffset(0.9);
   TOFCSCOverMinProf[i]->GetXaxis()->SetLabelSize(0.04);
   TOFCSCOverMinProf[i]->SetLineColor(Color[0]);
   TOFCSCOverMinProf[i]->SetFillColor(Color[0]);
   TOFCSCOverMinProf[i]->SetMarkerSize(0.4);
   TOFCSCOverMinProf[i]->SetMarkerStyle(Marker[0]);
   TOFCSCOverMinProf[i]->SetMarkerColor(Color[0]);
   TOFCSCOverMinProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_TOFCSCOverMin");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   VertexProf[i]->LabelsDeflate("X");
   VertexProf[i]->LabelsOption("av","X");
   VertexProf[i]->GetXaxis()->SetNdivisions(505);
   VertexProf[i]->SetTitle("");
   VertexProf[i]->SetStats(kFALSE);
   VertexProf[i]->GetXaxis()->SetTitle("");
   VertexProf[i]->GetYaxis()->SetTitle("1/#beta");
   VertexProf[i]->GetYaxis()->SetTitleOffset(0.9);
   VertexProf[i]->GetXaxis()->SetLabelSize(0.04);
   VertexProf[i]->SetLineColor(Color[0]);
   VertexProf[i]->SetFillColor(Color[0]);
   VertexProf[i]->SetMarkerSize(0.4);
   VertexProf[i]->SetMarkerStyle(Marker[0]);
   VertexProf[i]->SetMarkerColor(Color[0]);
   VertexProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_Vertex");
   delete c1;


   c1 = new TCanvas("c1","c1",600,600);
   VertexDTProf[i]->LabelsDeflate("X");
   VertexDTProf[i]->LabelsOption("av","X");
   VertexDTProf[i]->GetXaxis()->SetNdivisions(505);
   VertexDTProf[i]->SetTitle("");
   VertexDTProf[i]->SetStats(kFALSE);
   VertexDTProf[i]->GetXaxis()->SetTitle("");
   VertexDTProf[i]->GetYaxis()->SetTitle("1/#beta");
   VertexDTProf[i]->GetYaxis()->SetTitleOffset(0.9);
   VertexDTProf[i]->GetXaxis()->SetLabelSize(0.04);
   VertexDTProf[i]->SetLineColor(Color[0]);
   VertexDTProf[i]->SetFillColor(Color[0]);
   VertexDTProf[i]->SetMarkerSize(0.4);
   VertexDTProf[i]->SetMarkerStyle(Marker[0]);
   VertexDTProf[i]->SetMarkerColor(Color[0]);
   VertexDTProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_VertexDT");
   delete c1;

   c1 = new TCanvas("c1","c1",600,600);
   VertexCSCProf[i]->LabelsDeflate("X");
   VertexCSCProf[i]->LabelsOption("av","X");
   VertexCSCProf[i]->GetXaxis()->SetNdivisions(505);
   VertexCSCProf[i]->SetTitle("");
   VertexCSCProf[i]->SetStats(kFALSE);
   VertexCSCProf[i]->GetXaxis()->SetTitle("");
   VertexCSCProf[i]->GetYaxis()->SetTitle("1/#beta");
   VertexCSCProf[i]->GetYaxis()->SetTitleOffset(0.9);
   VertexCSCProf[i]->GetXaxis()->SetLabelSize(0.04);
   VertexCSCProf[i]->SetLineColor(Color[0]);
   VertexCSCProf[i]->SetFillColor(Color[0]);
   VertexCSCProf[i]->SetMarkerSize(0.4);
   VertexCSCProf[i]->SetMarkerStyle(Marker[0]);
   VertexCSCProf[i]->SetMarkerColor(Color[0]);
   VertexCSCProf[i]->Draw("E1");
   c1->Modified();
   c1->SetGridx(true);
   DrawPreliminary(triggers[i], SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,OUTDIRNAME + triggers[i],"Profile_VertexCSC");
   delete c1



;
   }


   OutputHisto->Write();
   OutputHisto->Close();  
}



bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold)
{
   unsigned int filterIndex = trEv.filterIndex(InputPath);
   //if(filterIndex<trEv.sizeFilters())printf("SELECTED INDEX =%i --> %s    XXX   %s\n",filterIndex,trEv.filterTag(filterIndex).label().c_str(), trEv.filterTag(filterIndex).process().c_str());
         
   if (filterIndex<trEv.sizeFilters()){
      const trigger::Vids& VIDS(trEv.filterIds(filterIndex));
      const trigger::Keys& KEYS(trEv.filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      const trigger::TriggerObjectCollection& TOC(trEv.getObjects());


      if(!averageThreshold){
         int NObjectAboveThresholdObserved = 0;
         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            if(TO.pt()> NewThreshold) NObjectAboveThresholdObserved++;
   	    //cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<< endl;
         }          
         if(NObjectAboveThresholdObserved>=NObjectAboveThreshold)return true;

      }else{
         std::vector<double> ObjPt;

         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            ObjPt.push_back(TO.pt());
            //cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<< endl;
         }  
         if((int)(ObjPt.size())<NObjectAboveThreshold)return false;
         std::sort(ObjPt.begin(), ObjPt.end());
         
         double Average = 0;
         for(int i=0; i<NObjectAboveThreshold;i++){
            Average+= ObjPt[ObjPt.size()-1-i];            
         }Average/=NObjectAboveThreshold;
	 //cout << "AVERAGE = " << Average << endl;
         
         if(Average>NewThreshold)return true;                  
      }
   }
   return false;
}


