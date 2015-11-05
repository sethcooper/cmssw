
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




void StabilityCheck(string DIRNAME="COMPILE", string OUTDIRNAME="pictures", string JobIndexStr="0", string NJobsStr="1")
{
  printf("DIRNAME = %s\n", DIRNAME.c_str());
  if(DIRNAME=="COMPILE") return;
  OUTDIRNAME+="/";

  int JobIndex;  sscanf(JobIndexStr.c_str(),"%d",&JobIndex);
  int NJobs;     sscanf(NJobsStr   .c_str(),"%d",&NJobs);

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
   printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n");
   printf("Run on the following samples:\n");
   for(unsigned int s=0;s<samples.size();s++){samples[s].print();}
   printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n\n");

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

   TH1D** NVert = new TH1D*[triggers.size()];;
   TH1D** dEdx = new TH1D*[triggers.size()];;
   TH1D** dEdxM = new TH1D*[triggers.size()];;
   TH1D** dEdxMS = new TH1D*[triggers.size()];;
   TH1D** dEdxMP = new TH1D*[triggers.size()];;
   TH1D** dEdxMSC = new TH1D*[triggers.size()];;
   TH1D** dEdxMPC = new TH1D*[triggers.size()];;
   TH1D** dEdxMSF = new TH1D*[triggers.size()];;
   TH1D** dEdxMPF = new TH1D*[triggers.size()];;
   TH1D** Pt   = new TH1D*[triggers.size()];;
   TH1D** TOFAOD   = new TH1D*[triggers.size()];;
   TH1D** TOFAODDT   = new TH1D*[triggers.size()];;
   TH1D** TOFAODCSC   = new TH1D*[triggers.size()];;
   TH1D** TOF   = new TH1D*[triggers.size()];;
   TH1D** TOFDT   = new TH1D*[triggers.size()];;
   TH1D** TOFCSC   = new TH1D*[triggers.size()];;
   TH1D** VertexAOD   = new TH1D*[triggers.size()];;
   TH1D** VertexAODDT   = new TH1D*[triggers.size()];;
   TH1D** VertexAODCSC   = new TH1D*[triggers.size()];;
   TH1D** Vertex   = new TH1D*[triggers.size()];;
   TH1D** VertexDT   = new TH1D*[triggers.size()];;
   TH1D** VertexCSC = new TH1D*[triggers.size()];;

 
   char OutputFileName[1024];  sprintf(OutputFileName, "%s/Histos_%i.root", OUTDIRNAME.c_str(), JobIndex);
   TFile* OutputHisto = new TFile(OutputFileName,"RECREATE");
   TypeMode      = 0;

   std::unordered_map<unsigned int,double> TrackerGains;
   double dEdxSF [2];
   dEdxSF [0] = 1.00000;
   dEdxSF [1] = 1.21836;
   TH3F* dEdxTemplates = loadDeDxTemplate("../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true);

//   LoadDeDxCalibration(TrackerGains, "../../../data/Data13TeVGains.root"); 
 
   moduleGeom::loadGeometry("../../../data/CMS_GeomTree.root");
   muonTimingCalculator tofCalculator;
   tofCalculator.loadTimeOffset("../../../data/MuonTimeOffset.txt");
   unsigned int CurrentRun = 0;

   FILE * gainsTXT = fopen ("../DeDxStudy/gains.txt", "r");
   char GainsFile [19];
   vector <string> GainsFiles;
   for (int i = 0; i < 33; i++){
       fscanf (gainsTXT, "%s %*s %*s", GainsFile);
       GainsFiles.push_back(string(GainsFile));
   }

   LoadDeDxCalibration(TrackerGains, "../DeDxStudy/Gains/"+GainsFiles[0]+".root"); 

   fwlite::ChainEvent ev(DataFileName);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Looping on Tree              :");

   int NEvents = ev.size() / NJobs;
   int FirstEvent = JobIndex * NEvents;
   int TreeStep = NEvents/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t e=FirstEvent;e<FirstEvent+NEvents;e++){
      ev.to(e); 
      if(e%TreeStep==0){printf(".");fflush(stdout);}

      //if run change, update conditions
      if(CurrentRun != ev.eventAuxiliary().run()){
         CurrentRun = ev.eventAuxiliary().run();
         tofCalculator.setRun(CurrentRun);

         reloadGainsFile (TrackerGains, "../DeDxStudy/Gains/", GainsFiles, CurrentRun);

         TDirectory* dir = OutputHisto;
         char DIRECTORY[2048]; sprintf(DIRECTORY,"%6i",ev.eventAuxiliary().run());
         TDirectory::AddDirectory(kTRUE);
         TH1::AddDirectory(kTRUE);
         dir = OutputHisto->mkdir(DIRECTORY, DIRECTORY);
         dir->cd();

         for(unsigned int i=0;i<triggers.size();i++){
            NVert[i] = new TH1D((triggers[i] + "NVert"    ).c_str(), "NVert"  , 100, 0.0, 100);
            Pt  [i]  = new TH1D((triggers[i] + "Pt"       ).c_str(), "Pt"     ,1000, 0.0,1000);

            dEdx[i]    = new TH1D((triggers[i] + "dEdx"   ).c_str(), "dEdx"   , 100, 0.0, 5.0);
            dEdxM[i]   = new TH1D((triggers[i] + "dEdxM"  ).c_str(), "dEdxM"  , 100, 0.0, 5.0);
            dEdxMS[i]  = new TH1D((triggers[i] + "dEdxMS" ).c_str(), "dEdxMS" , 100, 0.0, 5.0);
            dEdxMP[i]  = new TH1D((triggers[i] + "dEdxMP" ).c_str(), "dEdxMP" , 100, 0.0, 5.0);
            dEdxMSC[i] = new TH1D((triggers[i] + "dEdxMSC").c_str(), "dEdxMSC", 100, 0.0, 5.0);
            dEdxMPC[i] = new TH1D((triggers[i] + "dEdxMPC").c_str(), "dEdxMPC", 100, 0.0, 5.0);
            dEdxMSF[i] = new TH1D((triggers[i] + "dEdxMSF").c_str(), "dEdxMSF", 100, 0.0, 5.0);
            dEdxMPF[i] = new TH1D((triggers[i] + "dEdxMPF").c_str(), "dEdxMPF", 100, 0.0, 5.0);

            TOFAOD   [i] = new TH1D((triggers[i] + "TOFAOD"  ).c_str(), "TOFAOD"      , 100, -1.0, 3.0);
            TOFAODDT [i] = new TH1D((triggers[i] + "TOFAODDT"  ).c_str(), "TOFAODDT"  , 100, -1.0, 3.0);
            TOFAODCSC[i] = new TH1D((triggers[i] + "TOFAODCSC"  ).c_str(), "TOFAODCSC", 100, -1.0, 3.0);

            TOF      [i] = new TH1D((triggers[i] + "TOF"  ).c_str(), "TOF"            , 100, -1.0, 3.0);
            TOFDT    [i] = new TH1D((triggers[i] + "TOFDT"  ).c_str(), "TOFDT"        , 100, -1.0, 3.0);
            TOFCSC   [i] = new TH1D((triggers[i] + "TOFCSC"  ).c_str(), "TOFCSC"      , 100, -1.0, 3.0);

            VertexAOD   [i] = new TH1D((triggers[i] + "VertexAOD"  ).c_str(), "VertexAOD"      , 100, -10.0, 10.0);
            VertexAODDT [i] = new TH1D((triggers[i] + "VertexAODDT"  ).c_str(), "VertexAODDT"  , 100, -10.0, 10.0);
            VertexAODCSC[i] = new TH1D((triggers[i] + "VertexAODCSC"  ).c_str(), "VertexAODCSC", 100, -10.0, 10.0);

            Vertex      [i] = new TH1D((triggers[i] + "Vertex"  ).c_str(), "Vertex"            , 100, -10.0, 10.0);
            VertexDT    [i] = new TH1D((triggers[i] + "VertexDT"  ).c_str(), "VertexDT"        , 100, -10.0, 10.0);
            VertexCSC   [i] = new TH1D((triggers[i] + "VertexCSC"  ).c_str(), "VertexCSC"      , 100, -10.0, 10.0);
         }
      }

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

//         printf("AOD %i %i %iOR%i %i\n", tofaod!=NULL?1:0, tofaod->nDof()>=GlobalMinNDOF?1:0,  dttofaod->nDof()>=GlobalMinNDOFDT?1:0 , csctofaod->nDof()>=GlobalMinNDOFCSC?1:0, tofaod->inverseBetaErr()<=GlobalMaxTOFErr?1:0  );
//         printf("OTF %i %i %iOR%i %i\n", tof   !=NULL?1:0, tof   ->nDof()>=GlobalMinNDOF?1:0,  dttof   ->nDof()>=GlobalMinNDOFDT?1:0 , csctof   ->nDof()>=GlobalMinNDOFCSC?1:0, tof   ->inverseBetaErr()<=GlobalMaxTOFErr?1:0  );

         for(unsigned int i=0;i<triggers.size();i++){
            if(!PassingTrigger(ev,triggers[i])){continue;}

            NVert[i]->Fill(vertexColl.size()); 
             

            if(tofaod && tofaod->nDof()>=GlobalMinNDOF && (dttofaod->nDof()>=GlobalMinNDOFDT || csctofaod->nDof()>=GlobalMinNDOFCSC) && tofaod->inverseBetaErr()<=GlobalMaxTOFErr && fabs(dttofaod->inverseBeta()-1)<50){
               TOFAOD[i]->Fill(tofaod->inverseBeta());
               if(dttofaod->nDof()>=GlobalMinNDOFDT) TOFAODDT[i]->Fill(dttofaod->inverseBeta());
               if(csctofaod->nDof()>=GlobalMinNDOFCSC) TOFAODCSC[i]->Fill(csctofaod->inverseBeta());
               VertexAOD[i]->Fill(tofaod->timeAtIpInOut());
               if(dttofaod->nDof()>=GlobalMinNDOFDT) VertexAODDT[i]->Fill(dttofaod->timeAtIpInOut());
               if(csctofaod->nDof()>=GlobalMinNDOFCSC) VertexAODCSC[i]->Fill(csctofaod->timeAtIpInOut());
            }

            if(tof && tof->nDof()>=GlobalMinNDOF && (dttof->nDof()>=GlobalMinNDOFDT || csctof->nDof()>=GlobalMinNDOFCSC) && tof->inverseBetaErr()<=GlobalMaxTOFErr && fabs(dttof->inverseBeta()-1)<50){
               TOF[i]->Fill(tof->inverseBeta());
               if(dttof->nDof()>=GlobalMinNDOFDT) TOFDT[i]->Fill(dttof->inverseBeta());
               if(csctof->nDof()>=GlobalMinNDOFCSC) TOFCSC[i]->Fill(csctof->inverseBeta());
               Vertex[i]->Fill(tof->timeAtIpInOut());
               if(dttof->nDof()>=GlobalMinNDOFDT) VertexDT[i]->Fill(dttof->timeAtIpInOut());
               if(csctof->nDof()>=GlobalMinNDOFCSC) VertexCSC[i]->Fill(csctof->timeAtIpInOut());
            }  

/*            if(i==0 && tof && tof->nDof()>=GlobalMinNDOF && (dttof->nDof()>=GlobalMinNDOFDT || csctof->nDof()>=GlobalMinNDOFCSC) && tof->inverseBetaErr()<=GlobalMaxTOFErr && fabs(dttof->inverseBeta()-1)>50){
               printf("Large values %f+-%f vs %f+-%f (aod)\n", dttof->inverseBeta(), dttof->inverseBetaErr(), dttofaod->inverseBeta(), dttofaod->inverseBetaErr());
               printf("%i vs %i  (Min=%f)\n", dttof->nDof(), dttofaod->nDof(), GlobalMinNDOFDT);
               const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
               const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
               tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, 1, true ); //apply T0 correction on data but not on signal MC
            }
*/

            dEdx[i]->Fill(dedxSObj.dEdx());
            dEdxM[i]->Fill(dedxMObj.dEdx());
            dEdxMS[i]->Fill(dedxMSObj.dEdx());
            dEdxMP[i]->Fill(dedxMPObj.dEdx());
            if(fabs(track->eta())<0.5){
            dEdxMSC[i]->Fill(dedxMSObj.dEdx());
            dEdxMPC[i]->Fill(dedxMPObj.dEdx());
            }
            if(fabs(track->eta())>1.5){
            dEdxMSF[i]->Fill(dedxMSObj.dEdx());
            dEdxMPF[i]->Fill(dedxMPObj.dEdx());
            }
            Pt[i]->Fill(hscp.trackRef()->pt());
         }
      }
   }printf("\n");

/*
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

   TOFAODDT2D[i]->LabelsDeflate("X");
   TOFAODDT2D[i]->LabelsOption("av","X");
   TOFDT2D[i]->LabelsDeflate("X");
   TOFDT2D[i]->LabelsOption("av","X");


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
*/

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


