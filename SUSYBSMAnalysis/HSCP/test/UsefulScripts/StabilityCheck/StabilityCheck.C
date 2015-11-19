
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

   int highestPtGoodVertex = -1;
   int goodVerts=0;
   double dzMin=10000;
   for(unsigned int i=0;i<vertexColl.size();i++){
     if(vertexColl[i].isFake() || fabs(vertexColl[i].z())>24 || vertexColl[i].position().rho()>2 || vertexColl[i].ndof()<=4)continue; //only consider good vertex
     goodVerts++;

//     if(highestPtGoodVertex<0)highestPtGoodVertex = i;
     if(fabs(track->dz (vertexColl[i].position())) < fabs(dzMin) ){
         dzMin = fabs(track->dz (vertexColl[i].position()));
         highestPtGoodVertex = i;
//       dz  = track->dz (vertexColl[i].position());
//       dxy = track->dxy(vertexColl[i].position());
     }
   }if(highestPtGoodVertex<0)highestPtGoodVertex=0;

   double dz  = track->dz (vertexColl[highestPtGoodVertex].position());
   double dxy = track->dxy(vertexColl[highestPtGoodVertex].position());

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
   int sampleIdStart, sampleIdEnd; sscanf(JobIndexStr.c_str(),"%d",&sampleIdStart); sampleIdEnd=sampleIdStart;
   keepOnlyTheXtoYSamples(samples,sampleIdStart,sampleIdEnd);
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
   TH1D** dEdxHitStripAOD= new TH1D*[triggers.size()];;
   TH1D** dEdxHitPixelAOD= new TH1D*[triggers.size()];;
   TH1D** dEdxMin1AOD= new TH1D*[triggers.size()];;
   TH1D** dEdxMin2AOD= new TH1D*[triggers.size()];;
   TH1D** dEdxMin3AOD= new TH1D*[triggers.size()];;
   TH1D** dEdxMin4AOD= new TH1D*[triggers.size()];;
   TH1D** dEdxAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMTAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMSAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMPAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMSCAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMPCAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMSFAOD = new TH1D*[triggers.size()];;
   TH1D** dEdxMPFAOD = new TH1D*[triggers.size()];;

   TH1D** dEdxHitStrip= new TH1D*[triggers.size()];;
   TH1D** dEdxHitPixel= new TH1D*[triggers.size()];;
   TH1D** dEdxMin1= new TH1D*[triggers.size()];;
   TH1D** dEdxMin2= new TH1D*[triggers.size()];;
   TH1D** dEdxMin3= new TH1D*[triggers.size()];;
   TH1D** dEdxMin4= new TH1D*[triggers.size()];;
   TH1D** dEdx = new TH1D*[triggers.size()];;
   TH1D** dEdxMT = new TH1D*[triggers.size()];;
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

   dedxGainCorrector trackerCorrector;
   double dEdxSF [2];
   dEdxSF [0] = 1.00000;
   dEdxSF [1] = 1.21836;
   TH3F* dEdxTemplates = loadDeDxTemplate("../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true);

   if(true){    trackerCorrector.LoadDeDxCalibration("../../../data/Data13TeVGains_v2.root");
   }else{ trackerCorrector.TrackerGains = NULL; //FIXME check gain for MC
   }
 
   moduleGeom::loadGeometry("../../../data/CMS_GeomTree.root");
   muonTimingCalculator tofCalculator;
   tofCalculator.loadTimeOffset("../../../data/MuonTimeOffset.txt");
   unsigned int CurrentRun = 0;

   fwlite::ChainEvent ev(DataFileName);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Looping on Tree              :");

   int NEvents = ev.size();// / NJobs;
   int FirstEvent = 0;//JobIndex * NEvents;
   int TreeStep = NEvents/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t e=FirstEvent;e<FirstEvent+NEvents;e++){
      ev.to(e); 
      if(e%TreeStep==0){printf(".");fflush(stdout);}

      //if run change, update conditions
      if(CurrentRun != ev.eventAuxiliary().run()){
         CurrentRun = ev.eventAuxiliary().run();
         tofCalculator.setRun(CurrentRun);
         trackerCorrector.setRun(CurrentRun);

         TDirectory* dir = OutputHisto;
         char DIRECTORY[2048]; sprintf(DIRECTORY,"%6i",ev.eventAuxiliary().run());
         TDirectory::AddDirectory(kTRUE);
         TH1::AddDirectory(kTRUE);
         dir = (TDirectory*)OutputHisto->Get(DIRECTORY);
         if(dir==NULL){
            dir = OutputHisto->mkdir(DIRECTORY, DIRECTORY);
            dir->cd();

            for(unsigned int i=0;i<triggers.size();i++){
               NVert[i] = new TH1D((triggers[i] + "NVert"    ).c_str(), "NVert"  , 100, 0.0, 100);
               Pt  [i]  = new TH1D((triggers[i] + "Pt"       ).c_str(), "Pt"     ,1000, 0.0,1000);

               dEdxHitStripAOD[i] = new TH1D((triggers[i] + "dEdxHitStripAOD").c_str(), "dEdxHitStripAOD", 400, 0.0,20.0);
               dEdxHitPixelAOD[i] = new TH1D((triggers[i] + "dEdxHitPixelAOD").c_str(), "dEdxHitPixelAOD", 400, 0.0,20.0);
               dEdxMin1AOD[i] = new TH1D((triggers[i] + "dEdxMin1AOD").c_str(), "dEdxMin1", 200, 0.0,10.0);
               dEdxMin2AOD[i] = new TH1D((triggers[i] + "dEdxMin2AOD").c_str(), "dEdxMin2", 200, 0.0,10.0);
               dEdxMin3AOD[i] = new TH1D((triggers[i] + "dEdxMin3AOD").c_str(), "dEdxMin3", 200, 0.0,10.0);
               dEdxMin4AOD[i] = new TH1D((triggers[i] + "dEdxMin4AOD").c_str(), "dEdxMin4", 200, 0.0,10.0);

               dEdxAOD[i]    = new TH1D((triggers[i] + "dEdxAOD"   ).c_str(), "dEdxAOD"   , 100, 0.0, 1.0);
               dEdxMTAOD[i]  = new TH1D((triggers[i] + "dEdxMTAOD" ).c_str(), "dEdxMTAOD" , 200, 0.0,10.0);
               dEdxMAOD[i]   = new TH1D((triggers[i] + "dEdxMAOD"  ).c_str(), "dEdxMAOD"  , 200, 0.0,10.0);
               dEdxMSAOD[i]  = new TH1D((triggers[i] + "dEdxMSAOD" ).c_str(), "dEdxMSAOD" , 200, 0.0,10.0);
               dEdxMPAOD[i]  = new TH1D((triggers[i] + "dEdxMPAOD" ).c_str(), "dEdxMPAOD" , 200, 0.0,10.0);
               dEdxMSCAOD[i] = new TH1D((triggers[i] + "dEdxMSCAOD").c_str(), "dEdxMSCAOD", 200, 0.0,10.0);
               dEdxMPCAOD[i] = new TH1D((triggers[i] + "dEdxMPCAOD").c_str(), "dEdxMPCAOD", 200, 0.0,10.0);
               dEdxMSFAOD[i] = new TH1D((triggers[i] + "dEdxMSFAOD").c_str(), "dEdxMSFAOD", 200, 0.0,10.0);
               dEdxMPFAOD[i] = new TH1D((triggers[i] + "dEdxMPFAOD").c_str(), "dEdxMPFAOD", 200, 0.0,10.0);


               dEdxHitStrip[i] = new TH1D((triggers[i] + "dEdxHitStrip").c_str(), "dEdxHitStrip", 400, 0.0,20.0);
               dEdxHitPixel[i] = new TH1D((triggers[i] + "dEdxHitPixel").c_str(), "dEdxHitPixel", 400, 0.0,20.0);
               dEdxMin1[i] = new TH1D((triggers[i] + "dEdxMin1").c_str(), "dEdxMin1", 200, 0.0,10.0);
               dEdxMin2[i] = new TH1D((triggers[i] + "dEdxMin2").c_str(), "dEdxMin2", 200, 0.0,10.0);
               dEdxMin3[i] = new TH1D((triggers[i] + "dEdxMin3").c_str(), "dEdxMin3", 200, 0.0,10.0);
               dEdxMin4[i] = new TH1D((triggers[i] + "dEdxMin4").c_str(), "dEdxMin4", 200, 0.0,10.0);

               dEdx[i]    = new TH1D((triggers[i] + "dEdx"   ).c_str(), "dEdx"   , 100, 0.0, 1.0);
               dEdxMT[i]  = new TH1D((triggers[i] + "dEdxMT" ).c_str(), "dEdxMT" , 200, 0.0,10.0);
               dEdxM[i]   = new TH1D((triggers[i] + "dEdxM"  ).c_str(), "dEdxM"  , 200, 0.0,10.0);
               dEdxMS[i]  = new TH1D((triggers[i] + "dEdxMS" ).c_str(), "dEdxMS" , 200, 0.0,10.0);
               dEdxMP[i]  = new TH1D((triggers[i] + "dEdxMP" ).c_str(), "dEdxMP" , 200, 0.0,10.0);
               dEdxMSC[i] = new TH1D((triggers[i] + "dEdxMSC").c_str(), "dEdxMSC", 200, 0.0,10.0);
               dEdxMPC[i] = new TH1D((triggers[i] + "dEdxMPC").c_str(), "dEdxMPC", 200, 0.0,10.0);
               dEdxMSF[i] = new TH1D((triggers[i] + "dEdxMSF").c_str(), "dEdxMSF", 200, 0.0,10.0);
               dEdxMPF[i] = new TH1D((triggers[i] + "dEdxMPF").c_str(), "dEdxMPF", 200, 0.0,10.0);

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
         }else{
            for(unsigned int i=0;i<triggers.size();i++){
               NVert[i] = (TH1D*)dir->Get((triggers[i] + "NVert").c_str());
               Pt  [i]  = (TH1D*)dir->Get((triggers[i] + "Pt").c_str());

               dEdxHitStripAOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxHitStripAOD").c_str());
               dEdxHitPixelAOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxHitPixelAOD").c_str());
               dEdxMin1AOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin1AOD").c_str());
               dEdxMin2AOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin2AOD").c_str());
               dEdxMin3AOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin3AOD").c_str());
               dEdxMin4AOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin4AOD").c_str());

               dEdxAOD[i]    = (TH1D*)dir->Get((triggers[i] + "dEdxAOD").c_str());
               dEdxMTAOD[i]   = (TH1D*)dir->Get((triggers[i] + "dEdxMTAOD").c_str());
               dEdxMAOD[i]   = (TH1D*)dir->Get((triggers[i] + "dEdxMAOD").c_str());
               dEdxMSAOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMSAOD").c_str());
               dEdxMPAOD[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMPAOD").c_str());
               dEdxMSCAOD[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMSCAOD").c_str());
               dEdxMPCAOD[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMPCAOD").c_str());
               dEdxMSFAOD[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMSFAOD").c_str());
               dEdxMPFAOD[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMPFAOD").c_str());


               dEdxHitStrip[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxHitStrip").c_str());
               dEdxHitPixel[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxHitPixel").c_str());
               dEdxMin1[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin1").c_str());
               dEdxMin2[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin2").c_str());
               dEdxMin3[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin3").c_str());
               dEdxMin4[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMin4").c_str());

               dEdx[i]    = (TH1D*)dir->Get((triggers[i] + "dEdx").c_str());
               dEdxMT[i]   = (TH1D*)dir->Get((triggers[i] + "dEdxMT").c_str());
               dEdxM[i]   = (TH1D*)dir->Get((triggers[i] + "dEdxM").c_str());
               dEdxMS[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMS").c_str());
               dEdxMP[i]  = (TH1D*)dir->Get((triggers[i] + "dEdxMP").c_str());
               dEdxMSC[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMSC").c_str());
               dEdxMPC[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMPC").c_str());
               dEdxMSF[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMSF").c_str());
               dEdxMPF[i] = (TH1D*)dir->Get((triggers[i] + "dEdxMPF").c_str());

               TOFAOD   [i] = (TH1D*)dir->Get((triggers[i] + "TOFAOD").c_str());
               TOFAODDT [i] = (TH1D*)dir->Get((triggers[i] + "TOFAODDT").c_str());
               TOFAODCSC[i] = (TH1D*)dir->Get((triggers[i] + "TOFAODCSC").c_str());

               TOF      [i] = (TH1D*)dir->Get((triggers[i] + "TOF").c_str());
               TOFDT    [i] = (TH1D*)dir->Get((triggers[i] + "TOFDT").c_str());
               TOFCSC   [i] = (TH1D*)dir->Get((triggers[i] + "TOFCSC").c_str());

               VertexAOD   [i] = (TH1D*)dir->Get((triggers[i] + "VertexAOD").c_str());
               VertexAODDT [i] = (TH1D*)dir->Get((triggers[i] + "VertexAODDT").c_str());
               VertexAODCSC[i] = (TH1D*)dir->Get((triggers[i] + "VertexAODCSC").c_str());

               Vertex      [i] = (TH1D*)dir->Get((triggers[i] + "Vertex").c_str());
               VertexDT    [i] = (TH1D*)dir->Get((triggers[i] + "VertexDT").c_str());
               VertexCSC   [i] = (TH1D*)dir->Get((triggers[i] + "VertexCSC").c_str());
            }
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
         DeDxData dedxMin1AODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, NULL, true, true, 99, false, 1, 0.1);
         DeDxData dedxMin2AODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, NULL, true, true, 99, false, 1, 0.2);
         DeDxData dedxMin3AODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, NULL, true, true, 99, false, 1, 0.3);
         DeDxData dedxMin4AODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, NULL, true, true, 99, false, 1, 0.4);
         DeDxData dedxSAODObj = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode==5, false, NULL, true, true, 99, false, 1);
         DeDxData dedxMAODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, NULL, true, true, 99, false, 1);
         DeDxData dedxMTAODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , true, NULL, true, true, 99, false, 1);
         DeDxData dedxMSAODObj = computedEdx(dedxHits, dEdxSF, NULL,          false,useClusterCleaning, false      , false, NULL, true, true, 99, false, 1, 0.0, dEdxHitStripAOD[2]);
         DeDxData dedxMPAODObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, NULL, false, true, 99, false, 1, 0.0, dEdxHitPixelAOD[2]);


         DeDxData dedxMin1Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.1);
         DeDxData dedxMin2Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.2);
         DeDxData dedxMin3Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.3);
         DeDxData dedxMin4Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.4);
         DeDxData dedxSObj = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode==5, false, trackerCorrector.TrackerGains, true, true, 99, false, 1);
         DeDxData dedxMObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1);
         DeDxData dedxMTObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , true, trackerCorrector.TrackerGains, true, true, 99, false, 1);
         DeDxData dedxMSObj = computedEdx(dedxHits, dEdxSF, NULL,          false,useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.0, dEdxHitStrip[2]);
         DeDxData dedxMPObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, false, true, 99, false, 1, 0.0, dEdxHitPixel[2]);

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


            dEdxMin1AOD[i]->Fill(dedxMin1AODObj.dEdx());
            dEdxMin2AOD[i]->Fill(dedxMin2AODObj.dEdx());
            dEdxMin3AOD[i]->Fill(dedxMin3AODObj.dEdx());
            dEdxMin4AOD[i]->Fill(dedxMin4AODObj.dEdx());
            dEdxAOD[i]->Fill(dedxSAODObj.dEdx());
            dEdxMTAOD[i]->Fill(dedxMTAODObj.dEdx());
            dEdxMAOD[i]->Fill(dedxMAODObj.dEdx());
            dEdxMSAOD[i]->Fill(dedxMSAODObj.dEdx());
            dEdxMPAOD[i]->Fill(dedxMPAODObj.dEdx());
            if(fabs(track->eta())<0.5){
            dEdxMSCAOD[i]->Fill(dedxMSAODObj.dEdx());
            dEdxMPCAOD[i]->Fill(dedxMPAODObj.dEdx());
            }
            if(fabs(track->eta())>1.5){
            dEdxMSFAOD[i]->Fill(dedxMSAODObj.dEdx());
            dEdxMPFAOD[i]->Fill(dedxMPAODObj.dEdx());
            }

            dEdxMin1[i]->Fill(dedxMin1Obj.dEdx());
            dEdxMin2[i]->Fill(dedxMin2Obj.dEdx());
            dEdxMin3[i]->Fill(dedxMin3Obj.dEdx());
            dEdxMin4[i]->Fill(dedxMin4Obj.dEdx());
            dEdx[i]->Fill(dedxSObj.dEdx());
            dEdxMT[i]->Fill(dedxMTObj.dEdx());
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


