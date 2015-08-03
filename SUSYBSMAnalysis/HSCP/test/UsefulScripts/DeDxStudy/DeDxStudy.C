
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
#include "TCutG.h"
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
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;


#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"


#endif

const double P_Min               = 1;
const double P_Max               = 15;
const int    P_NBins             = 14  ;
const double Path_Min            = 0.2 ;
const double Path_Max            = 1.6 ;
const int    Path_NBins          = 42  ;
const double Charge_Min          = 0   ;
const double Charge_Max          = 5000;
const int    Charge_NBins        = 500 ;

struct dEdxStudyObj
{
   string Name;
   bool isDiscrim;
   bool isEstim;
   bool isHit;

   bool usePixel;
   bool useStrip;

   TH3D* Charge_Vs_Path;
   TH1D* HdedxMIP;
   TH1D* HdedxSIG;
   TH2D* HdedxVsPHSCP;
   TH2D* HdedxVsP;
   TH2D* HdedxVsQP;
   TProfile2D* HdedxVsP_NS;
   TProfile* HdedxVsPProfile;
   TProfile* HdedxVsEtaProfile;
   TH2D* HdedxVsEta;
   TProfile* HNOSVsEtaProfile;
   TProfile* HNOMVsEtaProfile;
   TProfile* HNOMSVsEtaProfile;
   TH1D* HMass;
   TH1D* HP;
   TH1D* HHit; 
   TProfile* Charge_Vs_FS[15];

   TH3F* dEdxTemplates = NULL;
   std::unordered_map<unsigned int,double>* TrackerGains = NULL;

   dEdxStudyObj(string Name_, int type_, int subdet_, TH3F* dEdxTemplates_=NULL, std::unordered_map<unsigned int,double>* TrackerGains_=NULL){
      Name = Name_;

      if     (type_==0){ isHit=true;  isEstim= false; isDiscrim = false;}
      else if(type_==1){ isHit=false; isEstim= true;  isDiscrim = false;}
      else if(type_==2){ isHit=false; isEstim= false; isDiscrim = true; }
      else             { isHit=false; isEstim= false; isDiscrim = false;}

           if(subdet_==1){ usePixel = true;  useStrip = false;}
      else if(subdet_==2){ usePixel = false; useStrip = true; }
      else               { usePixel = true;  useStrip = true; }

      dEdxTemplates = dEdxTemplates_;
      TrackerGains  = TrackerGains_;

      string HistoName;
      //HitLevel plot
      if(isHit){ 
         HistoName = Name + "_ChargeVsPath";      Charge_Vs_Path        = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_Hit";               HHit                  = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20); 
         for(unsigned int g=0;g<15;g++){
            char Id[255]; sprintf(Id, "%02i", g);
            HistoName = Name + "_ChargeVsFS"+Id;    Charge_Vs_FS[g]       = new TProfile( HistoName.c_str(), HistoName.c_str(),  800, 0, 800);
         }
      }

      //Track Level plots
      if(isEstim || isDiscrim){
         HistoName = Name + "_MIP";               HdedxMIP              = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, isDiscrim?1.0:20);
         HistoName = Name + "_SIG";               HdedxSIG              = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, isDiscrim?1.0:20);
         HistoName = Name + "_dedxVsPHSCP";       HdedxVsPHSCP          = new TH2D(      HistoName.c_str(), HistoName.c_str(), 3000, 0, 2000,1500,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP";           HdedxVsP              = new TH2D(      HistoName.c_str(), HistoName.c_str(), 3000, 0, 30,1500,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsQP";          HdedxVsQP             = new TH2D(      HistoName.c_str(), HistoName.c_str(), 6000, -30, 30,1500,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_NS";        HdedxVsP_NS           = new TProfile2D(HistoName.c_str(), HistoName.c_str(), 3000, 0, 30,1500,0, isDiscrim?1.0:15);
         HistoName = Name + "_Profile";           HdedxVsPProfile       = new TProfile(  HistoName.c_str(), HistoName.c_str(),  100, 0,100);
         HistoName = Name + "_Eta";               HdedxVsEtaProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),  100,-3,  3);
         HistoName = Name + "_Eta2D";             HdedxVsEta            = new TH2D(      HistoName.c_str(), HistoName.c_str(),  100,-3,  3, 1000,0, isDiscrim?1.0:5);
         HistoName = Name + "_NOS";               HNOSVsEtaProfile      = new TProfile(  HistoName.c_str(), HistoName.c_str(),  100,-3,  3);
         HistoName = Name + "_NOM";               HNOMVsEtaProfile      = new TProfile(  HistoName.c_str(), HistoName.c_str(),  100,-3,  3);
         HistoName = Name + "_NOMS";              HNOMSVsEtaProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),  100,-3,  3);
         HistoName = Name + "_P";                 HP                    = new TH1D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 100);  
      }

      //estimator plot only
      if(isEstim){
         HistoName = Name + "_Mass";              HMass                 = new TH1D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10);
      }

   }
};

void DeDxStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root")
{
  if(DIRNAME=="COMPILE") return;

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

   bool isData = INPUT.find("MC")!=string::npos ? false : true;
   std::vector<string> FileName;
   if(INPUT.find(".root")<std::string::npos){
      char* pch=strtok(&INPUT[0],",");
      while (pch!=NULL){
         FileName.push_back(pch);    
         pch=strtok(NULL,",");
      }
   }else{
      string SampleId = INPUT;
      InitBaseDirectory();
      GetSampleDefinition(samples , DIRNAME+"/../../AnalysisCode/Analysis_Samples.txt");
      stSample& sample = samples[JobIdToIndex(SampleId, samples)];
      isData = (sample.Type==0);
      GetInputFiles(sample, BaseDirectory, FileName, 0);
   }

   TH3F* dEdxTemplates    = NULL;
   TH3F* dEdxTemplatesInc = NULL;
   if(isData){   //FIXME update template on data directory
	 dEdxSF           = 1.0;
         dEdxTemplates    = loadDeDxTemplate(DIRNAME + "/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", true);
         dEdxTemplatesInc = loadDeDxTemplate(DIRNAME + "/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", false);
   }else{
	 dEdxSF           = 1.09708;
         dEdxTemplates    = loadDeDxTemplate(DIRNAME + "/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", true);
         dEdxTemplatesInc = loadDeDxTemplate(DIRNAME + "/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", false); 
   }

   std::unordered_map<unsigned int,double> TrackerGains;
   LoadDeDxCalibration(TrackerGains, DIRNAME+"/../../../data/Data13TeVGains.root");

   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");  //File must be opened before the histogram are created

   std::vector<dEdxStudyObj*> results;
   results.push_back(new dEdxStudyObj("hit_PO"      , 0, 1, NULL            , NULL) );
   results.push_back(new dEdxStudyObj("hit_SO_raw"  , 0, 2, NULL            , NULL) );
   results.push_back(new dEdxStudyObj("hit_SO"      , 0, 2, NULL            , &TrackerGains) );
   results.push_back(new dEdxStudyObj("hit_SP"      , 0, 3, NULL            , &TrackerGains) );
   results.push_back(new dEdxStudyObj("harm2_PO_raw", 1, 1, NULL            , NULL) );
   results.push_back(new dEdxStudyObj("harm2_SO_raw", 1, 2, NULL            , NULL) );
   results.push_back(new dEdxStudyObj("harm2_SP_raw", 1, 3, NULL            , NULL) );
   results.push_back(new dEdxStudyObj("harm2_SO"    , 1, 2, NULL            , &TrackerGains) );
   results.push_back(new dEdxStudyObj("harm2_SP"    , 1, 3, NULL            , &TrackerGains) );
   results.push_back(new dEdxStudyObj("Ias_SO_inc"  , 2, 2, dEdxTemplatesInc, NULL) );
   results.push_back(new dEdxStudyObj("Ias_SO"      , 2, 2, dEdxTemplates   , NULL) );

   fwlite::ChainEvent ev(FileName);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Looping on Tree              :");
   int TreeStep = ev.size()/50;if(TreeStep==0)TreeStep=1;
   for(Long64_t e=0;e<ev.size();e++){
      ev.to(e); 
      if(e%TreeStep==0){printf(".");fflush(stdout);}

      fwlite::Handle<DeDxHitInfoAss> dedxCollH;
      dedxCollH.getByLabel(ev, "dedxHitInfo");
      if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");continue;}

      fwlite::Handle< std::vector<reco::Track> > trackCollHandle;
      trackCollHandle.getByLabel(ev,"RefitterForDeDx");
      if(!trackCollHandle.isValid()){
         trackCollHandle.getByLabel(ev,"generalTracks");
	    if (!trackCollHandle.isValid()){
	       printf("Invalid trackCollHandle\n");
	       continue;
	    }
      }

      for(unsigned int c=0;c<trackCollHandle->size();c++){
         //basic track quality cuts
         reco::TrackRef track = reco::TrackRef( trackCollHandle.product(), c );
         if(track.isNull())continue;
         if(track->chi2()/track->ndof()>5 )continue;  //WAS >1
         if(track->found()<8)continue;

         //load dEdx informations
         const DeDxHitInfo* dedxHits = NULL;
         DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
         if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
         if(!dedxHits)continue;

         //hit level dEdx information (only done for MIPs)
         if(track->p() > 5){
            for(unsigned int h=0;h<dedxHits->size();h++){
                DetId detid(dedxHits->detId(h));
                double scaleFactor = dEdxSF;
                double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
                double ChargeOverPathlength = scaleFactor*Norm*dedxHits->charge(h)/dedxHits->pathlength(h);

                int moduleGeometry = 1; //0 will be used for the pixels since strips geom start at 1
                if(detid.subdetId()>=3){ SiStripDetId SSdetId(detid); moduleGeometry = SSdetId.moduleGeometry();}

                for(unsigned int R=0;R<results.size();R++){ 
                   if(!results[R]->isHit) continue; //only consider results related to hit info here
                   if(!results[R]->usePixel and detid.subdetId() <3)continue; // skip pixels
                   if(!results[R]->useStrip and detid.subdetId()>=3)continue; // skip strips
         
                   results[R]->Charge_Vs_Path->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor*dedxHits->charge(h)/(dedxHits->pathlength(h)*10)); 
                   if(detid.subdetId()>=3)results[R]->Charge_Vs_FS[moduleGeometry]->Fill(dedxHits->stripCluster(h)->firstStrip(),  dedxHits->charge(h)); 
                   results[R]->HHit->Fill(ChargeOverPathlength);
                }
             }
          }

          for(unsigned int R=0;R<results.size();R++){ 
             if(!results[R]->isEstim and !results[R]->isDiscrim ) continue; //only consider results related to estimator/discriminator variables here

             DeDxData* dedxObj = computedEdx(dedxHits, dEdxSF, results[R]->dEdxTemplates,             results[R]->usePixel, useClusterCleaning, false, false, results[R]->TrackerGains, results[R]->useStrip );

             results[R]->HdedxVsP    ->Fill(track->p(), dedxObj->dEdx() );
             results[R]->HdedxVsQP   ->Fill(track->p()*track->charge(), dedxObj->dEdx() );
             results[R]->HdedxVsP_NS ->Fill(track->p(), dedxObj->dEdx(), dedxObj->numberOfSaturatedMeasurements() );
             results[R]->HdedxVsPHSCP->Fill(track->pt(), dedxObj->dEdx());

             if(track->pt()>10 && track->pt()<45 && dedxObj->numberOfMeasurements()>=(results[R]->useStrip?7:3) ){
               results[R]->HdedxVsEtaProfile->Fill(track->eta(), dedxObj->dEdx() );
               results[R]->HdedxVsEta->Fill(track->eta(), dedxObj->dEdx() );
               results[R]->HNOMVsEtaProfile->Fill(track->eta(),dedxObj->numberOfMeasurements() );
               results[R]->HNOSVsEtaProfile->Fill(track->eta(),dedxObj->numberOfSaturatedMeasurements() );
               results[R]->HNOMSVsEtaProfile->Fill(track->eta(),dedxObj->numberOfMeasurements() - dedxObj->numberOfSaturatedMeasurements() );
             }

             if(fabs(track->eta())>2.1)continue;
             if((int)dedxObj->numberOfMeasurements()<(results[R]->useStrip?10:3))continue;
             if(track->p()>5 && track->p()<40){
                results[R]->HdedxMIP->Fill(dedxObj->dEdx());
                results[R]->HP->Fill(track->p());
             } else if (track->pt() > 40)
                results[R]->HdedxSIG->Fill(dedxObj->dEdx());

             if(fabs(track->eta())<0.4)results[R]->HdedxVsPProfile->Fill(track->p(), dedxObj->dEdx() );

             if(results[R]->isEstim){  //mass can only be computed for dEdx estimators
                double Mass = GetMass(track->p(),dedxObj->dEdx(), false);
                if(dedxObj->dEdx()>4.0 && track->p()<3.0){
                   results[R]->HMass->Fill(Mass);
                }
             }
          }
      }
   }printf("\n");

   OutputHisto->Write();
   OutputHisto->Close();  
}
