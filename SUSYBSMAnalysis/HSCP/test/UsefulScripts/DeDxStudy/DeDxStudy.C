
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
   TH3D* Charge_Vs_Path;
   TH1D* HdedxMIP;
   TH2D* HdedxVsPHSCP;
   TH2D* HdedxVsP;
   TH2D* HdedxVsPM;
   TH2D* HdedxVsQP;
   TProfile* HdedxVsPProfile;
   TProfile* HdedxVsEtaProfile;
   TH2D* HdedxVsEta;
   TProfile* HNOSVsEtaProfile;
   TProfile* HNOMVsEtaProfile;
   TProfile* HNOMSVsEtaProfile;
   TH1D* HMass;
   TH1D* HP;
   TH1D* HHit; 
   TH2D* HIasVsP;
   TH2D* HIasVsPM;
   TH1D* HIasMIP;

   dEdxStudyObj(string saveName){
      Name = saveName;

      bool isDiscrim = false;  if(saveName.find("Ias")!=std::string::npos)isDiscrim = true;


      string HistoName;
      HistoName = saveName + "_ChargeVsPath"; Charge_Vs_Path        = new TH3D(    HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
      HistoName = saveName + "_MIP";          HdedxMIP              = new TH1D(    HistoName.c_str(), HistoName.c_str(),  200, 0, isDiscrim?1.0:20);
      HistoName = saveName + "_dedxVsPHSCP";  HdedxVsPHSCP          = new TH2D(    HistoName.c_str(), HistoName.c_str(), 3000, 0, 2000,1500,0, isDiscrim?1.0:15);
      HistoName = saveName + "_dedxVsP";      HdedxVsP              = new TH2D(    HistoName.c_str(), HistoName.c_str(), 3000, 0, 30,1500,0, isDiscrim?1.0:15);
      HistoName = saveName + "_dedxVsPM";     HdedxVsPM             = new TH2D(    HistoName.c_str(), HistoName.c_str(), 3000, 0, 30,1500,0, isDiscrim?1.0:15);
      HistoName = saveName + "_dedxVsQP";     HdedxVsQP             = new TH2D(    HistoName.c_str(), HistoName.c_str(), 6000, -30, 30,1500,0, isDiscrim?1.0:25);
      HistoName = saveName + "_Profile";      HdedxVsPProfile       = new TProfile(HistoName.c_str(), HistoName.c_str(), 100, 0,100);
      HistoName = saveName + "_Eta";          HdedxVsEtaProfile     = new TProfile(HistoName.c_str(), HistoName.c_str(), 100,-3,  3);
      HistoName = saveName + "_Eta2D";        HdedxVsEta            = new TH2D    (HistoName.c_str(), HistoName.c_str(), 100,-3,  3, 1000,0, isDiscrim?1.0:5);
      HistoName = saveName + "_NOS";          HNOSVsEtaProfile      = new TProfile(HistoName.c_str(), HistoName.c_str(), 100,-3,  3);
      HistoName = saveName + "_NOM";          HNOMVsEtaProfile      = new TProfile(HistoName.c_str(), HistoName.c_str(), 100,-3,  3);
      HistoName = saveName + "_NOMS";         HNOMSVsEtaProfile     = new TProfile(HistoName.c_str(), HistoName.c_str(), 100,-3,  3);
      HistoName = saveName + "_Mass";         HMass                 = new TH1D(    HistoName.c_str(), HistoName.c_str(),  500, 0, 10);
      HistoName = saveName + "_P";            HP                    = new TH1D(    HistoName.c_str(), HistoName.c_str(),  500, 0, 100);
      HistoName = saveName + "_Hit";          HHit                  = new TH1D(    HistoName.c_str(), HistoName.c_str(),  200, 0, 20); 
   }
};






void DeDxStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root")
{
  if(DIRNAME=="COMPILE") return;

   system("mkdir -p pictures");

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

   std::vector<string> FileName;
   bool isData = true;
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
   fwlite::ChainEvent ev(FileName);

   TH3F* dEdxTemplates = NULL;
   if(isData){    dEdxTemplates = loadDeDxTemplate(DIRNAME + "/../../../data/Data7TeV_Deco_SiStripDeDxMip_3D_Rcd.root");
   }else{         dEdxTemplates = loadDeDxTemplate(DIRNAME + "/../../../data/MC7TeV_Deco_SiStripDeDxMip_3D_Rcd.root");
   }


   string studies[] = {"harm2", "trunc40", "harm2_raw", "trunc40_raw", "Ias"};
   unsigned int Nstudies = sizeof(studies)/sizeof(string);
   dEdxStudyObj** results = new dEdxStudyObj*[Nstudies];
   for(unsigned int R=0;R<Nstudies;R++){ results[R] = new dEdxStudyObj(studies[R]); }

   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");

   std::unordered_map<unsigned int,double> TrackerGains;
   LoadDeDxCalibration(TrackerGains, DIRNAME+"/Gains.root");

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
          reco::TrackRef track = reco::TrackRef( trackCollHandle.product(), c );

          if(track->chi2()/track->ndof()>5 )continue;  //WAS >1
          if(track->found()<8)continue;


          const DeDxHitInfo* dedxHits = NULL;
          if(TypeMode!=3 && !track.isNull()) {
             DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
             if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
          }
          if(!dedxHits)continue;


         for(unsigned int h=0;h<dedxHits->size();h++){
             DetId detid(dedxHits->detId(h));
             if(detid.subdetId()<3)continue; // skip pixels
               
             double scaleFactor = 1.0;
             double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
             double ChargeOverPathlength = scaleFactor*Norm*dedxHits->charge(h)/dedxHits->pathlength(h);

             SiStripDetId SSdetId(detid);

             
    	     if (track->p() > 5){
                 for(unsigned int R=0;R<Nstudies;R++){ results[R]->Charge_Vs_Path->Fill (SSdetId.moduleGeometry(), dedxHits->pathlength(h)*10, dedxHits->charge(h)/(dedxHits->pathlength(h)*10)); }
             }
             for(unsigned int R=0;R<Nstudies;R++){ results[R]->HHit->Fill(ChargeOverPathlength); }
          }


          for(unsigned int R=0;R<Nstudies;R++){ 
             DeDxData* dedxObj = NULL;
             if      (studies[R]=="harm2"      ){ dedxObj = computedEdx(dedxHits, dEdxSF, NULL,          false, useClusterCleaning, false, false, &TrackerGains );  
             }else if(studies[R]=="harm2_raw"  ){ dedxObj = computedEdx(dedxHits, dEdxSF, NULL,          false, useClusterCleaning, false, false, NULL          ); 
             }else if(studies[R]=="trunc40"    ){ dedxObj = computedEdx(dedxHits, dEdxSF, NULL,          false, useClusterCleaning, true,  false, &TrackerGains );  
             }else if(studies[R]=="trunc40_raw"){ dedxObj = computedEdx(dedxHits, dEdxSF, NULL,          false, useClusterCleaning, true,  false, NULL          );  
             }else if(studies[R]=="Ias"        ){ dedxObj = computedEdx(dedxHits, dEdxSF, dEdxTemplates, false, useClusterCleaning, false, false );                 
             }else{ printf("Unknown case: %s!  Exit here.\n", studies[R].c_str());  exit(0);
             }

             results[R]->HdedxVsPHSCP->Fill(track->pt(), dedxObj->dEdx());

             if(track->pt()>20 && track->pt()<40 && dedxObj->numberOfMeasurements()>6 ){
               results[R]->HdedxVsEtaProfile->Fill(track->eta(), dedxObj->dEdx() );
               results[R]->HdedxVsEta->Fill(track->eta(), dedxObj->dEdx() );
               results[R]->HNOMVsEtaProfile->Fill(track->eta(),dedxObj->numberOfMeasurements() );
               results[R]->HNOSVsEtaProfile->Fill(track->eta(),dedxObj->numberOfSaturatedMeasurements() );
               results[R]->HNOMSVsEtaProfile->Fill(track->eta(),dedxObj->numberOfMeasurements() - dedxObj->numberOfSaturatedMeasurements() );
             }

             if(fabs(track->eta())>2.1)continue;
             if((int)dedxObj->numberOfMeasurements()<10)continue;
             if(track->p()>5 && track->p()<40){
                results[R]->HdedxMIP->Fill(dedxObj->dEdx());
                results[R]->HP->Fill(track->p());
             }
             results[R]->HdedxVsP ->Fill(track->p(), dedxObj->dEdx() );
             results[R]->HdedxVsQP->Fill(track->p()*track->charge(), dedxObj->dEdx() );

             if(fabs(track->eta())<0.4)results[R]->HdedxVsPProfile->Fill(track->p(), dedxObj->dEdx() );
             double Mass = GetMass(track->p(),dedxObj->dEdx(), false);
             if(dedxObj->dEdx()>4.0 && track->p()<3.0){
                results[R]->HMass->Fill(Mass);
                if(isnan((float)Mass) || Mass<0.94-0.3 || Mass>0.94+0.3)continue;
                results[R]->HdedxVsPM ->Fill(track->p(), dedxObj->dEdx() );
             }
          }
      }
   }printf("\n");
   OutputHisto->Write();
   OutputHisto->Close();  
}
