
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
      FileName.push_back(INPUT);    
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

   //system("mkdir -p pictures/");
   double P_Min               = 1;
   double P_Max               = 15;
   int    P_NBins             = 14  ;
   double Path_Min            = 0.2 ;
   double Path_Max            = 1.6 ;
   int    Path_NBins          = 42  ;
   double Charge_Min          = 0   ;
   double Charge_Max          = 5000;
   int    Charge_NBins        = 500 ;

   TH3D* Charge_Vs_Path = new TH3D ("Charge_Vs_Path"     , "Charge_Vs_Path" , P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);

   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");

   for(int LOOP=0;LOOP<=1;LOOP++){
      string saveName = "harm2";
      if(LOOP==1) saveName = "trunc40";

      TH1D* HdedxMIP          = new TH1D(    (saveName + "_MIP"    ).c_str(), "MIP"    ,  200, 0, 20);
      TH2D* HdedxVsPHSCP      = new TH2D(    (saveName + "_dedxVsPHSCP").c_str(), "dedxVsPHSCP", 3000, 0, 2000,1500,0,15);
      TH2D* HdedxVsP          = new TH2D(    (saveName + "_dedxVsP").c_str(), "dedxVsP", 3000, 0, 30,1500,0,15);
      TH2D* HdedxVsPM         = new TH2D(    (saveName + "_dedxVsPM").c_str(), "dedxVsPM", 3000, 0, 30,1500,0,15);
      TH2D* HdedxVsQP         = new TH2D(    (saveName + "_dedxVsQP").c_str(), "dedxVsQP", 6000, -30, 30,1500,0,25);
      TProfile* HdedxVsPProfile   = new TProfile((saveName + "_Profile").c_str(), "Profile",  100, 0,100);
      TProfile* HdedxVsEtaProfile = new TProfile((saveName + "_Eta"    ).c_str(), "Eta"    ,  100,-3,  3);
      TH2D* HdedxVsEta        = new TH2D    ((saveName + "_Eta2D"  ).c_str(), "Eta"    ,  100,-3,  3, 1000,0,5);
      TProfile* HNOSVsEtaProfile  = new TProfile((saveName + "_NOS"    ).c_str(), "NOS"    ,  100,-3,  3);
      TProfile* HNOMVsEtaProfile  = new TProfile((saveName + "_NOM"    ).c_str(), "NOM"    ,  100,-3,  3);
      TProfile* HNOMSVsEtaProfile = new TProfile((saveName + "_NOMS"    ).c_str(), "NOMS"    ,  100,-3,  3);
      TH1D* HMass             = new TH1D(    (saveName + "_Mass"   ).c_str(), "Mass"   ,  500, 0, 10);
      TH1D* HP                = new TH1D(    (saveName + "_P"      ).c_str(), "P"      ,  500, 0, 100);

      TH1D* HHit          = new TH1D(    (saveName + "_Hit"      ).c_str(), "P"      ,  200, 0, 20);
 
      TH2D* HIasVsP          = new TH2D(    (saveName + "_IasVsP").c_str(), "IasVsP", 3000, 0, 30,1500,0,1);
      TH2D* HIasVsPM         = new TH2D(    (saveName + "_IasVsPM").c_str(), "IasVsPM", 3000, 0, 30,1500,0,1);
      TH1D* HIasMIP          = new TH1D(    (saveName + "_IasMIP"    ).c_str(), "IasMIP"    ,  1000, 0, 1);



      printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
      printf("Looping on Tree              :");
      int TreeStep = ev.size()/50;if(TreeStep==0)TreeStep=1;
      for(Long64_t e=0;e<ev.size();e++){
   //      if(e>10)break;
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
		if (track->p() > 5) Charge_Vs_Path->Fill (SSdetId.moduleGeometry,
				dedxHits->pathlength(h),
				dedxHits->charge(h)/dedxHits->pathlength(h));
                HHit->Fill(ChargeOverPathlength);
             }

             DeDxData* dedxSObj = computedEdx(dedxHits, dEdxSF, dEdxTemplates, false, useClusterCleaning, false, false );
             DeDxData* dedxMObj = computedEdx(dedxHits, dEdxSF, NULL,          false, useClusterCleaning, false, LOOP==0?false:true  );

             HdedxVsPHSCP->Fill(track->pt(), dedxMObj->dEdx());

             if(track->pt()>20 && track->pt()<40 && dedxMObj->numberOfMeasurements()>6 ){
               HdedxVsEtaProfile->Fill(track->eta(), dedxMObj->dEdx() );
               HdedxVsEta->Fill(track->eta(), dedxMObj->dEdx() );
               HNOMVsEtaProfile->Fill(track->eta(),dedxMObj->numberOfMeasurements() );
               HNOSVsEtaProfile->Fill(track->eta(),dedxMObj->numberOfSaturatedMeasurements() );
               HNOMSVsEtaProfile->Fill(track->eta(),dedxMObj->numberOfMeasurements() - dedxMObj->numberOfSaturatedMeasurements() );
             }


             if(fabs(track->eta())>2.1)continue;
             if((int)dedxMObj->numberOfMeasurements()<10)continue;
             if(track->p()>5 && track->p()<40){
                HdedxMIP->Fill(dedxMObj->dEdx());
                HP->Fill(track->p());
                HIasMIP->Fill(dedxSObj->dEdx());
             }
             HdedxVsP ->Fill(track->p(), dedxMObj->dEdx() );
             HdedxVsQP->Fill(track->p()*track->charge(), dedxMObj->dEdx() );
             HIasVsP ->Fill(track->p(), dedxSObj->dEdx() );

             if(fabs(track->eta())<0.4)HdedxVsPProfile->Fill(track->p(), dedxMObj->dEdx() );
             double Mass = GetMass(track->p(),dedxMObj->dEdx(), false);
             if(dedxMObj->dEdx()>4.0 && track->p()<3.0){
                HMass->Fill(Mass);
                if(isnan((float)Mass) || Mass<0.94-0.3 || Mass>0.94+0.3)continue;
                HdedxVsPM ->Fill(track->p(), dedxMObj->dEdx() );
                HIasVsPM ->Fill(track->p(), dedxSObj->dEdx() );
             }



         }
      }printf("\n");


      HdedxMIP->Scale(1.0/HdedxMIP->Integral() );
      HHit->Scale(1.0/HHit->Integral() );

   }

   Charge_Vs_Path->SaveAs(("ChargeVsPath"+OUTPUT).c_str());
   OutputHisto->Write();
   OutputHisto->Close();  
}
