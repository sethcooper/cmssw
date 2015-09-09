// Original Author:  Loic Quertenmont
 

#include <exception>
#include <vector>
#include <unordered_map>

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
#include "TVector3.h"

namespace reco    { class Vertex; class Track; class Muon; class MuonTimeExtra;   }
namespace fwlite  { class ChainEvent;}
namespace edm     { class TriggerResultsByName; class InputTag;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCIndexer.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

using namespace fwlite;
using namespace reco;
using namespace std;
using namespace edm;

#include "../../AnalysisCode/Analysis_TOFUtility.h"

#endif


bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl);


std::map<string, TH1D*> HChamber_Timing;
TH1D* getHisto(string name){
   if(HChamber_Timing.find(name) == HChamber_Timing.end()){  HChamber_Timing[name.c_str()] = new TH1D(name.c_str(), name.c_str(), 800, -100, 100);   }
   return HChamber_Timing[name];
}

void MuonTimingStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root", string DIRECTORY="")
{
  if(DIRNAME=="COMPILE") return;

   std::vector<string> FileName;
   if(INPUT.find(".root")<std::string::npos){
      char* pch=strtok(&INPUT[0],",");
      while (pch!=NULL){
         FileName.push_back(pch);    
         pch=strtok(NULL,",");
      }
   }else{
      printf("There is no .root file in the input file list\nStop the job here\n");
      exit(0);
   }
 
   moduleGeom::loadGeometry(DIRNAME+"/../../../data/CMS_GeomTree.root");
   muonTimingCalculator tofCalculator;
   tofCalculator.loadTimeOffset(DIRNAME+"/../../../data/MuonTimeOffset.txt");

   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");  //File must be opened before the histogram are created
   string HistoName; 

   HistoName = "CSC_iBeta_AOD";  TH1D* HCSC_iBetaAOD   = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "CSC_iBeta_FLY0"; TH1D* HCSC_iBetaFLY0  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "CSC_iBeta_FLY1"; TH1D* HCSC_iBetaFLY1  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "CSC_iBeta_FLY2"; TH1D* HCSC_iBetaFLY2  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);

   HistoName = "DT_iBeta_AOD";  TH1D* HDT_iBetaAOD   = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "DT_iBeta_FLY0"; TH1D* HDT_iBetaFLY0  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "DT_iBeta_FLY1"; TH1D* HDT_iBetaFLY1  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "DT_iBeta_FLY2"; TH1D* HDT_iBetaFLY2  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);

   HistoName = "iBeta_AOD";  TH1D* H_iBetaAOD   = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "iBeta_FLY0"; TH1D* H_iBetaFLY0  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "iBeta_FLY1"; TH1D* H_iBetaFLY1  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);
   HistoName = "iBeta_FLY2"; TH1D* H_iBetaFLY2  = new TH1D(HistoName.c_str(), HistoName.c_str(), 150, 0.5, 2.0);


   TDirectory* dir = OutputHisto;
   if(DIRECTORY!=""){
      TDirectory::AddDirectory(kTRUE);
      TH1::AddDirectory(kTRUE);
      dir = OutputHisto->mkdir(DIRECTORY.c_str(), DIRECTORY.c_str());
      dir->cd();
   }

   HistoName = "CSC_Timing";  TH1D* HCSC_Timing   = new TH1D(HistoName.c_str(), HistoName.c_str(), 800, -100, 100);
   HistoName = "DT_Timing";   TH1D* HDT_Timing    = new TH1D(HistoName.c_str(), HistoName.c_str(), 800, -100, 100);

   unsigned int CurrentRun = 0;
   printf("Progressing Bar           :0%%       20%%       40%%       60%%       80%%       100%%\n");
   for(unsigned int f=0;f<FileName.size();f++){
     TFile* file = TFile::Open(FileName[f].c_str() );
     fwlite::Event ev(file);

     dir->cd(); //make sure all new histograms are saved in the output file

     printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)FileName.size());
     int treeStep(std::max(1,int(ev.size()/50))), iev=0;
     for(ev.toBegin(); !ev.atEnd(); ++ev){ iev++;
         if(iev%treeStep==0){printf(".");fflush(stdout);}

         if(CurrentRun != ev.eventAuxiliary().run()){
            CurrentRun = ev.eventAuxiliary().run();
            tofCalculator.setRun(CurrentRun);
         }

         fwlite::Handle < std::vector<reco::Muon> > muonCollHandle;
         muonCollHandle.getByLabel(ev, "ALCARECOMuAlCalIsolatedMu", "SelectedMuons");
         if(!muonCollHandle.isValid()){
            muonCollHandle.getByLabel(ev, "muons");
            if(!muonCollHandle.isValid()){printf("Muon Collection not found!\n"); continue;}
         }

         fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
         TOFCSCCollH.getByLabel(ev, "muons","csc");
         if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");return;}

         fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
         TOFDTCollH.getByLabel(ev, "muons","dt");
         if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");return;}

         fwlite::Handle<MuonTimeExtraMap> TOFCollH;
         TOFCollH.getByLabel(ev, "muons","combined");
         if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}

         fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
         CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
         if(!CSCSegmentCollHandle.isValid()){printf("CSC Segment Collection not found!\n"); continue;}
         const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;

         fwlite::Handle<DTRecSegment4DCollection> DTSegmentCollHandle;
         DTSegmentCollHandle.getByLabel(ev, "dt4DSegments");
         if(!DTSegmentCollHandle.isValid()){printf("DT Segment Collection not found!\n"); continue;}
         const DTRecSegment4DCollection& DTCSegmentColl = *DTSegmentCollHandle;

         fwlite::Handle < std::vector<reco::Vertex> > vertexCollHandle;
         vertexCollHandle.getByLabel(ev, "offlinePrimaryVertices");
         if(!vertexCollHandle.isValid()){printf("Vertex Collection not found!\n"); continue;}
         const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;

         for(unsigned int c=0;c<muonCollHandle->size();c++){
            //basic track quality cuts
            reco::MuonRef muon = reco::MuonRef( muonCollHandle.product(), c );
            if(muon.isNull())continue;
            if(!muon->isGlobalMuon())continue;
            if(muon->pt()<10)continue;
            if(isCompatibleWithCosmic(muon->track(), vertexColl))continue;


            //PLOT SEGMENT T0            
            vector<const CSCSegment*>& cscSegs = tofCalculator.matchCSC(*muon->standAloneMuon(), CSCSegmentColl);

            for(unsigned int ic=0;ic<cscSegs.size();ic++){  
               HCSC_Timing->Fill(cscSegs[ic]->time());
               getHisto(moduleGeom::getChamberName(cscSegs[ic]->cscDetId().rawId()))->Fill(cscSegs[ic]->time());
               getHisto(moduleGeom::getStationName(cscSegs[ic]->cscDetId().rawId()))->Fill(cscSegs[ic]->time());
            }

            vector<const DTRecSegment4D*>& dtSegs = tofCalculator.matchDT(*muon->standAloneMuon(), DTCSegmentColl);
            for(unsigned int id=0;id<dtSegs.size();id++){
               for(int phi=0;phi<2;phi++){
                  const DTRecSegment2D* segm=NULL;
                  if(phi) segm = dynamic_cast<const DTRecSegment2D*>(dtSegs[id]->phiSegment()); 
                  else segm = dynamic_cast<const DTRecSegment2D*>(dtSegs[id]->zSegment());
                  if(!segm)continue;
 
                  HDT_Timing->Fill(segm->t0());
                  getHisto(moduleGeom::getChamberName(segm->geographicalId()))->Fill(segm->t0());
                  getHisto(moduleGeom::getStationName(segm->geographicalId()))->Fill(segm->t0());
               }
            }

            //PLOT MUON iBETA
            const reco::MuonTimeExtra* csctof = &TOFCSCCollH->get(muon.key());
            const reco::MuonTimeExtra* dttof = &TOFDTCollH->get(muon.key());
            const reco::MuonTimeExtra* tof = &TOFCollH->get(muon.key());

            HDT_iBetaAOD ->Fill(dttof->inverseBeta());
            HCSC_iBetaAOD->Fill(csctof->inverseBeta());
            H_iBetaAOD->Fill(tof->inverseBeta());

            tofCalculator.tmSeq.clear();
            HDT_iBetaFLY0 ->Fill(tofCalculator.iBetaFromDT (muon, 0));
            HCSC_iBetaFLY0->Fill(tofCalculator.iBetaFromCSC(muon, 0));
            reco::MuonTimeExtra tof0 = tofCalculator.fillTimeFromMeasurements();
            H_iBetaFLY0->Fill(tof0.inverseBeta());

            tofCalculator.tmSeq.clear();
            HDT_iBetaFLY1 ->Fill(tofCalculator.iBetaFromDT (muon, 1));
            HCSC_iBetaFLY1->Fill(tofCalculator.iBetaFromCSC(muon, 1));
            reco::MuonTimeExtra tof1 = tofCalculator.fillTimeFromMeasurements();
            H_iBetaFLY1->Fill(tof1.inverseBeta());

            tofCalculator.tmSeq.clear();
            HDT_iBetaFLY2 ->Fill(tofCalculator.iBetaFromDT (muon, 2));
            HCSC_iBetaFLY2->Fill(tofCalculator.iBetaFromCSC(muon, 2));
            reco::MuonTimeExtra tof2 = tofCalculator.fillTimeFromMeasurements();
            H_iBetaFLY2->Fill(tof2.inverseBeta());
 
            printf("%f vs %f  (%i vs %i\n", tof->inverseBeta(), tof2.inverseBeta(), tof->nDof(), tof2.nDof());
            for (unsigned int i=0;i<tofCalculator.tmSeq.size();i++){
               printf("%6.2E ", tofCalculator.tmSeq[i].weightInvBeta);
            }printf("\n");


         }//muon
      }printf("\n");
      delete file;
   }

   OutputHisto->Write();
   OutputHisto->Close();  
}

bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl){
   for (unsigned int vertex_i=0;vertex_i<vertexColl.size();vertex_i++){
      if(fabs(track->dz (vertexColl[vertex_i].position())) < 0.5 && fabs(track->dxy(vertexColl[vertex_i].position())) < 0.2)return false;
   }
   return true;
}











