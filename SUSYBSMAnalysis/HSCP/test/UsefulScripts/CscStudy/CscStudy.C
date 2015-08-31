
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
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;


#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"

#endif

//Global variables for Segment matcher
//         DTsegments = cms.InputTag("dt4DSegments"),
double DTradius      = 0.01;
bool   TightMatchDT  = false;
bool   cscTightMatch = true;


bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl);
vector<const CSCSegment*> matchCSC(const reco::Track& muon, const CSCSegmentCollection& CSCSegmentColl);

void CscStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root")
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

   bool isData   = !(INPUT.find("MC")!=string::npos);
   bool isSignal = false;
   bool removeCosmics = true;
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
      isData   = (sample.Type==0);
      isSignal = (sample.Type==2);
      GetInputFiles(sample, BaseDirectory, FileName, 0);
   }


   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");  //File must be opened before the histogram are created

   string HistoName; 
   HistoName = "CSC_Timing";  TH1D* HCSC_Timing   = new TH1D(HistoName.c_str(), HistoName.c_str(), 2000, -100, 100);


   printf("Progressing Bar           :0%%       20%%       40%%       60%%       80%%       100%%\n");
   for(unsigned int f=0;f<FileName.size();f++){
     TFile* file = TFile::Open(FileName[f].c_str() );
     fwlite::Event ev(file);
     printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)FileName.size());
     int treeStep(ev.size()/50), iev=0;
     for(ev.toBegin(); !ev.atEnd(); ++ev){ iev++;
         if(iev%treeStep==0){printf(".");fflush(stdout);}

         fwlite::Handle < std::vector<reco::Muon> > muonCollHandle;
         muonCollHandle.getByLabel(ev, "muons");
         if(!muonCollHandle.isValid()){printf("Muon Collection not found!\n"); continue;}
         const std::vector<reco::Muon>& muonColl = *muonCollHandle;

         fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
         CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
         if(!CSCSegmentCollHandle.isValid()){printf("CSC Segment Collection not found!\n"); continue;}
         const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;

         fwlite::Handle < std::vector<reco::Vertex> > vertexCollHandle;
         vertexCollHandle.getByLabel(ev, "offlinePrimaryVertices");
         if(!vertexCollHandle.isValid()){printf("Vertex Collection not found!\n"); continue;}
         const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
         if(vertexColl.size()<1){printf("NO VERTICES\n"); continue;}

         for(unsigned int c=0;c<muonCollHandle->size();c++){
            //basic track quality cuts
            reco::MuonRef muon = reco::MuonRef( muonCollHandle.product(), c );
            if(muon.isNull())continue;
            if(!muon->isGlobalMuon())continue;
            if(muon->pt()<20)continue;
            if(isCompatibleWithCosmic(muon->track(), vertexColl))continue;

            vector<const CSCSegment*> cscSegs = matchCSC(*muon->combinedMuon(), CSCSegmentColl);
            for(unsigned int ic=0;ic<cscSegs.size();ic++){  
//               printf("Segment %i --> Time = %f\n", ic, cscSegs[ic]->time());
               HCSC_Timing->Fill(cscSegs[ic]->time());
            }


         }
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



//Fwlite version of the code found in https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_6_X/RecoMuon/TrackingTools/src/MuonSegmentMatcher.cc
vector<const CSCSegment*> matchCSC(const reco::Track& muon, const CSCSegmentCollection& CSCSegmentColl)
{
  vector<const CSCSegment*> pointerToCSCSegments;

  double matchRatioCSC=0;
  double CSCXCut = 0.001;
  double CSCYCut = 0.001;
  double countMuonCSCHits = 0;

  for(CSCSegmentCollection::const_iterator segmentCSC = CSCSegmentColl.begin(); segmentCSC != CSCSegmentColl.end(); segmentCSC++) {
    double CSCcountAgreeingHits=0;
    if ( !segmentCSC->isValid()) continue; 

    const vector<CSCRecHit2D>& CSCRechits2D = segmentCSC->specificRecHits();
    countMuonCSCHits = 0;
    CSCDetId myChamber((*segmentCSC).geographicalId().rawId());

    bool segments = false;

    for(trackingRecHit_iterator hitC = muon.recHitsBegin(); hitC != muon.recHitsEnd(); ++hitC) {
      if (!(*hitC)->isValid()) continue; 
      if ( (*hitC)->geographicalId().det() != DetId::Muon ) continue; 
      if ( (*hitC)->geographicalId().subdetId() != MuonSubdetId::CSC ) continue;
      if (!(*hitC)->isValid()) continue;
      if ( (*hitC)->recHits().size()>1) segments = true;

      //DETECTOR CONSTRUCTION
      DetId id = (*hitC)->geographicalId();
      CSCDetId cscDetIdHit(id.rawId());

      if (segments) {
	if(!(myChamber.rawId()==cscDetIdHit.rawId())) continue; 

        // and compare the local positions
        LocalPoint positionLocalCSC = (*hitC)->localPosition();
	LocalPoint segLocalCSC = segmentCSC->localPosition();
	if ((fabs(positionLocalCSC.x()-segLocalCSC.x())<CSCXCut) && 
	    (fabs(positionLocalCSC.y()-segLocalCSC.y())<CSCYCut)) 
	  pointerToCSCSegments.push_back(&(*segmentCSC)); 
        continue;
      }

      if(!(cscDetIdHit.ring()==myChamber.ring())) continue;
      if(!(cscDetIdHit.station()==myChamber.station())) continue;
      if(!(cscDetIdHit.endcap()==myChamber.endcap())) continue;
      if(!(cscDetIdHit.chamber()==myChamber.chamber())) continue;

      countMuonCSCHits++;

      LocalPoint positionLocalCSC = (*hitC)->localPosition();
	
      for (vector<CSCRecHit2D>::const_iterator hiti=CSCRechits2D.begin(); hiti!=CSCRechits2D.end(); hiti++) {

	if ( !hiti->isValid()) continue; 
	CSCDetId cscDetId((hiti->geographicalId()).rawId());
		
	if ((*hitC)->geographicalId().rawId()!=(hiti->geographicalId()).rawId()) continue;

	LocalPoint segLocalCSC = hiti->localPosition();
	//		cout<<"Layer Id (MuonHit) =  "<<cscDetIdHit<<" Muon Local Position (det frame) "<<positionLocalCSC <<endl;
	//		cout<<"Layer Id  (CSCHit) =  "<<cscDetId<<"  Hit Local Position (det frame) "<<segLocalCSC <<endl;
	if((fabs(positionLocalCSC.x()-segLocalCSC.x())<CSCXCut) && 
	   (fabs(positionLocalCSC.y()-segLocalCSC.y())<CSCYCut)) {
	  CSCcountAgreeingHits++;
	  //		  cout << "   Matched." << endl;
	}  
      }//End 2D rechit iteration
    }//End muon hit iteration
    
    matchRatioCSC = countMuonCSCHits == 0 ? 0 : CSCcountAgreeingHits/countMuonCSCHits;
		
    if ((matchRatioCSC>0.9) && ((countMuonCSCHits>1) || !cscTightMatch)) pointerToCSCSegments.push_back(&(*segmentCSC));

  } //End CSC Segment Iteration 

  return pointerToCSCSegments;
}


