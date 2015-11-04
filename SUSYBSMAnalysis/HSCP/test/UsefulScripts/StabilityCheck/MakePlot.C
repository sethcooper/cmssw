
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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TCutG.h"


#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"


std::map<unsigned int, double> RunToIntLumi;

bool LoadLumiToRun()
{
   float TotalIntLuminosity = 0;

   FILE* pFile = fopen("out.txt","r");
   if(!pFile){
      printf("Not Found: %s\n","out.txt");
      return false;
   }

   unsigned int Run; float IntLumi;
   unsigned int DeliveredLs; double DeliveredLumi;
   char Line[2048], Tmp1[2048], Tmp2[2048], Tmp3[2048];
   while ( ! feof (pFile) ){
     fscanf(pFile,"%s\n",Line);
//     printf("%s\n",Line);
     for(unsigned int i=0;Line[i]!='\0';i++){if(Line[i]==',')Line[i]=' ';} 
     sscanf(Line,"%d %s %s %s %f\n",&Run,Tmp1,Tmp2,Tmp3,&IntLumi);
     TotalIntLuminosity+= IntLumi/1000000.0;
//     printf("%6i --> %f/pb   (%s | %s | %s)\n",Run,TotalIntLuminosity,Tmp1,Tmp2,Tmp3);

     RunToIntLumi[Run] = TotalIntLuminosity;
   }
   fclose(pFile);
   return true;
}


TGraph* ConvertFromRunToIntLumi(TProfile* Object, const char* DrawOption, string YLabel, double YRange_Min=3.1, double YRange_Max=3.7){
   TGraphErrors* graph = new TGraphErrors(Object->GetXaxis()->GetNbins());
   for(int i=1;i<Object->GetXaxis()->GetNbins()+1;i++){
      int RunNumber;
      sscanf(Object->GetXaxis()->GetBinLabel(i),"%d",&RunNumber);
      graph->SetPoint(i-1, RunToIntLumi[RunNumber], Object->GetBinContent(i));
      graph->SetPointError(i-1, 0.0*RunToIntLumi[RunNumber], Object->GetBinError(i));
   }
   graph->Draw(DrawOption);
   graph->SetTitle("");
   graph->GetYaxis()->SetTitle(Object->GetYaxis()->GetTitle());
   graph->GetYaxis()->SetTitleOffset(1.10);
   graph->GetXaxis()->SetTitle("Int. Luminosity (/pb)");
   graph->GetYaxis()->SetTitle(YLabel.c_str());
   graph->SetMarkerColor(Object->GetMarkerColor());
   graph->SetMarkerStyle(Object->GetMarkerStyle());
   graph->GetXaxis()->SetNdivisions(510);
   if(YRange_Min!=YRange_Max)graph->GetYaxis()->SetRangeUser(YRange_Min,YRange_Max);
   return graph;
}

void MakedEdxPlot()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.09);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;

   TFile* InputFile = new TFile("pictures/Histos.root");

   std::vector<string> runList;
   TList* ObjList = InputFile->GetListOfKeys();
   for(int i=0;i<ObjList->GetSize();i++){
      TObject* tmp = GetObjectFromPath(InputFile,ObjList->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TDirectory")){
         runList.push_back(ObjList->At(i)->GetName());
         printf("Add a new run: %s\n", ObjList->At(i)->GetName() );
      }
      delete tmp;
   }

   TProfile* SingleMu_PtProf           = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50PtProf");      
   TProfile* SingleMu_dEdxProf         = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxProf");   
   TProfile* SingleMu_dEdxMProf        = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMSProf       = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMSProf");
   TProfile* SingleMu_dEdxMPProf       = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMPProf");
   TProfile* SingleMu_dEdxMSCProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMSCProf");
   TProfile* SingleMu_dEdxMPCProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMPCProf");
   TProfile* SingleMu_dEdxMSFProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMSFProf");
   TProfile* SingleMu_dEdxMPFProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMPFProf");

   TProfile* SingleMu_NVertProf        = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50NVertProf");

   SingleMu_NVertProf->LabelsDeflate("X");
   SingleMu_NVertProf->LabelsOption("av","X");

/*
   TFile* InputFileLumi166380 = new TFile("pictures/HistosLumi166380.root");
   TFile* InputFileLumi166512 = new TFile("pictures/HistosLumi166512.root");
   TFile* InputFileLumi167807 = new TFile("pictures/HistosLumi167807.root");
   TFile* InputFileLumi167898 = new TFile("pictures/HistosLumi167898.root");

   TProfile* SingleMu_dEdxMProfLumi166380         = (TProfile*)GetObjectFromPath(InputFileLumi166380, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMProfLumi166512         = (TProfile*)GetObjectFromPath(InputFileLumi166512, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMProfLumi167807         = (TProfile*)GetObjectFromPath(InputFileLumi167807, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMProfLumi167898         = (TProfile*)GetObjectFromPath(InputFileLumi167898, "HLT_Mu50dEdxMProf");
*/

   if(LoadLumiToRun()){
      TLegend* leg;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graph =  ConvertFromRunToIntLumi(SingleMu_dEdxMProf  , "A*", "I_{h} (MeV/cm)");
      TGraph* graphS = ConvertFromRunToIntLumi(SingleMu_dEdxMSProf, "*" , "I_{h} (MeV/cm)");
      TGraph* graphP = ConvertFromRunToIntLumi(SingleMu_dEdxMPProf, "*" , "I_{h} (MeV/cm)");
      graphS->SetMarkerColor(2);    graphS->SetMarkerStyle(26);
      graphP->SetMarkerColor(4);    graphP->SetMarkerStyle(32);


      TF1* myfunc = new TF1("Fitgraph" ,"pol1",250,5000);  graph ->Fit(myfunc ,"QN","",250,5000); myfunc ->SetLineWidth(2); myfunc ->SetLineColor(graph ->GetMarkerColor()); myfunc ->Draw("same");
      TF1* myfuncS= new TF1("FitgraphS","pol1",250,5000);  graphS->Fit(myfuncS,"QN","",250,5000); myfuncS->SetLineWidth(2); myfuncS->SetLineColor(graphS->GetMarkerColor()); myfuncS->Draw("same");
      TF1* myfuncP= new TF1("FitgraphP","pol1",250,5000);  graphP->Fit(myfuncP,"QN","",250,5000); myfuncP->SetLineWidth(2); myfuncP->SetLineColor(graphP->GetMarkerColor()); myfuncP->Draw("same");
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip+Pixel)", myfunc ->GetChisquare()/ myfunc ->GetNDF(), myfunc ->GetParameter(0),myfunc ->GetParError(0),myfunc ->GetParameter(1),myfunc ->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip)"      , myfuncS->GetChisquare()/ myfuncS->GetNDF(), myfuncS->GetParameter(0),myfuncS->GetParError(0),myfuncS->GetParameter(1),myfuncS->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Pixel)"      , myfuncP->GetChisquare()/ myfuncP->GetNDF(), myfuncP->GetParameter(0),myfuncP->GetParError(0),myfuncP->GetParameter(1),myfuncP->GetParError(1));
      leg = new TLegend(0.79,0.92,0.79-0.20,0.92 - 3*0.05);     leg->SetFillColor(0);     leg->SetBorderSize(0);
      leg->AddEntry(graph, "dE/dx (Strip+Pixel)" ,"P");
      leg->AddEntry(graphS, "dE/dx (Strip)" ,"P");
      leg->AddEntry(graphP, "dE/dx (Pixel)" ,"P");
      leg->Draw();
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_dEdxM");
      delete c1;  delete leg;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphSC = ConvertFromRunToIntLumi(SingleMu_dEdxMSCProf, "A*", "I_{h} (MeV/cm)");
      TGraph* graphSF = ConvertFromRunToIntLumi(SingleMu_dEdxMSFProf, "*" , "I_{h} (MeV/cm)");
      graphSC->SetMarkerColor(2);    graphSC->SetMarkerStyle(26);
      graphSF->SetMarkerColor(4);    graphSF->SetMarkerStyle(32);
      TF1* myfuncSC= new TF1("FitgraphSC","pol1",250,5000);  graphSC->Fit(myfuncSC,"QN","",250,5000); myfuncSC->SetLineWidth(2); myfuncSC->SetLineColor(graphSC->GetMarkerColor()); myfuncSC->Draw("same");
      TF1* myfuncSF= new TF1("FitgraphSF","pol1",250,5000);  graphSF->Fit(myfuncSF,"QN","",250,5000); myfuncSF->SetLineWidth(2); myfuncSF->SetLineColor(graphSF->GetMarkerColor()); myfuncSF->Draw("same");
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip) |eta|<0.5", myfuncSC->GetChisquare()/ myfuncSC->GetNDF(), myfuncSC->GetParameter(0),myfuncSC->GetParError(0),myfuncSC->GetParameter(1),myfuncSC->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip) |eta|>1.5", myfuncSF->GetChisquare()/ myfuncSF->GetNDF(), myfuncSF->GetParameter(0),myfuncSF->GetParError(0),myfuncSF->GetParameter(1),myfuncSF->GetParError(1));
      leg = new TLegend(0.79,0.92,0.79-0.20,0.92 - 3*0.05);     leg->SetFillColor(0);     leg->SetBorderSize(0);
      leg->AddEntry(graphSC, "dE/dx (Strip) |#eta|<0.5" ,"P");
      leg->AddEntry(graphSF, "dE/dx (Strip) |#eta|>1.5"  ,"P");
      leg->Draw();
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_dEdxMS");
      delete c1; delete leg;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphPC = ConvertFromRunToIntLumi(SingleMu_dEdxMPCProf, "A*", "I_{h} (MeV/cm)");
      TGraph* graphPF = ConvertFromRunToIntLumi(SingleMu_dEdxMPFProf, "*" , "I_{h} (MeV/cm)");
      graphPC->SetMarkerColor(2);    graphPC->SetMarkerStyle(26);
      graphPF->SetMarkerColor(4);    graphPF->SetMarkerStyle(32);
      TF1* myfuncPC= new TF1("FitgraphPC","pol1",250,5000);  graphPC->Fit(myfuncPC,"QN","",250,5000); myfuncPC->SetLineWidth(2); myfuncPC->SetLineColor(graphPC->GetMarkerColor()); myfuncPC->Draw("same");
      TF1* myfuncPF= new TF1("FitgraphPF","pol1",250,5000);  graphPF->Fit(myfuncPF,"QN","",250,5000); myfuncPF->SetLineWidth(2); myfuncPF->SetLineColor(graphPF->GetMarkerColor()); myfuncPF->Draw("same");
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Pixel) |eta|<0.5", myfuncPC->GetChisquare()/ myfuncPC->GetNDF(), myfuncPC->GetParameter(0),myfuncPC->GetParError(0),myfuncPC->GetParameter(1),myfuncPC->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Pixel) |eta|>1.5", myfuncPF->GetChisquare()/ myfuncPF->GetNDF(), myfuncPF->GetParameter(0),myfuncPF->GetParError(0),myfuncPF->GetParameter(1),myfuncPF->GetParError(1));
      leg = new TLegend(0.79,0.92,0.79-0.20,0.92 - 3*0.05);     leg->SetFillColor(0);     leg->SetBorderSize(0);
      leg->AddEntry(graphPC, "dE/dx (Pixel) |#eta|<0.5" ,"P");
      leg->AddEntry(graphPF, "dE/dx (Pixel) |#eta|>1.5" ,"P");
      leg->Draw();
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_dEdxMP");
      delete c1; delete leg;




      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphNV = ConvertFromRunToIntLumi(SingleMu_NVertProf, "A*" , "<#Reco Vertices>",0,0);
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_Vert");
      delete c1;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphpT = ConvertFromRunToIntLumi(SingleMu_PtProf, "A*" , "<p_{T}> (GeV/c)",0,0);
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_pT");
      delete c1;


   }else{
      printf("TEST TEST TEST\n");
   }


   for(int i=0;i<SingleMu_PtProf->GetXaxis()->GetNbins();i++){
//      if((i+3)%4==0)continue;
      SingleMu_PtProf->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxProf->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMProf->GetXaxis()->SetBinLabel(i,"");
      SingleMu_NVertProf->GetXaxis()->SetBinLabel(i,"");
   }  


   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_NVertProf;                 legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "<#Reco Vertices>", 0,0, 0,0);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_NVert");
   delete c1;

 

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_PtProf;                    legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "p_{T} (GeV/c)", 0,0, 0,0);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_Pt");
   delete c1;


   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_dEdxProf;                  legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "I_{as}", 0,0, 0.02,0.06);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_dEdx");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_dEdxMProf;                  legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "I_{h}", 0,0, 3.2,3.4);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_dEdxM");
   delete c1;

/*
   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_dEdxMProfLumi166380;       legend.push_back("SingleMu50 - Run166380");
   Histos[1] = SingleMu_dEdxMProfLumi166512;       legend.push_back("SingleMu50 - Run166512");
   Histos[2] = SingleMu_dEdxMProfLumi167807;       legend.push_back("SingleMu50 - Run167807");
   Histos[3] = SingleMu_dEdxMProfLumi167898;       legend.push_back("SingleMu50 - Run167898");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Lumi", "I_{h}", 0,0, 3.2,3.4);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_dEdxMRun");
   delete c1;
*/


}

void GetMeanAndRMS(TFile* InputFile, string path, double& mean, double& rms, bool gaussianFit=false){
   TH1D* histo = (TH1D*)GetObjectFromPath(InputFile, path);
   if(!histo){
      mean = -999;
      rms  = -999;
   }else{
      mean = histo->GetMean();
      rms  = histo->GetMeanError();

      if(histo->GetEntries()>10 && gaussianFit){
         rms  = histo->GetRMS(); 
         TF1* fit = new TF1("MyFit","gaus", mean-2*rms, mean+2*rms);
         fit->SetParameters(histo->GetEntries(),  histo->GetMean(), histo->GetRMS());
         histo->Fit(fit, "0QR"); 
         mean = fit->GetParameter(1);
         rms  = fit->GetParError(1);
         delete fit;
      }
      delete histo;
   }
}

TGraphErrors* getStabilityGraph(std::vector<string>& runList, TFile* InputFile, string HistoName, bool gaussianFit=false){
   TGraphErrors* graph = new TGraphErrors(runList.size()); 
   graph->SetName(HistoName.c_str());

   double mean, rms;
   for(unsigned int r=0;r<runList.size();r++){
      char name[1024]; sprintf(name, "%s/%s", runList[r].c_str(), HistoName.c_str());
      GetMeanAndRMS(InputFile, name, mean, rms); graph->SetPoint(r, r+0.5, mean);   graph->SetPointError(r, 0, rms);
   }
   graph->SetLineWidth(2);
   graph->SetLineColor(1);
   graph->SetFillColor(1);
   graph->SetMarkerSize(0.5);
   graph->SetMarkerStyle(20);
   graph->SetMarkerColor(1);
   return graph;
}

void MakePlot()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.07);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(55);  //Rainbow --> deifne the color of all plots
   gStyle->SetNdivisions(505);

   system("mkdir -p pictures");


   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;

   TFile* InputFile = new TFile("pictures/Histos.root");
   std::vector<string> runList;
   TList* ObjList = InputFile->GetListOfKeys();
   for(int i=0;i<ObjList->GetSize();i++){
      TObject* tmp = GetObjectFromPath(InputFile,ObjList->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TDirectory")){
         runList.push_back(ObjList->At(i)->GetName());
         printf("Add a new run: %s\n", ObjList->At(i)->GetName() );
      }
      delete tmp;
   }
   unsigned int N= runList.size();

   TH1D* frame = new TH1D("frame", "frame", N, 0, N);
   frame->SetTitle("");
   frame->SetStats(kFALSE);
   frame->GetXaxis()->SetNdivisions(505);
   frame->GetXaxis()->SetTitle("");
   frame->GetYaxis()->SetTitleOffset(0.95);
   for(unsigned int r=0;r<N;r++){frame->GetXaxis()->SetBinLabel(r+1, ((r+0)%2==0)?runList[r].c_str():"");}  //plot only a label every 2




   TGraphErrors* Any_Pt                = getStabilityGraph(runList, InputFile, "AnyPt");
   TGraphErrors* SingleMu_Pt           = getStabilityGraph(runList, InputFile, "HLT_Mu50Pt");
   TGraphErrors* PFMet_Pt              = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedPt");

   TGraphErrors* Any_dEdx              = getStabilityGraph(runList, InputFile, "AnydEdx");
   TGraphErrors* SingleMu_dEdx         = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdx");
   TGraphErrors* PFMet_dEdx            = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdx");

   TGraphErrors* Any_dEdxM             = getStabilityGraph(runList, InputFile, "AnydEdxM");
   TGraphErrors* SingleMu_dEdxM        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxM");
   TGraphErrors* PFMet_dEdxM           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxM");

   TGraphErrors* Any_dEdxMS             = getStabilityGraph(runList, InputFile, "AnydEdxMS");
   TGraphErrors* SingleMu_dEdxMS        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxMS");
   TGraphErrors* PFMet_dEdxMS           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxMS");

   TGraphErrors* Any_dEdxMP             = getStabilityGraph(runList, InputFile, "AnydEdxMP");
   TGraphErrors* SingleMu_dEdxMP        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxMP");
   TGraphErrors* PFMet_dEdxMP           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxMP");

   TGraphErrors* Any_dEdxMSC             = getStabilityGraph(runList, InputFile, "AnydEdxMSC");
   TGraphErrors* SingleMu_dEdxMSC        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxMSC");
   TGraphErrors* PFMet_dEdxMSC           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxMSC");

   TGraphErrors* Any_dEdxMPC             = getStabilityGraph(runList, InputFile, "AnydEdxMPC");
   TGraphErrors* SingleMu_dEdxMPC        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxMPC");
   TGraphErrors* PFMet_dEdxMPC           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxMPC");

   TGraphErrors* Any_dEdxMSF             = getStabilityGraph(runList, InputFile, "AnydEdxMSF");
   TGraphErrors* SingleMu_dEdxMSF        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxMSF");
   TGraphErrors* PFMet_dEdxMSF           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxMSF");

   TGraphErrors* Any_dEdxMPF             = getStabilityGraph(runList, InputFile, "AnydEdxMPF");
   TGraphErrors* SingleMu_dEdxMPF        = getStabilityGraph(runList, InputFile, "HLT_Mu50dEdxMPF");
   TGraphErrors* PFMet_dEdxMPF           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleaneddEdxMPF");


   TGraphErrors* Any_TOFAOD               = getStabilityGraph(runList, InputFile, "AnyTOFAOD", true);
   TGraphErrors* SingleMu_TOFAOD          = getStabilityGraph(runList, InputFile, "HLT_Mu50TOFAOD", true);
   TGraphErrors* PFMet_TOFAOD             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedTOFAOD", true);

   TGraphErrors* Any_TOFAODDT             = getStabilityGraph(runList, InputFile, "AnyTOFAODDT", true);
   TGraphErrors* SingleMu_TOFAODDT        = getStabilityGraph(runList, InputFile, "HLT_Mu50TOFAODDT", true);
   TGraphErrors* PFMet_TOFAODDT           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedTOFAODDT", true);

   TGraphErrors* Any_TOFAODCSC            = getStabilityGraph(runList, InputFile, "AnyTOFAODCSC", true);
   TGraphErrors* SingleMu_TOFAODCSC       = getStabilityGraph(runList, InputFile, "HLT_Mu50TOFAODCSC", true);
   TGraphErrors* PFMet_TOFAODCSC          = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedTOFAODCSC", true);

   TGraphErrors* Any_VertexAOD               = getStabilityGraph(runList, InputFile, "AnyVertexAOD", true);
   TGraphErrors* SingleMu_VertexAOD          = getStabilityGraph(runList, InputFile, "HLT_Mu50VertexAOD", true);
   TGraphErrors* PFMet_VertexAOD             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedVertexAOD", true);

   TGraphErrors* Any_VertexAODDT             = getStabilityGraph(runList, InputFile, "AnyVertexAODDT", true);
   TGraphErrors* SingleMu_VertexAODDT        = getStabilityGraph(runList, InputFile, "HLT_Mu50VertexAODDT", true);
   TGraphErrors* PFMet_VertexAODDT           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedVertexAODDT", true);

   TGraphErrors* Any_VertexAODCSC            = getStabilityGraph(runList, InputFile, "AnyVertexAODCSC", true);
   TGraphErrors* SingleMu_VertexAODCSC       = getStabilityGraph(runList, InputFile, "HLT_Mu50VertexAODCSC", true);
   TGraphErrors* PFMet_VertexAODCSC          = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedVertexAODCSC", true);

   TGraphErrors* Any_TOF               = getStabilityGraph(runList, InputFile, "AnyTOF", true);
   TGraphErrors* SingleMu_TOF          = getStabilityGraph(runList, InputFile, "HLT_Mu50TOF", true);
   TGraphErrors* PFMet_TOF             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedTOF", true);

   TGraphErrors* Any_TOFDT             = getStabilityGraph(runList, InputFile, "AnyTOFDT", true);
   TGraphErrors* SingleMu_TOFDT        = getStabilityGraph(runList, InputFile, "HLT_Mu50TOFDT", true);
   TGraphErrors* PFMet_TOFDT           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedTOFDT", true);

   TGraphErrors* Any_TOFCSC            = getStabilityGraph(runList, InputFile, "AnyTOFCSC", true);
   TGraphErrors* SingleMu_TOFCSC       = getStabilityGraph(runList, InputFile, "HLT_Mu50TOFCSC", true);
   TGraphErrors* PFMet_TOFCSC          = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedTOFCSC", true);

   TGraphErrors* Any_Vertex               = getStabilityGraph(runList, InputFile, "AnyVertex", true);
   TGraphErrors* SingleMu_Vertex          = getStabilityGraph(runList, InputFile, "HLT_Mu50Vertex", true);
   TGraphErrors* PFMet_Vertex             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedVertex", true);

   TGraphErrors* Any_VertexDT             = getStabilityGraph(runList, InputFile, "AnyVertexDT", true);
   TGraphErrors* SingleMu_VertexDT        = getStabilityGraph(runList, InputFile, "HLT_Mu50VertexDT", true);
   TGraphErrors* PFMet_VertexDT           = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedVertexDT", true);

   TGraphErrors* Any_VertexCSC            = getStabilityGraph(runList, InputFile, "AnyVertexCSC", true);
   TGraphErrors* SingleMu_VertexCSC       = getStabilityGraph(runList, InputFile, "HLT_Mu50VertexCSC", true);
   TGraphErrors* PFMet_VertexCSC          = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedVertexCSC", true);

//   TGraphErrors* Any_HdEdx                 = getStabilityGraph(runList, InputFile, "AnyHdEdx");
//   TGraphErrors* SingleMu_HdEdx          = getStabilityGraph(runList, InputFile, "HLT_Mu50HdEdx");
//   TGraphErrors* PFMet_HdEdx             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedHdEdx");

//   TGraphErrors* Any_HPt                 = getStabilityGraph(runList, InputFile, "AnyHPt");
//   TGraphErrors* SingleMu_HPt          = getStabilityGraph(runList, InputFile, "HLT_Mu50HPt");
//   TGraphErrors* PFMet_HPt             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedHPt");

//   TGraphErrors* Any_HTOFAOD                 = getStabilityGraph(runList, InputFile, "AnyHTOFAOD");
//   TGraphErrors* SingleMu_HTOFAOD          = getStabilityGraph(runList, InputFile, "HLT_Mu50HTOFAOD");
//   TGraphErrors* PFMet_HTOFAOD             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedHTOFAOD");

//   TGraphErrors* Any_HTOF                 = getStabilityGraph(runList, InputFile, "AnyHTOF");
//   TGraphErrors* SingleMu_HTOF          = getStabilityGraph(runList, InputFile, "HLT_Mu50HTOF");
//   TGraphErrors* PFMet_HTOF             = getStabilityGraph(runList, InputFile, "HLT_PFMET170_NoiseCleanedHTOF");


/*
   for(int i=0;i<SingleMu_Pt->GetXaxis()->GetNbins();i++){
//      if((i+3)%12==0)continue; 
      if((i+2)%4==0)continue; 
//      continue;// //uncomment to show less label that bins
      Any_Pt         ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_Pt    ->GetXaxis()->SetBinLabel(i,"");
      PFMet_Pt       ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdx       ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdx  ->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdx     ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxM      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxM ->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxM    ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxMS      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMS ->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxMS    ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxMP      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMP ->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxMP    ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxMSC     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMSC->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxMSC   ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxMPC     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMPC->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxMPC   ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxMSF     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMSF->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxMSF   ->GetXaxis()->SetBinLabel(i,"");

      Any_dEdxMPF     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMPF->GetXaxis()->SetBinLabel(i,"");
      PFMet_dEdxMPF   ->GetXaxis()->SetBinLabel(i,"");

      Any_TOFAOD        ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_TOFAOD   ->GetXaxis()->SetBinLabel(i,"");
      PFMet_TOFAOD      ->GetXaxis()->SetBinLabel(i,"");

      Any_TOFAODDT      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_TOFAODDT ->GetXaxis()->SetBinLabel(i,"");
      PFMet_TOFAODDT    ->GetXaxis()->SetBinLabel(i,"");

      Any_TOFAODCSC     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_TOFAODCSC->GetXaxis()->SetBinLabel(i,"");
      PFMet_TOFAODCSC   ->GetXaxis()->SetBinLabel(i,"");

      Any_VertexAOD        ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_VertexAOD   ->GetXaxis()->SetBinLabel(i,"");
      PFMet_VertexAOD      ->GetXaxis()->SetBinLabel(i,"");

      Any_VertexAODDT      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_VertexAODDT ->GetXaxis()->SetBinLabel(i,"");
      PFMet_VertexAODDT    ->GetXaxis()->SetBinLabel(i,"");

      Any_VertexAODCSC     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_VertexAODCSC->GetXaxis()->SetBinLabel(i,"");
      PFMet_VertexAODCSC   ->GetXaxis()->SetBinLabel(i,"");

      Any_TOF        ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_TOF   ->GetXaxis()->SetBinLabel(i,"");
      PFMet_TOF      ->GetXaxis()->SetBinLabel(i,"");

      Any_TOFDT      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_TOFDT ->GetXaxis()->SetBinLabel(i,"");
      PFMet_TOFDT    ->GetXaxis()->SetBinLabel(i,"");

      Any_TOFCSC     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_TOFCSC->GetXaxis()->SetBinLabel(i,"");
      PFMet_TOFCSC   ->GetXaxis()->SetBinLabel(i,"");

      Any_Vertex        ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_Vertex   ->GetXaxis()->SetBinLabel(i,"");
      PFMet_Vertex      ->GetXaxis()->SetBinLabel(i,"");

      Any_VertexDT      ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_VertexDT ->GetXaxis()->SetBinLabel(i,"");
      PFMet_VertexDT    ->GetXaxis()->SetBinLabel(i,"");

      Any_VertexCSC     ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_VertexCSC->GetXaxis()->SetBinLabel(i,"");
      PFMet_VertexCSC   ->GetXaxis()->SetBinLabel(i,"");

      Any_HdEdx          ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_HdEdx     ->GetXaxis()->SetBinLabel(i,"");
      PFMet_HdEdx        ->GetXaxis()->SetBinLabel(i,"");

      Any_HPt            ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_HPt       ->GetXaxis()->SetBinLabel(i,"");
      PFMet_HPt          ->GetXaxis()->SetBinLabel(i,"");

      Any_HTOFAOD        ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_HTOFAOD   ->GetXaxis()->SetBinLabel(i,"");
      PFMet_HTOFAOD      ->GetXaxis()->SetBinLabel(i,"");

      Any_HTOF           ->GetXaxis()->SetBinLabel(i,"");
      SingleMu_HTOF      ->GetXaxis()->SetBinLabel(i,"");
      PFMet_HTOF         ->GetXaxis()->SetBinLabel(i,"");
   }  
*/

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("p_{T} (GeV/c)");   frame->SetMinimum(0.0);   frame->SetMaximum(150.0);  frame->Draw("AXIS");
   SingleMu_Pt->Draw("0 P same");;                  legend.push_back("SingleMu50");
   //for(unsigned int i=0;i<legend.size();i++){((Tile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_Pt");
   delete c1;


   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{as}");   frame->SetMinimum(0.0);   frame->SetMaximum(0.05);  frame->Draw("AXIS");
   //Any_dEdx->Draw("0 P same");                       legend.push_back("Any");
   SingleMu_dEdx->Draw("0 P same");                  legend.push_back("SingleMu50");
   //PFMet_dEdx->Draw("0 P same");                     legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdx");
   delete c1;


   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h}");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxM->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxM->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxM->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxM");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h} S");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxMS->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxMS->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxMS->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxMS");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h} P");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxMP->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxMP->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxMP->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxMP");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h} SC");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxMSC->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxMSC->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxMSC->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxMSC");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h} PC");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxMPC->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxMPC->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxMPC->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxMPC");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h} SF");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxMSF->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxMSF->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxMSF->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxMSF");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("I_{h} PF");   frame->SetMinimum(2.5);   frame->SetMaximum(3.7);  frame->Draw("AXIS");
   //Any_dEdxMPF->Draw("0 P same");                      legend.push_back("Any");
   SingleMu_dEdxMPF->Draw("0 P same");                 legend.push_back("SingleMu50");
   //PFMet_dEdxMPF->Draw("0 P same");                    legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_dEdxMPF");
   delete c1;


   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF} (AOD)");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   //Any_TOFAOD->Draw("0 P same");                        legend.push_back("Any");   
   SingleMu_TOFAOD->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_TOFAOD->Draw("0 P same");                      legend.push_back("PFMET170");
   
   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOFAOD");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF_DT} (AOD)");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   //Any_TOFAODDT->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_TOFAODDT->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_TOFAODDT->Draw("0 P same");                      legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOFAODDT");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF_CSC} (AOD)");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   //Any_TOFAODCSC->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_TOFAODCSC->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_TOFAODCSC->Draw("0 P same");                      legend.push_back("PFMET170");
      
   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOFAODCSC");
   delete c1;



   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF} (AOD)");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   SingleMu_TOFAOD->Draw("0 P same");                   legend.push_back("DT+CSC");
   SingleMu_TOFAODDT->Draw("0 P same");                 legend.push_back("DT");
   SingleMu_TOFAODCSC->Draw("0 P same");                legend.push_back("CSC");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P", 0.29, 0.52, 0.20, 0.05);
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOFAOD_Combined");
   delete c1;



   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("Vertex time [ns] (AOD)");   frame->SetMinimum(-2.0);   frame->SetMaximum(2.0);  frame->Draw("AXIS");
   //Any_VertexAOD->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_VertexAOD->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_VertexAOD->Draw("0 P same");                      legend.push_back("PFMET170");
   
   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_VertexAOD");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("Vertex time DT [ns] (AOD)");   frame->SetMinimum(-2.0);   frame->SetMaximum(2.0);  frame->Draw("AXIS");
   //Any_VertexAODDT->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_VertexAODDT->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_VertexAODDT->Draw("0 P same");                      legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_VertexAODDT");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("Vertex time CSC [ns] (AOD)");   frame->SetMinimum(-2.0);   frame->SetMaximum(2.0);  frame->Draw("AXIS");
   //Any_VertexAODCSC->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_VertexAODCSC->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_VertexAODCSC->Draw("0 P same");                      legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_VertexAODCSC");
   delete c1;


   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF}");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   //Any_TOF->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_TOF->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_TOF->Draw("0 P same");                      legend.push_back("PFMET170");
   
   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOF");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF_DT}");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   //Any_TOFDT->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_TOFDT->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_TOFDT->Draw("0 P same");                      legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOFDT");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF_CSC}");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   //Any_TOFCSC->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_TOFCSC->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_TOFCSC->Draw("0 P same");                      legend.push_back("PFMET170");
      
   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOFCSC");
   delete c1;



   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("1/#beta_{TOF}");   frame->SetMinimum(0.85);   frame->SetMaximum(1.15);  frame->Draw("AXIS");
   SingleMu_TOF->Draw("0 P same");                   legend.push_back("DT+CSC");
   SingleMu_TOFDT->Draw("0 P same");                 legend.push_back("DT");
   SingleMu_TOFCSC->Draw("0 P same");                legend.push_back("CSC");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P", 0.29, 0.52, 0.20, 0.05);
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_TOF_Combined");
   delete c1;



   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("Vertex time [ns]");   frame->SetMinimum(-5.0);   frame->SetMaximum(5.0);  frame->Draw("AXIS");
   //Any_Vertex->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_Vertex->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_Vertex->Draw("0 P same");                      legend.push_back("PFMET170");
   
   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_Vertex");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("Vertex time DT [ns]");   frame->SetMinimum(-5.0);   frame->SetMaximum(5.0);  frame->Draw("AXIS");
   //Any_VertexDT->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_VertexDT->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_VertexDT->Draw("0 P same");                      legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_VertexDT");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   frame->GetYaxis()->SetTitle("Vertex time CSC [ns]");   frame->SetMinimum(-5.0);   frame->SetMaximum(5.0);  frame->Draw("AXIS");
   //Any_VertexCSC->Draw("0 P same");                        legend.push_back("Any");
   SingleMu_VertexCSC->Draw("0 P same");                   legend.push_back("SingleMu50");
   //PFMet_VertexCSC->Draw("0 P same");                      legend.push_back("PFMET170");

   //for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   //DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","Summary_Profile_VertexCSC");
   delete c1;

//   MakedEdxPlot();
}


