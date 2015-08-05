
#include <exception>
#include <vector>
#include <fstream>

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
#include "TGraph.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCutG.h"
#include "../../AnalysisCode/tdrstyle.C"
#include "../../AnalysisCode/Analysis_PlotFunction.h"

using namespace std;


void getScaleFactor(TFile* InputFile1, string OutName, string ObjName1, TFile* InputFile2, string ObjName2="EMPTY");
void ExtractConstants(TH2D* input, int FileIndex=0);
void DrawComparisons (TFile* InputFile1, TFile* InputFile2=NULL, string ObjName1="Ias_SO", string ObjName2="Ias_SO_inc");

//const double K = 2.4496; //Truncated40
//const double C = 2.2364; //Truncated40

//const double K = 2.4236; //Truncated40
//const double C = 2.6474; //Truncated40

//const double K = 2.24; //Truncated40
//const double C = 2.72; //Truncated40

//const double K = 2.37; //Truncated40
//const double C = 2.55; //Truncated40

//const double K = 2.63; //Truncated40
//const double C = 2.39; //Truncated40

//2012 constants
//const double K = 2.529; //Harm2
//const double C = 2.772; //Harm2

//2015B prompt constants -- one for Each File
double K [2] = {2.779, 2.779}; double Kerr [2] = {0.001, 0.001};   //Harm2
double C [2] = {2.879, 2.779}; double Cerr [2] = {0.001, 0.001};   //Harm2

double GetMass(double P, double I, int FileIndex=0){
   return sqrt((I-C[FileIndex])/K[FileIndex])*P;
}

TF1* GetMassLine(double M, bool left=false, int FileIndex=0)
{  
   double BetaMax = 0.9;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));

   double BetaMin = 0.2;
   double PMin = sqrt((BetaMin*BetaMin*M*M)/(1-BetaMin*BetaMin));
   
   if(left){PMax*=-1; PMin*=-1;}

   TF1* MassLine = new TF1("MassLine","[2] + ([0]*[0]*[1])/(x*x)", PMin, PMax);
   MassLine->SetParName  (0,"M");
   MassLine->SetParName  (1,"K");
   MassLine->SetParName  (2,"C");
   MassLine->SetParameter(0, M);
   MassLine->SetParameter(1, K[FileIndex]);
   MassLine->SetParameter(2, C[FileIndex]);
   MassLine->SetLineWidth(2);
   return MassLine;
}

void MakePlot(string INPUT, string INPUT2="EMPTY")
{
   size_t firstToken = INPUT.find_first_of ("_"),
          lastToken  = INPUT.find_first_of (".") - firstToken;
   string SaveName (INPUT, firstToken, lastToken),
          SaveName2 ("");
   if (INPUT2 != "EMPTY"){
      firstToken = INPUT2.find_first_of ("_");
      lastToken  = INPUT2.find_first_of (".") - firstToken;
      SaveName2  = INPUT2.substr (firstToken, lastToken);
   }

   string SaveDir ("pictures" + SaveName + SaveName2 + "/");
   system(("mkdir -p " + SaveDir).c_str());
   system("mkdir -p fit");

   string SaveSuffix = INPUT2=="EMPTY" ? "" : SaveName;

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.125);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
//   gStyle->SetPalette(51); 
   gStyle->SetPalette(1); 
   gStyle->SetNdivisions(510,"X");

   TF1* PionLine = GetMassLine(0.140);
   PionLine->SetLineColor(1);
   PionLine->SetLineWidth(2);

   TF1* KaonLine = GetMassLine(0.494);
   KaonLine->SetLineColor(1);
   KaonLine->SetLineWidth(2);

   TF1* ProtonLine = GetMassLine(0.938);
   ProtonLine->SetLineColor(1);
   ProtonLine->SetLineWidth(2);

   TF1* DeuteronLine = GetMassLine(1.88);
   DeuteronLine->SetLineColor(1);
   DeuteronLine->SetLineWidth(2);

   TF1* TritonLine = GetMassLine(2.80);
   TritonLine->SetLineColor(1);
   TritonLine->SetLineWidth(2);

   TF1* ProtonLineFit = GetMassLine(0.938);
   ProtonLineFit->SetLineColor(2);
   ProtonLineFit->SetLineWidth(2);
   ProtonLineFit->SetRange(0.6,1.2);

   TF1* KaonLineLeft = GetMassLine(0.494, true);
   KaonLineLeft->SetLineColor(1);
   KaonLineLeft->SetLineWidth(2);

   TF1* ProtonLineLeft = GetMassLine(0.938, true);
   ProtonLineLeft->SetLineColor(1);
   ProtonLineLeft->SetLineWidth(2);

   TF1* DeuteronLineLeft = GetMassLine(1.88, true);
   DeuteronLineLeft->SetLineColor(1);
   DeuteronLineLeft->SetLineWidth(2);

   TF1* TritonLineLeft = GetMassLine(2.80, true);
   TritonLineLeft->SetLineColor(1);
   TritonLineLeft->SetLineWidth(2);


   TFile* InputFile  = new TFile(INPUT .c_str());
   TFile* InputFile2 = INPUT2.find("EMPTY")!=string::npos ? NULL : new TFile(INPUT2.c_str());


   std::vector<string> ObjName;
   ObjName.push_back("hit_SP");
   ObjName.push_back("hit_SP_in");
   ObjName.push_back("harm2_SO");
   ObjName.push_back("harm2_SP");
   ObjName.push_back("harm2_SO_in");
   ObjName.push_back("harm2_PO_raw"); // FIXME does not fit well
   ObjName.push_back("harm2_SO_raw"); // FIXME does not fit well
   ObjName.push_back("harm2_SP_raw"); // FIXME does not fit well
   ObjName.push_back("Ias_SO_inc");
   ObjName.push_back("Ias_SO");

//   DrawComparisons (InputFile, InputFile2);
//   DrawComparisons (InputFile, InputFile2, "harm2_SO_raw", "harm2_SO");
//   DrawComparisons (InputFile, InputFile2, "Ias_SO", "Ias_SO_inc");

   ofstream ExtractConstantsReport, ExtractConstantsReport2;
   ExtractConstantsReport.open ((SaveDir + "ConstantsReport" + SaveName + ".txt").c_str(), ofstream::out);
   if (InputFile2) ExtractConstantsReport2.open ((SaveDir + "ConstantsReport" + SaveName2 + ".txt").c_str(), ofstream::out);

   for(unsigned int i=0;i<ObjName.size();i++){
      TH3F*       dEdxTemplate       = (TH3F*)      GetObjectFromPath(InputFile, (ObjName[i] + "_ChargeVsPath"      ).c_str() );
      TH1D*       HdedxMIP           = (TH1D*)      GetObjectFromPath(InputFile, (ObjName[i] + "_MIP"               ).c_str() );
      TH1D*       HdedxSIG           = (TH1D*)      GetObjectFromPath(InputFile, (ObjName[i] + "_SIG"               ).c_str() );
      TH1D*       HMass              = (TH1D*)      GetObjectFromPath(InputFile, (ObjName[i] + "_Mass"              ).c_str() );
      TH2D*       HdedxVsP           = (TH2D*)      GetObjectFromPath(InputFile, (ObjName[i] + "_dedxVsP"           ).c_str() );
      TH2D*       HdedxVsQP          = (TH2D*)      GetObjectFromPath(InputFile, (ObjName[i] + "_dedxVsQP"          ).c_str() );
      TProfile*   HdedxVsPProfile    = (TProfile*)  GetObjectFromPath(InputFile, (ObjName[i] + "_Profile"           ).c_str() );
      TH2D*       HdedxVsEta         = (TH2D*)      GetObjectFromPath(InputFile, (ObjName[i] + "_Eta2D"             ).c_str() );
      TProfile*   HdedxVsEtaProfile  = (TProfile*)  GetObjectFromPath(InputFile, (ObjName[i] + "_Eta"               ).c_str() );
      TProfile2D* HdedxVsP_NS        = (TProfile2D*)GetObjectFromPath(InputFile, (ObjName[i] + "_dedxVsP_NS"        ).c_str() );
      TH3F*       dEdxTemplate2      = NULL;
      TH1D*       HdedxMIP2          = NULL;
      TH1D*       HdedxSIG2          = NULL;
      TH1D*       HMass2             = NULL;
      TH2D*       HdedxVsP2          = NULL;
      TH2D*       HdedxVsQP2         = NULL;
      TProfile*   HdedxVsPProfile2   = NULL;
      TH2D*       HdedxVsEta2        = NULL;
      TProfile*   HdedxVsEtaProfile2 = NULL;
      TProfile2D* HdedxVsP_NS2       = NULL;

      if (InputFile2) {
         dEdxTemplate2      = (TH3F*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_ChargeVsPath"  ).c_str() );
         HdedxMIP2          = (TH1D*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_MIP"           ).c_str() );
         HdedxSIG2          = (TH1D*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_SIG"           ).c_str() );
         HMass2             = (TH1D*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_Mass"          ).c_str() );
         HdedxVsP2          = (TH2D*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_dedxVsP"       ).c_str() );
         HdedxVsQP2         = (TH2D*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_dedxVsQP"      ).c_str() );
         HdedxVsPProfile2   = (TProfile*)  GetObjectFromPath(InputFile2, (ObjName[i] + "_Profile"       ).c_str() );
         HdedxVsEta2        = (TH2D*)      GetObjectFromPath(InputFile2, (ObjName[i] + "_Eta2D"         ).c_str() );
         HdedxVsEtaProfile2 = (TProfile*)  GetObjectFromPath(InputFile2, (ObjName[i] + "_Eta"           ).c_str() );
         HdedxVsP_NS2       = (TProfile2D*)GetObjectFromPath(InputFile2, (ObjName[i] + "_dedxVsP_NS"    ).c_str() );
      }

      if (ObjName[i].find("hit_SP")!=string::npos){
         dEdxTemplate->SetName("Charge_Vs_Path");
         dEdxTemplate->SaveAs (("dEdxTemplate_" + ObjName[i] + SaveName + ".root").c_str());

         if (InputFile2){
            dEdxTemplate2->SetName("Charge_Vs_Path");
            dEdxTemplate2->SaveAs (("dEdxTemplate_" + ObjName[i] + SaveName2 + ".root").c_str());
         }

         continue;
      }

      ExtractConstants(HdedxVsP);
      ExtractConstantsReport << ObjName[i] << " ... K = " << K[0] << " +- " << Kerr[0] << "\t"
                                           << " ... C = " << C[0] << " +- " << Cerr[0] << endl;

      TPaveText* T = new TPaveText(0.05, 0.995, 0.95, 0.945, "NDC");
      T->SetTextFont(43);  //give the font size in pixel (instead of fraction)
      T->SetTextSize(21);  //font size
      T->SetBorderSize(0);
      T->SetFillColor(0);
      T->SetFillStyle(0);
      T->SetTextAlign(22);
      T->AddText("#bf{CMS} Preliminary   -   2.74 pb^{-1}   -    #sqrt{s} = 13 TeV");

      std::cout << "TESTA\n";
      TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
      c1->SetLogz(true);
      HdedxVsP->SetStats(kFALSE);
      HdedxVsP->GetXaxis()->SetTitle("momentum (GeV/c)");
      HdedxVsP->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxVsP->SetAxisRange(0,5,"X");
      HdedxVsP->SetAxisRange(0,15,"Y");
      HdedxVsP->Draw("COLZ");

      PionLine->Draw("same");
      KaonLine->Draw("same");
      ProtonLine->Draw("same");
      DeuteronLine->Draw("same");
//      TritonLine->Draw("same");
      ProtonLineFit->Draw("same");
      T->Draw("same");
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_dedxVsP", true);
//      c1->SaveAs((SaveDir + ObjName[i] + SaveSuffix + "_dedxVsP.C").c_str());
//      c1->SaveAs((SaveDir + ObjName[i] + SaveSuffix + "_dedxVsP.pdf").c_str());
      delete c1;


      c1 = new TCanvas("c1", "c1", 600,600);
      c1->SetLogz(true);
      HdedxVsQP->SetStats(kFALSE);
      HdedxVsQP->GetXaxis()->SetTitle("charge * momentum (GeV/c)");
      HdedxVsQP->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxVsQP->SetAxisRange(-5,5,"X");
      HdedxVsQP->SetAxisRange(0,15,"Y");
      HdedxVsQP->Draw("COLZ");

//      KaonLine->Draw("same");
//      ProtonLine->Draw("same");
//      DeuteronLine->Draw("same");
//      TritonLine->Draw("same");
//      ProtonLineFit->Draw("same");

//      KaonLineLeft->Draw("same");
//      ProtonLineLeft->Draw("same");
//      DeuteronLineLeft->Draw("same");
//      TritonLineLeft->Draw("same");

      
      T->Draw("same");     
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_dedxVsQP", true);
//      c1->SaveAs((SaveDir + ObjName[i] + SaveSuffix + "_dedxVsQP.C"  ).c_str());
//      c1->SaveAs((SaveDir + ObjName[i] + SaveSuffix + "_dedxVsQP.pdf").c_str());
      delete c1;

      c1 = new TCanvas("c1", "c1", 600,600);
      c1->SetLogz(true);
      HdedxVsP_NS->SetStats(kFALSE);
      HdedxVsP_NS->GetXaxis()->SetTitle("track momentum (GeV/c)");
      HdedxVsP_NS->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxVsP_NS->GetZaxis()->SetTitle("Number of saturated strips");
      HdedxVsP_NS->SetAxisRange(0.0, 5.0, "X");
      HdedxVsP_NS->Draw("COLZ");
      PionLine->Draw("same");
      KaonLine->Draw("same");
      ProtonLine->Draw("same");
      DeuteronLine->Draw("same");
//      TritonLine->Draw("same");
      ProtonLineFit->Draw("same");
      T->Draw("same");
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_dedxVsP_NS");
      delete c1;

      c1 = new TCanvas("c1", "c1", 600,600);
      c1->SetLogz(true);
      HdedxVsEta->SetStats(kFALSE);
      HdedxVsEta->GetXaxis()->SetTitle("Eta");
      HdedxVsEta->GetYaxis()->SetTitle("I_{as}");
      HdedxVsEta->SetAxisRange(-2.1,2.1,"X");
      HdedxVsEta->Draw("COLZ");
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_Eta2D");
      delete c1;

      if (InputFile2) {
         ExtractConstants(HdedxVsP, 1);
         ExtractConstantsReport2 << ObjName[i] << " ... K = " << K[1] << " +- " << Kerr[1] << "\t"
                                               << " ... C = " << C[1] << " +- " << Cerr[1] << endl;

         TF1* PionLine2 = GetMassLine(0.140, false, 1);
         PionLine2->SetLineColor(1);
         PionLine2->SetLineWidth(2);

         TF1* KaonLine2 = GetMassLine(0.494, false, 1);
         KaonLine2->SetLineColor(1);
         KaonLine2->SetLineWidth(2);

         TF1* ProtonLine2 = GetMassLine(0.938, false, 1);
         ProtonLine2->SetLineColor(1);
         ProtonLine2->SetLineWidth(2);

         TF1* DeuteronLine2 = GetMassLine(1.88, false, 1);
         DeuteronLine2->SetLineColor(1);
         DeuteronLine2->SetLineWidth(2);

         TF1* TritonLine2 = GetMassLine(2.80, false, 1);
         TritonLine2->SetLineColor(1);
         TritonLine2->SetLineWidth(2);

         TF1* ProtonLineFit2 = GetMassLine(0.938, false, 1);
         ProtonLineFit2->SetLineColor(2);
         ProtonLineFit2->SetLineWidth(2);
         ProtonLineFit2->SetRange(0.6,1.2);

         TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         HdedxVsP2->SetStats(kFALSE);
         HdedxVsP2->GetXaxis()->SetTitle("momentum (GeV/c)");
         HdedxVsP2->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
         HdedxVsP2->SetAxisRange(0,5,"X");
         HdedxVsP2->SetAxisRange(0,15,"Y");
         HdedxVsP2->Draw("COLZ");
    
         PionLine2->Draw("same");
         KaonLine2->Draw("same");
         ProtonLine2->Draw("same");
         DeuteronLine2->Draw("same");
//         TritonLine->Draw("same");
         ProtonLineFit2->Draw("same");
         T->Draw("same");
         SaveCanvas(c1, SaveDir, ObjName[i] + SaveName2 + "_dedxVsP", true);
//         c1->SaveAs((SaveDir + ObjName[i] + SaveName2 + "_dedxVsP.C").c_str());
//         c1->SaveAs((SaveDir + ObjName[i] + SaveName2 + "_dedxVsP.pdf").c_str());
         delete c1;


         c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         HdedxVsQP->SetStats(kFALSE);
         HdedxVsQP->GetXaxis()->SetTitle("charge * momentum (GeV/c)");
         HdedxVsQP->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
         HdedxVsQP->SetAxisRange(-5,5,"X");
         HdedxVsQP->SetAxisRange(0,15,"Y");
         HdedxVsQP->Draw("COLZ");
       
//         KaonLine->Draw("same");
//         ProtonLine->Draw("same");
//         DeuteronLine->Draw("same");
//         TritonLine->Draw("same");
//         ProtonLineFit->Draw("same");

//         KaonLineLeft->Draw("same");
//         ProtonLineLeft->Draw("same");
//         DeuteronLineLeft->Draw("same");
//         TritonLineLeft->Draw("same");
      
         
         T->Draw("same");     
         SaveCanvas(c1, SaveDir, ObjName[i] + SaveName2 + "_dedxVsQP", true);
//         c1->SaveAs((SaveDir + ObjName[i] + SaveName2 + "_dedxVsQP.C"  ).c_str());
//         c1->SaveAs((SaveDir + ObjName[i] + SaveName2 + "_dedxVsQP.pdf").c_str());
         delete c1;

         c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         HdedxVsP_NS2->SetStats(kFALSE);
         HdedxVsP_NS2->GetXaxis()->SetTitle("track momentum (GeV/c)");
         HdedxVsP_NS2->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
         HdedxVsP_NS2->GetZaxis()->SetTitle("Number of saturated strips");
         HdedxVsP_NS2->SetAxisRange(0.0, 5.0, "X");
         HdedxVsP_NS2->Draw("COLZ");
         PionLine2->Draw("same");
         KaonLine2->Draw("same");
         ProtonLine2->Draw("same");
         DeuteronLine2->Draw("same");
//         TritonLine->Draw("same");
         ProtonLineFit2->Draw("same");
         T->Draw("same");
         SaveCanvas(c1, SaveDir, ObjName[i] + SaveName2 + "_dedxVsP_NS");
         delete c1;

         c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         HdedxVsEta2->SetStats(kFALSE);
         HdedxVsEta2->GetXaxis()->SetTitle("Eta");
         HdedxVsEta2->GetYaxis()->SetTitle("I_{as}");
         HdedxVsEta2->SetAxisRange(-2.1,2.1,"X");
         HdedxVsEta2->Draw("COLZ");
         SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_Eta2D");
         delete c1;
      }

      std::cout << "TESTB\n";

      SaveSuffix = InputFile2 ? SaveName + SaveName2 : "" ;

      c1 = new TCanvas("c1", "c1", 600,600);
      TLegend* leg = new TLegend (0.50, 0.80, 0.80, 0.90);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
      HdedxVsPProfile->SetStats(kFALSE);
      HdedxVsPProfile->SetAxisRange(2.5,5,"Y");
      HdedxVsPProfile->GetXaxis()->SetTitle("track momentum (GeV/c)");
      HdedxVsPProfile->GetYaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxVsPProfile->SetLineColor  (kBlack);
      HdedxVsPProfile->SetMarkerColor(kBlack);
      HdedxVsPProfile->SetMarkerStyle(20);
      HdedxVsPProfile->Draw("");
      if (InputFile2) {
         HdedxVsPProfile2->SetLineColor  (kBlue);
         HdedxVsPProfile2->SetMarkerColor(kBlue);
         HdedxVsPProfile2->SetMarkerStyle(22);
         HdedxVsPProfile2->Draw("same");
         leg->AddEntry (HdedxVsPProfile , SaveName .c_str(), "LP");
         leg->AddEntry (HdedxVsPProfile2, SaveName2.c_str(), "LP");
         leg->Draw();
      }
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_Profile");
      delete leg;
      delete c1;

      c1 = new TCanvas("c1", "c1", 600,600);
      leg = new TLegend (0.50, 0.80, 0.80, 0.90);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
      HdedxVsEtaProfile->SetStats(kFALSE);
      HdedxVsEtaProfile->GetXaxis()->SetTitle("#eta");
      HdedxVsEtaProfile->GetYaxis()->SetTitle("I_{as}");
      HdedxVsEtaProfile->SetAxisRange(-2.1,2.1,"X");
      HdedxVsEtaProfile->SetLineColor  (kBlack);
      HdedxVsEtaProfile->SetMarkerColor(kBlack);
      HdedxVsEtaProfile->SetMarkerStyle(20);
      HdedxVsEtaProfile->Draw("");
      if (InputFile2) {
         HdedxVsEtaProfile2->SetLineColor  (kBlue);
         HdedxVsEtaProfile2->SetMarkerColor(kBlue);
         HdedxVsEtaProfile2->SetMarkerStyle(22);
         HdedxVsEtaProfile2->Draw("same");
         leg->AddEntry (HdedxVsEtaProfile , SaveName .c_str(), "LP");
         leg->AddEntry (HdedxVsEtaProfile2, SaveName2.c_str(), "LP");
         leg->Draw();
      }
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_HdedxVsEtaProfile");
      delete leg;
      delete c1;

      c1 = new TCanvas("c1", "c1", 600,600);
      leg = new TLegend (0.50, 0.80, 0.80, 0.90);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
      c1->SetLogy(true);
      c1->SetGridx(true);
      HdedxMIP->SetStats(kFALSE);
      HdedxMIP->GetXaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxMIP->GetYaxis()->SetTitle("fraction of tracks");
      HdedxMIP->SetAxisRange(0,5,"X");
      HdedxMIP->SetLineColor(kBlack);
      HdedxMIP->SetLineWidth(2);
      HdedxMIP->Scale (1.0/HdedxMIP->Integral());
      HdedxMIP->Draw("");
      if (InputFile2) {
         HdedxMIP2->SetLineColor(kBlue);
         HdedxMIP2->SetLineWidth(2);
         HdedxMIP2->Scale (1.0/HdedxMIP2->Integral());
         HdedxMIP2->Draw("same");
         leg->AddEntry (HdedxMIP , SaveName .c_str(), "L");
         leg->AddEntry (HdedxMIP2, SaveName2.c_str(), "L");
         leg->Draw();
      }
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_MIP", true);
      delete leg;
      delete c1;

      c1 = new TCanvas("c1", "c1", 600,600);
      leg = new TLegend (0.50, 0.80, 0.80, 0.90);
		leg->SetFillColor(0);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
      c1->SetLogy(true);
      c1->SetGridx(true);
      HdedxSIG->SetStats(kFALSE);
      HdedxSIG->GetXaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxSIG->GetYaxis()->SetTitle("fraction of tracks");
      HdedxSIG->SetAxisRange(0,5,"X");
      HdedxSIG->SetLineColor(kBlack);
      HdedxSIG->SetLineWidth(2);
      HdedxSIG->Scale (1.0/HdedxSIG->Integral());
      HdedxSIG->Draw("");
      if (InputFile2) {
         HdedxSIG2->SetLineColor(kBlue);
         HdedxSIG2->SetLineWidth(2);
         HdedxSIG2->Scale (1.0/HdedxMIP2->Integral());
         HdedxSIG2->Draw("same");
         leg->AddEntry (HdedxSIG , SaveName .c_str(), "L");
         leg->AddEntry (HdedxSIG2, SaveName2.c_str(), "L");
         leg->Draw();
      }
      SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix + "_SIG", true);
      delete leg;
      delete c1;
/*
      c1 = new TCanvas("c1", "c1", 600,600);
      c1->SetLogy(true);
      c1->SetGridx(true);
      TLegend* leg = new TLegend (0.10, 0.80, 0.40, 0.90);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      HdedxSIG->SetStats(kFALSE);
      HdedxSIG->GetXaxis()->SetTitle(ObjName[i].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
      HdedxSIG->GetYaxis()->SetTitle("number of tracks");
      HdedxSIG->SetAxisRange(0,5,"X");
      HdedxMIP->SetLineColor (kBlue);
      HdedxMIP->Scale(1.0/HdedxMIP->Integral());
      HdedxSIG->Scale(1.0/HdedxSIG->Integral());
      leg->AddEntry (HdedxMIP, "Background", "L");
      leg->AddEntry (HdedxSIG, "Signal"    , "L");
      HdedxSIG->Draw("hist");
      HdedxMIP->Draw("same");
      leg->Draw();
      SaveCanvas(c1, SaveDir, ObjName[i] + "_SIGvsMIP", true);
      delete leg;
      delete c1;
*/
      std::cout << "TESTC\n";

      if (ObjName[i].find("harm2")!=std::string::npos){
            c1 = new TCanvas("c1", "c1", 600,600);
            leg = new TLegend (0.50, 0.80, 0.80, 0.90);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            c1->SetLogy(true);
            c1->SetGridx(true);
            HMass->Reset();
            for(int x=1;x<=HdedxVsP->GetNbinsX();x++){
               if(HdedxVsP->GetXaxis()->GetBinCenter(x)>3.0)continue;
               for(int y=1;y<=HdedxVsP->GetNbinsY();y++){
                  if(HdedxVsP->GetYaxis()->GetBinCenter(y)<5.0)continue;
                  HMass->Fill(GetMass (HdedxVsP->GetXaxis()->GetBinCenter(x),HdedxVsP->GetYaxis()->GetBinCenter(y)), HdedxVsP->GetBinContent(x,y));
               }
            }
            HMass->SetStats(kFALSE);
            HMass->GetXaxis()->SetTitle("Mass (GeV)");
            HMass->GetYaxis()->SetTitle("fraction of tracks");
            HMass->SetAxisRange(0,5,"X");
            HMass->Scale(1.0/HMass->Integral());
            HMass->SetLineColor  (kBlack);
            HMass->SetMarkerColor(kBlack);
            HMass->SetMarkerStyle(20);
            HMass->Draw("");
            if (InputFile2){
               HMass2->Reset();
               for(int x=1;x<=HdedxVsP->GetNbinsX();x++){
                  if(HdedxVsP2->GetXaxis()->GetBinCenter(x)>3.0)continue;
                  for(int y=1;y<=HdedxVsP2->GetNbinsY();y++){
                     if(HdedxVsP2->GetYaxis()->GetBinCenter(y)<5.0)continue;
                     HMass2->Fill(GetMass (HdedxVsP2->GetXaxis()->GetBinCenter(x),HdedxVsP2->GetYaxis()->GetBinCenter(y), 1), HdedxVsP2->GetBinContent(x,y));
                  }
               }
               HMass2->Scale(1.0/HMass2->Integral());
               HMass2->SetLineColor  (kBlue);
               HMass2->SetMarkerColor(kBlue);
               HMass2->SetMarkerStyle(22);
               HMass2->Draw("same");
               leg->AddEntry (HMass , SaveName .c_str(), "P");
               leg->AddEntry (HMass2, SaveName2.c_str(), "P");
               leg->Draw ();
            }
      
            double lineStart = InputFile2 ? std::min (HMass2->GetMinimum(), HMass->GetMinimum()) : HMass->GetMinimum(),
                   lineEnd   = InputFile2 ? std::max (HMass2->GetMaximum(), HMass->GetMaximum()) : HMass->GetMaximum();
            TLine* lineKaon = new TLine(0.493667, lineStart, 0.493667, lineEnd);
            lineKaon->SetLineWidth(2);
            lineKaon->SetLineStyle(2);
            lineKaon->SetLineColor(9);
            TLine* lineProton = new TLine(0.938272, lineStart, 0.938272, lineEnd);
            lineProton->SetLineWidth(2);
            lineProton->SetLineStyle(2);
            lineProton->SetLineColor(9);
            TLine* lineDeuteron = new TLine(1.88, lineStart, 1.88, lineEnd);
            lineDeuteron->SetLineWidth(2);
            lineDeuteron->SetLineStyle(2);
            lineDeuteron->SetLineColor(9);
            TLine* lineTriton = new TLine(2.80, lineStart, 2.80, lineEnd);
            lineTriton->SetLineWidth(2);
            lineTriton->SetLineStyle(2);
            lineTriton->SetLineColor(9);
      
            lineKaon->Draw("same");
            lineProton->Draw("same");
            lineDeuteron->Draw("same");
//            lineTriton->Draw("same");
            SaveCanvas(c1, SaveDir, ObjName[i] + SaveSuffix +  "_Mass", true);
            delete leg;
            delete c1;
      } else continue;
      
      std::cout << "TESTD\n";

   }
   ExtractConstantsReport.close();
   if (InputFile2) {
      ExtractConstantsReport2.close();

      vector <string> ObjNames;
      ObjNames.push_back ("Ias_SO");
      ObjNames.push_back ("Ias_SO_inc");
      ObjNames.push_back ("harm2_SO");
      ObjNames.push_back ("harm2_SO_in");
      ObjNames.push_back ("harm2_SP");

      TCanvas* c1   = new TCanvas ("c1", "c1", 600,600); 
      TLegend* leg  = new TLegend (0.50, 0.70, 0.80, 0.90);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      c1->SetLogy (true);
      TH1D h;
      h.GetXaxis()->SetTitle("signal efficiency");
      h.GetYaxis()->SetTitle("background efficiency");
      h.SetStats(0);
      h.Draw();
      TGraph** ROC = new TGraph* [ObjNames.size()];
      for (size_t NameIndex = 0; NameIndex < ObjNames.size(); NameIndex++)
      {
         TH1D* HdedxSIG1 = (TH1D*) GetObjectFromPath(InputFile , (ObjNames[NameIndex] + "_SIG").c_str() );
         TH1D* HdedxSIG2 = (TH1D*) GetObjectFromPath(InputFile2, (ObjNames[NameIndex] + "_SIG").c_str() );
         ROC[NameIndex]  = new TGraph(HdedxSIG1->GetNbinsX());

         double fullBkg  = HdedxSIG1->Integral(),
                fullSig  = HdedxSIG2->Integral();
         for (unsigned int cut_i = 1; cut_i <= HdedxSIG1->GetNbinsX(); cut_i++)
            ROC[NameIndex]->SetPoint (cut_i-1, 1 - HdedxSIG2->Integral(1, cut_i)/fullSig, 1 - HdedxSIG1->Integral(1, cut_i)/fullBkg);

	      ROC[NameIndex]->SetLineColor   (NameIndex);
	      ROC[NameIndex]->SetLineWidth   (2);
	      ROC[NameIndex]->Draw("same");

         leg->AddEntry (ROC[NameIndex], ObjNames[NameIndex].c_str(), "L");
         HdedxSIG1->~TH1D(); HdedxSIG2->~TH1D();
      }
      leg->Draw();
      SaveCanvas(c1, SaveDir, "Comparison_ROC");
      for (size_t NameIndex = 0; NameIndex < ObjNames.size(); NameIndex++)
         delete ROC[NameIndex];
      delete ROC;
      delete leg;
      delete c1;
   }


//   getScaleFactor(InputFile, "All_Rescale", "harm2_SO", InputFile2);
//   if (InputFile2) getScaleFactor(InputFile, "All_Rescale_" + SaveName + SaveName2, "harm2_SO_raw", InputFile2);
//   getScaleFactor(InputFile, "All_Rescale", "harm2_SO_raw", NULL, "harm2_PO_raw");
}



void getScaleFactor(TFile* InputFile1, string OutName, string ObjName1, TFile* InputFile2, string ObjName2){
   TProfile*   HdedxVsPProfile1;
   TProfile*   HdedxVsPProfile2;
   TH1D*       HdedxMIP1;
   TH1D*       HdedxMIP2;

   if (InputFile2!=NULL){
      HdedxVsPProfile1 = (TProfile*)GetObjectFromPath(InputFile1, (ObjName1 + "_Profile"  ).c_str() );
      HdedxVsPProfile2 = (TProfile*)GetObjectFromPath(InputFile2, (ObjName1 + "_Profile"  ).c_str() );

      HdedxMIP1        = (TProfile*)GetObjectFromPath(InputFile1, (ObjName1 + "_MIP"  ).c_str() );
      HdedxMIP2        = (TProfile*)GetObjectFromPath(InputFile2, (ObjName1 + "_MIP"  ).c_str() );
   } else if (ObjName2.find("EMPTY")==string::npos) {
      HdedxVsPProfile1 = (TProfile*)GetObjectFromPath(InputFile1, (ObjName1 + "_Profile"  ).c_str() );
      HdedxVsPProfile2 = (TProfile*)GetObjectFromPath(InputFile1, (ObjName2 + "_Profile"  ).c_str() );

      HdedxMIP1        = (TProfile*)GetObjectFromPath(InputFile1, (ObjName1 + "_MIP"  ).c_str() );
      HdedxMIP2        = (TProfile*)GetObjectFromPath(InputFile1, (ObjName2 + "_MIP"  ).c_str() );
   } else return;

   TF1* mygausMIP = new TF1("mygausMIP","gaus", 1, 5);
   HdedxMIP1->Fit("mygausMIP","Q0","");
   double peakMIP1  = mygausMIP->GetParameter(1);
   HdedxMIP2->Fit("mygausMIP","Q0","");
   double peakMIP2  = mygausMIP->GetParameter(1);

   std::cout << "SCALE FACTOR WITH MIP     = " << peakMIP1/peakMIP2 << endl;

   TH1D* Chi2Dist = new TH1D("Chi2Dist","Chi2Dist",300, 0.9 ,1.15);

   double Minimum = 999999;
   double AbsGain = -1;

   for(int i=1;i<=Chi2Dist->GetNbinsX();i++){
      double ScaleFactor = Chi2Dist->GetXaxis()->GetBinCenter(i);
      TProfile* Rescaled = (TProfile*)HdedxVsPProfile2->Clone("Cloned");
      Rescaled->Scale(ScaleFactor);
      double Dist = 0;
      double Error = 0;
      for(int x=1;x<=HdedxVsPProfile1->GetNbinsX();x++){
         double Momentum = HdedxVsPProfile1->GetXaxis()->GetBinCenter(x);
         if(Momentum<5)continue;//|| Momentum>20)continue;
         if(HdedxVsPProfile1->GetBinError(x)<=0)continue;
         Dist += pow(HdedxVsPProfile1->GetBinContent(x) - Rescaled->GetBinContent(x),2) / std::max(1E-8,pow(HdedxVsPProfile1->GetBinError(x),2));
         Error += pow(HdedxVsPProfile1->GetBinError(x),2);
      }
      Dist *= Error;

      if(Dist<Minimum){Minimum=Dist;AbsGain=ScaleFactor;}

      //std::cout << "Rescale = " << ScaleFactor << " --> SquareDist = " << Dist << endl;
      Chi2Dist->Fill(ScaleFactor,Dist);
      delete Rescaled;
   }

   std::cout << "SCALE FACTOR WITH PROFILE = " << AbsGain << endl;

   TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
   HdedxMIP2->SetStats(kFALSE);
   HdedxMIP2->SetAxisRange(0,10,"X");
   HdedxMIP2->GetXaxis()->SetNdivisions(516);
   HdedxMIP2->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
   HdedxMIP2->GetYaxis()->SetTitle("Tracks");
   HdedxMIP2->SetLineColor(1);
   HdedxMIP2->Draw("");
   TH1D* HdedxMIP3 = (TH1D*)HdedxMIP2->Clone("aaa");
   HdedxMIP3->SetLineColor(8);
   HdedxMIP3->GetXaxis()->Set(HdedxMIP3->GetXaxis()->GetNbins(), HdedxMIP3->GetXaxis()->GetXmin()*2.0, HdedxMIP3->GetXaxis()->GetXmax()*(peakMIP1/peakMIP2) );
   HdedxMIP3->Draw("same");
   TH1D* HdedxMIP4 = (TH1D*)HdedxMIP2->Clone("bbb");
   HdedxMIP4->SetLineColor(4);
   HdedxMIP4->GetXaxis()->Set(HdedxMIP4->GetXaxis()->GetNbins(), HdedxMIP4->GetXaxis()->GetXmin()*2.0, HdedxMIP4->GetXaxis()->GetXmax()*(AbsGain) );
   HdedxMIP4->Draw("same");
   HdedxMIP1->SetLineColor(2);
   HdedxMIP1->Draw("same");
   c1->SetLogy(true);
   c1->SetGridx(true); 
   SaveCanvas(c1, "pictures/", OutName + "_MIP");
   delete c1;


   c1 = new TCanvas("c1", "c1", 600,600);
   Chi2Dist->SetStats(kFALSE);
   Chi2Dist->GetXaxis()->SetNdivisions(504);
   Chi2Dist->GetXaxis()->SetTitle("Rescale Factor");
   Chi2Dist->GetYaxis()->SetTitle("Weighted Square Distance");
   Chi2Dist->Draw("");
   c1->SetLogy(true);
   c1->SetGridx(true); 
   SaveCanvas(c1, "pictures/", OutName + "_Dist");
   delete c1;

   c1 = new TCanvas("c1", "c1", 600,600);
   HdedxVsPProfile1->SetStats(kFALSE);
   HdedxVsPProfile1->SetAxisRange(5,50,"X");
   HdedxVsPProfile1->SetAxisRange(2.5,3.5,"Y");
   HdedxVsPProfile1->GetXaxis()->SetTitle("track momentum (GeV/c)");
   HdedxVsPProfile1->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
   HdedxVsPProfile1->SetMarkerColor(2);
   HdedxVsPProfile1->Draw("");

   HdedxVsPProfile2->SetMarkerColor(1);
   HdedxVsPProfile2->Draw("same");
   TProfile* HdedxVsPProfile3 = (TProfile*)HdedxVsPProfile2->Clone("abc");
   HdedxVsPProfile3->SetMarkerColor(8);
   HdedxVsPProfile3->Scale(peakMIP1/peakMIP2);
   HdedxVsPProfile3->Draw("same");
   TProfile* HdedxVsPProfile4 = (TProfile*)HdedxVsPProfile2->Clone("afs");
   HdedxVsPProfile4->SetMarkerColor(4);
   HdedxVsPProfile4->Scale(AbsGain);
   HdedxVsPProfile4->Draw("same");

   SaveCanvas(c1, "pictures/", OutName + "_Profile");
   delete c1;
}



void ExtractConstants(TH2D* input, int FileIndex){
       double MinRange = 0.60;
       double MaxRange = 1.20;
       char buffer[2048];
       bool hasConverged = false;

       for(unsigned int loop=0;loop<5 and !hasConverged; loop++){
	      TH2D* inputnew = (TH2D*)input->Clone("tempTH2D");
	      inputnew->Rebin2D(5,10);
	      for(int x=1;x<=inputnew->GetNbinsX();x++){
	      for(int y=1;y<=inputnew->GetNbinsY();y++){
		double Mass = GetMass(inputnew->GetXaxis()->GetBinCenter(x),inputnew->GetYaxis()->GetBinCenter(y));
		if(isnan(float(Mass)) || Mass<0.94-0.3 || Mass>0.94+0.3)inputnew->SetBinContent(x,y,0);        
	      }}

	      TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
	      c1->SetLogz(true);
	      inputnew->SetStats(kFALSE);
	      inputnew->GetXaxis()->SetTitle("track momentum (GeV/c)");
	      inputnew->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      inputnew->SetAxisRange(0,5,"X");
	      inputnew->SetAxisRange(0,15,"Y");
	      inputnew->Draw("COLZ");

	//      KaonLine->Draw("same");
	//      ProtonLine->Draw("same");
	//      DeuteronLine->Draw("same");
	//      TritonLine->Draw("same");
	      SaveCanvas(c1, "fit/", "dedxVsP");
	      delete c1;


	       TH1D* FitResult = new TH1D("FitResult"       , "FitResult"      ,inputnew->GetXaxis()->GetNbins(),inputnew->GetXaxis()->GetXmin(),inputnew->GetXaxis()->GetXmax());
	       FitResult->SetTitle("");
	       FitResult->SetStats(kFALSE);  
	       FitResult->GetXaxis()->SetTitle("P [GeV/c]");
	       FitResult->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
	       FitResult->GetYaxis()->SetTitleOffset(1.20);
	       FitResult->Reset();     

	       for(int x=1;x<inputnew->GetXaxis()->FindBin(5);x++){
		  double P       = inputnew->GetXaxis()->GetBinCenter(x);
	    
		  TH1D* Projection = (TH1D*)(inputnew->ProjectionY("proj",x,x))->Clone();
		  if(Projection->Integral()<100)continue;
		  Projection->SetAxisRange(0.1,25,"X");
		  Projection->Sumw2();
		  Projection->Scale(1.0/Projection->Integral());

		  TF1* mygaus = new TF1("mygaus","gaus", 2.5, 15);
		  Projection->Fit("mygaus","Q0 RME");
		  double chiFromFit  = (mygaus->GetChisquare())/(mygaus->GetNDF());
		  FitResult->SetBinContent(x, mygaus->GetParameter(1));
		  FitResult->SetBinError  (x, mygaus->GetParError (1));
		  mygaus->SetLineColor(2);
		  mygaus->SetLineWidth(2);

		  c1  = new TCanvas("canvas", "canvas", 600,600);
		  Projection->Draw();
		  Projection->SetTitle("");
		  Projection->SetStats(kFALSE);
		  Projection->GetXaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
		  Projection->GetYaxis()->SetTitle("#Entries");
		  Projection->GetYaxis()->SetTitleOffset(1.30);
		  Projection->SetAxisRange(1E-5,1.0,"Y");

		  mygaus->Draw("same");


		  TPaveText* stt = new TPaveText(0.55,0.82,0.79,0.92, "NDC");
		  stt->SetFillColor(0);
		  stt->SetTextAlign(31);
		  sprintf(buffer,"Proton  #mu:%5.1fMeV/cm",mygaus->GetParameter(1));      stt->AddText(buffer);
		  sprintf(buffer,"Proton  #sigma:%5.1fMeV/cm",mygaus->GetParameter(2));      stt->AddText(buffer);
		  stt->Draw("same");

		  //std::cout << "P = " << P << "  --> Proton dE/dx = " << mygaus->GetParameter(1) << endl;

		  c1->SetLogy(true);
		  sprintf(buffer,"%sProjectionFit_P%03i_%03i","fit/",(int)(100*FitResult->GetXaxis()->GetBinLowEdge(x)),(int)(100*FitResult->GetXaxis()->GetBinUpEdge(x)) );
		  if(P>=MinRange && P<=MaxRange){SaveCanvas(c1,"./",buffer);}
		  delete c1;
                  delete Projection;
                  delete mygaus;
                  delete stt;
	       }
	       c1  = new TCanvas("canvas", "canvas", 600,600);
	       FitResult->SetAxisRange(0,2.5,"X");
	       FitResult->SetAxisRange(0,15,"Y");
	       FitResult->Draw("");

	       TLine* line1 = new TLine(MinRange, FitResult->GetMinimum(), MinRange, FitResult->GetMaximum());
	       line1->SetLineWidth(2);
	       line1->SetLineStyle(2);
	       line1->Draw();

	       TLine* line2 = new TLine(MaxRange, FitResult->GetMinimum(), MaxRange, FitResult->GetMaximum());
	       line2->SetLineWidth(2);
	       line2->SetLineStyle(2);
	       line2->Draw();

	       //   TF1* myfit = new TF1("myfit","[1]+(pow(0.93827,2) + x*x)/([0]*x*x)", MinRange, MaxRange);
	       TF1* myfit = new TF1("myfit","[0]*pow(0.93827/x,2) + [1]", MinRange, MaxRange);
	       myfit->SetParName  (0,"K");
	       myfit->SetParName  (1,"C");
	       myfit->SetParameter(0, 2.7);
	       myfit->SetParameter(1, 2.7);
	       myfit->SetParLimits(0, 1.00,4.0);
	       myfit->SetParLimits(1, 1.00,4.0);
	       myfit->SetLineWidth(2);
	       myfit->SetLineColor(2);
	       FitResult->Fit("myfit", "M R E I 0");
	       myfit->SetRange(MinRange,MaxRange);
	       myfit->Draw("same");

	       double prevConstants [] = {K[FileIndex], Kerr[FileIndex], C[FileIndex], Cerr[FileIndex]};
	       K   [FileIndex] = myfit->GetParameter(0);
	       C   [FileIndex] = myfit->GetParameter(1);
	       Kerr[FileIndex] = myfit->GetParError(0);
	       Cerr[FileIndex] = myfit->GetParError(1);

	       printf("K Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n", prevConstants[0], prevConstants[1], K[FileIndex], Kerr[FileIndex], 100.0*(K[FileIndex]-prevConstants[0])/K[FileIndex]);
	       printf("C Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n", prevConstants[2], prevConstants[3], C[FileIndex], Cerr[FileIndex], 100.0*(C[FileIndex]-prevConstants[2])/C[FileIndex]);

          if(std::max(fabs(100.0*(K[FileIndex]-prevConstants[0])/K[FileIndex]), fabs(100.0*(C[FileIndex]-prevConstants[2])/C[FileIndex]))<1.0)hasConverged=true;  //<1% variation of the constant --> converged


	       TPaveText* st = new TPaveText(0.40,0.78,0.79,0.89, "NDC");
	       st->SetFillColor(0);
	       sprintf(buffer,"K = %4.3f +- %6.4f",myfit->GetParameter(0), myfit->GetParError(0));
	       st->AddText(buffer);
	       sprintf(buffer,"C = %4.3f +- %6.4f",myfit->GetParameter(1), myfit->GetParError(1));
	       st->AddText(buffer);
	       st->Draw("same");
	       sprintf(buffer,"%sFit","fit/");
	       SaveCanvas(c1,"./",buffer);              
	       delete c1;

          delete line1;
          delete line2;
          delete myfit;
          delete FitResult;
          delete inputnew;
       }
}

void DrawComparisons (TFile* InputFile1, TFile* InputFile2, string ObjName1, string ObjName2){
	TProfile*   HdedxVsEtaProfile1  = (TProfile*)  GetObjectFromPath(InputFile1, (ObjName1 + "_Eta" ).c_str() );
	TProfile*   HdedxVsEtaProfile2  = (TProfile*)  GetObjectFromPath(InputFile1, (ObjName2 + "_Eta" ).c_str() );
	TH1D*       HdedxMIP1           = (TH1D*)      GetObjectFromPath(InputFile1, (ObjName1 + "_MIP" ).c_str() );
	TH1D*       HdedxMIP2           = (TH1D*)      GetObjectFromPath(InputFile1, (ObjName2 + "_MIP" ).c_str() );
	
	TCanvas* c1  = new TCanvas("c1", "c1", 600,600);
	TLegend* leg = new TLegend(0.50, 0.80, 0.80, 0.90);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry (HdedxVsEtaProfile1, ObjName1.c_str(), "P");
	leg->AddEntry (HdedxVsEtaProfile2, ObjName2.c_str(), "P");
	HdedxVsEtaProfile1->SetStats(kFALSE);
	HdedxVsEtaProfile2->SetMarkerStyle(23);
	HdedxVsEtaProfile1->SetMarkerColor(kBlack);
	HdedxVsEtaProfile2->SetMarkerColor(kBlue);
	HdedxVsEtaProfile1->GetXaxis()->SetTitle("#eta");
	HdedxVsEtaProfile1->GetYaxis()->SetTitle(ObjName1.find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
	HdedxVsEtaProfile1->Draw("");
	HdedxVsEtaProfile2->Draw("same");
	leg->Draw();
	SaveCanvas(c1, "pictures/", "Comparison_"+ObjName1+"_"+ObjName2+"_HdedxVsEtaProfile");
	delete leg;
	delete c1;

	c1 = new TCanvas("c1", "c1", 600,600);
	leg = new TLegend (0.50, 0.80, 0.80, 0.90);
	c1->SetLogy(true);
	c1->SetGridx(true);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	HdedxMIP1->SetStats(kFALSE);
	HdedxMIP1->GetXaxis()->SetTitle(ObjName1.find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
	HdedxMIP1->GetYaxis()->SetTitle("fraction of tracks");
	HdedxMIP1->SetAxisRange(0,5,"X");
	HdedxMIP1->SetLineColor (kBlack);
	HdedxMIP2->SetLineColor (kBlue);
	HdedxMIP2->Scale(1.0/HdedxMIP2->Integral());
	HdedxMIP1->Scale(1.0/HdedxMIP1->Integral());
	leg->AddEntry (HdedxMIP1, ObjName1.c_str(), "L");
	leg->AddEntry (HdedxMIP2, ObjName2.c_str(), "L");
	HdedxMIP1->Draw("hist");
	HdedxMIP2->Draw("same");
	leg->Draw();
	SaveCanvas(c1, "pictures/", "Comparison_" + ObjName1 + "_" + ObjName2 + "_MIP", true);
	delete leg;
	delete c1;
}
