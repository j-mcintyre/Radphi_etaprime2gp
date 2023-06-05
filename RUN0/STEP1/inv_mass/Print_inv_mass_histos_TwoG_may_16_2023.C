// Author: J. McIntyre
// 16 May 2023
// This file prints png files of selected histograms and
// replaces selected parameters (i.e. title, axis name, etc.)

#include <TProofServ.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>
#include <Math/MinimizerOptions.h>
#include <TStyle.h>
#include <TMath.h>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TH2.h>
#include <math.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TSystem.h>
#include <TImage.h>

#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLatex.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TROOT.h>

#include <TLegend.h>

//using namespace std;       // Removes the need for std::

TH1D *hcut;
TH1D *hcut0;
TH1D *hcut1;
TH1D *hcut2;
TH1D *hcut3;
TH1D *hcut4;
TH1D *hcut5;
TH1D *hcut6;
TH1D *hcut7;
TH1D *hcut8;

TH1D *h_ncoin;
TH1D *h_dt;
TH1D *h2g_all;
TH1D *h_dcotime;
TH1D *h_dcotime2;
TH1D *t_2gamma;
TH1D *M_2gamma;
TH1D *TChan;

TH1D *M2g_E0_[3];
TH1D *M2g_E1_[3];
TH1D *M2g_E2_[3];
TH1D *M2g_E3_[3];

TH1D *Allt_E0;
TH1D *Allt_E1;
TH1D *Allt_E2;
TH1D *Allt_E3;

void Print_inv_mass_histos_TwoG_may_16_2023()
{
   /*
   double norm2g = 1.0;
   double norm3g = norm2 / norm3;
   double norm4g = norm2 / norm4;
   double norm5g = norm2 / norm5;
   double norm6g = norm2 / norm6;
   double norm7g = norm2 / norm7;

   //hcut2->GetXaxis()->SetRangeUser(0.55, 1.0);

   // Setting canvas size for (2,1) pad configuration
   int xsplit_canvas = 800;
   int ysplit_canvas = 800;

   // Setting canvas size for (2,1) pad configuration
   int xsplit3_canvas = 800;
   int ysplit3_canvas = 1000;

   // Setting canvas size for single pad configuration
   int x_canvas = 800;
   int y_canvas = 600;
   */

   //TString title;
   TString name0;
   TString name1;
   TString name2;
   TString name3;
   TString name4;
   TString name5;
   TString name6;
   TString name7;
   TString name8;
   TString name9;
   TString name10;
   TString name11;
   TString name12;
   TString name13;
   TString name14;
   TString name15;
   TString name16;

   TFile *rootf1 = new TFile("TwoG_may_16_2023Analysis.root");
   hcut = (TH1D*)rootf1->Get("hcut");
   hcut0 = (TH1D*)rootf1->Get("hcut0");
   hcut1 = (TH1D*)rootf1->Get("hcut1");
   hcut2 = (TH1D*)rootf1->Get("hcut2");
   hcut3 = (TH1D*)rootf1->Get("hcut3");
   hcut4 = (TH1D*)rootf1->Get("hcut4");
   hcut5 = (TH1D*)rootf1->Get("hcut5");
   hcut6 = (TH1D*)rootf1->Get("hcut6");
   hcut7 = (TH1D*)rootf1->Get("hcut7");
   hcut8 = (TH1D*)rootf1->Get("hcut8");

   h_ncoin = (TH1D*)rootf1->Get("h_ncoin");
   h_dt = (TH1D*)rootf1->Get("h_dt");
   h2g_all = (TH1D*)rootf1->Get("h2g_all");
   h_dcotime = (TH1D*)rootf1->Get("h_dcotime");
   h_dcotime2 = (TH1D*)rootf1->Get("h_dcotime2");
   t_2gamma = (TH1D*)rootf1->Get("t_2gamma");
   M_2gamma = (TH1D*)rootf1->Get("M_2gamma");
   TChan = (TH1D*)rootf1->Get("TChan");

   Allt_E0 = (TH1D*)rootf1->Get("Allt_E0");
   Allt_E1 = (TH1D*)rootf1->Get("Allt_E1");
   Allt_E2 = (TH1D*)rootf1->Get("Allt_E2");
   Allt_E3 = (TH1D*)rootf1->Get("Allt_E3");

   for (int i = 0; i < 3; i++) {
      name0.Form("M2g_E0_%d", i);
      M2g_E0_[i] = (TH1D*)rootf1->Get(name0);

      name1.Form("M2g_E1_%d", i);
      M2g_E1_[i] = (TH1D*)rootf1->Get(name1);

      name2.Form("M2g_E2_%d", i);
      M2g_E2_[i] = (TH1D*)rootf1->Get(name2);

      name3.Form("M2g_E3_%d", i);
      M2g_E3_[i] = (TH1D*)rootf1->Get(name3);
   }

   TString PName;
   TCanvas *c_100 = new TCanvas("c_100","", 800, 600);
   c_100->cd();
   gStyle->SetOptStat(11);
   //h->SetTitle("Varying Energy Rescaler");
   //h->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
   //h->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
   //TPaveStats *st_100 = (TPaveStats*)hcut->FindObject("stats");
   //h->GetXaxis()->SetTitleOffset(1.2);
   //h->GetYaxis()->SetTitleOffset(1.3);
   //hcut1->SetMinimum(0.0);
   hcut->GetXaxis()->SetRangeUser(0.0, 1.2);
   hcut->SetLineColor(1);
   hcut->SetFillColor(17);
   hcut->Draw("hist");
   PName.Form("M2g_hcut.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_hcut_LogScale.png");
   hcut->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   Allt_E0->GetXaxis()->SetRangeUser(0.0, 1.2);
   Allt_E0->SetLineColor(1);
   Allt_E0->SetFillColor(17);
   Allt_E0->Draw("hist");
   PName.Form("M2g_Allt_E0.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_Allt_E0_LogScale.png");
   Allt_E0->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   Allt_E1->GetXaxis()->SetRangeUser(0.0, 1.2);
   Allt_E1->SetLineColor(1);
   Allt_E1->SetFillColor(17);
   Allt_E1->Draw("hist");
   PName.Form("M2g_Allt_E1.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_Allt_E1_LogScale.png");
   Allt_E1->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   Allt_E2->GetXaxis()->SetRangeUser(0.0, 1.2);
   Allt_E2->SetLineColor(1);
   Allt_E2->SetFillColor(17);
   Allt_E2->Draw("hist");
   PName.Form("M2g_Allt_E2.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_Allt_E2_LogScale.png");
   Allt_E2->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   Allt_E3->GetXaxis()->SetRangeUser(0.0, 1.2);
   Allt_E3->SetLineColor(1);
   Allt_E3->SetFillColor(17);
   Allt_E3->Draw("hist");
   PName.Form("M2g_Allt_E3.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_Allt_E3_LogScale.png");
   Allt_E3->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   h_ncoin->SetLineColor(1);
   h_ncoin->SetFillColor(17);
   h_ncoin->Draw("hist");
   PName.Form("M2g_h_ncoin.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_h_ncoin_LogScale.png");
   h_ncoin->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   h_dt->SetLineColor(1);
   h_dt->SetFillColor(17);
   h_dt->Draw("hist");
   PName.Form("M2g_h_dt.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_h_dt_LogScale.png");
   h_dt->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   h2g_all->SetLineColor(1);
   h2g_all->SetFillColor(17);
   h2g_all->Draw("hist");
   int xx = h2g_all->GetMaximumBin();
   double MaxY = h2g_all->GetBinContent(xx);
   int NBins = h2g_all->GetNbinsX();
   int x0 = 0.1 * NBins;
   int xf = 0.9 * NBins;
   double xx0 = h2g_all->GetBinCenter(x0);
   double xxf = h2g_all->GetBinCenter(xf);
   TLine *zeroline4 = new TLine(xx0, MaxY, xxf, MaxY);
   zeroline4->SetLineColor(2);
   zeroline4->SetLineStyle(7);
   zeroline4->SetLineWidth(1);
   zeroline4->Draw("same");
   PName.Form("M2g_h2g_all.png");
   c_100->Print(PName);

   //
   h_dcotime->SetLineColor(1);
   h_dcotime->SetFillColor(17);
   h_dcotime->Draw("hist");
   xx = h_dcotime->GetMaximumBin();
   MaxY = h_dcotime->GetBinContent(xx);
   NBins = h_dcotime->GetNbinsX();
   x0 = 0.1 * NBins;
   xf = 0.9 * NBins;
   xx0 = h_dcotime->GetBinCenter(x0);
   xxf = h_dcotime->GetBinCenter(xf);
   TLine *zeroline1 = new TLine(xx0, MaxY, xxf, MaxY);
   zeroline1->SetLineColor(2);
   zeroline1->SetLineStyle(7);
   zeroline1->SetLineWidth(1);
   zeroline1->Draw("same");
   PName.Form("M2g_h_dcotime.png");
   c_100->Print(PName);


   //
   h_dcotime2->SetLineColor(1);
   h_dcotime2->SetFillColor(17);
   h_dcotime2->Draw("hist");
   xx = h_dcotime2->GetMaximumBin();
   MaxY = h_dcotime2->GetBinContent(xx);
   NBins = h_dcotime2->GetNbinsX();
   x0 = 0.1 * NBins;
   xf = 0.9 * NBins;
   xx0 = h_dcotime2->GetBinCenter(x0);
   xxf = h_dcotime2->GetBinCenter(xf);
   TLine *zeroline2 = new TLine(xx0, MaxY, xxf, MaxY);
   zeroline2->SetLineColor(2);
   zeroline2->SetLineStyle(7);
   zeroline2->SetLineWidth(1);
   zeroline2->Draw("same");
   PName.Form("M2g_h_dcotime2.png");
   c_100->Print(PName);

   //
   t_2gamma->SetLineColor(1);
   t_2gamma->SetFillColor(17);
   t_2gamma->Draw("hist");
   PName.Form("M2g_t_2gamma.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_t_2gamma_LogScale.png");
   t_2gamma->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   M_2gamma->GetXaxis()->SetRangeUser(0.0, 1.2);
   M_2gamma->SetLineColor(1);
   M_2gamma->SetFillColor(17);
   M_2gamma->Draw("hist");
   PName.Form("M2g_M_2gamma.png");
   c_100->Print(PName);
   c_100->SetLogy(1);
   PName.Form("M2g_M_2gamma_LogScale.png");
   M_2gamma->Draw("hist");
   c_100->Update();
   c_100->Print(PName);
   c_100->SetLogy(0);
   c_100->Update();

   //
   TChan->SetLineColor(1);
   TChan->SetFillColor(17);
   TChan->Draw("hist");
   PName.Form("M2g_TChan.png");
   c_100->Print(PName);

   //
   char Ename[30];
   char Etitle[30];
   char EtitleLog[30];
   for (int j = 0; j < 9; j++) {
      sprintf(Ename,"%s%d","hcut",j);
      TH1D *h = (TH1D*)gDirectory->Get(Ename);
      h->GetXaxis()->SetRangeUser(0.0, 1.2);
      h->SetLineColor(1);
      h->SetFillColor(17);

      // Draw histogram with logarithmic scale
      c_100->SetLogy(1);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG with logarithmic scale
      sprintf(EtitleLog,"%s%d%s","M2g_hcut",j,"_LogScale.png");
      TImage *r6 = TImage::Create();
      r6->FromPad(c_100);
      r6->WriteImage(EtitleLog);

      // Reset logarithmic scale and draw histogram again
      c_100->SetLogy(0);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG without logarithmic scale
      sprintf(Etitle,"%s%d%s","M2g_hcut",j,".png");
      TImage *i6 = TImage::Create();
      i6->FromPad(c_100);
      i6->WriteImage(Etitle);
   }

   //
   for (int k = 0; k < 3; k++) {
      sprintf(Ename,"%s%d","M2g_E0_",k);
      TH1D *h = (TH1D*)gDirectory->Get(Ename);
      h->GetXaxis()->SetRangeUser(0.0, 1.2);
      h->SetLineColor(1);
      h->SetFillColor(17);

      // Draw histogram with logarithmic scale
      c_100->SetLogy(1);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG with logarithmic scale
      sprintf(EtitleLog,"%s%d%s","M2g_E0_",k,"_LogScale.png");
      TImage *r6 = TImage::Create();
      r6->FromPad(c_100);
      r6->WriteImage(EtitleLog);

      // Reset logarithmic scale and draw histogram again
      c_100->SetLogy(0);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG without logarithmic scale
      sprintf(Etitle,"%s%d%s","M2g_E0_",k,".png");
      TImage *i6 = TImage::Create();
      i6->FromPad(c_100);
      i6->WriteImage(Etitle);
   }


   //
   for (int k = 0; k < 3; k++) {
      sprintf(Ename,"%s%d","M2g_E1_",k);
      TH1D *h = (TH1D*)gDirectory->Get(Ename);
      h->GetXaxis()->SetRangeUser(0.0, 1.2);
      h->SetLineColor(1);
      h->SetFillColor(17);

      // Draw histogram with logarithmic scale
      c_100->SetLogy(1);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG with logarithmic scale
      sprintf(EtitleLog,"%s%d%s","M2g_E1_",k,"_LogScale.png");
      TImage *r6 = TImage::Create();
      r6->FromPad(c_100);
      r6->WriteImage(EtitleLog);

      // Reset logarithmic scale and draw histogram again
      c_100->SetLogy(0);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG without logarithmic scale
      sprintf(Etitle,"%s%d%s","M2g_E1_",k,".png");
      TImage *i6 = TImage::Create();
      i6->FromPad(c_100);
      i6->WriteImage(Etitle);
   }

   //
   for (int k = 0; k < 3; k++) {
      sprintf(Ename,"%s%d","M2g_E2_",k);
      TH1D *h = (TH1D*)gDirectory->Get(Ename);
      h->GetXaxis()->SetRangeUser(0.0, 1.2);
      h->SetLineColor(1);
      h->SetFillColor(17);

      // Draw histogram with logarithmic scale
      c_100->SetLogy(1);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG with logarithmic scale
      sprintf(EtitleLog,"%s%d%s","M2g_E2_",k,"_LogScale.png");
      TImage *r6 = TImage::Create();
      r6->FromPad(c_100);
      r6->WriteImage(EtitleLog);

      // Reset logarithmic scale and draw histogram again
      c_100->SetLogy(0);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG without logarithmic scale
      sprintf(Etitle,"%s%d%s","M2g_E2_",k,".png");
      TImage *i6 = TImage::Create();
      i6->FromPad(c_100);
      i6->WriteImage(Etitle);
   }

   //
   for (int k = 0; k < 3; k++) {
      sprintf(Ename,"%s%d","M2g_E3_",k);
      TH1D *h = (TH1D*)gDirectory->Get(Ename);
      h->GetXaxis()->SetRangeUser(0.0, 1.2);
      h->SetLineColor(1);
      h->SetFillColor(17);

      // Draw histogram with logarithmic scale
      c_100->SetLogy(1);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG with logarithmic scale
      sprintf(EtitleLog,"%s%d%s","M2g_E3_",k,"_LogScale.png");
      TImage *r6 = TImage::Create();
      r6->FromPad(c_100);
      r6->WriteImage(EtitleLog);

      // Reset logarithmic scale and draw histogram again
      c_100->SetLogy(0);
      h->Draw("hist");
      c_100->Update();

      // Save histogram as PNG without logarithmic scale
      sprintf(Etitle,"%s%d%s","M2g_E3_",k,".png");
      TImage *i6 = TImage::Create();
      i6->FromPad(c_100);
      i6->WriteImage(Etitle);
   }


#if 0
   TCanvas *c_908 = new TCanvas("c_908","", 800, 600);
   c_908->cd();
   //hmass->SetTitle("");
   hmass->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
   hmass->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
   gStyle->SetOptStat(11);
   TPaveStats *st_908 = (TPaveStats*)hmass->FindObject("stats");
   hmass->SetLineColor(1);
   hmass->SetFillColor(17);
   hmass->GetXaxis()->SetTitleOffset(1.2);
   hmass->GetYaxis()->SetTitleOffset(1.5);
   hmass->GetXaxis()->SetRangeUser(0.0, 1.2);
   hmass->SetMinimum(0.0);
   hmass->Draw("hist");
   c_908->Update();
   TImage *img_908 = TImage::Create();
   img_908->FromPad(c_908);
   img_908->WriteImage("h2g_hmass.png");
   delete c_908;
   delete img_908;
#endif

   delete c_100;

   // Organize folder info
   gSystem->Exec("mkdir -p histos/Cuts");
   gSystem->Exec("mkdir -p histos/E_Binned/t_Binned");
   gSystem->Exec("mkdir -p histos/E_Binned/All_t");
   gSystem->Exec("mkdir -p histos/Cuts/LogPlots");
   gSystem->Exec("mkdir -p histos/E_Binned/t_Binned/LogPlots");
   gSystem->Exec("mkdir -p histos/E_Binned/All_t/LogPlots");

   gSystem->Exec("mv M2g_hcut* histos/Cuts/.");
   gSystem->Exec("mv histos/Cuts/*_LogScale.png histos/Cuts/LogPlots/.");
   gSystem->Exec("mv M2g_E* histos/E_Binned/t_Binned/.");
   gSystem->Exec("mv histos/E_Binned/t_Binned/*_LogScale.png histos/E_Binned/t_Binned/LogPlots/.");
   gSystem->Exec("mv M2g_Allt_E* histos/E_Binned/All_t/.");
   gSystem->Exec("mv histos/E_Binned/All_t/*_LogScale.png histos/E_Binned/All_t/LogPlots/.");
   gSystem->Exec("mv *.png histos/.");

   std::cout << "Output file has been written" << std::endl; 

}