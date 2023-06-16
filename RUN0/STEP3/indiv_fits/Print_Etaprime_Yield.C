// Author: J. McIntyre
// 16 May 2023
// This file generates new combined histograms and prints png files of selected histograms.

#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TH2.h>
#include <math.h>
#include <TImage.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLatex.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TLegend.h>

//using namespace std;       // Removes the need for std::

void Print_Etaprime_Yield()
{
   // Gather selected hisograms
   TFile *fStep1 = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP3/indiv_fits/TwoG_may_16_2023_STEP3.root");
   TFile *fStep2 = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP3/indiv_fits/min_40MeV_Fit_width/TwoG_may_16_2023_STEP3_40MeV.root");
   TFile *fStep3 = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP3/indiv_fits/min_60MeV_Fit_width/TwoG_may_16_2023_STEP3_60MeV.root");
   TH1D *h45 = (TH1D*)fStep1->Get("Etap_yield");
   TH1D *h40 = (TH1D*)fStep2->Get("Etap_yield");
   TH1D *h60 = (TH1D*)fStep3->Get("Etap_yield");

   TString PName;

   // Add new TH1D to ROOT file for reference in later analysis work
   //TFile *yieldFile = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP3/indiv_fits/TwoG_may_16_2023_STEP3.root", "UPDATE");

   // Add new ROOT file
   //TFile *yieldFile = new TFile("TwoG_may_16_2023_STEP3_Print.root", "NEW");

   TH1D *Etap_yield_comparison = new TH1D("Etap_yield_comparison", "#eta' yield for fits with various minimum parameter width limits;x_{_{0-1 }}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }0) ,  x_{_{1-2}}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }1) ,  ... , x_{_{11-12}}= (E_{#gamma_{Bin}}=_{ }3, |t|_{_{Bin}}=_{ }2);#eta' yield per (E_{#gamma}, |t|) Bin", 60, 0, 12);
   //Etap_yieldEtap_yield_comparisonSumw2();

   TH1D *Etap_yield_comparison1 = new TH1D("Etap_yield_comparison1", "#eta' yield for fits with various minimum parameter width limits;x_{_{0-1 }}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }0) ,  x_{_{1-2}}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }1) ,  ... , x_{_{11-12}}= (E_{#gamma_{Bin}}=_{ }3, |t|_{_{Bin}}=_{ }2);#eta' yield per (E_{#gamma}, |t|) Bin", 60, 0, 12);
   //Etap_yieldEtap_yield_comparisonSumw2();

   TH1D *Etap_yield_comparison2 = new TH1D("Etap_yield_comparison2", "#eta' yield for fits with various minimum parameter width limits;x_{_{0-1 }}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }0) ,  x_{_{1-2}}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }1) ,  ... , x_{_{11-12}}= (E_{#gamma_{Bin}}=_{ }3, |t|_{_{Bin}}=_{ }2);#eta' yield per (E_{#gamma}, |t|) Bin", 60, 0, 12);
   //Etap_yieldEtap_yield_comparisonSumw2();

   TH1D *Etap_yield_comparison3 = new TH1D("Etap_yield_comparison3", "#eta' yield for fits with various minimum parameter width limits;x_{_{0-1 }}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }0) ,  x_{_{1-2}}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }1) ,  ... , x_{_{11-12}}= (E_{#gamma_{Bin}}=_{ }3, |t|_{_{Bin}}=_{ }2);#eta' yield per (E_{#gamma}, |t|) Bin", 60, 0, 12);
   //Etap_yieldEtap_yield_comparisonSumw2();

   for (int i = 0; i < 12; i++) {
      int j = i + 1;
      int abq = 2 + i * 5;
      int bbq = 3 + i * 5;
      int cbq = 4 + i * 5;
      double yh40 = h40->GetBinContent(j);
      double eh40 = h40->GetBinError(j);
      double yh45 = h45->GetBinContent(j);
      double eh45 = h45->GetBinError(j);
      double yh60 = h60->GetBinContent(j);
      double eh60 = h60->GetBinError(j);
      if (yh40 > 5.0){
         Etap_yield_comparison1->SetBinContent(abq, yh40);
         Etap_yield_comparison1->SetBinError(abq, eh40);
         Etap_yield_comparison->SetBinContent(abq, yh40);
         Etap_yield_comparison->SetBinError(abq, eh40);
      }
      if (yh45 > 5.0){
         Etap_yield_comparison2->SetBinContent(bbq, yh45);
         Etap_yield_comparison2->SetBinError(bbq, eh45);
         Etap_yield_comparison->SetBinContent(bbq, yh45);
         Etap_yield_comparison->SetBinError(bbq, eh45);
      }
      if (yh60 > 5.0){
         Etap_yield_comparison3->SetBinContent(cbq, yh60);
         Etap_yield_comparison3->SetBinError(cbq, eh60);
         Etap_yield_comparison->SetBinContent(cbq, yh60);
         Etap_yield_comparison->SetBinError(cbq, eh60);
      }
   }

   TCanvas *c_101 = new TCanvas("c_101", "Etaprime Yield for various fit widths", 800, 600);
   c_101->cd();
   gStyle->SetOptStat(0);
   Etap_yield_comparison1->SetLineColor(1);
   Etap_yield_comparison2->SetLineColor(2);
   Etap_yield_comparison3->SetLineColor(kBlue);
   Etap_yield_comparison1->SetMaximum(4000);
   Etap_yield_comparison1->GetXaxis()->SetTitleOffset(1.3);
   Etap_yield_comparison1->GetYaxis()->SetTitleOffset(1.2);
   Etap_yield_comparison1->Draw();
   Etap_yield_comparison2->Draw("same");
   Etap_yield_comparison3->Draw("same");

   TLegend *legnd1 = new TLegend(0.1, 0.1, 0.4, 0.25);
   legnd1->SetTextSize(0.028);
   //legnd1->AddEntry((TObject*)0, "#eta' Yield", "");
   legnd1->AddEntry(Etap_yield_comparison1, "40 MeV/c^{2} min. fit width", "l");
   legnd1->AddEntry(Etap_yield_comparison2, "45 MeV/c^{2} min. fit width", "l");
   legnd1->AddEntry(Etap_yield_comparison3, "60 MeV/c^{2} min. fit width", "l");
   c_101->SetLogy(1);
   legnd1->Draw("same");
   c_101->Update();
   PName.Form("Etap_Yields_LogScale.png");
   c_101->Print(PName);

   c_101->SetLogy(0);
   Etap_yield_comparison1->GetXaxis()->SetTitleOffset(1.3);
   Etap_yield_comparison1->GetYaxis()->SetTitleOffset(1.3);
   Etap_yield_comparison1->SetMinimum(0.0);
   Etap_yield_comparison1->Draw();
   Etap_yield_comparison2->Draw("same");
   Etap_yield_comparison3->Draw("same");

   TLegend *legnd2 = new TLegend(0.1, 0.75, 0.4, 0.9);
   legnd2->SetTextSize(0.028);
   //legnd2->AddEntry((TObject*)0, "#eta' Yield", "");
   legnd2->AddEntry(Etap_yield_comparison1, "40 MeV/c^{2} min. fit width", "l");
   legnd2->AddEntry(Etap_yield_comparison2, "45 MeV/c^{2} min. fit width", "l");
   legnd2->AddEntry(Etap_yield_comparison3, "60 MeV/c^{2} min. fit width", "l");
   legnd2->Draw("same");

   c_101->Update();
   Etap_yield_comparison2->Draw("same");
   PName.Form("Etap_Yields.png");
   c_101->Print(PName);

   //gFile = yieldFile;
   //gDirectory->WriteObject(Etap_yield_comparison, "Etap_yield_comparison");
   //yieldFile->Close();

   /*
   TCanvas *c_100 = new TCanvas("c_100", "Etaprime Yield for various fit widths", 800, 600);
   c_100->cd();
   gStyle->SetOptStat(0);
   h40->SetTitle("#eta' yield for gaussian fits with various minimum width limits");
   h40->GetXaxis()->SetTitle("x_{_{0-1 }}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }0) ,  x_{_{1-2}}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }1) ,  ... , x_{_{11-12}}= (E_{#gamma_{Bin}}=_{ }3, |t|_{_{Bin}}=_{ }2)");
   h40->GetYaxis()->SetTitle("#eta' yield per (E_{#gamma}, |t|) Bin");
   //TPaveStats *st_100 = (TPaveStats*)h40->FindObject("stats");
   h40->GetXaxis()->SetTitleOffset(1.3);
   h40->GetYaxis()->SetTitleOffset(1.2);
   c_100->Update();
   h40->SetMaximum(4000);
   //h40->GetXaxis()->SetRangeUser(0.0, 1.2);
   h45->SetLineColor(1);
   h40->SetLineColor(2);
   h60->SetLineColor(kGreen);
   h40->Draw();
   h45->Draw("same");
   h60->Draw("same");

   TLegend *legnd3 = new TLegend(0.1, 0.1, 0.4, 0.25);
   legnd3->SetTextSize(0.028);
   //legnd3->AddEntry((TObject*)0, "#eta' Yield", "");
   legnd3->AddEntry(h40, "40 MeV/c^{2} min. fit width", "l");
   legnd3->AddEntry(h45, "45 MeV/c^{2} min. fit width", "l");
   legnd3->AddEntry(h60, "60 MeV/c^{2} min. fit width", "l");
   legnd3->Draw("same");

   c_100->SetLogy(1);
   PName.Form("Etap_Yields_LogScale.png");
   c_100->Update();
   h60->Draw("same");
   c_100->Print(PName);

   c_100->SetLogy(0);
   c_100->Update();
   h40->GetXaxis()->SetTitleOffset(1.3);
   h40->GetYaxis()->SetTitleOffset(1.3);
   c_100->Update();
   h40->Draw();
   h45->Draw("same");
   h60->Draw("same");
   TLegend *legnd4 = new TLegend(0.1, 0.75, 0.4, 0.9);
   legnd4->SetTextSize(0.028);
   //legnd4->AddEntry((TObject*)0, "#eta' Yield", "");
   legnd4->AddEntry(h40, "40 MeV/c^{2} min. fit width", "l");
   legnd4->AddEntry(h45, "45 MeV/c^{2} min. fit width", "l");
   legnd4->AddEntry(h60, "60 MeV/c^{2} min. fit width", "l");
   legnd4->Draw("same");
   h40->SetMinimum(0.0);
   c_100->Update();
   h60->Draw("same");
   PName.Form("Etap_Yields.png");
   c_100->Print(PName);
   */

   //delete c_100;
   delete c_101;

   std::cout << "Output file has been written \n";


   /*
   // Organize folder info
   gSystem->Exec("mkdir -p histos/Cuts");

   gSystem->Exec("mv M2g_hcut* histos/Cuts/.");

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
   */
}
