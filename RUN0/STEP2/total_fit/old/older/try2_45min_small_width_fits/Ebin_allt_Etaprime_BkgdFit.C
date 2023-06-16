///////////////////////////////////////////////////
// Author: J. McIntyre
// 16 May 2023
///////////////////////////////////////////////////
// This file is used for fitting the overall |t| 
// (within the range of interest) histogram with 
// basic cleanup cuts and tag weighting applied.
// The non-binned (in |t|) histogram is used to 
// determine the background model that will be used
// to fit the background in each individual |t| binned
// histogram. A factor coefficient is added to the model
// for scaling based on each bin's value.
//
// The background is modeled as a pol3.
// 
// To change this file modify the following:
// ------------------------------------------------
// This files name     --> Ebin_allt_Etaprime_BkgdFit
// Real data root file --> TwoG_may_16_2023Analysis.root
// Real data histogram --> M_2gamma
// ------------------------------------------------

/*
   ////////////////////////////////////////////////
   TLine *zeroline3 = new TLine(0.74, 0, 1.18, 0);
   zeroline3->SetLineColor(1);
   zeroline3->SetLineStyle(10);
   zeroline3->SetLineWidth(2);
   zeroline3->Draw("hist");
   ////////////////////////////////////////////////
*/

// ------------------------------------------------
// DRAWING OPTIONS
//----------------
// gStyle->SetOptStat(0);   // No Legend
// gStyle->SetOptTitle(0);   // No title
// histo->SetLineWidth(5);   // Line thickness
// histo->SetLineStyle(5);   // Line type (1=solid)
// histo->SetLineColor(kBlue);   // Line color
// histo->SetFillColor(kBlue);   // Filling color to histo (38 = light blue)
// histo->SetFillColorAlpha(kBlue,0.25);   // Suppose to be semi-transparent
// ------------------------------------------------

// Regarding #include, using <> tells the compiler to search for the .h file in the system include directories,
// while using "" has it search for the .h file in the current directory and then in the system include directories.
#include <iostream>
#include <TF1.h>
#include <cmath>
#include <TLegend.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>
#include <Math/MinimizerOptions.h>
#include <Math/Minimizer.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <fstream>
#include <math.h>
#include <TMath.h>
#include <string>
#include <TString.h>
#include <stdlib.h>
#include <Math/Functor.h>
#include <TSystem.h>
#include <TROOT.h>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <stdio.h>
#include <TPad.h>
#include <TText.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include <TAxis.h>
#include <TObjString.h>
#include <TList.h>
#include <Math/RootFinderAlgorithms.h>
#include <Math/RootFinder.h>
#include <Math/GSLMultiRootFinder.h>
#include <Math/WrappedMultiTF1.h>
#include <TKey.h>
#include <TGraphPainter.h>
#include <TAttFill.h>

/*
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TH2.h>
#include <TF2.h>
#include <TGraphSmooth.h>
#include <TSpline.h>
#include <TArrow.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TApplication.h>
#include <TRint.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TMarker.h>
#include "TRandom3.h"
*/

// using namespace std;

#define PI 3.14159265359
#define sigmaRANGE 2.25

TH1D *remass;
TH1D *Secondhisto;
TH1D *Addhisto;
TH1D *newt;
TF1  *fsig1;
TF1  *fsig2;
TF1  *fsig3;
TF1  *fsig4;
TF1  *fback;
TF1  *ftot;
TH1D *SG_Pi;
TH1D *SG_Eta;
TH1D *SG_Omega;
TF1  *fPi_SG; 
TF1  *fEta_SG; 
TF1  *fOmega_SG;

/////////////////////////////////////////
// User defined equations for the fits //
/////////////////////////////////////////

// Pi_0 Fit
double fsignal1(double *x, double *par) {
   double xx = x[0];
   double gaus1 = pow(par[0], 2) * exp(-0.5 * pow((xx - par[1])/par[2], 2)) +
                  pow(par[3], 2) * exp(-0.5 * pow((xx - par[4])/par[5], 2));
   return gaus1;
}

// Eta Fit
double fsignal2(double *x, double *par) {
   double xx = x[0];
   double gaus2 = pow(par[0], 2) * exp(-0.5 * pow((xx - par[1])/par[2], 2)) +
                  pow(par[3], 2) * exp(-0.5 * pow((xx - par[4])/par[5], 2));

   return gaus2;
}

// Leakage Fit
double fsignal3(double *x, double *par) {
   double xx = x[0];
   double gaus3 = pow(par[0], 2) * exp(-0.5 * pow((xx - par[1])/par[2], 2)) +
                  pow(par[3], 2) * exp(-0.5 * pow((xx - par[4])/par[5], 2));
   return gaus3;
}

// Etap Fit, not enough data for a double gaussian to work well
double fsignal4(double *x, double *par) {
   double xx = x[0];
   double gaus4 = pow(par[0], 2) * exp(-0.5 * pow((xx - par[1])/par[2], 2));
   return gaus4;
}

// Background Fit (pol3 * Scaler)
double fbackground(double *x, double *par) {
   double xx = x[0];
   double poly3 = (par[0] + par[1]*xx + par[2]*pow(xx, 2) + par[3]*pow(xx, 3)) * (sqrt(pow(par[4], 2)));
   return poly3;
}


double ftotal(double *x, double *par) {
  return  fsignal1(x, par) + fsignal2(x, &par[6]) + fsignal3(x, &par[12]) + fsignal4(x, &par[18]) + fbackground(x, &par[21]);
  }


void  Ebin_allt_Etaprime_BkgdFit() {

   // Root file with new histograms for reference in later analysis work
  	TFile *MyFile = new TFile("TwoG_may_16_2023_STEP2.root", "NEW");


   // Lower edge of |t| bins
	TH1D *hBinDW = new TH1D("hBinDW", "#eta' tbin[Bin No.][0] value;Bin No.;|t| value [GeV^2/c^{4}]", 3, 0, 3);
   hBinDW->Sumw2();

   // Upper edge of |t| bins
	TH1D *hBinUP = new TH1D("hBinUP", "#eta' tbin[Bin No.][1] value;Bin No.;|t| value [GeV^2/c^{4}]", 3, 0, 3);
   hBinUP->Sumw2();

   // Selecting the number of tagger energy bins
	TH1D *EChanBINS = new TH1D("EChanBINS", "Number of Tagger Energy Bins;;Photon Energy Bins", 1, 0.0, 1.0);
   EChanBINS->Sumw2();

   //////////////////////////////////////////////////////////////////////////
   // CHANGE BIN VALUE FOR TH1D BASED ON NUMBER OF TAGGER ENERGY BINS USED //
   //////////////////////////////////////////////////////////////////////////
	TH1D *EChanRange = new TH1D("EChanRange", "Tagger channels used for Photon energy bins;E_{ph} Bin = int[(X-value / 2) - 0.5];Tagger Channel", 8, 0.0, 8.0);
   EChanRange->Sumw2();

   //////////////////////////////////////////////////////////////////////////
   // CHANGE BIN VALUE FOR TH1D BASED ON NUMBER OF TAGGER ENERGY BINS USED //
   //////////////////////////////////////////////////////////////////////////
	TH1D *EChanValue = new TH1D("EChanValue", "Photon Energy Bins;E_{ph} Bin = int[(X-value / 2) - 0.5];Photon Energy [GeV]", 8, 0.0, 8.0);
   EChanValue->Sumw2();

   ////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////
   // ALL SUBSEQUENT ANALYSIS TSELECTORS WILL REFERENCE THE HISTOGRAMS MADE  //
   // HERE FOR WHAT |t| BINS TO USE.                                         //
   ////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////
   double BinM2g[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

	// Two histos with |t| bin ranges
	for (int c1 = 0; c1 < 3; c1++) {
		int xbin = c1 + 1;
		hBinDW->SetBinContent(xbin, BinM2g[c1][0]);
		hBinUP->SetBinContent(xbin, BinM2g[c1][1]);
   }

	gFile = MyFile;
	gDirectory->WriteObject(hBinDW, "hBinDW");
	gDirectory->WriteObject(hBinUP, "hBinUP");

   ///////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////
   // BELOW IS WHERE YOU CHANGE THE TAGGER ENERGY BINS USED IN ALL FURTHER ANALYSIS //
   ///////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////

   // High energy end of bin in GeV
   // Echannel# energy value is border between tagger channel # & (# - 1)
   double Echannel0  = 5.385;
   double Echannel1  = 5.335;
   double Echannel2  = 5.280;
   double Echannel3  = 5.220;
   double Echannel4  = 5.165;
   double Echannel5  = 5.115;
   double Echannel6  = 5.065;
   double Echannel7  = 5.015;
   double Echannel8  = 4.965;
   double Echannel9  = 4.915;
   double Echannel10 = 4.860;
   double Echannel11 = 4.800;
   double Echannel12 = 4.740;
   double Echannel13 = 4.675;
   double Echannel14 = 4.615;
   double Echannel15 = 4.570;
   double Echannel16 = 4.530;
   double Echannel17 = 4.485;
   double Echannel18 = 4.435;
   double Echannel19 = 4.385;

   //////////////////////////////////////////////
   // CHANGE NUMBER OF PHOTON ENERGY BINS HERE //
   //////////////////////////////////////////////
   int NumEChan = 4;
   // Setting Tagger Energy Range
   EChanBINS->SetBinContent(1, NumEChan);

	gFile = MyFile;
	gDirectory->WriteObject(EChanBINS, "EChanBINS");

   // Tagger Energy Bin E0 = channel 0 - channel 4
   double Range0a = Echannel0;
   double Range0b = Echannel5;
   int TC_Range0a = 0;
   int TC_Range0b = 4;

   // Tagger Energy Bin E1 = channel 5 - channel 9
   double Range1a = Echannel5;
   double Range1b = Echannel10;
   int TC_Range1a = 5;
   int TC_Range1b = 9;

   // Tagger Energy Bin E2 = channel 10 - channel 13
   double Range2a = Echannel10;
   double Range2b = Echannel14;
   int TC_Range2a = 10;
   int TC_Range2b = 13;

   // Tagger Energy Bin E0 = channel 14 - channel 18
   double Range3a = Echannel14;
   double Range3b = Echannel19;
   int TC_Range3a = 14;
   int TC_Range3b = 19; // One channel beyond tagger range

   // Tagger channels used for photon energy binning
   EChanRange->SetBinContent(1, TC_Range0a);
   EChanRange->SetBinContent(2, TC_Range0b);
   EChanRange->SetBinContent(3, TC_Range1a);
   EChanRange->SetBinContent(4, TC_Range1b);
   EChanRange->SetBinContent(5, TC_Range2a);
   EChanRange->SetBinContent(6, TC_Range2b);
   EChanRange->SetBinContent(7, TC_Range3a);
   EChanRange->SetBinContent(8, TC_Range3b);

	gFile = MyFile;
	gDirectory->WriteObject(EChanRange, "EChanRange");

   // Photon energy binning
   EChanValue->SetBinContent(1, Range0a);
   EChanValue->SetBinContent(2, Range0b);
   EChanValue->SetBinContent(3, Range1a);
   EChanValue->SetBinContent(4, Range1b);
   EChanValue->SetBinContent(5, Range2a);
   EChanValue->SetBinContent(6, Range2b);
   EChanValue->SetBinContent(7, Range3a);
   EChanValue->SetBinContent(8, Range3b);

	gFile = MyFile;
	gDirectory->WriteObject(EChanValue, "EChanValue");
   MyFile->Close();

   // Real data root file
   TFile *rootf = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP1/inv_mass/TwoG_may_16_2023Analysis.root");

   // Generate TH1D that will store the mass values to use for sideband subtraction regions
   TString NEWName;
	NEWName.Form("Mass Sideband Subtraction Range for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range0b, Range0a);
   TH1D *SBRANGE = new TH1D("SBRANGE", NEWName, 6, 0.0, 6.0);
   //SBRANGE->Sumw2();

   NEWName.Form("Mass Sideband Subtraction Range for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range1b, Range1a);
   TH1D *SBRANGE2 = new TH1D("SBRANGE2", NEWName, 6, 0.0, 6.0);
   //SBRANGE2->Sumw2();

   NEWName.Form("Mass Sideband Subtraction Range for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range2b, Range2a);
   TH1D *SBRANGE3 = new TH1D("SBRANGE3", NEWName, 6, 0.0, 6.0);
   //SBRANGE3->Sumw2();

   NEWName.Form("Mass Sideband Subtraction Range for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range3b, Range3a);
   TH1D *SBRANGE4 = new TH1D("SBRANGE4", NEWName, 6, 0.0, 6.0);
   //SBRANGE4->Sumw2();

   // TCanvas setup for capturing png files
   TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
   c1->Divide(2,2);
   TCanvas *c2 = new TCanvas("c2", "c2", 1200, 900);
   c2->Divide(2,2);
   TCanvas *c3 = new TCanvas("c3", "c3", 1200, 900);
   c3->Divide(2,2);
   TCanvas *c10 = new TCanvas("c10", "c10", 800, 600);
   TCanvas *c11 = new TCanvas("c11", "c11", 800, 600);
   TCanvas *c123 = new TCanvas("c123", "c123", 800, 600);
   TCanvas *c321 = new TCanvas("c321", "c321", 800, 600);
   TCanvas *c555 = new TCanvas("c555", "c555", 800, 600);

   // Variables used to set fit limits.
   // Defined initially outside the loop,
   // set to zero within the loop at the start of each iteration.
   double uparlim1 = 0.0;
   double dparlim1 = 0.0;
   double uparlim4 = 0.0;
   double dparlim4 = 0.0;
   double uparlim7 = 0.0;
   double dparlim7 = 0.0;
   double uparlim10 = 0.0;
   double dparlim10 = 0.0;
   double uparlim13 = 0.0;
   double dparlim13 = 0.0;
   double uparlim16 = 0.0;
   double dparlim16 = 0.0;
   double uparlim19 = 0.0;
   double dparlim19 = 0.0;
   double uparlim2 = 0.0;
   double dparlim2 = 0.0;
   double uparlim5 = 0.0;
   double dparlim5 = 0.0;
   double uparlim8 = 0.0;
   double dparlim8 = 0.0;
   double uparlim11 = 0.0;
   double dparlim11 = 0.0;
   double uparlim14 = 0.0;
   double dparlim14 = 0.0;
   double uparlim17 = 0.0;
   double dparlim17 = 0.0;
   double uparlim20 = 0.0;
   double dparlim20 = 0.0;

   // Declaring variable for fit convergence (outside loop),
   // so that TCanvases aren't deleted if a fit fails.
   // Undeleted histos on screen will allows easier troubleshooting.   
   int conv_fit = 1;

   TString htitle;
   int totalhisto;
   for (int Echannel = 0; Echannel < NumEChan; Echannel++) {
      // Left from previous version which included addind |t| binned histos together.
      // For naming png's. Can be replaced with Echannel.
      totalhisto = Echannel; 

      // Set variable, for fit limits, to zero
      uparlim1 = 0.0;
      dparlim1 = 0.0;
      uparlim4 = 0.0;
      dparlim4 = 0.0;
      uparlim7 = 0.0;
      dparlim7 = 0.0;
      uparlim10 = 0.0;
      dparlim10 = 0.0;
      uparlim13 = 0.0;
      dparlim13 = 0.0;
      uparlim16 = 0.0;
      dparlim16 = 0.0;
      uparlim19 = 0.0;
      dparlim19 = 0.0;
      uparlim2 = 0.0;
      dparlim2 = 0.0;
      uparlim5 = 0.0;
      dparlim5 = 0.0;
      uparlim8 = 0.0;
      dparlim8 = 0.0;
      uparlim11 = 0.0;
      dparlim11 = 0.0;
      uparlim14 = 0.0;
      dparlim14 = 0.0;
      uparlim17 = 0.0;
      dparlim17 = 0.0;
      uparlim20 = 0.0;
      dparlim20 = 0.0;

      htitle.Form("Allt_E%d", Echannel);
      remass = (TH1D*)rootf->Get(htitle);

		// Bin width for real data
		double binw = remass->GetBinWidth(1);  // In GeV/c^2
      int binwMeV = binw * 1000;             // In MeV/c^2

		// Histograms maximum heights
		double norm = remass->GetMaximum();

		// Cloning MassBinned histogram for use in different TPads
		TH1D *newhfinal2_1  = (TH1D*)remass->Clone("newhfinal2_1");
		TH1D *newhfinal2_2  = (TH1D*)remass->Clone("newhfinal2_2");
		TH1D *newhfinal2_3  = (TH1D*)remass->Clone("newhfinal2_3");
		TH1D *newhfinal2_4  = (TH1D*)remass->Clone("newhfinal2_4");
		TH1D *newhfinal2_1a = (TH1D*)remass->Clone("newhfinal2_1a");
		TH1D *newhfinal2_2a = (TH1D*)remass->Clone("newhfinal2_2a");
		TH1D *newhfinal2_3a = (TH1D*)remass->Clone("newhfinal2_3a");
		TH1D *newhfinal2_1b = (TH1D*)remass->Clone("newhfinal2_1b");
		TH1D *newhfinal2_1c = (TH1D*)remass->Clone("newhfinal2_1c");
		TH1D *newhfinal2_2b = (TH1D*)remass->Clone("newhfinal2_2b");
		TH1D *newhfinal2_2c = (TH1D*)remass->Clone("newhfinal2_2c");
		TH1D *newhfinal2_3b = (TH1D*)remass->Clone("newhfinal2_3b");
		TH1D *newhfinal2_3c = (TH1D*)remass->Clone("newhfinal2_3c");
		TH1D *newhfinal2_5  = (TH1D*)remass->Clone("newhfinal2_5");
		TH1D *newhfinal2_5a = (TH1D*)remass->Clone("newhfinal2_5a");
		TH1D *newhfinal2_5b = (TH1D*)remass->Clone("newhfinal2_5b");
		TH1D *newhfinal2_6  = (TH1D*)remass->Clone("newhfinal2_6");
		TH1D *newhfinal2_6b = (TH1D*)remass->Clone("newhfinal2_6b");
		TH1D *newhfinal2_6c = (TH1D*)remass->Clone("newhfinal2_6c");
		TH1D *newhfinal2_7  = (TH1D*)remass->Clone("newhfinal2_7");
		TH1D *newhfinal2_8  = (TH1D*)remass->Clone("newhfinal2_8");
		TH1D *newhfinal2_9  = (TH1D*)remass->Clone("newhfinal2_9");
		TH1D *newhfinal2_123 = (TH1D*)remass->Clone("newhfinal2_123");
		TH1D *newhfinal2_321 = (TH1D*)remass->Clone("newhfinal2_321");
		TH1D *newhfinal2_555 = (TH1D*)remass->Clone("newhfinal2_555");

      // Clones for filling a plot with colored sections
		TH1D *hfillClone1 = (TH1D*)remass->Clone("hfillClone1");
		TH1D *hfillClone2 = (TH1D*)remass->Clone("hfillClone2");
		TH1D *hfillClone3 = (TH1D*)remass->Clone("hfillClone3");
		TH1D *hfillClone4 = (TH1D*)remass->Clone("hfillClone4");
		TH1D *hfillClone5 = (TH1D*)remass->Clone("hfillClone5");
		TH1D *hfillClone6 = (TH1D*)remass->Clone("hfillClone6");
		TH1D *hfillClone7 = (TH1D*)remass->Clone("hfillClone7");
		TH1D *hfillClone8 = (TH1D*)remass->Clone("hfillClone8");
		TH1D *hfillClone9 = (TH1D*)remass->Clone("hfillClone9");

      // Used to fit a Single Gaussian in the range of the Pi0, Eta, and Omega
      // These parameters will help with the fitting as the invariant mass walks in |t|
      SG_Pi     = (TH1D*)remass->Clone("SG_Pi");
      SG_Eta    = (TH1D*)remass->Clone("SG_Eta");
      SG_Omega  = (TH1D*)remass->Clone("SG_Omega");
      fPi_SG    = new TF1("fPi_SG", "gaus", 0.1, 0.2); 
      fEta_SG   = new TF1("fEta_SG", "gaus", 0.48, 0.6); 
      fOmega_SG = new TF1("fOmega_SG", "gaus", 0.7, 0.8); 
      double PiPar[3];
      double EtaPar[3];
      double OmegaPar[3];
      SG_Pi->Fit(fPi_SG, "R");
      SG_Eta->Fit(fEta_SG, "R");
      SG_Omega->Fit(fOmega_SG, "R");
      fPi_SG->GetParameters(&PiPar[0]);
      fEta_SG->GetParameters(&EtaPar[0]);
      fOmega_SG->GetParameters(&OmegaPar[0]);

		// Signals of Interest & the 3 gamma leakage
		fsig1 = new TF1("fsig1", fsignal1, 0.0, 0.3, 6);
		fsig2 = new TF1("fsig2", fsignal2, 0.3, 0.8, 6);
		fsig3 = new TF1("fsig3", fsignal3, 0.32, 1.15, 6);
		fsig4 = new TF1("fsig4", fsignal4, 0.70, 1.20, 3);

		// Background. The 3 gamma leakage (a.k.a. background) is kept separate (fsig3).
		fback = new TF1("fback", fbackground, 0.0, 1.28, 5);   // 3rd Order Polynomial

		// Total Signal
		ftot = new TF1("ftot", ftotal, 0.0, 1.28, 26);

		// Starting parameters for the fits
      /* Manual Input (estimate from real data histo using Excel or the like) */
		double BkgdPar[6];
		BkgdPar[0] = 1.0;
		BkgdPar[1] = 1.0;
		BkgdPar[2] = 1.0;
		BkgdPar[3] = 1.0;
		BkgdPar[4] = 1.0;     // In individual fits (next analysis step) this parameter will be a scaling factor 

		double par[26] = {sqrt(norm * 0.639753), 0.1369, 0.0148, 
								sqrt(norm * 0.373938), 0.1389, 0.0256, 
								sqrt(norm * 0.0680778), 0.5412, 0.0314, 
								sqrt(norm * 0.0225913), 0.555, 0.046, 
								sqrt(norm * 0.00461673), 0.744, 0.054, 
								sqrt(norm * 0.00433187), 0.721, 0.072, 
								sqrt(norm * 0.00119216), 0.9574, 0.0514,
                        BkgdPar[0], BkgdPar[1], BkgdPar[2], BkgdPar[3], BkgdPar[4]};

		ftot->SetParameters(par);

		// Helps the fits not to run off
		for (int i = 0; i < 26; i++) {
			ftot->FixParameter(i, par[i]);
		}

		remass->Fit("ftot", "RB");

      // Do not release par[25]. 
      // This is a scaling factor for use in individual |t| binned fits (next step).
		for (int i = 0; i < 25; i++) {
			ftot->ReleaseParameter(i);
		}

      // Used to prevent STATUS=CALL LIMIT. Gives fitter more tries.
      ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

      if (Echannel == 2) {
         // Signal #1 - Pi0 mean limits for the double gaussian
         // Gaussian 1
         uparlim1 = PiPar[1] + 0.040;
         dparlim1 = PiPar[1] - 0.020;
         // Gaussian 2
         uparlim4 = PiPar[1] + 0.040;
         dparlim4 = PiPar[1] - 0.020;

         // Signal #1 - Pi0 std dev limits for the double gaussian
         // Gaussian 1
         uparlim2 = 0.085;
         dparlim2 = 0.010;
         // Gaussian 2
         uparlim5 = 0.085;
         dparlim5 = 0.010;
      }
      else {
         // Signal #1 - Pi0 mean limits for the double gaussian
         // Gaussian 1
         uparlim1 = PiPar[1] + 0.030;
         dparlim1 = PiPar[1] - 0.020;
         // Gaussian 2
         uparlim4 = PiPar[1] + 0.030;
         dparlim4 = PiPar[1] - 0.020;

         // Signal #1 - Pi0 std dev limits for the double gaussian
         // Gaussian 1
         uparlim2 = 0.075;
         dparlim2 = 0.010;
         // Gaussian 2
         uparlim5 = 0.075;
         dparlim5 = 0.010;
      }

      // Signal #2 - Eta mean limits for the double gaussian
      // Gaussian 1
      uparlim7 = EtaPar[1] + 0.045;
      dparlim7 = EtaPar[1] - 0.010;
      // Gaussian 2
      uparlim10 = EtaPar[1] + 0.045;
      dparlim10 = EtaPar[1] - 0.010;

      // Signal #2 - Eta std dev limits for the double gaussian
      // Gaussian 1
      uparlim8 = 0.080;
      dparlim8 = 0.010;
      // Gaussian 2
      uparlim11 = 0.080;
      dparlim11 = 0.010;

      if (Echannel == 2) {
         // Signal #3 - Omega Leakage mean limits for the double gaussian
         // Gaussian 1
         uparlim13 = 0.755;
         dparlim13 = 0.730;
         // Gaussian 2
         uparlim16 = 0.755;
         dparlim16 = 0.730;
      }
      else {
         // Signal #3 - Omega Leakage mean limits for the double gaussian
         // Gaussian 1
         uparlim13 = 0.745;
         dparlim13 = 0.730;
         // Gaussian 2
         uparlim16 = 0.745;
         dparlim16 = 0.730;
      }

      // Signal #3 - Omega Leakage std dev limits for the double gaussian
      // Gaussian 1
      uparlim14 = 0.075;
      dparlim14 = 0.035;
      // Gaussian 2
      uparlim17 = 0.075;
      dparlim17 = 0.035;

      // Signal #4 - Etap mean limits for the single gaussian
      // Single Gaussian
      uparlim19 = 0.990;
      dparlim19 = 0.940;

      // Signal #4 - Etap std dev limit for the single gaussian
      // Single Gaussian
      uparlim20 = 0.080;
      dparlim20 = 0.045;

      ///////////////////////////////////////
      /* Setting up the limits for the fit */
      ///////////////////////////////////////
      // Signal #1 - Pi0 mean limits for the double gaussian
      ftot->SetParLimits(1, dparlim1, uparlim1);
      ftot->SetParLimits(4, dparlim4, uparlim4);

      // Signal #2 - Eta mean limits for the double gaussian
      ftot->SetParLimits(7, dparlim7, uparlim7);
      ftot->SetParLimits(10, dparlim10, uparlim10);

      // Signal #3 - Omega Leakage mean limits for the double gaussian
      ftot->SetParLimits(13, dparlim13, uparlim13);
      ftot->SetParLimits(16, dparlim16, uparlim16);

      // Signal #4 - Etap mean limits for the single gaussian
      ftot->SetParLimits(19, dparlim19, uparlim19);

      // Signal #1 - Pi0 std dev limits for the double gaussian
      ftot->SetParLimits(2, dparlim2, uparlim2);
      ftot->SetParLimits(5, dparlim5, uparlim5);

      // Signal #2 - Eta std dev limits for the double gaussian
      ftot->SetParLimits(8, dparlim8, uparlim8);
      ftot->SetParLimits(11, dparlim11, uparlim11);

      // Signal #3 - Omega Leakage std dev limits for the double gaussian
      ftot->SetParLimits(14, dparlim14, uparlim14);
      ftot->SetParLimits(17, dparlim17, uparlim17);

      // Signal #4 - Etap std dev limit for the single gaussian
      ftot->SetParLimits(20, dparlim20, uparlim20);

		// Set Parameter Names
		ftot->SetParName(0,  "Pi0_____peak_1");
		ftot->SetParName(1,  "Pi0_____mean_1");
		ftot->SetParName(2,  "Pi0____sigma_1");
		ftot->SetParName(3,  "Pi0_____peak_2");
		ftot->SetParName(4,  "Pi0_____mean_2");
		ftot->SetParName(5,  "Pi0____sigma_2");

		ftot->SetParName(6,  "Eta_____peak_1");
		ftot->SetParName(7,  "Eta_____mean_1");
		ftot->SetParName(8,  "Eta____sigma_1");
		ftot->SetParName(9,  "Eta_____peak_2");
		ftot->SetParName(10, "Eta_____mean_2");
		ftot->SetParName(11, "Eta____sigma_2");

		ftot->SetParName(12, "3G_leak_peak_1");
		ftot->SetParName(13, "3G_leak_mean_1");
		ftot->SetParName(14, "3G_leak_sigma1");
		ftot->SetParName(15, "3G_leak_peak_2");
		ftot->SetParName(16, "3G_leak_mean_2");
		ftot->SetParName(17, "3G_leak_sigma2");

		ftot->SetParName(18, "Etap____peak_1");
		ftot->SetParName(19, "Etap____mean_1");
		ftot->SetParName(20, "Etap___sigma_1");
	 
		ftot->SetParName(21, "p0____________");
		ftot->SetParName(22, "p1____________");
		ftot->SetParName(23, "p2____________");
		ftot->SetParName(24, "p3____________");
		ftot->SetParName(25, "Bkgd Scaler___");

		////////////////////////////////////////////////////
		// Finding the actual fit and setting line colors //
		// Added "newpar" since I thought I might need it //
		////////////////////////////////////////////////////
      auto fitResult = remass->Fit("ftot", "SBRM");

      auto covMatrix = fitResult->GetCovarianceMatrix();
      std::cout << "Covariance matrix from the fit ";
      //covMatrix.Print();
      int conv_fit = fitResult->IsValid();
      double *p = ftot->GetParameters();

      std::cout << "**********************" << std::endl;
      std::cout << conv_fit << std::endl;
      std::cout << "**********************" << std::endl;

      // Passing in updated fit parameters
		double newpar[26];
		ftot->GetParameters(newpar);
		ftot->SetLineColor(kRed);
		fsig1->SetParameters(newpar);
		fsig1->SetLineColor(kBlue);
		fsig2->SetParameters(&newpar[6]);
		fsig2->SetLineColor(kBlue);
		fsig3->SetParameters(&newpar[12]);
		fsig3->SetLineColor(kMagenta);
		fsig4->SetParameters(&newpar[18]);
		fsig4->SetLineColor(1);
		fback->SetParameters(&newpar[21]);
		fback->SetLineColor(kGreen);

      // Location of the fits' maximum
		double GausMax1 = fsig1->GetMaximumX();	
		double GausMax2 = fsig2->GetMaximumX();
      double GausMax3 = fsig3->GetMaximumX();
		double GausMax4 = fsig4->GetMaximumX();

		///////////////////////////////////////////
		//////////////// Real Data ////////////////
		///////////////////////////////////////////
		double peak_sig1   = ftot->GetParameter(0);    // Pi0
		double peak2_sig1  = ftot->GetParameter(3);
		double peak_sig2   = ftot->GetParameter(6);    // Eta
		double peak2_sig2  = ftot->GetParameter(9);

		double peak_sig3   = ftot->GetParameter(12);   // Leakage "signal"
		double peak2_sig3  = ftot->GetParameter(15);   // Leakage "signal"

		double peak_sig4   = ftot->GetParameter(18);   // Eta prime

		double mean_sig1   = ftot->GetParameter(1);    // Pi0
		double mean2_sig1  = ftot->GetParameter(4);
		double mean_sig2   = ftot->GetParameter(7);    // Eta
		double mean2_sig2  = ftot->GetParameter(10);

		double mean_sig3   = ftot->GetParameter(13);   // Leakage "signal"
		double mean2_sig3  = ftot->GetParameter(16);   // Leakage "signal"

		double mean_sig4   = ftot->GetParameter(19);   // Eta prime

		double sigma_sig1  = ftot->GetParameter(2);    // Pi0
		double sigma2_sig1 = ftot->GetParameter(5);
		double sigma_sig2  = ftot->GetParameter(8);    // Eta
		double sigma2_sig2 = ftot->GetParameter(11);

		double sigma_sig3  = ftot->GetParameter(14);   // Leakage "signal"
		double sigma2_sig3 = ftot->GetParameter(17);   // Leakage "signal"

		double sigma_sig4  = ftot->GetParameter(20);   // Eta prime

      // Weighted average of the fitting double gaussians' mean & sigma
      double w1AvgPi = peak_sig1  / (peak_sig1 + peak2_sig1);
      double w2AvgPi = peak2_sig1 / (peak_sig1 + peak2_sig1);
      double mu_Pi = (w1AvgPi * mean_sig1) + (w2AvgPi * mean2_sig1);
      double sigma_Pi = sqrt(w1AvgPi * (pow(sigma_sig1, 2.0) + pow((mean_sig1 - mu_Pi), 2.0)) + w2AvgPi * (pow(sigma2_sig1, 2.0) + pow((mean2_sig1 - mu_Pi), 2.0)));

      double w1AvgEta = peak_sig2  / (peak_sig2 + peak2_sig2);
      double w2AvgEta = peak2_sig2 / (peak_sig2 + peak2_sig2);
      double mu_Eta = (w1AvgEta * mean_sig2) + (w2AvgEta * mean2_sig2);
      double sigma_Eta = sqrt(w1AvgEta * (pow(sigma_sig2, 2.0) + pow((mean_sig2 - mu_Eta), 2.0)) + w2AvgEta * (pow(sigma2_sig2, 2.0) + pow((mean2_sig2 - mu_Eta), 2.0)));

      double w1AvgOmega = peak_sig3  / (peak_sig3 + peak2_sig3);
      double w2AvgOmega = peak2_sig3 / (peak_sig3 + peak2_sig3);
      double mu_Omega = (w1AvgOmega * mean_sig3) + (w2AvgOmega * mean2_sig3);
      double sigma_Omega = sqrt(w1AvgOmega * (pow(sigma_sig3, 2.0) + pow((mean_sig3 - mu_Omega), 2.0)) + w2AvgOmega * (pow(sigma2_sig3, 2.0) + pow((mean2_sig3 - mu_Omega), 2.0)));

      double TOTALPi   = fsig1->Integral(0.05, 0.30) / binw;
      double TOTALEta  = fsig2->Integral(0.30, 0.80) / binw;
      double TOTALEtap = fsig4->Integral(0.70, 1.20) / binw;

      // Calculating the ranges for accidental subtraction weights // 
      // Center of Neutral Pion Signal Region
		double mx1 = mu_Pi;     // Weighted average of double gaussian fit
      /*
      // Alternate center for signal region
		double mx1 = GausMax1;  // Maximum of the fit fsig1
      */

      // sigmawidth = half of signal region rounded to nearest MeV
      // (1000[MeV/GeV conversion] * selected number of sigma width * weighted average std. dev. of double gaussian fit)
      int sigmawidth = lround(1000 * sigmaRANGE * sigma_Pi);

      // See if it is close to the next bin. If so add a bin.
      int AddBin = 0;
      if (sigmawidth % binwMeV >= 3) {
         AddBin = 1;
      }

      // Get the number of bins equivalent to half the signal region width
      sigmawidth /= binwMeV;
      int sigmaBin = sigmawidth + AddBin;

      // Bin value for the double gaussian fit weighted average mean location
      int pimean = remass->FindBin(mx1);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*                                    Invariant Mass Increases -->                                    */
      /* SB1 Region = (PiLower, PiLow) || Signal Region = (Pidown, Piup) || SB2 Region = (PiHigh, PiHigher) */
      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      // X-values (GeV/c^2) to use in integrals of a TF1 fit (don't forget to divide by bin width of real histo)
      // This makes the signal region 1 bin larger than needed
		double Piup   = remass->GetXaxis()->GetBinUpEdge(pimean + sigmaBin);
		double Pidown = remass->GetXaxis()->GetBinLowEdge(pimean - sigmaBin);

      // Far end of SB Regions
      // Make SB1 Region width = sigmaBin
		double PiLower = remass->GetXaxis()->GetBinLowEdge(pimean - 2 * sigmaBin - 1);
      bool HitLowerBound = false;
      int AddSBbin = 0;
      if (PiLower <= 0.060) {
         PiLower = 0.060;
         HitLowerBound = true;
         int LB = remass->FindBin(0.060);
         // No. of Bin difference btw 60 MeV/c^2 Bin and Lower Bound Bin
         // How many bins needed to be added to SB2, so that width of SB1 + SB2 = Signal Region
         AddSBbin = LB - (pimean - 2 * sigmaBin - 1);
      }
      // Since, Signal Region = (2 * sigmaBin) + 1 bin (central bin)
      // Make SB2 Region width = sigmaBin + 1 Bin
      // This way width of SB1 + SB2 = Signal Region
      // Remember there is a 1 Bin gap between Sideband and Signal Regions
      // AddSBbin is added in case SB1 is narrower due to the cutoff at 50 MeV/c^2
		double PiHigher = remass->GetXaxis()->GetBinUpEdge(pimean + 2 * sigmaBin + 2 + AddSBbin);

      // Inner Edge of SB Regions
		double PiHigh = remass->GetXaxis()->GetBinUpEdge(pimean + sigmaBin + 1);  // So that signal and sideband regions do not overlap
		double PiLow  = remass->GetXaxis()->GetBinLowEdge(pimean - sigmaBin - 1); // So that signal and sideband regions do not overlap

      std::cout << std::endl;
      if (HitLowerBound) {
         std::cout << "Lower Bound of 50 MeV/c^2 was reached" << std::endl;
      }
      std::cout << "PiLower = " << PiLower * 1000 << std::endl;
      std::cout << "PiLow = " << PiLow * 1000 << std::endl;
      std::cout << "Pidown = " << Pidown * 1000 << std::endl;
      std::cout << "Piup = " << Piup * 1000 << std::endl;
      std::cout << "PiHigh = " << PiHigh * 1000 << std::endl;
      std::cout << "PiHigher = " << PiHigher * 1000 << std::endl;
      std::cout << std::endl;


      // Calculating the ranges for accidental subtraction weights Eta Regions// 
      // Center of Eta Signal Region
		double mx2 = mu_Eta;     // Weighted average of double gaussian fit
      /*
      // Alternate center for signal region
		double mx2 = GausMax2;  // Maximum of the fit fsig2
      */

      // sigmawidth = half of signal region rounded to nearest MeV
      // (1000[MeV/GeV conversion] * selected number of sigma width * weighted average std. dev. of double gaussian fit)
      int sigmawidthEta = lround(1000 * sigmaRANGE * sigma_Eta);

      // See if it is close to the next bin
      int AddBinEta = 0;
      if (sigmawidthEta % binwMeV >= 3) {
         AddBinEta = 1;
      }

      // Get the bin value from the real data histogram
      sigmawidthEta /= binwMeV;
      int sigmaBinEta = sigmawidthEta + AddBinEta;

      // Bin value for the double gaussian fit weighted average mean location
      int etamean = remass->FindBin(mx2);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*                                       Invariant Mass Increases -->                                       */
      /* SB1 Region = (EtaLower, EtaLow) || Signal Region = (Etadown, Etaup) || SB2 Region = (EtaHigh, EtaHigher) */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // X-values (GeV/c^2) to use in integrals of a TF1 fit (don't forget to divide by bin width of real histo)
      // This makes the signal region 1 bin larger than needed
		double Etaup   = remass->GetXaxis()->GetBinUpEdge(etamean + sigmaBinEta);
		double Etadown = remass->GetXaxis()->GetBinLowEdge(etamean - sigmaBinEta);

      // Far end of SB Regions
      // Make SB1 Region width = sigmaBin
		double EtaLower  = remass->GetXaxis()->GetBinLowEdge(etamean - 2 * sigmaBinEta - 1);

      // Since, Signal Region = (2 * sigmaBin) + 1 bin (central bin)
      // Make SB2 Region width = sigmaBin + 1 Bin
      // This way width of SB1 + SB2 = Signal Region
      // Remember there is a 1 Bin gap between Sideband and Signal Regions
		double EtaHigher = remass->GetXaxis()->GetBinUpEdge(etamean + 2 * sigmaBinEta + 2);

      // Inner Edge of SB Regions
		double EtaHigh = remass->GetXaxis()->GetBinUpEdge(etamean + sigmaBinEta + 1);  // So that signal and sideband regions do not overlap
		double EtaLow  = remass->GetXaxis()->GetBinLowEdge(etamean - sigmaBinEta - 1); // So that signal and sideband regions do not overlap

      std::cout << std::endl;
      std::cout << "EtaLower = " << EtaLower * 1000 << std::endl;
      std::cout << "EtaLow = " << EtaLow * 1000 << std::endl;
      std::cout << "Etadown = " << Etadown * 1000 << std::endl;
      std::cout << "Etaup = " << Etaup * 1000 << std::endl;
      std::cout << "EtaHigh = " << EtaHigh * 1000 << std::endl;
      std::cout << "EtaHigher = " << EtaHigher * 1000 << std::endl;
      std::cout << std::endl;


      // Calculating the ranges for accidental subtraction weights Etap Regions// 
      // Center of Etap Signal Region
		double mx4 = mean_sig4;  // Maximum of the fit fsig4

      // sigmawidth = half of signal region rounded to nearest MeV
      // (1000[MeV/GeV conversion] * selected number of sigma width * weighted average std. dev. of double gaussian fit)
      int sigmawidthEtap = lround(1000 * sigmaRANGE * sigma_sig4);

      // See if it is close to the next bin
      int AddBinEtap = 0;
      if (sigmawidthEtap % binwMeV >= 3) {
         AddBinEtap = 1;
      }

      // Get the bin value from the real data histogram
      sigmawidthEtap /= binwMeV;
      int sigmaBinEtap = sigmawidthEtap + AddBinEtap;

      // Bin value for the double gaussian fit weighted average mean location
      int etapmean = remass->FindBin(mx4);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /*                                       Invariant Mass Increases -->                                       */
      /* SB1 Region = (EtapLower, EtapLow) || Signal Region = (Etapdown, Etapup) || SB2 Region = (EtapHigh, EtapHigher) */
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // X-values (GeV/c^2) to use in integrals of a TF1 fit (don't forget to divide by bin width of real histo)
      // This makes the signal region 1 bin larger than needed
		double Etapup   = remass->GetXaxis()->GetBinUpEdge(etapmean + sigmaBinEtap);
		double Etapdown = remass->GetXaxis()->GetBinLowEdge(etapmean - sigmaBinEtap);

      // Far end of SB Regions
      // Make SB1 Region width = sigmaBin
		double EtapLower  = remass->GetXaxis()->GetBinLowEdge(etapmean - 2 * sigmaBinEtap - 1);

      // Since, Signal Region = (2 * sigmaBin) + 1 bin (central bin)
      // Make SB2 Region width = sigmaBin + 1 Bin
      // This way width of SB1 + SB2 = Signal Region
      // Remember there is a 1 Bin gap between Sideband and Signal Regions
		double EtapHigher = remass->GetXaxis()->GetBinUpEdge(etapmean + 2 * sigmaBinEtap + 2);

      // Inner Edge of SB Regions
		double EtapHigh = remass->GetXaxis()->GetBinUpEdge(etapmean + sigmaBinEtap + 1);  // So that signal and sideband regions do not overlap
		double EtapLow  = remass->GetXaxis()->GetBinLowEdge(etapmean - sigmaBinEtap - 1); // So that signal and sideband regions do not overlap

      std::cout << std::endl;
      std::cout << "EtapLower = " << EtapLower * 1000 << std::endl;
      std::cout << "EtapLow = " << EtapLow * 1000 << std::endl;
      std::cout << "Etapdown = " << Etapdown * 1000 << std::endl;
      std::cout << "Etapup = " << Etapup * 1000 << std::endl;
      std::cout << "EtapHigh = " << EtapHigh * 1000 << std::endl;
      std::cout << "EtapHigher = " << EtapHigher * 1000 << std::endl;
      std::cout << std::endl;

      gFile = MyFile;
      TString hdirName;
      // Root file with new histograms for reference in later analysis work
      TFile *MyFile = new TFile("TwoG_may_16_2023_STEP2.root", "UPDATE");

		if (Echannel == 0) {
         SBRANGE->SetBinContent(1, EtapLower);
         SBRANGE->SetBinContent(2, EtapLow);
         SBRANGE->SetBinContent(3, Etapdown);
         SBRANGE->SetBinContent(4, Etapup);
         SBRANGE->SetBinContent(5, EtapHigh);
         SBRANGE->SetBinContent(6, EtapHigher);
         
         hdirName.Form("SB_Regions_allt_E%d", Echannel);
         gDirectory->WriteObject(SBRANGE, hdirName);
		}

		else if (Echannel == 1) {
         SBRANGE2->SetBinContent(1, EtapLower);
         SBRANGE2->SetBinContent(2, EtapLow);
         SBRANGE2->SetBinContent(3, Etapdown);
         SBRANGE2->SetBinContent(4, Etapup);
         SBRANGE2->SetBinContent(5, EtapHigh);
         SBRANGE2->SetBinContent(6, EtapHigher);

         hdirName.Form("SB_Regions_allt_E%d", Echannel);
         gDirectory->WriteObject(SBRANGE2, hdirName);
		}

		else if (Echannel == 2) {
         SBRANGE3->SetBinContent(1, EtapLower);
         SBRANGE3->SetBinContent(2, EtapLow);
         SBRANGE3->SetBinContent(3, Etapdown);
         SBRANGE3->SetBinContent(4, Etapup);
         SBRANGE3->SetBinContent(5, EtapHigh);
         SBRANGE3->SetBinContent(6, EtapHigher);

         hdirName.Form("SB_Regions_allt_E%d", Echannel);
         gDirectory->WriteObject(SBRANGE3, hdirName);
		}

		else {
         SBRANGE4->SetBinContent(1, EtapLower);
         SBRANGE4->SetBinContent(2, EtapLow);
         SBRANGE4->SetBinContent(3, Etapdown);
         SBRANGE4->SetBinContent(4, Etapup);
         SBRANGE4->SetBinContent(5, EtapHigh);
         SBRANGE4->SetBinContent(6, EtapHigher);

         hdirName.Form("SB_Regions_allt_E%d", Echannel);
         gDirectory->WriteObject(SBRANGE4, hdirName);
		}

      // Signal Region
      double PercentPiSig = fsig1->Integral(Pidown, Piup) / binw;
      PercentPiSig /= TOTALPi;

      double PercentEtaSig = fsig2->Integral(Etadown, Etaup) / binw;
      PercentEtaSig /= TOTALEta;

      double PercentEtapSig = fsig4->Integral(Etapdown, Etapup) / binw;
      PercentEtapSig /= TOTALEtap;

      // Lower Sideband Region
      double PercentPiSB1 = fsig1->Integral(PiLower, PiLow) / binw;
      PercentPiSB1 /= TOTALPi;

      double PercentEtaSB1 = fsig2->Integral(EtaLower, EtaLow) / binw;
      PercentEtaSB1 /= TOTALEta;

      double PercentEtapSB1 = fsig4->Integral(EtapLower, EtapLow) / binw;
      PercentEtapSB1 /= TOTALEtap;

      // Upper Sideband Region
      double PercentPiSB2 = fsig1->Integral(PiHigh, PiHigher) / binw;
      PercentPiSB2 /= TOTALPi;

      double PercentEtaSB2 = fsig2->Integral(EtaHigh, EtaHigher) / binw;
      PercentEtaSB2 /= TOTALEta;

      double PercentEtapSB2 = fsig4->Integral(EtapHigh, EtapHigher) / binw;
      PercentEtapSB2 /= TOTALEtap;

		// Pi0 Search histos for determination of signal and sideband widths for SB subtraction
		c123->cd();
		gStyle->SetOptStat(0);
		newhfinal2_123->Draw("hist");
		newhfinal2_123->SetMinimum(0.0);
		newhfinal2_123->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_123->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_123->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_123->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_123->SetTitle("#pi^{0} Mass Range");
		newhfinal2_123->SetFillColor(18);
		newhfinal2_123->GetXaxis()->SetRangeUser(0.04, 0.3);
		fback->Draw("same");
		ftot->Draw("same");

		TLine *zline123 = new TLine(Pidown, 0.0, Pidown, 0.5*norm);
		zline123->SetLineColor(1);
		zline123->SetLineStyle(9);
		zline123->SetLineWidth(2);
		zline123->Draw("same");

		TLine *zline123a = new TLine(Piup, 0.0, Piup, 0.5*norm);
		zline123a->SetLineColor(1);
		zline123a->SetLineStyle(9);
		zline123a->SetLineWidth(2);
		zline123a->Draw("same");

		TLine *zline123b = new TLine(PiLow, 0.0, PiLow, 0.25*norm);
		zline123b->SetLineColor(6);
		zline123b->SetLineStyle(9);
		zline123b->SetLineWidth(2);
		zline123b->Draw("same");

		TLine *zline123c = new TLine(PiHigh, 0.0, PiHigh, 0.25*norm);
		zline123c->SetLineColor(6);
		zline123c->SetLineStyle(9);
		zline123c->SetLineWidth(2);
		zline123c->Draw("same");

		TLine *zline123x = new TLine(PiLower, 0.0, PiLower, 0.25*norm);
		zline123x->SetLineColor(6);
		zline123x->SetLineStyle(9);
		zline123x->SetLineWidth(2);
		zline123x->Draw("same");

		TLine *zline123z = new TLine(PiHigher, 0.0, PiHigher, 0.25*norm);
		zline123z->SetLineColor(6);
		zline123z->SetLineStyle(9);
		zline123z->SetLineWidth(2);
		zline123z->Draw("same");

		TLine *zline123d = new TLine(mu_Pi, 0.0, mu_Pi, 1.025*norm);
		zline123d->SetLineColor(4);
		zline123d->SetLineStyle(3);
		zline123d->SetLineWidth(3);
		zline123d->Draw("same");

		TString name123;
		name123.Form("SB1 #rightarrow %.1f to %.1f MeV/c^{2}", PiLower*1000, PiLow*1000);

		TString name123sig;
		name123sig.Form("Signal #rightarrow %.1f to %.1f MeV/c^{2}", Pidown*1000, Piup*1000);

		TString name123a;
		name123a.Form("SB2 #rightarrow %.1f to %.1f MeV/c^{2}", PiHigh*1000, PiHigher*1000);

		TString name123b;
		name123b.Form("SB1 = %.2f%% of #pi^{0} Signal", PercentPiSB1 * 100);

		TString name123d;
		name123d.Form("SB2 = %.2f%% of #pi^{0} Signal", PercentPiSB2 * 100);

		TString name123e;
		name123e.Form("Sig. Region (%.2f#sigma) = %.1f%% of #pi^{0}", sigmaRANGE, PercentPiSig * 100);

		TString name123f;
		name123f.Form("#pi^{0} Fit Mean = %.1f MeV/c^{2}", mu_Pi * 1000);

		TString name123g;
		name123g.Form("#pi^{0} Fit #sigma = %.1f MeV/c^{2}", sigma_Pi * 1000);

		TLegend *legnd123 = new TLegend(0.50,0.55,0.90,0.90);
		legnd123->SetTextSize(0.026);
		legnd123->AddEntry(newhfinal2_123,"Real Data","f");
		legnd123->AddEntry(ftot,"Signal & Background Fit","l");
		legnd123->AddEntry(fback,"Background Fit","l");
		legnd123->AddEntry((TObject*)0, name123f, "");
		legnd123->AddEntry((TObject*)0, name123g, "");
		legnd123->AddEntry((TObject*)0, name123, "");
		legnd123->AddEntry((TObject*)0, name123sig, "");
		legnd123->AddEntry((TObject*)0, name123a, "");
		legnd123->AddEntry((TObject*)0, name123b, "");
		legnd123->AddEntry((TObject*)0, name123e, "");
		legnd123->AddEntry((TObject*)0, name123d, "");
		legnd123->Draw("same");
      TString xc123;
      xc123.Form("PiFit_%d.png", Echannel);
		c123->Print(xc123);

		// Eta Search
		double normEta = fsig2->Eval(GausMax2);
		c321->cd();
		gStyle->SetOptStat(0);
		newhfinal2_321->Draw("hist");
		newhfinal2_321->SetMinimum(0.0);
		newhfinal2_321->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_321->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_321->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_321->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_321->SetTitle("#eta Mass Range");
		newhfinal2_321->SetFillColor(18);
		newhfinal2_321->GetXaxis()->SetRangeUser(0.28, 0.85);
		fback->Draw("same");
      fsig2->Draw("same");
      fsig3->Draw("same");
		ftot->Draw("same");

		TLine *zline321 = new TLine(Etadown, 0.0, Etadown, 0.6*normEta);
		zline321->SetLineColor(1);
		zline321->SetLineStyle(9);
		zline321->SetLineWidth(2);
		zline321->Draw("same");

		TLine *zline321a = new TLine(Etaup, 0.0, Etaup, 0.6*normEta);
		zline321a->SetLineColor(1);
		zline321a->SetLineStyle(9);
		zline321a->SetLineWidth(2);
		zline321a->Draw("same");

		TLine *zline321b = new TLine(EtaLower, 0.0, EtaLower, 0.28*normEta);
		zline321b->SetLineColor(6);
		zline321b->SetLineStyle(9);
		zline321b->SetLineWidth(2);
		zline321b->Draw("same");

		TLine *zline321c = new TLine(EtaHigher, 0.0, EtaHigher, 0.28*normEta);
		zline321c->SetLineColor(6);
		zline321c->SetLineStyle(9);
		zline321c->SetLineWidth(2);
		zline321c->Draw("same");

		TLine *zline321x = new TLine(EtaLow, 0.0, EtaLow, 0.28*normEta);
		zline321x->SetLineColor(6);
		zline321x->SetLineStyle(9);
		zline321x->SetLineWidth(2);
		zline321x->Draw("same");

		TLine *zline321z = new TLine(EtaHigh, 0.0, EtaHigh, 0.28*normEta);
		zline321z->SetLineColor(6);
		zline321z->SetLineStyle(9);
		zline321z->SetLineWidth(2);
		zline321z->Draw("same");

		TLine *zline321d = new TLine(mu_Eta, 0.0, mu_Eta, 1.05*normEta);
		zline321d->SetLineColor(4);
		zline321d->SetLineStyle(3);
		zline321d->SetLineWidth(3);
		zline321d->Draw("same");

		TString name321;
		name321.Form("SB1 #rightarrow %.1f to %.1f MeV/c^{2}", EtaLower*1000, EtaLow*1000);

		TString name321sig;
		name321sig.Form("Signal #rightarrow %.1f to %.1f MeV/c^{2}", Etadown*1000, Etaup*1000);

		TString name321a;
		name321a.Form("SB2 #rightarrow %.1f to %.1f MeV/c^{2}", EtaHigh*1000, EtaHigher*1000);

		TString name321b;
		name321b.Form("SB1 = %.2f%% of #eta Signal", PercentEtaSB1 * 100);

		TString name321c;
		name321c.Form("#eta Fit Mean = %.1f MeV/c^{2}", mu_Eta * 1000);

		TString name321d;
		name321d.Form("SB2 = %.2f%% of #eta Signal", PercentEtaSB2 * 100);

		TString name321e;
		name321e.Form("Sig. Region (%.2f#sigma) = %.1f%% of #eta", sigmaRANGE, PercentEtaSig * 100);

		TString name321g;
		name321g.Form("#eta Fit #sigma = %.1f MeV/c^{2}", sigma_Eta * 1000);

		TLegend *legnd321 = new TLegend(0.55,0.56,0.90,0.90);
		legnd321->SetTextSize(0.026);
		legnd321->AddEntry(newhfinal2_321,"Real Data","f");
		legnd321->AddEntry(fsig2,"#eta Signal Fit","l");
		legnd321->AddEntry(ftot,"Signal & Background Fit","l");
		legnd321->AddEntry(fback,"Background Fit","l");
		legnd321->AddEntry((TObject*)0, name321c, "");
		legnd321->AddEntry((TObject*)0, name321g, "");
		legnd321->AddEntry((TObject*)0, name321, "");
		legnd321->AddEntry((TObject*)0, name321sig, "");
		legnd321->AddEntry((TObject*)0, name321a, "");
		legnd321->AddEntry((TObject*)0, name321b, "");
		legnd321->AddEntry((TObject*)0, name321e, "");
		legnd321->AddEntry((TObject*)0, name321d, "");
		legnd321->Draw("same");
      TString xc321;
      xc321.Form("EtaFit_%d.png", Echannel);
		c321->Print(xc321);

		// Etap Search
		double normEtap = fsig4->Eval(mean_sig4);
		double EtapPlot = ftot->Eval(0.800);
		c555->cd();
		gStyle->SetOptStat(0);
		newhfinal2_555->Draw("hist");
		newhfinal2_555->SetMaximum(EtapPlot);
		newhfinal2_555->SetMinimum(0.0);
		newhfinal2_555->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_555->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_555->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_555->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_555->SetTitle("#eta' Mass Range");
		newhfinal2_555->SetFillColor(18);
		newhfinal2_555->GetXaxis()->SetRangeUser(0.74, 1.18);
		fback->Draw("same");
      fsig3->Draw("same");
      fsig4->Draw("same");
		ftot->Draw("same");

		TLine *zline555 = new TLine(Etapdown, 0.0, Etapdown, 1.0*normEtap);
		zline555->SetLineColor(1);
		zline555->SetLineStyle(9);
		zline555->SetLineWidth(2);
		zline555->Draw("same");

		TLine *zline555a = new TLine(Etapup, 0.0, Etapup, 1.0*normEtap);
		zline555a->SetLineColor(1);
		zline555a->SetLineStyle(9);
		zline555a->SetLineWidth(2);
		zline555a->Draw("same");

		TLine *zline555b = new TLine(EtapLower, 0.0, EtapLower, 0.75*normEtap);
		zline555b->SetLineColor(6);
		zline555b->SetLineStyle(9);
		zline555b->SetLineWidth(2);
		zline555b->Draw("same");

		TLine *zline555c = new TLine(EtapHigher, 0.0, EtapHigher, 0.75*normEtap);
		zline555c->SetLineColor(6);
		zline555c->SetLineStyle(9);
		zline555c->SetLineWidth(2);
		zline555c->Draw("same");

		TLine *zline555x = new TLine(EtapLow, 0.0, EtapLow, 0.75*normEtap);
		zline555x->SetLineColor(6);
		zline555x->SetLineStyle(9);
		zline555x->SetLineWidth(2);
		zline555x->Draw("same");

		TLine *zline555z = new TLine(EtapHigh, 0.0, EtapHigh, 0.75*normEtap);
		zline555z->SetLineColor(6);
		zline555z->SetLineStyle(9);
		zline555z->SetLineWidth(2);
		zline555z->Draw("same");

		TLine *zline555d = new TLine(mean_sig4, 0.0, mean_sig4, 1.25*normEtap);
		zline555d->SetLineColor(4);
		zline555d->SetLineStyle(3);
		zline555d->SetLineWidth(3);
		zline555d->Draw("same");

		TString name555;
		name555.Form("SB1 #rightarrow %.1f to %.1f MeV/c^{2}", EtapLower*1000, EtapLow*1000);

		TString name555sig;
		name555sig.Form("Signal #rightarrow %.1f to %.1f MeV/c^{2}", Etapdown*1000, Etapup*1000);

		TString name555a;
		name555a.Form("SB2 #rightarrow %.1f to %.1f MeV/c^{2}", EtapHigh*1000, EtapHigher*1000);

		TString name555b;
		name555b.Form("SB1 = %.2f%% of #eta' Signal", PercentEtapSB1 * 100);

		TString name555c;
		name555c.Form("#eta' Fit Mean = %.1f MeV/c^{2}", mean_sig4 * 1000);

		TString name555d;
		name555d.Form("SB2 = %.2f%% of #eta' Signal", PercentEtapSB2 * 100);

		TString name555e;
		name555e.Form("Sig. Region (%.2f#sigma) = %.1f%% of #eta'", sigmaRANGE, PercentEtapSig * 100);

		TString name555g;
		name555g.Form("#eta' Fit #sigma = %.1f MeV/c^{2}", sigma_sig4 * 1000);

		TLegend *legnd555 = new TLegend(0.50,0.54,0.90,0.90);
		legnd555->SetTextSize(0.026);
		legnd555->AddEntry(newhfinal2_555,"Real Data","f");
		legnd555->AddEntry(fsig4,"#eta' Signal Fit","l");
		legnd555->AddEntry(ftot,"Signal & Background Fit","l");
		legnd555->AddEntry(fback,"Background Fit","l");
		legnd555->AddEntry((TObject*)0, name555c, "");
		legnd555->AddEntry((TObject*)0, name555g, "");
		legnd555->AddEntry((TObject*)0, name555, "");
		legnd555->AddEntry((TObject*)0, name555sig, "");
		legnd555->AddEntry((TObject*)0, name555a, "");
		legnd555->AddEntry((TObject*)0, name555b, "");
		legnd555->AddEntry((TObject*)0, name555e, "");
		legnd555->AddEntry((TObject*)0, name555d, "");
		legnd555->Draw("same");
      TString xc555;
      xc555.Form("EtapFit_%d.png", Echannel);
		c555->Print(xc555);

      // Extra calculation, not really needed
		// Signal and background in SIGNAL REGIONS /////////////
		double counts_Pi   = ftot->Integral(Pidown, Piup) / binw;
		double counts_Eta  = ftot->Integral(Etadown, Etaup) / binw;
		double counts_Etap = ftot->Integral(Etapdown, Etapup) / binw;


		/* Signal only in SIGNAL REGION */
      /* S_0 */
		double Sigcount_Pi   = fsig1->Integral(Pidown, Piup) / binw;
		double Sigcount_Eta  = fsig2->Integral(Etadown, Etaup) / binw;
		double Sigcount_Etap = fsig4->Integral(Etapdown, Etapup) / binw;


		/* Sideband 1 Region (left side). Background only. */
      /* B_1 */
      // Neutral Pion SB1 Region
      // No other signals show up this low in invariant mass
		double sb1count_Pi   = fback->Integral(PiLower, PiLow) / binw;

      // Eta SB1 Region
		double ab1sb1        = fback->Integral(EtaLower, EtaLow) / binw;
		double ab2sb1        = fsig1->Integral(EtaLower, EtaLow) / binw;
		double ab3sb1        = fsig3->Integral(EtaLower, EtaLow) / binw;
		double ab4sb1        = fsig4->Integral(EtaLower, EtaLow) / binw;
		double sb1count_Eta  = ab1sb1 + ab2sb1 + ab3sb1 + ab4sb1;

      // Eta prime SB1 Region
      // No Pi0 signal (fsig1) in this region
		double ac1sb1        = fback->Integral(EtapLower, EtapLow) / binw;
		double ac2sb1        = fsig2->Integral(EtapLower, EtapLow) / binw;
		double ac3sb1        = fsig3->Integral(EtapLower, EtapLow) / binw;
		double sb1count_Etap = ac1sb1 + ac2sb1 + ac3sb1;


		/* Sideband 2 Region (right side). Background only. */
      /* B_2 */
      // Neutral Pion SB2 Region
		double aa1sb2        = fback->Integral(PiHigh, PiHigher) / binw;
		double aa2sb2        = fsig2->Integral(PiHigh, PiHigher) / binw;
		double aa3sb2        = fsig3->Integral(PiHigh, PiHigher) / binw;
		double sb2count_Pi   = aa1sb2 + aa2sb2 + aa3sb2;

      // Eta SB2 Region
		double ab1sb2        = fback->Integral(EtaHigh, EtaHigher) / binw;
		double ab2sb2        = fsig1->Integral(EtaHigh, EtaHigher) / binw;
		double ab3sb2        = fsig3->Integral(EtaHigh, EtaHigher) / binw;
		double ab4sb2        = fsig4->Integral(EtaHigh, EtaHigher) / binw;
		double sb2count_Eta  = ab1sb2 + ab2sb2 + ab3sb2 + ab4sb2;

      // Eta prime SB2 Region
		double ac1sb2        = fback->Integral(EtapHigh, EtapHigher) / binw;
		double ac2sb2        = fsig2->Integral(EtapHigh, EtapHigher) / binw;
		double ac3sb2        = fsig3->Integral(EtapHigh, EtapHigher) / binw;
		double sb2count_Etap = ac1sb2 + ac2sb2 + ac3sb2;


		/* Background only in SIGNAL REGION */
      /* B_0 */
      // Neutral Pion Signal Region
      // No other signals show up this low in invariant mass
		double Bg_Pi        = fback->Integral(Pidown, Piup) / binw;
		double BgEta_Pi     = fback->Integral(Pidown, Piup) / binw;
		double BgLeak_Pi    = fback->Integral(Pidown, Piup) / binw;
		double Bgcount_Pi   = Bg_Pi + BgEta_Pi + BgLeak_Pi;

      // Eta Signal Region
		double BgPi_Eta     = fsig1->Integral(Etadown, Etaup) / binw;
		double Bgleak_Eta   = fsig3->Integral(Etadown, Etaup) / binw;
		double BgEtap_Eta   = fsig4->Integral(Etadown, Etaup) / binw;
		double BgBack_Eta   = fback->Integral(Etadown, Etaup) / binw;
		double Bgcount_Eta  = BgPi_Eta + Bgleak_Eta + BgEtap_Eta + BgBack_Eta;

      // Eta prime Signal Region
		double BgEta_Etap   = fsig2->Integral(Etapdown, Etapup) / binw;
		double Bgleak_Etap  = fsig3->Integral(Etapdown, Etapup) / binw;
		double BgBack_Etap  = fback->Integral(Etapdown, Etapup) / binw;
		double Bgcount_Etap = BgEta_Etap + Bgleak_Etap + BgBack_Etap;


		/* Sideband 1 Region (left side). Signal only. */
      /* S_1 */
      // Neutral Pion Signal Region
		double SigSb1_Pi   = fsig1->Integral(PiLower, PiLow) / binw;

      // Eta Signal Region
		double SigSb1_Eta  = fsig2->Integral(EtaLower, EtaLow) / binw;

      // Eta prime Signal Region
		double SigSb1_Etap = fsig4->Integral(EtapLower, EtapLow) / binw;


		/* Sideband 2 Region (right side). Signal only. */
      /* S_2 */
      // Neutral Pion Signal Region
		double SigSb2_Pi   = fsig1->Integral(PiHigh, PiHigher) / binw;

      // Eta Signal Region
		double SigSb2_Eta  = fsig2->Integral(EtaHigh, EtaHigher) / binw;

      // Eta prime Signal Region
		double SigSb2_Etap = fsig4->Integral(EtapHigh, EtapHigher) / binw;


		/* RATIOS USED TO CALCULATE ACCIDENTAL SUBTRACTION WEIGHTS */
      // Neutral Pion Signal Region
		double R1_Pi   = SigSb1_Pi / Sigcount_Pi;       // R1 = s1/s0
		double R2_Pi   = SigSb2_Pi / Sigcount_Pi;       // R2 = s2/s0
		double C1_Pi   = sb1count_Pi / Bgcount_Pi;      // C1 = b1/b0
		double C2_Pi   = sb2count_Pi / Bgcount_Pi;      // C2 = b2/b0

      // Eta Signal Region
		double R1_Eta  = SigSb1_Eta / Sigcount_Eta;     // R1 = s1/s0
		double R2_Eta  = SigSb2_Eta / Sigcount_Eta;     // R2 = s2/s0
		double C1_Eta  = sb1count_Eta / Bgcount_Eta;    // C1 = b1/b0
		double C2_Eta  = sb2count_Eta / Bgcount_Eta;    // C2 = b2/b0

      // Eta prime Signal Region
		double R1_Etap = SigSb1_Etap / Sigcount_Etap;   // R1 = s1/s0
		double R2_Etap = SigSb2_Etap / Sigcount_Etap;   // R2 = s2/s0
		double C1_Etap = sb1count_Etap / Bgcount_Etap;  // C1 = b1/b0
		double C2_Etap = sb2count_Etap / Bgcount_Etap;  // C2 = b2/b0

      ////////////////////////////////////////////////////
      /*          EQUAL SIDEBAND WEIGHT METHOD          */
		////// Sideband subtraction weighting factors //////
		// w0 = signal                                    //
		// w0 = (1 + R1 + R2) / (1 - (R1 + R2)/(C1 + C2)) //
		////////////////////////////////////////////////////
		// w1 = sb left || w2 = sb right                  //
		// w1 = w2 = -(w0 / (C1 + C2))                    //
		////////////////////////////////////////////////////

      ////////////////////////////////////////////////////
      /*                ALTERNATE METHOD                */
      /*         UNEQUAL SIDEBAND WEIGHT METHOD         */
		////// Sideband subtraction weighting factors //////
		// w0 = signal                                    //
		// w0 = (1 + R1 + R2) /                           //
      //      (1 - ((C1*R1 + C2*R2) / (C1*C1 + C2*C2))) //
		////////////////////////////////////////////////////
		// w1 = sb left                                   //
		// w1 = -w0 * (C1 / (C1*C1 + C2*C2))              //
      //                                                //
		// w2 = sb right                                  //
		// w1 = -w0 * (C1 / (C1*C1 + C2*C2))              //
		////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////
      /* Use equal sideband weight method for Pi0              */
      /* Use unequal sideband weight method for Eta & Etaprime */
      ///////////////////////////////////////////////////////////

		/////////////////// For Pi ///////////////////
		// w0 for Pi0
		double w_pi_sig = (1.0 + R1_Pi + R2_Pi) / (1.0 - ((R1_Pi + R2_Pi) / (C1_Pi + C2_Pi)));
		// w1 = w2 for Pi0
		double w_pi_sb  = -(w_pi_sig / (C1_Pi + C2_Pi));

		/////////////////// For Eta ///////////////////
		// w0 for Eta
		double w_eta_sig = (1 + R1_Eta + R2_Eta) / (1 - ((C1_Eta * R1_Eta + C2_Eta * R2_Eta) / (C1_Eta * C1_Eta + C2_Eta * C2_Eta)));
		// w1 for Eta
		double w_eta_sb1  = -(w_eta_sig) * (C1_Eta / (C1_Eta * C1_Eta + C2_Eta * C2_Eta));
		// w2 for Eta
		double w_eta_sb2  = -(w_eta_sig) * (C2_Eta / (C1_Eta * C1_Eta + C2_Eta * C2_Eta));

		/////////////////// For Etap ///////////////////
		// w0 for Etap
		double w_etap_sig = (1 + R1_Etap + R2_Etap) / (1 - ((C1_Etap * R1_Etap + C2_Etap * R2_Etap) / (C1_Etap * C1_Etap + C2_Etap * C2_Etap)));
		// w1 for Etap
		double w_etap_sb1  = -(w_etap_sig) * (C1_Etap / (C1_Etap * C1_Etap + C2_Etap * C2_Etap));
		// w2 for Etap
		double w_etap_sb2  = -(w_etap_sig) * (C2_Etap / (C1_Etap * C1_Etap + C2_Etap * C2_Etap));
	 

		// Ratio of Signal-to-Background for the +/-2sigma signal region 
		double ratio_pi   = Sigcount_Pi   / Bgcount_Pi;
		double ratio_eta  = Sigcount_Eta  / Bgcount_Eta;
		double ratio_etap = Sigcount_Etap / Bgcount_Etap;

		/////////////////////////////////////////////////
		/////////////// Fit Summary /////////////////////
		/////////////////////////////////////////////////
		// Pi0 Search
		c1->cd(1);
		gStyle->SetOptStat(0);
		newhfinal2_1->Draw("hist");
		newhfinal2_1->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_1->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_1->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_1->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_1->SetTitle("#pi^{0} Mass Range");
		newhfinal2_1->SetFillColor(18);
		newhfinal2_1->GetXaxis()->SetRangeUser(0.03, 0.26);
		fback->Draw("same");
		ftot->Draw("same");

		TLegend *legnd5 = new TLegend(0.6, 0.68, 0.89, 0.89);
		legnd5->SetTextSize(0.028);
		legnd5->AddEntry(newhfinal2_1,"Real Data","f");
		legnd5->AddEntry(ftot,"Signal & Background Fit","l");
		legnd5->AddEntry(fback,"Background Fit","l");
		legnd5->Draw("same");

		c1->cd(2);
		newhfinal2_1b->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_1b->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_1b->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_1b->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_1b->SetFillColor(18);
		newhfinal2_1b->SetMinimum(0.0);
		newhfinal2_1b->GetXaxis()->SetRangeUser(0.03, 0.26);

      TString NEWtitle;
		if (Echannel == 0) {
	      NEWtitle.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range0b, Range0a);
		}
		else if (Echannel == 1) {
	      NEWtitle.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range1b, Range1a);
		}
		else if (Echannel == 2) {
	      NEWtitle.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range2b, Range2a);
		}
		else {
         NEWtitle.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range3b, Range3a);
		}

		newhfinal2_1b->SetTitle(NEWtitle);
		newhfinal2_1b->Draw("hist");
		fback->Draw("same");
		fsig1->Draw("same");

		TLine *zline1232 = new TLine(Pidown, 0.0, Pidown, 0.5*norm);
		zline1232->SetLineColor(1);
		zline1232->SetLineStyle(9);
		zline1232->SetLineWidth(2);
		zline1232->Draw("same");

		TLine *zline1232a = new TLine(Piup, 0.0, Piup, 0.5*norm);
		zline1232a->SetLineColor(1);
		zline1232a->SetLineStyle(9);
		zline1232a->SetLineWidth(2);
		zline1232a->Draw("same");

		TLine *zline1232b = new TLine(PiLow, 0.0, PiLow, 0.25*norm);
		zline1232b->SetLineColor(6);
		zline1232b->SetLineStyle(9);
		zline1232b->SetLineWidth(2);
		zline1232b->Draw("same");

		TLine *zline1232c = new TLine(PiHigh, 0.0, PiHigh, 0.25*norm);
		zline1232c->SetLineColor(6);
		zline1232c->SetLineStyle(9);
		zline1232c->SetLineWidth(2);
		zline1232c->Draw("same");

		TLine *zline1232x = new TLine(PiLower, 0.0, PiLower, 0.25*norm);
		zline1232x->SetLineColor(6);
		zline1232x->SetLineStyle(9);
		zline1232x->SetLineWidth(2);
		zline1232x->Draw("same");

		TLine *zline1232z = new TLine(PiHigher, 0.0, PiHigher, 0.25*norm);
		zline1232z->SetLineColor(6);
		zline1232z->SetLineStyle(9);
		zline1232z->SetLineWidth(2);
		zline1232z->Draw("same");

		TLine *zline1232d = new TLine(mu_Pi, 0.0, mu_Pi, 1.025*norm);
		zline1232d->SetLineColor(4);
		zline1232d->SetLineStyle(3);
		zline1232d->SetLineWidth(3);
		zline1232d->Draw("same");

		TLegend *legnd5b = new TLegend(0.6,0.58,0.89,0.89);
		legnd5b->SetTextSize(0.028);
		legnd5b->AddEntry(newhfinal2_1b,"Real Data","f");
		legnd5b->AddEntry(fsig1,"#pi^{0} Signal Fit","l");
		legnd5b->AddEntry(fback,"Background Fit","l");
      TString Leg5b;
      Leg5b.Form("%.2f#sigma Signal Region", sigmaRANGE);
		legnd5b->AddEntry(zline1232, Leg5b, "l");
		legnd5b->Draw("same");

		c1->cd(3);
		newhfinal2_1a->GetXaxis()->SetRangeUser(0.03, 0.26);
		newhfinal2_1a->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
		newhfinal2_1a->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_1a->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_1a->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_1a->SetTitle("Background in #pi^{0} Region");
		fback->SetFillColor(8);
		newhfinal2_1a->Add(newhfinal2_5, -2);
		double mini1a = fback->Eval(0.24);
		mini1a += 200;
		newhfinal2_1a->SetMaximum(mini1a);
		mini1a *= -0.08;
		newhfinal2_1a->SetMinimum(mini1a);
		newhfinal2_1a->Draw("hist");
		fback->Draw("same");

		TLine *zeroline = new TLine(0.03, 0, 0.26, 0);
		zeroline->SetLineColor(1);
		zeroline->SetLineStyle(10);
		zeroline->SetLineWidth(2);
		zeroline->Draw("same");

		c1->cd(4);
		gStyle->SetOptFit(1100);
		newhfinal2_1c->GetXaxis()->SetRangeUser(0.03, 0.26);
		newhfinal2_1c->SetTitle("#pi^{0} Signal Fit");
		newhfinal2_1c->SetLineColor(0);
		newhfinal2_1c->Draw("hist");
		fsig1->Draw("same");

		TLine *zerolineD2DD = new TLine(GausMax1, 0, GausMax1, 1.03*norm);
		zerolineD2DD->SetLineColor(6);
		zerolineD2DD->SetLineStyle(9);
		zerolineD2DD->SetLineWidth(2);
		zerolineD2DD->Draw("same");

		double pidiff = (abs(((GausMax1*1000)-134.9766)/134.9766))*100;
		TString titlepi;
		TString namepi;
		namepi.Form("#pi^{0} Fit Peak = %.1f MeV/c^{2}", GausMax1*1000);
		titlepi.Form("Difference = %.2f %%", pidiff);

		TLegend *legnd5a = new TLegend(0.6, 0.68, 0.89, 0.89);
		legnd5a->SetTextSize(0.028);
		legnd5a->AddEntry(fsig1,"#pi^{0} Signal Fit","l");
		legnd5a->AddEntry(zerolineD2DD, namepi, "l");
		legnd5a->AddEntry((TObject*)0, "PDG Value = 135.0 MeV/c^{2}", "");
		legnd5a->AddEntry((TObject*)0, titlepi, "");
		legnd5a->Draw("same");
		char histo1[20];
		sprintf(histo1,"%s%d%s","Pi_mass", totalhisto, ".png");
		c1->Print(histo1);

		///////////////////////////////////////////////////////////////////
		// Eta Search
		c2->cd(1);
		double normE = fsig2->Eval(GausMax2);
		gStyle->SetOptStat(0);
		newhfinal2_2->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_2->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_2->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_2->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_2->SetTitle("#eta Mass Range");
		newhfinal2_2->SetFillColor(18);
		newhfinal2_2->SetMinimum(0.0);
		newhfinal2_2->GetXaxis()->SetRangeUser(0.3, 0.8);
		newhfinal2_2->Draw("hist");
      fsig3->Draw("same");
		fback->Draw("same");
		ftot->Draw("same");

		TLegend *legnd4 = new TLegend(0.60, 0.745, 0.89, 0.89);
		legnd4->SetTextSize(0.028);
		legnd4->AddEntry(newhfinal2_2,"Real Data","f");
		legnd4->AddEntry(ftot,"Signal & Background Fit","l");
		legnd4->AddEntry(fback,"Background Fit","l");
		legnd4->Draw("same");

		c2->cd(2);
		newhfinal2_2b->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
		newhfinal2_2b->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_2b->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_2b->GetYaxis()->SetTitleOffset(1.5);

      TString title;
		if (Echannel == 0) {
	      title.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range0b, Range0a);
		}
		else if (Echannel == 1) {
	      title.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range1b, Range1a);
		}
		else if (Echannel == 2) {
	      title.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range2b, Range2a);
		}
		else {
         title.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range3b, Range3a);
		}

      newhfinal2_2b->SetTitle(title);
		newhfinal2_2b->SetFillColor(18);
		newhfinal2_2b->SetMinimum(0.0);
		newhfinal2_2b->GetXaxis()->SetRangeUser(0.3, 0.8);
		newhfinal2_2b->Draw("hist");
		fsig2->Draw("same");
      fsig3->Draw("same");
		fback->Draw("same");

		TLine *zline3212 = new TLine(Etadown, 0.0, Etadown, 0.6*normEta);
		zline3212->SetLineColor(1);
		zline3212->SetLineStyle(9);
		zline3212->SetLineWidth(2);
		zline3212->Draw("same");

		TLine *zline3212a = new TLine(Etaup, 0.0, Etaup, 0.6*normEta);
		zline3212a->SetLineColor(1);
		zline3212a->SetLineStyle(9);
		zline3212a->SetLineWidth(2);
		zline3212a->Draw("same");

		TLine *zline3212b = new TLine(EtaLower, 0.0, EtaLower, 0.28*normEta);
		zline3212b->SetLineColor(6);
		zline3212b->SetLineStyle(9);
		zline3212b->SetLineWidth(2);
		zline3212b->Draw("same");

		TLine *zline3212c = new TLine(EtaHigher, 0.0, EtaHigher, 0.28*normEta);
		zline3212c->SetLineColor(6);
		zline3212c->SetLineStyle(9);
		zline3212c->SetLineWidth(2);
		zline3212c->Draw("same");

		TLine *zline3212x = new TLine(EtaLow, 0.0, EtaLow, 0.28*normEta);
		zline3212x->SetLineColor(6);
		zline3212x->SetLineStyle(9);
		zline3212x->SetLineWidth(2);
		zline3212x->Draw("same");

		TLine *zline3212z = new TLine(EtaHigh, 0.0, EtaHigh, 0.28*normEta);
		zline3212z->SetLineColor(6);
		zline3212z->SetLineStyle(9);
		zline3212z->SetLineWidth(2);
		zline3212z->Draw("same");

		TLine *zline3212d = new TLine(mu_Eta, 0.0, mu_Eta, 1.05*normEta);
		zline3212d->SetLineColor(4);
		zline3212d->SetLineStyle(3);
		zline3212d->SetLineWidth(3);
		zline3212d->Draw("same");

		TLegend *legnd5c = new TLegend(0.60,0.58,0.89,0.89);
		legnd5c->SetTextSize(0.028);
		legnd5c->AddEntry(newhfinal2_2b,"Real Data","f");
		legnd5c->AddEntry(fsig1,"#eta Signal Fit","l");
		legnd5c->AddEntry(fback,"Background Fit","l");
      TString Leg5c;
      Leg5c.Form("%.2f#sigma Signal Region", sigmaRANGE);
		legnd5c->AddEntry(zline3212, Leg5c, "l");
		legnd5c->Draw("same");

		c2->cd(3);
		newhfinal2_2a->GetXaxis()->SetRangeUser(0.3, 0.8);
		newhfinal2_2a->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
		newhfinal2_2a->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_2a->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_2a->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_2a->SetTitle("Background in #eta Region");
		fback->SetFillColor(8);
		newhfinal2_2a->Add(newhfinal2_5a, -2);
		double mini2 = fsig3->Eval(0.75);
		mini2 += 200;
		newhfinal2_2a->SetMaximum(mini2*1.1);
		mini2 *= -0.08;
		newhfinal2_2a->SetMinimum(mini2);
		newhfinal2_2a->Draw("hist");
      fsig3->Draw("same");
		fback->Draw("same");

		TLine *zeroline4 = new TLine(0.3, 0, 0.8, 0);
		zeroline4->SetLineColor(1);
		zeroline4->SetLineStyle(10);
		zeroline4->SetLineWidth(2);
		zeroline4->Draw("same");

		c2->cd(4);
		gStyle->SetOptFit(1100);
		newhfinal2_2c->Add(newhfinal2_7, 300);
		newhfinal2_2c->SetTitle("#eta Signal Fit");
		newhfinal2_2c->Draw("hist");
		newhfinal2_2c->GetXaxis()->SetRangeUser(0.3, 0.8);
		newhfinal2_2c->SetTitle("#eta Signal Fit");
		newhfinal2_2c->SetMinimum(0.0);
		newhfinal2_2c->SetMaximum(normE*1.05);
		fsig2->Draw("same");

		TLine *zerolineD2D = new TLine(GausMax2, 0, GausMax2, 1.03*normE);
		zerolineD2D->SetLineColor(6);
		zerolineD2D->SetLineStyle(9);
		zerolineD2D->SetLineWidth(2);
		zerolineD2D->Draw("same");

		double etadiff = (abs(((GausMax2*1000) - 547.862) / 547.862))*100;
		TString titleeta;
		TString nameeta;
		nameeta.Form("#eta Fit Peak = %.1f MeV/c^{2}", GausMax2*1000);
		titleeta.Form("Difference = %.2f %%", etadiff);

		TLegend *legndACC = new TLegend(0.6, 0.68, 0.89, 0.89);
		legndACC->SetTextSize(0.028);
		legndACC->AddEntry(fsig2,"#eta Signal Fit","l");
		legndACC->AddEntry(zerolineD2D, nameeta, "l");
		legndACC->AddEntry((TObject*)0, "PDG Value = 547.9 MeV/c^{2}", "");
		legndACC->AddEntry((TObject*)0, titleeta, "");
		legndACC->Draw("same");

		sprintf(histo1,"%s%d%s","Eta_mass", totalhisto, ".png");
		c2->Print(histo1);

		///////////////////////////////////////////////////////////////////
		// Etap Search
		c3->cd(1);
		gStyle->SetOptStat(0);
		newhfinal2_3->Draw("hist");
		newhfinal2_3->SetMinimum(0.0);
		newhfinal2_3->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
		newhfinal2_3->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_3->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_3->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_3->SetTitle("#eta' Mass Range");
		newhfinal2_3->SetFillColor(18);
		newhfinal2_3->SetLineColor(1);
		newhfinal2_3->GetXaxis()->SetRangeUser(0.74, 1.18);
		fback->Draw("same");
      fsig3->Draw("same");
		ftot->Draw("same");

		TLegend *legnd = new TLegend(0.6, 0.68, 0.89, 0.89);
		legnd->SetTextSize(0.028);
		legnd->AddEntry(newhfinal2_3, "Real Data", "f");
		legnd->AddEntry(ftot, "Signal & Background Fit", "l");
		legnd->AddEntry(fback, "Background Fit", "l");
		legnd->Draw("same");

		c3->cd(2);
		double normE3 = ftot->Eval(mean_sig4);
		double mini3 = -0.08 * normE3;
		newhfinal2_3b->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
		newhfinal2_3b->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_3b->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_3b->GetYaxis()->SetTitleOffset(1.5);

      TString TitleX;
		if (Echannel == 0) {
	      TitleX.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range0b, Range0a);
		}
		else if (Echannel == 1) {
	      TitleX.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range1b, Range1a);
		}
		else if (Echannel == 2) {
	      TitleX.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range2b, Range2a);
		}
		else {
         TitleX.Form("%.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4}", Range3b, Range3a);
		}

      newhfinal2_3b->SetTitle(TitleX);
		newhfinal2_3b->GetXaxis()->SetRangeUser(0.74, 1.18);
		newhfinal2_3b->SetFillColor(18);
		newhfinal2_3b->SetMinimum(0.0);
		newhfinal2_3b->Draw("hist");
		fback->Draw("same");
      fsig3->Draw("same");
		fsig4->Draw("same");

		TLine *zerolineA21 = new TLine(Etapdown, 0.0, Etapdown, 1.1*normEtap);
		zerolineA21->SetLineColor(kBlue);
		zerolineA21->SetLineStyle(9);
		zerolineA21->SetLineWidth(2);
		zerolineA21->Draw("same");

		TLine *zerolineB2 = new TLine(Etapup, 0.0, Etapup, 1.1*normEtap);
		zerolineB2->SetLineColor(kBlue);
		zerolineB2->SetLineStyle(9);
		zerolineB2->SetLineWidth(2);
		zerolineB2->Draw("same");

		TLine *zerolineC21 = new TLine(EtapLower, 0.0, EtapLower, 0.9*normEtap);
		zerolineC21->SetLineColor(2);
		zerolineC21->SetLineStyle(9);
		zerolineC21->SetLineWidth(2);
		zerolineC21->Draw("same");

		TLine *zerolineD21 = new TLine(EtapHigher, 0.0, EtapHigher, 0.9*normEtap);
		zerolineD21->SetLineColor(2);
		zerolineD21->SetLineStyle(9);
		zerolineD21->SetLineWidth(2);
		zerolineD21->Draw("same");

		TLine *zerolineE21 = new TLine(EtapLow, 0.0, EtapLow, 0.9*normEtap);
		zerolineE21->SetLineColor(2);
		zerolineE21->SetLineStyle(9);
		zerolineE21->SetLineWidth(2);
		zerolineE21->Draw("same");

		TLine *zerolineF21 = new TLine(EtapHigh, 0.0, EtapHigh, 0.9*normEtap);
		zerolineF21->SetLineColor(2);
		zerolineF21->SetLineStyle(9);
		zerolineF21->SetLineWidth(2);
		zerolineF21->Draw("same");

		TLine *zerolineG21 = new TLine(mean_sig4, 0.0, mean_sig4, 1.35*normEtap);
		zerolineG21->SetLineColor(4);
		zerolineG21->SetLineStyle(3);
		zerolineG21->SetLineWidth(3);
		zerolineG21->Draw("same");

		TLegend *legndz = new TLegend(0.36, 0.74, 0.64, 0.89);
		legndz->SetTextSize(0.028);
		legndz->AddEntry(newhfinal2_3b, "Real Data", "f");
		legndz->AddEntry(fsig4,"#eta' Signal Fit","l");
		legndz->AddEntry(fback,"Background Fit","l");
      TString LegZ;
      LegZ.Form("%.2f#sigma Signal Region", sigmaRANGE);
		legndz->AddEntry(zerolineA21, LegZ, "l");
		legndz->Draw("same");

		c3->cd(3);
		newhfinal2_3a->GetXaxis()->SetRangeUser(0.74, 1.18);
		newhfinal2_3a->GetXaxis()->SetRangeUser(0.74, 1.18);
		newhfinal2_3a->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
		newhfinal2_3a->GetXaxis()->SetTitleOffset(1.2);
		newhfinal2_3a->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
		newhfinal2_3a->GetYaxis()->SetTitleOffset(1.5);
		newhfinal2_3a->SetTitle("Background in #eta' Region");
		fback->SetFillColor(8);
		newhfinal2_3a->Add(newhfinal2_5b, -2);
		double mini3a = fsig3->Eval(0.75);
		newhfinal2_3a->SetMaximum(mini3a*1.1);
		mini3a *= -0.08;
		newhfinal2_3a->SetMinimum(mini3a);
		newhfinal2_3a->Draw("hist");
		fback->Draw("same");
      fsig3->Draw("same");

		TLine *zeroline32 = new TLine(0.70, 0.0, 1.18, 0.0);
		zeroline32->SetLineColor(1);
		zeroline32->SetLineStyle(10);
		zeroline32->SetLineWidth(2);
		zeroline32->Draw("same");

		c3->cd(4);
		gStyle->SetOptFit(1100);
		newhfinal2_3c->Add(newhfinal2_6, 10);
		newhfinal2_3c->Draw("hist");
		newhfinal2_3c->GetXaxis()->SetRangeUser(0.74, 1.18);
		newhfinal2_3c->SetTitle("#eta' Signal Fit");
		newhfinal2_3c->SetMinimum(0.0);
		double gryf = fsig4->Eval(mean_sig4);
		newhfinal2_3c->SetMaximum(gryf*1.05);
		fsig4->Draw("same");

		TLine *zerolineD2 = new TLine(GausMax4, 0.0, GausMax4, 1.03*gryf);
		zerolineD2->SetLineColor(6);
		zerolineD2->SetLineStyle(9);
		zerolineD2->SetLineWidth(2);
		zerolineD2->Draw("same");

		double etapdiff = (abs(((GausMax4*1000) - 957.78) / 957.78))*100;
		TString titleetap;
		TString nameetap;
		nameetap.Form("#eta' Fit Peak = %.1f MeV/c^{2}", GausMax4*1000);
		titleetap.Form("Difference = %.2f %%", etapdiff);

		TLegend *legndAC = new TLegend(0.6, 0.68, 0.89, 0.89);
		legndAC->SetTextSize(0.028);
		legndAC->AddEntry(fsig4, "#eta' Signal Fit", "l");
		legndAC->AddEntry(zerolineD2, nameetap, "l");
		legndAC->AddEntry((TObject*)0, "PDG Value = 957.8 MeV/c^{2}", "");
		legndAC->AddEntry((TObject*)0, titleetap, "");
		legndAC->Draw("same");

		sprintf(histo1,"%s%d%s","Etap_mass", totalhisto, ".png");
		c3->Print(histo1);

      c10->cd();
      gStyle->SetOptFit(0);
      c10->SetLogy(1);
      sprintf(histo1,"%s%d%s","mass_", totalhisto, ".png");
      newhfinal2_9->Draw();
      newhfinal2_9->GetXaxis()->SetRangeUser(0.0, 1.2);
      newhfinal2_9->SetTitle(TitleX);
      // Pi0 SB1
      int P1bin1 = hfillClone1->FindBin(PiLower);
      int P1bin2 = hfillClone1->FindBin(PiLow);

      // Pi0 Signal
      int P2bin1 = hfillClone2->FindBin(Pidown);
      int P2bin2 = hfillClone2->FindBin(Piup);

      // Pi0 SB2
      int P3bin1 = hfillClone3->FindBin(PiHigh);
      int P3bin2 = hfillClone3->FindBin(PiHigher);

      // Eta SB1
      int Eta1bin1 = hfillClone4->FindBin(EtaLower);
      int Eta1bin2 = hfillClone4->FindBin(EtaLow);

      // Eta Signal
      int Eta2bin1 = hfillClone5->FindBin(Etadown);
      int Eta2bin2 = hfillClone5->FindBin(Etaup);

      // Eta SB2
      int Eta3bin1 = hfillClone6->FindBin(EtaHigh);
      int Eta3bin2 = hfillClone6->FindBin(EtaHigher);

      // Etap SB1
      int Etap1bin1 = hfillClone7->FindBin(EtapLower);
      int Etap1bin2 = hfillClone7->FindBin(EtapLow);

      // Etap Signal
      int Etap2bin1 = hfillClone8->FindBin(Etapdown);
      int Etap2bin2 = hfillClone8->FindBin(Etapup);

      // Etap SB2
      int Etap3bin1 = hfillClone9->FindBin(EtapHigh);
      int Etap3bin2 = hfillClone9->FindBin(EtapHigher);

      // Set Fill Color
      hfillClone1->SetFillColor(20);
      hfillClone2->SetFillColor(38);
      hfillClone3->SetFillColor(20);

      hfillClone4->SetFillColor(20);
      hfillClone5->SetFillColor(38);
      hfillClone6->SetFillColor(20);

      hfillClone7->SetFillColor(20);
      hfillClone8->SetFillColor(38);
      hfillClone9->SetFillColor(20);

      // Set Line Color
      hfillClone1->SetLineColor(0);
      hfillClone2->SetLineColor(0);
      hfillClone3->SetLineColor(0);
      hfillClone4->SetLineColor(0);
      hfillClone5->SetLineColor(0);
      hfillClone6->SetLineColor(0);
      hfillClone7->SetLineColor(0);
      hfillClone8->SetLineColor(0);
      hfillClone9->SetLineColor(0);

      // Adjust range of colored section to show
      hfillClone1->GetXaxis()->SetRange(P1bin1, P1bin2);        // bin range
      hfillClone2->GetXaxis()->SetRange(P2bin1, P2bin2);        // bin range
      hfillClone3->GetXaxis()->SetRange(P3bin1, P3bin2);        // bin range
      hfillClone4->GetXaxis()->SetRange(Eta1bin1, Eta1bin2);    // bin range
      hfillClone5->GetXaxis()->SetRange(Eta2bin1, Eta2bin2);    // bin range
      hfillClone6->GetXaxis()->SetRange(Eta3bin1, Eta3bin2);    // bin range
      hfillClone7->GetXaxis()->SetRange(Etap1bin1, Etap1bin2);  // bin range
      hfillClone8->GetXaxis()->SetRange(Etap2bin1, Etap2bin2);  // bin range
      hfillClone9->GetXaxis()->SetRange(Etap3bin1, Etap3bin2);  // bin range

      // Draw Colored histo sections onto full-range histo
      // Need to use "hist" in order to use fill color
      hfillClone1->Draw("hist same");
      hfillClone2->Draw("hist same");
      hfillClone3->Draw("hist same");
      hfillClone4->Draw("hist same");
      hfillClone5->Draw("hist same");
      hfillClone6->Draw("hist same");
      hfillClone7->Draw("hist same");
      hfillClone8->Draw("hist same");
      hfillClone9->Draw("hist same");
      ftot->Draw("same");
      fback->Draw("same");
      c10->Print(histo1);
      c10->SetLogy(0);

      // Plot of all fits
      c11->cd();
      // Used to look closer at the background fit
      // 1.1* gives an easier to read plot
      newhfinal2_8->GetXaxis()->SetRangeUser(0.6, 0.7);
      double TotV = newhfinal2_8->GetMaximum();
      double TotMax  =  1.05 * TotV;
      double TotMin  = -0.05 * TotV;

      newhfinal2_8->GetXaxis()->SetRangeUser(0.0, 1.2);
      newhfinal2_8->Draw();
      newhfinal2_8->SetMaximum(TotMax);
      newhfinal2_8->SetMinimum(TotMin);
      newhfinal2_8->SetTitle(TitleX);
      fback->Draw("same");
      newhfinal2_8->SetLineColor(18);
      fsig2->SetLineColor(kYellow);
      fsig4->SetLineColor(7);
      TLine *zline = new TLine(0.0, 0.0, 1.2, 0.0);
      zline->SetLineColor(46);
      zline->SetLineStyle(9);
      zline->SetLineWidth(1);
      zline->Draw("same");
      fsig1->Draw("same");
      fsig2->Draw("same");
      fsig3->Draw("same");
      fsig4->Draw("same");
      ftot->Draw("same");
      TString histo10;
      histo10.Form("M2g_E%d_AllSignals.png", Echannel);
      c11->Print(histo10);

		////////////////////////////////////////////////////////////////// 

		TF1 *fitresultTot = remass->GetFunction("ftot");

		double chisqTot = fitresultTot->GetChisquare();
		double ndofTot = fitresultTot->GetNDF();
		double probTot = fitresultTot->GetProb();
		int nbins = remass->GetNbinsX();

		std::cout << std::endl;
		std::cout << "norm -> " << norm << std::endl;
		std::cout << "bin width -> " << binw << std::endl;

		std::cout << std::endl;
		std::cout << "chisqTot -> " << chisqTot << std::endl;
		std::cout << "ndofTot -> " << ndofTot << std::endl;
		std::cout << "fit probabilityTot -> " << probTot << std::endl; 


      TString BkNEW;
		if (Echannel == 0) {
	      BkNEW.Form("Polynomial Background Parameters for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range0b, Range0a);
		}
		else if (Echannel == 1) {
	      BkNEW.Form("Polynomial Background Parameters for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range1b, Range1a);
		}
		else if (Echannel == 2) {
	      BkNEW.Form("Polynomial Background Parameters for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range2b, Range2a);
		}
		else {
         BkNEW.Form("Polynomial Background Parameters for: %.3f #geq E_{#gamma} < %.3f GeV  &  0.1 #leq |t| #leq 1.050 GeV^{2}/c^{4};[ SB1 | SIGNAL | SB2 ] Regions;M_{2#gamma} [GeV/c^{2}]", Range3b, Range3a);
		}

      TH1D *BkgdRANGE = new TH1D("BkgdRANGE", BkNEW, 5, 0.0, 5.0);
      //BkgdRANGE->Sumw2();

      // BkgdPar[0]
      double bkpar0 = ftot->GetParameter(21);
      // BkgdPar[1]
      double bkpar1 = ftot->GetParameter(22);
      // BkgdPar[2]
      double bkpar2 = ftot->GetParameter(23);
      // BkgdPar[3]
      double bkpar3 = ftot->GetParameter(24);
      // BkgdPar[4]
      double bkpar4 = ftot->GetParameter(25);

      BkgdRANGE->SetBinContent(1, bkpar0);
      BkgdRANGE->SetBinContent(2, bkpar1);
      BkgdRANGE->SetBinContent(3, bkpar2);
      BkgdRANGE->SetBinContent(4, bkpar3);
      BkgdRANGE->SetBinContent(5, bkpar4);

      TString Bkhtitle;
      Bkhtitle.Form("Bkgd_parameters_allt_E%d", Echannel);

		gFile = MyFile;
		gDirectory->WriteObject(BkgdRANGE, Bkhtitle);


		ofstream output;
		output.open("2g_E_binned_Etaprime_fit_allt_TwoG_may_16_2023.txt", std::ios_base::app);
		//// Or you can simply use the following ////
		// ofstream output ("2gamma_results.txt");
		if (output.is_open()) {
			output << "Histogram = " << htitle << std::endl;
         if (Echannel == 0) {
            output << "Tagger Channel Range = " << totalhisto << "  -->  E_photon = " << Range0b << " to " << Range0a << " GeV" << std::endl;
         }
         else if (Echannel == 1) {
            output << "Tagger Channel Range = " << totalhisto << "  -->  E_photon = " << Range1b << " to " << Range1a << " GeV" << std::endl;
         }
         else if (Echannel == 2) {
            output << "Tagger Channel Range = " << totalhisto << "  -->  E_photon = " << Range2b << " to " << Range2a << " GeV" << std::endl;
         }
         else {
            output << "Tagger Channel Range = " << totalhisto << "  -->  E_photon = " << Range3b << " to " << Range3a << " GeV" << std::endl;
         }
			output << "|t| = 0.1 to 1.050 GeV^2/c^4 -- Using XPi_bins[3][2]" << std::endl;
			output << std::endl;

         TString statname;
         if (conv_fit == 1) {
            statname.Form("Converged");
         }
         else {
            statname.Form("Failed");
         }
         output << "Status = " << statname << std::endl;
         output << std::endl;
         output << "chisqTot -> " << chisqTot << std::endl;
         output << "ndofTot -> " << ndofTot << std::endl;
         output << "fit probabilityTot -> " << probTot << std::endl;
         output << std::endl;

			output << "Pi0 fit peak           = " << GausMax1*1000 << " MeV/c^2"  << std::endl;
         output << "PiPar = " << PiPar[1] << std::endl;
			output << "Eta fit peak           = " << GausMax2*1000 << " MeV/c^2"  << std::endl;
         output << "EtaPar = " << EtaPar[1] << std::endl;
			output << "Omega Leakage fit peak = " << GausMax3*1000 << " MeV/c^2"  << std::endl;
         output << "OmegaPar = " << OmegaPar[1] << std::endl;
			output << "Etap fit peak          = " << GausMax4*1000 << " MeV/c^2"  << std::endl;
			output << std::endl;

         output << "Pi0_____mean_1 -> " << ftot->GetParameter(1);
         if (((mean_sig1 - 0.00000001) < dparlim1) || ((mean_sig1 + 0.00000001) > uparlim1)) {
            output << "  ****";
         }
         output << std::endl;

         output << "Pi0____sigma_1 -> " << ftot->GetParameter(2);
         if (((sigma_sig1 - 0.00000001) < dparlim2) || ((sigma_sig1 + 0.00000001) > uparlim2)) {
            output << "  ****";
         }
         output << std::endl;

         output << "Pi0_____mean_2 -> " << ftot->GetParameter(4);
         if (((mean2_sig1 - 0.00000001) < dparlim4) || ((mean2_sig1 + 0.00000001) > uparlim4)) {
            output << "  ****";
         }
         output << std::endl;

         output << "Pi0____sigma_2 -> " << ftot->GetParameter(5);
         if (((sigma2_sig1 - 0.00000001) < dparlim5) || ((sigma2_sig1 + 0.00000001) > uparlim5)) {
            output << "  ****";
         }
         output << std::endl;
         output << std::endl;

         output << "Eta_____mean_1 -> " << ftot->GetParameter(7);
         if (((mean_sig2 - 0.00000001) < dparlim7) || ((mean_sig2 + 0.00000001) > uparlim7)) {
            output << "  ****";
         }
         output << std::endl;

         output << "Eta____sigma_1 -> " << ftot->GetParameter(8);
         if (((sigma_sig2 - 0.00000001) < dparlim8) || ((sigma_sig2 + 0.00000001) > uparlim8)) {
            output << "  ****";
         }
         output << std::endl;

         output << "Eta_____mean_2 -> " << ftot->GetParameter(10);
         if (((mean2_sig2 - 0.00000001) < dparlim10) || ((mean2_sig2 + 0.00000001) > uparlim10)) {
            output << "  ****";
         }
         output << std::endl;

         output << "Eta____sigma_2 -> " << ftot->GetParameter(11);
         if (((sigma2_sig2 - 0.00000001) < dparlim11) || ((sigma2_sig2 + 0.00000001) > uparlim11)) {
            output << "  ****";
         }
         output << std::endl;
         output << std::endl;

         output << "3G_leak_mean_1 -> " << ftot->GetParameter(13);
         if (((mean_sig3 - 0.00000001) < dparlim13) || ((mean_sig3 + 0.00000001) > uparlim13)) {
            output << "  ****";
         }
         output << std::endl;
         output << "3G_leak_sigma1 -> " << ftot->GetParameter(14);
         if (((sigma_sig3 - 0.00000001) < dparlim14) || ((sigma_sig3 + 0.00000001) > uparlim14)) {
            output << "  ****";
         }
         output << std::endl;
         output << "3G_leak_mean_2 -> " << ftot->GetParameter(16);
         if (((mean2_sig3 - 0.00000001) < dparlim16) || ((mean2_sig3 + 0.00000001) > uparlim16)) {
            output << "  ****";
         }
         output << std::endl;
         output << "3G_leak_sigma2 -> " << ftot->GetParameter(17);
         if (((sigma2_sig3 - 0.00000001) < dparlim17) || ((sigma2_sig3 + 0.00000001) > uparlim17)) {
            output << "  ****";
         }
         output << std::endl;
         output << std::endl;

         output << "Etap____mean_1 -> " << ftot->GetParameter(19);
         if (((mean_sig4 - 0.00000001) < dparlim19) || ((mean_sig4 + 0.00000001) > uparlim19)) {
            output << "  ****";
         }
         output << std::endl;
         output << "Etap___sigma_1 -> " << ftot->GetParameter(20);
         if (((sigma_sig4 - 0.00000001) < dparlim20) || ((sigma_sig4 + 0.00000001) > uparlim20)) {
            output << "  ****";
         }
         output << std::endl;
         output << std::endl;
         output << "p0     -> " << ftot->GetParameter(21) << std::endl;
         output << "p1     -> " << ftot->GetParameter(22) << std::endl;
         output << "p2     -> " << ftot->GetParameter(23) << std::endl;
         output << "p3     -> " << ftot->GetParameter(24) << std::endl;

			output << std::endl;
			output << "   // Energy Rescaling Factor" << std::endl;
			output << "ER = 0.97;" << std::endl; 
			output << "0.1 <= |t| <= 1.050 GeV^2/c^4" << std::endl;
			output << std::endl;

			output << "   // Tagging Scaling Factor" << std::endl;
			output << "TSF = 0.935198 +/- 0.00249013;" << std::endl; 
			output << std::endl;

			output << "bin width -> " << binw << std::endl;
			output << std::endl;

         output << "EtapLower  = " << EtapLower << std::endl;
         output << "EtapLow    = " << EtapLow << std::endl;
         output << "Etapdown   = " << Etapdown << std::endl;
         output << "Etapup     = " << Etapup << std::endl;
         output << "EtapHigh   = " << EtapHigh << std::endl;
         output << "EtapHigher = " << EtapHigher << std::endl;

			output << "_____________________________________________" << std::endl;
			output << std::endl;
		}
		else {
			std::cout << "UNABLE TO OPEN OUTPUT FILE." << std::endl;
		}
		output.close();

      MyFile->Close();

   }

   if (conv_fit == 1) {
      std::cout << "All Fits Finished" << std::endl;
      delete c1;
      delete c2;
      delete c3;
      delete c10;
      delete c11;
      delete c123;
      delete c321;
      delete c555;

      // Organize folder info
      gSystem->Exec("mkdir -p histos/Pi0");
      gSystem->Exec("mkdir -p histos/Eta");
      gSystem->Exec("mkdir -p histos/Etap");
      gSystem->Exec("mkdir -p histos/LogPlot");
      gSystem->Exec("mkdir -p histos/AllSignals");
      gSystem->Exec("mkdir -p histos/Fits");

      gSystem->Exec("mv PiFit_*.png histos/Fits/.");
      gSystem->Exec("mv EtaFit_*.png histos/Fits/.");
      gSystem->Exec("mv EtapFit_*.png histos/Fits/.");
      gSystem->Exec("mv Pi_*.png histos/Pi0/.");
      gSystem->Exec("mv Eta_*.png histos/Eta/.");
      gSystem->Exec("mv Etap_*.png histos/Etap/.");
      gSystem->Exec("mv mass_*.png histos/LogPlot/.");
      gSystem->Exec("mv *_AllSignals.png histos/AllSignals/.");
   }
}
