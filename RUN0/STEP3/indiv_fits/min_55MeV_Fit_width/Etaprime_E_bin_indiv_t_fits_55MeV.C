///////////////////////////////////////////////////
// Author: J. McIntyre
// 16 May 2023
///////////////////////////////////////////////////
// This file is used for fitting the individual |t| & E_photon
// binned histograms with basic cuts and tag weighting applied.
// The non-binned (in |t|) histogram was used to determine the 
// background model (STEP2). This will be used to fit the 
// background in each individual |t| & E_photon binned histogram. 
// A scaling factor coefficient is added to the fit model from 
// STEP2 to scale the background fit in STEP3. The background 
// parameters determined in STEP2, for photon energy binned over 
// all |t|, are fixed in STEP3.
//
// The background is modeled as a pol3 plus a double gaussian
// for the Omega 3 gamma leakage.
//
// ------------------------------------------------

/*
   //----------------------
   // Add line to histogram
   //----------------------
   TLine *zeroline3 = new TLine(0.74, 0, 1.18, 0);
   zeroline3->SetLineColor(1);
   zeroline3->SetLineStyle(10);
   zeroline3->SetLineWidth(2);
   zeroline3->Draw("hist");

   //----------------
   // Drawing Options
   //----------------
   gStyle->SetOptStat(0);                 // No Legend
   gStyle->SetOptTitle(0);                // No title
   histo->SetLineWidth(5);                // Line thickness
   histo->SetLineStyle(5);                // Line type (1=solid)
   histo->SetLineColor(kBlue);            // Line color
   histo->SetFillColor(kBlue);            // Filling color to histo (38 = light blue)
   histo->SetFillColorAlpha(kBlue,0.25);  // Suppose to be semi-transparent
   //----------------
*/

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
#include <TMinuit.h>


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
TH1D *EChanPARS;
TH1D *EEtaRange;
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


void Etaprime_E_bin_indiv_t_fits_55MeV() {

   // Create a Root file to later fill with new histograms for reference in further analysis work
   TFile *StarterFile = new TFile("TwoG_may_16_2023_STEP3_55MeV.root", "NEW");
   StarterFile->Close();

   // File from STEP2 containing information for calculations & fitting
   TFile *fStep2 = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP2/total_fit/min_55MeV_Fit_width/TwoG_may_16_2023_STEP2_55MeV.root");

   // Importing |t| bin values
   TH1D *hBinDW  = (TH1D*)fStep2->Get("hBinDW");
   TH1D *hBinUP  = (TH1D*)fStep2->Get("hBinUP");

   int tBinNum = hBinUP->GetNbinsX();
   double Etap_bins[tBinNum][2];
   for (int b1 = 0; b1 < tBinNum; b1++) {
      int bs = b1 + 1;
      Etap_bins[b1][0] = hBinDW->GetBinContent(bs);
      Etap_bins[b1][1] = hBinUP->GetBinContent(bs);

      // For troubleshooting
      //std::cout << "t bin[0] = " << Etap_bins[b1][0] << "\n";
      //std::cout << "t bin[1] = " << Etap_bins[b1][1] << "\n";
   }

   // Import number of photon energy bins used
   TH1D *EChanBINS = (TH1D*)fStep2->Get("EChanBINS");
   int NumEChan = EChanBINS->GetBinContent(1);
   // For troubleshooting
   //std::cout << "NumEChan = " << NumEChan << "\n";

   // Import number of polynomial parameters used 
   // We know there will always be at least one photon energy bin (e.g. "E0" histogram)
   TH1D *PolParam = (TH1D*)fStep2->Get("Bkgd_parameters_allt_E0");
   int PolParNum = PolParam->GetNbinsX();
   // For troubleshooting
   //std::cout << "PolParNum = " << PolParNum << "\n";

   // Import polynomial parameters used to fit all |t| histos in STEP2
   // Values will be fixed for fitting in this step
   TString echanname;
   double EbinBkgdPar[NumEChan][PolParNum];
   for (int i = 0; i < NumEChan; i++) {
      echanname.Form("Bkgd_parameters_allt_E%d", i);
      EChanPARS = (TH1D*)fStep2->Get(echanname);
      for (int j = 0; j < PolParNum; j++) {
         EbinBkgdPar[i][j] = EChanPARS->GetBinContent(j + 1);
         // For troubleshooting
         //std::cout << "EbinBkgdPar = " << EbinBkgdPar[i][j] << "\n";
      }
   }

   // Import values for Neutral Pion Sideband Subtraction Regions
   // I use the sideband subtraction regions determined from fitting
   // over all |t| in each tagger range in STEP2
   TString rname;
   double RangeUp[NumEChan];
   double RangeDown[NumEChan];
   double RangeLow[NumEChan];
   double RangeLower[NumEChan];
   double RangeHigh[NumEChan];
   double RangeHigher[NumEChan];

   for (int k = 0; k < NumEChan; k++) {
      rname.Form("SB_Regions_allt_E%d", k);
      EEtaRange = (TH1D*)fStep2->Get(rname);
      RangeLower[k]  = EEtaRange->GetBinContent(1);
      RangeLow[k]    = EEtaRange->GetBinContent(2);
      RangeDown[k]   = EEtaRange->GetBinContent(3);
      RangeUp[k]     = EEtaRange->GetBinContent(4);
      RangeHigh[k]   = EEtaRange->GetBinContent(5);
      RangeHigher[k] = EEtaRange->GetBinContent(6);

      // For troubleshooting
      //std::cout << "E channel = " << k << "\n";
      //std::cout << "RangeLower = " << RangeLower[k] << "\n";
      //std::cout << "RangeLow = " << RangeLow[k] << "\n";
      //std::cout << "RangeDown = " << RangeDown[k] << "\n";
      //std::cout << "RangeUp = " << RangeUp[k] << "\n";
      //std::cout << "RangeHigh = " << RangeHigh[k] << "\n";
      //std::cout << "RangeHigher = " << RangeHigher[k] << "\n";
   }

   // Import values for Tagger Energy Bins/Channels
   double RangeA[NumEChan];
   double RangeB[NumEChan];
   int TC_RangeA[NumEChan];
   int TC_RangeB[NumEChan];

   TH1D *EChanRange = (TH1D*)fStep2->Get("EChanRange");
   TH1D *EChanValue = (TH1D*)fStep2->Get("EChanValue");

   for (int m = 0; m < NumEChan; m++) {
      int zz = m * 2 + 1;
		// Tagger channels used for photon energy binning
      TC_RangeA[m] = EChanRange->GetBinContent(zz);
		// Photon energy value used for photon energy binning
      RangeA[m]    = EChanValue->GetBinContent(zz);

      zz += 1;
      TC_RangeB[m] = EChanRange->GetBinContent(zz);
		RangeB[m]    = EChanValue->GetBinContent(zz);
   }

   // All values imported from STEP2 root file, close TFile
   fStep2->Close();

   // Array used to store Etaprime fit yields
   int histoCount = NumEChan * tBinNum;
   double EtapYield[histoCount];
   double EtapYieldError[histoCount];

   // Real data root file from STEP1 after cleanup cuts and accidental subtraction
   TFile *rootf = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP1/inv_mass/TwoG_may_16_2023Analysis.root");

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

   // Declaring variable for fit convergence (outside loop),
   // so that TCanvases aren't deleted if a fit fails.
   // Undeleted histos on screen will allows easier troubleshooting.   
   int conv_fit = 1;

   // Variables used to set fit limits.
   // Defined initially outside the nested loops,
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

   // A variable used to only print certain text
   // once in 2g_Etaprime_Results_55MeV.txt
   int reset = 0;

   // Setup an array of weights to be entered into a TH1D histogram 
   // at the end of each Echannel loop. 
   double w_etap_sigW[tBinNum];
   double w0_errW[tBinNum];
   double w1_etap_sbW[tBinNum];
   double w1_errW[tBinNum];
   double w2_etap_sbW[tBinNum];
   double w2_errW[tBinNum];

   // Troubleshooting
   double wR1_err[tBinNum];
   double wR2_err[tBinNum];
   double wC1_err[tBinNum];
   double wC2_err[tBinNum];
   double wR1_Etap[tBinNum];
   double wR2_Etap[tBinNum];
   double wC1_Etap[tBinNum];
   double wC2_Etap[tBinNum];
   double wiEtap_s0[tBinNum];
   double wEtapGaus_s0[tBinNum];
   double wiEtap_s1[tBinNum];
   double wEtapGaus_s1[tBinNum];
   double wiEtap_s2[tBinNum];
   double wEtapGaus_s2[tBinNum];
   double wiEtap_b0[tBinNum];
   double wEtapGaus_b0[tBinNum];
   double wiEtap_b1[tBinNum];
   double wEtapGaus_b1[tBinNum];
   double wiEtap_b2[tBinNum];
   double wEtapGaus_b2[tBinNum];

   TString htitle;
   int totalhisto = 0;

   // Boundaries used in calculating Sideband Subtraction Weights (set later) 
   double XEtapUp     = 0.0;
   double XEtapDown   = 0.0;
   double XEtapLow    = 0.0;
   double XEtapLower  = 0.0;
   double XEtapHigh   = 0.0;
   double XEtapHigher = 0.0;

   // For troubleshooting
   int Echannel_count = 0;
   int tchannel_count = 0;

   // USED TO IGNORE FITS WITH NO ENTRIES
   bool ZeroFit = false;
   double ZeroFitArray[tBinNum];
   double ZFA[histoCount];

   for (int Echannel = 0; Echannel < NumEChan; Echannel++) {
      // For troubleshooting
      Echannel_count = Echannel;

      // Background polynomial parameters determined 
      // in STEP2 for each photon energy bin
		double BkgdPar[5];
      for (int ppar = 0; ppar < PolParNum; ppar++) {
         BkgdPar[ppar] = EbinBkgdPar[Echannel][ppar];
      }

      // Use Sideband Subtraction boundary values from STEP2 (over all |t| for a photon energy bin)
		XEtapUp     = RangeUp[Echannel];
		XEtapDown   = RangeDown[Echannel];
		XEtapLow    = RangeLow[Echannel];
		XEtapLower  = RangeLower[Echannel];
		XEtapHigh   = RangeHigh[Echannel];
		XEtapHigher = RangeHigher[Echannel];

      for (int tchannel = 0; tchannel < tBinNum; tchannel++) {
         // USED TO IGNORE FITS WITH NO ENTRIES
         ZeroFit = false;

         // For troubleshooting
         tchannel_count = tchannel;
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

         totalhisto = (Echannel * tBinNum) + tchannel;
         htitle.Form("M2g_E%d_%d", Echannel, tchannel);
         remass = (TH1D*)rootf->Get(htitle);

         std::cout << "*********** Starting ***********" << "\n";
         std::cout << "          " << htitle << "          " << "\n";
         std::cout << "********************************" << "\n";

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

         // Switch to Minuit 
         ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");


         int bv0 = 0;
         double binval0 = 0.0;
         double xvalue = 0.0;
         bool Adjustment = false;
         // Way to fix background to a data point, when the fitter needs help
         if (Echannel == 0 && tchannel == 0) {
            // E0_t0
            xvalue = 1.11; // GeV/c^2
            bv0 = xvalue / binw;
            binval0 = remass->GetBinContent(bv0);
            Adjustment = true;
         }

         else if (Echannel == 1 && tchannel == 0) {
            // E1_t0; 1100 MeV/c^2
            xvalue = 1.1; // GeV/c^2
            bv0 = xvalue / binw;
            binval0 = remass->GetBinContent(bv0);
            Adjustment = true;
         }

         else if (Echannel == 2 && tchannel == 0) {
            // E1_t2; 935 MeV/c^2
            xvalue = 1.135; // GeV/c^2
            bv0 = xvalue / binw;
            binval0 = remass->GetBinContent(bv0);
            Adjustment = true;
         }

         else if (Echannel == 3 && tchannel == 0) {
            // E2_t0; 1120 MeV/c^2
            xvalue = 1.12; // GeV/c^2
            bv0 = xvalue / binw;
            binval0 = remass->GetBinContent(bv0);
            Adjustment = true;
         }
         /*
         else if (Echannel == 3 && tchannel == 0) {
            // E3_t0; 1150 MeV/c^2
            xvalue = 1.155; // GeV/c^2
            bv0 = xvalue / binw;
            binval0 = remass->GetBinContent(bv0);
            Adjustment = true;
         }

         else if (Echannel == 3 && tchannel == 2) {
            // E3_t2; 980 MeV/c^2
            xvalue = 0.98; // GeV/c^2
            bv0 = xvalue / binw;
            binval0 = remass->GetBinContent(bv0);
            Adjustment = true;
         }
         */
         // Background value at xvalue
         double bkground = BkgdPar[0] +
                           BkgdPar[1] * xvalue +
                           BkgdPar[2] * pow(xvalue, 2.0) +
                           BkgdPar[3] * pow(xvalue, 3.0);

         double bkpoly3 = std::abs(binval0 / bkground);

         if (!Adjustment) {
            bkpoly3 = BkgdPar[4];
         }


         // Used to fit a Single Gaussian in the range of the Pi0, Eta, and Omega
         // These parameters will help with the fitting as the invariant mass walks in |t|
			TH1D *SG_Pi    = (TH1D*)remass->Clone("SG_Pi");
			TH1D *SG_Eta   = (TH1D*)remass->Clone("SG_Eta");
			TH1D *SG_Omega = (TH1D*)remass->Clone("SG_Omega");
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
			fsig1 = new TF1("fsig1", fsignal1, 0.0, 0.3, 6);      // fsignal1, 0.05, 0.3, 6)
			fsig2 = new TF1("fsig2", fsignal2, 0.3, 0.8, 6);      // fsignal2, 0.3, 0.8, 6); 
			fsig3 = new TF1("fsig3", fsignal3, 0.32, 1.15, 6);    // Leakage // fsignal3, 0.5, 1.3, 6);
			fsig4 = new TF1("fsig4", fsignal4, 0.70, 1.20, 3);    // fsignal4, 0.72, 1.3, 3);

			// Background. The 3 gamma leakage (a.k.a. background) is kept separate (fsig3).
			fback = new TF1("fback", fbackground, 0.0, 1.28, 5);   // 3th Order Polynomial

			// Total Signal
			ftot = new TF1("ftot", ftotal, 0.0, 1.28, 26);

         double par[26];
         TFitResultPtr fitResult;
         int BestLoop = 0;
         double fudgeFactor = 0.002;

         while (BestLoop != 1) {
            fudgeFactor += 0.001;
            double parC[26] = {sqrt(norm * 0.639753), (PiPar[1] - fudgeFactor), (PiPar[2] - fudgeFactor),
                              sqrt(norm * 0.373938), (PiPar[1] + fudgeFactor), (PiPar[2] + fudgeFactor),
                              sqrt(norm * 0.0680778), (EtaPar[1] - fudgeFactor), (EtaPar[2] - fudgeFactor),
                              sqrt(norm * 0.0225913), (EtaPar[1] + fudgeFactor), (EtaPar[2] + fudgeFactor),
                              sqrt(norm * 0.00461673), 0.728, 0.045,  // Avoid using Single Gaussian values since Omega disappears at high |t|
                              sqrt(norm * 0.00433187), 0.742, 0.055,
                              sqrt(norm * 0.00119216), 0.9574, 0.0514,
                              BkgdPar[0], BkgdPar[1], BkgdPar[2], BkgdPar[3], bkpoly3};

            for (int i = 0; i < 26; i++) {
               par[i] = parC[i];
            }

            ftot->SetParameters(par);

            // Helps the fits not to run off
            for (int i = 0; i < 26; i++) {
               ftot->FixParameter(i, par[i]);
            }

            remass->Fit("ftot", "RB");

            // PAR[21] THROUGH PAR[24] ARE IMPORTED FROM PREVIOUS STEP (ALL |t| FIT)
            for (int i = 0; i < 21; i++) {
               ftot->ReleaseParameter(i);
            }

            if (!Adjustment) {
               // Bkgd Scaling Factor
               ftot->ReleaseParameter(25);
            }


            // Used to prevent STATUS=CALL LIMIT. Gives fitter more tries.
            ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

            ////////////////////////////////////////////////////////////////////////
            // Total Signal - Limits for means and std dev of the double gaussian //
            ////////////////////////////////////////////////////////////////////////

            //--------------------------------------------------------------------//
            /* If a fit fails, modify the appropriate "if statement" below to fix */
            //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

            if (tchannel >= 18) {
               // Signal #1 - Pi0 mean limits for the double gaussian
               // Gaussian 1
               uparlim1 = PiPar[1] + 0.040;
               dparlim1 = PiPar[1] - 0.010;
               // Gaussian 2
               uparlim4 = PiPar[1] + 0.040;
               dparlim4 = PiPar[1] - 0.010;

               // Signal #1 - Pi0 std dev limits for the double gaussian
               // Gaussian 1
               uparlim2 = 0.085;
               dparlim2 = 0.010;
               // Gaussian 2
               uparlim5 = 0.085;
               dparlim5 = 0.010;
            }
            else if (tchannel <= 1) {
               // Signal #1 - Pi0 mean limits for the double gaussian
               // Gaussian 1
               uparlim1 = PiPar[1] + 0.028;
               dparlim1 = PiPar[1] - 0.008;
               // Gaussian 2
               uparlim4 = PiPar[1] + 0.028;
               dparlim4 = PiPar[1] - 0.008;

               // Signal #1 - Pi0 std dev limits for the double gaussian
               // Gaussian 1
               uparlim2 = 0.035;
               dparlim2 = 0.010;
               // Gaussian 2
               uparlim5 = 0.035;
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
               uparlim2 = 0.055;
               dparlim2 = 0.010;
               // Gaussian 2
               uparlim5 = 0.055;
               dparlim5 = 0.010;
            }

            if (Echannel == 0 && tchannel == 1) {
               // Signal #2 - Eta mean limits for the double gaussian
               // Gaussian 1
               uparlim7 = EtaPar[1] + 0.045;
               dparlim7 = EtaPar[1] - 0.010;
               // Gaussian 2
               uparlim10 = EtaPar[1] + 0.045;
               dparlim10 = EtaPar[1] - 0.010;

               // Signal #2 - Eta std dev limits for the double gaussian
               // Gaussian 1
               uparlim8 = 0.075;
               dparlim8 = 0.010;
               // Gaussian 2
               uparlim11 = 0.075;
               dparlim11 = 0.010;
            }
            else if (Echannel == 2 && tchannel == 17) {
               // Signal #2 - Eta mean limits for the double gaussian
               // Gaussian 1
               uparlim7 = EtaPar[1] + 0.045;
               dparlim7 = EtaPar[1] - 0.010;
               // Gaussian 2
               uparlim10 = EtaPar[1] + 0.045;
               dparlim10 = EtaPar[1] - 0.010;

               // Signal #2 - Eta std dev limits for the double gaussian
               // Gaussian 1
               uparlim8 = 0.075;
               dparlim8 = 0.025;
               // Gaussian 2
               uparlim11 = 0.075;
               dparlim11 = 0.025;
            }
            else {
               // Signal #2 - Eta mean limits for the double gaussian
               // Gaussian 1
               uparlim7 = EtaPar[1] + 0.050;
               dparlim7 = EtaPar[1] - 0.015;
               // Gaussian 2
               uparlim10 = EtaPar[1] + 0.050;
               dparlim10 = EtaPar[1] - 0.015;

               // Signal #2 - Eta std dev limits for the double gaussian
               // Gaussian 1
               uparlim8 = 0.075;
               dparlim8 = 0.015;
               // Gaussian 2
               uparlim11 = 0.075;
               dparlim11 = 0.015;
            }

            // Signal #3 - Omega Leakage mean limits for the double gaussian
            // Gaussian 1
            uparlim13 = 0.758;
            dparlim13 = 0.708;
            // Gaussian 2
            uparlim16 = 0.758;
            dparlim16 = 0.708;

            // Signal #3 - Omega Leakage std dev limits for the double gaussian
            // Gaussian 1
            uparlim14 = 0.075;
            dparlim14 = 0.025;
            // Gaussian 2
            uparlim17 = 0.075;
            dparlim17 = 0.025;

            // Signal #4 - Etap mean limits for the single gaussian
            // Single Gaussian
            uparlim19 = 0.990;
            dparlim19 = 0.940;

            // Signal #4 - Etap std dev limit for the single gaussian
            // Single Gaussian
            uparlim20 = 0.08;
            dparlim20 = 0.055;

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

            fitResult = remass->Fit("ftot", "SBRM");
            BestLoop = fitResult->IsValid();
         }

         auto covMatrix = fitResult->GetCovarianceMatrix();
         //std::cout << "Covariance matrix from the fit ";
         //covMatrix.Print();
         conv_fit = fitResult->IsValid();
         double *p = ftot->GetParameters();
         const double *p_err = ftot->GetParErrors();

         std::cout << "**********************" << "\n";
         std::cout << conv_fit << "\n";
         std::cout << "**********************" << "\n";

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
			double newpar[26];
			double newerror[26];
			ftot->GetParameters(newpar);
			ftot->SetLineColor(kRed);
         for (int pe = 0; pe < 26; pe++) {
            newerror[pe] = ftot->GetParError(pe);
         }
			fsig1->SetParameters(newpar);
			fsig1->SetParErrors(newerror);
			fsig1->SetLineColor(kBlue);
			fsig2->SetParameters(&newpar[6]);
			fsig2->SetParErrors(&newerror[6]);
			fsig2->SetLineColor(kBlue);
			fsig3->SetParameters(&newpar[12]);
			fsig3->SetParErrors(&newerror[12]);
			fsig3->SetLineColor(kMagenta);
			fsig4->SetParameters(&newpar[18]);
			fsig4->SetParErrors(&newerror[18]);
			fsig4->SetLineColor(1);
			fback->SetParameters(&newpar[21]);
			fback->SetParErrors(&newerror[21]);
			fback->SetLineColor(kGreen);


			/* Etaprime Signal Yield */
         // These constants make it easier to read the equations below.
         // upperX & lowerX is the range that we fit Etap
         double lowerX = 0.70;
         double upperX = 1.20;

         // Storing Etap fit yield, so it canbe added to a TH1D at the end of this macro
         double YIELD = fsig4->Integral(lowerX, upperX) / binw;
         int hfit = Echannel * tBinNum + tchannel;
         EtapYield[hfit] = YIELD;

         // Do not skip filling histogram at end of macro
         ZeroFitArray[tchannel] = 0.0;
         ZFA[hfit] = 0.0;

         // Ignores fits with events with < 5 entries (basically a null fit)
         if (YIELD <= 5.0) {
            ZeroFit = true;
            YIELD = 0.0;
            // Skip filling histogram at end of macro due to a null fit
            ZeroFitArray[tchannel] = 1.0;
            ZFA[hfit] = 1.0;
         }

         ////////////////////////////////////////////////////
         /* Integral of a Gaussian from lowerX to upperX */
         ////////////////////////////////////////////////////
         double iEtap_yield = -(sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - upperX) / (sqrt(2.0) * p[20]))) +
                               (sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - lowerX) / (sqrt(2.0) * p[20])));
         iEtap_yield /= binw;

         // Single gaussian parameters
         // Partial derivitive of "iEtap_yield" wrt p[18]
         double g18yield =  -sqrt(2.0 * PI) * p[18] * p[20] * TMath::Erf((p[19] - upperX) / (sqrt(2.0) * p[20])) +
                             sqrt(2.0 * PI) * p[18] * p[20] * TMath::Erf((p[19] - lowerX) / (sqrt(2.0) * p[20]));

         // Partial derivitive of "iEtap_yield" wrt p[19]
         double g19yield = pow(p[18], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[19] - lowerX), 2.0) / (2.0 * pow(p[20], 2.0))))) - 
                                              TMath::Power(TMath::E(), (-(pow((p[19] - upperX), 2.0) / (2.0 * pow(p[20], 2.0))))));

         // Partial derivitive of "iEtap_yield" wrt p[20]
         double g20yield = pow(p[18], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[19] - lowerX) / (sqrt(2.0) * p[20]))) -
                                              (sqrt(PI / 2.0) * TMath::Erf((p[19] - upperX) / (sqrt(2.0) * p[20]))) -
                                              ((p[19] - lowerX) / p[20]) * TMath::Power(TMath::E(), (-(pow((p[19] - lowerX), 2.0)) / (2.0 * pow(p[20], 2.0)))) +
                                              ((p[19] - upperX) / p[20]) * TMath::Power(TMath::E(), (-(pow((p[19] - upperX), 2.0)) / (2.0 * pow(p[20], 2.0)))));


         double EtapYield_Error = sqrt(std::fabs(g18yield * g18yield * covMatrix(18, 18) + g19yield * g19yield * covMatrix(19, 19) + g20yield * g20yield * covMatrix(20, 20) +
                                                 2.0 * g18yield * g19yield * covMatrix(18, 19) + 2.0 * g18yield * g20yield * covMatrix(18, 20) +
                                                 2.0 * g19yield * g20yield * covMatrix(19, 20)));
         EtapYield_Error /= binw;

         // Storing Etap fit yield error, so it canbe added to a TH1D at the end of this macro
         EtapYieldError[hfit] = EtapYield_Error;


			//Output file of Etaprime yield
			ofstream CheckYield;
         TString hName;
         hName.Form("E%d_t%d", Echannel, tchannel);
         CheckYield.open("Etaprime_yield_check_55MeV.txt", std::ios_base::app);
         if (CheckYield.is_open()) {
            CheckYield << "---- " << hName << " ----" << "\n";
            CheckYield << "----   Yield   ----" << "\n";
            CheckYield << "Fit Integral   : " << YIELD << "\n";
            CheckYield << "Manual Integral: " << iEtap_yield << " +/- " << EtapYield_Error << "\n";
            CheckYield << "\n";
         }
         else {
            std::cout << "Unable to open fit yield." << "\n";
         }
         CheckYield.close();

			//Output file of Etaprime yield
			ofstream FitYield;
         TString HistName;
         HistName.Form("E%d_t%d", Echannel, tchannel);
         FitYield.open("Etaprime_yield_fit_width_55MeV_min.txt", std::ios_base::app);
         if (FitYield.is_open()) {
            FitYield << HistName << " = " << YIELD << "\n";
         }
         else {
            std::cout << "Unable to open fit yield." << "\n";
         }
         FitYield.close();


         // Location of the fits' maximum
			double GausMax1 = fsig1->GetMaximumX();
			double GausMax2 = fsig2->GetMaximumX();
			double GausMax3 = fsig3->GetMaximumX();
			double GausMax4 = fsig4->GetMaximumX();

         // Helps eliminate plotting fits where the peak value is < 1.0 events
         double PeakVal  = fsig1->GetMaximum();
         double PeakVal2 = fsig2->GetMaximum();
         double PeakVal3 = fsig4->GetMaximum();

         ///////////////////////////////////////////
         //////////////// Real Data ////////////////
         ///////////////////////////////////////////
			double peak_sig1   = std::abs(ftot->GetParameter(0));    // Pi0
			double peak2_sig1  = std::abs(ftot->GetParameter(3));
			double mean_sig1   = ftot->GetParameter(1);
			double mean2_sig1  = ftot->GetParameter(4);
			double sigma_sig1  = ftot->GetParameter(2);
			double sigma2_sig1 = ftot->GetParameter(5);

			double peak_sig2   = std::abs(ftot->GetParameter(6));    // Eta
			double peak2_sig2  = std::abs(ftot->GetParameter(9));
			double mean_sig2   = ftot->GetParameter(7);
			double mean2_sig2  = ftot->GetParameter(10);
			double sigma_sig2  = ftot->GetParameter(8);
			double sigma2_sig2 = ftot->GetParameter(11);

			double peak_sig3   = std::abs(ftot->GetParameter(12));   // Omega Leakage "signal"
			double peak2_sig3  = std::abs(ftot->GetParameter(15));
			double mean_sig3   = ftot->GetParameter(13);
			double mean2_sig3  = ftot->GetParameter(16);
			double sigma_sig3  = ftot->GetParameter(14);
			double sigma2_sig3 = ftot->GetParameter(17);

			double peak_sig4   = std::abs(ftot->GetParameter(18));   // Eta prime
			double mean_sig4   = ftot->GetParameter(19);
			double sigma_sig4  = ftot->GetParameter(20);

         // Weighted average of the fitting double gaussians' mean & sigma
         double w1AvgPi  = peak_sig1  / (peak_sig1 + peak2_sig1);
         double w2AvgPi  = peak2_sig1 / (peak_sig1 + peak2_sig1);
         double mu_Pi    = (w1AvgPi * mean_sig1) + (w2AvgPi * mean2_sig1);
         double sigma_Pi = sqrt(w1AvgPi * (pow(sigma_sig1, 2.0) + pow((mean_sig1 - mu_Pi), 2.0)) + w2AvgPi * (pow(sigma2_sig1, 2.0) + pow((mean2_sig1 - mu_Pi), 2.0)));

         double w1AvgEta  = peak_sig2  / (peak_sig2 + peak2_sig2);
         double w2AvgEta  = peak2_sig2 / (peak_sig2 + peak2_sig2);
         double mu_Eta    = (w1AvgEta * mean_sig2) + (w2AvgEta * mean2_sig2);
         double sigma_Eta = sqrt(w1AvgEta * (pow(sigma_sig2, 2.0) + pow((mean_sig2 - mu_Eta), 2.0)) + w2AvgEta * (pow(sigma2_sig2, 2.0) + pow((mean2_sig2 - mu_Eta), 2.0)));

         double w1AvgOmega  = peak_sig3  / (peak_sig3 + peak2_sig3);
         double w2AvgOmega  = peak2_sig3 / (peak_sig3 + peak2_sig3);
         double mu_Omega    = (w1AvgOmega * mean_sig3) + (w2AvgOmega * mean2_sig3);
         double sigma_Omega = sqrt(w1AvgOmega * (pow(sigma_sig3, 2.0) + pow((mean_sig3 - mu_Omega), 2.0)) + w2AvgOmega * (pow(sigma2_sig3, 2.0) + pow((mean2_sig3 - mu_Omega), 2.0)));

         double TOTALPi   = fsig1->Integral(0.05, 0.30) / binw;
         double TOTALEta  = fsig2->Integral(0.30, 0.80) / binw;
         double TOTALEtap = fsig4->Integral(0.70, 1.20) / binw;

         /*
         // For troubleshooting
         std::cout << "XEtapLower = " << XEtapLower * 1000 << "\n";
         std::cout << "XEtapLow = " << XEtapLow * 1000 << "\n";
         std::cout << "XEtapDown = " << XEtapDown * 1000 << "\n";
         std::cout << "XEtapUp = " << XEtapUp * 1000 << "\n";
         std::cout << "XEtapHigh = " << XEtapHigh * 1000 << "\n";
         std::cout << "XEtapHigher = " << XEtapHigher * 1000 << "\n";
         std::cout << "\n";
         */


			/* Signal only in SIGNAL REGION */
			/* S_0 */
			double Sigcount_Etap = fsig4->Integral(XEtapDown, XEtapUp) / binw;

         // These constants make it easier to read the equations below.
         double highX = XEtapUp;
         double lowX  = XEtapDown;
         ////////////////////////////////////////////////////
         /* Integral of a Gaussian from lowX to highX */
         ////////////////////////////////////////////////////
         double iEtap_s0 = -(sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))) +
                            (sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20])));
         iEtap_s0 /= binw;

         // First gaussian parameters
         // Partial derivitive of "iEtap_s0" wrt p[18]
         double g18s0 =  -sqrt(2.0 * PI) * p[18] * p[20] * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20])) +
                          sqrt(2.0 * PI) * p[18] * p[20] * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20]));

         // Partial derivitive of "iEtap_s0" wrt p[19]
         double g19s0 = pow(p[18], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[19] - lowX), 2.0)  / (2.0 * pow(p[20], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[19] - highX), 2.0) / (2.0 * pow(p[20], 2.0))))));

         // Partial derivitive of "iEtap_s0" wrt p[20]
         double g20s0 = pow(p[18], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))) -
                                           ((p[19] - lowX)  / p[20]) * TMath::Power(TMath::E(), (-(pow((p[19] - lowX), 2.0))  / (2.0 * pow(p[20], 2.0)))) +
                                           ((p[19] - highX) / p[20]) * TMath::Power(TMath::E(), (-(pow((p[19] - highX), 2.0)) / (2.0 * pow(p[20], 2.0)))));


         double EtapGaus_s0 = sqrt(std::fabs(g18s0 * g18s0 * covMatrix(18, 18) + g19s0 * g19s0 * covMatrix(19, 19) + g20s0 * g20s0 * covMatrix(20, 20) +
                                               2.0 * g18s0 * g19s0 * covMatrix(18, 19) + 2.0 * g18s0 * g20s0 * covMatrix(18, 20) +
                                               2.0 * g19s0 * g20s0 * covMatrix(19, 20)));
         EtapGaus_s0 /= binw;

         std::cout << "****************************************************************************" << "\n";
         std::cout << "****************************************************************************" << "\n";
         std::cout << "Eta' Integral = " << Sigcount_Etap << " vs. calc -> " << iEtap_s0 << "\n";
         std::cout << "Error = " << EtapGaus_s0 << "\n";
         std::cout << "****************************************************************************" << "\n";
         std::cout << "****************************************************************************" << "\n";

			///////////////////////////////////////////////////////////
			// Just used for plotting fits onto histograms           //
			// Separate calculation are performed for Pi0 & Eta      //
			// in other TSelector analysis.                          //
			//                                                       //
			//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

         // Calculating the ranges for accidental subtraction weights Pi0 Regions// 
         // Center of Pi Signal Region
         double mx1 = mu_Pi;     // Weighted average of double gaussian fit
         /*
         // Alternate center for signal region
         double mx1 = GausMax1;  // Maximum of the fit fsig1
         */

         // sigmawidth = half of signal region rounded to nearest MeV
         // (1000[MeV/GeV conversion] * selected number of sigma width * weighted average std. dev. of double gaussian fit)
         int sigmawidthPi = lround(1000 * sigmaRANGE * sigma_Pi);

         // See if it is close to the next bin
         int AddBinPi = 0;
         if (sigmawidthPi % binwMeV >= 3) {
            AddBinPi = 1;
         }

         // Get the bin value from the real data histogram
         sigmawidthPi /= binwMeV;
         int sigmaBinPi = sigmawidthPi + AddBinPi;

         // Bin value for the double gaussian fit weighted average mean location
         int pimean = remass->FindBin(mx1);

         //////////////////////////////////////////////////////////////////////////////////////////////////////////////
         /*                                       Invariant Mass Increases -->                                       */
         /*    SB1 Region = (PiLower, PiLow) || Signal Region = (Pidown, Piup) || SB2 Region = (PiHigh, PiHigher)    */
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////
         // X-values (GeV/c^2) to use in integrals of a TF1 fit (don't forget to divide by bin width of real histo)
         // This makes the signal region 1 bin larger than needed
         double Piup   = remass->GetXaxis()->GetBinUpEdge(pimean + sigmaBinPi);
         double Pidown = remass->GetXaxis()->GetBinLowEdge(pimean - sigmaBinPi);

         // Far end of SB Regions
         // Make SB1 Region width = sigmaBin
         double PiLower  = remass->GetXaxis()->GetBinLowEdge(pimean - 2 * sigmaBinPi - 1);

         // Since, Signal Region = (2 * sigmaBin) + 1 bin (central bin)
         // Make SB2 Region width = sigmaBin + 1 Bin
         // This way width of SB1 + SB2 = Signal Region
         // Remember there is a 1 Bin gap between Sideband and Signal Regions
         double PiHigher = remass->GetXaxis()->GetBinUpEdge(pimean + 2 * sigmaBinPi + 2);

         // Inner Edge of SB Regions
         double PiHigh = remass->GetXaxis()->GetBinUpEdge(pimean + sigmaBinPi + 1);  // So that signal and sideband regions do not overlap
         double PiLow  = remass->GetXaxis()->GetBinLowEdge(pimean - sigmaBinPi - 1); // So that signal and sideband regions do not overlap

         /*
         std::cout << "\n";
         std::cout << "PiLower = " << PiLower * 1000 << "\n";
         std::cout << "PiLow = " << PiLow * 1000 << "\n";
         std::cout << "Pidown = " << Pidown * 1000 << "\n";
         std::cout << "Piup = " << Piup * 1000 << "\n";
         std::cout << "PiHigh = " << PiHigh * 1000 << "\n";
         std::cout << "PiHigher = " << PiHigher * 1000 << "\n";
         std::cout << "\n";
         */

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

         /*
         std::cout << "\n";
         std::cout << "EtaLower = " << EtaLower * 1000 << "\n";
         std::cout << "EtaLow = " << EtaLow * 1000 << "\n";
         std::cout << "Etadown = " << Etadown * 1000 << "\n";
         std::cout << "Etaup = " << Etaup * 1000 << "\n";
         std::cout << "EtaHigh = " << EtaHigh * 1000 << "\n";
         std::cout << "EtaHigher = " << EtaHigher * 1000 << "\n";
         std::cout << "\n";
         */

			//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
			//                                                       //
			// Just used for plotting fits onto histograms           //
			// Separate calculation are performed for Eta & Etaprime //
			// in other TSelector analysis.                          //
			///////////////////////////////////////////////////////////

         // Signal Region
         double PercentPiSig = fsig1->Integral(Pidown, Piup) / binw;
         PercentPiSig /= TOTALPi;

         double PercentEtaSig = fsig2->Integral(Etadown, Etaup) / binw;
         PercentEtaSig /= TOTALEta;

         double PercentEtapSig = fsig4->Integral(XEtapDown, XEtapUp) / binw;
         PercentEtapSig /= TOTALEtap;


         // Lower Sideband Region
         double PercentPiSB1 = fsig1->Integral(PiLower, PiLow) / binw;
         PercentPiSB1 /= TOTALPi;

         double PercentEtaSB1 = fsig2->Integral(EtaLower, EtaLow) / binw;
         PercentEtaSB1 /= TOTALEta;

         double PercentEtapSB1 = fsig4->Integral(XEtapLower, XEtapLow) / binw;
         PercentEtapSB1 /= TOTALEtap;


         // Upper Sideband Region
         double PercentPiSB2 = fsig1->Integral(PiHigh, PiHigher) / binw;
         PercentPiSB2 /= TOTALPi;

         double PercentEtaSB2 = fsig2->Integral(EtaHigh, EtaHigher) / binw;
         PercentEtaSB2 /= TOTALEta;

         double PercentEtapSB2 = fsig4->Integral(XEtapHigh, XEtapHigher) / binw;
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
         newhfinal2_123->GetXaxis()->SetRangeUser(0.0, 0.25);
         double PMax = newhfinal2_123->GetMaximum();
         newhfinal2_123->GetXaxis()->SetRangeUser(0.04, 0.3);
         fback->Draw("same");
         ftot->Draw("same");

         TLine *zline123 = new TLine(Pidown, 0.0, Pidown, 0.5*PMax);
         zline123->SetLineColor(1);
         zline123->SetLineStyle(9);
         zline123->SetLineWidth(2);
         zline123->Draw("same");

         TLine *zline123a = new TLine(Piup, 0.0, Piup, 0.5*PMax);
         zline123a->SetLineColor(1);
         zline123a->SetLineStyle(9);
         zline123a->SetLineWidth(2);
         zline123a->Draw("same");

         TLine *zline123b = new TLine(PiLow, 0.0, PiLow, 0.25*PMax);
         zline123b->SetLineColor(6);
         zline123b->SetLineStyle(9);
         zline123b->SetLineWidth(2);
         zline123b->Draw("same");

         TLine *zline123c = new TLine(PiHigh, 0.0, PiHigh, 0.25*PMax);
         zline123c->SetLineColor(6);
         zline123c->SetLineStyle(9);
         zline123c->SetLineWidth(2);
         zline123c->Draw("same");

         TLine *zline123x = new TLine(PiLower, 0.0, PiLower, 0.25*PMax);
         zline123x->SetLineColor(6);
         zline123x->SetLineStyle(9);
         zline123x->SetLineWidth(2);
         zline123x->Draw("same");

         TLine *zline123z = new TLine(PiHigher, 0.0, PiHigher, 0.25*PMax);
         zline123z->SetLineColor(6);
         zline123z->SetLineStyle(9);
         zline123z->SetLineWidth(2);
         zline123z->Draw("same");

         TLine *zline123d = new TLine(mu_Pi, 0.0, mu_Pi, 1.025*PMax);
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
         xc123.Form("PiFit_E%d_t%d.png", Echannel, tchannel);
         c123->Print(xc123);

         // Eta Search
         double normEta = fsig2->Eval(mu_Eta);
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
         xc321.Form("EtaFit_E%d_t%d.png", Echannel, tchannel);
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

         TLine *zline555 = new TLine(XEtapDown, 0.0, XEtapDown, 1.0*normEtap);
         zline555->SetLineColor(1);
         zline555->SetLineStyle(9);
         zline555->SetLineWidth(2);
         zline555->Draw("same");

         TLine *zline555a = new TLine(XEtapUp, 0.0, XEtapUp, 1.0*normEtap);
         zline555a->SetLineColor(1);
         zline555a->SetLineStyle(9);
         zline555a->SetLineWidth(2);
         zline555a->Draw("same");

         TLine *zline555b = new TLine(XEtapLower, 0.0, XEtapLower, 0.75*normEtap);
         zline555b->SetLineColor(6);
         zline555b->SetLineStyle(9);
         zline555b->SetLineWidth(2);
         zline555b->Draw("same");

         TLine *zline555c = new TLine(XEtapHigher, 0.0, XEtapHigher, 0.75*normEtap);
         zline555c->SetLineColor(6);
         zline555c->SetLineStyle(9);
         zline555c->SetLineWidth(2);
         zline555c->Draw("same");

         TLine *zline555x = new TLine(XEtapLow, 0.0, XEtapLow, 0.75*normEtap);
         zline555x->SetLineColor(6);
         zline555x->SetLineStyle(9);
         zline555x->SetLineWidth(2);
         zline555x->Draw("same");

         TLine *zline555z = new TLine(XEtapHigh, 0.0, XEtapHigh, 0.75*normEtap);
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
         name555.Form("SB1 #rightarrow %.1f to %.1f MeV/c^{2}", XEtapLower*1000, XEtapLow*1000);

         TString name555sig;
         name555sig.Form("Signal #rightarrow %.1f to %.1f MeV/c^{2}", XEtapDown*1000, XEtapUp*1000);

         TString name555a;
         name555a.Form("SB2 #rightarrow %.1f to %.1f MeV/c^{2}", XEtapHigh*1000, XEtapHigher*1000);

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
         xc555.Form("EtapFit_E%d_t%d.png", Echannel, tchannel);
         c555->Print(xc555);

         /* Signal only in SIGNAL REGION */
         // Signal count in their respective regions
			double Sigcount_Pi  = fsig1->Integral(Pidown, Piup) / binw;
			double Sigcount_Eta = fsig2->Integral(Etadown, Etaup) / binw;


			/* Sideband 1 Region (left side). Background only. */
			/* B_1 */
			// Etap SB1 Region
			double ab1sb1        = fback->Integral(XEtapLower, XEtapLow) / binw;
			double ab3sb1        = fsig2->Integral(XEtapLower, XEtapLow) / binw;
			double ab4sb1        = fsig3->Integral(XEtapLower, XEtapLow) / binw;
			double sb1count_Etap  = ab1sb1 + ab3sb1 + ab4sb1;

			// Etap SB1 Region
         highX = XEtapLow;
         lowX  = XEtapLower;
         /////////////////////////////////
         /* Integral from lowX to highX */
         /////////////////////////////////
			double iEtap_b1 = sqrt(pow(p[25], 2.0)) * ( p[21] * (highX - lowX) +
                                                    (p[22] / 2.0) * (pow(highX, 2.0) - pow(lowX, 2.0)) +
                                                    (p[23] / 3.0) * (pow(highX, 3.0) - pow(lowX, 3.0)) +
                                                    (p[24] / 4.0) * (pow(highX, 4.0) - pow(lowX, 4.0)) ) +

                           (sqrt(PI / 2.0) * pow(p[6], 2.0) * p[8]  * TMath::Erf((p[7]  - lowX)  / (sqrt(2.0) * p[8]))) +
                           (sqrt(PI / 2.0) * pow(p[9], 2.0) * p[11] * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]))) -
                           (sqrt(PI / 2.0) * pow(p[6], 2.0) * p[8]  * TMath::Erf((p[7]  - highX) / (sqrt(2.0) * p[8]))) -
                           (sqrt(PI / 2.0) * pow(p[9], 2.0) * p[11] * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))) -

                           (sqrt(PI / 2.0) * pow(p[12], 2.0) * p[14] * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))) -
                           (sqrt(PI / 2.0) * pow(p[15], 2.0) * p[17] * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))) +
                           (sqrt(PI / 2.0) * pow(p[12], 2.0) * p[14] * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]))) +
                           (sqrt(PI / 2.0) * pow(p[15], 2.0) * p[17] * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17])));

         iEtap_b1 /= binw;

			double g21b1 = sqrt(pow(p[25], 2.0)) * (highX - lowX);

			double g22b1 = sqrt(pow(p[25], 2.0)) * (pow(highX, 2.0) - pow(lowX, 2.0)) / 2.0;

			double g23b1 = sqrt(pow(p[25], 2.0)) * (pow(highX, 3.0) - pow(lowX, 3.0)) / 3.0;

			double g24b1 = sqrt(pow(p[25], 2.0)) * (pow(highX, 4.0) - pow(lowX, 4.0)) / 4.0;

			double g25b1 = ( p[21] * (highX - lowX) +
                         (p[22] / 2.0) * (pow(highX, 2.0) - pow(lowX, 2.0)) +
                         (p[23] / 3.0) * (pow(highX, 3.0) - pow(lowX, 3.0)) +
                         (p[24] / 4.0) * (pow(highX, 4.0) - pow(lowX, 4.0)) );

         // Eta First gaussian parameters
         // Partial derivitive wrt p[6]
         double g6b1 =  -sqrt(2.0 * PI) * p[6] * p[8] * TMath::Erf((p[7] - highX) / (sqrt(2.0) * p[8]))
                       + sqrt(2.0 * PI) * p[6] * p[8] * TMath::Erf((p[7] - lowX)  / (sqrt(2.0) * p[8]));

         // Partial derivitive wrt p[7]
         double g7b1 = pow(p[6], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[7] - lowX), 2.0)  / (2.0 * pow(p[8], 2.0))))) - 
                                         TMath::Power(TMath::E(), (-(pow((p[7] - highX), 2.0) / (2.0 * pow(p[8], 2.0))))));

         // Partial derivitive wrt p[8]
         double g8b1 = pow(p[6], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[7] - lowX)  / (sqrt(2.0) * p[8]))) -
                                         (sqrt(PI / 2.0) * TMath::Erf((p[7] - highX) / (sqrt(2.0) * p[8]))) -
                                         ((p[7] - lowX)  / p[8]) * TMath::Power(TMath::E(), (-(pow((p[7] - lowX), 2.0))  / (2.0 * pow(p[8], 2.0)))) +
                                         ((p[7] - highX) / p[8]) * TMath::Power(TMath::E(), (-(pow((p[7] - highX), 2.0)) / (2.0 * pow(p[8], 2.0)))));

         // Eta Second gaussian parameters
         // Partial derivitive wrt p[9]
         double g9b1 =  -sqrt(2.0 * PI) * p[9] * p[11] * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))
                       + sqrt(2.0 * PI) * p[9] * p[11] * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]));

         // Partial derivitive wrt p[10]
         double g10b1 = pow(p[9], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[10] - lowX), 2.0)  / (2.0 * pow(p[11], 2.0))))) - 
                                          TMath::Power(TMath::E(), (-(pow((p[10] - highX), 2.0) / (2.0 * pow(p[11], 2.0))))));

         // Partial derivitive wrt p[11]
         double g11b1 = pow(p[9], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]))) -
                                          (sqrt(PI / 2.0) * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))) -
                                          ((p[10] - lowX)  / p[11]) * TMath::Power(TMath::E(), (-(pow((p[10] - lowX), 2.0))  / (2.0 * pow(p[11], 2.0)))) +
                                          ((p[10] - highX) / p[11]) * TMath::Power(TMath::E(), (-(pow((p[10] - highX), 2.0)) / (2.0 * pow(p[11], 2.0)))));

         // Omega First gaussian parameters
         // Partial derivitive wrt p[12]
         double g12b1 =  -sqrt(2.0 * PI) * p[12] * p[14] * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))
                        + sqrt(2.0 * PI) * p[12] * p[14] * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]));

         // Partial derivitive wrt p[13]
         double g13b1 = pow(p[12], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[13] - lowX), 2.0)  / (2.0 * pow(p[14], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[13] - highX), 2.0) / (2.0 * pow(p[14], 2.0))))));

         // Partial derivitive wrt p[14]
         double g14b1 = pow(p[12], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))) -
                                           ((p[13] - lowX)  / p[14]) * TMath::Power(TMath::E(), (-(pow((p[13] - lowX), 2.0))  / (2.0 * pow(p[14], 2.0)))) +
                                           ((p[13] - highX) / p[14]) * TMath::Power(TMath::E(), (-(pow((p[13] - highX), 2.0)) / (2.0 * pow(p[14], 2.0)))));

         // Omega Second gaussian parameters
         // Partial derivitive wrt p[15]
         double g15b1 =  -sqrt(2.0 * PI) * p[15] * p[17] * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))
                        + sqrt(2.0 * PI) * p[15] * p[17] * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17]));

         // Partial derivitive wrt p[16]
         double g16b1 = pow(p[15], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[16] - lowX), 2.0)  / (2.0 * pow(p[17], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[16] - highX), 2.0) / (2.0 * pow(p[17], 2.0))))));

         // Partial derivitive wrt p[17]
         double g17b1 = pow(p[15], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))) -
                                           ((p[16] - lowX)  / p[17]) * TMath::Power(TMath::E(), (-(pow((p[16] - lowX), 2.0))  / (2.0 * pow(p[17], 2.0)))) +
                                           ((p[16] - highX) / p[17]) * TMath::Power(TMath::E(), (-(pow((p[16] - highX), 2.0)) / (2.0 * pow(p[17], 2.0)))));


			double EtapGaus_b1 = sqrt(std::fabs(g21b1 * g21b1 * covMatrix(21, 21) + g22b1 * g22b1 * covMatrix(22, 22) + g23b1 * g23b1 * covMatrix(23, 23) +
                                             g24b1 * g24b1 * covMatrix(24, 24) + g25b1 * g25b1 * covMatrix(25, 25) +

                                             2.0 * g21b1 * g22b1 * covMatrix(21, 22) + 2.0 * g21b1 * g23b1 * covMatrix(21, 23) +
                                             2.0 * g21b1 * g24b1 * covMatrix(21, 24) + 2.0 * g21b1 * g25b1 * covMatrix(21, 25) +

                                             2.0 * g22b1 * g23b1 * covMatrix(22, 23) + 2.0 * g22b1 * g24b1 * covMatrix(22, 24) +
                                             2.0 * g22b1 * g25b1 * covMatrix(22, 25) +

                                             2.0 * g23b1 * g24b1 * covMatrix(23, 24) + 2.0 * g23b1 * g25b1 * covMatrix(23, 25) +

                                             2.0 * g24b1 * g25b1 * covMatrix(24, 25)) +


                                   std::fabs(g6b1 * g6b1 * covMatrix(6, 6) + g7b1 * g7b1 * covMatrix(7, 7) + g8b1 * g8b1 * covMatrix(8, 8) +
                                             g9b1 * g9b1 * covMatrix(9, 9) + g10b1 * g10b1 * covMatrix(10, 10) + g11b1 * g11b1 * covMatrix(11, 11) +

                                             2.0 * g6b1 * g7b1 * covMatrix(6, 7) + 2.0 * g6b1 * g8b1 * covMatrix(6, 8) +
                                             2.0 * g6b1 * g9b1 * covMatrix(6, 9) + 2.0 * g6b1 * g10b1 * covMatrix(6, 10) + 2.0 * g6b1 * g11b1 * covMatrix(6, 11) +

                                             2.0 * g7b1 * g8b1 * covMatrix(7, 8) + 2.0 * g7b1 * g9b1 * covMatrix(7, 9) +
                                             2.0 * g7b1 * g10b1 * covMatrix(7, 10) + 2.0 * g7b1 * g11b1 * covMatrix(7, 11) +

                                             2.0 * g8b1 * g9b1 * covMatrix(8, 9) + 2.0 * g8b1 * g10b1 * covMatrix(8, 10) +
                                             2.0 * g8b1 * g11b1 * covMatrix(8, 11) +

                                             2.0 * g9b1 * g10b1 * covMatrix(9, 10) + 2.0 * g9b1 * g11b1 * covMatrix(9, 11) +

                                             2.0 * g10b1 * g11b1 * covMatrix(10, 11)) +


                                   std::fabs(g12b1 * g12b1 * covMatrix(12, 12) + g13b1 * g13b1 * covMatrix(13, 13) + g14b1 * g14b1 * covMatrix(14, 14) +
                                             g15b1 * g15b1 * covMatrix(15, 15) + g16b1 * g16b1 * covMatrix(16, 16) + g17b1 * g17b1 * covMatrix(17, 17) +

                                             2.0 * g12b1 * g13b1 * covMatrix(12, 13) + 2.0 * g12b1 * g14b1 * covMatrix(12, 14) +
                                             2.0 * g12b1 * g15b1 * covMatrix(12, 15) + 2.0 * g12b1 * g16b1 * covMatrix(12, 16) + 2.0 * g12b1 * g17b1 * covMatrix(12, 17) +

                                             2.0 * g13b1 * g14b1 * covMatrix(13, 14) + 2.0 * g13b1 * g15b1 * covMatrix(13, 15) +
                                             2.0 * g13b1 * g16b1 * covMatrix(13, 16) + 2.0 * g13b1 * g17b1 * covMatrix(13, 17) +

                                             2.0 * g14b1 * g15b1 * covMatrix(14, 15) + 2.0 * g14b1 * g16b1 * covMatrix(14, 16) +
                                             2.0 * g14b1 * g17b1 * covMatrix(14, 17) +

                                             2.0 * g15b1 * g16b1 * covMatrix(15, 16) + 2.0 * g15b1 * g17b1 * covMatrix(15, 17) +

                                             2.0 * g16b1 * g17b1 * covMatrix(16, 17)));
         EtapGaus_b1 /= binw;

         std::cout << " Check B1 Region" << "\n";
         std::cout << "Integral = " << sb1count_Etap << "\n";
         std::cout << "Manual = " << iEtap_b1 << "\n";
         std::cout << "Error = " << EtapGaus_b1 << "\n";
         std::cout << "\n";


			/* Sideband 2 Region (right side). Background only. */
			/* B_2 */
			// Etap SB2 Region
			double aa1sb2 = fback->Integral(XEtapHigh, XEtapHigher) / binw;
			double aa3sb2 = fsig2->Integral(XEtapHigh, XEtapHigher) / binw;
			double aa4sb2 = fsig3->Integral(XEtapHigh, XEtapHigher) / binw;
			double sb2count_Etap = aa1sb2 + aa3sb2 + aa4sb2;

			// Etap SB2 Region
         highX = XEtapHigher;
         lowX  = XEtapHigh;
         //////////////////////////////////////////////////
         /* Integral of SB2 Backgound from lowX to highX */
         //////////////////////////////////////////////////
			double iEtap_b2 = sqrt(pow(p[25], 2.0)) * ( p[21] * (highX - lowX) +
                                                    (p[22] / 2.0) * (pow(highX, 2.0) - pow(lowX, 2.0)) +
                                                    (p[23] / 3.0) * (pow(highX, 3.0) - pow(lowX, 3.0)) +
                                                    (p[24] / 4.0) * (pow(highX, 4.0) - pow(lowX, 4.0)) ) +

                           (sqrt(PI / 2.0) * pow(p[6], 2.0) * p[8]  * TMath::Erf((p[7]  - lowX)  / (sqrt(2.0) * p[8]))) +
                           (sqrt(PI / 2.0) * pow(p[9], 2.0) * p[11] * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]))) -
                           (sqrt(PI / 2.0) * pow(p[6], 2.0) * p[8]  * TMath::Erf((p[7]  - highX) / (sqrt(2.0) * p[8]))) -
                           (sqrt(PI / 2.0) * pow(p[9], 2.0) * p[11] * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))) -

                           (sqrt(PI / 2.0) * pow(p[12], 2.0) * p[14] * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))) -
                           (sqrt(PI / 2.0) * pow(p[15], 2.0) * p[17] * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))) +
                           (sqrt(PI / 2.0) * pow(p[12], 2.0) * p[14] * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]))) +
                           (sqrt(PI / 2.0) * pow(p[15], 2.0) * p[17] * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17])));
         iEtap_b2 /= binw;

			double g21b2 = sqrt(pow(p[25], 2.0)) * (highX - lowX);

			double g22b2 = sqrt(pow(p[25], 2.0)) * (pow(highX, 2.0) - pow(lowX, 2.0)) / 2.0;

			double g23b2 = sqrt(pow(p[25], 2.0)) * (pow(highX, 3.0) - pow(lowX, 3.0)) / 3.0;

			double g24b2 = sqrt(pow(p[25], 2.0)) * (pow(highX, 4.0) - pow(lowX, 4.0)) / 4.0;

			double g25b2 = ( p[21] * (highX - lowX) +
                         (p[22] / 2.0) * (pow(highX, 2.0) - pow(lowX, 2.0)) +
                         (p[23] / 3.0) * (pow(highX, 3.0) - pow(lowX, 3.0)) +
                         (p[24] / 4.0) * (pow(highX, 4.0) - pow(lowX, 4.0)) );

         // Eta First gaussian parameters
         // Partial derivitive wrt p[6]
         double g6b2 =  -sqrt(2.0 * PI) * p[6] * p[8] * TMath::Erf((p[7] - highX) / (sqrt(2.0) * p[8]))
                       + sqrt(2.0 * PI) * p[6] * p[8] * TMath::Erf((p[7] - lowX)  / (sqrt(2.0) * p[8]));

         // Partial derivitive wrt p[7]
         double g7b2 = pow(p[6], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[7] - lowX), 2.0)  / (2.0 * pow(p[8], 2.0))))) - 
                                         TMath::Power(TMath::E(), (-(pow((p[7] - highX), 2.0) / (2.0 * pow(p[8], 2.0))))));

         // Partial derivitive wrt p[8]
         double g8b2 = pow(p[6], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[7] - lowX)  / (sqrt(2.0) * p[8]))) -
                                         (sqrt(PI / 2.0) * TMath::Erf((p[7] - highX) / (sqrt(2.0) * p[8]))) -
                                         ((p[7] - lowX)  / p[8]) * TMath::Power(TMath::E(), (-(pow((p[7] - lowX), 2.0))  / (2.0 * pow(p[8], 2.0)))) +
                                         ((p[7] - highX) / p[8]) * TMath::Power(TMath::E(), (-(pow((p[7] - highX), 2.0)) / (2.0 * pow(p[8], 2.0)))));

         // Eta Second gaussian parameters
         // Partial derivitive wrt p[9]
         double g9b2 =  -sqrt(2.0 * PI) * p[9] * p[11] * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))
                       + sqrt(2.0 * PI) * p[9] * p[11] * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]));

         // Partial derivitive wrt p[10]
         double g10b2 = pow(p[9], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[10] - lowX), 2.0)  / (2.0 * pow(p[11], 2.0))))) - 
                                          TMath::Power(TMath::E(), (-(pow((p[10] - highX), 2.0) / (2.0 * pow(p[11], 2.0))))));

         // Partial derivitive wrt p[11]
         double g11b2 = pow(p[9], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]))) -
                                          (sqrt(PI / 2.0) * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))) -
                                          ((p[10] - lowX)  / p[11]) * TMath::Power(TMath::E(), (-(pow((p[10] - lowX), 2.0))  / (2.0 * pow(p[11], 2.0)))) +
                                          ((p[10] - highX) / p[11]) * TMath::Power(TMath::E(), (-(pow((p[10] - highX), 2.0)) / (2.0 * pow(p[11], 2.0)))));

         // Omega First gaussian parameters
         // Partial derivitive wrt p[12]
         double g12b2 =  -sqrt(2.0 * PI) * p[12] * p[14] * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))
                        + sqrt(2.0 * PI) * p[12] * p[14] * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]));

         // Partial derivitive wrt p[13]
         double g13b2 = pow(p[12], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[13] - lowX), 2.0)  / (2.0 * pow(p[14], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[13] - highX), 2.0) / (2.0 * pow(p[14], 2.0))))));

         // Partial derivitive wrt p[14]
         double g14b2 = pow(p[12], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))) -
                                           ((p[13] - lowX)  / p[14]) * TMath::Power(TMath::E(), (-(pow((p[13] - lowX), 2.0))  / (2.0 * pow(p[14], 2.0)))) +
                                           ((p[13] - highX) / p[14]) * TMath::Power(TMath::E(), (-(pow((p[13] - highX), 2.0)) / (2.0 * pow(p[14], 2.0)))));

         // Omega Second gaussian parameters
         // Partial derivitive wrt p[15]
         double g15b2 =  -sqrt(2.0 * PI) * p[15] * p[17] * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))
                        + sqrt(2.0 * PI) * p[15] * p[17] * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17]));

         // Partial derivitive wrt p[16]
         double g16b2 = pow(p[15], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[16] - lowX), 2.0)  / (2.0 * pow(p[17], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[16] - highX), 2.0) / (2.0 * pow(p[17], 2.0))))));

         // Partial derivitive wrt p[17]
         double g17b2 = pow(p[15], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))) -
                                           ((p[16] - lowX)  / p[17]) * TMath::Power(TMath::E(), (-(pow((p[16] - lowX), 2.0))  / (2.0 * pow(p[17], 2.0)))) +
                                           ((p[16] - highX) / p[17]) * TMath::Power(TMath::E(), (-(pow((p[16] - highX), 2.0)) / (2.0 * pow(p[17], 2.0)))));


			double EtapGaus_b2 = sqrt(std::fabs(g21b2 * g21b2 * covMatrix(21, 21) + g22b2 * g22b2 * covMatrix(22, 22) + g23b2 * g23b2 * covMatrix(23, 23) +
                                             g24b2 * g24b2 * covMatrix(24, 24) + g25b2 * g25b2 * covMatrix(25, 25) +

                                             2.0 * g21b2 * g22b2 * covMatrix(21, 22) + 2.0 * g21b2 * g23b2 * covMatrix(21, 23) +
                                             2.0 * g21b2 * g24b2 * covMatrix(21, 24) + 2.0 * g21b2 * g25b2 * covMatrix(21, 25) +

                                             2.0 * g22b2 * g23b2 * covMatrix(22, 23) + 2.0 * g22b2 * g24b2 * covMatrix(22, 24) +
                                             2.0 * g22b2 * g25b2 * covMatrix(22, 25) +

                                             2.0 * g23b2 * g24b2 * covMatrix(23, 24) + 2.0 * g23b2 * g25b2 * covMatrix(23, 25) +

                                             2.0 * g24b2 * g25b2 * covMatrix(24, 25)) +


                                   std::fabs(g6b2 * g6b2 * covMatrix(6, 6) + g7b2 * g7b2 * covMatrix(7, 7) + g8b2 * g8b2 * covMatrix(8, 8) +
                                             g9b2 * g9b2 * covMatrix(9, 9) + g10b2 * g10b2 * covMatrix(10, 10) + g11b2 * g11b2 * covMatrix(11, 11) +

                                             2.0 * g6b2 * g7b2 * covMatrix(6, 7) + 2.0 * g6b2 * g8b2 * covMatrix(6, 8) +
                                             2.0 * g6b2 * g9b2 * covMatrix(6, 9) + 2.0 * g6b2 * g10b2 * covMatrix(6, 10) + 2.0 * g6b2 * g11b2 * covMatrix(6, 11) +

                                             2.0 * g7b2 * g8b2 * covMatrix(7, 8) + 2.0 * g7b2 * g9b2 * covMatrix(7, 9) +
                                             2.0 * g7b2 * g10b2 * covMatrix(7, 10) + 2.0 * g7b2 * g11b2 * covMatrix(7, 11) +

                                             2.0 * g8b2 * g9b2 * covMatrix(8, 9) + 2.0 * g8b2 * g10b2 * covMatrix(8, 10) +
                                             2.0 * g8b2 * g11b2 * covMatrix(8, 11) +

                                             2.0 * g9b2 * g10b2 * covMatrix(9, 10) + 2.0 * g9b2 * g11b2 * covMatrix(9, 11) +

                                             2.0 * g10b2 * g11b2 * covMatrix(10, 11)) +


                                   std::fabs(g12b2 * g12b2 * covMatrix(12, 12) + g13b2 * g13b2 * covMatrix(13, 13) + g14b2 * g14b2 * covMatrix(14, 14) +
                                             g15b2 * g15b2 * covMatrix(15, 15) + g16b2 * g16b2 * covMatrix(16, 16) + g17b2 * g17b2 * covMatrix(17, 17) +

                                             2.0 * g12b2 * g13b2 * covMatrix(12, 13) + 2.0 * g12b2 * g14b2 * covMatrix(12, 14) +
                                             2.0 * g12b2 * g15b2 * covMatrix(12, 15) + 2.0 * g12b2 * g16b2 * covMatrix(12, 16) + 2.0 * g12b2 * g17b2 * covMatrix(12, 17) +

                                             2.0 * g13b2 * g14b2 * covMatrix(13, 14) + 2.0 * g13b2 * g15b2 * covMatrix(13, 15) +
                                             2.0 * g13b2 * g16b2 * covMatrix(13, 16) + 2.0 * g13b2 * g17b2 * covMatrix(13, 17) +

                                             2.0 * g14b2 * g15b2 * covMatrix(14, 15) + 2.0 * g14b2 * g16b2 * covMatrix(14, 16) +
                                             2.0 * g14b2 * g17b2 * covMatrix(14, 17) +

                                             2.0 * g15b2 * g16b2 * covMatrix(15, 16) + 2.0 * g15b2 * g17b2 * covMatrix(15, 17) +

                                             2.0 * g16b2 * g17b2 * covMatrix(16, 17)));
         EtapGaus_b2 /= binw;

         std::cout << " Check B2 Region" << "\n";
         std::cout << "Integral = " << sb2count_Etap << "\n";
         std::cout << "Manual = " << iEtap_b2 << "\n";
         std::cout << "Error = " << EtapGaus_b2 << "\n";
         std::cout << "\n";

         // vvvv Troubleshooting vvvv

         //Output file to check alternate method for error calculations
			ofstream FitErrorCalc;
         TString EName;
         EName.Form("E%d_t%d", Echannel, tchannel);
         FitErrorCalc.open("Etap_fit_status_55MeV.txt", std::ios_base::app);
         if (FitErrorCalc.is_open()) {
            FitErrorCalc << "___ " << EName << " ___" << "\n";
            FitErrorCalc << "B2               = " << iEtap_b2 << "\n";
            FitErrorCalc << "B2 CovM Error    = " << EtapGaus_b2 << "\n";
            FitErrorCalc << "Percent Error    = " << EtapGaus_b2 / sb2count_Etap * 100.0 << "% \n";
         }
         else {
            std::cout << "Unable to open alternate error calculation check file." << "\n";
         }
         FitErrorCalc.close();
         // ^^^^ Troubleshooting ^^^^

			/* Background only in SIGNAL REGION */
			/* B_0 */
			// Eta prime Signal Region
			double BgEta_Etap   = fsig2->Integral(XEtapDown, XEtapUp) / binw;
			double Bgleak_Etap  = fsig3->Integral(XEtapDown, XEtapUp) / binw;
			double BgBack_Etap  = fback->Integral(XEtapDown, XEtapUp) / binw;
			double Bgcount_Etap = BgEta_Etap + Bgleak_Etap + BgBack_Etap;

			// Etap Signal Region
         highX = XEtapUp;
         lowX  = XEtapDown;
         /////////////////////////////////
         /* Integral from lowX to highX */
         /////////////////////////////////
			double iEtap_b0 = sqrt(pow(p[25], 2.0)) * ( p[21] * (highX - lowX) +
                                                    (p[22] / 2.0) * (pow(highX, 2.0) - pow(lowX, 2.0)) +
                                                    (p[23] / 3.0) * (pow(highX, 3.0) - pow(lowX, 3.0)) +
                                                    (p[24] / 4.0) * (pow(highX, 4.0) - pow(lowX, 4.0)) ) +

                           (sqrt(PI / 2.0) * pow(p[0], 2.0) * p[2]  * TMath::Erf((p[1] - lowX)  / (sqrt(2.0) * p[2]))) +
                           (sqrt(PI / 2.0) * pow(p[3], 2.0) * p[5]  * TMath::Erf((p[4] - lowX)  / (sqrt(2.0) * p[5]))) -
                           (sqrt(PI / 2.0) * pow(p[0], 2.0) * p[2]  * TMath::Erf((p[1] - highX) / (sqrt(2.0) * p[2]))) -
                           (sqrt(PI / 2.0) * pow(p[3], 2.0) * p[5]  * TMath::Erf((p[4] - highX) / (sqrt(2.0) * p[5]))) +

                           (sqrt(PI / 2.0) * pow(p[6], 2.0) * p[8]  * TMath::Erf((p[7]  - lowX)  / (sqrt(2.0) * p[8]))) +
                           (sqrt(PI / 2.0) * pow(p[9], 2.0) * p[11] * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]))) -
                           (sqrt(PI / 2.0) * pow(p[6], 2.0) * p[8]  * TMath::Erf((p[7]  - highX) / (sqrt(2.0) * p[8]))) -
                           (sqrt(PI / 2.0) * pow(p[9], 2.0) * p[11] * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))) -

                           (sqrt(PI / 2.0) * pow(p[12], 2.0) * p[14] * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))) -
                           (sqrt(PI / 2.0) * pow(p[15], 2.0) * p[17] * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))) +
                           (sqrt(PI / 2.0) * pow(p[12], 2.0) * p[14] * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]))) +
                           (sqrt(PI / 2.0) * pow(p[15], 2.0) * p[17] * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17])));
         iEtap_b0 /= binw;

			double g21b0 = ( p[21] * (highX - lowX) +
                         (p[22] / 2.0) * (pow(highX, 2.0) - pow(lowX, 2.0)) +
                         (p[23] / 3.0) * (pow(highX, 3.0) - pow(lowX, 3.0)) +
                         (p[24] / 4.0) * (pow(highX, 4.0) - pow(lowX, 4.0)) );

			double g22b0 = sqrt(pow(p[25], 2.0)) * (highX - lowX);

			double g23b0 = sqrt(pow(p[25], 2.0)) * (pow(highX, 2.0) - pow(lowX, 2.0)) / 2.0;

			double g24b0 = sqrt(pow(p[25], 2.0)) * (pow(highX, 3.0) - pow(lowX, 3.0)) / 3.0;

			double g25b0 = sqrt(pow(p[25], 2.0)) * (pow(highX, 4.0) - pow(lowX, 4.0)) / 4.0;

         // Pi0 First gaussian parameters
         // Partial derivitive wrt p[0]
         double g0b0 =  -sqrt(2.0 * PI) * p[0] * p[2] * TMath::Erf((p[1] - highX) / (sqrt(2.0) * p[2]))
                       + sqrt(2.0 * PI) * p[0] * p[2] * TMath::Erf((p[1] - lowX)  / (sqrt(2.0) * p[2]));

         // Partial derivitive wrt p[1]
         double g1b0 = pow(p[0], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[1] - lowX), 2.0)  / (2.0 * pow(p[2], 2.0))))) - 
                                         TMath::Power(TMath::E(), (-(pow((p[1] - highX), 2.0) / (2.0 * pow(p[2], 2.0))))));

         // Partial derivitive wrt p[2]
         double g2b0 = pow(p[0], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[1] - lowX)  / (sqrt(2.0) * p[2]))) -
                                         (sqrt(PI / 2.0) * TMath::Erf((p[1] - highX) / (sqrt(2.0) * p[2]))) -
                                         ((p[1] - lowX)  / p[2]) * TMath::Power(TMath::E(), (-(pow((p[1] - lowX), 2.0))  / (2.0 * pow(p[2], 2.0)))) +
                                         ((p[1] - highX) / p[2]) * TMath::Power(TMath::E(), (-(pow((p[1] - highX), 2.0)) / (2.0 * pow(p[2], 2.0)))));

         // Pi0 Second gaussian parameters
         // Partial derivitive wrt p[3]
         double g3b0 =  -sqrt(2.0 * PI) * p[3] * p[5] * TMath::Erf((p[4] - highX) / (sqrt(2.0) * p[5]))
                       + sqrt(2.0 * PI) * p[3] * p[5] * TMath::Erf((p[4] - lowX)  / (sqrt(2.0) * p[5]));

         // Partial derivitive wrt p[4]
         double g4b0 = pow(p[3], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[4] - lowX), 2.0)  / (2.0 * pow(p[5], 2.0))))) - 
                                         TMath::Power(TMath::E(), (-(pow((p[4] - highX), 2.0) / (2.0 * pow(p[5], 2.0))))));

         // Partial derivitive wrt p[5]
         double g5b0 = pow(p[3], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[4] - lowX)  / (sqrt(2.0) * p[5]))) -
                                        (sqrt(PI / 2.0) * TMath::Erf((p[4] - highX) / (sqrt(2.0) * p[5]))) -
                                        ((p[4] - lowX)  / p[5]) * TMath::Power(TMath::E(), (-(pow((p[4] - lowX), 2.0))  / (2.0 * pow(p[5], 2.0)))) +
                                        ((p[4] - highX) / p[5]) * TMath::Power(TMath::E(), (-(pow((p[4] - highX), 2.0)) / (2.0 * pow(p[5], 2.0)))));

         // Eta First gaussian parameters
         // Partial derivitive wrt p[6]
         double g6b0 =  -sqrt(2.0 * PI) * p[6] * p[8] * TMath::Erf((p[7] - highX) / (sqrt(2.0) * p[8]))
                       + sqrt(2.0 * PI) * p[6] * p[8] * TMath::Erf((p[7] - lowX)  / (sqrt(2.0) * p[8]));

         // Partial derivitive wrt p[7]
         double g7b0 = pow(p[6], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[7] - lowX), 2.0)  / (2.0 * pow(p[8], 2.0))))) - 
                                         TMath::Power(TMath::E(), (-(pow((p[7] - highX), 2.0) / (2.0 * pow(p[8], 2.0))))));

         // Partial derivitive wrt p[8]
         double g8b0 = pow(p[6], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[7] - lowX)  / (sqrt(2.0) * p[8]))) -
                                         (sqrt(PI / 2.0) * TMath::Erf((p[7] - highX) / (sqrt(2.0) * p[8]))) -
                                         ((p[7] - lowX)  / p[8]) * TMath::Power(TMath::E(), (-(pow((p[7] - lowX), 2.0))  / (2.0 * pow(p[8], 2.0)))) +
                                         ((p[7] - highX) / p[8]) * TMath::Power(TMath::E(), (-(pow((p[7] - highX), 2.0)) / (2.0 * pow(p[8], 2.0)))));

         // Eta Second gaussian parameters
         // Partial derivitive wrt p[9]
         double g9b0 =  -sqrt(2.0 * PI) * p[9] * p[11] * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))
                       + sqrt(2.0 * PI) * p[9] * p[11] * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]));

         // Partial derivitive wrt p[10]
         double g10b0 = pow(p[9], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[10] - lowX), 2.0)  / (2.0 * pow(p[11], 2.0))))) - 
                                          TMath::Power(TMath::E(), (-(pow((p[10] - highX), 2.0) / (2.0 * pow(p[11], 2.0))))));

         // Partial derivitive wrt p[11]
         double g11b0 = pow(p[9], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[10] - lowX)  / (sqrt(2.0) * p[11]))) -
                                          (sqrt(PI / 2.0) * TMath::Erf((p[10] - highX) / (sqrt(2.0) * p[11]))) -
                                          ((p[10] - lowX)  / p[11]) * TMath::Power(TMath::E(), (-(pow((p[10] - lowX), 2.0))  / (2.0 * pow(p[11], 2.0)))) +
                                          ((p[10] - highX) / p[11]) * TMath::Power(TMath::E(), (-(pow((p[10] - highX), 2.0)) / (2.0 * pow(p[11], 2.0)))));

         // Omega First gaussian parameters
         // Partial derivitive wrt p[12]
         double g12b0 =  -sqrt(2.0 * PI) * p[12] * p[14] * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))
                        + sqrt(2.0 * PI) * p[12] * p[14] * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]));

         // Partial derivitive wrt p[13]
         double g13b0 = pow(p[12], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[13] - lowX), 2.0)  / (2.0 * pow(p[14], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[13] - highX), 2.0) / (2.0 * pow(p[14], 2.0))))));

         // Partial derivitive wrt p[14]
         double g14b0 = pow(p[12], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[13] - lowX)  / (sqrt(2.0) * p[14]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[13] - highX) / (sqrt(2.0) * p[14]))) -
                                           ((p[13] - lowX)  / p[14]) * TMath::Power(TMath::E(), (-(pow((p[13] - lowX), 2.0))  / (2.0 * pow(p[14], 2.0)))) +
                                           ((p[13] - highX) / p[14]) * TMath::Power(TMath::E(), (-(pow((p[13] - highX), 2.0)) / (2.0 * pow(p[14], 2.0)))));

         // Omega Second gaussian parameters
         // Partial derivitive wrt p[15]
         double g15b0 =  -sqrt(2.0 * PI) * p[15] * p[17] * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))
                        + sqrt(2.0 * PI) * p[15] * p[17] * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17]));

         // Partial derivitive wrt p[16]
         double g16b0 = pow(p[15], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[16] - lowX), 2.0)  / (2.0 * pow(p[17], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[16] - highX), 2.0) / (2.0 * pow(p[17], 2.0))))));

         // Partial derivitive wrt p[17]
         double g17b0 = pow(p[15], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[16] - lowX)  / (sqrt(2.0) * p[17]))) -
                                           (sqrt(PI / 2.0) * TMath::Erf((p[16] - highX) / (sqrt(2.0) * p[17]))) -
                                           ((p[16] - lowX)  / p[17]) * TMath::Power(TMath::E(), (-(pow((p[16] - lowX), 2.0))  / (2.0 * pow(p[17], 2.0)))) +
                                           ((p[16] - highX) / p[17]) * TMath::Power(TMath::E(), (-(pow((p[16] - highX), 2.0)) / (2.0 * pow(p[17], 2.0)))));


			double EtapGaus_b0 = sqrt(std::fabs(g21b0 * g21b0 * covMatrix(21, 21) + g22b0 * g22b0 * covMatrix(22, 22) + g23b0 * g23b0 * covMatrix(23, 23) +
                                             g24b0 * g24b0 * covMatrix(24, 24) + g25b0 * g25b0 * covMatrix(25, 25) +

                                             2.0 * g21b0 * g22b0 * covMatrix(21, 22) + 2.0 * g21b0 * g23b0 * covMatrix(21, 23) +
                                             2.0 * g21b0 * g24b0 * covMatrix(21, 24) + 2.0 * g21b0 * g25b0 * covMatrix(21, 25) +

                                             2.0 * g22b0 * g23b0 * covMatrix(22, 23) + 2.0 * g22b0 * g24b0 * covMatrix(22, 24) +
                                             2.0 * g22b0 * g25b0 * covMatrix(22, 25) +

                                             2.0 * g23b0 * g24b0 * covMatrix(23, 24) + 2.0 * g23b0 * g25b0 * covMatrix(23, 25) +

                                             2.0 * g24b0 * g25b0 * covMatrix(24, 25)) +

                                   std::fabs(g0b0 * g0b0 * covMatrix(0, 0) + g1b0 * g1b0 * covMatrix(1, 1) + g2b0 * g2b0 * covMatrix(2, 2) +
                                             g3b0 * g3b0 * covMatrix(3, 3) + g4b0 * g4b0 * covMatrix(4, 4) + g5b0 * g5b0 * covMatrix(5, 5) +

                                             2.0 * g0b0 * g1b0 * covMatrix(0, 1) + 2.0 * g0b0 * g2b0 * covMatrix(0, 2) +
                                             2.0 * g0b0 * g3b0 * covMatrix(0, 3) + 2.0 * g0b0 * g4b0 * covMatrix(0, 4) + 2.0 * g0b0 * g5b0 * covMatrix(0, 5) +

                                             2.0 * g1b0 * g2b0 * covMatrix(1, 2) + 2.0 * g1b0 * g3b0 * covMatrix(1, 3) +
                                             2.0 * g1b0 * g4b0 * covMatrix(1, 4) + 2.0 * g1b0 * g5b0 * covMatrix(1, 5) +

                                             2.0 * g2b0 * g3b0 * covMatrix(2, 3) + 2.0 * g2b0 * g4b0 * covMatrix(2, 4) +
                                             2.0 * g2b0 * g5b0 * covMatrix(2, 5) +

                                             2.0 * g3b0 * g4b0 * covMatrix(3, 4) + 2.0 * g3b0 * g5b0 * covMatrix(3, 5) +

                                             2.0 * g4b0 * g5b0 * covMatrix(4, 5)) +

                                   std::fabs(g6b0 * g6b0 * covMatrix(6, 6) + g7b0 * g7b0 * covMatrix(7, 7) + g8b0 * g8b0 * covMatrix(8, 8) +
                                             g9b0 * g9b0 * covMatrix(9, 9) + g10b0 * g10b0 * covMatrix(10, 10) + g11b0 * g11b0 * covMatrix(11, 11) +

                                             2.0 * g6b0 * g7b0 * covMatrix(6, 7) + 2.0 * g6b0 * g8b0 * covMatrix(6, 8) +
                                             2.0 * g6b0 * g9b0 * covMatrix(6, 9) + 2.0 * g6b0 * g10b0 * covMatrix(6, 10) + 2.0 * g6b0 * g11b0 * covMatrix(6, 11) +

                                             2.0 * g7b0 * g8b0 * covMatrix(7, 8) + 2.0 * g7b0 * g9b0 * covMatrix(7, 9) +
                                             2.0 * g7b0 * g10b0 * covMatrix(7, 10) + 2.0 * g7b0 * g11b0 * covMatrix(7, 11) +

                                             2.0 * g8b0 * g9b0 * covMatrix(8, 9) + 2.0 * g8b0 * g10b0 * covMatrix(8, 10) +
                                             2.0 * g8b0 * g11b0 * covMatrix(8, 11) +

                                             2.0 * g9b0 * g10b0 * covMatrix(9, 10) + 2.0 * g9b0 * g11b0 * covMatrix(9, 11) +

                                             2.0 * g10b0 * g11b0 * covMatrix(10, 11)) +

                                   std::fabs(g12b0 * g12b0 * covMatrix(12, 12) + g13b0 * g13b0 * covMatrix(13, 13) + g14b0 * g14b0 * covMatrix(14, 14) +
                                             g15b0 * g15b0 * covMatrix(15, 15) + g16b0 * g16b0 * covMatrix(16, 16) + g17b0 * g17b0 * covMatrix(17, 17) +

                                             2.0 * g12b0 * g13b0 * covMatrix(12, 13) + 2.0 * g12b0 * g14b0 * covMatrix(12, 14) +
                                             2.0 * g12b0 * g15b0 * covMatrix(12, 15) + 2.0 * g12b0 * g16b0 * covMatrix(12, 16) + 2.0 * g12b0 * g17b0 * covMatrix(12, 17) +

                                             2.0 * g13b0 * g14b0 * covMatrix(13, 14) + 2.0 * g13b0 * g15b0 * covMatrix(13, 15) +
                                             2.0 * g13b0 * g16b0 * covMatrix(13, 16) + 2.0 * g13b0 * g17b0 * covMatrix(13, 17) +

                                             2.0 * g14b0 * g15b0 * covMatrix(14, 15) + 2.0 * g14b0 * g16b0 * covMatrix(14, 16) +
                                             2.0 * g14b0 * g17b0 * covMatrix(14, 17) +

                                             2.0 * g15b0 * g16b0 * covMatrix(15, 16) + 2.0 * g15b0 * g17b0 * covMatrix(15, 17) +

                                             2.0 * g16b0 * g17b0 * covMatrix(16, 17)));
         EtapGaus_b0 /= binw;

         std::cout << " Check B0 Region" << "\n";
         std::cout << "Integral = " << Bgcount_Etap << "\n";
         std::cout << "Manual = " << iEtap_b0 << "\n";
         std::cout << "Error = " << EtapGaus_b0 << "\n";
         std::cout << "\n";


			/* Sideband 1 Region (left side). Signal only. */
			/* S_1 */
			// Etap Signal Region
			double SigSb1_Etap = fsig4->Integral(XEtapLower, XEtapLow) / binw;

			// Etap Signal SB1 Region
         highX = XEtapLow;
         lowX  = XEtapLower;
         ////////////////////////////////////////////////////
         /* Integral of Double Gaussian from lowX to highX */
         ////////////////////////////////////////////////////
         double iEtap_s1 =  -(sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))) +
                             (sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20])));
         iEtap_s1 /= binw;

         // First gaussian parameters
         // Partial derivitive of "ig" wrt p[18]
         double p0s1 =  -2.0 * sqrt(PI / 2.0) * p[18] * p[20] * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))
                       + 2.0 * sqrt(PI / 2.0) * p[18] * p[20] * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20]));

         // Partial derivitive of "ig" wrt p[19]
         double p1s1 =  pow(p[18], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[19] - lowX), 2.0)  / (2.0 * pow(p[20], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[19] - highX), 2.0) / (2.0 * pow(p[20], 2.0))))));

         // Partial derivitive of "ig" wrt p[20]
         double p2s1 =   pow(p[18], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20]))) -
                                            (sqrt(PI / 2.0) * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))) -
                                            (((p[19] - lowX)  * TMath::Power(TMath::E(), (-(pow((p[19] - lowX), 2.0))  / (2.0 * pow(p[20], 2.0))))) / p[20]) +
                                            (((p[19] - highX) * TMath::Power(TMath::E(), (-(pow((p[19] - highX), 2.0)) / (2.0 * pow(p[20], 2.0))))) / p[20]));

         double EtapGaus_s1 = sqrt(std::fabs(p0s1 * p0s1 * covMatrix(18, 18) + p1s1 * p1s1 * covMatrix(19, 19) + p2s1 * p2s1 * covMatrix(20, 20) +
                                             2.0 * p0s1 * p1s1 * covMatrix(18, 19) + 2.0 * p0s1 * p2s1 * covMatrix(18, 20) +
                                             2.0 * p1s1 * p2s1 * covMatrix(19, 20)));
         EtapGaus_s1 /= binw;

         std::cout << " Check S1 Region" << "\n";
         std::cout << "Integral = " << SigSb1_Etap << "\n";
         std::cout << "Manual = " << iEtap_s1 << "\n";
         std::cout << "Error = " << EtapGaus_s1 << "\n";
         std::cout << "\n";


			/* Sideband 2 Region (right side). Signal only. */
			/* S_2 */
			// Etap Signal SB2 Region
			double SigSb2_Etap = fsig4->Integral(XEtapHigh, XEtapHigher) / binw;

			// Etap Signal SB2 Region
         highX = XEtapHigher;
         lowX  = XEtapHigh;
         ////////////////////////////////////////////////////
         /* Integral of Double Gaussian from lowX to highX */
         ////////////////////////////////////////////////////
         double iEtap_s2 =  -(sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))) +
                             (sqrt(PI / 2.0) * pow(p[18], 2.0) * p[20] * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20])));
         iEtap_s2 /= binw;

         // First gaussian parameters
         // Partial derivitive of "ig" wrt p[18]
         double p0s2 =  -2.0 * sqrt(PI / 2.0) * p[18] * p[20] * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))
                       + 2.0 * sqrt(PI / 2.0) * p[18] * p[20] * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20]));

         // Partial derivitive of "ig" wrt p[19]
         double p1s2 =  pow(p[18], 2.0) * (TMath::Power(TMath::E(), (-(pow((p[19] - lowX), 2.0)  / (2.0 * pow(p[20], 2.0))))) - 
                                           TMath::Power(TMath::E(), (-(pow((p[19] - highX), 2.0) / (2.0 * pow(p[20], 2.0))))));

         // Partial derivitive of "ig" wrt p[20]
         double p2s2 =   pow(p[18], 2.0) * ((sqrt(PI / 2.0) * TMath::Erf((p[19] - lowX)  / (sqrt(2.0) * p[20]))) -
                                            (sqrt(PI / 2.0) * TMath::Erf((p[19] - highX) / (sqrt(2.0) * p[20]))) -
                                            (((p[19] - lowX)  * TMath::Power(TMath::E(), (-(pow((p[19] - lowX), 2.0))  / (2.0 * pow(p[20], 2.0))))) / p[20]) +
                                            (((p[19] - highX) * TMath::Power(TMath::E(), (-(pow((p[19] - highX), 2.0)) / (2.0 * pow(p[20], 2.0))))) / p[20]));


         double EtapGaus_s2 = sqrt(std::fabs(p0s2 * p0s2 * covMatrix(18, 18) + p1s2 * p1s2 * covMatrix(19, 19) + p2s2 * p2s2 * covMatrix(20, 20) +
                                             2.0 * p0s2 * p1s2 * covMatrix(18, 19) + 2.0 * p0s2 * p2s2 * covMatrix(18, 20) +
                                             2.0 * p1s2 * p2s2 * covMatrix(19, 20)));
         EtapGaus_s2 /= binw;

         std::cout << " Check S2 Region" << "\n";
         std::cout << "Integral = " << SigSb2_Etap << "\n";
         std::cout << "Manual = " << iEtap_s2 << "\n";
         std::cout << "Error = " << EtapGaus_s2 << "\n";
         std::cout << "\n";

			//Output file
			ofstream Calcstatus1;
         Calcstatus1.open("Etap_Calculations_Check_55MeV.txt", std::ios_base::app);
         if (Calcstatus1.is_open()) {
            if (tchannel == 0) {
               Calcstatus1 << "____________________________________" << "\n"; 
               Calcstatus1 << "XEtapUp     = " << XEtapUp * 1000.0 << " MeV/c^2" << "\n";
               Calcstatus1 << "XEtapDown   = " << XEtapDown * 1000.0 << " MeV/c^2" << "\n";
               Calcstatus1 << "XEtapLow    = " << XEtapLow * 1000.0 << " MeV/c^2" << "\n";
               Calcstatus1 << "XEtapLower  = " << XEtapLower * 1000.0 << " MeV/c^2" << "\n";
               Calcstatus1 << "XEtapHigh   = " << XEtapHigh * 1000.0 << " MeV/c^2" << "\n";
               Calcstatus1 << "XEtapHigher = " << XEtapHigher * 1000.0 << " MeV/c^2" << "\n";
            }
            Calcstatus1 << " Check S0 Region" << "\n";
            Calcstatus1 << "Integral = " << Sigcount_Etap << "\n";
            Calcstatus1 << "Manual   = " << iEtap_s0 << "\n";
            Calcstatus1 << "Error    = " << EtapGaus_s0 << "\n" << "\n";

            Calcstatus1 << " Check S1 Region" << "\n";
            Calcstatus1 << "Integral = " << SigSb1_Etap << "\n";
            Calcstatus1 << "Manual   = " << iEtap_s1 << "\n";
            Calcstatus1 << "Error    = " << EtapGaus_s1 << "\n" << "\n";

            Calcstatus1 << " Check S2 Region" << "\n";
            Calcstatus1 << "Integral = " << SigSb2_Etap << "\n";
            Calcstatus1 << "Manual   = " << iEtap_s2 << "\n";
            Calcstatus1 << "Error    = " << EtapGaus_s2 << "\n" << "\n";

            Calcstatus1 << " Check B0 Region" << "\n";
            Calcstatus1 << "Integral = " << Bgcount_Etap << "\n";
            Calcstatus1 << "Manual   = " << iEtap_b0 << "\n";
            Calcstatus1 << "Error    = " << EtapGaus_b0 << "\n" << "\n";

            Calcstatus1 << " Check B1 Region" << "\n";
            Calcstatus1 << "Integral = " << sb1count_Etap << "\n";
            Calcstatus1 << "Manual   = " << iEtap_b1 << "\n";
            Calcstatus1 << "Error    = " << EtapGaus_b1 << "\n" << "\n";

            Calcstatus1 << " Check B2 Region" << "\n";
            Calcstatus1 << "Integral = " << sb2count_Etap << "\n";
            Calcstatus1 << "Manual   = " << iEtap_b2 << "\n";
            Calcstatus1 << "Error    = " << EtapGaus_b2 << "\n" << "\n";
            Calcstatus1 << "---- " << "\n";
            Calcstatus1 << "---- " << "\n";
            Calcstatus1 << "\n";
         }
         else {
            Calcstatus1 << "Unable to open calculation check" << "\n";
         }
         Calcstatus1.close();


			/* RATIOS USED TO CALCULATE ACCIDENTAL SUBTRACTION WEIGHTS */
         // Eta Signal - Using manual calculations
			double R1_Etap   = iEtap_s1 / iEtap_s0;      // R1 = s1/s0
			double R2_Etap   = iEtap_s2 / iEtap_s0;      // R2 = s2/s0
			double C1_Etap   = iEtap_b1 / iEtap_b0;      // C1 = b1/b0
			double C2_Etap   = iEtap_b2 / iEtap_b0;      // C2 = b2/b0

			// Eta Signal - Using integrals
			double R1_Etap2  = SigSb1_Etap / Sigcount_Etap;     // R1 = s1/s0
			double R2_Etap2  = SigSb2_Etap / Sigcount_Etap;     // R2 = s2/s0
			double C1_Etap2  = sb1count_Etap / Bgcount_Etap;    // C1 = b1/b0
			double C2_Etap2  = sb2count_Etap / Bgcount_Etap;    // C2 = b2/b0


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
			// w2 = -w0 * (C2 / (C1*C1 + C2*C2))              //
			////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////
			//    Option of using equal sideband weight method, but we    //
			// use unequal sideband weight method for Pi0, Eta & Etaprime //
			////////////////////////////////////////////////////////////////

         /*
			/////////////////// For Eta, using equal sideband weights ///////////////////
			// w0 for Eta
			double w_etap_sig = (1.0 + R1_Etap + R2_Etap) / (1.0 - ((R1_Etap + R2_Etap) / (C1_Etap + C2_Etap)));
			// w1 = w2 for Eta
			double w_eta_sb  = -(w_etap_sig / (C1_Etap + C2_Etap));

         double w0_err = sqrt(std::fabs((R1_err * R1_err * (pow((C1_Etap + C2_Etap), 2.0)) * (pow((1.0 + C1_Etap + C2_Etap), 2.0)) +
                                         R2_err * R2_err * (pow((C1_Etap + C2_Etap), 2.0)) * (pow((1.0 + C1_Etap + C2_Etap), 2.0)) +
                                         C1_err * C1_err * (pow((R1_Etap + R2_Etap), 2.0)) * (pow((1.0 + R1_Etap + R2_Etap), 2.0)) +
                                         C2_err * C2_err * (pow((R1_Etap + R2_Etap), 2.0)) * (pow((1.0 + R1_Etap + R2_Etap), 2.0))) /
                                         (pow((R1_Etap + R2_Etap - C1_Etap - C2_Etap), 4.0))));

         double w1_err = sqrt(std::fabs(((R1_err * R1_err + R2_err * R2_err) * (pow((1.0 + C1_Etap + C2_Etap), 2.0)) +
                                         (C1_err * C1_err + C2_err * C2_err) * (pow((1.0 + R1_Etap + R2_Etap), 2.0))) /
                                        (pow((R1_Etap + R2_Etap - C1_Etap - C2_Etap), 4.0))));
         */

			/////////////////// For Etap, using manual calculations ///////////////////
			// w0 for Eta
			double w_etap_sig = (1.0 + R1_Etap + R2_Etap) / (1.0 - ((C1_Etap * R1_Etap + C2_Etap * R2_Etap) / (C1_Etap * C1_Etap + C2_Etap * C2_Etap)));
			// w1 for Eta
			double w_etap_sb1 = -(w_etap_sig) * (C1_Etap / (C1_Etap * C1_Etap + C2_Etap * C2_Etap));
			// w2 for Eta
			double w_etap_sb2 = -(w_etap_sig) * (C2_Etap / (C1_Etap * C1_Etap + C2_Etap * C2_Etap));


			/////////////////// For Etap, using integrals of the fit functions ///////////////////
			// w0 for Eta
			double w_etap2_sig = (1.0 + R1_Etap2 + R2_Etap2) / (1.0 - ((C1_Etap2 * R1_Etap2 + C2_Etap2 * R2_Etap2) / (C1_Etap2 * C1_Etap2 + C2_Etap2 * C2_Etap2)));
			// w1 for Eta
			double w_etap2_sb1 = -(w_etap2_sig) * (C1_Etap2 / (C1_Etap2 * C1_Etap2 + C2_Etap2 * C2_Etap2));
			// w2 for Eta
			double w_etap2_sb2 = -(w_etap2_sig) * (C2_Etap2 / (C1_Etap2 * C1_Etap2 + C2_Etap2 * C2_Etap2));


			/////////////////// For Eta, errors using manual calculations ///////////////////
         double R1_err = std::fabs(R1_Etap) * sqrt(pow((EtapGaus_s0 / iEtap_s0), 2.0) + pow((EtapGaus_s1 / iEtap_s1), 2.0));

         double R2_err = std::fabs(R2_Etap) * sqrt(pow((EtapGaus_s0 / iEtap_s0), 2.0) + pow((EtapGaus_s2 / iEtap_s2), 2.0));

         double C1_err = std::fabs(C1_Etap) * sqrt(pow((EtapGaus_b0 / iEtap_b0), 2.0) + pow((EtapGaus_b1 / iEtap_b1), 2.0));

         double C2_err = std::fabs(C2_Etap) * sqrt(pow((EtapGaus_b0 / iEtap_b0), 2.0) + pow((EtapGaus_b2 / iEtap_b2), 2.0));

			/////////////////// For Eta, errors using manual calculations ///////////////////
         double w0_err = sqrt(std::fabs((pow(C1_err, 2.0) * pow((2.0 * C1_Etap *(pow(C1_Etap, 2.0) - C1_Etap * R1_Etap + pow(C2_Etap, 2.0) - C2_Etap * R2_Etap) - (2.0 * C1_Etap - R1_Etap) * (pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0))), 2.0) * pow((R2_Etap + R1_Etap + 1.0), 2.0) +
                                         pow(C2_err, 2.0) * pow((2.0 * C2_Etap *(pow(C1_Etap, 2.0) - C1_Etap * R1_Etap + pow(C2_Etap, 2.0) - C2_Etap * R2_Etap) - (2.0 * C2_Etap - R2_Etap) * (pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0))), 2.0) * pow((R2_Etap + R1_Etap + 1.0), 2.0) +
                                         pow(R2_err, 2.0) * pow((pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0)), 2.0) * pow((pow(C1_Etap, 2.0) - (C1_Etap * R1_Etap) + pow(C2_Etap, 2.0) - C2_Etap * R2_Etap + C2_Etap * (R2_Etap + R1_Etap + 1.0)), 2.0) +
                                         pow(R1_err, 2.0) * pow((pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0)), 2.0) * pow((pow(C1_Etap, 2.0) - C1_Etap * R1_Etap + C1_Etap * (R2_Etap + R1_Etap + 1.0) + pow(C2_Etap, 2.0) - C2_Etap * R2_Etap), 2.0))
                                         / pow((pow(C1_Etap, 2.0) - C1_Etap * R1_Etap + pow(C2_Etap, 2.0) - C2_Etap * R2_Etap), 4.0)));

         double w1_err = sqrt(std::fabs((4.0 * pow(w_etap_sig, 2.0) * pow(C1_Etap, 2.0) * pow(C2_Etap, 2.0) * pow(C2_err, 2.0) + pow(w_etap_sig, 2.0) * pow(C1_err, 2.0) * pow((pow(C1_Etap, 2.0) - pow(C2_Etap, 2.0)), 2.0) + pow(C1_Etap, 2.0) * pow(w0_err, 2.0) * pow((pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0)), 2.0)) 
                                        / (pow((pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0)), 4.0))));

         double w2_err = sqrt(std::fabs((4.0 * pow(w_etap_sig, 2.0) * pow(C1_Etap, 2.0) * pow(C2_Etap, 2.0) * pow(C1_err, 2.0) + pow(w_etap_sig, 2.0) * pow(C2_err, 2.0) * pow((pow(C2_Etap, 2.0) - pow(C1_Etap, 2.0)), 2.0) + pow(C2_Etap, 2.0) * pow(w0_err, 2.0) * pow((pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0)), 2.0)) 
                                        / (pow((pow(C1_Etap, 2.0) + pow(C2_Etap, 2.0)), 4.0))));

         // Saving values for different |t| bins
         wR1_err[tchannel] = R1_err;
         wR2_err[tchannel] = R2_err;
         wC1_err[tchannel] = C1_err;
         wC2_err[tchannel] = C2_err;
         wR1_Etap[tchannel] = R1_Etap;
         wR2_Etap[tchannel] = R2_Etap;
         wC1_Etap[tchannel] = C1_Etap;
         wC2_Etap[tchannel] = C2_Etap;
         wiEtap_s0[tchannel] = iEtap_s0;
         wEtapGaus_s0[tchannel] = EtapGaus_s0;
         wiEtap_s1[tchannel] = iEtap_s1;
         wEtapGaus_s1[tchannel] = EtapGaus_s1;
         wiEtap_s2[tchannel] = iEtap_s2;
         wEtapGaus_s2[tchannel] = EtapGaus_s2;
         wiEtap_b0[tchannel] = iEtap_b0;
         wEtapGaus_b0[tchannel] = EtapGaus_b0;
         wiEtap_b1[tchannel] = iEtap_b1;
         wEtapGaus_b1[tchannel] = EtapGaus_b1;
         wiEtap_b2[tchannel] = iEtap_b2;
         wEtapGaus_b2[tchannel] = EtapGaus_b2;


         // Does not set a weight for null fits (< 5 mesons of interest in fit)
         if (!ZeroFit) {
            // Storing the |t| binned weights in order to add to a .root file for use in later STEP#'s
            w_etap_sigW[tchannel] = w_etap_sig;
            w0_errW[tchannel] = w0_err;

            w1_etap_sbW[tchannel] = w_etap_sb1;
            w1_errW[tchannel] = w1_err;

            w2_etap_sbW[tchannel] = w_etap_sb2;
            w2_errW[tchannel] = w2_err;
         }

			/////////////////////////////////////////////////
			/////////////// Fit Summary /////////////////////
			/////////////////////////////////////////////////
			// Pi0 Search
			c1->cd(1);
			gStyle->SetOptStat(0);
         newhfinal2_1->GetXaxis()->SetRangeUser(0.0, 0.25);
         PMax = newhfinal2_1->GetMaximum();
			newhfinal2_1->GetXaxis()->SetRangeUser(0.03, 0.26);
			newhfinal2_1->Draw("hist");
			newhfinal2_1->GetXaxis()->SetTitle("M_{2#gamma} [GeV/c^{2}]");
			newhfinal2_1->GetXaxis()->SetTitleOffset(1.2);
			newhfinal2_1->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
			newhfinal2_1->GetYaxis()->SetTitleOffset(1.5);
			newhfinal2_1->SetTitle("#pi^{0} Mass Range");
			newhfinal2_1->SetFillColor(18);
			fback->Draw("same");
			ftot->Draw("same");

			TLegend *legnd5 = new TLegend(0.6,0.68,0.89,0.89);
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

			newhfinal2_1b->Draw("hist");
			fback->Draw("same");
			fsig1->Draw("same");

			TLine *zerolineA = new TLine(Pidown, 0.0, Pidown, 0.5*PMax);
			zerolineA->SetLineColor(1);
			zerolineA->SetLineStyle(9);
			zerolineA->SetLineWidth(2);
			zerolineA->Draw("same");

			TLine *zerolineB = new TLine(Piup, 0.0, Piup, 0.5*PMax);
			zerolineB->SetLineColor(1);
			zerolineB->SetLineStyle(9);
			zerolineB->SetLineWidth(2);
			zerolineB->Draw("same");

			TLine *zerolineC = new TLine(PiLower, 0.0, PiLower, 0.25*PMax);
			zerolineC->SetLineColor(6);
			zerolineC->SetLineStyle(9);
			zerolineC->SetLineWidth(2);
			zerolineC->Draw("same");

			TLine *zerolineD = new TLine(PiHigher, 0.0, PiHigher, 0.25*PMax);
			zerolineD->SetLineColor(6);
			zerolineD->SetLineStyle(9);
			zerolineD->SetLineWidth(2);
			zerolineD->Draw("same");

         TLine *zline1232x = new TLine(PiLow, 0.0, PiLow, 0.25*PMax);
         zline1232x->SetLineColor(6);
         zline1232x->SetLineStyle(9);
         zline1232x->SetLineWidth(2);
         zline1232x->Draw("same");

         TLine *zline1232z = new TLine(PiHigh, 0.0, PiHigh, 0.25*PMax);
         zline1232z->SetLineColor(6);
         zline1232z->SetLineStyle(9);
         zline1232z->SetLineWidth(2);
         zline1232z->Draw("same");

         TLine *zline1232d = new TLine(mu_Pi, 0.0, mu_Pi, 1.025*PMax);
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
         legnd5b->AddEntry(zerolineA, Leg5b, "l");
         legnd5b->Draw("same");

         if (PeakVal >= 1.0) {
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

				TLine *zeroline = new TLine(0.03, 0.0, 0.26, 0.0);
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

				TLine *zerolineD2DD = new TLine(GausMax1, 0.0, GausMax1, 1.03*norm);
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
         }
         else {
				c1->cd(3);
				newhfinal2_1a->GetXaxis()->SetRangeUser(0.0, 0.001);
				newhfinal2_1a->SetTitle("Null Fit");
				newhfinal2_1a->Draw("hist");

				c1->cd(4);
				newhfinal2_1c->GetXaxis()->SetRangeUser(0.0, 0.001);
				newhfinal2_1c->SetTitle("Null Fit");
				newhfinal2_1c->Draw("hist");
         }

			char histo1[20];
			sprintf(histo1,"%s%d%s%d%s","Pi0_E", Echannel, "_t_", tchannel, ".png");
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

			TLegend *legnd5c = new TLegend(0.60,0.58,0.89,0.89); // (0.57,0.58,0.89,0.89);
			legnd5c->SetTextSize(0.028);
			legnd5c->AddEntry(newhfinal2_2b,"Real Data","f");
			legnd5c->AddEntry(fsig1,"#eta Signal Fit","l");
			legnd5c->AddEntry(fback,"Background Fit","l");
         TString Leg5c;
         Leg5c.Form("%.2f#sigma Signal Region", sigmaRANGE);
         legnd5c->AddEntry(zline3212, Leg5c, "l");
			legnd5c->Draw("same");

         if (PeakVal2 >= 1.0) {
				c2->cd(3);
				newhfinal2_2a->GetXaxis()->SetRangeUser(0.3, 0.8);
				newhfinal2_2a->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
				newhfinal2_2a->GetXaxis()->SetTitleOffset(1.2);
				newhfinal2_2a->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
				newhfinal2_2a->GetYaxis()->SetTitleOffset(1.5);
				newhfinal2_2a->SetTitle("Background in #eta Region");
				fback->SetFillColor(8);
				newhfinal2_2a->Add(newhfinal2_5a, -2);
				double mini2 = fback->Eval(0.3);
				mini2 += 200;
				newhfinal2_2a->SetMaximum(mini2*1.1);
				mini2 *= -0.08;
				newhfinal2_2a->SetMinimum(mini2);
				newhfinal2_2a->Draw("hist");
				fback->Draw("same");

				TLine *zeroline4 = new TLine(0.3, 0.0, 0.8, 0.0);
				zeroline4->SetLineColor(1);
				zeroline4->SetLineStyle(10);
				zeroline4->SetLineWidth(2);
				zeroline4->Draw("same");

				c2->cd(4);
				gStyle->SetOptFit(1100);
				newhfinal2_2c->Add(newhfinal2_7, 300);
				newhfinal2_2c->SetTitle("#eta Signal Fit");
				newhfinal2_2c->Draw("hist");
				newhfinal2_2c->GetXaxis()->SetRangeUser(0.30, 0.8);
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

         }
         else {
				c2->cd(3);
				newhfinal2_2a->GetXaxis()->SetRangeUser(0.0, 0.001);
				newhfinal2_2a->SetTitle("Null Fit");
				newhfinal2_2a->Draw("hist");

				c2->cd(4);
				newhfinal2_2c->GetXaxis()->SetRangeUser(0.0, 0.001);
				newhfinal2_2c->SetTitle("Null Fit");
				newhfinal2_2c->Draw("hist");
         }

			sprintf(histo1,"%s%d%s%d%s","Eta_E", Echannel, "_t_", tchannel, ".png");
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
			ftot->Draw("same");
			ftot->Draw("same");

			TLegend *legnd = new TLegend(0.6, 0.68, 0.89, 0.89);
			legnd->SetTextSize(0.028);
			legnd->AddEntry(newhfinal2_3, "Real Data", "f");
			legnd->AddEntry(ftot, "Signal & Background Fit", "l");
			legnd->AddEntry(fback, "Background Fit", "l");
			legnd->Draw("same");

			c3->cd(2);
			newhfinal2_3b->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
			newhfinal2_3b->GetXaxis()->SetTitleOffset(1.2);
			newhfinal2_3b->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
			newhfinal2_3b->GetYaxis()->SetTitleOffset(1.5);
			newhfinal2_3b->GetXaxis()->SetRangeUser(0.74, 1.18);
			newhfinal2_3b->SetFillColor(18);
			newhfinal2_3b->SetMinimum(0.0);
			newhfinal2_3b->Draw("hist");
			fback->Draw("same");
         fsig3->Draw("same");
			fsig4->Draw("same");

         TLine *zerolineA21 = new TLine(XEtapDown, 0.0, XEtapDown, 1.1*normEtap);
         zerolineA21->SetLineColor(kBlue);
         zerolineA21->SetLineStyle(9);
         zerolineA21->SetLineWidth(2);
         zerolineA21->Draw("same");

         TLine *zerolineB2 = new TLine(XEtapUp, 0.0, XEtapUp, 1.1*normEtap);
         zerolineB2->SetLineColor(kBlue);
         zerolineB2->SetLineStyle(9);
         zerolineB2->SetLineWidth(2);
         zerolineB2->Draw("same");

         TLine *zerolineC21 = new TLine(XEtapLower, 0.0, XEtapLower, 0.9*normEtap);
         zerolineC21->SetLineColor(2);
         zerolineC21->SetLineStyle(9);
         zerolineC21->SetLineWidth(2);
         zerolineC21->Draw("same");

         TLine *zerolineD21 = new TLine(XEtapHigher, 0.0, XEtapHigher, 0.9*normEtap);
         zerolineD21->SetLineColor(2);
         zerolineD21->SetLineStyle(9);
         zerolineD21->SetLineWidth(2);
         zerolineD21->Draw("same");

         TLine *zerolineE21 = new TLine(XEtapLow, 0.0, XEtapLow, 0.9*normEtap);
         zerolineE21->SetLineColor(2);
         zerolineE21->SetLineStyle(9);
         zerolineE21->SetLineWidth(2);
         zerolineE21->Draw("same");

         TLine *zerolineF21 = new TLine(XEtapHigh, 0.0, XEtapHigh, 0.9*normEtap);
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

         if (PeakVal3 >= 1.0) {
				c3->cd(3);
				newhfinal2_3a->GetXaxis()->SetRangeUser(0.74, 1.18);
				newhfinal2_3a->GetXaxis()->SetTitle("M_{2#gamma}[GeV/c^{2}]");
				newhfinal2_3a->GetXaxis()->SetTitleOffset(1.2);
				newhfinal2_3a->GetYaxis()->SetTitle("Events per 5 MeV/c^{2}");
				newhfinal2_3a->GetYaxis()->SetTitleOffset(1.5);
				newhfinal2_3a->SetTitle("Background in #eta' Region");
				fback->SetFillColor(8);
				newhfinal2_3a->Add(newhfinal2_5b, -2);
				double mini3a = fback->Eval(0.74);
				newhfinal2_3a->SetMaximum(mini3a*1.1);
				mini3a *= -0.08;
				newhfinal2_3a->SetMinimum(mini3a);
				newhfinal2_3a->Draw("hist");
				fback->Draw("same");

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

				TLine *zerolineD2A = new TLine(GausMax4, 0.0, GausMax4, 1.03*gryf);
				zerolineD2A->SetLineColor(6);
				zerolineD2A->SetLineStyle(9);
				zerolineD2A->SetLineWidth(2);
				zerolineD2A->Draw("same");

				double etapdiff = (abs(((GausMax4*1000) - 957.78) / 957.78))*100;
				TString titleetap;
				TString nameetap;
				nameetap.Form("#eta' Fit Peak = %.1f MeV/c^{2}", GausMax4*1000);
				titleetap.Form("Difference = %.2f %%", etapdiff);

				TLegend *legndAC = new TLegend(0.6, 0.68, 0.89, 0.89);
				legndAC->SetTextSize(0.028);
				legndAC->AddEntry(fsig4, "#eta' Signal Fit", "l");
				legndAC->AddEntry(zerolineD2A, nameetap, "l");
				legndAC->AddEntry((TObject*)0, "PDG Value = 957.8 MeV/c^{2}", "");
				legndAC->AddEntry((TObject*)0, titleetap, "");
				legndAC->Draw("same");
         }
         else {
				c3->cd(3);
				newhfinal2_3a->GetXaxis()->SetRangeUser(0.0, 0.001);
				newhfinal2_3a->SetTitle("Null Fit");
				newhfinal2_3a->Draw("hist");

				c3->cd(4);
				newhfinal2_3c->GetXaxis()->SetRangeUser(0.0, 0.001);
				newhfinal2_3c->SetTitle("Null Fit");
				newhfinal2_3c->Draw("hist");
         }

			sprintf(histo1,"%s%d%s%d%s","Etap_E", Echannel, "_t_", tchannel, ".png");
		   c3->Print(histo1);


         // Log plot of the total fit 
         TString TitleX;
         TitleX.Form("2#gamma Inv. Mass Fit in E_{#gamma} Range #rightarrow %d & |t| Range #rightarrow %d", Echannel, tchannel);

         TString histo10;
         histo10.Form("M2g_E%d_%d_LogPlot.png", Echannel, tchannel);
         c10->cd();
         gStyle->SetOptFit(0);
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
         int Etap1bin1 = hfillClone7->FindBin(XEtapLower);
         int Etap1bin2 = hfillClone7->FindBin(XEtapLow);

         // Etap Signal
         int Etap2bin1 = hfillClone8->FindBin(XEtapDown);
         int Etap2bin2 = hfillClone8->FindBin(XEtapUp);

         // Etap SB2
         int Etap3bin1 = hfillClone9->FindBin(XEtapHigh);
         int Etap3bin2 = hfillClone9->FindBin(XEtapHigher);

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
         c10->SetLogy(1);
         c10->Print(histo10);
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
         histo10.Form("M2g_E%d_%d_AllSignals.png", Echannel, tchannel);
         c11->Print(histo10);

			////////////////////////////////////////////////////////////////// 

			TF1 *fitresultTot = remass->GetFunction("ftot");

			double chisqTot = fitresultTot->GetChisquare();
			double ndofTot  = fitresultTot->GetNDF();
			double probTot  = fitresultTot->GetProb();
			int nbins = remass->GetNbinsX();

			std::cout << "\n";
			std::cout << "norm -> " << norm << "\n";
			std::cout << "bin width -> " << binw << "\n";

			std::cout << "\n";
			std::cout << "chisqTot           -> " << chisqTot << "\n";
			std::cout << "ndofTot            -> " << ndofTot << "\n";
			std::cout << "fit probabilityTot -> " << probTot << "\n"; 

			//Output file to ensure all fits converge
			ofstream Fitstatus1;
         TString StatusName;
         StatusName.Form("E%d_t%d", Echannel, tchannel);
         Fitstatus1.open("Etap_fit_status_55MeV.txt", std::ios_base::app);
         if (Fitstatus1.is_open()) {
            Fitstatus1 << StatusName << " = " << conv_fit << "\n";
            Fitstatus1 << "\n";
         }
         else {
            std::cout << "Unable to open fit status." << "\n";
         }
         Fitstatus1.close();

         // Output file with Fit Summary
         reset += 1;
			ofstream output;
			output.open("2g_Etap_Results_55MeV.txt", std::ios_base::app);
			if (output.is_open()) {
            if (reset == 1){
               output << "//-------------------------------//" << "\n";
               output << "// Energy Rescaling Factor       //" << "\n";
               output << "// ER = 0.97                     //" << "\n";
               output << "// Tagging Scaling Factor        //" << "\n";
               output << "// TSF = 0.935198 +/- 0.00249013 //" << "\n"; 
               output << "//-------------------------------//" << "\n";
               output << "\n";
            }

            output << "---- " << htitle << " ----" << "\n";
				output << "|t| bin = " << tchannel << "\n";

            int RAA = 0;
            int RBB = 0;
				TString title;
				if (Echannel == 0) {
					title.Form("%.3f < E_{gamma} < %.3f GeV", RangeB[0], RangeA[0]);
               RAA = TC_RangeA[0];
               RBB = TC_RangeB[0];
				}
				else if (Echannel == 1) {
					title.Form("%.3f < E_{gamma} < %.3f GeV", RangeB[1], RangeA[1]);
               RAA = TC_RangeA[1];
               RBB = TC_RangeB[1];
				}
				else if (Echannel == 2) {
					title.Form("%.3f < E_{gamma} < %.3f GeV", RangeB[2], RangeA[2]);
               RAA = TC_RangeA[2];
               RBB = TC_RangeB[2];
				}
				else {
					title.Form("%.3f < E_{gamma} < %.3f GeV", RangeB[3], RangeA[3]);
               RAA = TC_RangeA[3];
               RBB = TC_RangeB[3];
				}

				output << "|t| = " << Etap_bins[tchannel][0] << " to " << Etap_bins[tchannel][1] << " GeV^2" << "\n";
            output << "E bin = " << Echannel << ",  Tagger channel = " << RAA << " to " << RBB << "\n";
				output << title << "\n";
				output << "\n";
				output << "chisqTot           -> " << chisqTot << "\n";
				output << "ndofTot            -> " << ndofTot << "\n";
				output << "fit probabilityTot -> " << probTot << "\n";
				output << "\n";
				output << "Pi0 Fit's Peak           = " << GausMax1 * 1000.0 << " MeV/c^2" << "\n";
				output << "Eta Fit's Peak           = " << GausMax2 * 1000.0 << " MeV/c^2" << "\n";
				output << "Omega Leakage Fit's Peak = " << GausMax3 * 1000.0 << " MeV/c^2" << "\n";
				output << "Etaprime Fit's Peak      = " << GausMax4 * 1000.0 << " MeV/c^2" << "\n";
				output << "\n";

				TString Verg;
				if (conv_fit == 1) {
					Verg.Form("Converged");
				}
				else {
					Verg.Form("Failed");
				}

            output << "Status = " << Verg << "\n";
				output << "\n";
				output << "chisqTot           -> " << chisqTot << "\n";
				output << "ndofTot            -> " << ndofTot << "\n";
				output << "fit probabilityTot -> " << probTot << "\n";
				output << "\n";

				output << "Pi0_____mean_1 -> " << ftot->GetParameter(1);
            if (((mean_sig1 - 0.00000001) < dparlim1) || ((mean_sig1 + 0.00000001) > uparlim1)) {
				   output << "  ****";
            }
				output << "\n";

				output << "Pi0____sigma_1 -> " << ftot->GetParameter(2);
            if (((sigma_sig1 - 0.00000001) < dparlim2) || ((sigma_sig1 + 0.00000001) > uparlim2)) {
				   output << "  ****";
            }
				output << "\n";

				output << "Pi0_____mean_2 -> " << ftot->GetParameter(4);
            if (((mean2_sig1 - 0.00000001) < dparlim4) || ((mean2_sig1 + 0.00000001) > uparlim4)) {
				   output << "  ****";
            }
				output << "\n";

				output << "Pi0____sigma_2 -> " << ftot->GetParameter(5);
            if (((sigma2_sig1 - 0.00000001) < dparlim5) || ((sigma2_sig1 + 0.00000001) > uparlim5)) {
				   output << "  ****";
            }
				output << "\n";
				output << "\n";

				output << "Eta_____mean_1 -> " << ftot->GetParameter(7);
            if (((mean_sig2 - 0.00000001) < dparlim7) || ((mean_sig2 + 0.00000001) > uparlim7)) {
				   output << "  ****";
            }
				output << "\n";

				output << "Eta____sigma_1 -> " << ftot->GetParameter(8);
            if (((sigma_sig2 - 0.00000001) < dparlim8) || ((sigma_sig2 + 0.00000001) > uparlim8)) {
				   output << "  ****";
            }
				output << "\n";

				output << "Eta_____mean_2 -> " << ftot->GetParameter(10);
            if (((mean2_sig2 - 0.00000001) < dparlim10) || ((mean2_sig2 + 0.00000001) > uparlim10)) {
				   output << "  ****";
            }
				output << "\n";

				output << "Eta____sigma_2 -> " << ftot->GetParameter(11);
            if (((sigma2_sig2 - 0.00000001) < dparlim11) || ((sigma2_sig2 + 0.00000001) > uparlim11)) {
				   output << "  ****";
            }
				output << "\n";
				output << "\n";

				output << "3G_leak_mean_1 -> " << ftot->GetParameter(13);
            if (((mean_sig3 - 0.00000001) < dparlim13) || ((mean_sig3 + 0.00000001) > uparlim13)) {
				   output << "  ****";
            }
				output << "\n";
				output << "3G_leak_sigma1 -> " << ftot->GetParameter(14);
            if (((sigma_sig3 - 0.00000001) < dparlim14) || ((sigma_sig3 + 0.00000001) > uparlim14)) {
				   output << "  ****";
            }
				output << "\n";
				output << "3G_leak_mean_2 -> " << ftot->GetParameter(16);
            if (((mean2_sig3 - 0.00000001) < dparlim16) || ((mean2_sig3 + 0.00000001) > uparlim16)) {
				   output << "  ****";
            }
				output << "\n";
				output << "3G_leak_sigma2 -> " << ftot->GetParameter(17);
            if (((sigma2_sig3 - 0.00000001) < dparlim17) || ((sigma2_sig3 + 0.00000001) > uparlim17)) {
				   output << "  ****";
            }
				output << "\n";
				output << "\n";

				output << "Etap____mean_1 -> " << ftot->GetParameter(19);
            if (((mean_sig4 - 0.00000001) < dparlim19) || ((mean_sig4 + 0.00000001) > uparlim19)) {
				   output << "  ****";
            }
				output << "\n";
				output << "Etap___sigma_1 -> " << ftot->GetParameter(20);
            if (((sigma_sig4 - 0.00000001) < dparlim20) || ((sigma_sig4 + 0.00000001) > uparlim20)) {
				   output << "  ****";
            }
				output << "\n";
				output << "\n";
				output << "p0     -> " << ftot->GetParameter(21) << "\n";
				output << "p1     -> " << ftot->GetParameter(22) << "\n";
				output << "p2     -> " << ftot->GetParameter(23) << "\n";
				output << "p3     -> " << ftot->GetParameter(24) << "\n";
				output << "Scaler -> " << ftot->GetParameter(25) << "\n";
				output << "_____________________________________________" << "\n";
				output << "\n";
         }
			else {
				std::cout << "UNABLE TO OPEN OUTPUT FILE." << "\n";
			}
			output.close();

         std::cout << "*********** Ending ***********" << "\n";
         std::cout << "        " << htitle << "        " << "\n";
         std::cout << "******************************" << "\n";

         // This stops the fitting program if a fit does not converge.
         // It allows us to determine the first failed fit.
         //
         // Go to the "IF statements" at around line 356,
         // and modify "SetParLimits" to get the fit in question to converge.
         if (conv_fit != 1) {
            break;
         }
      }

      if (conv_fit != 1) {
         delete c1;
         delete c2;
         delete c3;
         delete c10;
         delete c11;
         delete c123;
         delete c321;
         delete c555;
         std::cout << "FITTING FAILED AT: E" << Echannel_count << " t" << tchannel_count << "\n";
         break;
      }

      // Root file with new histograms for reference in later analysis work
      TFile *MyFile = new TFile("TwoG_may_16_2023_STEP3_55MeV.root", "UPDATE");

      TString hist0Name = Form("h_W0_E%d", Echannel);
      // w0 weight and error
      // bin 1 = w0, bin 2 = error
      TH1D *hEtap_W0 = new TH1D(hist0Name, "#eta' Signal Region Sideband Subtraction Weight;Weight #pm Error [1^{st} Bin #pm 2^{nd} Bin];Value", 3, 0, 3);
      //hEtap_W0->Sumw2();

      // w1 weight and error
      // bin 1 = w1, bin 2 = error
      TString hist1Name = Form("h_W1_E%d", Echannel);
      TH1D *hEtap_W1 = new TH1D(hist1Name, "#eta' Lower Mass Region Sideband Subtraction Weight;Weight #pm Error [1^{st} Bin #pm 2^{nd} Bin];Value", 3, 0, 3);
      //hEtap_W1->Sumw2();

      // w2 weight and error
      // bin 1 = w2, bin 2 = error
      TString hist2Name = Form("h_W2_E%d", Echannel);
      TH1D *hEtap_W2 = new TH1D(hist2Name, "#eta' Higher Mass Region Sideband Subtraction Weight;Weight #pm Error [1^{st} Bin #pm 2^{nd} Bin];Value", 3, 0, 3);
      //hEtap_W2->Sumw2();

      for (int yurt = 0; yurt < tBinNum; yurt++) {
         int acq = yurt + 1;
         // Prevents a histogram entry for null fits
         if (ZeroFitArray[yurt] == 1.0) {
            continue;
         }
         hEtap_W0->SetBinContent(acq, w_etap_sigW[yurt]);
         hEtap_W0->SetBinError(acq, w0_errW[yurt]);

         hEtap_W1->SetBinContent(acq, w1_etap_sbW[yurt]);
         hEtap_W1->SetBinError(acq, w1_errW[yurt]);

         hEtap_W2->SetBinContent(acq, w2_etap_sbW[yurt]);
         hEtap_W2->SetBinError(acq, w2_errW[yurt]);
      }
      gFile = MyFile;
      gDirectory->WriteObject(hEtap_W0, hist0Name);
      gDirectory->WriteObject(hEtap_W1, hist1Name);
      gDirectory->WriteObject(hEtap_W2, hist2Name);
      MyFile->Close();


      //Output file
      TString WgtName;
      WgtName.Form("Errors_55MeV.txt");
      ofstream FitPar2;
      FitPar2.open(WgtName, std::ios_base::app);
      if (FitPar2.is_open()) {
         FitPar2 << "E channel = " << Echannel << "\n";
         FitPar2 << "_____________________________________" << "\n";
         for (int iterT = 0; iterT < tBinNum; iterT++) {
            FitPar2 << "*** |t| channel = " << iterT << " ***" << "\n";
            FitPar2 << "w_etap_sig  = " << w_etap_sigW[iterT] << "\n";
            FitPar2 << "w0_err      = " << w0_errW[iterT] << "\n";
            FitPar2 << "w_etap_sb1  = " << w1_etap_sbW[iterT] << "\n";
            FitPar2 << "w1_err      = " << w1_errW[iterT] << "\n";
            FitPar2 << "w_etap_sb2  = " << w2_etap_sbW[iterT] << "\n";
            FitPar2 << "w2_err      = " << w2_errW[iterT] << "\n";

            FitPar2 << "R1_Etap     = " << wR1_Etap[iterT] << "\n";
            FitPar2 << "R1_err      = " << wR1_err[iterT] << "\n";
            FitPar2 << "R2_Etap     = " << wR2_Etap[iterT] << "\n";
            FitPar2 << "R2_err      = " << wR2_err[iterT] << "\n";
            FitPar2 << "C1_Etap     = " << wC1_Etap[iterT] << "\n";
            FitPar2 << "C1_err      = " << wC1_err[iterT] << "\n";
            FitPar2 << "C2_Etap     = " << wC2_Etap[iterT] << "\n";
            FitPar2 << "C2_err      = " << wC2_err[iterT] << "\n";

            FitPar2 << "iEtap_s0    = " << wiEtap_s0[iterT] << "\n";
            FitPar2 << "EtapGaus_s0 = " << wEtapGaus_s0[iterT] << "\n";
            FitPar2 << "iEtap_s1    = " << wiEtap_s1[iterT] << "\n";
            FitPar2 << "EtapGaus_s1 = " << wEtapGaus_s1[iterT] << "\n";
            FitPar2 << "iEtap_s2    = " << wiEtap_s2[iterT] << "\n";
            FitPar2 << "EtapGaus_s2 = " << wEtapGaus_s2[iterT] << "\n";
            FitPar2 << "iEtap_b0    = " << wiEtap_b0[iterT] << "\n";
            FitPar2 << "EtapGaus_b0 = " << wEtapGaus_b0[iterT] << "\n";
            FitPar2 << "iEtap_b1    = " << wiEtap_b1[iterT] << "\n";
            FitPar2 << "EtapGaus_b1 = " << wEtapGaus_b1[iterT] << "\n";
            FitPar2 << "iEtap_b2    = " << wiEtap_b2[iterT] << "\n";
            FitPar2 << "EtapGaus_b2 = " << wEtapGaus_b2[iterT] << "\n";
            FitPar2 << "-------------------" << "\n";
            FitPar2 << "\n";
         }
         FitPar2 << "\n";
      }
   }

   // Root file with new histograms for reference in later analysis work
   TFile *yieldFile = new TFile("TwoG_may_16_2023_STEP3_55MeV.root", "UPDATE");

   TH1D *Etap_yield = new TH1D("Etap_yield", "#eta' yield for fits with a minimum parameter width limit of 40 MeV/c^{2};x_{_{0-1 }}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }0) ,  x_{_{1-2}}= (E_{#gamma_{Bin}}=_{ }0, |t|_{_{Bin}}=_{ }1) ,  ... , x_{_{11-12}}= (E_{#gamma_{Bin}}=_{ }3, |t|_{_{Bin}}=_{ }2);#eta' yield per (E_{#gamma}, |t|) Bin", 12, 0, 12);
   //Etap_yield->Sumw2();


   for (int k = 0; k < histoCount; k++) {
      int abq = k + 1;
      // Leave bin empty for null fits
      if (ZFA[k] == 1.0) {
         continue;
      }
      // Fill histo
      Etap_yield->SetBinContent(abq, EtapYield[k]);
      Etap_yield->SetBinError(abq, EtapYieldError[k]);
   }
   gFile = yieldFile;
   gDirectory->WriteObject(Etap_yield, "Etap_yield");
   yieldFile->Close();

   if (conv_fit == 1) {
      delete c1;
      delete c2;
      delete c3;
      delete c10;
      delete c11;
      delete c123;
      delete c321;
      delete c555;

      std::cout << "FITTING COMPLETED" << "\n";

      // Organize folder info
      gSystem->Exec("mkdir -p calc_check");
      gSystem->Exec("mkdir -p histos/Pi0");
      gSystem->Exec("mkdir -p histos/Eta");
      gSystem->Exec("mkdir -p histos/Etap");
      gSystem->Exec("mkdir -p histos/LogPlot");
      gSystem->Exec("mkdir -p histos/AllSignals");
      gSystem->Exec("mkdir -p FitSummary/Pi0");
      gSystem->Exec("mkdir -p FitSummary/Eta");
      gSystem->Exec("mkdir -p FitSummary/Etap");

      gSystem->Exec("mv Pi0_*.png histos/Pi0/.");
      gSystem->Exec("mv Eta_*.png histos/Eta/.");
      gSystem->Exec("mv Etap_*.png histos/Etap/.");
      gSystem->Exec("mv *_LogPlot.png histos/LogPlot/.");
      gSystem->Exec("mv *_AllSignals.png histos/AllSignals/.");
      gSystem->Exec("mv PiFit_E*.png FitSummary/Pi0/.");
      gSystem->Exec("mv EtaFit_E*.png FitSummary/Eta/.");
      gSystem->Exec("mv EtapFit_E*.png FitSummary/Etap/.");
      gSystem->Exec("mv *Results_55MeV.txt FitSummary/.");
      gSystem->Exec("mv *fit_status_55MeV.txt FitSummary/.");

      gSystem->Exec("mv Etap_Calculations_Check_55MeV.txt calc_check/.");
      gSystem->Exec("mv Etaprime_yield_check_55MeV.txt calc_check/.");
      gSystem->Exec("mv Errors_55MeV.txt calc_check/.");
      gSystem->Exec("mv Etaprime_yield_fit_width_55MeV_min.txt calc_check/.");
   }

   // Failed Fitting
   else {
      delete c1;
      delete c2;
      delete c3;
      delete c10;
      delete c11;
      delete c123;
      delete c321;
      delete c555;

      std::cout << "FITTING FAILED" << "\n";

      // Organize folder info
      gSystem->Exec("mkdir -p FAIL/TRASH");
      gSystem->Exec("mkdir -p FAIL/TRASH/AllSignals/.");
      gSystem->Exec("mkdir -p FAIL/TRASH/Pi0/.");
      gSystem->Exec("mkdir -p FAIL/TRASH/Eta/.");
      gSystem->Exec("mkdir -p FAIL/TRASH/Etap/.");

      gSystem->Exec("mv Pi0_*.png FAIL/TRASH/.");
      gSystem->Exec("mv Eta_*.png FAIL/TRASH/.");
      gSystem->Exec("mv Etap_*.png FAIL/TRASH/.");
      gSystem->Exec("mv *_LogPlot.png FAIL/TRASH/.");
      gSystem->Exec("mv *_AllSignals.png FAIL/TRASH/AllSignals/.");
      gSystem->Exec("mv PiFit_E*.png FAIL/TRASH/Pi0/.");
      gSystem->Exec("mv EtaFit_E*.png FAIL/TRASH/Eta/.");
      gSystem->Exec("mv EtapFit_E*.png FAIL/TRASH/Etap/.");
      gSystem->Exec("mv *Results_55MeV.txt FAIL/.");
      gSystem->Exec("mv *fit_status_55MeV.txt FAIL/TRASH/.");
      gSystem->Exec("mv *.txt FAIL/TRASH/.");
   }
}
