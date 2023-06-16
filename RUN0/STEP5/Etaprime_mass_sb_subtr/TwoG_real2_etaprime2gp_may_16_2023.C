#define TwoG_real2_etaprime2gp_may_16_2023_cxx

// Author: J. McIntyre
// May 16, 2023
// This file contains cuts on real data for a 2 photon sample.
////////////////////////////////////////////////////////////////////////
// This was written to look at the 2 gamma sample.                    //
// __________________________________________________________________ //
// __________________________________________________________________ //
//                                                                    //
// THIS FILE HAS UPDATES THAT:                                        //
// --------------------------                                         //
//    - USES THE LATEST ENERGY RESCALING FACTOR (ER) = 0.97           //
// __________________________________________________________________ //
//                                                                    //
// This file is a re-do of ER_real_etaprime2gp_oct_21_2021.C,         //
// /nfs/direct/annex/mcintyre/analysisIII/ER/ER_0_955/5MeV_mass_bins  //
// since the ER & ASF were updated and a valid cross section for pi0  //
// was plotted. Now this file is the start of a final analysis of the //
// two gamma sample for the pi0, eta, and etap cross sections.        //
//                                                                    //
// The new ER equals 0.97, as determined by the various attempts in   //
// folder: /nfs/direct/annex/mcintyre/analysisIX/TwoGamma/inv_mass    //
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*
   The following methods are defined in this file:
      Begin():        called every time a loop on the tree starts,
                      a convenient place to create your histograms.
      SlaveBegin():   called after Begin(), when on PROOF called only on the
                      slave servers.
      Process():      called for each event, in this function you decide what
                      to read and fill your histograms.
      SlaveTerminate: called at the end of the loop on the tree, when on PROOF
                      called only on the slave servers.
      Terminate():    called at the end of the loop on the tree,
                      a convenient place to draw/fit your histograms.
*/

// Regarding #include, using <> tells the compiler to search for the .h file in the system include directories,
// while using "" has it search for the .h file in the current directory and then in the system include directories.
#include "TwoG_real2_etaprime2gp_may_16_2023.h"
#include <iostream>
#include <TProofServ.h>
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
#include "TMinuit.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TSpline.h>
#include <TH2.h>

/*
#include <TF2.h>
#include <TGraphSmooth.h>
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

#define sqr(A) ((A)*(A))     // Define a function to square a number

#define PI 3.14159265359

#define fwdgamma 2           // Defines the number of photons in 
                             // the LGD for the events of interest
#define bgvgamma 2           // Defines how many photons are permitted in the BGV
                             // [bgvgamma - fwdgamma] = photons_in_BGV

#define ResonanceEner 0.300  // Cut on missing energy (e.g. resonance energy of the proton)(GeV)

#define ASF 0.935198         // Accidental Scaling Factor -- **pass-9-2020 dataset** Average ASF = 0.935198 +/- 0.00249013 
                             // (*old* pass-1-2020 -> 0.935144 +/- 0.00162298) 

#define ER 0.97              // Energy Rescaler -- Shifts invariant mass plot to 
                             // match well known mesons
                             // Previous tries: 0.955, 0.953, 0.9668

#define TYPECUT 1            // Cluster_cleanup cut

// using namespace std;         // Removes the need for std::

TFile *fHistFile;            // Set Pointer to the output file


void TwoG_real2_etaprime2gp_may_16_2023::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void TwoG_real2_etaprime2gp_may_16_2023::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();

   // Get previous analysis procedure values from STEP2 root file
   TFile *rootf = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP2/total_fit/TwoG_may_16_2023_STEP2.root");

   // |t| bin values
   hBinDW = (TH1D*)rootf->Get("hBinDW");
   hBinDW->SetDirectory(0);
   hBinUP = (TH1D*)rootf->Get("hBinUP");
   hBinUP->SetDirectory(0);

   // Information about photon energy binning
   EChanBINS = (TH1D*)rootf->Get("EChanBINS");
   EChanBINS->SetDirectory(0);
   EChanRange = (TH1D*)rootf->Get("EChanRange");
   EChanRange->SetDirectory(0);
   EChanValue = (TH1D*)rootf->Get("EChanValue");
   EChanValue->SetDirectory(0);

   // Mass values of the sideband subtraction regions
   SB_Regions_allt_E0 = (TH1D*)rootf->Get("SB_Regions_allt_E0");
   SB_Regions_allt_E0->SetDirectory(0);
   SB_Regions_allt_E1 = (TH1D*)rootf->Get("SB_Regions_allt_E1");
   SB_Regions_allt_E1->SetDirectory(0);
   SB_Regions_allt_E2 = (TH1D*)rootf->Get("SB_Regions_allt_E2");
   SB_Regions_allt_E2->SetDirectory(0);
   SB_Regions_allt_E3 = (TH1D*)rootf->Get("SB_Regions_allt_E3");
   SB_Regions_allt_E3->SetDirectory(0);

   rootf->Close();

   // Get previous analysis procedure TSplines used for
   // sideband subtraction weights from STEP4 root file
   TFile *rootS4f = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_etaprime2gp/RUN0/STEP4/Etaprime_mass_SBwgt/TwoG_may_16_2023_STEP4.root");

   // TSpline does not need SetDirectory(0) when TFile is closed
   // SIGNAL REGION
   w0_spline_E0 = (TSpline*)rootS4f->Get("w0_spline_E0");
   w0_spline_E1 = (TSpline*)rootS4f->Get("w0_spline_E1");
   w0_spline_E2 = (TSpline*)rootS4f->Get("w0_spline_E2");
   w0_spline_E3 = (TSpline*)rootS4f->Get("w0_spline_E3");
   // LOWER MASS SIDEBAND REGION
   w1_spline_E0 = (TSpline*)rootS4f->Get("w1_spline_E0");
   w1_spline_E1 = (TSpline*)rootS4f->Get("w1_spline_E1");
   w1_spline_E2 = (TSpline*)rootS4f->Get("w1_spline_E2");
   w1_spline_E3 = (TSpline*)rootS4f->Get("w1_spline_E3");
   //HIGHER MASS SIDEBAND REGION
   w2_spline_E0 = (TSpline*)rootS4f->Get("w2_spline_E0");
   w2_spline_E1 = (TSpline*)rootS4f->Get("w2_spline_E1");
   w2_spline_E2 = (TSpline*)rootS4f->Get("w2_spline_E2");
   w2_spline_E3 = (TSpline*)rootS4f->Get("w2_spline_E3");

   rootS4f->Close();


   // Histogram Definitions
   tEtap0 = new TH1D("tEtap0", "|t| from E_{#gamma} = 5.115 to 5.385 GeV for #eta';|t| [GeV^{2}/c^{4]];Events per 5 MeV^{2}", 380, 0.0, 1.9);
   tEtap0->Sumw2();
   GetOutputList()->Add(tEtap0);

   dphi_Etap0 = new TH1D("dphi_Etap0", "#Delta#phi from E_{#gamma} = 5.115 to 5.385 GeV for #eta' all |t|;#Delta#phi [degree];Events per 2 deg.", 180, -180.0, 180.0);
   dphi_Etap0->Sumw2();
   GetOutputList()->Add(dphi_Etap0);

   tEtap1 = new TH1D("tEtap1", "|t| from E_{#gamma} = 4.860 to 5.115 GeV for #eta';|t| [GeV^{2}/c^{4]];Events per 5 MeV^{2}", 380, 0.0, 1.9);
   tEtap1->Sumw2();
   GetOutputList()->Add(tEtap1);

   dphi_Etap1 = new TH1D("dphi_Etap1", "#Delta#phi from E_{#gamma} = 4.860 to 5.115 GeV for #eta' all |t|;#Delta#phi [degree];Events per 2 deg.", 180, -180.0, 180.0);
   dphi_Etap1->Sumw2();
   GetOutputList()->Add(dphi_Etap1);

   tEtap2 = new TH1D("tEtap2", "|t| from E_{#gamma} = 4.615 to 4.860 GeV for #eta';|t| [GeV^{2}/c^{4]];Events per 5 MeV^{2}", 380, 0.0, 1.9);
   tEtap2->Sumw2();
   GetOutputList()->Add(tEtap2);

   dphi_Etap2 = new TH1D("dphi_Etap2", "#Delta#phi from E_{#gamma} = 4.615 to 4.860 GeV for #eta' all |t|;#Delta#phi [degree];Events per 2 deg.", 180, -180.0, 180.0);
   dphi_Etap2->Sumw2();
   GetOutputList()->Add(dphi_Etap2);

   tEtap3 = new TH1D("tEtap3", "|t| from E_{#gamma} = 4.358 to 4.615 GeV for #eta';|t| [GeV^{2}/c^{4]];Events per 5 MeV^{2}", 380, 0.0, 1.9);
   tEtap3->Sumw2();
   GetOutputList()->Add(tEtap3);

   dphi_Etap3 = new TH1D("dphi_Etap3", "#Delta#phi from E_{#gamma} = 4.358 to 4.615 GeV for #eta' all |t|;#Delta#phi [degree];Events per 2 deg.", 180, -180.0, 180.0);
   dphi_Etap3->Sumw2();
   GetOutputList()->Add(dphi_Etap3);

   // |t| Bin Range Values
   // Ensure this matches "BinM2g[3][2]" line 219 of "Ebin_allt_Etap0_BkgdFit.C" in STEP2
   double XEtap_bins[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

   TString title;
   TString name;
	for (int i = 0; i < 3; i++) {
		name.Form("Etap_dphi_E0_t%d", i);
		title.Form("#Delta#phi from E_{#gamma} = 5.115 to 5.385 GeV for #eta' fit, %.3f <= |t| < %.3f;#Delta#phi [degree];Events per 2 deg.", XEtap_bins[i][0], XEtap_bins[i][1]);
		Etap_dphi_E0_t[i] = new TH1D(name, title, 180, -180.0, 180.0);
		Etap_dphi_E0_t[i]->Sumw2();
		GetOutputList()->Add(Etap_dphi_E0_t[i]);

		name.Form("Etap_dphi_E1_t%d", i);
		title.Form("#Delta#phi from E_{#gamma} = 4.860 to 5.115 GeV for #eta' fit, %.3f <= |t| < %.3f;#Delta#phi [degree];Events per 2 deg.", XEtap_bins[i][0], XEtap_bins[i][1]);
		Etap_dphi_E1_t[i] = new TH1D(name, title, 180, -180.0, 180.0);
		Etap_dphi_E1_t[i]->Sumw2();
		GetOutputList()->Add(Etap_dphi_E1_t[i]);

		name.Form("Etap_dphi_E2_t%d", i);
		title.Form("#Delta#phi from E_{#gamma} = 4.615 to 4.860 GeV for #eta' fit, %.3f <= |t| < %.3f;#Delta#phi [degree];Events per 2 deg.", XEtap_bins[i][0], XEtap_bins[i][1]);
		Etap_dphi_E2_t[i] = new TH1D(name, title, 180, -180.0, 180.0);
		Etap_dphi_E2_t[i]->Sumw2();
		GetOutputList()->Add(Etap_dphi_E2_t[i]);

		name.Form("Etap_dphi_E3_t%d", i);
		title.Form("#Delta#phi from E_{#gamma} = 4.358 to 4.615 GeV for #eta' fit, %.3f <= |t| < %.3f;#Delta#phi [degree];Events per 2 deg.", XEtap_bins[i][0], XEtap_bins[i][1]);
		Etap_dphi_E3_t[i] = new TH1D(name, title, 180, -180.0, 180.0);
		Etap_dphi_E3_t[i]->Sumw2();
		GetOutputList()->Add(Etap_dphi_E3_t[i]);
   }
}


Bool_t TwoG_real2_etaprime2gp_may_16_2023::Process(Long64_t entry)
{
   if (fChain == 0) {   // REMEMBER:
      return false;     // false is valid for ROOT 6 and above.
   }                    // kFALSE works for all ROOT, but not C++.

   GetEntry(entry);

   // Inputs from fitting
   /*
   // Get |t| bin ranges from two histos created in earlier analysis
   double Etap_bins[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};
   */

   // Get number of |t| bins used
   int tBinNum = hBinUP->GetNbinsX();

   // Get upper and lower limits of each |t| bin
   /* double Etap_bins[tBinNum][2]; */
   double Etap_bins[3][2];
   for (int b1 = 0; b1 < tBinNum; b1++) {
      int bs = b1 + 1;
      Etap_bins[b1][0] = hBinDW->GetBinContent(bs);
      Etap_bins[b1][1] = hBinUP->GetBinContent(bs);
   }

   // Number of tagger energy bins
   int NumEChan = EChanBINS->GetBinContent(1);

   // Eta Sideband Subtraction Regions
   // Import values for Eta meson Sideband Subtraction Regions
   // I use the sideband subtraction regions determined from fitting
   // over all |t| in each tagger range in STEP2
   double EtaUp[4];
   double EtaDown[4];
   double EtaLow[4];
   double EtaLower[4];
   double EtaHigh[4];
   double EtaHigher[4];
	EtaLower[0]  = SB_Regions_allt_E0->GetBinContent(1);
	EtaLow[0]    = SB_Regions_allt_E0->GetBinContent(2);
	EtaDown[0]   = SB_Regions_allt_E0->GetBinContent(3);
	EtaUp[0]     = SB_Regions_allt_E0->GetBinContent(4);
	EtaHigh[0]   = SB_Regions_allt_E0->GetBinContent(5);
	EtaHigher[0] = SB_Regions_allt_E0->GetBinContent(6);

	EtaLower[1]  = SB_Regions_allt_E1->GetBinContent(1);
	EtaLow[1]    = SB_Regions_allt_E1->GetBinContent(2);
	EtaDown[1]   = SB_Regions_allt_E1->GetBinContent(3);
	EtaUp[1]     = SB_Regions_allt_E1->GetBinContent(4);
	EtaHigh[1]   = SB_Regions_allt_E1->GetBinContent(5);
	EtaHigher[1] = SB_Regions_allt_E1->GetBinContent(6);

	EtaLower[2]  = SB_Regions_allt_E2->GetBinContent(1);
	EtaLow[2]    = SB_Regions_allt_E2->GetBinContent(2);
	EtaDown[2]   = SB_Regions_allt_E2->GetBinContent(3);
	EtaUp[2]     = SB_Regions_allt_E2->GetBinContent(4);
	EtaHigh[2]   = SB_Regions_allt_E2->GetBinContent(5);
	EtaHigher[2] = SB_Regions_allt_E2->GetBinContent(6);

	EtaLower[3]  = SB_Regions_allt_E3->GetBinContent(1);
	EtaLow[3]    = SB_Regions_allt_E3->GetBinContent(2);
	EtaDown[3]   = SB_Regions_allt_E3->GetBinContent(3);
	EtaUp[3]     = SB_Regions_allt_E3->GetBinContent(4);
	EtaHigh[3]   = SB_Regions_allt_E3->GetBinContent(5);
	EtaHigher[3] = SB_Regions_allt_E3->GetBinContent(6);

   // Import values for Tagger Energy Bins/Channels
   double RangeA[4];
   double RangeB[4];
   int TC_RangeA[4];
   int TC_RangeB[4];
	// Tagger channels used for photon energy binning
	TC_RangeA[0] = EChanRange->GetBinContent(1);
	TC_RangeB[0] = EChanRange->GetBinContent(2);
	// Photon energy value used for photon energy binning
	RangeA[0]    = EChanValue->GetBinContent(1);
	RangeB[0]    = EChanValue->GetBinContent(2);

	// Tagger channels used for photon energy binning
	TC_RangeA[1] = EChanRange->GetBinContent(3);
	TC_RangeB[1] = EChanRange->GetBinContent(4);
	// Photon energy value used for photon energy binning
	RangeA[1]    = EChanValue->GetBinContent(3);
	RangeB[1]    = EChanValue->GetBinContent(4);

	// Tagger channels used for photon energy binning
	TC_RangeA[2] = EChanRange->GetBinContent(5);
	TC_RangeB[2] = EChanRange->GetBinContent(6);
	// Photon energy value used for photon energy binning
	RangeA[2]    = EChanValue->GetBinContent(5);
	RangeB[2]    = EChanValue->GetBinContent(6);

	// Tagger channels used for photon energy binning
	TC_RangeA[3] = EChanRange->GetBinContent(7);
	TC_RangeB[3] = EChanRange->GetBinContent(8);
	// Photon energy value used for photon energy binning
	RangeA[3]    = EChanValue->GetBinContent(7);
	RangeB[3]    = EChanValue->GetBinContent(8);

   // Rescaling the LGD
   double EFRWD = ER * Efrwd;
   double PVECT[7][4];
   for (int ph = 0; ph < 7; ph++) {
      for (int fv = 0; fv < 4; fv++) {
         PVECT[ph][fv] = ER * pvect[ph][fv];
      }
   }

   // Setting TLorentzVectors for all photons
   TLorentzVector photon_0;
   TLorentzVector photon_1;
   photon_0.SetPxPyPzE(PVECT[0][1],PVECT[0][2],PVECT[0][3],PVECT[0][0]);
   photon_1.SetPxPyPzE(PVECT[1][1],PVECT[1][2],PVECT[1][3],PVECT[1][0]);

   //////////////////////////////////////////////////////////
   //////////// Determining the Invariant Masses ////////////
   //////////////////////////////////////////////////////////
   // 2 Gamma Invariant Mass (imass2)
   double imass2 = 0.0;
   imass2 = (photon_0 + photon_1).M();

   #if TYPECUT
   // Cut to look at cluster_cleanup turned on
   if (evtype > 10) {
      return false;
   }
   #endif

   ///////////////////////////////////////////////////////
   //////////// Performing the first few cuts ////////////
   ///////////////////////////////////////////////////////
   // Cut #1 - Number of showers reconstructed in the LGD
   if (nfrwd != fwdgamma) {
      return false;   
   }

   // Cut #2 - Number of reconstructed recoil particles
   // Exactly 1 pixel cluster in the BSD
   if (nrec != 1) {
      return false;
   }

   // Cut #3 - No showers unrelated to the recoil reconstructed in the BGD
   if (nphot != bgvgamma) {
      return false;
   }

   // Cut #4 - Individual shower energy in LGD
   // Individual shower energy in the LGD must be at least 0.05 GeV
   int cut4 = 0;
   for (int i = 0; i < nphot; i++) {
      if ((PVECT[i][0]) < 0.050) {
         cut4 += 1;
      }
   }
   if (cut4 != 0) {
      return false;
   }

   // Cut #5 - Total shower energy in the LGD
   // Total shower energy in the LGD must be at least 3.0 GeV
   // Use this if allowing a photon in BGD (e.g. nphot != nfrwd)
   // Otherwise i < nphot can be used.
   double tot_shower_e = 0.0;
   for (int i = 0; i < nfrwd; i++) {
      tot_shower_e += PVECT[i][0];
   }
   if (tot_shower_e < 3.0) {
      return false;
   }

   // Cut #6 - No coincident hits in the forward veto scintillators (CPV)
   // No coincidence between a hit in the CPV & a recoil particle
   int cpv_recoil_coin = 0;
   // Cycling through the hits in the CPV
   for (int i = 0; i < nhcpv; i++) {
      // Cycling through number of times there were a CPV hit [i]
      for (int j = 0; j < ntcpv[i]; j++) {
         // Leading edge times minus recoil proton time
         double dt = le[t1cpv[i] + j] - trec0;
         if (dt >= -3.0 && dt <= 3.0) {
            cpv_recoil_coin += 1;
         }
      }
   }
   if (cpv_recoil_coin != 0) {
      return false;
   }

   // Cut #7 - No low-energy showers around the beam hole in the LGD
   // Individual showers in the LGD that are < 0.5 GeV must satisfy
   // the following condition. This reduces the observed EM background
   // that peaks at low angles & low energies.
   int beam_hole_cut = 0;
   for (int i = 0; i < nphot; i++) {
      double lgd_xy = (sqrt((sqr(PVECT[i][1])) + (sqr(PVECT[i][2]))));
      double lgd_z = PVECT[i][3];
      // Polar angle in deg.
      double theta = ((atan2(lgd_xy, lgd_z))*(180/PI));
      // individual shower energy 
      double ishower_E = PVECT[i][0];
      double test_level = (1.033 - (0.13*theta));
      if ((ishower_E <= test_level) && (ishower_E < 0.5)) {
         beam_hole_cut += 1;
      }
   }
   if (beam_hole_cut != 0) {
      return false;
   }

   // Cut #8 - Suppresses proton resonances (cleans up the signal)
   int summed_eloss = 0;
   double E_missing = 0.0;
   for (int i = 0; i < ncoin; i++) {
      E_missing = (coenergy[i] - EFRWD);
      if (E_missing < ResonanceEner) {
         summed_eloss += 1;
      }
   }
   if (summed_eloss == 0) {
      return false;
   }

   ///////////////////// CALCULATION SECTION /////////////////////
   // Summed shower components of the LGD (lab frame)
   double e_lab  = 0.0;
   double px_lab = 0.0; 
   double py_lab = 0.0; 
   double pz_lab = 0.0; 
   for (int indv = 0; indv < nphot; indv++) {
      e_lab  += PVECT[indv][0];     
      px_lab += PVECT[indv][1];     
      py_lab += PVECT[indv][2];     
      pz_lab += PVECT[indv][3];     
   }

   // Summed LGD gammas' 4-vector in the lab frame
   TLorentzVector decay_lab;
   decay_lab.SetPxPyPzE(px_lab, py_lab, pz_lab, e_lab);

   // delta phi = azimuthal angular difference (off the decay plane)
   // between the meson of interest & the recoil proton.
   /* Delta phi, which should be zero, is used to suppress proton resonances. */
   double phi_rec = 0.0;
   double rec_phi = 0.0;
   double two_phi = 0.0;
   double del_phi = 0.0;

   phi_rec = phirec[0] * (180.0 / PI); 
   rec_phi = constrain_angle(phi_rec);   // Gives the angle - 180 deg.
   two_phi = decay_lab.Phi() * (180.0 / PI);
   del_phi = angle_diff(two_phi, rec_phi);

   // Filling t-bins
   for (int taghit = 0; taghit < ncoin; taghit++) {
      double tabs = 0.0;
      double dt = 0.0;
      double tag_wgt = 0.0;
      double Piwgt = 0.0;
      double Pidphiwgt = 0.0;
      double WSig = 0.0;
      double WSB1 = 0.0;
      double WSB2 = 0.0;
      int tagErange = 4;  // out of spec range to start, so I know if there is a coding error

      // Time of best recoil proton candidate minus tagger hit time of hit[taghit].
      dt = trec0 - cotime[taghit];

      // Setting tag weight (accidental subtraction weight) for each coincidence hit.
      if (dt >= -1.0 && dt <= 5.0) {
         tag_wgt = 1.0;
      }

      else if (dt > 5.0 && dt <= 11.0) {
         tag_wgt = -(1.0 / ASF);
      }

      else {
         tag_wgt = 0.0;
      }

      // Mandelstam variable |t|
      // tabs = -((imass2 * imass2) - (2 * coenergy[taghit] * e_lab) + (2 * coenergy[taghit] * pz_lab));
      tabs = -(sqr(imass2) - 2 * coenergy[taghit] * (e_lab - pz_lab));

      // Getting mass sideband subtraction weights weights
      /* IF THERE ARE MORE THAN 4 PHOTON ENERGY BINS, THEN YOU NEED MODIFY THIS SECTION */
		if (cochan[taghit] >= TC_RangeA[0] && cochan[taghit] <= TC_RangeB[0]) {
         WSig = w0_spline_E0->Eval(tabs);
         WSB1 = w1_spline_E0->Eval(tabs);
         WSB2 = w2_spline_E0->Eval(tabs);
         tagErange = 0;
		}
		else if (cochan[taghit] >= TC_RangeA[1] && cochan[taghit] <= TC_RangeB[1]) {
         WSig = w0_spline_E1->Eval(tabs);
         WSB1 = w1_spline_E1->Eval(tabs);
         WSB2 = w2_spline_E1->Eval(tabs);
         tagErange = 1;
		}
		else if (cochan[taghit] >= TC_RangeA[2] && cochan[taghit] <= TC_RangeB[2]) {
         WSig = w0_spline_E2->Eval(tabs);
         WSB1 = w1_spline_E2->Eval(tabs);
         WSB2 = w2_spline_E2->Eval(tabs);
         tagErange = 2;
		}
		else if (cochan[taghit] >= TC_RangeA[3] && cochan[taghit] <= TC_RangeB[3]) {
         WSig = w0_spline_E3->Eval(tabs);
         WSB1 = w1_spline_E3->Eval(tabs);
         WSB2 = w2_spline_E3->Eval(tabs);
         tagErange = 3;
		}

      // Since SB mass ranges is the same over all |t|, but not tagger energy ranges
		if (imass2 >= EtaLower[tagErange] && imass2 < EtaDown[tagErange]) {
			Piwgt = WSB1;
		}

		else if (imass2 >= EtaDown[tagErange] && imass2 <= EtaUp[tagErange]) {
			Piwgt = WSig;
		}

		else if (imass2 > EtaUp[tagErange] && imass2 <= EtaHigher[tagErange]) {
			Piwgt = WSB2;
		}

		else {
			Piwgt = 0.0;
		}

      // Filling histograms
      /* for (int xp = 0; xp < tBinNum; xp++) { */
      for (int xp = 0; xp < 3; xp++) {
         if (tabs >= Etap_bins[xp][0] && tabs <= Etap_bins[xp][1]) {
            /* IF THERE ARE MORE THAN 4 PHOTON ENERGY BINS, THEN YOU NEED MODIFY THIS SECTION */
            if (tagErange == 0) {
					tEtap0->Fill(tabs, (tag_wgt * Piwgt));
					dphi_Etap0->Fill(del_phi, (tag_wgt * Piwgt));
					Etap_dphi_E0_t[xp]->Fill(del_phi, (tag_wgt * Piwgt));
            }
				else if (tagErange == 1) {
					tEtap1->Fill(tabs, (tag_wgt * Piwgt));
					dphi_Etap1->Fill(del_phi, (tag_wgt * Piwgt));
					Etap_dphi_E1_t[xp]->Fill(del_phi, (tag_wgt * Piwgt));
				}
				else if (tagErange == 2) {
					tEtap2->Fill(tabs, (tag_wgt * Piwgt));
					dphi_Etap2->Fill(del_phi, (tag_wgt * Piwgt));
					Etap_dphi_E2_t[xp]->Fill(del_phi, (tag_wgt * Piwgt));
				}
				else if (tagErange == 3) {
					tEtap3->Fill(tabs, (tag_wgt * Piwgt));
					dphi_Etap3->Fill(del_phi, (tag_wgt * Piwgt));
					Etap_dphi_E3_t[xp]->Fill(del_phi, (tag_wgt * Piwgt));
            }
         }
      }
   }

   return kTRUE;
}

void TwoG_real2_etaprime2gp_may_16_2023::SlaveTerminate()
{
  
}

void TwoG_real2_etaprime2gp_may_16_2023::Terminate()
{
   // Write histograms to the output file
   fHistFile = new TFile("TwoG_real2_etaprime2gp_may_16_2023Analysis.root", "RECREATE");   // Declare output file
   std::cout << "pushing result histograms to a file..." << std::endl;
   std::cout << "my outputlist has " << GetOutputList()->GetEntries() << " entries." << std::endl;
   for (TIter iter = GetOutputList()->begin(); iter != GetOutputList()->end(); ++iter) {
      (*iter)->Write();
   }
   fHistFile->Close();
   std::cout << "Output file has been written" << std::endl;
}

//////////////////////////////////////////////////////

   double TwoG_real2_etaprime2gp_may_16_2023::constrain_angle(double x) {
      double anglex;
      x = remainder(x, 360.0);
      if (x < 0) {
         x += 360.0;
      }
      anglex = (x - 180.0);
      return anglex;
   }

   double TwoG_real2_etaprime2gp_may_16_2023::angle_diff(double angle1, double angle2) {
      double diff = angle1 - angle2;
      while (diff < -180.0) {
         diff += 360.0;
      }
      while (diff > 180.0) {
         diff -= 360.0;
      }
      return diff;
   }

   double TwoG_real2_etaprime2gp_may_16_2023::delta_energy(int n_coin, float tagger_energy[30], double particle_ener) {
      double del_ener = 0.0;
      for (int i = 0; i < n_coin; i++) {
         del_ener = tagger_energy[i] - particle_ener;
      }
      return del_ener;
   }
