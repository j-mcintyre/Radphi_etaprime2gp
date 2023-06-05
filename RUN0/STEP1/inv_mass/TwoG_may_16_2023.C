#define TwoG_may_16_2023_cxx

// Author: J. McIntyre
// May 16, 2023
// This file contains cuts on real data for a 2 photon final state sample.
//////////////////////////////////////////////////////////////////////////
// This was written to look at the 2 gamma sample.                      //
// ____________________________________________________________________ //
// ____________________________________________________________________ //
//                                                                      //
// THIS FILE HAS UPDATES THAT:                                          //
// --------------------------                                           //
//    - USES THE LATEST ENERGY RESCALING FACTOR (ER) = 0.97             //
// ____________________________________________________________________ //
//                                                                      //
// This file is a re-do of TwoG_jul_16_2021.C, which was a re-do of     //
// ER_real_pi02gp_oct_21_2021.C                                         //
// /nfs/direct/annex/mcintyre/analysisIII/ER/ER_0_955/5MeV_mass_bins    //
// since the ER & ASF were updated and a valid cross section for pi0    //
// was plotted. Now this file is the start of a final analysis of the   //
// two gamma sample for the pi0, eta, and etap cross sections.          //
//                                                                      //
// The new ER equals 0.97, as determined by the various attempts in     //
// folder: /nfs/direct/annex/mcintyre/analysisIX/TwoGamma/inv_mass      //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

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

#include "TwoG_may_16_2023.h"
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <math.h>
#include <cmath>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <TProofServ.h>
#include <TF1.h>
#include <vector>
#include <TH1D.h>
#include <TH1.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>

/*
#include <TChain.h>
#include <Math/MinimizerOptions.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TFile.h>
#include <fstream>
#include <string>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TLatex.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TLegend.h>
#include <sstream>
#include <TApplication.h>
#include <TRint.h>
#include <TDirectory.h>
#include <TAxis.h>
#include <TColor.h>
#include <TPad.h>
#include <TKey.h>
#include <TText.h>
#include <TMarker.h>
#include <TF2.h>
#include "TObjString.h"
#include "TRandom3.h"
#include "TGraphSmooth.h"
#include "TSpline.h"
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

#define TYPECUT 1            // Cluster_cleanup cut (1 = On, 0 = Off)

#define THESISPLOTS 1        // Thesis plots for dt & ASF

TFile *fHistFile;            // Set Pointer to the output file


void TwoG_may_16_2023::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}

void TwoG_may_16_2023::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
   // Histogram definition

   // No cuts, No cluster cleanup algorithm
   hcut = new TH1D("hcut", "M_{2#gamma}: no cleanup cuts, cluster cleanup algorithm not applied;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut->Sumw2();
   GetOutputList()->Add(hcut);

   // No cuts, Cluster cleanup algorithm applied
   hcut0 = new TH1D("hcut0", "M_{2#gamma}: no cleanup cuts, cluster cleanup algorithm applied;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut0->Sumw2();
   GetOutputList()->Add(hcut0);

   // Cut 1
   hcut1 = new TH1D("hcut1", "Two showers reconstructed in the LGD;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut1->Sumw2();
   GetOutputList()->Add(hcut1);

   // Cut 2
   hcut2 = new TH1D("hcut2", "Exactly 1 pixel cluster in the BSD;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut2->Sumw2();
   GetOutputList()->Add(hcut2);

   // Cut 3
   hcut3 = new TH1D("hcut3", "No showers unrelated to the recoil p reconstructed in the BGD;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut3->Sumw2();
   GetOutputList()->Add(hcut3);

   // Cut 4
   hcut4 = new TH1D("hcut4", "Individual shower energy in the LGD must be at least 0.05 GeV;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut4->Sumw2();
   GetOutputList()->Add(hcut4);

   // Cut 5
   hcut5 = new TH1D("hcut5", "Total shower energy in the LGD must be at least 3.0 GeV;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut5->Sumw2();
   GetOutputList()->Add(hcut5);

   // Cut 6
   hcut6 = new TH1D("hcut6", "No coincidence between a hit in the CPV & a recoil particle;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut6->Sumw2();
   GetOutputList()->Add(hcut6);

   // Cut 7
   hcut7 = new TH1D("hcut7", "No low-energy showers around the beam hole in the LGD;M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut7->Sumw2();
   GetOutputList()->Add(hcut7);

   // Cut 8
   hcut8 = new TH1D("hcut8", "Missing Energy Cut (suppresses proton resonances);M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   hcut8->Sumw2();
   GetOutputList()->Add(hcut8);

   // Coincidence hits after cuts 1 - 8
   h_ncoin = new TH1D("h_ncoin", "Number of tagger hits per event;Coincidence Hits;Number of 2#gamma Events", 124, -1.0, 30.0);
   h_ncoin->Sumw2();
   GetOutputList()->Add(h_ncoin);

   // Time difference between recoil proton and the i'th tagger hit
   // (Time of best recoil proton candidate minus tagger hit time of hit[i])
   h_dt = new TH1D("h_dt", "Time difference between recoil proton and tagger hit;t_{recoil particle} - t_{tagger hit} [ns];Events per 0.1 ns", 1900, -90.0, 100.0);
   h_dt->Sumw2();
   GetOutputList()->Add(h_dt);

   // Cuts 1-8, channel diff. > 2, Peak at -40 ns, tag_{j} btw -53 to -27 ns;cotime_{i} - cotime_{j} [ns]
   h2g_all = new TH1D("h2g_all", "#Deltat_{tag} for peak at -40 ns and range -53 to -27 ns;t_{accidental tag (i)} - t_{accidental tag (j)} [ns];Events per 0.1 ns", 320, -16.0, 16.0);
   h2g_all->Sumw2();
   GetOutputList()->Add(h2g_all);

   // t_{tag_{i}} - t_{tag_{j}}, Cuts 1-8, channel diff. > 2, Peak at -42 ns, tag_{j} btw -61 to -25 ns
   h_dcotime = new TH1D("h_dcotime", "#Deltat_{tag} for peak at -42 ns and range -61 to -25 ns;t_{accidental tag (i)} - t_{accidental tag (j)} [ns];Events per 0.1 ns", 420, -22.0, 20.0);
   h_dcotime->Sumw2();
   GetOutputList()->Add(h_dcotime);

   // t_{tag_{i}} - t_{tag_{j}}, Cuts 1-8, channel diff. > 2, Peak at -44 ns, tag_{j} btw -61 to -25 ns
   h_dcotime2 = new TH1D("h_dcotime2", "#Deltat_{tag} for peak at -44 ns and range -61 to -25 ns;t_{accidental tag (i)} - t_{accidental tag (j)} [ns];Events per 0.1 ns", 420, -20.0, 22.0);
   h_dcotime2->Sumw2();
   GetOutputList()->Add(h_dcotime2);

   // |t| for 2#gamma multiplicity
   t_2gamma = new TH1D("t_2gamma", "Mandelstam Variable |t|;|t| [GeV^{2}];Events per 5 MeV^{2}", 1000, 0.0, 5.0);
   t_2gamma->Sumw2();
   GetOutputList()->Add(t_2gamma);

   // M_{2#gamma} for 2#gamma final state (0.1 to 1.9 GeV^2)
   M_2gamma = new TH1D("M_2gamma", "M_{2#gamma} for 0.1 #leq |t| #geq 1.9 GeV^{2};M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", 280, 0.0, 1.4);
   M_2gamma->Sumw2();
   GetOutputList()->Add(M_2gamma);

   // Tagger Channels
   TChan = new TH1D("TChan", "Tagger Channel Events for a 2#gamma Final State;Tagger Channel No.;Events per channel", 84, -1.0, 20.0);
   TChan->Sumw2();
   GetOutputList()->Add(TChan);

   ////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////
   // CHANGE |t| BINS HERE AND AT THE BEGINNING OF Process IN THIS TSELECTOR //
   // ALL SUBSEQUENT ANALYSIS TSELECTORS WILL REFERENCE THE HISTOGRAMS MADE  //
   // IN STEP2 FOR WHAT |t| BINS TO USE.                                     //
   ////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////
   double M2g_bins[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

   //////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////
   // BELOW IS WHERE YOU CHANGE THE TAGGER ENERGY BINS USED IN ALL FURTHER ANALYSIS//
   //////////////////////////////////////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////
   // High energy end of bin in GeV
   // Tagchan# energy value is border between tagger channel # & (# - 1)
   double Tagchan0 = 5.385;
   double Tagchan1 = 5.335;
   double Tagchan2 = 5.280;
   double Tagchan3 = 5.220;
   double Tagchan4 = 5.165;
   double Tagchan5 = 5.115;
   double Tagchan6 = 5.065;
   double Tagchan7 = 5.015;
   double Tagchan8 = 4.965;
   double Tagchan9 = 4.915;
   double Tagchan10 = 4.860;
   double Tagchan11 = 4.800;
   double Tagchan12 = 4.740;
   double Tagchan13 = 4.675;
   double Tagchan14 = 4.615;
   double Tagchan15 = 4.570;
   double Tagchan16 = 4.530;
   double Tagchan17 = 4.485;
   double Tagchan18 = 4.435;
   double Tagchan19 = 4.385;

   // Tagger Energy Bin E0 = channel 0 thru channel 4
   double Range0a = Tagchan0;
   double Range0b = Tagchan5;

   // Tagger Energy Bin E1 = channel 5 thru channel 9
   double Range1a = Tagchan5;
   double Range1b = Tagchan10;

   // Tagger Energy Bin E2 = channel 10 thru channel 13
   double Range2a = Tagchan10;
   double Range2b = Tagchan14;

   // Tagger Energy Bin E0 = channel 14 thru channel 18
   double Range3a = Tagchan14;
   double Range3b = Tagchan19;

   //////////////////////////////////////////////////////////////////////////////////////
   // CHANGE ENERGY RANGE IN TITLE FOR TH1D BASED ON NUMBER OF TAGGER ENERGY BINS USED //
   //////////////////////////////////////////////////////////////////////////////////////
   TString title;
   TString name;

   for (int i = 0; i < 3; i++) {
      name.Form("M2g_E0_%d", i);
      title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for %.3f #leq |t| < %.3f;M_{2#gamma} [MeV/c^{2}];Events per 5 MeV/c^{2}", Range0b, Range0a, M2g_bins[i][0], M2g_bins[i][1]);
      M2g_E0_[i] = new TH1D(name, title, 280, 0.0, 1.4);
      M2g_E0_[i]->Sumw2();
      GetOutputList()->Add(M2g_E0_[i]);
   }

   //////////////////////////////////////////////////////////////////////////////////////
   // CHANGE ENERGY RANGE IN TITLE FOR TH1D BASED ON NUMBER OF TAGGER ENERGY BINS USED //
   //////////////////////////////////////////////////////////////////////////////////////
   for (int i = 0; i < 3; i++) {
      name.Form("M2g_E1_%d", i);
      title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for %.3f #leq |t| < %.3f;M_{2#gamma} [MeV/c^{2}];Events per 5 MeV/c^{2}", Range1b, Range1a, M2g_bins[i][0], M2g_bins[i][1]);
      M2g_E1_[i] = new TH1D(name, title, 280, 0.0, 1.4);
      M2g_E1_[i]->Sumw2();
      GetOutputList()->Add(M2g_E1_[i]);
   }

   //////////////////////////////////////////////////////////////////////////////////////
   // CHANGE ENERGY RANGE IN TITLE FOR TH1D BASED ON NUMBER OF TAGGER ENERGY BINS USED //
   //////////////////////////////////////////////////////////////////////////////////////
   for (int i = 0; i < 3; i++) {
      name.Form("M2g_E2_%d", i);
      title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for %.3f #leq |t| < %.3f;M_{2#gamma} [MeV/c^{2}];Events per 5 MeV/c^{2}", Range2b, Range2a, M2g_bins[i][0], M2g_bins[i][1]);
      M2g_E2_[i] = new TH1D(name, title, 280, 0.0, 1.4);
      M2g_E2_[i]->Sumw2();
      GetOutputList()->Add(M2g_E2_[i]);
   }

   //////////////////////////////////////////////////////////////////////////////////////
   // CHANGE ENERGY RANGE IN TITLE FOR TH1D BASED ON NUMBER OF TAGGER ENERGY BINS USED //
   //////////////////////////////////////////////////////////////////////////////////////
   for (int i = 0; i < 3; i++) {
      name.Form("M2g_E3_%d", i);
      title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for %.3f #leq |t| < %.3f;M_{2#gamma} [MeV/c^{2}];Events per 5 MeV/c^{2}", Range3b, Range3a, M2g_bins[i][0], M2g_bins[i][1]);
      M2g_E3_[i] = new TH1D(name, title, 280, 0.0, 1.4);
      M2g_E3_[i]->Sumw2();
      GetOutputList()->Add(M2g_E3_[i]);
   }

   // M_{2#gamma} for 2#gamma final state (0.1 to 1.05 GeV^2)
   title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for 0.1 #leq |t| #geq 1.9 GeV^{2};M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", Range0b, Range0a);
   Allt_E0 = new TH1D("Allt_E0", title, 280, 0.0, 1.4);
   Allt_E0->Sumw2();
   GetOutputList()->Add(Allt_E0);

   // M_{2#gamma} for 2#gamma final state (0.1 to 1.05 GeV^2)
   title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for 0.1 #leq |t| #geq 1.9 GeV^{2};M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", Range1b, Range1a);
   Allt_E1 = new TH1D("Allt_E1", title, 280, 0.0, 1.4);
   Allt_E1->Sumw2();
   GetOutputList()->Add(Allt_E1);

   // M_{2#gamma} for 2#gamma final state (0.1 to 1.05 GeV^2)
   title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for 0.1 #leq |t| #geq 1.9 GeV^{2};M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", Range2b, Range2a);

   Allt_E2 = new TH1D("Allt_E2", title, 280, 0.0, 1.4);
   Allt_E2->Sumw2();
   GetOutputList()->Add(Allt_E2);

   // M_{2#gamma} for 2#gamma final state (0.1 to 1.05 GeV^2)
   title.Form("M_{2#gamma} from E_{#gamma} = %.3f to %.3f GeV for 0.1 #leq |t| #geq 1.9 GeV^{2};M_{2#gamma} [GeV/c^{2}];Events per 5 MeV/c^{2}", Range3b, Range3a);
   Allt_E3 = new TH1D("Allt_E3", title, 280, 0.0, 1.4);
   Allt_E3->Sumw2();
   GetOutputList()->Add(Allt_E3);
}


Bool_t TwoG_may_16_2023::Process(Long64_t entry)
{
   if (fChain == 0) {   // REMEMBER:
      return false;     // false is valid for ROOT 6 and above.
   }                    // kFALSE works for all ROOT, but not C++.

   GetEntry(entry);

   // Rescaling the LGD
   double EFRWD = ER * Efrwd;
   double PVECT[7][4];
   for (int ph = 0; ph < 7; ph++) {          // I like to Postfix Increment, but for the "for loop" it doesn't matter
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

   // Histogram prior to cluster_cleanup
   hcut->Fill(imass2);

   // Cut to look at cluster_cleanup turned on
   #if TYPECUT
   if (evtype > 10) {
      return false;
   }
   #endif

   // Histogram after cluster_cleanup
   hcut0->Fill(imass2);

   ///////////////////////////////////////////////////////
   //////////// Performing the first few cuts ////////////
   ///////////////////////////////////////////////////////
   // Cut #1 - Number of showers reconstructed in the LGD
   if (nfrwd != fwdgamma) {
      return false;   
   }
   hcut1->Fill(imass2);

   // Cut #2 - Number of reconstructed recoil particles
   // Exactly 1 pixel cluster in the BSD
   if (nrec != 1) {
      return false;
   }
   hcut2->Fill(imass2);

   // Cut #3 - No showers unrelated to the recoil reconstructed in the BGD
   if (nphot != bgvgamma) {
      return false;
   }
   hcut3->Fill(imass2);

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
   hcut4->Fill(imass2);

   // Cut #5 - Total shower energy in the LGD
   // Total shower energy in the LGD must be at least 3.0 GeV
   // Use this code if allowing a photon in BGD (e.g. nphot != nfrwd)
   // Otherwise i < nphot can be used.
   double tot_shower_e = 0.0;
   for (int i = 0; i < nfrwd; i++) {
      tot_shower_e += PVECT[i][0];
   }
   if (tot_shower_e < 3.0) {
      return false;
   }
   hcut5->Fill(imass2);

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
   hcut6->Fill(imass2);

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
   hcut7->Fill(imass2);

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
   hcut8->Fill(imass2);

   // Histos showing the number of coincidence hits we are dealing with
   double ncoinplus = ncoin + 0.1;
   double ncoinminus = ncoin - 0.1;
   h_ncoin->Fill(ncoinplus);
   h_ncoin->Fill(ncoinminus);

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

   #if THESISPLOTS
   // Producing histogram used to determine ASF
   double dt_tag2 = 0.0;
   double dt2 = 0.0;
   double dt3 = 0.0;
   int delta_cochan2 = 0;

   for (int j = 0; j < ncoin; j++) {
      // reset variables to zero for reuse without bias
      dt_tag2 = 0.0;
      dt2 = 0.0;
      dt3 = 0.0;
      delta_cochan2 = 0;

      // Time difference between recoil proton and the j'th tagger hit
      // Time of best recoil proton candidate minus tagger hit time of hit[j].
      dt2 = trec0 - cotime[j];

      // Plot for thesis
      h_dt->Fill(dt2);
      
      if ((dt2 >= -41.0) && (dt2 <= -39.0)) {
         for (int k = 0; k < ncoin; k++) {
            delta_cochan2 = abs(cochan[j] - cochan[k]);
            if ((j != k) && (delta_cochan2 >= 3)) {
               dt_tag2 = cotime[j] - cotime[k];
               dt3 = trec0 - cotime[k];
               if ((dt3 >= -53.0) && (dt3 <= -27.0)) {
                  // Histogram used for determining ASF
                  h2g_all->Fill(dt_tag2);
               }
            }
         }
      }
   }

   // Alternate timing range for histogram used to determine ASF
   for (int l = 0; l < ncoin; l++) {
		// reset variables to zero for reuse without bias
		dt_tag2 = 0.0;
		dt2 = 0.0;
		dt3 = 0.0;
		delta_cochan2 = 0;

      // Time difference between recoil proton and the l'th tagger hit
      // Time of best recoil proton candidate minus tagger hit time of hit[l].
      dt2 = trec0 - cotime[l];

      if ((dt2 >= -43.0) && (dt2 <= -41.0)) {
         for (int m = 0; m < ncoin; m++) {
            delta_cochan2 = abs(cochan[l] - cochan[m]);
            if ((l != m) && (delta_cochan2 >= 3)) {
               dt_tag2 = cotime[l] - cotime[m];
               dt3 = trec0 - cotime[m];
               if ((dt3 >= -61.0) && (dt3 <= -25.0)) {
                  h_dcotime->Fill(dt_tag2);
               }
            }
         }
      }

		// reset variables to zero for reuse without bias
		dt_tag2 = 0.0;
		dt3 = 0.0;
		delta_cochan2 = 0;

      if ((dt2 >= -45.0) && (dt2 <= -43.0)) {
         for (int n = 0; n < ncoin; n++) {
            delta_cochan2 = abs(cochan[l] - cochan[n]);
            if ((l != n) && (delta_cochan2 >= 3)) {
               dt_tag2 = cotime[l] - cotime[n];
               dt3 = trec0 - cotime[n];
               if ((dt3 >= -61.0) && (dt3 <= -25.0)) {
                  h_dcotime2->Fill(dt_tag2);
               }
            }
         }
      }
   }
   #endif

   ////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   // CHANGE |t| BINS HERE AND SlaveBegin IN THIS TSELECTOR                      //
   // ALL SUBSEQUENT ANALYSIS TSELECTORS WILL REFERENCE THE HISTOGRAMS MADE      //
   // IN STEP2 FOR WHAT |t| BINS TO USE.                                         //
   ////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////
   // ALL FURTHER ANALYSIS WILL REFERENCE THE INFO HERE VIA HISTOS CREATED BELOW //
   // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv //

   // Setting |t| bins
   double BinM2g[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

   // High energy end of bin in GeV
   // Tagchan# energy value is border between tagger channel # & (# - 1)
   double Tagchan0 = 5.385;
   double Tagchan1 = 5.335;
   double Tagchan2 = 5.280;
   double Tagchan3 = 5.220;
   double Tagchan4 = 5.165;
   double Tagchan5 = 5.115;
   double Tagchan6 = 5.065;
   double Tagchan7 = 5.015;
   double Tagchan8 = 4.965;
   double Tagchan9 = 4.915;
   double Tagchan10 = 4.860;
   double Tagchan11 = 4.800;
   double Tagchan12 = 4.740;
   double Tagchan13 = 4.675;
   double Tagchan14 = 4.615;
   double Tagchan15 = 4.570;
   double Tagchan16 = 4.530;
   double Tagchan17 = 4.485;
   double Tagchan18 = 4.435;
   double Tagchan19 = 4.385;

   // Tagger Energy Bin E0 = channel 0 - channel 4
   double Range0a = Tagchan0;
   double Range0b = Tagchan5;
   int TC_Range0a = 0;
   int TC_Range0b = 4;

   // Tagger Energy Bin E1 = channel 5 - channel 9
   double Range1a = Tagchan5;
   double Range1b = Tagchan10;
   int TC_Range1a = 5;
   int TC_Range1b = 9;

   // Tagger Energy Bin E2 = channel 10 - channel 13
   double Range2a = Tagchan10;
   double Range2b = Tagchan14;
   int TC_Range2a = 10;
   int TC_Range2b = 13;

   // Tagger Energy Bin E0 = channel 14 - channel 18
   double Range3a = Tagchan14;
   double Range3b = Tagchan19;
   int TC_Range3a = 14;
   int TC_Range3b = 19; // One channel beyond tagger range

   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ //
   //      ALL FURTHER ANALYSIS WILL REFERENCE THE HISTOS CREATED BELOW      //
   ////////////////////////////////////////////////////////////////////////////

   for (int taghit = 0; taghit < ncoin; taghit++) {
      double tabs = 0.0;
      double dt = 0.0;
      double tag_wgt = 0.0;

      // Tagger channel histogram
      TChan->Fill(cochan[taghit]);

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
      // |t| histogram
      t_2gamma->Fill(tabs, tag_wgt);

      // Histograms used for fitting a background polynomial 
      // which will then be fixed (and just scaled) in the individual |t| bin fits 
      if (tabs >= BinM2g[0][0] && tabs <= BinM2g[2][1]) {
         M_2gamma->Fill(imass2, tag_wgt);
         if (cochan[taghit] >= TC_Range0a && cochan[taghit] <= TC_Range0b) {
            Allt_E0->Fill(imass2, tag_wgt);
         }
         else if (cochan[taghit] >= TC_Range1a && cochan[taghit] <= TC_Range1b) {
            Allt_E1->Fill(imass2, tag_wgt);
         }
         else if (cochan[taghit] >= TC_Range2a && cochan[taghit] <= TC_Range2b) {
            Allt_E2->Fill(imass2, tag_wgt);
         }
         else if (cochan[taghit] >= TC_Range3a && cochan[taghit] <= TC_Range3b) {
            Allt_E3->Fill(imass2, tag_wgt);
         }
      }

      // Histos used for individual |t| bin fitting of signal and 
      // background for mass sideband subtraction weight calculation
      for (int cycE = 0; cycE < 3; cycE++) {
         if (tabs >= BinM2g[cycE][0] && tabs < BinM2g[cycE][1]) {
            if (cochan[taghit] >= TC_Range0a && cochan[taghit] <= TC_Range0b) {
               M2g_E0_[cycE]->Fill(imass2, tag_wgt);
            }
            else if (cochan[taghit] >= TC_Range1a && cochan[taghit] <= TC_Range1b) {
               M2g_E1_[cycE]->Fill(imass2, tag_wgt);
            }
            else if (cochan[taghit] >= TC_Range2a && cochan[taghit] <= TC_Range2b) {
               M2g_E2_[cycE]->Fill(imass2, tag_wgt);
            }
            else if (cochan[taghit] >= TC_Range3a && cochan[taghit] <= TC_Range3b) {
               M2g_E3_[cycE]->Fill(imass2, tag_wgt);
            }
         }
      }
   }
   return kTRUE;
}

void TwoG_may_16_2023::SlaveTerminate()
{
  
}

void TwoG_may_16_2023::Terminate()
{
   // Write histograms to the output file
   fHistFile = new TFile("TwoG_may_16_2023Analysis.root", "RECREATE");   // Declare output file
   std::cout << "pushing result histograms to a file..." << std::endl;
   std::cout << "my outputlist has " << GetOutputList()->GetEntries() << " entries." << std::endl;
   for (TIter iter = GetOutputList()->begin(); iter != GetOutputList()->end(); ++iter) {
      (*iter)->Write();
   }
   fHistFile->Close();
   std::cout << "Output file has been written" << std::endl;
}

//////////////////////////////////////////////////////

   double TwoG_may_16_2023::constrain_angle(double x) {
      double anglex;
      x = remainder(x, 360.0);
      if (x < 0) {
         x += 360.0;
      }
      anglex = (x - 180.0);
      return anglex;
   }

   double TwoG_may_16_2023::angle_diff(double angle1, double angle2) {
      double diff = angle1 - angle2;
      while (diff < -180.0) {
         diff += 360.0;
      }
      while (diff > 180.0) {
         diff -= 360.0;
      }
      return diff;
   }

   double TwoG_may_16_2023::delta_energy(int n_coin, float tagger_energy[30], double particle_ener) {
      double del_ener = 0.0;
      for (int i = 0; i < n_coin; i++) {
         del_ener = tagger_energy[i] - particle_ener;
      }
   return del_ener;
   }