///////////////////////////////////////////////////////////
// Author: J. McIntyre
// 16 May 2023
///////////////////////////////////////////////////////////
// This file is used to run PROOF on Radphi dataset
///////////////////////////////////////////////////////////
//
//    To change this file modify the following:
//    -----------------------------------------
//    PROOF Service Used --> 1097 
//                           /* Jim's 1099    */
//                           /* Fridah's 1097 */
//    Radphi Dataset     --> pass-9-2020
//    This files name    --> proof_TwoG_may_16_2023
//    TSelector Name     --> TwoG_may_16_2023
//
//   ///////////////// For MC Data ///////////////////
//   MC reaction(s) of interest --> pi02gp_clean
//   MC reaction(s) of interest --> eta2gp_clean
//   MC reaction(s) of interest --> etaprime2gp_clean
//
//        "If more are needed update lines 69 & 76"
//
//            "Uncomment lines 59, 63, 74, 81"
//           "Comment out lines 58, 62, 73, 80"
//
//   //////////////// For Real Data //////////////////
//            "Uncomment lines 58, 62, 73, 80"
//           "Comment out lines 59, 63, 74, 81"
//
//
// --------------------------------------------------------
// Associated real data .root file
// --> TwoG_may_16_2023Analysis.root
// --------------------------------------------------------

#include <sys/stat.h>
#include <unistd.h>
#include <iostream>

#include <TROOT.h>
#include <TProof.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TMacro.h>
#include <TH2D.h>
#include <fstream>
#include <string>
#define DOING_PROOF 1
#define ASYNC_PROOF 1


TString proofserver("proof://stat31.phys.uconn.edu:1097/?N");

/* Can use nod 25, 26, 27, 28, or 29 for xrootdpath */
TString xrootdpath("root://nod26.phys.uconn.edu/Gluex/radphi/pass-9-2020");    // Real data
//TString xrootdpath("root://nod26.phys.uconn.edu/Gluex/radphi/sims_2020");    // MC data
//TString xrootdpath("/pnfs/phys.uconn.edu/data/Gluex/radphi/sims_2020");      // Alternate path

TString pnfspath("root://nod26.phys.uconn.edu/Gluex/radphi/pass-9-2020");      // Real data
//TString pnfspath("root://nod26.phys.uconn.edu/Gluex/radphi/sims_2020");      // MC data
//TString pnfspath("/pnfs/phys.uconn.edu/data/Gluex/radphi/pass-9-2020");      // Alternate path

TChain *fullchain()
{
   TChain *ch = new TChain("h1");
   int filecount = 0;
   void *datadir = gSystem->OpenDirectory(pnfspath);
   while (const char *fname = gSystem->GetDirEntry(datadir)) {
      std::string fnames(fname);
      if (fnames.find(".root") != fnames.npos) {  // For real data
      //if ((fnames.find("etaprime2gp_gen") != fnames.npos) && (fnames.find(".root") != fnames.npos)) {  // For MC data
         ch->Add(xrootdpath + "/" + TString(fnames));
         filecount += 1;
      }
   }
   std::cout << "\n";
   std::cout << "There are ==> " << filecount << " files being processed from: " << xrootdpath << "\n" ;
   //std::cout << "There are ==> " << filecount << " Monte Carlo files with name(s) - " << "etaprime2gp_gen" << " - being processed. \n" ;
   gSystem->FreeDirectory(datadir);
   return ch;
}


void proof_TwoG_may_16_2023(int nevents=2000000000)
{
   TChain *ch=fullchain();

#if DOING_PROOF
   if (gProof == 0) {
      TProof::Open(proofserver);
      gProof->SetParameter("PROOF_LookupOpt", "none");
      gProof->SetParameter("PROOF_MaxSlavesPerNode", (Int_t)108);
   }

   gProof->Load("TwoG_may_16_2023.C+O");
   ch->SetProof();
#else
   gROOT->ProcessLine(".L TwoG_may_16_2023.C+O");
#endif

   ch->Process("TwoG_may_16_2023", "", nevents);

#if DOING_PROOF
   if (gProof == 0) {
      gProof->Close();
   }
#endif
}
