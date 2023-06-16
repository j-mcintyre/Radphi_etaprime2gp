////////////////////////////
// Author: James McIntyre //
// Etap Analysis STEP #2  //
////////////////////////////
The C++ macro in STEP2 fits histograms from STEP1, prints PNGs, and creates a ROOT file with parameters
used in the subsequent analysis steps. The fits in this macro are used to determine the background used
for each photon energy bin (e.g. over all |t| for an E_gamma bin). They are also used to determine the 
mass sideband subtraction range to use for each photon energy bin when fitting individual |t| binned 
histograms in STEP3.


***********************************************************************************************
*********** This fitting macro sets the minimum Eta prime width limit at 45 MeV/c^2 ***********
***********************************************************************************************


The sideband subtraction range limits for each photon energy bin (Echannel) are listed in a TH1D in the STEP2 root file:
      EtapLower[Echannel]  = Left sideband lower mass limit
      EtapLow[Echannel]    = Left sideband upper mass limit
      Etapdown[Echannel]   = Signal region lower mass limit
      Etapup[Echannel]     = Signal region upper mass limit
      EtapHigh[Echannel]   = Right sideband lower mass limit
      EtapHigher[Echannel] = Right sideband upper mass limit


The fitted histograms (listed below) cover mandelstam |t| from 0.1 - 1.9 GeV^2 and are binned in
four photon energy ranges (E0 - E3) and 20 Mandelstam |t| bins (80 histograms total):
   htitle.Form("M2g_E%d_%d", Echannel, tchannel)
   TH1D M2g_E0_0
   TH1D M2g_E0_1
   TH1D M2g_E0_2
        ...
   TH1D M2g_E1_0
   TH1D M2g_E1_1
   TH1D M2g_E1_2
        ...


/////////////////
// |t| binning //
/////////////////
The |t| range and bins are listed on line 201 (BinM2g[20][2]) of STEP2 C++ macro(or there about).
To change |t| binning for the analysis you need to change "BinM2g" and "M2g_bins" in STEP2 macro 
(used to define histo names) and "BinM2g" (used in histogram filling) in STEP1 TSelector.
Current range:
   double BinM2g[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

The above |t| bin values are stored in the STEP2 ROOT file as:
   TH1D hBinDW    // Lower bound
   TH1D hBinUP    // Upper bound

They will be imported into STEP3 macro via:
      Etap_bins[tchannel][0] = hBinDW->GetBinContent(tchannel + 1);
      Etap_bins[tchannel][1] = hBinUP->GetBinContent(tchannel + 1);

////////////////////////////////////
// Changing Photon Energy Binning //
////////////////////////////////////
Number of Energy Bins (currently 4):
   TH1D EChanBINS  (STEP2)

Tagger Channels used in each Photon Energy Bin (Lower & Upper Channels are listed):
   TH1D EChanRange  (STEP2)

Histogram with corresponding photon energies (GeV) for the range listed above("EChanRange"):
   TH1D EChanValue  (STEP2)


For reference, the Hall B tagger channel energies (GeV) used are:
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