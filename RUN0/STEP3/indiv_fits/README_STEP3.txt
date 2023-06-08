////////////////////////////
// Author: James McIntyre //
// Etap Analysis STEP #3  //
////////////////////////////
The C++ macro in STEP3 fits histograms from STEP1, prints PNGs, and creates a ROOT file with parameters
used in the subsequent analysis steps.

The sideband subtraction range limits for each photon energy bin (Echannel) are imported fro STEP2:
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


////////////////////////////////////////
// |t| binning is imported from STEP2 //
////////////////////////////////////////
The |t| range and bins are listed on line 201 (BinM2g[20][2]) of STEP2 C++ macro(or there about).
To change |t| binning for the analysis you need to change "BinM2g" and "M2g_bins" in STEP2 macro 
(used to define histo names) and "BinM2g" (used in histogram filling) in STEP1 TSelector.
Current range:
   double BinM2g[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

The above |t| bin values are stored in the STEP2 ROOT file as:
   TH1D hBinDW    // Lower bound
   TH1D hBinUP    // Upper bound

They are imported into STEP3 macro via:
      Etap_bins[tchannel][0] = hBinDW->GetBinContent(tchannel + 1);
      Etap_bins[tchannel][1] = hBinUP->GetBinContent(tchannel + 1);

////////////////////////////////////
// Changing Photon Energy Binning //
////////////////////////////////////
Number of Energy Bins (currently 4):
   TH1D EChanBINS  (from STEP2)
Imported as:
   int NumEChan = EChanBINS->GetBinContent(1);

Tagger Channels used in each Photon Energy Bin (Lower & Upper Channels are listed):
   TH1D EChanRange  (from STEP2)
Imported as:
   TC_RangeA[Echannel]
   TC_RangeB[Echannel]

Histogram with corresponding photon energies (GeV) for the range listed above("EChanRange"):
   TH1D EChanValue  (from STEP2)
Imported as:
   RangeA[Echannel]
   RangeB[Echannel]


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

/////////////////////
// Fitting Process //
/////////////////////
Invariant mass peaks for each signal change with photon energy and |t| 
(particularly at high and low |t|). The following fitting steps are really
useful for individual |t| bin fits (STEP3), but are useful here too.
1) Create and fit two TF1 single gaussian fits, one for the range of the Pi0 and 
   one for the Eta. Save the mean values as a double.
2) Create TF1 total fit (ftot) from a combination of double gaussians (Pi0, Eta, and Omega),
   single gaussian (Eta prime), and a third order polynomial with a scaling factor (pol3 * Scaler).
3) Fix fit parameters to values "somewhat" around what they should be, except pol3 & scaler set them to 1.
   Too close values can cause the fitter to fail.
4) Run "fixed parameter" fit with options "RB".
5) Release all parameters, except the Scaler. Keep that set to 1.
6) Set parameter limits (SetParLimits). Use the mean values +/- some value for the Pi0 and 
   Eta mean limits. Set widths to a fairly wide range. Keep the widths for a single double gaussian 
   the same range.
7) Fit the histogram of interest with options "SBRM".
8) Save pol3 parameters in a ROOT file, along with the other values mentioned at the beginning of this file.
9) Print fitting histograms, create directories, and save/organize PNG's in the created subdirectories.
10)Once the macro completes a visual inspection of the fit histograms saved as PNG’s is required. Look for inappropriate fit. 

/////////////////////
// Troubleshooting //
/////////////////////
If the fitting completes, but gives an odd fit:
1) Write an IF statement for that particular histogram (i.e. Echannel == 1 && tchannel == 7) and alter the fit limits.

If the fits fails try the following (in the order listed):
1) Remove "M" from the options in step (7) above, then try to run the macro again.
   The fitter sometimes gets stuck minimizing the fit with this option.
2) Try changing the initial fixed fit parameter values that are set in step 3 above.
   Values too close to the correct value can cause the fitter to fail. Use an IF Statement
   to change the initial values for the failed histogram only. Sometimes there is only one 
   histogram that is the problem child. This step usually fixes those fittings that fail 
   even though the fit looks really good (e.g. initial parameters are too close and 
   the fitter gets stuck while trying to minimize the fit).
3) If the previous step fails, add an IF statement and implement different parameter limits
   in step (6) above. Base your changes off of the output text file created from your
   previous fit (2g_E_binned_Etaprime_fit_allt_*.txt). The IF statement is added to change the
   parameter limits for the failed histogram and those there after. Limit the number of elements
   you change at once, several attempts may be required.
4) If all of the above fails to fix things, have a look at the failed TH1D. The single gaussian fit
   for Pi0 or Eta might give a mean that's way off (i.e. from a signal too weak or too broad).
   Also, check the text file "2g_E_binned_Etaprime_fit_allt_*.txt" for the actual fitted values.
   If these means are way off then your SetParLimit values will be out of range of the signal.

////////////
// OUTPUT //
////////////
STEP3's macro produces a root file containing the sideband weights and their errors for analysis use in further steps:
TFile:
   TwoG_may_16_2023_STEP3.root

Histograms containing analysis information:
	TH1D h_W0_E3  -> Signal Region Weight
	TH1D h_W0_E3  -> Signal Region Error

	TH1D h_W1_E3  -> Left (lower mass) Sideband Region Weight
	TH1D h_W1_E3  -> Left (lower mass) Sideband Region Error

	TH1D h_W2_E3  -> Right (higher mass) Sideband Region Weight
	TH1D h_W2_E3  -> Right (higher mass) Sideband Region Error




