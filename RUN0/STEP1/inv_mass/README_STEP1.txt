////////////////////////////
// Author: James McIntyre //
// Etap Analysis STEP #1  //
////////////////////////////
The TSelector of STEP1 applies various cuts to the 2 gamma dataset (2578 ROOT files), using PROOF,
to reduce background. Then the main histograms of interest are produced and saved to a ROOT file
with the name <TSelector name>Analysis.root. Additional histograms are also saved to this ROOT file,
as mentioned below, for use in future STEPs or for plots to be placed in my thesis. Subsequent STEPs
will further reduce the background using sideband subtraction and then determine the Etap cross-section.

|t| Range
   double BinM2g[3][2] = {{0.100, 0.300}, {0.300, 0.650}, {0.650, 1.050}};

Histograms that contain the invariant mass spectrum at different stages of the cleanup cut process.
The histogram titles explain them in detail.
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

Some histograms for use in my thesis. If interested, the histogram titles describe it well without 
looking at the code itself.
   TH1D *h_ncoin;
   TH1D *h_dt;
   TH1D *h2g_all;
   TH1D *h_dcotime;
   TH1D *h_dcotime2;
   TH1D *t_2gamma;
   TH1D *M_2gamma;
   TH1D *M_2gtEta;
   TH1D *M_2gtEtap;
   TH1D *TChan;

These histograms are of two gamma invariant mass binned in photon energy (E#, e.g. tagger channel ranges)
over the full range of |t| (e.g. tbin[first][0] to tbin[last][1]). These histograms are fitted to 
determine the background shape (polynomial). These histograms are used before the ones below (individual
|t| bins) so that the background polynomial background parameters can be fixed and a scaling factor is added.
The scaling factor is used to determine the background fit in the individual |t| bin fits.
  NOTE: Radphi background is treated as:
        (1) A third order polynomial for each photon energy bin (which only scales wrt |t|) 
           AND
        (2) A double gaussian around 50 MeV below the Omega mass (Omega signal leakage from 3 gamma final state)
        The Omega leakage is treated like a signal (double gaussian fit) and allowed to vary slightly
        in width and walk in mean. These are restricted since omega is so close to the eta signal and the fitting
        program wants to combine the two at higher |t|. Essentially, as |t| increases the omega leakage decreases
        and its mass lowers. It looks like it's trying to hide under the eta. For this reason we must limit the 
        omega fit's walk in mass and question the eta fitting at higher |t|.

   TH1D *Allt_E0;
   TH1D *Allt_E1;
   TH1D *Allt_E2;
   TH1D *Allt_E3;

Histograms for etaprime meson. These contain tag weighted (accidental subtracted) plots of
two gamma invariant mass for each tagger range (E#, e.g. photon energy) in individual |t| bins. These are
used to fit signal & background for mass sideband subtraction calculations (used to get fill weights for 
subsequent TSelector use).
   TH1D *M2g_E0_[20];
   TH1D *M2g_E1_[20];
   TH1D *M2g_E2_[20];
   TH1D *M2g_E3_[20];
