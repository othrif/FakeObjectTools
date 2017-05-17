#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <iostream>
#include <algorithm>
#include <sstream>
using namespace std;

void ytRelative_difference_of_efficiency(float pt_low = 10., float pt_up = 200., bool rebin = true, bool eta_less_than_two = false)
{
    stringstream ss_pt_low, ss_pt_up;
    ss_pt_low << pt_low;
    ss_pt_up << pt_up;
    string pt_range = ss_pt_low.str() + " GeV < pT < " + ss_pt_up.str() + " GeV";

    TString path = "/Users/ytshen/Documents/Working/OU/HEP/my_codes/SS_v50/Results/20170310/";

    TFile *data_Zee = TFile::Open(path + "hist-RLE-merged-data-elec.root");
    TFile *data_Zmumu = TFile::Open(path + "hist-RLE-merged-data-muon.root");
    TFile *mc_ttbar_elec = TFile::Open(path + "RLE_MC_ttbar_electron/hist-20170310.root");
    TFile *mc_ttbar_muon = TFile::Open(path + "RLE_MC_ttbar_muon/hist-20170310.root");
    TFile *mc_GG_ttn1_elec = TFile::Open(path + "RLE_MC_GG_ttn1_electron_20170314/hist-20170314.root");
    TFile *mc_GG_ttn1_muon = TFile::Open(path + "RLE_MC_GG_ttn1_muon_20170314/hist-20170314.root");

    // Get 3dim histogrmas
    TH3F *h_baseline_pt_dRjet_mll_Zee = (TH3F *)data_Zee->Get("h_baseline_pt_dRjet_mll");
    TH3F *h_baseline_pt_dRjet_mll_Zmumu = (TH3F *)data_Zmumu->Get("h_baseline_pt_dRjet_mll");
    TH3F *h_baseline_pt_dRjet_mll_ttbar_elec = (TH3F *)mc_ttbar_elec->Get("h_baseline_pt_dRjet_mll");
    TH3F *h_baseline_pt_dRjet_mll_ttbar_muon = (TH3F *)mc_ttbar_muon->Get("h_baseline_pt_dRjet_mll");
    TH3F *h_baseline_pt_dRjet_mll_GG_ttn1_elec = (TH3F *)mc_GG_ttn1_elec->Get("h_baseline_pt_dRjet_mll");
    TH3F *h_baseline_pt_dRjet_mll_GG_ttn1_muon = (TH3F *)mc_GG_ttn1_muon->Get("h_baseline_pt_dRjet_mll");

    h_baseline_pt_dRjet_mll_Zee->Sumw2();
    h_baseline_pt_dRjet_mll_Zmumu->Sumw2();
    h_baseline_pt_dRjet_mll_ttbar_elec->Sumw2();
    h_baseline_pt_dRjet_mll_ttbar_muon->Sumw2();
    h_baseline_pt_dRjet_mll_GG_ttn1_elec->Sumw2();
    h_baseline_pt_dRjet_mll_GG_ttn1_muon->Sumw2();

    TH3F *h_signal_pt_dRjet_mll_Zee = (TH3F *)data_Zee->Get("h_signal_pt_dRjet_mll");
    TH3F *h_signal_pt_dRjet_mll_Zmumu = (TH3F *)data_Zmumu->Get("h_signal_pt_dRjet_mll");
    TH3F *h_signal_pt_dRjet_mll_ttbar_elec = (TH3F *)mc_ttbar_elec->Get("h_signal_pt_dRjet_mll");
    TH3F *h_signal_pt_dRjet_mll_ttbar_muon = (TH3F *)mc_ttbar_muon->Get("h_signal_pt_dRjet_mll");
    TH3F *h_signal_pt_dRjet_mll_GG_ttn1_elec = (TH3F *)mc_GG_ttn1_elec->Get("h_signal_pt_dRjet_mll");
    TH3F *h_signal_pt_dRjet_mll_GG_ttn1_muon = (TH3F *)mc_GG_ttn1_muon->Get("h_signal_pt_dRjet_mll");

    h_signal_pt_dRjet_mll_Zee->Sumw2();
    h_signal_pt_dRjet_mll_Zmumu->Sumw2();
    h_signal_pt_dRjet_mll_ttbar_elec->Sumw2();
    h_signal_pt_dRjet_mll_ttbar_muon->Sumw2();
    h_signal_pt_dRjet_mll_GG_ttn1_elec->Sumw2();
    h_signal_pt_dRjet_mll_GG_ttn1_muon->Sumw2();

    int pt_bin_low = h_baseline_pt_dRjet_mll_Zmumu->GetXaxis()->FindBin(pt_low + 0.01);
    int pt_bin_up = h_baseline_pt_dRjet_mll_Zmumu->GetXaxis()->FindBin(pt_up - 0.01);

    // Project onto dRjet axis
    TH1F *h_baseline_dRjet_Zee = 0, *h_baseline_dRjet_ttbar_elec= 0, *h_baseline_dRjet_GG_ttn1_elec = 0;
    TH1F *h_signal_dRjet_Zee = 0, *h_signal_dRjet_ttbar_elec = 0, *h_signal_dRjet_GG_ttn1_elec = 0;

    if (!eta_less_than_two) { // electron uses all |eta| range
        h_baseline_dRjet_Zee = (TH1F *)h_baseline_pt_dRjet_mll_Zee->ProjectionY("h_baseline_dRjet_Zee", pt_bin_low, pt_bin_up)->Clone();
        h_baseline_dRjet_ttbar_elec = (TH1F *)h_baseline_pt_dRjet_mll_ttbar_elec->ProjectionY("h_baseline_dRjet_ttbar_elec", pt_bin_low, pt_bin_up)->Clone();
        h_baseline_dRjet_GG_ttn1_elec = (TH1F *)h_baseline_pt_dRjet_mll_GG_ttn1_elec->ProjectionY("h_baseline_dRjet_GG_ttn1_elec", pt_bin_low, pt_bin_up)->Clone();

        h_signal_dRjet_Zee = (TH1F *)h_signal_pt_dRjet_mll_Zee->ProjectionY("h_signal_dRjet_Zee", pt_bin_low, pt_bin_up)->Clone();
        h_signal_dRjet_ttbar_elec = (TH1F *)h_signal_pt_dRjet_mll_ttbar_elec->ProjectionY("h_signal_dRjet_ttbar_elec", pt_bin_low, pt_bin_up)->Clone();
        h_signal_dRjet_GG_ttn1_elec = (TH1F *)h_signal_pt_dRjet_mll_GG_ttn1_elec->ProjectionY("h_signal_dRjet_GG_ttn1_elec", pt_bin_low, pt_bin_up)->Clone();
    }
    else { // electron uses 0 < |eta| < 2.0
        TH3F *h_baseline_pt_eta_dRjet_Zee = (TH3F *)data_Zee->Get("h_baseline_pt_eta_dRjet");
        TH3F *h_baseline_pt_eta_dRjet_ttbar_elec = (TH3F *)mc_ttbar_elec->Get("h_baseline_pt_eta_dRjet");
        TH3F *h_baseline_pt_eta_dRjet_GG_ttn1_elec = (TH3F *)mc_GG_ttn1_elec->Get("h_baseline_pt_eta_dRjet");

        h_baseline_pt_eta_dRjet_Zee->Sumw2();
        h_baseline_pt_eta_dRjet_ttbar_elec->Sumw2();
        h_baseline_pt_eta_dRjet_GG_ttn1_elec->Sumw2();

        TH3F *h_signal_pt_eta_dRjet_Zee = (TH3F *)data_Zee->Get("h_signal_pt_eta_dRjet");
        TH3F *h_signal_pt_eta_dRjet_ttbar_elec = (TH3F *)mc_ttbar_elec->Get("h_signal_pt_eta_dRjet");
        TH3F *h_signal_pt_eta_dRjet_GG_ttn1_elec = (TH3F *)mc_GG_ttn1_elec->Get("h_signal_pt_eta_dRjet");

        h_signal_pt_eta_dRjet_Zee->Sumw2();
        h_signal_pt_eta_dRjet_ttbar_elec->Sumw2();
        h_signal_pt_eta_dRjet_GG_ttn1_elec->Sumw2();

        // For electron, we use 0 < |eta| < 2 only.
        int eta_bin_low = h_baseline_pt_eta_dRjet_Zee->GetYaxis()->FindBin(0 + 0.01);
        int eta_bin_up = h_baseline_pt_eta_dRjet_Zee->GetYaxis()->FindBin(2 - 0.01);

        h_baseline_dRjet_Zee = (TH1F *)h_baseline_pt_eta_dRjet_Zee->ProjectionZ("h_baseline_dRjet_Zee", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)->Clone();
        h_baseline_dRjet_ttbar_elec = (TH1F *)h_baseline_pt_eta_dRjet_ttbar_elec->ProjectionZ("h_baseline_dRjet_ttbar_elec", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)->Clone();
        h_baseline_dRjet_GG_ttn1_elec = (TH1F *)h_baseline_pt_eta_dRjet_GG_ttn1_elec->ProjectionZ("h_baseline_dRjet_GG_ttn1_elec", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)->Clone();

        h_signal_dRjet_Zee = (TH1F *)h_signal_pt_eta_dRjet_Zee->ProjectionZ("h_signal_dRjet_Zee", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)->Clone();
        h_signal_dRjet_ttbar_elec = (TH1F *)h_signal_pt_eta_dRjet_ttbar_elec->ProjectionZ("h_signal_dRjet_ttbar_elec", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)->Clone();
        h_signal_dRjet_GG_ttn1_elec = (TH1F *)h_signal_pt_eta_dRjet_GG_ttn1_elec->ProjectionZ("h_signal_dRjet_GG_ttn1_elec", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)->Clone();
    }

    // muon always use all |eta| range.
    TH1F *h_baseline_dRjet_Zmumu = (TH1F *)h_baseline_pt_dRjet_mll_Zmumu->ProjectionY("h_baseline_dRjet_Zmumu", pt_bin_low, pt_bin_up)->Clone();
    TH1F *h_baseline_dRjet_ttbar_muon = (TH1F *)h_baseline_pt_dRjet_mll_ttbar_muon->ProjectionY("h_baseline_dRjet_ttbar_muon", pt_bin_low, pt_bin_up)->Clone();
    TH1F *h_baseline_dRjet_GG_ttn1_muon = (TH1F *)h_baseline_pt_dRjet_mll_GG_ttn1_muon->ProjectionY("h_baseline_dRjet_GG_ttn1_muon", pt_bin_low, pt_bin_up)->Clone();

    TH1F *h_signal_dRjet_Zmumu = (TH1F *)h_signal_pt_dRjet_mll_Zmumu->ProjectionY("h_signal_dRjet_Zmumu", pt_bin_low, pt_bin_up)->Clone();
    TH1F *h_signal_dRjet_ttbar_muon = (TH1F *)h_signal_pt_dRjet_mll_ttbar_muon->ProjectionY("h_signal_dRjet_ttbar_muon", pt_bin_low, pt_bin_up)->Clone();
    TH1F *h_signal_dRjet_GG_ttn1_muon = (TH1F *)h_signal_pt_dRjet_mll_GG_ttn1_muon->ProjectionY("h_signal_dRjet_GG_ttn1_muon", pt_bin_low, pt_bin_up)->Clone();

    // Rebinned dR
    //double xbins_rebinned[7] = {0., 0.2, 0.4, 0.6, 1., 2., 4.};
    //double xbins_rebinned[9] = {0., 0.1, 0.2, 0.3, 0.4, 0.6, 1., 2., 4.};
    double xbins_rebinned[11] = {0., 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.6, 1., 2., 4.};
    int number_of_rebinned_xbins = sizeof(xbins_rebinned) / sizeof(double) - 1;
    //cout << "number_of_rebinned_xbins=" << number_of_rebinned_xbins << endl;

    TH1F *h_baseline_dRjet_rebinned_Zee = (TH1F *)h_baseline_dRjet_Zee->Rebin(number_of_rebinned_xbins, "h_baseline_dRjet_rebinned_Zee", xbins_rebinned);
    TH1F *h_baseline_dRjet_rebinned_ttbar_elec = (TH1F *)h_baseline_dRjet_ttbar_elec->Rebin(number_of_rebinned_xbins, "h_baseline_dRjet_rebinned_ttbar_elec", xbins_rebinned);
    TH1F *h_baseline_dRjet_rebinned_GG_ttn1_elec = (TH1F *)h_baseline_dRjet_GG_ttn1_elec->Rebin(number_of_rebinned_xbins, "h_baseline_dRjet_rebinned_GG_ttn1_elec", xbins_rebinned);

    TH1F *h_signal_dRjet_rebinned_Zee = (TH1F *)h_signal_dRjet_Zee->Rebin(number_of_rebinned_xbins, "h_signal_dRjet_rebinned_Zee", xbins_rebinned);
    TH1F *h_signal_dRjet_rebinned_ttbar_elec = (TH1F *)h_signal_dRjet_ttbar_elec->Rebin(number_of_rebinned_xbins, "h_signal_dRjet_rebinned_ttbar_elec", xbins_rebinned);
    TH1F *h_signal_dRjet_rebinned_GG_ttn1_elec = (TH1F *)h_signal_dRjet_GG_ttn1_elec->Rebin(number_of_rebinned_xbins, "h_signal_dRjet_rebinned_GG_ttn1_elec", xbins_rebinned);

    TH1F *h_baseline_dRjet_rebinned_Zmumu = (TH1F *)h_baseline_dRjet_Zmumu->Rebin(number_of_rebinned_xbins, "h_baseline_dRjet_rebinned_Zmumu", xbins_rebinned);
    TH1F *h_baseline_dRjet_rebinned_ttbar_muon = (TH1F *)h_baseline_dRjet_ttbar_muon->Rebin(number_of_rebinned_xbins, "h_baseline_dRjet_rebinned_ttbar_muon", xbins_rebinned);
    TH1F *h_baseline_dRjet_rebinned_GG_ttn1_muon = (TH1F *)h_baseline_dRjet_GG_ttn1_muon->Rebin(number_of_rebinned_xbins, "h_baseline_dRjet_rebinned_GG_ttn1_muon", xbins_rebinned);

    TH1F *h_signal_dRjet_rebinned_Zmumu = (TH1F *)h_signal_dRjet_Zmumu->Rebin(number_of_rebinned_xbins, "h_signal_dRjet_rebinned_Zmumu", xbins_rebinned);
    TH1F *h_signal_dRjet_rebinned_ttbar_muon = (TH1F *)h_signal_dRjet_ttbar_muon->Rebin(number_of_rebinned_xbins, "h_signal_dRjet_rebinned_ttbar_muon", xbins_rebinned);
    TH1F *h_signal_dRjet_rebinned_GG_ttn1_muon = (TH1F *)h_signal_dRjet_GG_ttn1_muon->Rebin(number_of_rebinned_xbins, "h_signal_dRjet_rebinned_GG_ttn1_muon", xbins_rebinned);

    // Calculate efficiencies
    TH1F *h_eff_dRjet_Zee = 0, *h_eff_dRjet_Zmumu = 0,
         *h_eff_dRjet_ttbar_elec = 0, *h_eff_dRjet_ttbar_muon = 0,
         *h_eff_dRjet_GG_ttn1_elec = 0, *h_eff_dRjet_GG_ttn1_muon = 0;
    TH1F *h_baseline_dRjet_denominator_Zee = 0, *h_baseline_dRjet_denominator_Zmumu = 0,
         *h_baseline_dRjet_denominator_ttbar_elec = 0, *h_baseline_dRjet_denominator_ttbar_muon = 0,
         *h_baseline_dRjet_denominator_GG_ttn1_elec = 0, *h_baseline_dRjet_denominator_GG_ttn1_muon = 0;

    if (rebin) {
        h_eff_dRjet_Zee = (TH1F *)h_signal_dRjet_rebinned_Zee->Clone();
        h_eff_dRjet_Zmumu = (TH1F *)h_signal_dRjet_rebinned_Zmumu->Clone();
        h_eff_dRjet_ttbar_elec = (TH1F *)h_signal_dRjet_rebinned_ttbar_elec->Clone();
        h_eff_dRjet_ttbar_muon = (TH1F *)h_signal_dRjet_rebinned_ttbar_muon->Clone();
        h_eff_dRjet_GG_ttn1_elec = (TH1F *)h_signal_dRjet_rebinned_GG_ttn1_elec->Clone();
        h_eff_dRjet_GG_ttn1_muon = (TH1F *)h_signal_dRjet_rebinned_GG_ttn1_muon->Clone();

        h_baseline_dRjet_denominator_Zee = h_baseline_dRjet_rebinned_Zee;
        h_baseline_dRjet_denominator_Zmumu = h_baseline_dRjet_rebinned_Zmumu;
        h_baseline_dRjet_denominator_ttbar_elec = h_baseline_dRjet_rebinned_ttbar_elec;
        h_baseline_dRjet_denominator_ttbar_muon = h_baseline_dRjet_rebinned_ttbar_muon;
        h_baseline_dRjet_denominator_GG_ttn1_elec = h_baseline_dRjet_rebinned_GG_ttn1_elec;
        h_baseline_dRjet_denominator_GG_ttn1_muon = h_baseline_dRjet_rebinned_GG_ttn1_muon;

        // use binomial error because signal and baseline are correlated.
        h_eff_dRjet_Zee->Divide(h_signal_dRjet_rebinned_Zee, h_baseline_dRjet_rebinned_Zee, 1, 1, "B");
        h_eff_dRjet_Zmumu->Divide(h_signal_dRjet_rebinned_Zmumu, h_baseline_dRjet_rebinned_Zmumu, 1, 1, "B");
        h_eff_dRjet_ttbar_elec->Divide(h_signal_dRjet_rebinned_ttbar_elec, h_baseline_dRjet_rebinned_ttbar_elec, 1, 1, "B");
        h_eff_dRjet_ttbar_muon->Divide(h_signal_dRjet_rebinned_ttbar_muon, h_baseline_dRjet_rebinned_ttbar_muon, 1, 1, "B");
        h_eff_dRjet_GG_ttn1_elec->Divide(h_signal_dRjet_rebinned_GG_ttn1_elec, h_baseline_dRjet_rebinned_GG_ttn1_elec, 1, 1, "B");
        h_eff_dRjet_GG_ttn1_muon->Divide(h_signal_dRjet_rebinned_GG_ttn1_muon, h_baseline_dRjet_rebinned_GG_ttn1_muon, 1, 1, "B");
    }
    else {
        h_eff_dRjet_Zee = (TH1F *)h_signal_dRjet_Zee->Clone();
        h_eff_dRjet_Zmumu = (TH1F *)h_signal_dRjet_Zmumu->Clone();
        h_eff_dRjet_ttbar_elec = (TH1F *)h_signal_dRjet_ttbar_elec->Clone();
        h_eff_dRjet_ttbar_muon = (TH1F *)h_signal_dRjet_ttbar_muon->Clone();
        h_eff_dRjet_GG_ttn1_elec = (TH1F *)h_signal_dRjet_GG_ttn1_elec->Clone();
        h_eff_dRjet_GG_ttn1_muon = (TH1F *)h_signal_dRjet_GG_ttn1_muon->Clone();

        h_baseline_dRjet_denominator_Zee = h_baseline_dRjet_Zee;
        h_baseline_dRjet_denominator_Zmumu = h_baseline_dRjet_Zmumu;
        h_baseline_dRjet_denominator_ttbar_elec = h_baseline_dRjet_ttbar_elec;
        h_baseline_dRjet_denominator_ttbar_muon = h_baseline_dRjet_ttbar_muon;
        h_baseline_dRjet_denominator_GG_ttn1_elec = h_baseline_dRjet_GG_ttn1_elec;
        h_baseline_dRjet_denominator_GG_ttn1_muon = h_baseline_dRjet_GG_ttn1_muon;

        // use binomial error because signal and baseline are correlated.
        h_eff_dRjet_Zee->Divide(h_signal_dRjet_Zee, h_baseline_dRjet_Zee, 1, 1, "B");
        h_eff_dRjet_Zmumu->Divide(h_signal_dRjet_Zmumu, h_baseline_dRjet_Zmumu, 1, 1, "B");
        h_eff_dRjet_ttbar_elec->Divide(h_signal_dRjet_ttbar_elec, h_baseline_dRjet_ttbar_elec, 1, 1, "B");
        h_eff_dRjet_ttbar_muon->Divide(h_signal_dRjet_ttbar_muon, h_baseline_dRjet_ttbar_muon, 1, 1, "B");
        h_eff_dRjet_GG_ttn1_elec->Divide(h_signal_dRjet_GG_ttn1_elec, h_baseline_dRjet_GG_ttn1_elec, 1, 1, "B");
        h_eff_dRjet_GG_ttn1_muon->Divide(h_signal_dRjet_GG_ttn1_muon, h_baseline_dRjet_GG_ttn1_muon, 1, 1, "B");
    }

    // Calculate the average efficiencies of Zee and Zmumu
    // Use N(signal)/N(baseline) to calculate the average
    int Nbins_baseline_Zee = h_baseline_dRjet_Zee->GetXaxis()->GetNbins();
    float Number_of_baseline_Zee = 0, Number_of_signal_Zee = 0;
    for (int i = 1; i <= Nbins_baseline_Zee + 1; i++) { // consider the overflow bins
        Number_of_baseline_Zee += h_baseline_dRjet_Zee->GetBinContent(i);
        Number_of_signal_Zee += h_signal_dRjet_Zee->GetBinContent(i);
    }
    float average_eff_Zee = Number_of_signal_Zee / Number_of_baseline_Zee;
    cout << "average_eff_Zee=" << average_eff_Zee << endl;
    // Fill the average efficiency histogram
    TH1F *h_average_eff_Zee = (TH1F *)h_eff_dRjet_Zee->Clone();
    h_average_eff_Zee->Sumw2();
    h_average_eff_Zee->Reset();
    for (int i = 1; i <= Nbins_baseline_Zee + 1; i++) { // consider the overflow bins
        h_average_eff_Zee->SetBinContent(i, average_eff_Zee);
    }

    int Nbins_baseline_Zmumu = h_baseline_dRjet_Zmumu->GetXaxis()->GetNbins();
    float Number_of_baseline_Zmumu = 0, Number_of_signal_Zmumu = 0;
    for (int i = 1; i <= Nbins_baseline_Zmumu + 1; i++) { // consider the overflow bins
        Number_of_baseline_Zmumu += h_baseline_dRjet_Zmumu->GetBinContent(i);
        Number_of_signal_Zmumu += h_signal_dRjet_Zmumu->GetBinContent(i);
    }
    float average_eff_Zmumu = Number_of_signal_Zmumu / Number_of_baseline_Zmumu;
    cout << "average_eff_Zmumu=" << average_eff_Zmumu << endl;
    // Fill the average efficiency histogram
    TH1F *h_average_eff_Zmumu = (TH1F *)h_eff_dRjet_Zmumu->Clone();
    h_average_eff_Zmumu->Sumw2();
    h_average_eff_Zmumu->Reset();
    for (int i = 1; i <= Nbins_baseline_Zmumu + 1; i ++) { // consider the overflow bins
        h_average_eff_Zmumu->SetBinContent(i, average_eff_Zmumu);
    }

    TH1F *h_average_eff_rebinned_Zee = 0, *h_average_eff_rebinned_Zmumu = 0;
    if (rebin) {
        h_average_eff_rebinned_Zee = (TH1F *)h_average_eff_Zee->Rebin(number_of_rebinned_xbins, "h_average_eff_rebinned_Zee", xbins_rebinned);
        h_average_eff_rebinned_Zmumu = (TH1F *)h_average_eff_Zmumu->Rebin(number_of_rebinned_xbins, "h_average_eff_rebinned_Zmumu", xbins_rebinned);
    }

    // Calculate the relative difference of the efficiencies
    TH1F *h_efficiency_difference_ttbar_elec = (TH1F *)h_eff_dRjet_ttbar_elec->Clone();
    h_efficiency_difference_ttbar_elec->Add(h_average_eff_Zee, -1);

    TH1F *h_relative_difference_eff_ttbar_elec = (TH1F *)h_efficiency_difference_ttbar_elec->Clone();
    h_relative_difference_eff_ttbar_elec->Reset();
    for (int i = 1; i <= h_efficiency_difference_ttbar_elec->GetXaxis()->GetNbins(); i++) {
        float bin_content = h_efficiency_difference_ttbar_elec->GetBinContent(i);
        float bin_error = h_efficiency_difference_ttbar_elec->GetBinError(i);
        // Make the absolute value
        if (bin_content < 0.) {
            h_relative_difference_eff_ttbar_elec->SetBinContent(i, fabs(bin_content));
            h_relative_difference_eff_ttbar_elec->SetBinError(i, bin_error);
        }
        else {
            h_relative_difference_eff_ttbar_elec->SetBinContent(i, bin_content);
            h_relative_difference_eff_ttbar_elec->SetBinError(i, bin_error);
        }
    }
    h_relative_difference_eff_ttbar_elec->Divide(h_average_eff_Zee);

    cout << "Relative difference of efficiency ttbar_elec:" << endl;
    cout << "Nbins = " << h_relative_difference_eff_ttbar_elec->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_relative_difference_eff_ttbar_elec->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_relative_difference_eff_ttbar_elec->GetBinContent(i);
        float bin_error = h_relative_difference_eff_ttbar_elec->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
/*
    double xbins_ttbar_elec[7] = {0., 0.2, 0.4, 0.6, 1., 2., 4.};
    TH1F *h_rebin_ttbar_elec = (TH1F *)h_relative_difference_eff_ttbar_elec->Rebin(6, "h_rebin_ttbar_elec", xbins_ttbar_elec);

    cout << "Rebinned Relative difference of efficiency ttbar_elec:" << endl;
    cout << "Nbins = " << h_rebin_ttbar_elec->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_rebin_ttbar_elec->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_rebin_ttbar_elec->GetBinContent(i);
        float bin_error = h_rebin_ttbar_elec->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
*/
    TH1F *h_efficiency_difference_ttbar_muon = (TH1F *)h_eff_dRjet_ttbar_muon->Clone();
    h_efficiency_difference_ttbar_muon->Add(h_average_eff_Zmumu, -1);

    TH1F *h_relative_difference_eff_ttbar_muon = (TH1F *)h_efficiency_difference_ttbar_muon->Clone();
    h_relative_difference_eff_ttbar_muon->Reset();
    for (int i = 1; i <= h_efficiency_difference_ttbar_muon->GetXaxis()->GetNbins(); i++) {
        float bin_content = h_efficiency_difference_ttbar_muon->GetBinContent(i);
        float bin_error = h_efficiency_difference_ttbar_muon->GetBinError(i);
        // Make the absolute value
        if (bin_content < 0.) {
            h_relative_difference_eff_ttbar_muon->SetBinContent(i, fabs(bin_content));
            h_relative_difference_eff_ttbar_muon->SetBinError(i, bin_error);
        }
        else {
            h_relative_difference_eff_ttbar_muon->SetBinContent(i, bin_content);
            h_relative_difference_eff_ttbar_muon->SetBinError(i, bin_error);
        }
    }
    h_relative_difference_eff_ttbar_muon->Divide(h_average_eff_Zmumu);

    cout << "Relative difference of efficiency ttbar_muon:" << endl;
    cout << "Nbins = " << h_relative_difference_eff_ttbar_muon->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_relative_difference_eff_ttbar_muon->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_relative_difference_eff_ttbar_muon->GetBinContent(i);
        float bin_error = h_relative_difference_eff_ttbar_muon->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
/*
	double xbins_ttbar_muon[7] = {0., 0.2, 0.4, 0.6, 1., 2., 4.};
	TH1F *h_rebin_ttbar_muon = (TH1F *)h_relative_difference_eff_ttbar_muon->Rebin(6, "h_rebin_ttbar_muon", xbins_ttbar_muon);

	cout << "Rebinned Relative difference of efficiency ttbar_muon:" << endl;
	cout << "Nbins = " << h_rebin_ttbar_muon->GetXaxis()->GetNbins() << endl;
	for (int i = 0; i <= h_rebin_ttbar_muon->GetXaxis()->GetNbins() + 1; i++) {
		float bin_content = h_rebin_ttbar_muon->GetBinContent(i);
		float bin_error = h_rebin_ttbar_muon->GetBinError(i);
		cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
	}
	cout << endl;
*/
    TH1F *h_efficiency_difference_GG_ttn1_elec = (TH1F *)h_eff_dRjet_GG_ttn1_elec->Clone();
    h_efficiency_difference_GG_ttn1_elec->Add(h_average_eff_Zee, -1);

    TH1F *h_relative_difference_eff_GG_ttn1_elec = (TH1F *)h_efficiency_difference_GG_ttn1_elec->Clone();
    h_relative_difference_eff_GG_ttn1_elec->Reset();
    for (int i = 1; i <= h_efficiency_difference_GG_ttn1_elec->GetXaxis()->GetNbins(); i++) {
        float bin_content = h_efficiency_difference_GG_ttn1_elec->GetBinContent(i);
        float bin_error = h_efficiency_difference_GG_ttn1_elec->GetBinError(i);
        if (bin_content < 0.) {
            h_relative_difference_eff_GG_ttn1_elec->SetBinContent(i, fabs(bin_content));
            h_relative_difference_eff_GG_ttn1_elec->SetBinError(i, bin_error);
        }
        else {
            h_relative_difference_eff_GG_ttn1_elec->SetBinContent(i, bin_content);
            h_relative_difference_eff_GG_ttn1_elec->SetBinError(i, bin_error);
        }
    }
    h_relative_difference_eff_GG_ttn1_elec->Divide(h_average_eff_Zee);

    cout << "Relative difference of efficiency GG_ttn1_elec:" << endl;
    cout << "Nbins = " << h_relative_difference_eff_GG_ttn1_elec->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_relative_difference_eff_GG_ttn1_elec->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_relative_difference_eff_GG_ttn1_elec->GetBinContent(i);
        float bin_error = h_relative_difference_eff_GG_ttn1_elec->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
/*
    double xbins_GG_ttn1_elec[7] = {0., 0.2, 0.4, 0.6, 1., 2., 4.};
    TH1F *h_rebin_GG_ttn1_elec = (TH1F *)h_relative_difference_eff_GG_ttn1_elec->Rebin(6, "h_rebin_GG_ttn1_elec", xbins_GG_ttn1_elec);

    cout << "Rebinned Relative difference of efficiency GG_ttn1_elec:" << endl;
    cout << "Nbins = " << h_rebin_GG_ttn1_elec->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_rebin_GG_ttn1_elec->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_rebin_GG_ttn1_elec->GetBinContent(i);
        float bin_error = h_rebin_GG_ttn1_elec->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
*/
    TH1F *h_efficiency_difference_GG_ttn1_muon = (TH1F *)h_eff_dRjet_GG_ttn1_muon->Clone();
    h_efficiency_difference_GG_ttn1_muon->Add(h_average_eff_Zmumu, -1);

    TH1F *h_relative_difference_eff_GG_ttn1_muon = (TH1F *)h_efficiency_difference_GG_ttn1_muon->Clone();
    h_relative_difference_eff_GG_ttn1_muon->Reset();
    for (int i = 1; i <= h_efficiency_difference_GG_ttn1_muon->GetXaxis()->GetNbins(); i++) {
        float bin_content = h_efficiency_difference_GG_ttn1_muon->GetBinContent(i);
        float bin_error = h_efficiency_difference_GG_ttn1_muon->GetBinError(i);
        if (bin_content < 0.) {
            h_relative_difference_eff_GG_ttn1_muon->SetBinContent(i, fabs(bin_content));
            h_relative_difference_eff_GG_ttn1_muon->SetBinError(i, bin_error);
        }
        else {
            h_relative_difference_eff_GG_ttn1_muon->SetBinContent(i, bin_content);
            h_relative_difference_eff_GG_ttn1_muon->SetBinError(i, bin_error);
        }
    }
    h_relative_difference_eff_GG_ttn1_muon->Divide(h_average_eff_Zmumu);

    cout << "Relative difference of efficiency GG_ttn1_muon:" << endl;
    cout << "Nbins = " << h_relative_difference_eff_GG_ttn1_muon->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_relative_difference_eff_GG_ttn1_muon->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_relative_difference_eff_GG_ttn1_muon->GetBinContent(i);
        float bin_error = h_relative_difference_eff_GG_ttn1_muon->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
/*
    double xbins_GG_ttn1_muon[7] = {0., 0.2, 0.4, 0.6, 1., 2., 4.};
    TH1F *h_rebin_GG_ttn1_muon = (TH1F *)h_relative_difference_eff_GG_ttn1_muon->Rebin(6, "h_rebin_GG_ttn1_muon", xbins_GG_ttn1_muon);

    cout << "Rebinned Relative difference of efficiency GG_ttn1_elec:" << endl;
    cout << "Nbins = " << h_rebin_GG_ttn1_muon->GetXaxis()->GetNbins() << endl;
    for (int i = 0; i <= h_rebin_GG_ttn1_muon->GetXaxis()->GetNbins() + 1; i++) {
        float bin_content = h_rebin_GG_ttn1_muon->GetBinContent(i);
        float bin_error = h_rebin_GG_ttn1_muon->GetBinError(i);
        cout << "bin=" << i << ": " << bin_content << " $pm$ " << bin_error << endl;
    }
    cout << endl;
*/
    // make plots electron
    TCanvas *efficiencies_ee = new TCanvas("efficiencies_ee", "", 600, 600);
    efficiencies_ee->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_eff_dRjet_Zee->SetTitle("");
    h_eff_dRjet_Zee->SetXTitle("#Delta R(e, jet)");
    h_eff_dRjet_Zee->SetYTitle("Electron Real Efficiencies");
    h_eff_dRjet_Zee->GetYaxis()->SetTitleOffset(1.5);
    h_eff_dRjet_Zee->SetMaximum(1.2);
    h_eff_dRjet_Zee->SetMinimum(0.);
    h_eff_dRjet_Zee->SetMarkerStyle(kFullCircle);
    h_eff_dRjet_Zee->SetMarkerColor(kBlack);
    h_eff_dRjet_Zee->SetLineColor(kBlack);
    h_eff_dRjet_Zee->Draw("E0");
    if (rebin) {
        h_average_eff_rebinned_Zee->SetMarkerStyle(kFullCircle);
        h_average_eff_rebinned_Zee->SetMarkerColor(kBlack);
        h_average_eff_rebinned_Zee->SetLineColor(kBlack);
        h_average_eff_rebinned_Zee->SetLineStyle(2);
        h_average_eff_rebinned_Zee->Draw("same");
    }
    else {
        h_average_eff_Zee->SetMarkerStyle(kFullCircle);
        h_average_eff_Zee->SetMarkerColor(kBlack);
        h_average_eff_Zee->SetLineColor(kBlack);
        h_average_eff_Zee->SetLineStyle(2);
        h_average_eff_Zee->Draw("same");
    }
    h_eff_dRjet_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_eff_dRjet_ttbar_elec->SetMarkerColor(kMagenta);
    h_eff_dRjet_ttbar_elec->SetLineColor(kMagenta);
    h_eff_dRjet_ttbar_elec->Draw("E0,same");
    h_eff_dRjet_GG_ttn1_elec->SetMarkerStyle(kFullCross);
    h_eff_dRjet_GG_ttn1_elec->SetMarkerColor(kOrange);
    h_eff_dRjet_GG_ttn1_elec->SetLineColor(kOrange);
    h_eff_dRjet_GG_ttn1_elec->Draw("E0,same");
    TPaveText *text1 = new TPaveText(0.3, 0.4, 0.6, 0.45, "NDC");
    text1->AddText(pt_range.c_str());
    text1->SetTextSize(0.03);
    text1->Draw("same");
    TLegend *leg1 = new TLegend(0.3, 0.2, 0.6, 0.4);
    leg1->AddEntry(h_eff_dRjet_Zee, "Z#rightarrow ee", "lp");
    if (rebin)
        leg1->AddEntry(h_average_eff_rebinned_Zee, "Z#rightarrow ee (average, before rebinning)", "l");
    else
        leg1->AddEntry(h_average_eff_Zee, "Z#rightarrow ee (average, before rebinning)", "l");
    leg1->AddEntry(h_eff_dRjet_ttbar_elec, "t#bar{t}", "lp");
    leg1->AddEntry(h_eff_dRjet_GG_ttn1_elec, "#tilde{g}#rightarrow t#bar{t}", "lp");
    leg1->SetTextSize(0.03);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");
    string filename_efficiencies_ee;
    if (eta_less_than_two)
        filename_efficiencies_ee = "real_efficiency_electron_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + "_eta020.pdf";
    else
        filename_efficiencies_ee = "real_efficiency_electron_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + ".pdf";
    efficiencies_ee->SaveAs(filename_efficiencies_ee.c_str(), "pdf");

    TCanvas *efficiency_difference_ee = new TCanvas("efficiency_difference_ee", "", 600, 600);
    efficiency_difference_ee->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_efficiency_difference_ttbar_elec->SetTitle("");
    h_efficiency_difference_ttbar_elec->SetXTitle("#Delta R(e, jet)");
    h_efficiency_difference_ttbar_elec->SetYTitle("|#Delta Efficiencies|");
    h_efficiency_difference_ttbar_elec->GetYaxis()->SetTitleOffset(1.5);
    h_efficiency_difference_ttbar_elec->SetMaximum(1.);
    h_efficiency_difference_ttbar_elec->SetMinimum(-1.);
    h_efficiency_difference_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_efficiency_difference_ttbar_elec->SetMarkerColor(kMagenta);
    h_efficiency_difference_ttbar_elec->SetLineColor(kMagenta);
    h_efficiency_difference_ttbar_elec->Draw();
    h_efficiency_difference_GG_ttn1_elec->SetMarkerStyle(kFullCross);
    h_efficiency_difference_GG_ttn1_elec->SetMarkerColor(kOrange);
    h_efficiency_difference_GG_ttn1_elec->SetLineColor(kOrange);
    h_efficiency_difference_GG_ttn1_elec->Draw("same");
    TLine *line2 = new TLine(0, 0, 4, 0);
    line2->SetLineColor(kBlack);
    line2->SetLineStyle(2);
    line2->Draw("same");
    TPaveText *text2 = new TPaveText(0.6, 0.8, 0.8, 0.85, "NDC");
    text2->AddText(pt_range.c_str());
    text2->SetTextSize(0.03);
    text2->Draw("same");
    TLegend *leg2 = new TLegend(0.6, 0.7, 0.8, 0.8);
    leg2->AddEntry(h_efficiency_difference_ttbar_elec, "t#bar{t}", "lp");
    leg2->AddEntry(h_efficiency_difference_GG_ttn1_elec, "#tilde{g}#rightarrow t#bar{t}", "lp");
    leg2->SetTextSize(0.03);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->Draw("same");
    string filename_efficiency_difference_ee;
    if (eta_less_than_two)
        filename_efficiency_difference_ee = "real_efficiency_difference_electron_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + "_eta020.pdf";
    else
        filename_efficiency_difference_ee = "real_efficiency_difference_electron_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + ".pdf";
    efficiency_difference_ee->SaveAs(filename_efficiency_difference_ee.c_str(), "pdf");

    TCanvas *relative_difference_eff_ee = new TCanvas("relative_difference_eff_ee", "", 600, 600);
    relative_difference_eff_ee->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_relative_difference_eff_ttbar_elec->Scale(100.);
    h_relative_difference_eff_ttbar_elec->SetTitle("");
    h_relative_difference_eff_ttbar_elec->SetXTitle("#Delta R(e, jet)");
    h_relative_difference_eff_ttbar_elec->SetYTitle("Relative difference of efficiencies (%)");
    h_relative_difference_eff_ttbar_elec->GetYaxis()->SetTitleOffset(1.5);
    h_relative_difference_eff_ttbar_elec->SetMaximum(101.);
    h_relative_difference_eff_ttbar_elec->SetMinimum(0.);
    h_relative_difference_eff_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_relative_difference_eff_ttbar_elec->SetMarkerColor(kMagenta);
    h_relative_difference_eff_ttbar_elec->SetLineColor(kMagenta);
    h_relative_difference_eff_ttbar_elec->Draw();
    //h_rebin_ttbar_elec->SetMarkerStyle(kOpenDiamond);
    //h_rebin_ttbar_elec->SetMarkerColor(kMagenta+3);
    //h_rebin_ttbar_elec->SetLineColor(kMagenta+3);
    //h_rebin_ttbar_elec->Draw("same");
    h_relative_difference_eff_GG_ttn1_elec->Scale(100.);
    h_relative_difference_eff_GG_ttn1_elec->SetMarkerStyle(kFullCross);
    h_relative_difference_eff_GG_ttn1_elec->SetMarkerColor(kOrange);
    h_relative_difference_eff_GG_ttn1_elec->SetLineColor(kOrange);
    h_relative_difference_eff_GG_ttn1_elec->Draw("same");
    //h_rebin_GG_ttn1_elec->SetMarkerStyle(kOpenCross);
    //h_rebin_GG_ttn1_elec->SetMarkerColor(kOrange+3);
    //h_rebin_GG_ttn1_elec->SetLineColor(kOrange+3);
    //h_rebin_GG_ttn1_elec->Draw("same");
    TPaveText *text3 = new TPaveText(0.6, 0.8, 0.8, 0.85, "NDC");
    text3->AddText(pt_range.c_str());
    text3->SetTextSize(0.03);
    text3->Draw("same");
    TLegend *leg3 = new TLegend(0.6, 0.7, 0.8, 0.8);
    leg3->AddEntry(h_relative_difference_eff_ttbar_elec, "t#bar{t}", "lp");
    leg3->AddEntry(h_relative_difference_eff_GG_ttn1_elec, "#tilde{g}#rightarrow t#bar{t}", "lp");
    //leg3->AddEntry(h_rebin_ttbar_elec, "t#bar{t} (Rebinned)", "lp");
    //leg3->AddEntry(h_rebin_GG_ttn1_elec, "#tilde{g}#rightarrow t#bar{t} (Rebinned)", "lp");
    leg3->SetTextSize(0.03);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(0);
    leg3->Draw("same");
    string filename_relative_difference_eff_ee;
    if (eta_less_than_two)
        filename_relative_difference_eff_ee = "real_efficiency_relative_difference_electron_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + "_eta020.pdf";
    else
        filename_relative_difference_eff_ee = "real_efficiency_relative_difference_electron_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + ".pdf";
    relative_difference_eff_ee->SaveAs(filename_relative_difference_eff_ee.c_str(), "pdf");

    // make plots muon
    TCanvas *efficiencies_mumu = new TCanvas("efficiencies_mumu", "", 600, 600);
    efficiencies_mumu->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_eff_dRjet_Zmumu->SetTitle("");
    h_eff_dRjet_Zmumu->SetXTitle("#Delta R(#mu, jet)");
    h_eff_dRjet_Zmumu->SetYTitle("Muon Real Efficiencies");
    h_eff_dRjet_Zmumu->GetYaxis()->SetTitleOffset(1.5);
    h_eff_dRjet_Zmumu->SetMaximum(1.1);
    h_eff_dRjet_Zmumu->SetMinimum(0.);
    h_eff_dRjet_Zmumu->SetMarkerStyle(kFullCircle);
    h_eff_dRjet_Zmumu->SetMarkerColor(kBlack);
    h_eff_dRjet_Zmumu->SetLineColor(kBlack);
    h_eff_dRjet_Zmumu->Draw("E0");
    if (rebin) {
        h_average_eff_rebinned_Zmumu->SetMarkerStyle(kFullCircle);
        h_average_eff_rebinned_Zmumu->SetMarkerColor(kBlack);
        h_average_eff_rebinned_Zmumu->SetLineColor(kBlack);
        h_average_eff_rebinned_Zmumu->SetLineStyle(2);
        h_average_eff_rebinned_Zmumu->Draw("same");
    }
    else {
        h_average_eff_Zmumu->SetMarkerStyle(kFullCircle);
        h_average_eff_Zmumu->SetMarkerColor(kBlack);
        h_average_eff_Zmumu->SetLineColor(kBlack);
        h_average_eff_Zmumu->SetLineStyle(2);
        h_average_eff_Zmumu->Draw("same");
    }
    h_eff_dRjet_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_eff_dRjet_ttbar_muon->SetMarkerColor(kMagenta);
    h_eff_dRjet_ttbar_muon->SetLineColor(kMagenta);
    h_eff_dRjet_ttbar_muon->Draw("E0,same");
    h_eff_dRjet_GG_ttn1_muon->SetMarkerStyle(kFullCross);
    h_eff_dRjet_GG_ttn1_muon->SetMarkerColor(kOrange);
    h_eff_dRjet_GG_ttn1_muon->SetLineColor(kOrange);
    h_eff_dRjet_GG_ttn1_muon->Draw("E0,same");
    TPaveText *text4 = new TPaveText(0.3, 0.4, 0.6, 0.45, "NDC");
    text4->AddText(pt_range.c_str());
    text4->SetTextSize(0.03);
    text4->Draw("same");
    TLegend *leg4 = new TLegend(0.3, 0.2, 0.6, 0.4);
    leg4->AddEntry(h_eff_dRjet_Zmumu, "Z#rightarrow #mu#mu", "lp");
    if (rebin)
        leg4->AddEntry(h_average_eff_rebinned_Zmumu, "Z#rightarrow #mu#mu (average, before rebinning)", "l");
    else
        leg4->AddEntry(h_average_eff_Zmumu, "Z#rightarrow #mu#mu (average, before rebinning)", "l");
    leg4->AddEntry(h_eff_dRjet_ttbar_muon, "t#bar{t}", "lp");
    leg4->AddEntry(h_eff_dRjet_GG_ttn1_muon, "#tilde{g}#rightarrow t#bar{t}", "lp");
    leg4->SetTextSize(0.03);
    leg4->SetBorderSize(0);
    leg4->SetFillColor(0);
    leg4->SetFillStyle(0);
    leg4->Draw("same");
    string filename_efficiencies_mumu = "real_efficiency_muon_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + ".pdf";
    efficiencies_mumu->SaveAs(filename_efficiencies_mumu.c_str(), "pdf");

    TCanvas *efficiency_difference_mumu = new TCanvas("efficiency_difference_mumu", "", 600, 600);
    efficiency_difference_mumu->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_efficiency_difference_ttbar_muon->SetTitle("");
    h_efficiency_difference_ttbar_muon->SetXTitle("#Delta R(#mu, jet)");
    h_efficiency_difference_ttbar_muon->SetYTitle("|#Delta Efficiencies|");
    h_efficiency_difference_ttbar_muon->GetYaxis()->SetTitleOffset(1.5);
    h_efficiency_difference_ttbar_muon->SetMaximum(1.);
    h_efficiency_difference_ttbar_muon->SetMinimum(-1.);
    h_efficiency_difference_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_efficiency_difference_ttbar_muon->SetMarkerColor(kMagenta);
    h_efficiency_difference_ttbar_muon->SetLineColor(kMagenta);
    h_efficiency_difference_ttbar_muon->Draw();
    h_efficiency_difference_GG_ttn1_muon->SetMarkerStyle(kFullCross);
    h_efficiency_difference_GG_ttn1_muon->SetMarkerColor(kOrange);
    h_efficiency_difference_GG_ttn1_muon->SetLineColor(kOrange);
    h_efficiency_difference_GG_ttn1_muon->Draw("same");
    TLine *line5 = new TLine(0, 0, 4, 0);
    line5->SetLineColor(kBlack);
    line5->SetLineStyle(2);
    line5->Draw("same");
    TPaveText *text5 = new TPaveText(0.6, 0.8, 0.8, 0.85, "NDC");
    text5->AddText(pt_range.c_str());
    text5->SetTextSize(0.03);
    text5->Draw("same");
    TLegend *leg5 = new TLegend(0.6, 0.7, 0.8, 0.8);
    leg5->AddEntry(h_efficiency_difference_ttbar_muon, "t#bar{t}", "lp");
    leg5->AddEntry(h_efficiency_difference_GG_ttn1_muon, "#tilde{g}#rightarrow t#bar{t}", "lp");
    leg5->SetTextSize(0.03);
    leg5->SetBorderSize(0);
    leg5->SetFillColor(0);
    leg5->SetFillStyle(0);
    leg5->Draw("same");
    string filename_efficiency_difference_mumu = "real_efficiency_difference_muon_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + ".pdf";
    efficiency_difference_mumu->SaveAs(filename_efficiency_difference_mumu.c_str(), "pdf");

    TCanvas *relative_difference_eff_mumu = new TCanvas("relative_difference_eff_mumu", "", 600, 600);
    relative_difference_eff_mumu->SetLeftMargin(0.12);
    gStyle->SetOptStat(0);
    h_relative_difference_eff_ttbar_muon->Scale(100.);
    h_relative_difference_eff_ttbar_muon->SetTitle("");
    h_relative_difference_eff_ttbar_muon->SetXTitle("#Delta R(#mu, jet)");
    h_relative_difference_eff_ttbar_muon->SetYTitle("Relative difference of efficiencies (%)");
    h_relative_difference_eff_ttbar_muon->GetYaxis()->SetTitleOffset(1.5);
    h_relative_difference_eff_ttbar_muon->SetMaximum(101.);
    h_relative_difference_eff_ttbar_muon->SetMinimum(0.);
    h_relative_difference_eff_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_relative_difference_eff_ttbar_muon->SetMarkerColor(kMagenta);
    h_relative_difference_eff_ttbar_muon->SetLineColor(kMagenta);
    h_relative_difference_eff_ttbar_muon->Draw();
    //h_rebin_ttbar_muon->SetMarkerStyle(kOpenDiamond);
    //h_rebin_ttbar_muon->SetMarkerColor(kMagenta+3);
    //h_rebin_ttbar_muon->SetLineColor(kMagenta+3);
    //h_rebin_ttbar_muon->Draw("same");
    h_relative_difference_eff_GG_ttn1_muon->Scale(100.);
    h_relative_difference_eff_GG_ttn1_muon->SetMarkerStyle(kFullCross);
    h_relative_difference_eff_GG_ttn1_muon->SetMarkerColor(kOrange);
    h_relative_difference_eff_GG_ttn1_muon->SetLineColor(kOrange);
    h_relative_difference_eff_GG_ttn1_muon->Draw("same");
    //h_rebin_GG_ttn1_muon->SetMarkerStyle(kOpenCross);
    //h_rebin_GG_ttn1_muon->SetMarkerColor(kOrange+3);
    //h_rebin_GG_ttn1_muon->SetLineColor(kOrange+3);
    //h_rebin_GG_ttn1_muon->Draw("same");
    TPaveText *text6 = new TPaveText(0.6, 0.8, 0.8, 0.85, "NDC");
    text6->AddText(pt_range.c_str());
    text6->SetTextSize(0.03);
    text6->Draw("same");
    TLegend *leg6 = new TLegend(0.6, 0.7, 0.8, 0.8);
    leg6->AddEntry(h_relative_difference_eff_ttbar_muon, "t#bar{t}", "lp");
    leg6->AddEntry(h_relative_difference_eff_GG_ttn1_muon, "#tilde{g}#rightarrow t#bar{t}", "lp");
    //leg6->AddEntry(h_rebin_ttbar_muon, "t#bar{t} (Rebinned)", "lp");
    //leg6->AddEntry(h_rebin_GG_ttn1_muon, "#tilde{g}#rightarrow t#bar{t} (Rebinned)", "lp");
    leg6->SetTextSize(0.03);
    leg6->SetBorderSize(0);
    leg6->SetFillColor(0);
    leg6->SetFillStyle(0);
    leg6->Draw("same");
    string filename_relative_difference_eff_mumu = "real_efficiency_relative_difference_muon_pt" + ss_pt_low.str() + "_" + ss_pt_up.str() + ".pdf";
    relative_difference_eff_mumu->SaveAs(filename_relative_difference_eff_mumu.c_str(), "pdf");
}
