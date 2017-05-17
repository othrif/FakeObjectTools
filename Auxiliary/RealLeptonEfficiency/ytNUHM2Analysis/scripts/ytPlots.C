#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

//#include "AtlasStyle/AtlasStyle.C"
//#include "AtlasStyle/AtlasLabels.C"
//#include "AtlasStyle/AtlasUtils.C"

#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

// TString path = "../../Results/20170222/";
TString path = "/Users/ytshen/Documents/Working/OU/HEP/my_codes/SS_v50/Results/20170310/";

TString electron = "hist-RLE-merged-data-elec.root";
TString muon = "hist-RLE-merged-data-muon.root";

TString electron_tag_trigger_matched = "hist-RLE-merged-data-elec_tag_trigger_matched.root";
TString muon_tag_trigger_matched = "hist-RLE-merged-data-muon_tag_trigger_matched.root";

TString electron_dilepton_trigger = "hist-RLE-merged-data-elec_dilepton_trigger.root";
TString muon_dilepton_trigger = "hist-RLE-merged-data-muon_dilepton_trigger.root";

TString electron_dilepton_trigger_tag_trigger_matched = "hist-RLE-merged-data-elec_dilepton_trigger_tag_trigger_matched.root";
TString muon_dilepton_trigger_tag_trigger_matched = "hist-RLE-merged-data-muon_dilepton_trigger_tag_trigger_matched.root";

TString Zee_root = "RLE_MC_Zee/hist-20170310.root";
TString Zmumu_root = "RLE_MC_Zmumu/hist-20170310.root";

TString Zee_truth_root = "RLE_MC_Zee_truth_match/hist-20170310.root";
TString Zmumu_truth_root = "RLE_MC_Zmumu_truth_match/hist-20170310.root";

TString Zee_TP_truth_root = "RLE_MC_Zee_TandP_truth_match/hist-20170310.root";
TString Zmumu_TP_truth_root = "RLE_MC_Zmumu_TandP_truth_match/hist-20170310.root";

TString ttbar_elec_root = "RLE_MC_ttbar_electron/hist-20170310.root";
TString ttbar_muon_root = "RLE_MC_ttbar_muon/hist-20170310.root";

TString Gtt_elec_root = "RLE_MC_GG_ttn1_electron_20170310/hist-20170310.root";
TString Gtt_muon_root = "RLE_MC_GG_ttn1_muon_20170310/hist-20170310.root";

TString Gtt_elec_all_signal_root = "RLE_MC_GG_ttn1_electron_all_signal/hist-20170118.root";
TString Gtt_elec_Medium_only_root = "RLE_MC_GG_ttn1_electron_medium_only/hist-20170118.root";
TString Gtt_elec_track_iso_only_root = "RLE_MC_GG_ttn1_electron_Track_iso_only/hist-20170118.root";
TString Gtt_elec_calo_iso_only_root = "RLE_MC_GG_ttn1_electron_Calo_iso_only/hist-20170118.root";
TString Gtt_elec_z0_only_root = "RLE_MC_GG_ttn1_electron_z0_cut_only/hist-20170118.root";
TString Gtt_elec_eta_only_root = "RLE_MC_GG_ttn1_electron_eta_cut_only/hist-20170118.root";

TString Zee_tag_trigger_matched_root = "RLE_MC_Zee_tag_trigger_matched/hist-20170310.root";
TString Zee_dilepton_trigger_root = "RLE_MC_Zee_dilepton_trigger/hist-20170310.root";
TString Zee_dilepton_trigger_tag_trigger_matched_root = "RLE_MC_Zee_dilepton_trigger_tag_trigger_matched/hist-20170310.root";

TString Zmumu_tag_trigger_matched_root = "RLE_MC_Zmumu_tag_trigger_matched/hist-20170310.root";
TString Zmumu_dilepton_trigger_root = "RLE_MC_Zmumu_dilepton_trigger/hist-20170310.root";
TString Zmumu_dilepton_trigger_tag_trigger_matched_root = "RLE_MC_Zmumu_dilepton_trigger_tag_trigger_matched/hist-20170310.root";

//
// Mll plots
//

// Make baseline level and signal level mll plots or ratio plots for pt_bin_low < pt < pt_bin_up.
// Using pt_bin_low = 0 and pt_bin_up = -1 means the full pt range.
// Default overlap data, Tag and Probe, and truth match results.
// Set truth_match = true to turn on truth match on the plot.
void yt_mll_plots()
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);

    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TH3F *data_elec_baseline_pt_eta_mll = (TH3F *)data_elec->Get("h_baseline_pt_eta_mll");
    TH3F *data_muon_baseline_pt_eta_mll = (TH3F *)data_muon->Get("h_baseline_pt_eta_mll");
    TH3F *data_elec_signal_pt_eta_mll = (TH3F *)data_elec->Get("h_signal_pt_eta_mll");
    TH3F *data_muon_signal_pt_eta_mll = (TH3F *)data_muon->Get("h_signal_pt_eta_mll");
    TH3F *Zee_TandP_baseline_pt_eta_mll = (TH3F *)Zee_TandP->Get("h_baseline_pt_eta_mll");
    TH3F *Zmumu_TandP_baseline_pt_eta_mll = (TH3F *)Zmumu_TandP->Get("h_baseline_pt_eta_mll");

    // Project all bins to Z axis
    TH1D *data_elec_baseline_mll = (TH1D *)data_elec_baseline_pt_eta_mll->ProjectionZ("data_elec_baseline_mll");
    TH1D *data_muon_baseline_mll = (TH1D *)data_muon_baseline_pt_eta_mll->ProjectionZ("data_muon_baseline_mll");
    TH1D *data_elec_signal_mll = (TH1D *)data_elec_signal_pt_eta_mll->ProjectionZ("data_elec_signal_mll");
    TH1D *data_muon_signal_mll = (TH1D *)data_muon_signal_pt_eta_mll->ProjectionZ("data_muon_signal_mll");
    TH1D *Zee_TandP_mll = (TH1D *)Zee_TandP_baseline_pt_eta_mll->ProjectionZ("Zee_TandP_mll");
    TH1D *Zmumu_TandP_mll = (TH1D *)Zmumu_TandP_baseline_pt_eta_mll->ProjectionZ("Zmumu_TandP_mll");

    double elec_max = max(data_elec_baseline_mll->GetMaximum(), data_elec_signal_mll->GetMaximum());
    double muon_max = max(data_muon_baseline_mll->GetMaximum(), data_muon_signal_mll->GetMaximum());

    // Create a dummy histogram in order to fill color between 80 < mll < 100
    TH1F *h_80_100 = new TH1F("h_80_100", "h_80_100", 1, 80, 100);
    h_80_100->SetFillColor(kYellow);
    h_80_100->SetFillColorAlpha(kYellow, 0.35);
    h_80_100->SetFillStyle(1001); // solid
    h_80_100->SetLineColor(kYellow);

    TCanvas *c_ee = new TCanvas("c_ee", "M_{ee}", 600, 600);
    c_ee->SetGrid();
    c_ee->SetLogy();
    c_ee->SetTickx();
    c_ee->SetTicky();
    c_ee->SetLeftMargin(0.12); //12 per cent of pad for left margin, default is 10%

    data_elec_baseline_mll->SetTitle("");
    data_elec_baseline_mll->SetXTitle("M_{ee} [GeV]");
    data_elec_baseline_mll->SetYTitle("Events / 2 GeV");
    data_elec_baseline_mll->SetYTitle("Number of leptons / 2 GeV / 36.5 fb^{-1}");
    data_elec_baseline_mll->GetYaxis()->SetTitleOffset(1.5);
    data_elec_baseline_mll->SetStats(kFALSE);
    data_elec_baseline_mll->SetMaximum(elec_max * 10);
    data_elec_baseline_mll->SetMinimum(1000.);
    data_elec_baseline_mll->SetLineColor(kBlack);
    data_elec_baseline_mll->SetMarkerColor(kBlack);
    data_elec_baseline_mll->SetMarkerStyle(kFullCircle);
    data_elec_baseline_mll->Rebin(2);
    data_elec_baseline_mll->Draw();

    data_elec_signal_mll->SetLineColor(kRed);
    data_elec_signal_mll->SetMarkerColor(kRed);
    data_elec_signal_mll->SetMarkerStyle(kFullTriangleUp);
    data_elec_signal_mll->Rebin(2);
    data_elec_signal_mll->Draw("same");

    h_80_100->Fill(90, elec_max * 10);
    h_80_100->Draw("hist,same");

    Zee_TandP_mll->SetLineColor(kAzure);
    Zee_TandP_mll->SetFillColor(kAzure);
    //Zee_TandP_mll->SetFillStyle(1001);
    Zee_TandP_mll->SetFillStyle(3004);
    Zee_TandP_mll->Rebin(2);
    Zee_TandP_mll->Draw("hist,same");

    TLegend *leg_ee = new TLegend(0.48, 0.75, 0.9, 0.9);
    //leg_ee->SetHeader("Data 2015 + 2016: 36.5 fb^{-1}", "C");
    leg_ee->SetHeader("Data 2015 + 2016: 36.5 fb^{-1}");
    leg_ee->AddEntry(data_elec_baseline_mll, "pass baseline cut", "p");
    leg_ee->AddEntry(data_elec_signal_mll, "pass signal cut", "p");
    leg_ee->AddEntry(Zee_TandP_mll, "Z #rightarrow ee, pass baseline cut", "l");
    //leg_ee->SetBorderSize(0);
    leg_ee->SetTextFont(42);
    leg_ee->SetTextSize(0.03);
    //leg_ee->SetFillColor(0);
    //leg_ee->SetFillStyle(0);
    leg_ee->Draw("same");

    c_ee->SaveAs("baseline_and_signal_mll_distribution_Mee.pdf", "pdf");

    TCanvas *c_mumu = new TCanvas("c_mumu", "M_{#mu#mu}", 600, 600);
    c_mumu->SetGrid();
    c_mumu->SetLogy();
    c_mumu->SetTickx();
    c_mumu->SetTicky();
    c_mumu->SetLeftMargin(0.12); //12 per cent of pad for left margin, default is 10%

    data_muon_baseline_mll->SetTitle("");
    data_muon_baseline_mll->SetXTitle("M_{#mu#mu} [GeV]");
    data_muon_baseline_mll->SetYTitle("Events / 2 GeV");
    data_muon_baseline_mll->SetYTitle("Number of leptons / 2 GeV / 36.5 fb^{-1}");
    data_muon_baseline_mll->GetYaxis()->SetTitleOffset(1.5);
    data_muon_baseline_mll->SetStats(kFALSE);
    data_muon_baseline_mll->SetMaximum(muon_max * 10);
    data_muon_baseline_mll->SetMinimum(1000.);
    data_muon_baseline_mll->SetLineColor(kBlack);
    data_muon_baseline_mll->SetMarkerColor(kBlack);
    data_muon_baseline_mll->SetMarkerStyle(kFullCircle);
    data_muon_baseline_mll->Rebin(2);
    data_muon_baseline_mll->Draw();

    data_muon_signal_mll->SetLineColor(kRed);
    data_muon_signal_mll->SetMarkerColor(kRed);
    data_muon_signal_mll->SetMarkerStyle(kFullTriangleUp);
    data_muon_signal_mll->Rebin(2);
    data_muon_signal_mll->Draw("same");

    h_80_100->Reset();
    h_80_100->Fill(90, muon_max * 10);
    h_80_100->Draw("hist,same");

    Zmumu_TandP_mll->SetLineColor(kAzure);
    Zmumu_TandP_mll->SetFillColor(kAzure);
    //Zmumu_TandP_mll->SetFillStyle(1001);
    Zmumu_TandP_mll->SetFillStyle(3004);
    Zmumu_TandP_mll->Rebin(2);
    Zmumu_TandP_mll->Draw("hist,same");

    TLegend *leg_mumu = new TLegend(0.48, 0.75, 0.9, 0.9);
    //leg_mumu->SetHeader("Data 2015 + 2016: 36.5 fb^{-1}", "C");
    leg_mumu->SetHeader("Data 2015 + 2016: 36.5 fb^{-1}");
    leg_mumu->AddEntry(data_muon_baseline_mll, "pass baseline cut", "p");
    leg_mumu->AddEntry(data_muon_signal_mll, "pass signal cut", "p");
    leg_mumu->AddEntry(Zmumu_TandP_mll, "Z #rightarrow #mu#mu, pass baseline cut", "l");
    //leg_mumu->SetBorderSize(0);
    leg_mumu->SetTextFont(42);
    leg_mumu->SetTextSize(0.03);
    //leg_mumu->SetFillColor(0);
    //leg_mumu->SetFillStyle(0);
    leg_mumu->Draw("same");

    c_mumu->SaveAs("baseline_and_signal_mll_distribution_Mmumu.pdf", "pdf");
}

void yt_baseline_mll_plots(int pt_bin_low = 0, int pt_bin_up = -1, bool truth_match = false)
{
    // default pt bins are:
    // 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200

    bool debug = false;

    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *Zee_truth;
    TFile *Zmumu_truth;
    if (truth_match) {
        Zee_truth = TFile::Open(path + Zee_truth_root);
        Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    }

    if (debug) {
        cout << data_elec << endl;
        cout << data_muon << endl;
        cout << Zee_TandP << endl;
        cout << Zmumu_TandP << endl;
        if (truth_match) {
            cout << Zee_truth << endl;
            cout << Zmumu_truth << endl;
        }
    }

    TH3F *data_elec_baseline_pt_eta_mll = (TH3F *)data_elec->Get("h_baseline_pt_eta_mll");
    TH3F *data_muon_baseline_pt_eta_mll = (TH3F *)data_muon->Get("h_baseline_pt_eta_mll");
    TH3F *Zee_TandP_baseline_pt_eta_mll = (TH3F *)Zee_TandP->Get("h_baseline_pt_eta_mll");
    TH3F *Zmumu_TandP_baseline_pt_eta_mll = (TH3F *)Zmumu_TandP->Get("h_baseline_pt_eta_mll");

    TH3F *Zee_truth_baseline_pt_eta_mll;
    TH3F *Zmumu_truth_baseline_pt_eta_mll;
    if (truth_match) {
        Zee_truth_baseline_pt_eta_mll = (TH3F *)Zee_truth->Get("h_baseline_pt_eta_mll");
        Zmumu_truth_baseline_pt_eta_mll = (TH3F *)Zmumu_truth->Get("h_baseline_pt_eta_mll");
    }
    if (debug) {
        cout << data_elec_baseline_pt_eta_mll << endl;
        cout << data_muon_baseline_pt_eta_mll << endl;
        cout << Zee_TandP_baseline_pt_eta_mll << endl;
        cout << Zmumu_TandP_baseline_pt_eta_mll << endl;
        if (truth_match) {
            cout << Zee_truth_baseline_pt_eta_mll << endl;
            cout << Zmumu_truth_baseline_pt_eta_mll << endl;
        }
    }

    TH1D *data_elec_mll = (TH1D *)data_elec_baseline_pt_eta_mll->ProjectionZ("data_elec_mll", pt_bin_low, pt_bin_up);
    TH1D *data_muon_mll = (TH1D *)data_muon_baseline_pt_eta_mll->ProjectionZ("data_muon_mll", pt_bin_low, pt_bin_up);
    TH1D *Zee_TandP_mll = (TH1D *)Zee_TandP_baseline_pt_eta_mll->ProjectionZ("Zee_TandP_mll", pt_bin_low, pt_bin_up);
    TH1D *Zmumu_TandP_mll = (TH1D *)Zmumu_TandP_baseline_pt_eta_mll->ProjectionZ("Zmumu_TandP_mll", pt_bin_low, pt_bin_up);

    TH1D *Zee_truth_mll;
    TH1D *Zmumu_truth_mll;
    if (truth_match) {
        Zee_truth_mll = (TH1D *)Zee_truth_baseline_pt_eta_mll->ProjectionZ("Zee_truth_mll", pt_bin_low, pt_bin_up);
        Zmumu_truth_mll = (TH1D *)Zmumu_truth_baseline_pt_eta_mll->ProjectionZ("Zmumu_truth_mll", pt_bin_low, pt_bin_up);
    }
    if (debug) {
        cout << data_elec_mll << endl;
        cout << data_muon_mll << endl;
        cout << Zee_TandP_mll << endl;
        cout << Zmumu_TandP_mll << endl;
        if (truth_match) {
            cout << Zee_truth_mll << endl;
            cout << Zmumu_truth_mll << endl;
        }
    }

    double pt_bin_low_value;
    double pt_bin_up_value;
    if (pt_bin_low != 0 && pt_bin_up != -1) {
        pt_bin_low_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinLowEdge(pt_bin_low);
        pt_bin_up_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinUpEdge(pt_bin_up);
    }
    else {
        pt_bin_low_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinLowEdge(1);
        int nbins = data_elec_baseline_pt_eta_mll->GetXaxis()->GetNbins();
        pt_bin_up_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinUpEdge(nbins - 1);
    }
    // Convert double to string using stringstream
    stringstream sstream_pt_low, sstream_pt_up;
    sstream_pt_low << pt_bin_low_value;
    sstream_pt_up << pt_bin_up_value;

    TCanvas *Mee_plot = new TCanvas("Mee_plot", "Baseline level M_{ee} plot", 500, 500);
    Mee_plot->SetLeftMargin(0.12);

    data_elec_mll->Rebin(2);
    Zee_TandP_mll->Rebin(2);

    double data_elec_value = data_elec_mll->GetMaximum();
    double Zee_TandP_max_value = Zee_TandP_mll->GetMaximum();

    data_elec_mll->SetMarkerColor(kBlack);
    data_elec_mll->SetMarkerStyle(kFullCircle);
    data_elec_mll->SetLineColor(kBlack);
    data_elec_mll->SetTitle("");
    data_elec_mll->SetXTitle("M_{ee} [GeV]");
    data_elec_mll->SetYTitle("Events / 2 GeV");
    data_elec_mll->GetYaxis()->SetTitleOffset(1.5);
    data_elec_mll->SetMaximum(max(data_elec_value, Zee_TandP_max_value) * 1.05);
    data_elec_mll->SetMinimum(0.);
    data_elec_mll->SetStats(kFALSE);
    data_elec_mll->Draw();
    //Zee_TandP_mll->SetMarkerColor(kRed);
    //Zee_TandP_mll->SetMarkerStyle(kFullSquare);
    Zee_TandP_mll->SetLineColor(kRed);
    Zee_TandP_mll->SetFillColor(kRed);
    Zee_TandP_mll->SetFillStyle(3004);
    Zee_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zee_truth_mll->SetMarkerColor(kBlue);
        //Zee_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zee_truth_mll->SetLineColor(kBlue);
        Zee_truth_mll->SetFillColor(kBlue);
        Zee_truth_mll->SetFillStyle(3005);
        Zee_truth_mll->Rebin(2);
        Zee_truth_mll->Draw("hist,same");
    }

    TPaveText *text1 = new TPaveText(0.5, 0.9, 0.9, 0.8, "NDC");
    text1->SetTextSize(0.04);
    text1->SetBorderSize(0);
    text1->SetFillStyle(0);
    text1->SetFillColor(0);
    text1->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text1->Draw("same");

    TLegend *leg1 = new TLegend(0.5, 0.6, 0.9, 0.8);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetFillColor(0);
    leg1->AddEntry(data_elec_mll, "Data");
    leg1->AddEntry(Zee_TandP_mll, "Zee Tag and Probe", "f");
    if (truth_match) {
        leg1->AddEntry(Zee_truth_mll, "Zee truth match", "f");
    }
    leg1->Draw("same");

    string filename1 = "baseline_level_Mee_pt" + sstream_pt_low.str() + sstream_pt_up.str() + ".pdf";
    Mee_plot->SaveAs(filename1.c_str(), "pdf");

    TCanvas *Mmumu_plot = new TCanvas("Mmumu_plot", "Baseline level M_{#mu#mu} plot", 500, 500);
    Mmumu_plot->SetLeftMargin(0.12);

    data_muon_mll->Rebin(2);
    Zmumu_TandP_mll->Rebin(2);

    double data_muon_value = data_muon_mll->GetMaximum();
    double Zmumu_TandP_max_value = Zmumu_TandP_mll->GetMaximum();

    data_muon_mll->SetMarkerColor(kBlack);
    data_muon_mll->SetMarkerStyle(kFullCircle);
    data_muon_mll->SetLineColor(kBlack);
    data_muon_mll->SetTitle("");
    data_muon_mll->SetXTitle("M_{#mu#mu} [GeV]");
    data_muon_mll->SetYTitle("Events / 2 GeV");
    data_muon_mll->GetYaxis()->SetTitleOffset(1.5);
    data_muon_mll->SetMaximum(max(data_muon_value, Zmumu_TandP_max_value) * 1.05);
    data_muon_mll->SetMinimum(0.);
    data_muon_mll->SetStats(kFALSE);
    data_muon_mll->Draw();
    //Zmumu_TandP_mll->SetMarkerColor(kRed);
    //Zmumu_TandP_mll->SetMarkerStyle(kFullSquare);
    Zmumu_TandP_mll->SetLineColor(kRed);
    Zmumu_TandP_mll->SetFillColor(kRed);
    Zmumu_TandP_mll->SetFillStyle(3004);
    Zmumu_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zmumu_truth_mll->SetMarkerColor(kBlue);
        //Zmumu_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zmumu_truth_mll->SetLineColor(kBlue);
        Zmumu_truth_mll->SetFillColor(kBlue);
        Zmumu_truth_mll->SetFillStyle(3005);
        Zmumu_truth_mll->Rebin(2);
        Zmumu_truth_mll->Draw("hist,same");
    }

    TPaveText *text2 = new TPaveText(0.5, 0.9, 0.9, 0.8, "NDC");
    text2->SetTextSize(0.04);
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->SetFillColor(0);
    text2->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text2->Draw("same");

    TLegend *leg2 = new TLegend(0.5, 0.6, 0.9, 0.8);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetFillColor(0);
    leg2->AddEntry(data_muon_mll, "Data");
    leg2->AddEntry(Zmumu_TandP_mll, "Z#mu#mu Tag and Probe", "f");
    if (truth_match) {
        leg2->AddEntry(Zmumu_truth_mll, "Z#mu#mu truth match", "f");
    }
    leg2->Draw("same");

    string filename2 = "baseline_level_Mmumu_pt" + sstream_pt_low.str() + sstream_pt_up.str() + ".pdf";
    Mmumu_plot->SaveAs(filename2.c_str(), "pdf");
}

void yt_baseline_mll_ratio_plots(int pt_bin_low = 0, int pt_bin_up = -1, bool truth_match = false, bool norm = false)
{
    // default pt bins are:
    // 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200

    bool debug = false;

    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *Zee_truth;
    TFile *Zmumu_truth;
    if (truth_match) {
        Zee_truth = TFile::Open(path + Zee_truth_root);
        Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    }

    if (debug) {
        cout << data_elec << endl;
        cout << data_muon << endl;
        cout << Zee_TandP << endl;
        cout << Zmumu_TandP << endl;
        if (truth_match) {
            cout << Zee_truth << endl;
            cout << Zmumu_truth << endl;
        }
    }

    TH3F *data_elec_baseline_pt_eta_mll = (TH3F *)data_elec->Get("h_baseline_pt_eta_mll");
    TH3F *data_muon_baseline_pt_eta_mll = (TH3F *)data_muon->Get("h_baseline_pt_eta_mll");
    TH3F *Zee_TandP_baseline_pt_eta_mll = (TH3F *)Zee_TandP->Get("h_baseline_pt_eta_mll");
    TH3F *Zmumu_TandP_baseline_pt_eta_mll = (TH3F *)Zmumu_TandP->Get("h_baseline_pt_eta_mll");

    TH3F *Zee_truth_baseline_pt_eta_mll;
    TH3F *Zmumu_truth_baseline_pt_eta_mll;
    if (truth_match) {
        Zee_truth_baseline_pt_eta_mll = (TH3F *)Zee_truth->Get("h_baseline_pt_eta_mll");
        Zmumu_truth_baseline_pt_eta_mll = (TH3F *)Zmumu_truth->Get("h_baseline_pt_eta_mll");
    }
    if (debug) {
        cout << data_elec_baseline_pt_eta_mll << endl;
        cout << data_muon_baseline_pt_eta_mll << endl;
        cout << Zee_TandP_baseline_pt_eta_mll << endl;
        cout << Zmumu_TandP_baseline_pt_eta_mll << endl;
        if (truth_match) {
            cout << Zee_truth_baseline_pt_eta_mll << endl;
            cout << Zmumu_truth_baseline_pt_eta_mll << endl;
        }
    }

    TH1D *data_elec_mll = (TH1D *)data_elec_baseline_pt_eta_mll->ProjectionZ("data_elec_mll", pt_bin_low, pt_bin_up)->Clone();
    TH1D *data_muon_mll = (TH1D *)data_muon_baseline_pt_eta_mll->ProjectionZ("data_muon_mll", pt_bin_low, pt_bin_up)->Clone();
    TH1D *Zee_TandP_mll = (TH1D *)Zee_TandP_baseline_pt_eta_mll->ProjectionZ("Zee_TandP_mll", pt_bin_low, pt_bin_up)->Clone();
    TH1D *Zmumu_TandP_mll = (TH1D *)Zmumu_TandP_baseline_pt_eta_mll->ProjectionZ("Zmumu_TandP_mll", pt_bin_low, pt_bin_up)->Clone();

    TH1D *Zee_truth_mll;
    TH1D *Zmumu_truth_mll;
    if (truth_match) {
        Zee_truth_mll = (TH1D *)Zee_truth_baseline_pt_eta_mll->ProjectionZ("Zee_truth_mll", pt_bin_low, pt_bin_up);
        Zmumu_truth_mll = (TH1D *)Zmumu_truth_baseline_pt_eta_mll->ProjectionZ("Zmumu_truth_mll", pt_bin_low, pt_bin_up);
    }
    if (debug) {
        cout << data_elec_mll << endl;
        cout << data_muon_mll << endl;
        cout << Zee_TandP_mll << endl;
        cout << Zmumu_TandP_mll << endl;
        if (truth_match) {
            cout << Zee_truth_mll << endl;
            cout << Zmumu_truth_mll << endl;
        }
    }

    double pt_bin_low_value;
    double pt_bin_up_value;
    if (pt_bin_low != 0 && pt_bin_up != -1) {
        pt_bin_low_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinLowEdge(pt_bin_low);
        pt_bin_up_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinUpEdge(pt_bin_up);
    }
    else {
        pt_bin_low_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinLowEdge(1);
        int nbins = data_elec_baseline_pt_eta_mll->GetXaxis()->GetNbins();
        pt_bin_up_value = data_elec_baseline_pt_eta_mll->GetXaxis()->GetBinUpEdge(nbins - 1);
    }
    // Convert double to string using stringstream
    stringstream sstream_pt_low, sstream_pt_up;
    sstream_pt_low << pt_bin_low_value;
    sstream_pt_up << pt_bin_up_value;

    //
    // Normalize MC to data using a Gaussian fit of the Z peak (85 < mll < 95)
    //
    data_elec_mll->Fit("gaus", "0", "", 85, 95);
    TF1 *func_ee = data_elec_mll->GetFunction("gaus");
    //cout << func_ee->Integral(85., 95.) << endl;  // function integral uses values as arguments.
    //cout << data_elec_mll->Integral(data_elec_mll->GetXaxis()->FindBin(85. + 0.01), data_elec_mll->GetXaxis()->FindBin(95. - 0.01)) << endl; // histogram integral uses bins as argument.
    double data_elec_mll_peak_area = func_ee->Integral(85., 95.);

    data_muon_mll->Fit("gaus", "0", "", 85, 95);
    TF1 *func_mumu = data_muon_mll->GetFunction("gaus");
    //cout << func_mumu->Integral(85., 95.) << endl;  // function integral uses values as arguments.
    //cout << data_muon_mll->Integral(data_muon_mll->GetXaxis()->FindBin(85. + 0.01), data_muon_mll->GetXaxis()->FindBin(95. - 0.01)) << endl; // histogram integral uses bins as argument.
    double data_muon_mll_peak_area = func_mumu->Integral(85., 95.);

    if (norm) {
        double Zee_TandP_mll_peak_area = Zee_TandP_mll->Integral(Zee_TandP_mll->GetXaxis()->FindBin(85. + 0.01), Zee_TandP_mll->GetXaxis()->FindBin(95. - 0.01));
        double Zmumu_TandP_mll_peak_area = Zmumu_TandP_mll->Integral(Zmumu_TandP_mll->GetXaxis()->FindBin(85. + 0.01), Zmumu_TandP_mll->GetXaxis()->FindBin(95. - 0.01));
        Zee_TandP_mll->Scale(data_elec_mll_peak_area / Zee_TandP_mll_peak_area);
        Zmumu_TandP_mll->Scale(data_muon_mll_peak_area / Zmumu_TandP_mll_peak_area);
        if (truth_match) {
            double Zee_truth_mll_peak_area = Zee_truth_mll->Integral(Zee_truth_mll->GetXaxis()->FindBin(85. + 0.01), Zee_truth_mll->GetXaxis()->FindBin(95. - 0.01));
            double Zmumu_truth_mll_peak_area = Zmumu_truth_mll->Integral(Zmumu_truth_mll->GetXaxis()->FindBin(85. + 0.01), Zmumu_truth_mll->GetXaxis()->FindBin(95. - 0.01));
            Zee_truth_mll->Scale(data_elec_mll_peak_area / Zee_truth_mll_peak_area);
            Zmumu_truth_mll->Scale(data_muon_mll_peak_area / Zmumu_truth_mll_peak_area);
        }
    }

    //
    // For Mee
    //
    TCanvas *Mee_plot = new TCanvas("Mee_plot", "Baseline level M_{ee} ratio plot", 600, 600);

    //Upper plot will be in pad1
    TPad *pad1_Mee = new TPad("pad1_Mee", "pad1_Mee", 0, 0.35, 1, 1.0);
    pad1_Mee->SetBottomMargin(0); // Upper and lower plot are joined
    pad1_Mee->SetRightMargin(0.08);
    pad1_Mee->SetLeftMargin(0.12);
    //pad1_Mee->SetGridy(); // grid lines
    //pad1_Mee->SetLogx();
    pad1_Mee->Draw();

    // lower plot will be in pad
    TPad *pad2_Mee = new TPad("pad2_Mee", "pad2_Mee", 0, 0.05, 1, 0.35);
    pad2_Mee->SetTopMargin(0);
    pad2_Mee->SetBottomMargin(0.3);
    pad2_Mee->SetRightMargin(0.08);
    pad2_Mee->SetLeftMargin(0.12);
    pad2_Mee->SetGridy(); // grid lines
    //pad2_Mee->SetLogx();
    pad2_Mee->Draw();

    pad1_Mee->cd(); // pad1 becomes the current pad
    //pad1_Mee->SetFrameLineWidth(2);

    // Draw curve here

    data_elec_mll->Rebin(2);
    Zee_TandP_mll->Rebin(2);
    Zee_truth_mll->Rebin(2);

    double data_elec_value = data_elec_mll->GetMaximum();
    double Zee_TandP_max_value = Zee_TandP_mll->GetMaximum();

    data_elec_mll->SetMarkerColor(kBlack);
    data_elec_mll->SetMarkerStyle(kFullCircle);
    data_elec_mll->SetLineColor(kBlack);
    data_elec_mll->SetTitle("");
    data_elec_mll->SetYTitle("Events / 2 GeV");
    data_elec_mll->GetYaxis()->SetTitleOffset(1.5);
    data_elec_mll->SetMaximum(max(data_elec_value, Zee_TandP_max_value) * 1.05);
    data_elec_mll->SetMinimum(0.1);
    data_elec_mll->SetStats(kFALSE);
    data_elec_mll->Draw();
    //Zee_TandP_mll->SetMarkerColor(kRed);
    //Zee_TandP_mll->SetMarkerStyle(kFullSquare);
    Zee_TandP_mll->SetLineColor(kRed);
    Zee_TandP_mll->SetFillColor(kRed);
    Zee_TandP_mll->SetFillStyle(3004);
    Zee_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zee_truth_mll->SetMarkerColor(kBlue);
        //Zee_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zee_truth_mll->SetLineColor(kBlue);
        Zee_truth_mll->SetFillColor(kBlue);
        Zee_truth_mll->SetFillStyle(3005);
        Zee_truth_mll->Draw("hist,same");
    }

    TPaveText *text1 = new TPaveText(0.5, 0.8, 0.9, 0.89, "NDC");
    text1->SetTextSize(0.05);
    text1->SetBorderSize(0);
    text1->SetFillStyle(0);
    text1->SetFillColor(0);
    text1->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text1->Draw("same");

    TLegend *leg1 = new TLegend(0.55, 0.6, 0.9, 0.8);
    leg1->AddEntry(data_elec_mll, "Data");
    leg1->AddEntry(Zee_TandP_mll, "Zee Tag and Probe", "f");
    if (truth_match) {
        leg1->AddEntry(Zee_truth_mll, "Zee truth match", "f");
    }
    leg1->SetBorderSize(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.05);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");

    pad2_Mee->cd(); // pad2 becomes the current pad

    TH1F *frame_Mee = pad2_Mee->DrawFrame(60, 0.8, 150, 1.19);
    frame_Mee->GetXaxis()->SetNdivisions(510);
    frame_Mee->GetYaxis()->SetNdivisions(405);
    frame_Mee->SetLineWidth(1);
    frame_Mee->SetXTitle("M_{ll} [GeV]");
    frame_Mee->GetXaxis()->SetTitleSize(20);
    frame_Mee->GetXaxis()->SetTitleFont(47);
    frame_Mee->GetXaxis()->SetTitleOffset(3.0);
    frame_Mee->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mee->GetXaxis()->SetLabelSize(15);
    frame_Mee->SetYTitle("Data / MC");
    frame_Mee->GetYaxis()->SetTitleSize(17);
    frame_Mee->GetYaxis()->SetTitleFont(43);
    frame_Mee->GetYaxis()->SetTitleOffset(1.5);
    frame_Mee->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mee->GetYaxis()->SetLabelSize(16);
    frame_Mee->Draw();

    TH1F *ratio_Mee_TandP = (TH1F *)data_elec_mll->Clone();
    ratio_Mee_TandP->Sumw2();
    ratio_Mee_TandP->Divide(Zee_TandP_mll);
    ratio_Mee_TandP->SetMarkerColor(kRed);
    ratio_Mee_TandP->SetMarkerStyle(kFullTriangleUp);
    ratio_Mee_TandP->SetLineColor(kRed);
    ratio_Mee_TandP->Draw("same");

    if (truth_match) {
        TH1F *ratio_Mee_truth = (TH1F *)data_elec_mll->Clone();
        ratio_Mee_truth->Sumw2();
        ratio_Mee_truth->Divide(Zee_truth_mll);
        ratio_Mee_truth->SetMarkerColor(kBlue);
        ratio_Mee_truth->SetMarkerStyle(kFullTriangleDown);
        ratio_Mee_truth->SetLineColor(kBlue);
        ratio_Mee_truth->Draw("same");
    }

    string filename1;
    if (!norm)
        filename1 = "baseline_level_Mee_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot.pdf";
    else
        filename1 = "baseline_level_Mee_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot_MC_normalized.pdf";
    Mee_plot->SaveAs(filename1.c_str(), "pdf");

    //
    // For Mmumu
    //
    TCanvas *Mmumu_plot = new TCanvas("Mmumu_plot", "Baseline level M_{#mu#mu} ratio plot", 600, 600);

    //Upper plot will be in pad1
    TPad *pad1_Mmumu = new TPad("pad1_Mmumu", "pad1_Mmumu", 0, 0.35, 1, 1.0);
    pad1_Mmumu->SetBottomMargin(0); // Upper and lower plot are joined
    pad1_Mmumu->SetRightMargin(0.08);
    pad1_Mmumu->SetLeftMargin(0.12);
    //pad1_Mmumu->SetGridy(); // grid lines
    //pad1_Mmumu->SetLogx();
    pad1_Mmumu->Draw();

    // lower plot will be in pad
    TPad *pad2_Mmumu = new TPad("pad2_Mmumu", "pad2_Mmumu", 0, 0.05, 1, 0.35);
    pad2_Mmumu->SetTopMargin(0);
    pad2_Mmumu->SetBottomMargin(0.3);
    pad2_Mmumu->SetRightMargin(0.08);
    pad2_Mmumu->SetLeftMargin(0.12);
    pad2_Mmumu->SetGridy(); // grid lines
    //pad2_Mmumu->SetLogx();
    pad2_Mmumu->Draw();

    pad1_Mmumu->cd(); // pad1 becomes the current pad
    //pad1_Mmumu->SetFrameLineWidth(2);

    // Draw curve here

    data_muon_mll->Rebin(2);
    Zmumu_TandP_mll->Rebin(2);
    Zmumu_truth_mll->Rebin(2);

    double data_muon_value = data_muon_mll->GetMaximum();
    double Zmumu_TandP_max_value = Zmumu_TandP_mll->GetMaximum();

    data_muon_mll->SetMarkerColor(kBlack);
    data_muon_mll->SetMarkerStyle(kFullCircle);
    data_muon_mll->SetLineColor(kBlack);
    data_muon_mll->SetTitle("");
    data_muon_mll->SetYTitle("Events / 2 GeV");
    data_muon_mll->GetYaxis()->SetTitleOffset(1.5);
    data_muon_mll->SetMaximum(max(data_muon_value, Zmumu_TandP_max_value) * 1.05);
    data_muon_mll->SetMinimum(0.1);
    data_muon_mll->SetStats(kFALSE);
    data_muon_mll->Draw();
    //Zmumu_TandP_mll->SetMarkerColor(kRed);
    //Zmumu_TandP_mll->SetMarkerStyle(kFullSquare);
    Zmumu_TandP_mll->SetLineColor(kRed);
    Zmumu_TandP_mll->SetFillColor(kRed);
    Zmumu_TandP_mll->SetFillStyle(3004);
    Zmumu_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zmumu_truth_mll->SetMarkerColor(kBlue);
        //Zmumu_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zmumu_truth_mll->SetLineColor(kBlue);
        Zmumu_truth_mll->SetFillColor(kBlue);
        Zmumu_truth_mll->SetFillStyle(3005);
        Zmumu_truth_mll->Draw("hist,same");
    }

    TPaveText *text2 = new TPaveText(0.5, 0.8, 0.9, 0.89, "NDC");
    text2->SetTextSize(0.05);
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->SetFillColor(0);
    text2->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text2->Draw("same");

    TLegend *leg2 = new TLegend(0.55, 0.6, 0.9, 0.8);
    leg2->AddEntry(data_muon_mll, "Data");
    leg2->AddEntry(Zmumu_TandP_mll, "Z#mu#mu Tag and Probe", "f");
    if (truth_match) {
        leg2->AddEntry(Zmumu_truth_mll, "Z#mu#mu truth match", "f");
    }
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.05);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->Draw("same");

    pad2_Mmumu->cd(); // pad2 becomes the current pad

    TH1F *frame_Mmumu = pad2_Mmumu->DrawFrame(60, 0.8, 150, 1.19);
    frame_Mmumu->GetXaxis()->SetNdivisions(510);
    frame_Mmumu->GetYaxis()->SetNdivisions(405);
    frame_Mmumu->SetLineWidth(1);
    frame_Mmumu->SetXTitle("M_{ll} [GeV]");
    frame_Mmumu->GetXaxis()->SetTitleSize(20);
    frame_Mmumu->GetXaxis()->SetTitleFont(47);
    frame_Mmumu->GetXaxis()->SetTitleOffset(3.0);
    frame_Mmumu->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mmumu->GetXaxis()->SetLabelSize(15);
    frame_Mmumu->SetYTitle("Data / MC");
    frame_Mmumu->GetYaxis()->SetTitleSize(17);
    frame_Mmumu->GetYaxis()->SetTitleFont(43);
    frame_Mmumu->GetYaxis()->SetTitleOffset(1.5);
    frame_Mmumu->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mmumu->GetYaxis()->SetLabelSize(16);
    frame_Mmumu->Draw();

    TH1F *ratio_Mmumu_TandP = (TH1F *)data_muon_mll->Clone();
    ratio_Mmumu_TandP->Sumw2();
    ratio_Mmumu_TandP->Divide(Zmumu_TandP_mll);
    ratio_Mmumu_TandP->SetMarkerColor(kRed);
    ratio_Mmumu_TandP->SetMarkerStyle(kFullTriangleUp);
    ratio_Mmumu_TandP->SetLineColor(kRed);
    ratio_Mmumu_TandP->Draw("same");

    if (truth_match) {
        TH1F *ratio_Mmumu_truth = (TH1F *)data_muon_mll->Clone();
        ratio_Mmumu_truth->Sumw2();
        ratio_Mmumu_truth->Divide(Zmumu_truth_mll);
        ratio_Mmumu_truth->SetMarkerColor(kBlue);
        ratio_Mmumu_truth->SetMarkerStyle(kFullTriangleDown);
        ratio_Mmumu_truth->SetLineColor(kBlue);
        ratio_Mmumu_truth->Draw("same");
    }

    string filename2;
    if (!norm)
        filename2 = "baseline_level_Mmumu_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot.pdf";
    else
        filename2 = "baseline_level_Mmumu_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot_MC_normalized.pdf";
    Mmumu_plot->SaveAs(filename2.c_str(), "pdf");
}

void yt_signal_mll_plots(int pt_bin_low = 0, int pt_bin_up = -1, bool truth_match = false)
{
    // default pt bins are:
    // 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200

    bool debug = false;

    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *Zee_truth;
    TFile *Zmumu_truth;
    if (truth_match) {
        Zee_truth = TFile::Open(path + Zee_truth_root);
        Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    }

    if (debug) {
        cout << data_elec << endl;
        cout << data_muon << endl;
        cout << Zee_TandP << endl;
        cout << Zmumu_TandP << endl;
        if (truth_match) {
            cout << Zee_truth << endl;
            cout << Zmumu_truth << endl;
        }
    }

    TH3F *data_elec_signal_pt_eta_mll = (TH3F *)data_elec->Get("h_signal_pt_eta_mll");
    TH3F *data_muon_signal_pt_eta_mll = (TH3F *)data_muon->Get("h_signal_pt_eta_mll");
    TH3F *Zee_TandP_signal_pt_eta_mll = (TH3F *)Zee_TandP->Get("h_signal_pt_eta_mll");
    TH3F *Zmumu_TandP_signal_pt_eta_mll = (TH3F *)Zmumu_TandP->Get("h_signal_pt_eta_mll");

    TH3F *Zee_truth_signal_pt_eta_mll;
    TH3F *Zmumu_truth_signal_pt_eta_mll;
    if (truth_match) {
        Zee_truth_signal_pt_eta_mll = (TH3F *)Zee_truth->Get("h_signal_pt_eta_mll");
        Zmumu_truth_signal_pt_eta_mll = (TH3F *)Zmumu_truth->Get("h_signal_pt_eta_mll");
    }
    if (debug) {
        cout << data_elec_signal_pt_eta_mll << endl;
        cout << data_muon_signal_pt_eta_mll << endl;
        cout << Zee_TandP_signal_pt_eta_mll << endl;
        cout << Zmumu_TandP_signal_pt_eta_mll << endl;
        if (truth_match) {
            cout << Zee_truth_signal_pt_eta_mll << endl;
            cout << Zmumu_truth_signal_pt_eta_mll << endl;
        }
    }

    TH1D *data_elec_mll = (TH1D *)data_elec_signal_pt_eta_mll->ProjectionZ("data_elec_mll", pt_bin_low, pt_bin_up);
    TH1D *data_muon_mll = (TH1D *)data_muon_signal_pt_eta_mll->ProjectionZ("data_muon_mll", pt_bin_low, pt_bin_up);
    TH1D *Zee_TandP_mll = (TH1D *)Zee_TandP_signal_pt_eta_mll->ProjectionZ("Zee_TandP_mll", pt_bin_low, pt_bin_up);
    TH1D *Zmumu_TandP_mll = (TH1D *)Zmumu_TandP_signal_pt_eta_mll->ProjectionZ("Zmumu_TandP_mll", pt_bin_low, pt_bin_up);

    TH1D *Zee_truth_mll;
    TH1D *Zmumu_truth_mll;
    if (truth_match) {
        Zee_truth_mll = (TH1D *)Zee_truth_signal_pt_eta_mll->ProjectionZ("Zee_truth_mll", pt_bin_low, pt_bin_up);
        Zmumu_truth_mll = (TH1D *)Zmumu_truth_signal_pt_eta_mll->ProjectionZ("Zmumu_truth_mll", pt_bin_low, pt_bin_up);
    }
    if (debug) {
        cout << data_elec_mll << endl;
        cout << data_muon_mll << endl;
        cout << Zee_TandP_mll << endl;
        cout << Zmumu_TandP_mll << endl;
        if (truth_match) {
            cout << Zee_truth_mll << endl;
            cout << Zmumu_truth_mll << endl;
        }
    }

    double pt_bin_low_value;
    double pt_bin_up_value;
    if (pt_bin_low != 0 && pt_bin_up != -1) {
        pt_bin_low_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinLowEdge(pt_bin_low);
        pt_bin_up_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinUpEdge(pt_bin_up);
    }
    else {
        pt_bin_low_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinLowEdge(1);
        int nbins = data_elec_signal_pt_eta_mll->GetXaxis()->GetNbins();
        pt_bin_up_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinUpEdge(nbins - 1);
    }
    // Convert double to string using stringstream
    stringstream sstream_pt_low, sstream_pt_up;
    sstream_pt_low << pt_bin_low_value;
    sstream_pt_up << pt_bin_up_value;

    TCanvas *Mee_plot = new TCanvas("Mee_plot", "Signal level M_{ee} plot", 500, 500);
    Mee_plot->SetLeftMargin(0.12);

    data_elec_mll->Rebin(2);
    Zee_TandP_mll->Rebin(2);

    double data_elec_value = data_elec_mll->GetMaximum();
    double Zee_TandP_max_value = Zee_TandP_mll->GetMaximum();

    data_elec_mll->SetMarkerColor(kBlack);
    data_elec_mll->SetMarkerStyle(kFullCircle);
    data_elec_mll->SetLineColor(kBlack);
    data_elec_mll->SetTitle("");
    data_elec_mll->SetXTitle("M_{ee} [GeV]");
    data_elec_mll->SetYTitle("Events / 2 GeV");
    data_elec_mll->GetYaxis()->SetTitleOffset(1.5);
    data_elec_mll->SetMaximum(max(data_elec_value, Zee_TandP_max_value) * 1.05);
    data_elec_mll->SetMinimum(0.);
    data_elec_mll->SetStats(kFALSE);
    data_elec_mll->Draw();
    //Zee_TandP_mll->SetMarkerColor(kRed);
    //Zee_TandP_mll->SetMarkerStyle(kFullSquare);
    Zee_TandP_mll->SetLineColor(kRed);
    Zee_TandP_mll->SetFillColor(kRed);
    Zee_TandP_mll->SetFillStyle(3004);
    Zee_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zee_truth_mll->SetMarkerColor(kBlue);
        //Zee_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zee_truth_mll->SetLineColor(kBlue);
        Zee_truth_mll->SetFillColor(kBlue);
        Zee_truth_mll->SetFillStyle(3005);
        Zee_truth_mll->Rebin(2);
        Zee_truth_mll->Draw("hist,same");
    }

    TPaveText *text1 = new TPaveText(0.5, 0.9, 0.9, 0.8, "NDC");
    text1->SetTextSize(0.04);
    text1->SetBorderSize(0);
    text1->SetFillStyle(0);
    text1->SetFillColor(0);
    text1->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text1->Draw("same");

    TLegend *leg1 = new TLegend(0.5, 0.6, 0.9, 0.8);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetFillColor(0);
    leg1->AddEntry(data_elec_mll, "Data");
    leg1->AddEntry(Zee_TandP_mll, "Zee Tag and Probe", "f");
    if (truth_match) {
        leg1->AddEntry(Zee_truth_mll, "Zee truth match", "f");
    }
    leg1->Draw("same");

    string filename1 = "signal_level_Mee_pt" + sstream_pt_low.str() + sstream_pt_up.str() + ".pdf";
    Mee_plot->SaveAs(filename1.c_str(), "pdf");

    TCanvas *Mmumu_plot = new TCanvas("Mmumu_plot", "Signal level M_{#mu#mu} plot", 500, 500);
    Mmumu_plot->SetLeftMargin(0.12);

    data_muon_mll->Rebin(2);
    Zmumu_TandP_mll->Rebin(2);

    double data_muon_value = data_muon_mll->GetMaximum();
    double Zmumu_TandP_max_value = Zmumu_TandP_mll->GetMaximum();

    data_muon_mll->SetMarkerColor(kBlack);
    data_muon_mll->SetMarkerStyle(kFullCircle);
    data_muon_mll->SetLineColor(kBlack);
    data_muon_mll->SetTitle("");
    data_muon_mll->SetXTitle("M_{#mu#mu} [GeV]");
    data_muon_mll->SetYTitle("Events / 1 GeV");
    data_muon_mll->GetYaxis()->SetTitleOffset(1.5);
    data_muon_mll->SetMaximum(max(data_muon_value, Zmumu_TandP_max_value) * 1.05);
    data_muon_mll->SetMinimum(0.);
    data_muon_mll->SetStats(kFALSE);
    data_muon_mll->Draw();
    //Zmumu_TandP_mll->SetMarkerColor(kRed);
    //Zmumu_TandP_mll->SetMarkerStyle(kFullSquare);
    Zmumu_TandP_mll->SetLineColor(kRed);
    Zmumu_TandP_mll->SetFillColor(kRed);
    Zmumu_TandP_mll->SetFillStyle(3004);
    Zmumu_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zmumu_truth_mll->SetMarkerColor(kBlue);
        //Zmumu_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zmumu_truth_mll->SetLineColor(kBlue);
        Zmumu_truth_mll->SetFillColor(kBlue);
        Zmumu_truth_mll->SetFillStyle(3005);
        Zmumu_truth_mll->Rebin(2);
        Zmumu_truth_mll->Draw("hist,same");
    }

    TPaveText *text2 = new TPaveText(0.5, 0.9, 0.9, 0.8, "NDC");
    text2->SetTextSize(0.04);
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->SetFillColor(0);
    text2->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text2->Draw("same");

    TLegend *leg2 = new TLegend(0.5, 0.6, 0.9, 0.8);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetFillColor(0);
    leg2->AddEntry(data_muon_mll, "Data");
    leg2->AddEntry(Zmumu_TandP_mll, "Z#mu#mu Tag and Probe", "f");
    if (truth_match) {
        leg2->AddEntry(Zmumu_truth_mll, "Z#mu#mu truth match", "f");
    }
    leg2->Draw("same");

    string filename2 = "signal_level_Mmumu_pt" + sstream_pt_low.str() + sstream_pt_up.str() + ".pdf";
    Mmumu_plot->SaveAs(filename2.c_str(), "pdf");
}

void yt_signal_mll_ratio_plots(int pt_bin_low = 0, int pt_bin_up = -1, bool truth_match = false, bool norm = false)
{
    // default pt bins are:
    // 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200

    bool debug = false;

    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);
    
    TFile *Zee_truth;
    TFile *Zmumu_truth;
    if (truth_match) {
        Zee_truth = TFile::Open(path + Zee_truth_root);
        Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    }

    if (debug) {
        cout << data_elec << endl;
        cout << data_muon << endl;
        cout << Zee_TandP << endl;
        cout << Zmumu_TandP << endl;
        if (truth_match) {
            cout << Zee_truth << endl;
            cout << Zmumu_truth << endl;
        }
    }

    TH3F *data_elec_signal_pt_eta_mll = (TH3F *)data_elec->Get("h_signal_pt_eta_mll");
    TH3F *data_muon_signal_pt_eta_mll = (TH3F *)data_muon->Get("h_signal_pt_eta_mll");
    TH3F *Zee_TandP_signal_pt_eta_mll = (TH3F *)Zee_TandP->Get("h_signal_pt_eta_mll");
    TH3F *Zmumu_TandP_signal_pt_eta_mll = (TH3F *)Zmumu_TandP->Get("h_signal_pt_eta_mll");

    TH3F *Zee_truth_signal_pt_eta_mll;
    TH3F *Zmumu_truth_signal_pt_eta_mll;
    if (truth_match) {
        Zee_truth_signal_pt_eta_mll = (TH3F *)Zee_truth->Get("h_signal_pt_eta_mll");
        Zmumu_truth_signal_pt_eta_mll = (TH3F *)Zmumu_truth->Get("h_signal_pt_eta_mll");
    }
    if (debug) {
        cout << data_elec_signal_pt_eta_mll << endl;
        cout << data_muon_signal_pt_eta_mll << endl;
        cout << Zee_TandP_signal_pt_eta_mll << endl;
        cout << Zmumu_TandP_signal_pt_eta_mll << endl;
        if (truth_match) {
            cout << Zee_truth_signal_pt_eta_mll << endl;
            cout << Zmumu_truth_signal_pt_eta_mll << endl;
        }
    }

    TH1D *data_elec_mll = (TH1D *)data_elec_signal_pt_eta_mll->ProjectionZ("data_elec_mll", pt_bin_low, pt_bin_up);
    TH1D *data_muon_mll = (TH1D *)data_muon_signal_pt_eta_mll->ProjectionZ("data_muon_mll", pt_bin_low, pt_bin_up);
    TH1D *Zee_TandP_mll = (TH1D *)Zee_TandP_signal_pt_eta_mll->ProjectionZ("Zee_TandP_mll", pt_bin_low, pt_bin_up);
    TH1D *Zmumu_TandP_mll = (TH1D *)Zmumu_TandP_signal_pt_eta_mll->ProjectionZ("Zmumu_TandP_mll", pt_bin_low, pt_bin_up);

    TH1D *Zee_truth_mll;
    TH1D *Zmumu_truth_mll;
    if (truth_match) {
        Zee_truth_mll = (TH1D *)Zee_truth_signal_pt_eta_mll->ProjectionZ("Zee_truth_mll", pt_bin_low, pt_bin_up);
        Zmumu_truth_mll = (TH1D *)Zmumu_truth_signal_pt_eta_mll->ProjectionZ("Zmumu_truth_mll", pt_bin_low, pt_bin_up);
    }
    if (debug) {
        cout << data_elec_mll << endl;
        cout << data_muon_mll << endl;
        cout << Zee_TandP_mll << endl;
        cout << Zmumu_TandP_mll << endl;
        if (truth_match) {
            cout << Zee_truth_mll << endl;
            cout << Zmumu_truth_mll << endl;
        }
    }

    double pt_bin_low_value;
    double pt_bin_up_value;
    if (pt_bin_low != 0 && pt_bin_up != -1) {
        pt_bin_low_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinLowEdge(pt_bin_low);
        pt_bin_up_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinUpEdge(pt_bin_up);
    }
    else {
        pt_bin_low_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinLowEdge(1);
        int nbins = data_elec_signal_pt_eta_mll->GetXaxis()->GetNbins();
        pt_bin_up_value = data_elec_signal_pt_eta_mll->GetXaxis()->GetBinUpEdge(nbins - 1);
    }
    // Convert double to string using stringstream
    stringstream sstream_pt_low, sstream_pt_up;
    sstream_pt_low << pt_bin_low_value;
    sstream_pt_up << pt_bin_up_value;

    //
    // Normalize MC to data using a Gaussian fit of the Z peak (85 < mll < 95)
    //
    data_elec_mll->Fit("gaus", "0", "", 85, 95);
    TF1 *func_ee = data_elec_mll->GetFunction("gaus");
    //cout << func_ee->Integral(85., 95.) << endl;  // function integral uses values as arguments.
    //cout << data_elec_mll->Integral(data_elec_mll->GetXaxis()->FindBin(85. + 0.01), data_elec_mll->GetXaxis()->FindBin(95. - 0.01)) << endl; // histogram integral uses bins as argument.
    double data_elec_mll_peak_area = func_ee->Integral(85., 95.);

    data_muon_mll->Fit("gaus", "0", "", 85, 95);
    TF1 *func_mumu = data_muon_mll->GetFunction("gaus");
    //cout << func_mumu->Integral(85., 95.) << endl;  // function integral uses values as arguments.
    //cout << data_muon_mll->Integral(data_muon_mll->GetXaxis()->FindBin(85. + 0.01), data_muon_mll->GetXaxis()->FindBin(95. - 0.01)) << endl; // histogram integral uses bins as argument.
    double data_muon_mll_peak_area = func_mumu->Integral(85., 95.);

    if (norm) {
        double Zee_TandP_mll_peak_area = Zee_TandP_mll->Integral(Zee_TandP_mll->GetXaxis()->FindBin(85. + 0.01), Zee_TandP_mll->GetXaxis()->FindBin(95. - 0.01));
        double Zmumu_TandP_mll_peak_area = Zmumu_TandP_mll->Integral(Zmumu_TandP_mll->GetXaxis()->FindBin(85. + 0.01), Zmumu_TandP_mll->GetXaxis()->FindBin(95. - 0.01));
        Zee_TandP_mll->Scale(data_elec_mll_peak_area / Zee_TandP_mll_peak_area);
        Zmumu_TandP_mll->Scale(data_muon_mll_peak_area / Zmumu_TandP_mll_peak_area);
        if (truth_match) {
            double Zee_truth_mll_peak_area = Zee_truth_mll->Integral(Zee_truth_mll->GetXaxis()->FindBin(85. + 0.01), Zee_truth_mll->GetXaxis()->FindBin(95. - 0.01));
            double Zmumu_truth_mll_peak_area = Zmumu_truth_mll->Integral(Zmumu_truth_mll->GetXaxis()->FindBin(85. + 0.01), Zmumu_truth_mll->GetXaxis()->FindBin(95. - 0.01));
            Zee_truth_mll->Scale(data_elec_mll_peak_area / Zee_truth_mll_peak_area);
            Zmumu_truth_mll->Scale(data_muon_mll_peak_area / Zmumu_truth_mll_peak_area);
        }
    }

    //
    // For Mee
    //
    TCanvas *Mee_plot = new TCanvas("Mee_plot", "Signal level M_{ee} ratio plot", 600, 600);

    //Upper plot will be in pad1
    TPad *pad1_Mee = new TPad("pad1_Mee", "pad1_Mee", 0, 0.35, 1, 1.0);
    pad1_Mee->SetBottomMargin(0); // Upper and lower plot are joined
    pad1_Mee->SetRightMargin(0.08);
    pad1_Mee->SetLeftMargin(0.12);
    //pad1_Mee->SetGridy(); // grid lines
    //pad1_Mee->SetLogx();
    pad1_Mee->Draw();

    // lower plot will be in pad
    TPad *pad2_Mee = new TPad("pad2_Mee", "pad2_Mee", 0, 0.05, 1, 0.35);
    pad2_Mee->SetTopMargin(0);
    pad2_Mee->SetBottomMargin(0.3);
    pad2_Mee->SetRightMargin(0.08);
    pad2_Mee->SetLeftMargin(0.12);
    pad2_Mee->SetGridy(); // grid lines
    //pad2_Mee->SetLogx();
    pad2_Mee->Draw();

    pad1_Mee->cd(); // pad1 becomes the current pad
    //pad1_Mee->SetFrameLineWidth(2);

    // Draw curve here

    data_elec_mll->Rebin(2);
    Zee_TandP_mll->Rebin(2);
    Zee_truth_mll->Rebin(2);

    double data_elec_value = data_elec_mll->GetMaximum();
    double Zee_TandP_max_value = Zee_TandP_mll->GetMaximum();

    data_elec_mll->SetMarkerColor(kBlack);
    data_elec_mll->SetMarkerStyle(kFullCircle);
    data_elec_mll->SetLineColor(kBlack);
    data_elec_mll->SetTitle("");
    data_elec_mll->SetYTitle("Events / 2 GeV");
    data_elec_mll->GetYaxis()->SetTitleOffset(1.5);
    data_elec_mll->SetMaximum(max(data_elec_value, Zee_TandP_max_value) * 1.05);
    data_elec_mll->SetMinimum(0.1);
    data_elec_mll->SetStats(kFALSE);
    data_elec_mll->Draw();
    //Zee_TandP_mll->SetMarkerColor(kRed);
    //Zee_TandP_mll->SetMarkerStyle(kFullSquare);
    Zee_TandP_mll->SetLineColor(kRed);
    Zee_TandP_mll->SetFillColor(kRed);
    Zee_TandP_mll->SetFillStyle(3004);
    Zee_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zee_truth_mll->SetMarkerColor(kBlue);
        //Zee_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zee_truth_mll->SetLineColor(kBlue);
        Zee_truth_mll->SetFillColor(kBlue);
        Zee_truth_mll->SetFillStyle(3005);
        Zee_truth_mll->Draw("hist,same");
    }

    TPaveText *text1 = new TPaveText(0.5, 0.8, 0.9, 0.89, "NDC");
    text1->SetTextSize(0.05);
    text1->SetBorderSize(0);
    text1->SetFillStyle(0);
    text1->SetFillColor(0);
    text1->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text1->Draw("same");

    TLegend *leg1 = new TLegend(0.55, 0.6, 0.9, 0.8);
    leg1->AddEntry(data_elec_mll, "Data");
    leg1->AddEntry(Zee_TandP_mll, "Zee Tag and Probe", "f");
    if (truth_match) {
        leg1->AddEntry(Zee_truth_mll, "Zee truth match", "f");
    }
    leg1->SetBorderSize(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.05);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");

    pad2_Mee->cd(); // pad2 becomes the current pad

    TH1F *frame_Mee = pad2_Mee->DrawFrame(60, 0.8, 150, 1.19);
    frame_Mee->GetXaxis()->SetNdivisions(510);
    frame_Mee->GetYaxis()->SetNdivisions(405);
    frame_Mee->SetLineWidth(1);
    frame_Mee->SetXTitle("M_{ll} [GeV]");
    frame_Mee->GetXaxis()->SetTitleSize(20);
    frame_Mee->GetXaxis()->SetTitleFont(47);
    frame_Mee->GetXaxis()->SetTitleOffset(3.0);
    frame_Mee->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mee->GetXaxis()->SetLabelSize(15);
    frame_Mee->SetYTitle("Data / MC");
    frame_Mee->GetYaxis()->SetTitleSize(17);
    frame_Mee->GetYaxis()->SetTitleFont(43);
    frame_Mee->GetYaxis()->SetTitleOffset(1.5);
    frame_Mee->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mee->GetYaxis()->SetLabelSize(16);
    frame_Mee->Draw();

    TH1F *ratio_Mee_TandP = (TH1F *)data_elec_mll->Clone();
    ratio_Mee_TandP->Sumw2();
    ratio_Mee_TandP->Divide(Zee_TandP_mll);
    ratio_Mee_TandP->SetMarkerColor(kRed);
    ratio_Mee_TandP->SetMarkerStyle(kFullTriangleUp);
    ratio_Mee_TandP->SetLineColor(kRed);
    ratio_Mee_TandP->Draw("same");

    if (truth_match) {
        TH1F *ratio_Mee_truth = (TH1F *)data_elec_mll->Clone();
        ratio_Mee_truth->Sumw2();
        ratio_Mee_truth->Divide(Zee_truth_mll);
        ratio_Mee_truth->SetMarkerColor(kBlue);
        ratio_Mee_truth->SetMarkerStyle(kFullTriangleDown);
        ratio_Mee_truth->SetLineColor(kBlue);
        ratio_Mee_truth->Draw("same");
    }

    string filename1;
    if (!norm)
        filename1 = "signal_level_Mee_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot.pdf";
    else
        filename1 = "signal_level_Mee_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot_MC_normalized.pdf";
    Mee_plot->SaveAs(filename1.c_str(), "pdf");

    //
    // For Mmumu
    //
    TCanvas *Mmumu_plot = new TCanvas("Mmumu_plot", "Signal level M_{#mu#mu} ratio plot", 600, 600);

    //Upper plot will be in pad1
    TPad *pad1_Mmumu = new TPad("pad1_Mmumu", "pad1_Mmumu", 0, 0.35, 1, 1.0);
    pad1_Mmumu->SetBottomMargin(0); // Upper and lower plot are joined
    pad1_Mmumu->SetRightMargin(0.08);
    pad1_Mmumu->SetLeftMargin(0.12);
    //pad1_Mmumu->SetGridy(); // grid lines
    //pad1_Mmumu->SetLogx();
    pad1_Mmumu->Draw();

    // lower plot will be in pad
    TPad *pad2_Mmumu = new TPad("pad2_Mmumu", "pad2_Mmumu", 0, 0.05, 1, 0.35);
    pad2_Mmumu->SetTopMargin(0);
    pad2_Mmumu->SetBottomMargin(0.3);
    pad2_Mmumu->SetRightMargin(0.08);
    pad2_Mmumu->SetLeftMargin(0.12);
    pad2_Mmumu->SetGridy(); // grid lines
    //pad2_Mmumu->SetLogx();
    pad2_Mmumu->Draw();

    pad1_Mmumu->cd(); // pad1 becomes the current pad
    //pad1_Mmumu->SetFrameLineWidth(2);

    // Draw curve here

    data_muon_mll->Rebin(2);
    Zmumu_TandP_mll->Rebin(2);
    Zmumu_truth_mll->Rebin(2);

    double data_muon_value = data_muon_mll->GetMaximum();
    double Zmumu_TandP_max_value = Zmumu_TandP_mll->GetMaximum();

    data_muon_mll->SetMarkerColor(kBlack);
    data_muon_mll->SetMarkerStyle(kFullCircle);
    data_muon_mll->SetLineColor(kBlack);
    data_muon_mll->SetTitle("");
    data_muon_mll->SetYTitle("Events / 2 GeV");
    data_muon_mll->GetYaxis()->SetTitleOffset(1.5);
    data_muon_mll->SetMaximum(max(data_muon_value, Zmumu_TandP_max_value) * 1.05);
    data_muon_mll->SetMinimum(0.1);
    data_muon_mll->SetStats(kFALSE);
    data_muon_mll->Draw();
    //Zmumu_TandP_mll->SetMarkerColor(kRed);
    //Zmumu_TandP_mll->SetMarkerStyle(kFullSquare);
    Zmumu_TandP_mll->SetLineColor(kRed);
    Zmumu_TandP_mll->SetFillColor(kRed);
    Zmumu_TandP_mll->SetFillStyle(3004);
    Zmumu_TandP_mll->Draw("hist,same");

    if (truth_match) {
        //Zmumu_truth_mll->SetMarkerColor(kBlue);
        //Zmumu_truth_mll->SetMarkerStyle(kFullTriangleUp);
        Zmumu_truth_mll->SetLineColor(kBlue);
        Zmumu_truth_mll->SetFillColor(kBlue);
        Zmumu_truth_mll->SetFillStyle(3005);
        Zmumu_truth_mll->Draw("hist,same");
    }

    TPaveText *text2 = new TPaveText(0.5, 0.8, 0.9, 0.89, "NDC");
    text2->SetTextSize(0.05);
    text2->SetBorderSize(0);
    text2->SetFillStyle(0);
    text2->SetFillColor(0);
    text2->AddText(string(sstream_pt_low.str() + " GeV < p_{T} <" + sstream_pt_up.str() + " GeV").c_str());
    text2->Draw("same");

    TLegend *leg2 = new TLegend(0.55, 0.6, 0.9, 0.8);
    leg2->AddEntry(data_muon_mll, "Data");
    leg2->AddEntry(Zmumu_TandP_mll, "Z#mu#mu Tag and Probe", "f");
    if (truth_match) {
        leg2->AddEntry(Zmumu_truth_mll, "Z#mu#mu truth match", "f");
    }
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.05);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->Draw("same");

    pad2_Mmumu->cd(); // pad2 becomes the current pad

    TH1F *frame_Mmumu = pad2_Mmumu->DrawFrame(60, 0.8, 150, 1.19);
    frame_Mmumu->GetXaxis()->SetNdivisions(510);
    frame_Mmumu->GetYaxis()->SetNdivisions(405);
    frame_Mmumu->SetLineWidth(1);
    frame_Mmumu->SetXTitle("M_{ll} [GeV]");
    frame_Mmumu->GetXaxis()->SetTitleSize(20);
    frame_Mmumu->GetXaxis()->SetTitleFont(47);
    frame_Mmumu->GetXaxis()->SetTitleOffset(3.0);
    frame_Mmumu->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mmumu->GetXaxis()->SetLabelSize(15);
    frame_Mmumu->SetYTitle("Data / MC");
    frame_Mmumu->GetYaxis()->SetTitleSize(17);
    frame_Mmumu->GetYaxis()->SetTitleFont(43);
    frame_Mmumu->GetYaxis()->SetTitleOffset(1.5);
    frame_Mmumu->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_Mmumu->GetYaxis()->SetLabelSize(16);
    frame_Mmumu->Draw();

    TH1F *ratio_Mmumu_TandP = (TH1F *)data_muon_mll->Clone();
    ratio_Mmumu_TandP->Sumw2();
    ratio_Mmumu_TandP->Divide(Zmumu_TandP_mll);
    ratio_Mmumu_TandP->SetMarkerColor(kRed);
    ratio_Mmumu_TandP->SetMarkerStyle(kFullTriangleUp);
    ratio_Mmumu_TandP->SetLineColor(kRed);
    ratio_Mmumu_TandP->Draw("same");

    if (truth_match) {
        TH1F *ratio_Mmumu_truth = (TH1F *)data_muon_mll->Clone();
        ratio_Mmumu_truth->Sumw2();
        ratio_Mmumu_truth->Divide(Zmumu_truth_mll);
        ratio_Mmumu_truth->SetMarkerColor(kBlue);
        ratio_Mmumu_truth->SetMarkerStyle(kFullTriangleDown);
        ratio_Mmumu_truth->SetLineColor(kBlue);
        ratio_Mmumu_truth->Draw("same");
    }

    string filename2;
    if (!norm)
        filename2 = "signal_level_Mmumu_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot.pdf";
    else
        filename2 = "signal_level_Mmumu_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_ratio_plot_MC_normalized.pdf";
    Mmumu_plot->SaveAs(filename2.c_str(), "pdf");
}



//
// Efficiency plots
//

void yt_make_cut_efficiency_plot(TString filename, TString lepton) // Use 1-dim histograms
{
    TFile *file = TFile::Open(filename);
/*
    // Only work for v44 ntuple
    TH1F *h_eff_pt = (TH1F *)file->Get("h_eff_pt");

    TH1F *h_cut_eff_LooseAndBLayerLLH_to_MediumLLH;
    TH1F *h_cut_eff_CaloIso;
    if (lepton == "electron") {
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH = (TH1F *)file->Get("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH");
        h_cut_eff_CaloIso = (TH1F *)file->Get("h_cut_eff_CaloIso");
    } 
    TH1F *h_cut_eff_TrackIso = (TH1F *)file->Get("h_cut_eff_TrackIso");
    TH1F *h_cut_eff_z0 = (TH1F *)file->Get("h_cut_eff_z0");
    TH1F *h_cut_eff_sigd0;
    if (lepton == "muon") {
        h_cut_eff_sigd0 = (TH1F *)file->Get("h_cut_eff_sigd0");
    }
*/
    TH1F *h_baseline_pt = (TH1F *)file->Get("h_baseline_pt");
    TH1F *h_signal_pt = (TH1F *)file->Get("h_signal_pt");
    TH1F *h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt;
    TH1F *h_baseline_CaloIso_cut_pt;
    if (lepton == "electron") {
        h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt = (TH1F *)file->Get("h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt");
        h_baseline_CaloIso_cut_pt = (TH1F *)file->Get("h_baseline_CaloIso_cut_pt");
    }
    TH1F *h_baseline_TrackIso_cut_pt = (TH1F *)file->Get("h_baseline_TrackIso_cut_pt");
    TH1F *h_baseline_z0_cut_pt = (TH1F *)file->Get("h_baseline_z0_cut_pt");
    TH1F *h_baseline_sigd0_cut_pt;
    if (lepton == "muon") {
        h_baseline_sigd0_cut_pt = (TH1F *)file->Get("h_baseline_sigd0_cut_pt");
    }

    // Get the histogram structures
    TH1F *h_eff_pt = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_LooseAndBLayerLLH_to_MediumLLH = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_CaloIso = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_TrackIso = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_z0 = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_sigd0 = (TH1F *)h_baseline_pt->Clone();

    // Reset the contents
    h_eff_pt->Reset();
    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Reset();
    h_cut_eff_CaloIso->Reset();
    h_cut_eff_TrackIso->Reset();
    h_cut_eff_z0->Reset();
    h_cut_eff_sigd0->Reset();

    h_eff_pt->SetName("h_eff_pt");
    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetName("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH");
    h_cut_eff_CaloIso->SetName("h_cut_eff_CaloIso");
    h_cut_eff_TrackIso->SetName("h_cut_eff_TrackIso");
    h_cut_eff_z0->SetName("h_cut_eff_z0");
    h_cut_eff_sigd0->SetName("h_cut_eff_sigd0");

    // Calculate the cut efficiency
    // option="B": Binomial errors are computed
    h_eff_pt->Divide(h_signal_pt, h_baseline_pt, 1, 1, "B");
    if (lepton == "electron") {
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Divide(h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt, h_baseline_pt,1, 1, "B");
        h_cut_eff_CaloIso->Divide(h_baseline_CaloIso_cut_pt, h_baseline_pt,1, 1, "B");
    }
    h_cut_eff_TrackIso->Divide(h_baseline_TrackIso_cut_pt, h_baseline_pt,1, 1, "B");
    h_cut_eff_z0->Divide(h_baseline_z0_cut_pt, h_baseline_pt,1, 1, "B");
    if (lepton == "muon") {
        h_cut_eff_sigd0->Divide(h_baseline_sigd0_cut_pt, h_baseline_pt,1, 1, "B");
    }

    //
    // The following part are the same as ytRealLeptonEfficiency_Data::make_cut_efficiency_plot()
    //

    TCanvas *cut_efficiency_plot = new TCanvas("cut_efficiency", "Cut efficiency comparison", 500, 500);
    cut_efficiency_plot->SetLeftMargin(0.12);
    cut_efficiency_plot->cd();
    //gStyle->SetOptStat(0);

    h_eff_pt->GetYaxis()->SetTitleOffset(1.5);
    h_eff_pt->SetLineColor(kBlack);
    h_eff_pt->SetMarkerColor(kBlack);
    h_eff_pt->SetMarkerStyle(kFullCircle);
    h_eff_pt->SetStats(kFALSE);
    h_eff_pt->SetMaximum(1.1);
    h_eff_pt->SetMinimum(0.4);
    h_eff_pt->SetTitle("");
    h_eff_pt->Draw();

    if (lepton == "electron") {
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetLineColor(kRed);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetMarkerColor(kRed);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetMarkerStyle(kFullSquare);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Draw("same");

        h_cut_eff_CaloIso->SetLineColor(kGreen);
        h_cut_eff_CaloIso->SetMarkerColor(kGreen);
        h_cut_eff_CaloIso->SetMarkerStyle(kFullTriangleUp);
        h_cut_eff_CaloIso->Draw("same");
    }

    h_cut_eff_TrackIso->SetLineColor(kBlue);
    h_cut_eff_TrackIso->SetMarkerColor(kBlue);
    h_cut_eff_TrackIso->SetMarkerStyle(kFullTriangleDown);
    h_cut_eff_TrackIso->Draw("same");

    h_cut_eff_z0->SetLineColor(kMagenta);
    h_cut_eff_z0->SetMarkerColor(kMagenta);
    h_cut_eff_z0->SetMarkerStyle(kFullDiamond);
    h_cut_eff_z0->Draw("same");

    if (lepton == "muon") {
        h_cut_eff_sigd0->SetLineColor(kCyan);
        h_cut_eff_sigd0->SetMarkerColor(kCyan);
        h_cut_eff_sigd0->SetMarkerStyle(kFullCross);
        h_cut_eff_sigd0->Draw("same");
    }

    float X_max = h_eff_pt->GetXaxis()->GetXmax();
    float X_min = h_eff_pt->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.45, 0.2, 0.8, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry("h_eff_pt", "Real Efficiency");
    if (lepton == "electron") {
        leg->AddEntry("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH", "LooseLH to MediumLH");
        leg->AddEntry("h_cut_eff_CaloIso", "Calo isolation");
    }
    leg->AddEntry("h_cut_eff_TrackIso", "Track isolation");
    leg->AddEntry("h_cut_eff_z0", "|z_{0} sin(#theta)| < 0.5");
    if (lepton == "muon") {
        leg->AddEntry("h_cut_eff_sigd0", "|d_{0}/#sigma_{d0}| < 3");
    }
    leg->Draw();

    string output_filename = "cut_efficiency_" + static_cast<string>(lepton) + ".pdf";
    cut_efficiency_plot->SaveAs(output_filename.c_str(), "pdf");
}

void yt_make_cut_efficiency_plot_2(TString filename, TString lepton) // Use 3-dim histograms
{
    TFile *file = TFile::Open(filename);

    TH1F *h_baseline_pt = (TH1F *)file->Get("h_baseline_pt");
    TH1F *h_signal_pt = (TH1F *)file->Get("h_signal_pt");

    TH3F *h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll; // electron only
    TH3F *h_baseline_CaloIso_cut_pt_eta_mll; // electron only
    if (lepton == "electron") {
        h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll = (TH3F *)file->Get("h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll");
        h_baseline_CaloIso_cut_pt_eta_mll = (TH3F *)file->Get("h_baseline_CaloIso_cut_pt_eta_mll");
    }
    TH3F *h_baseline_TrackIso_cut_pt_eta_mll = (TH3F *)file->Get("h_baseline_TrackIso_cut_pt_eta_mll");
    TH3F *h_baseline_z0_cut_pt_eta_mll = (TH3F *)file->Get("h_baseline_z0_cut_pt_eta_mll");
    TH3F *h_baseline_sigd0_cut_pt_eta_mll; // muon only
    if (lepton == "muon") {
        h_baseline_sigd0_cut_pt_eta_mll = (TH3F *)file->Get("h_baseline_sigd0_cut_pt_eta_mll");
    }

    TH1D *h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt;
    TH1D *h_baseline_CaloIso_cut_pt;
    TH1D *h_baseline_TrackIso_cut_pt;
    TH1D *h_baseline_z0_cut_pt;
    TH1D *h_baseline_sigd0_cut_pt;

    if (lepton == "electron") {
        // For electron case, we only use |eta| < 2.0 to calculate the efficiency
        int eta_low_bin = h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll->GetYaxis()->FindBin(0+0.01); // should be 1
        int eta_up_bin = h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll->GetYaxis()->FindBin(2-0.01); // should be 8

        h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt = (TH1D *)h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll->ProjectionX("", eta_low_bin, eta_up_bin)->Clone();
        h_baseline_CaloIso_cut_pt = (TH1D *)h_baseline_CaloIso_cut_pt_eta_mll->ProjectionX("", eta_low_bin, eta_up_bin)->Clone();
        h_baseline_TrackIso_cut_pt = (TH1D *)h_baseline_TrackIso_cut_pt_eta_mll->ProjectionX("", eta_low_bin, eta_up_bin)->Clone();
        h_baseline_z0_cut_pt = (TH1D *)h_baseline_z0_cut_pt_eta_mll->ProjectionX("", eta_low_bin, eta_up_bin)->Clone();
        h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt->SetName("h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt");
        h_baseline_CaloIso_cut_pt->SetName("h_baseline_CaloIso_cut_pt");
        h_baseline_TrackIso_cut_pt->SetName("h_baseline_TrackIso_cut_pt");
        h_baseline_z0_cut_pt->SetName("h_baseline_z0_cut_pt");
    }
    if (lepton == "muon") {
        h_baseline_TrackIso_cut_pt = (TH1D *)h_baseline_TrackIso_cut_pt_eta_mll->ProjectionX("")->Clone();
        h_baseline_z0_cut_pt = (TH1D *)h_baseline_z0_cut_pt_eta_mll->ProjectionX("")->Clone();
        h_baseline_sigd0_cut_pt = (TH1D *)h_baseline_sigd0_cut_pt_eta_mll->ProjectionX("")->Clone();
        h_baseline_TrackIso_cut_pt->SetName("h_baseline_TrackIso_cut_pt");
        h_baseline_z0_cut_pt->SetName("h_baseline_z0_cut_pt");
        h_baseline_sigd0_cut_pt->SetName("h_baseline_sigd0_cut_pt");
    }

    // Get the histogram structures
    TH1F *h_eff_pt = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_LooseAndBLayerLLH_to_MediumLLH = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_CaloIso = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_TrackIso = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_z0 = (TH1F *)h_baseline_pt->Clone();
    TH1F *h_cut_eff_sigd0 = (TH1F *)h_baseline_pt->Clone();

    // Reset the contents
    h_eff_pt->Reset();
    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Reset();
    h_cut_eff_CaloIso->Reset();
    h_cut_eff_TrackIso->Reset();
    h_cut_eff_z0->Reset();
    h_cut_eff_sigd0->Reset();

    h_eff_pt->SetName("h_eff_pt");
    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetName("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH");
    h_cut_eff_CaloIso->SetName("h_cut_eff_CaloIso");
    h_cut_eff_TrackIso->SetName("h_cut_eff_TrackIso");
    h_cut_eff_z0->SetName("h_cut_eff_z0");
    h_cut_eff_sigd0->SetName("h_cut_eff_sigd0");

    // Calculate the cut efficiency
    // option="B": Binomial errors are computed
    h_eff_pt->Divide(h_signal_pt, h_baseline_pt, 1, 1, "B");
    if (lepton == "electron") {
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Divide(h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt, h_baseline_pt,1, 1, "B");
        h_cut_eff_CaloIso->Divide(h_baseline_CaloIso_cut_pt, h_baseline_pt,1, 1, "B");
    }
    h_cut_eff_TrackIso->Divide(h_baseline_TrackIso_cut_pt, h_baseline_pt,1, 1, "B");
    h_cut_eff_z0->Divide(h_baseline_z0_cut_pt, h_baseline_pt,1, 1, "B");
    if (lepton == "muon") {
        h_cut_eff_sigd0->Divide(h_baseline_sigd0_cut_pt, h_baseline_pt,1, 1, "B");
    }

    //
    // The following part are the same as ytRealLeptonEfficiency_Data::make_cut_efficiency_plot()
    //

    TCanvas *cut_efficiency_plot = new TCanvas("cut_efficiency", "Cut efficiency comparison", 500, 500);
    cut_efficiency_plot->SetLeftMargin(0.12);
    cut_efficiency_plot->cd();
    //gStyle->SetOptStat(0);

    h_eff_pt->SetLineColor(kBlack);
    h_eff_pt->SetMarkerColor(kBlack);
    h_eff_pt->SetMarkerStyle(kFullCircle);
    h_eff_pt->SetStats(kFALSE);
    h_eff_pt->SetMaximum(1.1);
    h_eff_pt->SetMinimum(0.4);
    h_eff_pt->SetTitle("");
    h_eff_pt->GetYaxis()->SetTitleOffset(1.5);
    h_eff_pt->Draw();

    if (lepton == "electron") {
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetLineColor(kRed);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetMarkerColor(kRed);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetMarkerStyle(kFullSquare);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Draw("same");

        h_cut_eff_CaloIso->SetLineColor(kGreen);
        h_cut_eff_CaloIso->SetMarkerColor(kGreen);
        h_cut_eff_CaloIso->SetMarkerStyle(kFullTriangleUp);
        h_cut_eff_CaloIso->Draw("same");
    }

    h_cut_eff_TrackIso->SetLineColor(kBlue);
    h_cut_eff_TrackIso->SetMarkerColor(kBlue);
    h_cut_eff_TrackIso->SetMarkerStyle(kFullTriangleDown);
    h_cut_eff_TrackIso->Draw("same");

    h_cut_eff_z0->SetLineColor(kMagenta);
    h_cut_eff_z0->SetMarkerColor(kMagenta);
    h_cut_eff_z0->SetMarkerStyle(kFullDiamond);
    h_cut_eff_z0->Draw("same");

    if (lepton == "muon") {
        h_cut_eff_sigd0->SetLineColor(kCyan);
        h_cut_eff_sigd0->SetMarkerColor(kCyan);
        h_cut_eff_sigd0->SetMarkerStyle(kFullCross);
        h_cut_eff_sigd0->Draw("same");
    }

    float X_max = h_eff_pt->GetXaxis()->GetXmax();
    float X_min = h_eff_pt->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.45, 0.2, 0.8, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry("h_eff_pt", "Real Efficiency");
    if (lepton == "electron") {
        leg->AddEntry("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH", "LooseLH to MediumLH");
        leg->AddEntry("h_cut_eff_CaloIso", "Calo isolation");
    }
    leg->AddEntry("h_cut_eff_TrackIso", "Track isolation");
    leg->AddEntry("h_cut_eff_z0", "|z_{0} sin(#theta)| < 0.5");
    if (lepton == "muon") {
        leg->AddEntry("h_cut_eff_sigd0", "|d_{0}/#sigma_{d0}| < 3");
    }
    leg->Draw();

    string output_filename = "cut_efficiency_" + static_cast<string>(lepton) + ".pdf";
    cut_efficiency_plot->SaveAs(output_filename.c_str(), "pdf");
}

void yt_make_electron_real_efficiency_plot() // Use 2-dim histograms
{
    // TFile *file = TFile::Open(path + electron);
    TFile *file = TFile::Open(path + electron_tag_trigger_matched);
    TH2F *h_baseline_pt_eta = (TH2F *)file->Get("h_baseline_pt_eta");
    TH2F *h_signal_pt_eta = (TH2F *)file->Get("h_signal_pt_eta");

    TH1F *h_baseline_eta_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_1", 1, 3); // 0 < |eta| < 0.8
    TH1F *h_baseline_eta_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_2", 4, 5); // 0.8 < |eta| < 1.37
    TH1F *h_baseline_eta_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_3", 7, 8); // 1.52 < |eta| < 2

    TH1F *h_signal_eta_range_1 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_1", 1, 3); // 0 < |eta| < 0.8
    TH1F *h_signal_eta_range_2 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_2", 4, 5); // 0.8 < |eta| < 1.37
    TH1F *h_signal_eta_range_3 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_3", 7, 8); // 1.52 < |eta| < 2

    //
    // The following part are the same as ytRealLeptonEfficiency_Data::make_electron_real_efficiency_plot()
    //
/*
    // For electrons, there are 3 eta ranges: [0, 0.8], [0.8, 1.37], [1.52, 2]
    // Corresponding to eta bins: [1, 3], [4, 5], and [7, 8]
    // In order to get the correct bin, we have to use lower_edge + 0.01 and upper_edge - 0.01 to find the bin.
    // This is because FindBin() sometimes assign the lower_edge (upper_edge) to the correct_bin-1 (correct_bin+1).
    TH1F *h_baseline_eta_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_1", h_baseline_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(0.8 - 0.01));
    TH1F *h_baseline_eta_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_2", h_baseline_pt_eta->GetYaxis()->FindBin(0.8 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(1.37 - 0.01));
    TH1F *h_baseline_eta_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_3", h_baseline_pt_eta->GetYaxis()->FindBin(1.52 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(2.0 - 0.01));

    TH1F *h_signal_eta_range_1 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_1", h_signal_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(0.8 - 0.01));
    TH1F *h_signal_eta_range_2 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_2", h_signal_pt_eta->GetYaxis()->FindBin(0.8 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(1.37 - 0.01));
    TH1F *h_signal_eta_range_3 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_3", h_signal_pt_eta->GetYaxis()->FindBin(1.52 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(2.0 - 0.01));
*/
    TH1F *ratio_1 = (TH1F *)h_signal_eta_range_1->Clone();
    TH1F *ratio_2 = (TH1F *)h_signal_eta_range_2->Clone();
    TH1F *ratio_3 = (TH1F *)h_signal_eta_range_3->Clone();

    ratio_1->Reset();
    ratio_2->Reset();
    ratio_3->Reset();

    ratio_1->SetName("ratio_1");
    ratio_2->SetName("ratio_2");
    ratio_3->SetName("ratio_3");

    // ratio_1->Divide(h_baseline_eta_range_1);
    // ratio_2->Divide(h_baseline_eta_range_2);
    // ratio_3->Divide(h_baseline_eta_range_3);

    ratio_1->Divide(h_signal_eta_range_1, h_baseline_eta_range_1, 1, 1, "B");
    ratio_2->Divide(h_signal_eta_range_2, h_baseline_eta_range_2, 1, 1, "B");
    ratio_3->Divide(h_signal_eta_range_3, h_baseline_eta_range_3, 1, 1, "B");

    TCanvas *real_efficiency_plot = new TCanvas("real_efficiency_plot", "Electron Real Efficiency", 500, 500);
    real_efficiency_plot->SetLeftMargin(0.12);
    real_efficiency_plot->cd();

    ratio_1->SetLineColor(kBlack);
    ratio_1->SetMarkerColor(kBlack);
    ratio_1->SetMarkerStyle(kFullCircle);
    ratio_1->SetXTitle("p_{T} [GeV]");
    ratio_1->SetYTitle("Electron Real Efficiency");
    ratio_1->GetYaxis()->SetTitleOffset(1.5);
    ratio_1->SetStats(kFALSE);
    ratio_1->SetMaximum(1.1);
    ratio_1->SetMinimum(0.5);
    ratio_1->SetTitle("");
    ratio_1->Draw();

    ratio_2->SetLineColor(kRed);
    ratio_2->SetMarkerColor(kRed);
    ratio_2->SetMarkerStyle(kFullSquare);
    ratio_2->Draw("same");

    ratio_3->SetLineColor(kBlue);
    ratio_3->SetMarkerColor(kBlue);
    ratio_3->SetMarkerStyle(kFullTriangleUp);
    ratio_3->Draw("same");

    float X_max = ratio_1->GetXaxis()->GetXmax();
    float X_min = ratio_1->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.4, 0.2, 0.8, 0.5);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry("ratio_1", "0 < |#eta| < 0.8");
    leg->AddEntry("ratio_2", "0.8 < |#eta| < 1.37");
    leg->AddEntry("ratio_3", "1.52 < |#eta| < 2.0");
    leg->Draw();

    real_efficiency_plot->SaveAs("real_efficiency_electron_before_bkg_subtraction.pdf", "pdf");
}

void yt_make_muon_real_efficiency_plot() // Use 2-dim histograms
{
    // TFile *file = TFile::Open(path + muon);
    TFile *file = TFile::Open(path + muon_tag_trigger_matched);
    TH2F *h_baseline_pt_eta = (TH2F *)file->Get("h_baseline_pt_eta");
    TH2F *h_signal_pt_eta = (TH2F *)file->Get("h_signal_pt_eta");

    TH1F *h_baseline_eta_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_1", 1, 1); // 0 < |eta| < 0.6
    TH1F *h_baseline_eta_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_2", 2, 2); // 0.6 < |eta| < 1.2
    TH1F *h_baseline_eta_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_3", 3, 3); // 1.2 < |eta| < 1.8
    TH1F *h_baseline_eta_range_4 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_4", 4, 4); // 1.8 < |eta| < 2.5

    TH1F *h_signal_eta_range_1 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_1", 1, 1); // 0 < |eta| < 0.6
    TH1F *h_signal_eta_range_2 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_2", 2, 2); // 0.6 < |eta| < 1.2
    TH1F *h_signal_eta_range_3 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_3", 3, 3); // 1.2 < |eta| < 1.8
    TH1F *h_signal_eta_range_4 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_4", 4, 4); // 1.8 < |eta| < 2.5

    //
    // The following part are the same as ytRealLeptonEfficiency_Data::make_muon_real_efficiency_plot()
    //
/*
    // For muons, there are 4 eta ranges: [0, 0.6], [0.6, 1.2], [1.2, 1.8], [1.8, 2.5]
    // Corresponding to eta bins: [1, 1], [2, 2], [3, 3], and [4, 4]
    // In order to get the correct bin, we have to use lower_edge + 0.01 and upper_edge - 0.01 to find the bin.
    // This is because FindBin() sometimes assign the lower_edge (upper_edge) to the correct_bin-1 (correct_bin+1).
    TH1F *h_baseline_eta_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_1", h_baseline_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(0.6 - 0.01));
    TH1F *h_baseline_eta_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_2", h_baseline_pt_eta->GetYaxis()->FindBin(0.6 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(1.2 - 0.01));
    TH1F *h_baseline_eta_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_3", h_baseline_pt_eta->GetYaxis()->FindBin(1.2 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(1.8 - 0.01));
    TH1F *h_baseline_eta_range_4 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_4", h_baseline_pt_eta->GetYaxis()->FindBin(1.8 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(2.5 - 0.01));

    TH1F *h_signal_eta_range_1 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_1", h_signal_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(0.6 - 0.01));
    TH1F *h_signal_eta_range_2 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_2", h_signal_pt_eta->GetYaxis()->FindBin(0.6 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(1.2 - 0.01));
    TH1F *h_signal_eta_range_3 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_3", h_signal_pt_eta->GetYaxis()->FindBin(1.2 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(1.8 - 0.01));
    TH1F *h_signal_eta_range_4 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_4", h_signal_pt_eta->GetYaxis()->FindBin(1.8 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(2.5 - 0.01));
*/
    TH1F *ratio_1 = (TH1F *)h_signal_eta_range_1->Clone();
    TH1F *ratio_2 = (TH1F *)h_signal_eta_range_2->Clone();
    TH1F *ratio_3 = (TH1F *)h_signal_eta_range_3->Clone();
    TH1F *ratio_4 = (TH1F *)h_signal_eta_range_4->Clone();

    ratio_1->Reset();
    ratio_2->Reset();
    ratio_3->Reset();
    ratio_4->Reset();

    ratio_1->SetName("ratio_1");
    ratio_2->SetName("ratio_2");
    ratio_3->SetName("ratio_3");
    ratio_4->SetName("ratio_4");

    ratio_1->Divide(h_signal_eta_range_1, h_baseline_eta_range_1, 1, 1, "B");
    ratio_2->Divide(h_signal_eta_range_2, h_baseline_eta_range_2, 1, 1, "B");
    ratio_3->Divide(h_signal_eta_range_3, h_baseline_eta_range_3, 1, 1, "B");
    ratio_4->Divide(h_signal_eta_range_4, h_baseline_eta_range_4, 1, 1, "B");

    TCanvas *real_efficiency_plot = new TCanvas("real_efficiency_plot", "Muon Real Efficiency", 500, 500);
    real_efficiency_plot->SetLeftMargin(0.12);
    real_efficiency_plot->cd();

    ratio_1->SetLineColor(kBlack);
    ratio_1->SetMarkerColor(kBlack);
    ratio_1->SetMarkerStyle(kFullCircle);
    ratio_1->SetXTitle("p_{T} [GeV]");
    ratio_1->SetYTitle("Muon Real Efficiency");
    ratio_1->GetYaxis()->SetTitleOffset(1.5);
    ratio_1->SetStats(kFALSE);
    ratio_1->SetMaximum(1.1);
    ratio_1->SetMinimum(0.5);
    ratio_1->SetTitle("");
    ratio_1->Draw();

    ratio_2->SetLineColor(kRed);
    ratio_2->SetMarkerColor(kRed);
    ratio_2->SetMarkerStyle(kFullSquare);
    ratio_2->Draw("same");

    ratio_3->SetLineColor(kBlue);
    ratio_3->SetMarkerColor(kBlue);
    ratio_3->SetMarkerStyle(kFullTriangleUp);
    ratio_3->Draw("same");

    ratio_4->SetLineColor(kGreen);
    ratio_4->SetMarkerColor(kGreen);
    ratio_4->SetMarkerStyle(kFullTriangleDown);
    ratio_4->Draw("same");

    float X_max = ratio_1->GetXaxis()->GetXmax();
    float X_min = ratio_1->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.4, 0.2, 0.8, 0.6);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry("ratio_1", "0 < |#eta| < 0.6");
    leg->AddEntry("ratio_2", "0.6 < |#eta| < 1.2");
    leg->AddEntry("ratio_3", "1.2 < |#eta| < 1.8");
    leg->AddEntry("ratio_4", "1.8 < |#eta| < 2.5");
    leg->Draw();

    real_efficiency_plot->SaveAs("real_efficiency_muon.pdf", "pdf");
}

void yt_make_real_efficiency_plots(bool truth_match = false, bool ttbar = false, bool Gtt = false) // Use 1-dim histograms
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *Zee_truth;
    TFile *Zmumu_truth;
    if (truth_match) {
        Zee_truth = TFile::Open(path + Zee_truth_root);
        Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    }

    TFile *ttbar_elec;
    TFile *ttbar_muon;
    if (ttbar) {
        ttbar_elec = TFile::Open(path + ttbar_elec_root);
        ttbar_muon = TFile::Open(path + ttbar_muon_root);
    }

    TFile *Gtt_elec;
    TFile *Gtt_muon;
    if (Gtt) {
        Gtt_elec = TFile::Open(path + Gtt_elec_root);
        Gtt_muon = TFile::Open(path + Gtt_muon_root);
    }

    TH1F *h_elec_baseline_pt = (TH1F *)data_elec->Get("h_baseline_pt");
    TH1F *h_elec_baseline_eta = (TH1F *)data_elec->Get("h_baseline_eta");
    TH1F *h_elec_baseline_dRjet = (TH1F *)data_elec->Get("h_baseline_dRjet");

    TH1F *h_elec_signal_pt = (TH1F *)data_elec->Get("h_signal_pt");
    TH1F *h_elec_signal_eta = (TH1F *)data_elec->Get("h_signal_eta");
    TH1F *h_elec_signal_dRjet = (TH1F *)data_elec->Get("h_signal_dRjet");

    TH1F *h_muon_baseline_pt = (TH1F *)data_muon->Get("h_baseline_pt");
    TH1F *h_muon_baseline_eta = (TH1F *)data_muon->Get("h_baseline_eta");
    TH1F *h_muon_baseline_dRjet = (TH1F *)data_muon->Get("h_baseline_dRjet");

    TH1F *h_muon_signal_pt = (TH1F *)data_muon->Get("h_signal_pt");
    TH1F *h_muon_signal_eta = (TH1F *)data_muon->Get("h_signal_eta");
    TH1F *h_muon_signal_dRjet = (TH1F *)data_muon->Get("h_signal_dRjet");

    static const float m_default_pt_bins[14] = {
        10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200
    };
    static const float m_electron_abs_eta_bins[11] = {
        0., 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47
    };
    static const float m_muon_coarse_eta_bins[5] = {
        0, 0.6, 1.2, 1.8, 2.5
    };
/*
    static const float m_default_deltaR_bins[21] = {
        0, 0.2, 0.4, 0.6, 0.8,
        1.0, 1.2, 1.4, 1.6, 1.8,
        2.0, 2.2, 2.4, 2.6, 2.8,
        3.0, 3.2, 3.4, 3.6, 3.8,
        4.0
    };
*/
    static const float m_fine_deltaR_bins[36] = {
        0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.2, 1.4, 1.6, 1.8,
        2.0, 2.2, 2.4, 2.6, 2.8,
        3.0, 3.2, 3.4, 3.6, 3.8,
        4.0
    };

    vector<float> m_pt_bins = vector<float> (m_default_pt_bins, m_default_pt_bins + sizeof(m_default_pt_bins) / sizeof(float));
    vector<float> m_elec_eta_bins = vector<float> (m_electron_abs_eta_bins, m_electron_abs_eta_bins + sizeof(m_electron_abs_eta_bins) / sizeof(float));
    vector<float> m_muon_eta_bins = vector<float> (m_muon_coarse_eta_bins, m_muon_coarse_eta_bins + sizeof(m_muon_coarse_eta_bins) / sizeof(float));
    //vector<float> m_deltaR_bins = vector<float> (m_default_deltaR_bins, m_default_deltaR_bins + sizeof(m_default_deltaR_bins) / sizeof(float));
    vector<float> m_deltaR_bins = vector<float> (m_fine_deltaR_bins, m_fine_deltaR_bins + sizeof(m_fine_deltaR_bins) / sizeof(float));

    int n_pt_bins = m_pt_bins.size() - 1;
    int n_elec_eta_bins = m_elec_eta_bins.size() - 1;
    int n_muon_eta_bins = m_muon_eta_bins.size() - 1;
    int n_deltaR_bins = m_deltaR_bins.size() - 1;

    TH1F *h_elec_eff_pt    = new TH1F("h_elec_eff_pt", "eff_pt;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    TH1F *h_elec_eff_eta   = new TH1F("h_elec_eff_eta", "eff_eta;|#eta|;Efficiency", n_elec_eta_bins, &m_elec_eta_bins[0]);
    TH1F *h_elec_eff_dRjet = new TH1F("h_elec_eff_dRjet", "eff_dRjet;dRjet;Efficiency", n_deltaR_bins, &m_deltaR_bins[0]);

    TH1F *h_muon_eff_pt    = new TH1F("h_muon_eff_pt", "eff_pt;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    TH1F *h_muon_eff_eta   = new TH1F("h_muon_eff_eta", "eff_eta;|#eta|;Efficiency", n_muon_eta_bins, &m_muon_eta_bins[0]);
    TH1F *h_muon_eff_dRjet = new TH1F("h_muon_eff_dRjet", "eff_dRjet;dRjet;Efficiency", n_deltaR_bins, &m_deltaR_bins[0]);

    // Calculate the real lepton efficiency
    h_elec_eff_pt->Divide(h_elec_signal_pt, h_elec_baseline_pt, 1, 1, "B");
    h_elec_eff_eta->Divide(h_elec_signal_eta, h_elec_baseline_eta, 1, 1, "B");
    h_elec_eff_dRjet->Divide(h_elec_signal_dRjet, h_elec_baseline_dRjet, 1, 1, "B");

    h_muon_eff_pt->Divide(h_muon_signal_pt, h_muon_baseline_pt, 1, 1, "B");
    h_muon_eff_eta->Divide(h_muon_signal_eta, h_muon_baseline_eta, 1, 1, "B");
    h_muon_eff_dRjet->Divide(h_muon_signal_dRjet, h_muon_baseline_dRjet, 1, 1, "B");
/*
    TH1F *h_elec_eff_pt = (TH1F *)data_elec->Get("h_eff_pt");
    TH1F *h_elec_eff_eta = (TH1F *)data_elec->Get("h_eff_eta");
    TH1F *h_elec_eff_dRjet = (TH1F *)data_elec->Get("h_eff_dRjet");

    TH1F *h_muon_eff_pt = (TH1F *)data_muon->Get("h_eff_pt");
    TH1F *h_muon_eff_eta = (TH1F *)data_muon->Get("h_eff_eta");
    TH1F *h_muon_eff_dRjet = (TH1F *)data_muon->Get("h_eff_dRjet");
*/
    TH1F *h_Zee_eff_pt = (TH1F *)Zee_TandP->Get("h_eff_pt");
    TH1F *h_Zee_eff_eta = (TH1F *)Zee_TandP->Get("h_eff_eta");
    TH1F *h_Zee_eff_dRjet = (TH1F *)Zee_TandP->Get("h_eff_dRjet");

    TH1F *h_Zmumu_eff_pt = (TH1F *)Zmumu_TandP->Get("h_eff_pt");
    TH1F *h_Zmumu_eff_eta = (TH1F *)Zmumu_TandP->Get("h_eff_eta");
    TH1F *h_Zmumu_eff_dRjet = (TH1F *)Zmumu_TandP->Get("h_eff_dRjet");

    TH1F *h_Zee_truth_eff_pt;
    TH1F *h_Zee_truth_eff_eta;
    TH1F *h_Zee_truth_eff_dRjet;

    TH1F *h_Zmumu_truth_eff_pt;
    TH1F *h_Zmumu_truth_eff_eta;
    TH1F *h_Zmumu_truth_eff_dRjet;

    if (truth_match) {
        h_Zee_truth_eff_pt = (TH1F *)Zee_truth->Get("h_eff_pt");
        h_Zee_truth_eff_eta = (TH1F *)Zee_truth->Get("h_eff_eta");
        h_Zee_truth_eff_dRjet = (TH1F *)Zee_truth->Get("h_eff_dRjet");

        h_Zmumu_truth_eff_pt = (TH1F *)Zmumu_truth->Get("h_eff_pt");
        h_Zmumu_truth_eff_eta = (TH1F *)Zmumu_truth->Get("h_eff_eta");
        h_Zmumu_truth_eff_dRjet = (TH1F *)Zmumu_truth->Get("h_eff_dRjet");
    }

    TH1F *h_ttbar_elec_eff_pt;
    TH1F *h_ttbar_elec_eff_eta;
    TH1F *h_ttbar_elec_eff_dRjet;

    TH1F *h_ttbar_muon_eff_pt;
    TH1F *h_ttbar_muon_eff_eta;
    TH1F *h_ttbar_muon_eff_dRjet;

    if (ttbar) {
        h_ttbar_elec_eff_pt = (TH1F *)ttbar_elec->Get("h_eff_pt");
        h_ttbar_elec_eff_eta = (TH1F *)ttbar_elec->Get("h_eff_eta");
        h_ttbar_elec_eff_dRjet = (TH1F *)ttbar_elec->Get("h_eff_dRjet");

        h_ttbar_muon_eff_pt = (TH1F *)ttbar_muon->Get("h_eff_pt");
        h_ttbar_muon_eff_eta = (TH1F *)ttbar_muon->Get("h_eff_eta");
        h_ttbar_muon_eff_dRjet = (TH1F *)ttbar_muon->Get("h_eff_dRjet");
    }

    TH1F *h_Gtt_elec_eff_pt;
    TH1F *h_Gtt_elec_eff_eta;
    TH1F *h_Gtt_elec_eff_dRjet;

    TH1F *h_Gtt_muon_eff_pt;
    TH1F *h_Gtt_muon_eff_eta;
    TH1F *h_Gtt_muon_eff_dRjet;

    if (Gtt) {
        h_Gtt_elec_eff_pt = (TH1F *)Gtt_elec->Get("h_eff_pt");
        h_Gtt_elec_eff_eta = (TH1F *)Gtt_elec->Get("h_eff_eta");
        h_Gtt_elec_eff_dRjet = (TH1F *)Gtt_elec->Get("h_eff_dRjet");

        h_Gtt_muon_eff_pt = (TH1F *)Gtt_muon->Get("h_eff_pt");
        h_Gtt_muon_eff_eta = (TH1F *)Gtt_muon->Get("h_eff_eta");
        h_Gtt_muon_eff_dRjet = (TH1F *)Gtt_muon->Get("h_eff_dRjet");
    }

    h_elec_eff_pt->SetName("h_elec_eff_pt");
    h_elec_eff_eta->SetName("h_elec_eff_eta");
    h_elec_eff_dRjet->SetName("h_elec_eff_dRjet");

    h_muon_eff_pt->SetName("h_muon_eff_pt");
    h_muon_eff_eta->SetName("h_muon_eff_eta");
    h_muon_eff_dRjet->SetName("h_muon_eff_dRjet");

    h_Zee_eff_pt->SetName("h_Zee_eff_pt");
    h_Zee_eff_eta->SetName("h_Zee_eff_eta");
    h_Zee_eff_dRjet->SetName("h_Zee_eff_dRjet");

    h_Zmumu_eff_pt->SetName("h_Zmumu_eff_pt");
    h_Zmumu_eff_eta->SetName("h_Zmumu_eff_eta");
    h_Zmumu_eff_dRjet->SetName("h_Zmumu_eff_dRjet");

    if (truth_match) {
        h_Zee_truth_eff_pt->SetName("h_Zee_truth_eff_pt");
        h_Zee_truth_eff_eta->SetName("h_Zee_truth_eff_eta");
        h_Zee_truth_eff_dRjet->SetName("h_Zee_truth_eff_dRjet");

        h_Zmumu_truth_eff_pt->SetName("h_Zmumu_truth_eff_pt");
        h_Zmumu_truth_eff_eta->SetName("h_Zmumu_truth_eff_eta");
        h_Zmumu_truth_eff_dRjet->SetName("h_Zmumu_truth_eff_dRjet");
    }

    if (ttbar) {
        h_ttbar_elec_eff_pt->SetName("h_ttbar_elec_eff_pt");
        h_ttbar_elec_eff_eta->SetName("h_ttbar_elec_eff_eta");
        h_ttbar_elec_eff_dRjet->SetName("h_ttbar_elec_eff_dRjet");

        h_ttbar_muon_eff_pt->SetName("h_ttbar_muon_eff_pt");
        h_ttbar_muon_eff_eta->SetName("h_ttbar_muon_eff_eta");
        h_ttbar_muon_eff_dRjet->SetName("h_ttbar_muon_eff_dRjet");
    }

    if (Gtt) {
        h_Gtt_elec_eff_pt->SetName("h_Gtt_elec_eff_pt");
        h_Gtt_elec_eff_eta->SetName("h_Gtt_elec_eff_eta");
        h_Gtt_elec_eff_dRjet->SetName("h_Gtt_elec_eff_dRjet");

        h_Gtt_muon_eff_pt->SetName("h_Gtt_muon_eff_pt");
        h_Gtt_muon_eff_eta->SetName("h_Gtt_muon_eff_eta");
        h_Gtt_muon_eff_dRjet->SetName("h_Gtt_muon_eff_dRjet");
    }

    TCanvas *real_efficiency_plot = new TCanvas("real_efficiency_plot", "Real Efficiency", 1500, 1000);
    real_efficiency_plot->Divide(3, 2);

    real_efficiency_plot->cd(1);
    gPad->SetLeftMargin(0.12);
    h_elec_eff_pt->SetMarkerColor(kBlack);
    h_elec_eff_pt->SetMarkerStyle(kFullCircle);
    h_elec_eff_pt->SetLineColor(kBlack);
    h_elec_eff_pt->SetMaximum(1.1);
    h_elec_eff_pt->SetMinimum(0);
    h_elec_eff_pt->SetTitle("");
    h_elec_eff_pt->SetXTitle("p_{T} [GeV]");
    h_elec_eff_pt->SetYTitle("Real electron efficiency");
    h_elec_eff_pt->GetYaxis()->SetTitleOffset(1.5);
    h_elec_eff_pt->SetStats(kFALSE);
    h_elec_eff_pt->Draw();

    h_Zee_eff_pt->SetMarkerColor(kRed);
    h_Zee_eff_pt->SetMarkerStyle(kFullSquare);
    h_Zee_eff_pt->SetLineColor(kRed);
    h_Zee_eff_pt->Draw("same");

    if (truth_match) {
        h_Zee_truth_eff_pt->SetMarkerColor(kBlue);
        h_Zee_truth_eff_pt->SetMarkerStyle(kFullTriangleUp);
        h_Zee_truth_eff_pt->SetLineColor(kBlue);
        h_Zee_truth_eff_pt->Draw("same");
    }

    if (ttbar) {
        h_ttbar_elec_eff_pt->SetMarkerColor(kMagenta);
        h_ttbar_elec_eff_pt->SetMarkerStyle(kFullDiamond);
        h_ttbar_elec_eff_pt->SetLineColor(kMagenta);
        h_ttbar_elec_eff_pt->Draw("same");
    }

    if (Gtt) {
        h_Gtt_elec_eff_pt->SetMarkerColor(kOrange);
        h_Gtt_elec_eff_pt->SetMarkerStyle(kFullCross);
        h_Gtt_elec_eff_pt->SetLineColor(kOrange);
        h_Gtt_elec_eff_pt->Draw("same");
    }

    TLegend *leg1 = new TLegend(0.25, 0.2, 0.7, 0.5);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.05);
    leg1->AddEntry("h_elec_eff_pt", "Data");
    leg1->AddEntry("h_Zee_eff_pt", "Z #rightarrow ee T&P MC");
    if (truth_match) {
        leg1->AddEntry("h_Zee_truth_eff_pt", "Z #rightarrow ee truth matched MC");
    }
    if (ttbar) {
        leg1->AddEntry("h_ttbar_elec_eff_pt", "ttbar MC");
    }
    if (Gtt) {
        leg1->AddEntry("h_Gtt_elec_eff_pt", "Gtt MC");
    }
    leg1->Draw();

    real_efficiency_plot->cd(2);
    gPad->SetLeftMargin(0.12);
    h_elec_eff_eta->SetMarkerColor(kBlack);
    h_elec_eff_eta->SetMarkerStyle(kFullCircle);
    h_elec_eff_eta->SetLineColor(kBlack);
    h_elec_eff_eta->SetMaximum(1.1);
    h_elec_eff_eta->SetMinimum(0);
    h_elec_eff_eta->GetXaxis()->SetRangeUser(0, 2.0);
    h_elec_eff_eta->SetTitle("");
    h_elec_eff_eta->SetXTitle("|#eta|");
    h_elec_eff_eta->SetYTitle("Real electron efficiency");
    h_elec_eff_eta->GetYaxis()->SetTitleOffset(1.5);
    h_elec_eff_eta->SetStats(kFALSE);
    h_elec_eff_eta->Draw();

    h_Zee_eff_eta->SetMarkerColor(kRed);
    h_Zee_eff_eta->SetMarkerStyle(kFullSquare);
    h_Zee_eff_eta->SetLineColor(kRed);
    h_Zee_eff_eta->Draw("same");

    if (truth_match) {
        h_Zee_truth_eff_eta->SetMarkerColor(kBlue);
        h_Zee_truth_eff_eta->SetMarkerStyle(kFullTriangleUp);
        h_Zee_truth_eff_eta->SetLineColor(kBlue);
        h_Zee_truth_eff_eta->Draw("same");
    }

    if (ttbar) {
        h_ttbar_elec_eff_eta->SetMarkerColor(kMagenta);
        h_ttbar_elec_eff_eta->SetMarkerStyle(kFullDiamond);
        h_ttbar_elec_eff_eta->SetLineColor(kMagenta);
        h_ttbar_elec_eff_eta->Draw("same");
    }

    if (Gtt) {
        h_Gtt_elec_eff_eta->SetMarkerColor(kOrange);
        h_Gtt_elec_eff_eta->SetMarkerStyle(kFullCross);
        h_Gtt_elec_eff_eta->SetLineColor(kOrange);
        h_Gtt_elec_eff_eta->Draw("same");
    }

    TLegend *leg2 = new TLegend(0.25, 0.2, 0.7, 0.5);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.05);
    leg2->AddEntry("h_elec_eff_eta", "Data");
    leg2->AddEntry("h_Zee_eff_eta", "Z #rightarrow ee T&P MC");
    if (truth_match) {
        leg2->AddEntry("h_Zee_truth_eff_eta", "Z #rightarrow ee truth matched MC");
    }
    if (ttbar) {
        leg2->AddEntry("h_ttbar_elec_eff_eta", "ttbar MC");
    }
    if (Gtt) {
        leg2->AddEntry("h_Gtt_elec_eff_eta", "Gtt MC");
    }
    leg2->Draw();

    real_efficiency_plot->cd(3);
    gPad->SetLeftMargin(0.12);
    h_elec_eff_dRjet->SetMarkerColor(kBlack);
    h_elec_eff_dRjet->SetMarkerStyle(kFullCircle);
    h_elec_eff_dRjet->SetLineColor(kBlack);
    h_elec_eff_dRjet->SetMaximum(1.1);
    h_elec_eff_dRjet->SetMinimum(0);
    h_elec_eff_dRjet->SetTitle("");
    h_elec_eff_dRjet->SetXTitle("#Delta R(e, jet)");
    h_elec_eff_dRjet->SetYTitle("Real electron efficiency");
    h_elec_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
    h_elec_eff_dRjet->SetStats(kFALSE);
    h_elec_eff_dRjet->Draw();

    h_Zee_eff_dRjet->SetMarkerColor(kRed);
    h_Zee_eff_dRjet->SetMarkerStyle(kFullSquare);
    h_Zee_eff_dRjet->SetLineColor(kRed);
    h_Zee_eff_dRjet->Draw("same");

    if (truth_match) {
        h_Zee_truth_eff_dRjet->SetMarkerColor(kBlue);
        h_Zee_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
        h_Zee_truth_eff_dRjet->SetLineColor(kBlue);
        h_Zee_truth_eff_dRjet->Draw("same");
    }

    if (ttbar) {
        h_ttbar_elec_eff_dRjet->SetMarkerColor(kMagenta);
        h_ttbar_elec_eff_dRjet->SetMarkerStyle(kFullDiamond);
        h_ttbar_elec_eff_dRjet->SetLineColor(kMagenta);
        h_ttbar_elec_eff_dRjet->Draw("same");
    }

    if (Gtt) {
        h_Gtt_elec_eff_dRjet->SetMarkerColor(kOrange);
        h_Gtt_elec_eff_dRjet->SetMarkerStyle(kFullCross);
        h_Gtt_elec_eff_dRjet->SetLineColor(kOrange);
        h_Gtt_elec_eff_dRjet->Draw("same");
    }

    TLegend *leg3 = new TLegend(0.25, 0.2, 0.7, 0.5);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.05);
    leg3->AddEntry("h_elec_eff_dRjet", "Data");
    leg3->AddEntry("h_Zee_eff_dRjet", "Z #rightarrow ee T&P MC");
    if (truth_match) {
        leg3->AddEntry("h_Zee_truth_eff_dRjet", "Z #rightarrow ee truth matched MC");
    }
    if (ttbar) {
        leg3->AddEntry("h_ttbar_elec_eff_dRjet", "ttbar MC");
    }
    if (Gtt) {
        leg3->AddEntry("h_Gtt_elec_eff_dRjet", "Gtt MC");
    }
    leg3->Draw();

    real_efficiency_plot->cd(4);
    gPad->SetLeftMargin(0.12);
    h_muon_eff_pt->SetMarkerColor(kBlack);
    h_muon_eff_pt->SetMarkerStyle(kFullCircle);
    h_muon_eff_pt->SetLineColor(kBlack);
    h_muon_eff_pt->SetMaximum(1.1);
    h_muon_eff_pt->SetMinimum(0);
    h_muon_eff_pt->SetTitle("");
    h_muon_eff_pt->SetYTitle("Real muon efficiency");
    h_muon_eff_pt->GetYaxis()->SetTitleOffset(1.5);
    h_muon_eff_pt->SetStats(kFALSE);
    h_muon_eff_pt->Draw();

    h_Zmumu_eff_pt->SetMarkerColor(kRed);
    h_Zmumu_eff_pt->SetMarkerStyle(kFullSquare);
    h_Zmumu_eff_pt->SetLineColor(kRed);
    h_Zmumu_eff_pt->Draw("same");

    if (truth_match) {
        h_Zmumu_truth_eff_pt->SetMarkerColor(kBlue);
        h_Zmumu_truth_eff_pt->SetMarkerStyle(kFullTriangleUp);
        h_Zmumu_truth_eff_pt->SetLineColor(kBlue);
        h_Zmumu_truth_eff_pt->Draw("same");
    }

    if (ttbar) {
        h_ttbar_muon_eff_pt->SetMarkerColor(kMagenta);
        h_ttbar_muon_eff_pt->SetMarkerStyle(kFullDiamond);
        h_ttbar_muon_eff_pt->SetLineColor(kMagenta);
        h_ttbar_muon_eff_pt->Draw("same");
    }

    if (Gtt) {
        h_Gtt_muon_eff_pt->SetMarkerColor(kOrange);
        h_Gtt_muon_eff_pt->SetMarkerStyle(kFullCross);
        h_Gtt_muon_eff_pt->SetLineColor(kOrange);
        h_Gtt_muon_eff_pt->Draw("same");
    }

    TLegend *leg4 = new TLegend(0.25, 0.2, 0.7, 0.5);
    leg4->SetBorderSize(0);
    leg4->SetFillColor(0);
    leg4->SetTextFont(42);
    leg4->SetTextSize(0.05);
    leg4->AddEntry("h_muon_eff_pt", "Data");
    leg4->AddEntry("h_Zmumu_eff_pt", "Z #rightarrow #mu#mu T&P MC");
    if (truth_match) {
        leg4->AddEntry("h_Zmumu_truth_eff_pt", "Z #rightarrow #mu#mu truth matched MC");
    }
    if (ttbar) {
        leg4->AddEntry("h_ttbar_muon_eff_pt", "ttbar MC");
    }
    if (Gtt) {
        leg4->AddEntry("h_Gtt_muon_eff_pt", "Gtt MC");
    }
    leg4->Draw();

    real_efficiency_plot->cd(5);
    gPad->SetLeftMargin(0.12);
    h_muon_eff_eta->SetMarkerColor(kBlack);
    h_muon_eff_eta->SetMarkerStyle(kFullCircle);
    h_muon_eff_eta->SetLineColor(kBlack);
    h_muon_eff_eta->SetMaximum(1.1);
    h_muon_eff_eta->SetMinimum(0);
    h_muon_eff_eta->SetTitle("");
    h_muon_eff_eta->SetXTitle("|#eta|");
    h_muon_eff_eta->SetYTitle("Real muon efficiency");
    h_muon_eff_eta->GetYaxis()->SetTitleOffset(1.5);
    h_muon_eff_eta->SetStats(kFALSE);
    h_muon_eff_eta->Draw();

    h_Zmumu_eff_eta->SetMarkerColor(kRed);
    h_Zmumu_eff_eta->SetMarkerStyle(kFullSquare);
    h_Zmumu_eff_eta->SetLineColor(kRed);
    h_Zmumu_eff_eta->Draw("same");

    if (truth_match) {
        h_Zmumu_truth_eff_eta->SetMarkerColor(kBlue);
        h_Zmumu_truth_eff_eta->SetMarkerStyle(kFullTriangleUp);
        h_Zmumu_truth_eff_eta->SetLineColor(kBlue);
        h_Zmumu_truth_eff_eta->Draw("same");
    }

    if (ttbar) {
        h_ttbar_muon_eff_eta->SetMarkerColor(kMagenta);
        h_ttbar_muon_eff_eta->SetMarkerStyle(kFullDiamond);
        h_ttbar_muon_eff_eta->SetLineColor(kMagenta);
        h_ttbar_muon_eff_eta->Draw("same");
    }

    if (Gtt) {
        h_Gtt_muon_eff_eta->SetMarkerColor(kOrange);
        h_Gtt_muon_eff_eta->SetMarkerStyle(kFullCross);
        h_Gtt_muon_eff_eta->SetLineColor(kOrange);
        h_Gtt_muon_eff_eta->Draw("same");
    }

    TLegend *leg5 = new TLegend(0.25, 0.2, 0.7, 0.5);
    leg5->SetBorderSize(0);
    leg5->SetFillColor(0);
    leg5->SetTextFont(42);
    leg5->SetTextSize(0.05);
    leg5->AddEntry("h_muon_eff_eta", "Data");
    leg5->AddEntry("h_Zmumu_eff_eta", "Z #rightarrow #mu#mu T&P MC");
    if (truth_match) {
        leg5->AddEntry("h_Zmumu_truth_eff_eta", "Z #rightarrow #mu#mu truth matched MC");
    }
    if (ttbar) {
        leg5->AddEntry("h_ttbar_muon_eff_eta", "ttbar MC");
    }
    if (Gtt) {
        leg5->AddEntry("h_Gtt_muon_eff_eta", "Gtt MC");
    }
    leg5->Draw();

    real_efficiency_plot->cd(6);
    gPad->SetLeftMargin(0.12);
    h_muon_eff_dRjet->SetMarkerColor(kBlack);
    h_muon_eff_dRjet->SetMarkerStyle(kFullCircle);
    h_muon_eff_dRjet->SetLineColor(kBlack);
    h_muon_eff_dRjet->SetMaximum(1.1);
    h_muon_eff_dRjet->SetMinimum(0);
    h_muon_eff_dRjet->SetTitle("");
    h_muon_eff_dRjet->SetXTitle("#Delta R(#mu, jet)");
    h_muon_eff_dRjet->SetYTitle("Real muon efficiency");
    h_muon_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
    h_muon_eff_dRjet->SetStats(kFALSE);
    h_muon_eff_dRjet->Draw();

    h_Zmumu_eff_dRjet->SetMarkerColor(kRed);
    h_Zmumu_eff_dRjet->SetMarkerStyle(kFullSquare);
    h_Zmumu_eff_dRjet->SetLineColor(kRed);
    h_Zmumu_eff_dRjet->Draw("same");

    if (truth_match) {
        h_Zmumu_truth_eff_dRjet->SetMarkerColor(kBlue);
        h_Zmumu_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
        h_Zmumu_truth_eff_dRjet->SetLineColor(kBlue);
        h_Zmumu_truth_eff_dRjet->Draw("same");
    }

    if (ttbar) {
        h_ttbar_muon_eff_dRjet->SetMarkerColor(kMagenta);
        h_ttbar_muon_eff_dRjet->SetMarkerStyle(kFullDiamond);
        h_ttbar_muon_eff_dRjet->SetLineColor(kMagenta);
        h_ttbar_muon_eff_dRjet->Draw("same");
    }

    if (Gtt) {
        h_Gtt_muon_eff_dRjet->SetMarkerColor(kOrange);
        h_Gtt_muon_eff_dRjet->SetMarkerStyle(kFullCross);
        h_Gtt_muon_eff_dRjet->SetLineColor(kOrange);
        h_Gtt_muon_eff_dRjet->Draw("same");
    }

    TLegend *leg6 = new TLegend(0.25, 0.2, 0.7, 0.5);
    leg6->SetBorderSize(0);
    leg6->SetFillColor(0);
    leg6->SetTextFont(42);
    leg6->SetTextSize(0.05);
    leg6->AddEntry("h_muon_eff_dRjet", "Data");
    leg6->AddEntry("h_Zmumu_eff_dRjet", "Z #rightarrow #mu#mu T&P MC");
    if (truth_match) {
        leg6->AddEntry("h_Zmumu_truth_eff_dRjet", "Z #rightarrow #mu#mu truth matched MC");
    }
    if (ttbar) {
        leg6->AddEntry("h_ttbar_muon_eff_dRjet", "ttbar MC");
    }
    if (Gtt) {
        leg6->AddEntry("h_Gtt_muon_eff_dRjet", "Gtt MC");
    }
    leg6->Draw();

    real_efficiency_plot->SaveAs("real_efficiency.pdf", "pdf");
}

void yt_make_real_efficiency_ratio_plots(TString lepton, TString plot, bool truth_match = false, bool ttbar = false, bool Gtt = false) // Use 1-dim histograms
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *Zee_truth;
    TFile *Zmumu_truth;
    if (truth_match) {
        Zee_truth = TFile::Open(path + Zee_truth_root);
        Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    }

    TFile *ttbar_elec;
    TFile *ttbar_muon;
    if (ttbar) {
        ttbar_elec = TFile::Open(path + ttbar_elec_root);
        ttbar_muon = TFile::Open(path + ttbar_muon_root);
    }

    TFile *Gtt_elec;
    TFile *Gtt_muon;
    if (Gtt) {
        Gtt_elec = TFile::Open(path + Gtt_elec_root);
        Gtt_muon = TFile::Open(path + Gtt_muon_root);
    }

    TH1F *h_elec_baseline_pt = (TH1F *)data_elec->Get("h_baseline_pt");
    TH1F *h_elec_baseline_eta = (TH1F *)data_elec->Get("h_baseline_eta");
    TH1F *h_elec_baseline_dRjet = (TH1F *)data_elec->Get("h_baseline_dRjet");

    TH1F *h_elec_signal_pt = (TH1F *)data_elec->Get("h_signal_pt");
    TH1F *h_elec_signal_eta = (TH1F *)data_elec->Get("h_signal_eta");
    TH1F *h_elec_signal_dRjet = (TH1F *)data_elec->Get("h_signal_dRjet");

    TH1F *h_muon_baseline_pt = (TH1F *)data_muon->Get("h_baseline_pt");
    TH1F *h_muon_baseline_eta = (TH1F *)data_muon->Get("h_baseline_eta");
    TH1F *h_muon_baseline_dRjet = (TH1F *)data_muon->Get("h_baseline_dRjet");

    TH1F *h_muon_signal_pt = (TH1F *)data_muon->Get("h_signal_pt");
    TH1F *h_muon_signal_eta = (TH1F *)data_muon->Get("h_signal_eta");
    TH1F *h_muon_signal_dRjet = (TH1F *)data_muon->Get("h_signal_dRjet");

    static const float m_default_pt_bins[14] = {
        10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200
    };
    static const float m_electron_abs_eta_bins[11] = {
        0., 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47
    };
    static const float m_muon_coarse_eta_bins[5] = {
        0, 0.6, 1.2, 1.8, 2.5
    };
    static const float m_fine_deltaR_bins[36] = {
        0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.2, 1.4, 1.6, 1.8,
        2.0, 2.2, 2.4, 2.6, 2.8,
        3.0, 3.2, 3.4, 3.6, 3.8,
        4.0
    };
    vector<float> m_pt_bins = vector<float> (m_default_pt_bins, m_default_pt_bins + sizeof(m_default_pt_bins) / sizeof(float));
    vector<float> m_elec_eta_bins = vector<float> (m_electron_abs_eta_bins, m_electron_abs_eta_bins + sizeof(m_electron_abs_eta_bins) / sizeof(float));
    vector<float> m_muon_eta_bins = vector<float> (m_muon_coarse_eta_bins, m_muon_coarse_eta_bins + sizeof(m_muon_coarse_eta_bins) / sizeof(float));
    vector<float> m_deltaR_bins = vector<float> (m_fine_deltaR_bins, m_fine_deltaR_bins + sizeof(m_fine_deltaR_bins) / sizeof(float));

    int n_pt_bins = m_pt_bins.size() - 1;
    int n_elec_eta_bins = m_elec_eta_bins.size() - 1;
    int n_muon_eta_bins = m_muon_eta_bins.size() - 1;
    int n_deltaR_bins = m_deltaR_bins.size() - 1;

    TH1F *h_elec_eff_pt    = new TH1F("h_elec_eff_pt", "eff_pt;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    TH1F *h_elec_eff_eta   = new TH1F("h_elec_eff_eta", "eff_eta;|#eta|;Efficiency", n_elec_eta_bins, &m_elec_eta_bins[0]);
    TH1F *h_elec_eff_dRjet = new TH1F("h_elec_eff_dRjet", "eff_dRjet;dRjet;Efficiency", n_deltaR_bins, &m_deltaR_bins[0]);

    TH1F *h_muon_eff_pt    = new TH1F("h_muon_eff_pt", "eff_pt;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    TH1F *h_muon_eff_eta   = new TH1F("h_muon_eff_eta", "eff_eta;|#eta|;Efficiency", n_muon_eta_bins, &m_muon_eta_bins[0]);
    TH1F *h_muon_eff_dRjet = new TH1F("h_muon_eff_dRjet", "eff_dRjet;dRjet;Efficiency", n_deltaR_bins, &m_deltaR_bins[0]);

    // Calculate the real lepton efficiency
    h_elec_eff_pt->Divide(h_elec_signal_pt, h_elec_baseline_pt, 1, 1, "B");
    h_elec_eff_eta->Divide(h_elec_signal_eta, h_elec_baseline_eta, 1, 1, "B");
    h_elec_eff_dRjet->Divide(h_elec_signal_dRjet, h_elec_baseline_dRjet, 1, 1, "B");

    h_muon_eff_pt->Divide(h_muon_signal_pt, h_muon_baseline_pt, 1, 1, "B");
    h_muon_eff_eta->Divide(h_muon_signal_eta, h_muon_baseline_eta, 1, 1, "B");
    h_muon_eff_dRjet->Divide(h_muon_signal_dRjet, h_muon_baseline_dRjet, 1, 1, "B");

    TH1F *h_Zee_eff_pt = (TH1F *)Zee_TandP->Get("h_eff_pt");
    TH1F *h_Zee_eff_eta = (TH1F *)Zee_TandP->Get("h_eff_eta");
    TH1F *h_Zee_eff_dRjet = (TH1F *)Zee_TandP->Get("h_eff_dRjet");

    TH1F *h_Zmumu_eff_pt = (TH1F *)Zmumu_TandP->Get("h_eff_pt");
    TH1F *h_Zmumu_eff_eta = (TH1F *)Zmumu_TandP->Get("h_eff_eta");
    TH1F *h_Zmumu_eff_dRjet = (TH1F *)Zmumu_TandP->Get("h_eff_dRjet");

    TH1F *h_Zee_truth_eff_pt;
    TH1F *h_Zee_truth_eff_eta;
    TH1F *h_Zee_truth_eff_dRjet;

    TH1F *h_Zmumu_truth_eff_pt;
    TH1F *h_Zmumu_truth_eff_eta;
    TH1F *h_Zmumu_truth_eff_dRjet;

    if (truth_match) {
        h_Zee_truth_eff_pt = (TH1F *)Zee_truth->Get("h_eff_pt");
        h_Zee_truth_eff_eta = (TH1F *)Zee_truth->Get("h_eff_eta");
        h_Zee_truth_eff_dRjet = (TH1F *)Zee_truth->Get("h_eff_dRjet");

        h_Zmumu_truth_eff_pt = (TH1F *)Zmumu_truth->Get("h_eff_pt");
        h_Zmumu_truth_eff_eta = (TH1F *)Zmumu_truth->Get("h_eff_eta");
        h_Zmumu_truth_eff_dRjet = (TH1F *)Zmumu_truth->Get("h_eff_dRjet");
    }

    TH1F *h_ttbar_elec_eff_pt;
    TH1F *h_ttbar_elec_eff_eta;
    TH1F *h_ttbar_elec_eff_dRjet;

    TH1F *h_ttbar_muon_eff_pt;
    TH1F *h_ttbar_muon_eff_eta;
    TH1F *h_ttbar_muon_eff_dRjet;

    if (ttbar) {
        h_ttbar_elec_eff_pt = (TH1F *)ttbar_elec->Get("h_eff_pt");
        h_ttbar_elec_eff_eta = (TH1F *)ttbar_elec->Get("h_eff_eta");
        h_ttbar_elec_eff_dRjet = (TH1F *)ttbar_elec->Get("h_eff_dRjet");

        h_ttbar_muon_eff_pt = (TH1F *)ttbar_muon->Get("h_eff_pt");
        h_ttbar_muon_eff_eta = (TH1F *)ttbar_muon->Get("h_eff_eta");
        h_ttbar_muon_eff_dRjet = (TH1F *)ttbar_muon->Get("h_eff_dRjet");
    }

    TH1F *h_Gtt_elec_eff_pt;
    TH1F *h_Gtt_elec_eff_eta;
    TH1F *h_Gtt_elec_eff_dRjet;

    TH1F *h_Gtt_muon_eff_pt;
    TH1F *h_Gtt_muon_eff_eta;
    TH1F *h_Gtt_muon_eff_dRjet;

    if (Gtt) {
        h_Gtt_elec_eff_pt = (TH1F *)Gtt_elec->Get("h_eff_pt");
        h_Gtt_elec_eff_eta = (TH1F *)Gtt_elec->Get("h_eff_eta");
        h_Gtt_elec_eff_dRjet = (TH1F *)Gtt_elec->Get("h_eff_dRjet");

        h_Gtt_muon_eff_pt = (TH1F *)Gtt_muon->Get("h_eff_pt");
        h_Gtt_muon_eff_eta = (TH1F *)Gtt_muon->Get("h_eff_eta");
        h_Gtt_muon_eff_dRjet = (TH1F *)Gtt_muon->Get("h_eff_dRjet");
    }

    h_elec_eff_pt->SetName("h_elec_eff_pt");
    h_elec_eff_eta->SetName("h_elec_eff_eta");
    h_elec_eff_dRjet->SetName("h_elec_eff_dRjet");

    h_muon_eff_pt->SetName("h_muon_eff_pt");
    h_muon_eff_eta->SetName("h_muon_eff_eta");
    h_muon_eff_dRjet->SetName("h_muon_eff_dRjet");

    h_Zee_eff_pt->SetName("h_Zee_eff_pt");
    h_Zee_eff_eta->SetName("h_Zee_eff_eta");
    h_Zee_eff_dRjet->SetName("h_Zee_eff_dRjet");

    h_Zmumu_eff_pt->SetName("h_Zmumu_eff_pt");
    h_Zmumu_eff_eta->SetName("h_Zmumu_eff_eta");
    h_Zmumu_eff_dRjet->SetName("h_Zmumu_eff_dRjet");

    if (truth_match) {
        h_Zee_truth_eff_pt->SetName("h_Zee_truth_eff_pt");
        h_Zee_truth_eff_eta->SetName("h_Zee_truth_eff_eta");
        h_Zee_truth_eff_dRjet->SetName("h_Zee_truth_eff_dRjet");

        h_Zmumu_truth_eff_pt->SetName("h_Zmumu_truth_eff_pt");
        h_Zmumu_truth_eff_eta->SetName("h_Zmumu_truth_eff_eta");
        h_Zmumu_truth_eff_dRjet->SetName("h_Zmumu_truth_eff_dRjet");
    }

    if (ttbar) {
        h_ttbar_elec_eff_pt->SetName("h_ttbar_elec_eff_pt");
        h_ttbar_elec_eff_eta->SetName("h_ttbar_elec_eff_eta");
        h_ttbar_elec_eff_dRjet->SetName("h_ttbar_elec_eff_dRjet");

        h_ttbar_muon_eff_pt->SetName("h_ttbar_muon_eff_pt");
        h_ttbar_muon_eff_eta->SetName("h_ttbar_muon_eff_eta");
        h_ttbar_muon_eff_dRjet->SetName("h_ttbar_muon_eff_dRjet");
    }

    if (Gtt) {
        h_Gtt_elec_eff_pt->SetName("h_Gtt_elec_eff_pt");
        h_Gtt_elec_eff_eta->SetName("h_Gtt_elec_eff_eta");
        h_Gtt_elec_eff_dRjet->SetName("h_Gtt_elec_eff_dRjet");

        h_Gtt_muon_eff_pt->SetName("h_Gtt_muon_eff_pt");
        h_Gtt_muon_eff_eta->SetName("h_Gtt_muon_eff_eta");
        h_Gtt_muon_eff_dRjet->SetName("h_Gtt_muon_eff_dRjet");
    }

    TCanvas *real_efficiency_ratio_plot = new TCanvas("real_efficiency_ratio_plot", "Real Efficiency", 600, 600);

    //Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetRightMargin(0.08);
    pad1->SetLeftMargin(0.12);
    //pad1->SetGridy(); // grid lines
    //pad1->SetLogx();
    pad1->Draw();

    // lower plot will be in pad
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetRightMargin(0.08);
    pad2->SetLeftMargin(0.12);
    pad2->SetGridy(); // grid lines
    //pad2->SetLogx();
    pad2->Draw();

    pad1->cd(); // pad1 becomes the current pad
    //pad1->SetFrameLineWidth(2);

    // Draw curve here

    if (plot == "pt") {
        if (lepton == "electron") {
            h_elec_eff_pt->SetMarkerColor(kBlack);
            h_elec_eff_pt->SetMarkerStyle(kFullCircle);
            h_elec_eff_pt->SetLineColor(kBlack);
            h_elec_eff_pt->SetMaximum(1.1);
            h_elec_eff_pt->SetMinimum(0);
            h_elec_eff_pt->SetTitle("");
            h_elec_eff_pt->SetXTitle("p_{T} [GeV]");
            h_elec_eff_pt->SetYTitle("Real electron efficiency");
            h_elec_eff_pt->GetYaxis()->SetTitleOffset(1.5);
            h_elec_eff_pt->SetStats(kFALSE);
            h_elec_eff_pt->Draw();

            h_Zee_eff_pt->SetMarkerColor(kRed);
            h_Zee_eff_pt->SetMarkerStyle(kFullSquare);
            h_Zee_eff_pt->SetLineColor(kRed);
            h_Zee_eff_pt->Draw("same");

            if (truth_match) {
                h_Zee_truth_eff_pt->SetMarkerColor(kBlue);
                h_Zee_truth_eff_pt->SetMarkerStyle(kFullTriangleUp);
                h_Zee_truth_eff_pt->SetLineColor(kBlue);
                h_Zee_truth_eff_pt->Draw("same");
            }

            if (ttbar) {
                h_ttbar_elec_eff_pt->SetMarkerColor(kMagenta);
                h_ttbar_elec_eff_pt->SetMarkerStyle(kFullDiamond);
                h_ttbar_elec_eff_pt->SetLineColor(kMagenta);
                h_ttbar_elec_eff_pt->Draw("same");
            }

            if (Gtt) {
                h_Gtt_elec_eff_pt->SetMarkerColor(kOrange);
                h_Gtt_elec_eff_pt->SetMarkerStyle(kFullCross);
                h_Gtt_elec_eff_pt->SetLineColor(kOrange);
                h_Gtt_elec_eff_pt->Draw("same");
            }

            TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->AddEntry("h_elec_eff_pt", "Data");
            leg->AddEntry("h_Zee_eff_pt", "Z #rightarrow ee T&P MC");
            if (truth_match) {
                leg->AddEntry("h_Zee_truth_eff_pt", "Z #rightarrow ee truth matched MC");
            }
            if (ttbar) {
                leg->AddEntry("h_ttbar_elec_eff_pt", "ttbar MC");
            }
            if (Gtt) {
                leg->AddEntry("h_Gtt_elec_eff_pt", "Gtt MC");
            }
            leg->Draw();

            pad2->cd(); // pad2 becomes the current pad

            double x_min = h_elec_eff_pt->GetXaxis()->GetXmin();
            double x_max = h_elec_eff_pt->GetXaxis()->GetXmax();

            TH1F *frame = pad2->DrawFrame(x_min, 0.9, x_max, 1.09);
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetYaxis()->SetNdivisions(405);
            frame->SetLineWidth(1);
            frame->SetXTitle("p_{T} [GeV]");
            frame->GetXaxis()->SetTitleSize(20);
            frame->GetXaxis()->SetTitleFont(47);
            frame->GetXaxis()->SetTitleOffset(3.0);
            frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetXaxis()->SetLabelSize(15);
            frame->SetYTitle("Data / MC");
            frame->GetYaxis()->SetTitleSize(17);
            frame->GetYaxis()->SetTitleFont(43);
            frame->GetYaxis()->SetTitleOffset(1.5);
            frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetYaxis()->SetLabelSize(16);
            frame->Draw();

            TH1F *ratio_TandP = (TH1F *)h_elec_eff_pt->Clone();
            ratio_TandP->Sumw2();
            ratio_TandP->Divide(h_Zee_eff_pt);
            ratio_TandP->SetMarkerColor(kRed);
            ratio_TandP->SetMarkerStyle(kFullSquare);
            ratio_TandP->SetLineColor(kRed);
            ratio_TandP->Draw("same");

            if (truth_match) {
                TH1F *ratio_truth = (TH1F *)h_elec_eff_pt->Clone();
                ratio_truth->Sumw2();
                ratio_truth->Divide(h_Zee_truth_eff_pt);
                ratio_truth->SetMarkerColor(kBlue);
                ratio_truth->SetMarkerStyle(kFullTriangleUp);
                ratio_truth->SetLineColor(kBlue);
                ratio_truth->Draw("same");
            }

            if (ttbar) {
                TH1F *ratio_ttbar = (TH1F *)h_elec_eff_pt->Clone();
                ratio_ttbar->Sumw2();
                ratio_ttbar->Divide(h_ttbar_elec_eff_pt);
                ratio_ttbar->SetMarkerColor(kMagenta);
                ratio_ttbar->SetMarkerStyle(kFullDiamond);
                ratio_ttbar->SetLineColor(kMagenta);
                ratio_ttbar->Draw("same");
            }

            if (Gtt) {
                TH1F *ratio_gtt = (TH1F *)h_elec_eff_pt->Clone();
                ratio_gtt->Sumw2();
                ratio_gtt->Divide(h_Gtt_elec_eff_pt);
                ratio_gtt->SetMarkerColor(kOrange);
                ratio_gtt->SetMarkerStyle(kFullCross);
                ratio_gtt->SetLineColor(kOrange);
                ratio_gtt->Draw("same");
            }
        }
        else if (lepton == "muon") {
            h_muon_eff_pt->SetMarkerColor(kBlack);
            h_muon_eff_pt->SetMarkerStyle(kFullCircle);
            h_muon_eff_pt->SetLineColor(kBlack);
            h_muon_eff_pt->SetMaximum(1.1);
            h_muon_eff_pt->SetMinimum(0);
            h_muon_eff_pt->SetTitle("");
            h_muon_eff_pt->SetYTitle("Real muon efficiency");
            h_muon_eff_pt->GetYaxis()->SetTitleOffset(1.5);
            h_muon_eff_pt->SetStats(kFALSE);
            h_muon_eff_pt->Draw();

            h_Zmumu_eff_pt->SetMarkerColor(kRed);
            h_Zmumu_eff_pt->SetMarkerStyle(kFullSquare);
            h_Zmumu_eff_pt->SetLineColor(kRed);
            h_Zmumu_eff_pt->Draw("same");

            if (truth_match) {
                h_Zmumu_truth_eff_pt->SetMarkerColor(kBlue);
                h_Zmumu_truth_eff_pt->SetMarkerStyle(kFullTriangleUp);
                h_Zmumu_truth_eff_pt->SetLineColor(kBlue);
                h_Zmumu_truth_eff_pt->Draw("same");
            }

            if (ttbar) {
                h_ttbar_muon_eff_pt->SetMarkerColor(kMagenta);
                h_ttbar_muon_eff_pt->SetMarkerStyle(kFullDiamond);
                h_ttbar_muon_eff_pt->SetLineColor(kMagenta);
                h_ttbar_muon_eff_pt->Draw("same");
            }

            if (Gtt) {
                h_Gtt_muon_eff_pt->SetMarkerColor(kOrange);
                h_Gtt_muon_eff_pt->SetMarkerStyle(kFullCross);
                h_Gtt_muon_eff_pt->SetLineColor(kOrange);
                h_Gtt_muon_eff_pt->Draw("same");
            }

            TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->AddEntry("h_muon_eff_pt", "Data");
            leg->AddEntry("h_Zmumu_eff_pt", "Z #rightarrow #mu#mu T&P MC");
            if (truth_match) {
                leg->AddEntry("h_Zmumu_truth_eff_pt", "Z #rightarrow #mu#mu truth matched MC");
            }
            if (ttbar) {
                leg->AddEntry("h_ttbar_muon_eff_pt", "ttbar MC");
            }
            if (Gtt) {
                leg->AddEntry("h_Gtt_muon_eff_pt", "Gtt MC");
            }
            leg->Draw();

            pad2->cd(); // pad2 becomes the current pad

            double x_min = h_muon_eff_pt->GetXaxis()->GetXmin();
            double x_max = h_muon_eff_pt->GetXaxis()->GetXmax();

            TH1F *frame = pad2->DrawFrame(x_min, 0.9, x_max, 1.29);
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetYaxis()->SetNdivisions(405);
            frame->SetLineWidth(1);
            frame->SetXTitle("p_{T} [GeV]");
            frame->GetXaxis()->SetTitleSize(20);
            frame->GetXaxis()->SetTitleFont(47);
            frame->GetXaxis()->SetTitleOffset(3.0);
            frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetXaxis()->SetLabelSize(15);
            frame->SetYTitle("Data / MC");
            frame->GetYaxis()->SetTitleSize(17);
            frame->GetYaxis()->SetTitleFont(43);
            frame->GetYaxis()->SetTitleOffset(1.5);
            frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetYaxis()->SetLabelSize(16);
            frame->Draw();

            TH1F *ratio_TandP = (TH1F *)h_muon_eff_pt->Clone();
            ratio_TandP->Sumw2();
            ratio_TandP->Divide(h_Zmumu_eff_pt);
            ratio_TandP->SetMarkerColor(kRed);
            ratio_TandP->SetMarkerStyle(kFullSquare);
            ratio_TandP->SetLineColor(kRed);
            ratio_TandP->Draw("same");

            if (truth_match) {
                TH1F *ratio_truth = (TH1F *)h_muon_eff_pt->Clone();
                ratio_truth->Sumw2();
                ratio_truth->Divide(h_Zmumu_truth_eff_pt);
                ratio_truth->SetMarkerColor(kBlue);
                ratio_truth->SetMarkerStyle(kFullTriangleUp);
                ratio_truth->SetLineColor(kBlue);
                ratio_truth->Draw("same");
            }

            if (ttbar) {
                TH1F *ratio_ttbar = (TH1F *)h_muon_eff_pt->Clone();
                ratio_ttbar->Sumw2();
                ratio_ttbar->Divide(h_ttbar_muon_eff_pt);
                ratio_ttbar->SetMarkerColor(kMagenta);
                ratio_ttbar->SetMarkerStyle(kFullDiamond);
                ratio_ttbar->SetLineColor(kMagenta);
                ratio_ttbar->Draw("same");
            }

            if (Gtt) {
                TH1F *ratio_gtt = (TH1F *)h_muon_eff_pt->Clone();
                ratio_gtt->Sumw2();
                ratio_gtt->Divide(h_Gtt_muon_eff_pt);
                ratio_gtt->SetMarkerColor(kOrange);
                ratio_gtt->SetMarkerStyle(kFullCross);
                ratio_gtt->SetLineColor(kOrange);
                ratio_gtt->Draw("same");
            }
        }
    }
    else if (plot == "eta") {
        if (lepton == "electron") {
            h_elec_eff_eta->SetMarkerColor(kBlack);
            h_elec_eff_eta->SetMarkerStyle(kFullCircle);
            h_elec_eff_eta->SetLineColor(kBlack);
            h_elec_eff_eta->SetMaximum(1.1);
            h_elec_eff_eta->SetMinimum(0);
            h_elec_eff_eta->GetXaxis()->SetRangeUser(0, 2.0);
            h_elec_eff_eta->SetTitle("");
            h_elec_eff_eta->SetXTitle("|#eta|");
            h_elec_eff_eta->SetYTitle("Real electron efficiency");
            h_elec_eff_eta->GetYaxis()->SetTitleOffset(1.5);
            h_elec_eff_eta->SetStats(kFALSE);
            h_elec_eff_eta->Draw();

            h_Zee_eff_eta->SetMarkerColor(kRed);
            h_Zee_eff_eta->SetMarkerStyle(kFullSquare);
            h_Zee_eff_eta->SetLineColor(kRed);
            h_Zee_eff_eta->Draw("same");

            if (truth_match) {
                h_Zee_truth_eff_eta->SetMarkerColor(kBlue);
                h_Zee_truth_eff_eta->SetMarkerStyle(kFullTriangleUp);
                h_Zee_truth_eff_eta->SetLineColor(kBlue);
                h_Zee_truth_eff_eta->Draw("same");
            }

            if (ttbar) {
                h_ttbar_elec_eff_eta->SetMarkerColor(kMagenta);
                h_ttbar_elec_eff_eta->SetMarkerStyle(kFullDiamond);
                h_ttbar_elec_eff_eta->SetLineColor(kMagenta);
                h_ttbar_elec_eff_eta->Draw("same");
            }

            if (Gtt) {
                h_Gtt_elec_eff_eta->SetMarkerColor(kOrange);
                h_Gtt_elec_eff_eta->SetMarkerStyle(kFullCross);
                h_Gtt_elec_eff_eta->SetLineColor(kOrange);
                h_Gtt_elec_eff_eta->Draw("same");
            }

            TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->AddEntry("h_elec_eff_eta", "Data");
            leg->AddEntry("h_Zee_eff_eta", "Z #rightarrow ee T&P MC");
            if (truth_match) {
                leg->AddEntry("h_Zee_truth_eff_eta", "Z #rightarrow ee truth matched MC");
            }
            if (ttbar) {
                leg->AddEntry("h_ttbar_elec_eff_eta", "ttbar MC");
            }
            if (Gtt) {
                leg->AddEntry("h_Gtt_elec_eff_eta", "Gtt MC");
            }
            leg->Draw();

            pad2->cd(); // pad2 becomes the current pad

            double x_min = h_elec_eff_eta->GetXaxis()->GetXmin();
            double x_max = h_elec_eff_eta->GetXaxis()->GetXmax();

            TH1F *frame = pad2->DrawFrame(x_min, 0.9, 2.0, 1.09);
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetYaxis()->SetNdivisions(405);
            frame->SetLineWidth(1);
            frame->SetXTitle("|#eta|");
            frame->GetXaxis()->SetTitleSize(20);
            frame->GetXaxis()->SetTitleFont(47);
            frame->GetXaxis()->SetTitleOffset(3.0);
            frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetXaxis()->SetLabelSize(15);
            frame->SetYTitle("Data / MC");
            frame->GetYaxis()->SetTitleSize(17);
            frame->GetYaxis()->SetTitleFont(43);
            frame->GetYaxis()->SetTitleOffset(1.5);
            frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetYaxis()->SetLabelSize(16);
            frame->Draw();

            TH1F *ratio_TandP = (TH1F *)h_elec_eff_eta->Clone();
            ratio_TandP->Sumw2();
            ratio_TandP->Divide(h_Zee_eff_eta);
            ratio_TandP->SetMarkerColor(kRed);
            ratio_TandP->SetMarkerStyle(kFullSquare);
            ratio_TandP->SetLineColor(kRed);
            ratio_TandP->Draw("same");

            if (truth_match) {
                TH1F *ratio_truth = (TH1F *)h_elec_eff_eta->Clone();
                ratio_truth->Sumw2();
                ratio_truth->Divide(h_Zee_truth_eff_eta);
                ratio_truth->SetMarkerColor(kBlue);
                ratio_truth->SetMarkerStyle(kFullTriangleUp);
                ratio_truth->SetLineColor(kBlue);
                ratio_truth->Draw("same");
            }

            if (ttbar) {
                TH1F *ratio_ttbar = (TH1F *)h_elec_eff_eta->Clone();
                ratio_ttbar->Sumw2();
                ratio_ttbar->Divide(h_ttbar_elec_eff_eta);
                ratio_ttbar->SetMarkerColor(kMagenta);
                ratio_ttbar->SetMarkerStyle(kFullDiamond);
                ratio_ttbar->SetLineColor(kMagenta);
                ratio_ttbar->Draw("same");
            }

            if (Gtt) {
                TH1F *ratio_gtt = (TH1F *)h_elec_eff_eta->Clone();
                ratio_gtt->Sumw2();
                ratio_gtt->Divide(h_Gtt_elec_eff_eta);
                ratio_gtt->SetMarkerColor(kOrange);
                ratio_gtt->SetMarkerStyle(kFullCross);
                ratio_gtt->SetLineColor(kOrange);
                ratio_gtt->Draw("same");
            }
        }
        else if (lepton == "muon") {
            h_muon_eff_eta->SetMarkerColor(kBlack);
            h_muon_eff_eta->SetMarkerStyle(kFullCircle);
            h_muon_eff_eta->SetLineColor(kBlack);
            h_muon_eff_eta->SetMaximum(1.1);
            h_muon_eff_eta->SetMinimum(0);
            h_muon_eff_eta->SetTitle("");
            h_muon_eff_eta->SetXTitle("|#eta|");
            h_muon_eff_eta->SetYTitle("Real muon efficiency");
            h_muon_eff_eta->GetYaxis()->SetTitleOffset(1.5);
            h_muon_eff_eta->SetStats(kFALSE);
            h_muon_eff_eta->Draw();

            h_Zmumu_eff_eta->SetMarkerColor(kRed);
            h_Zmumu_eff_eta->SetMarkerStyle(kFullSquare);
            h_Zmumu_eff_eta->SetLineColor(kRed);
            h_Zmumu_eff_eta->Draw("same");

            if (truth_match) {
                h_Zmumu_truth_eff_eta->SetMarkerColor(kBlue);
                h_Zmumu_truth_eff_eta->SetMarkerStyle(kFullTriangleUp);
                h_Zmumu_truth_eff_eta->SetLineColor(kBlue);
                h_Zmumu_truth_eff_eta->Draw("same");
            }

            if (ttbar) {
                h_ttbar_muon_eff_eta->SetMarkerColor(kMagenta);
                h_ttbar_muon_eff_eta->SetMarkerStyle(kFullDiamond);
                h_ttbar_muon_eff_eta->SetLineColor(kMagenta);
                h_ttbar_muon_eff_eta->Draw("same");
            }

            if (Gtt) {
                h_Gtt_muon_eff_eta->SetMarkerColor(kOrange);
                h_Gtt_muon_eff_eta->SetMarkerStyle(kFullCross);
                h_Gtt_muon_eff_eta->SetLineColor(kOrange);
                h_Gtt_muon_eff_eta->Draw("same");
            }

            TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->AddEntry("h_muon_eff_eta", "Data");
            leg->AddEntry("h_Zmumu_eff_eta", "Z #rightarrow #mu#mu T&P MC");
            if (truth_match) {
                leg->AddEntry("h_Zmumu_truth_eff_eta", "Z #rightarrow #mu#mu truth matched MC");
            }
            if (ttbar) {
                leg->AddEntry("h_ttbar_muon_eff_eta", "ttbar MC");
            }
            if (Gtt) {
                leg->AddEntry("h_Gtt_muon_eff_eta", "Gtt MC");
            }
            leg->Draw();

            pad2->cd(); // pad2 becomes the current pad

            double x_min = h_Zmumu_eff_eta->GetXaxis()->GetXmin();
            double x_max = h_Zmumu_eff_eta->GetXaxis()->GetXmax();

            TH1F *frame = pad2->DrawFrame(x_min, 0.9, x_max, 1.09);
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetYaxis()->SetNdivisions(405);
            frame->SetLineWidth(1);
            frame->SetXTitle("|#eta|");
            frame->GetXaxis()->SetTitleSize(20);
            frame->GetXaxis()->SetTitleFont(47);
            frame->GetXaxis()->SetTitleOffset(3.0);
            frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetXaxis()->SetLabelSize(15);
            frame->SetYTitle("Data / MC");
            frame->GetYaxis()->SetTitleSize(17);
            frame->GetYaxis()->SetTitleFont(43);
            frame->GetYaxis()->SetTitleOffset(1.5);
            frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetYaxis()->SetLabelSize(16);
            frame->Draw();

            TH1F *ratio_TandP = (TH1F *)h_muon_eff_eta->Clone();
            ratio_TandP->Sumw2();
            ratio_TandP->Divide(h_Zmumu_eff_eta);
            ratio_TandP->SetMarkerColor(kRed);
            ratio_TandP->SetMarkerStyle(kFullSquare);
            ratio_TandP->SetLineColor(kRed);
            ratio_TandP->Draw("same");

            if (truth_match) {
                TH1F *ratio_truth = (TH1F *)h_muon_eff_eta->Clone();
                ratio_truth->Sumw2();
                ratio_truth->Divide(h_Zmumu_truth_eff_eta);
                ratio_truth->SetMarkerColor(kBlue);
                ratio_truth->SetMarkerStyle(kFullTriangleUp);
                ratio_truth->SetLineColor(kBlue);
                ratio_truth->Draw("same");
            }

            if (ttbar) {
                TH1F *ratio_ttbar = (TH1F *)h_muon_eff_eta->Clone();
                ratio_ttbar->Sumw2();
                ratio_ttbar->Divide(h_ttbar_muon_eff_eta);
                ratio_ttbar->SetMarkerColor(kMagenta);
                ratio_ttbar->SetMarkerStyle(kFullDiamond);
                ratio_ttbar->SetLineColor(kMagenta);
                ratio_ttbar->Draw("same");
            }

            if (Gtt) {
                TH1F *ratio_gtt = (TH1F *)h_muon_eff_eta->Clone();
                ratio_gtt->Sumw2();
                ratio_gtt->Divide(h_Gtt_muon_eff_eta);
                ratio_gtt->SetMarkerColor(kOrange);
                ratio_gtt->SetMarkerStyle(kFullCross);
                ratio_gtt->SetLineColor(kOrange);
                ratio_gtt->Draw("same");
            }
        }
    }
    else if (plot == "dRjet") {
        if (lepton == "electron") {
            h_elec_eff_dRjet->SetMarkerColor(kBlack);
            h_elec_eff_dRjet->SetMarkerStyle(kFullCircle);
            h_elec_eff_dRjet->SetLineColor(kBlack);
            h_elec_eff_dRjet->SetMaximum(1.1);
            h_elec_eff_dRjet->SetMinimum(0);
            h_elec_eff_dRjet->SetTitle("");
            h_elec_eff_dRjet->SetXTitle("#Delta R(e, jet)");
            h_elec_eff_dRjet->SetYTitle("Real electron efficiency");
            h_elec_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
            h_elec_eff_dRjet->SetStats(kFALSE);
            h_elec_eff_dRjet->Draw();

            h_Zee_eff_dRjet->SetMarkerColor(kRed);
            h_Zee_eff_dRjet->SetMarkerStyle(kFullSquare);
            h_Zee_eff_dRjet->SetLineColor(kRed);
            h_Zee_eff_dRjet->Draw("same");

            if (truth_match) {
                h_Zee_truth_eff_dRjet->SetMarkerColor(kBlue);
                h_Zee_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
                h_Zee_truth_eff_dRjet->SetLineColor(kBlue);
                h_Zee_truth_eff_dRjet->Draw("same");
            }

            if (ttbar) {
                h_ttbar_elec_eff_dRjet->SetMarkerColor(kMagenta);
                h_ttbar_elec_eff_dRjet->SetMarkerStyle(kFullDiamond);
                h_ttbar_elec_eff_dRjet->SetLineColor(kMagenta);
                h_ttbar_elec_eff_dRjet->Draw("same");
            }

            if (Gtt) {
                h_Gtt_elec_eff_dRjet->SetMarkerColor(kOrange);
                h_Gtt_elec_eff_dRjet->SetMarkerStyle(kFullCross);
                h_Gtt_elec_eff_dRjet->SetLineColor(kOrange);
                h_Gtt_elec_eff_dRjet->Draw("same");
            }

            TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->AddEntry("h_elec_eff_dRjet", "Data");
            leg->AddEntry("h_Zee_eff_dRjet", "Z #rightarrow ee T&P MC");
            if (truth_match) {
                leg->AddEntry("h_Zee_truth_eff_dRjet", "Z #rightarrow ee truth matched MC");
            }
            if (ttbar) {
                leg->AddEntry("h_ttbar_elec_eff_dRjet", "ttbar MC");
            }
            if (Gtt) {
                leg->AddEntry("h_Gtt_elec_eff_dRjet", "Gtt MC");
            }
            leg->Draw();

            pad2->cd(); // pad2 becomes the current pad

            double x_min = h_elec_eff_dRjet->GetXaxis()->GetXmin();
            double x_max = h_elec_eff_dRjet->GetXaxis()->GetXmax();

            TH1F *frame = pad2->DrawFrame(x_min, 0.5, x_max, 1.49);
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetYaxis()->SetNdivisions(405);
            frame->SetLineWidth(1);
            frame->SetXTitle("#Delta R(e, jet)");
            frame->GetXaxis()->SetTitleSize(20);
            frame->GetXaxis()->SetTitleFont(47);
            frame->GetXaxis()->SetTitleOffset(3.0);
            frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetXaxis()->SetLabelSize(15);
            frame->SetYTitle("Data / MC");
            frame->GetYaxis()->SetTitleSize(17);
            frame->GetYaxis()->SetTitleFont(43);
            frame->GetYaxis()->SetTitleOffset(1.5);
            frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetYaxis()->SetLabelSize(16);
            frame->Draw();

            TH1F *ratio_TandP = (TH1F *)h_elec_eff_dRjet->Clone();
            ratio_TandP->Sumw2();
            ratio_TandP->Divide(h_Zee_eff_dRjet);
            ratio_TandP->SetMarkerColor(kRed);
            ratio_TandP->SetMarkerStyle(kFullSquare);
            ratio_TandP->SetLineColor(kRed);
            ratio_TandP->Draw("same");

            if (truth_match) {
                TH1F *ratio_truth = (TH1F *)h_elec_eff_dRjet->Clone();
                ratio_truth->Sumw2();
                ratio_truth->Divide(h_Zee_truth_eff_dRjet);
                ratio_truth->SetMarkerColor(kBlue);
                ratio_truth->SetMarkerStyle(kFullTriangleUp);
                ratio_truth->SetLineColor(kBlue);
                ratio_truth->Draw("same");
            }

            if (ttbar) {
                TH1F *ratio_ttbar = (TH1F *)h_elec_eff_dRjet->Clone();
                ratio_ttbar->Sumw2();
                ratio_ttbar->Divide(h_ttbar_elec_eff_dRjet);
                ratio_ttbar->SetMarkerColor(kMagenta);
                ratio_ttbar->SetMarkerStyle(kFullDiamond);
                ratio_ttbar->SetLineColor(kMagenta);
                ratio_ttbar->Draw("same");
            }

            if (Gtt) {
                TH1F *ratio_gtt = (TH1F *)h_elec_eff_dRjet->Clone();
                ratio_gtt->Sumw2();
                ratio_gtt->Divide(h_Gtt_elec_eff_dRjet);
                ratio_gtt->SetMarkerColor(kOrange);
                ratio_gtt->SetMarkerStyle(kFullCross);
                ratio_gtt->SetLineColor(kOrange);
                ratio_gtt->Draw("same");
            }
        }
        else if (lepton == "muon") {
            h_muon_eff_dRjet->SetMarkerColor(kBlack);
            h_muon_eff_dRjet->SetMarkerStyle(kFullCircle);
            h_muon_eff_dRjet->SetLineColor(kBlack);
            h_muon_eff_dRjet->SetMaximum(1.1);
            h_muon_eff_dRjet->SetMinimum(0);
            h_muon_eff_dRjet->SetTitle("");
            h_muon_eff_dRjet->SetXTitle("#Delta R(#mu, jet)");
            h_muon_eff_dRjet->SetYTitle("Real muon efficiency");
            h_muon_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
            h_muon_eff_dRjet->SetStats(kFALSE);
            h_muon_eff_dRjet->Draw();

            h_Zmumu_eff_dRjet->SetMarkerColor(kRed);
            h_Zmumu_eff_dRjet->SetMarkerStyle(kFullSquare);
            h_Zmumu_eff_dRjet->SetLineColor(kRed);
            h_Zmumu_eff_dRjet->Draw("same");

            if (truth_match) {
                h_Zmumu_truth_eff_dRjet->SetMarkerColor(kBlue);
                h_Zmumu_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
                h_Zmumu_truth_eff_dRjet->SetLineColor(kBlue);
                h_Zmumu_truth_eff_dRjet->Draw("same");
            }

            if (ttbar) {
                h_ttbar_muon_eff_dRjet->SetMarkerColor(kMagenta);
                h_ttbar_muon_eff_dRjet->SetMarkerStyle(kFullDiamond);
                h_ttbar_muon_eff_dRjet->SetLineColor(kMagenta);
                h_ttbar_muon_eff_dRjet->Draw("same");
            }

            if (Gtt) {
                h_Gtt_muon_eff_dRjet->SetMarkerColor(kOrange);
                h_Gtt_muon_eff_dRjet->SetMarkerStyle(kFullCross);
                h_Gtt_muon_eff_dRjet->SetLineColor(kOrange);
                h_Gtt_muon_eff_dRjet->Draw("same");
            }

            TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.05);
            leg->AddEntry("h_muon_eff_dRjet", "Data");
            leg->AddEntry("h_Zmumu_eff_dRjet", "Z #rightarrow #mu#mu T&P MC");
            if (truth_match) {
                leg->AddEntry("h_Zmumu_truth_eff_dRjet", "Z #rightarrow #mu#mu truth matched MC");
            }
            if (ttbar) {
                leg->AddEntry("h_ttbar_muon_eff_dRjet", "ttbar MC");
            }
            if (Gtt) {
                leg->AddEntry("h_Gtt_muon_eff_dRjet", "Gtt MC");
            }
            leg->Draw();

            pad2->cd(); // pad2 becomes the current pad

            double x_min = h_muon_eff_dRjet->GetXaxis()->GetXmin();
            double x_max = h_muon_eff_dRjet->GetXaxis()->GetXmax();

            TH1F *frame = pad2->DrawFrame(x_min, 0.5, x_max, 1.49);
            frame->GetXaxis()->SetNdivisions(510);
            frame->GetYaxis()->SetNdivisions(405);
            frame->SetLineWidth(1);
            frame->SetXTitle("#Delta R(#mu, jet)");
            frame->GetXaxis()->SetTitleSize(20);
            frame->GetXaxis()->SetTitleFont(47);
            frame->GetXaxis()->SetTitleOffset(3.0);
            frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetXaxis()->SetLabelSize(15);
            frame->SetYTitle("Data / MC");
            frame->GetYaxis()->SetTitleSize(17);
            frame->GetYaxis()->SetTitleFont(43);
            frame->GetYaxis()->SetTitleOffset(1.5);
            frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
            frame->GetYaxis()->SetLabelSize(16);
            frame->Draw();

            TH1F *ratio_TandP = (TH1F *)h_muon_eff_dRjet->Clone();
            ratio_TandP->Sumw2();
            ratio_TandP->Divide(h_Zmumu_eff_dRjet);
            ratio_TandP->SetMarkerColor(kRed);
            ratio_TandP->SetMarkerStyle(kFullSquare);
            ratio_TandP->SetLineColor(kRed);
            ratio_TandP->Draw("same");

            if (truth_match) {
                TH1F *ratio_truth = (TH1F *)h_muon_eff_dRjet->Clone();
                ratio_truth->Sumw2();
                ratio_truth->Divide(h_Zmumu_truth_eff_dRjet);
                ratio_truth->SetMarkerColor(kBlue);
                ratio_truth->SetMarkerStyle(kFullTriangleUp);
                ratio_truth->SetLineColor(kBlue);
                ratio_truth->Draw("same");
            }

            if (ttbar) {
                TH1F *ratio_ttbar = (TH1F *)h_muon_eff_dRjet->Clone();
                ratio_ttbar->Sumw2();
                ratio_ttbar->Divide(h_ttbar_muon_eff_dRjet);
                ratio_ttbar->SetMarkerColor(kMagenta);
                ratio_ttbar->SetMarkerStyle(kFullDiamond);
                ratio_ttbar->SetLineColor(kMagenta);
                ratio_ttbar->Draw("same");
            }

            if (Gtt) {
                TH1F *ratio_gtt = (TH1F *)h_muon_eff_dRjet->Clone();
                ratio_gtt->Sumw2();
                ratio_gtt->Divide(h_Gtt_muon_eff_dRjet);
                ratio_gtt->SetMarkerColor(kOrange);
                ratio_gtt->SetMarkerStyle(kFullCross);
                ratio_gtt->SetLineColor(kOrange);
                ratio_gtt->Draw("same");
            }
        }
    }

    string filename = "real_efficiency_ratio_plot_" + static_cast<string>(lepton) + "_" + static_cast<string>(plot) + ".pdf";
    real_efficiency_ratio_plot->SaveAs(filename.c_str(), "pdf");
}

void yt_real_efficiency_with_trigger_eta_range(bool isData, TString lepton, int eta_low_bin, int eta_up_bin) // Use 2-dim histograms
{
    TFile *file_single_lepton_trigger,
          *file_single_lepton_trigger_tag_trigger_matched,
          *file_dilepton_trigger,
          *file_dilepton_trigger_tag_trigger_matched;

    if (isData) {
        if (lepton == "electron") {
            file_single_lepton_trigger = TFile::Open(path + electron);
            file_single_lepton_trigger_tag_trigger_matched = TFile::Open(path + electron_tag_trigger_matched);
            file_dilepton_trigger = TFile::Open(path + electron_dilepton_trigger);
            file_dilepton_trigger_tag_trigger_matched = TFile::Open(path + electron_dilepton_trigger_tag_trigger_matched);
        }
        else if (lepton == "muon") {
            file_single_lepton_trigger = TFile::Open(path + muon);
            file_single_lepton_trigger_tag_trigger_matched = TFile::Open(path + muon_tag_trigger_matched);
            file_dilepton_trigger = TFile::Open(path + muon_dilepton_trigger);
            file_dilepton_trigger_tag_trigger_matched = TFile::Open(path + muon_dilepton_trigger_tag_trigger_matched);
        }
    }
    else { // isMC
        if (lepton == "electron") {
            file_single_lepton_trigger = TFile::Open(path + Zee_root);
            file_single_lepton_trigger_tag_trigger_matched = TFile::Open(path + Zee_tag_trigger_matched_root);
            file_dilepton_trigger = TFile::Open(path + Zee_dilepton_trigger_root);
            file_dilepton_trigger_tag_trigger_matched = TFile::Open(path + Zee_dilepton_trigger_tag_trigger_matched_root);
        }
        else if (lepton == "muon") {
            file_single_lepton_trigger = TFile::Open(path + Zmumu_root);
            file_single_lepton_trigger_tag_trigger_matched = TFile::Open(path + Zmumu_tag_trigger_matched_root);
            file_dilepton_trigger = TFile::Open(path + Zmumu_dilepton_trigger_root);
            file_dilepton_trigger_tag_trigger_matched = TFile::Open(path + Zmumu_dilepton_trigger_tag_trigger_matched_root);
        }
    } 

    // Naming:
    // sin: single_lepton_trigger
    // sin_tag: single_lepton_trigger_tag_trigger_matched
    // di: dilepton_trigge
    // di_sin: dilepton_trigger_but_fail_single_lepton_trigger
    // di_tag: dilepton_trigger_tag_trigger_matched
    // di_sin_tag: dilepton_trigger_but_fail_single_lepton_trigger_tag_trigger_matched

    // Get the 2-dim histograms
    TH2F *h_sin_baseline_pt_eta = (TH2F *)file_single_lepton_trigger->Get("h_baseline_pt_eta");
    TH2F *h_sin_signal_pt_eta = (TH2F *)file_single_lepton_trigger->Get("h_signal_pt_eta");

    TH2F *h_sin_tag_baseline_pt_eta = (TH2F *)file_single_lepton_trigger_tag_trigger_matched->Get("h_baseline_pt_eta");
    TH2F *h_sin_tag_signal_pt_eta = (TH2F *)file_single_lepton_trigger_tag_trigger_matched->Get("h_signal_pt_eta");

    TH2F *h_di_baseline_pt_eta = (TH2F *)file_dilepton_trigger->Get("h_baseline_pt_eta");
    TH2F *h_di_signal_pt_eta = (TH2F *)file_dilepton_trigger->Get("h_signal_pt_eta");

    TH2F *h_di_tag_baseline_pt_eta = (TH2F *)file_dilepton_trigger_tag_trigger_matched->Get("h_baseline_pt_eta");
    TH2F *h_di_tag_signal_pt_eta = (TH2F *)file_dilepton_trigger_tag_trigger_matched->Get("h_signal_pt_eta");

    // Project histograms to 1-dim pt axis
    TH1F *h_sin_baseline = (TH1F *)h_sin_baseline_pt_eta->ProjectionX("h_sin_baseline", eta_low_bin, eta_up_bin);
    TH1F *h_sin_signal = (TH1F *)h_sin_signal_pt_eta->ProjectionX("h_sin_signal", eta_low_bin, eta_up_bin);

    TH1F *h_sin_tag_baseline = (TH1F *)h_sin_tag_baseline_pt_eta->ProjectionX("h_sin_tag_baseline", eta_low_bin, eta_up_bin);
    TH1F *h_sin_tag_signal = (TH1F *)h_sin_tag_signal_pt_eta->ProjectionX("h_sin_tag_signal", eta_low_bin, eta_up_bin);

    TH1F *h_di_baseline = (TH1F *)h_di_baseline_pt_eta->ProjectionX("h_di_baseline", eta_low_bin, eta_up_bin);
    TH1F *h_di_signal = (TH1F *)h_di_signal_pt_eta->ProjectionX("h_di_signal", eta_low_bin, eta_up_bin);

    TH1F *h_di_tag_baseline = (TH1F *)h_di_tag_baseline_pt_eta->ProjectionX("h_di_tag_baseline", eta_low_bin, eta_up_bin);
    TH1F *h_di_tag_signal = (TH1F *)h_di_tag_signal_pt_eta->ProjectionX("h_di_tag_signal", eta_low_bin, eta_up_bin);

    // Calculate the efficiency
    TH1F *eff_sin = (TH1F *)h_sin_signal->Clone();
    TH1F *eff_sin_tag = (TH1F *)h_sin_tag_signal->Clone();
    TH1F *eff_di = (TH1F *)h_di_signal->Clone();
    TH1F *eff_di_tag = (TH1F *)h_di_tag_signal->Clone();

    eff_sin->Sumw2();
    eff_sin_tag->Sumw2();
    eff_di->Sumw2();
    eff_di_tag->Sumw2();

    eff_sin->Reset();
    eff_sin_tag->Reset();
    eff_di->Reset();
    eff_di_tag->Reset();

    eff_sin->SetName("eff_sin");
    eff_sin_tag->SetName("eff_sin_tag");
    eff_di->SetName("eff_di");
    eff_di_tag->SetName("eff_di_tag");

    eff_sin->Divide(h_sin_signal, h_sin_baseline, 1, 1, "B");
    eff_sin_tag->Divide(h_sin_tag_signal, h_sin_tag_baseline, 1, 1, "B");
    eff_di->Divide(h_di_signal, h_di_baseline, 1, 1, "B");
    eff_di_tag->Divide(h_di_tag_signal, h_di_tag_baseline, 1, 1, "B");

    double eta_low_value;
    double eta_up_value;
    eta_low_value = h_sin_baseline_pt_eta->GetYaxis()->GetBinLowEdge(eta_low_bin);
    eta_up_value = h_sin_baseline_pt_eta->GetYaxis()->GetBinUpEdge(eta_up_bin);
    // Convert double to string using stringstream
    stringstream sstream_eta_low, sstream_eta_up;
    sstream_eta_low << eta_low_value;
    sstream_eta_up << eta_up_value;

    TCanvas *trigger_uncertainty = new TCanvas("trigger_uncertainty", "Electron Real Efficiency", 500, 500);

    //Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetRightMargin(0.08);
    pad1->SetLeftMargin(0.12);
    //pad1->SetGridy(); // grid lines
    //pad1->SetLogx();
    pad1->Draw();

    // lower plot will be in pad
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetRightMargin(0.08);
    pad2->SetLeftMargin(0.12);
    pad2->SetGridy(); // grid lines
    //pad2->SetLogx();
    pad2->Draw();

    pad1->cd(); // pad1 becomes the current pad
    //pad1->SetFrameLineWidth(2);

    // Draw curve here
    vector<double> vec_max;
    double max = 0.;

    vec_max.push_back(eff_sin->GetMaximum());
    vec_max.push_back(eff_sin_tag->GetMaximum());
    vec_max.push_back(eff_di->GetMaximum());
    vec_max.push_back(eff_di_tag->GetMaximum());
    sort(vec_max.begin(), vec_max.end());
    max = vec_max.back();

    eff_sin->SetMarkerColor(kBlack);
    eff_sin->SetMarkerStyle(kFullCircle);
    eff_sin->SetLineColor(kBlack);
    eff_sin->SetTitle("");
    eff_sin->SetXTitle("p_{T} [GeV]");
    eff_sin->SetYTitle("Events");
    eff_sin->GetYaxis()->SetTitleOffset(1.5);
    eff_sin->SetMaximum(max * 1.1);
    if (lepton == "electron")
        eff_sin->SetMinimum(0.51);
    else if (lepton == "muon")
        eff_sin->SetMinimum(0.71);
    eff_sin->SetStats(kFALSE);
    eff_sin->Draw();

    eff_sin_tag->SetMarkerColor(kRed);
    eff_sin_tag->SetMarkerStyle(kFullSquare);
    eff_sin_tag->SetLineColor(kRed);
    eff_sin_tag->Draw("E1, same");

    eff_di->SetMarkerColor(kGreen);
    eff_di->SetMarkerStyle(kFullTriangleUp);
    eff_di->SetLineColor(kGreen);
    eff_di->Draw("E1, same");

    eff_di_tag->SetMarkerColor(kBlue);
    eff_di_tag->SetMarkerStyle(kFullTriangleDown);
    eff_di_tag->SetLineColor(kBlue);
    eff_di_tag->Draw("E1, same");

    TLegend *leg1 = new TLegend(0.2, 0.05, 0.8, 0.4);
    leg1->AddEntry(eff_sin, "Single lepton trigger", "lp");
    leg1->AddEntry(eff_sin_tag, "Single lepton trigger + tag trigger matched", "lp");
    leg1->AddEntry(eff_di, "Dilepton trigger", "lp");
    leg1->AddEntry(eff_di_tag, "Dilepton trigger + tag trigger matched", "lp");
    leg1->SetBorderSize(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");

    TPaveText *text1 = new TPaveText(0.3, 0.45, 0.4, 0.55, "NDC");
    text1->SetTextSize(0.05);
    text1->SetBorderSize(0);
    text1->SetFillStyle(0);
    text1->SetFillColor(0);
    text1->AddText(string(sstream_eta_low.str() + " < |#eta| < " + sstream_eta_up.str()).c_str());
    text1->Draw("same");

    pad2->cd(); // pad2 becomes the current pad

    TH1F *frame;
    if (lepton == "electron")
        frame = pad2->DrawFrame(10, 0.95, 200, 1.049);
    if (lepton == "muon")
        frame = pad2->DrawFrame(10, 0.98, 200, 1.039);
    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetNdivisions(405);
    frame->SetLineWidth(1);
    frame->SetXTitle("p_{T} [GeV]");
    frame->GetXaxis()->SetTitleSize(20);
    frame->GetXaxis()->SetTitleFont(47);
    frame->GetXaxis()->SetTitleOffset(3.0);
    frame->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame->GetXaxis()->SetLabelSize(15);
    frame->SetYTitle("Ratio");
    frame->GetYaxis()->SetTitleSize(17);
    frame->GetYaxis()->SetTitleFont(43);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame->GetYaxis()->SetLabelSize(16);
    frame->Draw();

    TH1F *ratio_sin = (TH1F *)eff_sin->Clone();
    // TH1F *ratio_sin_tag = (TH1F *)eff_sin_tag->Clone();
    TH1F *ratio_di = (TH1F *)eff_di->Clone();
    TH1F *ratio_di_tag = (TH1F *)eff_di_tag->Clone();

    ratio_sin->Sumw2();
    // ratio_sin_tag->Sumw2();
    ratio_di->Sumw2();
    ratio_di_tag->Sumw2();

    ratio_sin->Divide(eff_sin_tag);
    ratio_di->Divide(eff_sin_tag);
    ratio_di_tag->Divide(eff_sin_tag);
    // ratio_sin_tag->Divide(eff_sin);
    // ratio_di->Divide(eff_sin);
    // ratio_di_tag->Divide(eff_sin);

    ratio_sin->SetMarkerColor(kBlack);
    ratio_sin->SetMarkerStyle(kFullCircle);
    ratio_sin->SetLineColor(kBlack);
    ratio_sin->Draw("same");
    // ratio_sin_tag->SetMarkerColor(kRed);
    // ratio_sin_tag->SetMarkerStyle(kFullSquare);
    // ratio_sin_tag->SetLineColor(kRed);
    // ratio_sin_tag->Draw("same");

    ratio_di->SetMarkerColor(kGreen);
    ratio_di->SetMarkerStyle(kFullTriangleUp);
    ratio_di->SetLineColor(kGreen);
    ratio_di->Draw("same");

    ratio_di_tag->SetMarkerColor(kBlue);
    ratio_di_tag->SetMarkerStyle(kFullTriangleDown);
    ratio_di_tag->SetLineColor(kBlue);
    ratio_di_tag->Draw("same");

    // Reset stringstream
    sstream_eta_low.str("");
    sstream_eta_low.clear();
    sstream_eta_up.str("");
    sstream_eta_up.clear();

    eta_low_value *= 100.;
    eta_up_value *= 100.;
    sstream_eta_low << eta_low_value;
    sstream_eta_up << eta_up_value;

    string filename = "trigger_uncertainty_eta" + sstream_eta_low.str() + sstream_eta_up.str() + "_ratio_plot.pdf";
    if (lepton == "electron")
        filename = "trigger_uncertainty_electron_eta" + sstream_eta_low.str() + sstream_eta_up.str() + "_ratio_plot.pdf";
    else if (lepton == "muon")
        filename = "trigger_uncertainty_muon_eta" + sstream_eta_low.str() + sstream_eta_up.str() + "_ratio_plot.pdf";
    trigger_uncertainty->SaveAs(filename.c_str());
}

void yt_real_efficiency_function_of_eta_with_different_pt()
{
    TFile *data_elec = TFile::Open(path + electron);

    TH2F *h_baseline_pt_eta = (TH2F *)data_elec->Get("h_baseline_pt_eta");
    TH2F *h_signal_pt_eta = (TH2F *)data_elec->Get("h_signal_pt_eta");

    int pt_range_1_bin_low = h_baseline_pt_eta->GetXaxis()->FindBin(10+0.1);
    int pt_range_1_bin_up = h_baseline_pt_eta->GetXaxis()->FindBin(40-0.1);
    int pt_range_2_bin_low = h_baseline_pt_eta->GetXaxis()->FindBin(40+0.1);
    int pt_range_2_bin_up = h_baseline_pt_eta->GetXaxis()->FindBin(80-0.1);
    int pt_range_3_bin_low = h_baseline_pt_eta->GetXaxis()->FindBin(80+0.1);
    int pt_range_3_bin_up = h_baseline_pt_eta->GetXaxis()->FindBin(200-0.1);

    // cout << "pt_range_1_bin_low=" << pt_range_1_bin_low << ", pt_range_1_bin_up=" << pt_range_1_bin_up << endl;
    // cout << "pt_range_2_bin_low=" << pt_range_2_bin_low << ", pt_range_2_bin_up=" << pt_range_2_bin_up << endl;
    // cout << "pt_range_3_bin_low=" << pt_range_3_bin_low << ", pt_range_3_bin_up=" << pt_range_3_bin_up << endl;

    TH1F *h_baseline_pt_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionY("h_baseline_pt_range_1", pt_range_1_bin_low, pt_range_1_bin_up); // 10 < pT < 40
    TH1F *h_baseline_pt_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionY("h_baseline_pt_range_2", pt_range_2_bin_low, pt_range_2_bin_up); // 40 < pT < 80
    TH1F *h_baseline_pt_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionY("h_baseline_pt_range_3", pt_range_3_bin_low, pt_range_3_bin_up); // 80 < pT < 200

    TH1F *h_signal_pt_range_1 = (TH1F *)h_signal_pt_eta->ProjectionY("h_signal_pt_range_1", pt_range_1_bin_low, pt_range_1_bin_up); // 10 < pT < 40
    TH1F *h_signal_pt_range_2 = (TH1F *)h_signal_pt_eta->ProjectionY("h_signal_pt_range_2", pt_range_2_bin_low, pt_range_2_bin_up); // 40 < pT < 80
    TH1F *h_signal_pt_range_3 = (TH1F *)h_signal_pt_eta->ProjectionY("h_signal_pt_range_3", pt_range_3_bin_low, pt_range_3_bin_up); // 80 < pT < 200

    TH1F *h_ratio_1 = (TH1F *)h_signal_pt_range_1->Clone();
    TH1F *h_ratio_2 = (TH1F *)h_signal_pt_range_2->Clone();
    TH1F *h_ratio_3 = (TH1F *)h_signal_pt_range_3->Clone();

    h_ratio_1->Sumw2();
    h_ratio_2->Sumw2();
    h_ratio_3->Sumw2();

    h_ratio_1->Reset();
    h_ratio_2->Reset();
    h_ratio_3->Reset();

    h_ratio_1->SetName("ratio_1");
    h_ratio_2->SetName("ratio_2");
    h_ratio_3->SetName("ratio_3");

    h_ratio_1->Divide(h_signal_pt_range_1, h_baseline_pt_range_1, 1, 1, "B");
    h_ratio_2->Divide(h_signal_pt_range_2, h_baseline_pt_range_2, 1, 1, "B");
    h_ratio_3->Divide(h_signal_pt_range_3, h_baseline_pt_range_3, 1, 1, "B");

    TCanvas *real_efficiency_plot = new TCanvas("real_efficiency_plot", "Electron Real Efficiency as function of |#eta|", 500, 500);
    real_efficiency_plot->SetLeftMargin(0.12);
    real_efficiency_plot->cd();

    h_ratio_1->SetLineColor(kBlack);
    h_ratio_1->SetMarkerColor(kBlack);
    h_ratio_1->SetMarkerStyle(kFullCircle);
    h_ratio_1->SetXTitle("|#eta|");
    h_ratio_1->SetYTitle("Electron Real Efficiency");
    h_ratio_1->GetXaxis()->SetRangeUser(0, 2.);
    h_ratio_1->GetYaxis()->SetTitleOffset(1.5);
    h_ratio_1->SetStats(kFALSE);
    h_ratio_1->SetMaximum(1.1);
    h_ratio_1->SetMinimum(0.5);
    h_ratio_1->SetTitle("");
    h_ratio_1->Draw();

    h_ratio_2->SetLineColor(kRed);
    h_ratio_2->SetMarkerColor(kRed);
    h_ratio_2->SetMarkerStyle(kFullSquare);
    h_ratio_2->Draw("same");

    h_ratio_3->SetLineColor(kBlue);
    h_ratio_3->SetMarkerColor(kBlue);
    h_ratio_3->SetMarkerStyle(kFullTriangleUp);
    h_ratio_3->Draw("same");

    // float X_max = h_ratio_1->GetXaxis()->GetXmax();
    float X_min = h_ratio_1->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., 2., 1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.2, 0.2, 0.6, 0.4);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry("ratio_1", "10 < p_{T} < 40");
    leg->AddEntry("ratio_2", "40 < p_{T} < 80");
    leg->AddEntry("ratio_3", "80 < p_{T} < 200");
    leg->Draw();

    real_efficiency_plot->SaveAs("real_efficiency_electron_function_of_eta_with_different_pt.pdf", "pdf");
}


void yt_Gtt_electron_real_efficiency_individual_cut_study()
{
    TFile *data_elec = TFile::Open(path + electron);
    TFile *Zee_TandP = TFile::Open(path + Zee_root);
    TFile *Zee_truth = TFile::Open(path + Zee_truth_root);

    TFile *Gtt_elec_all_signal = TFile::Open(path + Gtt_elec_all_signal_root);
    TFile *Gtt_elec_Medium_only = TFile::Open(path + Gtt_elec_Medium_only_root);
    TFile *Gtt_elec_track_iso_only = TFile::Open(path + Gtt_elec_track_iso_only_root);
    TFile *Gtt_elec_calo_iso_only = TFile::Open(path + Gtt_elec_calo_iso_only_root);
    TFile *Gtt_elec_z0_only = TFile::Open(path + Gtt_elec_z0_only_root);
    TFile *Gtt_elec_eta_only = TFile::Open(path + Gtt_elec_eta_only_root);
    TFile *Gtt_elec_all_cuts = TFile::Open(path + Gtt_elec_root);

    TH1F *h_elec_eff_dRjet = (TH1F *)data_elec->Get("h_eff_dRjet");
    TH1F *h_Zee_eff_dRjet = (TH1F *)Zee_TandP->Get("h_eff_dRjet");
    TH1F *h_Zee_truth_eff_dRjet = (TH1F *)Zee_truth->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_all_signal_dRjet = (TH1F *)Gtt_elec_all_signal->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_Medium_only_dRjet = (TH1F *)Gtt_elec_Medium_only->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_track_iso_only_dRjet = (TH1F *)Gtt_elec_track_iso_only->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_calo_iso_only_dRjet =  (TH1F *)Gtt_elec_calo_iso_only->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_z0_only_dRjet = (TH1F *)Gtt_elec_z0_only->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_eta_only_dRjet = (TH1F *)Gtt_elec_eta_only->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_all_cuts_dRjet = (TH1F *)Gtt_elec_all_cuts->Get("h_eff_dRjet");

    h_elec_eff_dRjet->SetName("h_elec_eff_dRjet");
    h_Zee_eff_dRjet->SetName("h_Zee_eff_dRjet");
    h_Zee_truth_eff_dRjet->SetName("h_Zee_truth_eff_dRjet");
    h_Gtt_elec_all_signal_dRjet->SetName("h_Gtt_elec_all_signal_dRjet");
    h_Gtt_elec_Medium_only_dRjet->SetName("h_Gtt_elec_Medium_only_dRjet");
    h_Gtt_elec_track_iso_only_dRjet->SetName("h_Gtt_elec_track_iso_only_dRjet");
    h_Gtt_elec_calo_iso_only_dRjet->SetName("h_Gtt_elec_calo_iso_only_dRjet");
    h_Gtt_elec_z0_only_dRjet->SetName("h_Gtt_elec_z0_only_dRjet");
    h_Gtt_elec_eta_only_dRjet->SetName("h_Gtt_elec_eta_only_dRjet");
    h_Gtt_elec_all_cuts_dRjet->SetName("h_Gtt_elec_all_cuts_dRjet");

    TCanvas *Gtt_electron_real_efficiency_study = new TCanvas("real_efficiency_plot", "Real Efficiency", 600, 600);
    Gtt_electron_real_efficiency_study->SetLeftMargin(0.12);

    h_elec_eff_dRjet->SetMarkerColor(kBlack);
    h_elec_eff_dRjet->SetMarkerStyle(kFullCircle);
    h_elec_eff_dRjet->SetLineColor(kBlack);
    h_elec_eff_dRjet->SetMaximum(1.1);
    h_elec_eff_dRjet->SetMinimum(0);
    h_elec_eff_dRjet->SetTitle("");
    h_elec_eff_dRjet->SetXTitle("#Delta R(e, jet)");
    h_elec_eff_dRjet->SetYTitle("Real electron efficiency");
    h_elec_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
    h_elec_eff_dRjet->SetStats(kFALSE);
    h_elec_eff_dRjet->Draw();

    h_Zee_eff_dRjet->SetMarkerColor(kRed);
    h_Zee_eff_dRjet->SetMarkerStyle(kFullSquare);
    h_Zee_eff_dRjet->SetLineColor(kRed);
    h_Zee_eff_dRjet->Draw("same");

    h_Zee_truth_eff_dRjet->SetMarkerColor(kBlue);
    h_Zee_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
    h_Zee_truth_eff_dRjet->SetLineColor(kBlue);
    h_Zee_truth_eff_dRjet->Draw("same");

    h_Gtt_elec_all_cuts_dRjet->SetMarkerColor(kOrange);
    h_Gtt_elec_all_cuts_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_all_cuts_dRjet->SetLineColor(kOrange);
    h_Gtt_elec_all_cuts_dRjet->Draw("same");

    h_Gtt_elec_all_signal_dRjet->SetMarkerColor(kOrange-3);
    h_Gtt_elec_all_signal_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_all_signal_dRjet->SetLineColor(kOrange-3);
    h_Gtt_elec_all_signal_dRjet->Draw("same");

    h_Gtt_elec_Medium_only_dRjet->SetMarkerColor(kOrange-6);
    h_Gtt_elec_Medium_only_dRjet->SetMarkerStyle(kOpenCross);
    h_Gtt_elec_Medium_only_dRjet->SetLineColor(kOrange-6);
    h_Gtt_elec_Medium_only_dRjet->Draw("same");

    h_Gtt_elec_track_iso_only_dRjet->SetMarkerColor(kOrange+3);
    h_Gtt_elec_track_iso_only_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_track_iso_only_dRjet->SetLineColor(kOrange+3);
    h_Gtt_elec_track_iso_only_dRjet->Draw("same");

    h_Gtt_elec_calo_iso_only_dRjet->SetMarkerColor(kOrange+6);
    h_Gtt_elec_calo_iso_only_dRjet->SetMarkerStyle(kOpenCross);
    h_Gtt_elec_calo_iso_only_dRjet->SetLineColor(kOrange+6);
    h_Gtt_elec_calo_iso_only_dRjet->Draw("same");

    h_Gtt_elec_z0_only_dRjet->SetMarkerColor(kOrange+9);
    h_Gtt_elec_z0_only_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_z0_only_dRjet->SetLineColor(kOrange+9);
    h_Gtt_elec_z0_only_dRjet->Draw("same");

    h_Gtt_elec_eta_only_dRjet->SetMarkerColor(kGreen);
    h_Gtt_elec_eta_only_dRjet->SetMarkerStyle(kOpenCross);
    h_Gtt_elec_eta_only_dRjet->SetLineColor(kGreen);
    h_Gtt_elec_eta_only_dRjet->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.2, 0.7, 0.5);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->AddEntry("h_elec_eff_dRjet", "Data");
    leg->AddEntry("h_Zee_eff_dRjet", "Z #rightarrow ee T&P MC");
    leg->AddEntry("h_Zee_truth_eff_dRjet", "Z #rightarrow ee truth matched MC");
    leg->AddEntry("h_Gtt_elec_all_cuts_dRjet", "Gtt MC (all cuts)");
    leg->AddEntry("h_Gtt_elec_all_signal_dRjet", "Gtt MC (all are signal)");
    leg->AddEntry("h_Gtt_elec_Medium_only_dRjet", "Gtt MC (Medium only)");
    leg->AddEntry("h_Gtt_elec_track_iso_only_dRjet", "Gtt MC (Track iso only)");
    leg->AddEntry("h_Gtt_elec_calo_iso_only_dRjet", "Gtt MC (Calo iso only)");
    leg->AddEntry("h_Gtt_elec_z0_only_dRjet", "Gtt MC (z0 cut only)");
    leg->AddEntry("h_Gtt_elec_eta_only_dRjet", "Gtt MC (#eta cut only)");
    leg->Draw();

    Gtt_electron_real_efficiency_study->SaveAs("real_efficiency_Gtt_elec_individual_cut_study.pdf", "pdf");
}

void yt_Gtt_electron_real_efficiency_combine_cuts_study()
{
    TFile *data_elec = TFile::Open(path + electron);
    TFile *Zee_TandP = TFile::Open(path + Zee_root);
    TFile *Zee_truth = TFile::Open(path + Zee_truth_root);
  
    TFile *Gtt_elec_Medium_and_calo = TFile::Open(path + "submitDir_MC_GG_ttn1_electron/hist-1025_Medium_and_calo.root");
    TFile *Gtt_elec_track_and_z0 = TFile::Open(path + "submitDir_MC_GG_ttn1_electron/hist-1025_track_and_z0.root");
    TFile *Gtt_elec_except_eta = TFile::Open(path + "submitDir_MC_GG_ttn1_electron/hist-1025_except_eta.root");
    TFile *Gtt_elec_eta_only = TFile::Open(path + "submitDir_MC_GG_ttn1_electron/hist-1025_eta_only.root");
    TFile *Gtt_elec_all_cuts = TFile::Open(path + Gtt_elec_root);

    TH1F *h_elec_eff_dRjet = (TH1F *)data_elec->Get("h_eff_dRjet");
    TH1F *h_Zee_eff_dRjet = (TH1F *)Zee_TandP->Get("h_eff_dRjet");
    TH1F *h_Zee_truth_eff_dRjet = (TH1F *)Zee_truth->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_Medium_and_calo_dRjet = (TH1F *)Gtt_elec_Medium_and_calo->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_track_and_z0_dRjet = (TH1F *)Gtt_elec_track_and_z0->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_except_eta_dRjet =  (TH1F *)Gtt_elec_except_eta->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_eta_only_dRjet = (TH1F *)Gtt_elec_eta_only->Get("h_eff_dRjet");
    TH1F *h_Gtt_elec_all_cuts_dRjet = (TH1F *)Gtt_elec_all_cuts->Get("h_eff_dRjet");

    h_elec_eff_dRjet->SetName("h_elec_eff_dRjet");
    h_Zee_eff_dRjet->SetName("h_Zee_eff_dRjet");
    h_Zee_truth_eff_dRjet->SetName("h_Zee_truth_eff_dRjet");
    h_Gtt_elec_Medium_and_calo_dRjet->SetName("h_Gtt_elec_Medium_and_calo_dRjet");
    h_Gtt_elec_track_and_z0_dRjet->SetName("h_Gtt_elec_track_and_z0_dRjet");
    h_Gtt_elec_except_eta_dRjet->SetName("h_Gtt_elec_except_eta_dRjet");
    h_Gtt_elec_eta_only_dRjet->SetName("h_Gtt_elec_eta_only_dRjet");
    h_Gtt_elec_all_cuts_dRjet->SetName("h_Gtt_elec_all_cuts_dRjet");

    TCanvas *Gtt_electron_real_efficiency_study = new TCanvas("real_efficiency_plot", "Real Efficiency", 600, 600);
    Gtt_electron_real_efficiency_study->SetLeftMargin(0.12);

    h_elec_eff_dRjet->SetMarkerColor(kBlack);
    h_elec_eff_dRjet->SetMarkerStyle(kFullCircle);
    h_elec_eff_dRjet->SetLineColor(kBlack);
    h_elec_eff_dRjet->SetMaximum(1.1);
    h_elec_eff_dRjet->SetMinimum(0);
    h_elec_eff_dRjet->SetTitle("");
    h_elec_eff_dRjet->SetXTitle("#Delta R(e, jet)");
    h_elec_eff_dRjet->SetYTitle("Real electron efficiency");
    h_elec_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
    h_elec_eff_dRjet->SetStats(kFALSE);
    h_elec_eff_dRjet->Draw();

    h_Zee_eff_dRjet->SetMarkerColor(kRed);
    h_Zee_eff_dRjet->SetMarkerStyle(kFullSquare);
    h_Zee_eff_dRjet->SetLineColor(kRed);
    h_Zee_eff_dRjet->Draw("same");

    h_Zee_truth_eff_dRjet->SetMarkerColor(kBlue);
    h_Zee_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
    h_Zee_truth_eff_dRjet->SetLineColor(kBlue);
    h_Zee_truth_eff_dRjet->Draw("same");

    h_Gtt_elec_all_cuts_dRjet->SetMarkerColor(kOrange);
    h_Gtt_elec_all_cuts_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_all_cuts_dRjet->SetLineColor(kOrange);
    h_Gtt_elec_all_cuts_dRjet->Draw("same");
/*
    h_Gtt_elec_Medium_and_calo_dRjet->SetMarkerColor(kOrange-3);
    h_Gtt_elec_Medium_and_calo_dRjet->SetMarkerStyle(kOpenCross);
    h_Gtt_elec_Medium_and_calo_dRjet->SetLineColor(kOrange-3);
    h_Gtt_elec_Medium_and_calo_dRjet->Draw("same");

    h_Gtt_elec_track_and_z0_dRjet->SetMarkerColor(kOrange-6);
    h_Gtt_elec_track_and_z0_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_track_and_z0_dRjet->SetLineColor(kOrange-6);
    h_Gtt_elec_track_and_z0_dRjet->Draw("same");
*/
    h_Gtt_elec_except_eta_dRjet->SetMarkerColor(kOrange+3);
    h_Gtt_elec_except_eta_dRjet->SetMarkerStyle(kOpenCross);
    h_Gtt_elec_except_eta_dRjet->SetLineColor(kOrange+3);
    h_Gtt_elec_except_eta_dRjet->Draw("same");

    h_Gtt_elec_eta_only_dRjet->SetMarkerColor(kGreen);
    h_Gtt_elec_eta_only_dRjet->SetMarkerStyle(kOpenCross);
    h_Gtt_elec_eta_only_dRjet->SetLineColor(kGreen);
    h_Gtt_elec_eta_only_dRjet->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.2, 0.7, 0.5);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->AddEntry("h_elec_eff_dRjet", "Data");
    leg->AddEntry("h_Zee_eff_dRjet", "Z #rightarrow ee T&P MC");
    leg->AddEntry("h_Zee_truth_eff_dRjet", "Z #rightarrow ee truth matched MC");
    leg->AddEntry("h_Gtt_elec_all_cuts_dRjet", "Gtt MC (all cuts)");
    //leg->AddEntry("h_Gtt_elec_Medium_and_calo_dRjet", "Gtt MC (Medium & calo)");
    //leg->AddEntry("h_Gtt_elec_track_and_z0_dRjet", "Gtt MC (Track & z_{0} only)");
    leg->AddEntry("h_Gtt_elec_except_eta_dRjet", "Gtt MC (except #eta)");
    leg->AddEntry("h_Gtt_elec_eta_only_dRjet", "Gtt MC (#eta cut only)");
    leg->Draw();

    Gtt_electron_real_efficiency_study->SaveAs("real_efficiency_Gtt_elec_combine_cuts_study.pdf", "pdf");
}

void yt_Gtt_electron_real_efficiency_comparison()
{
    TFile *Gtt_elec = TFile::Open(path + Gtt_elec_root);

    TH3F *h_Gtt_elec_baseline_pt_eta_dRjet = (TH3F *)Gtt_elec->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_Gtt_elec_signal_pt_eta_dRjet = (TH3F *)Gtt_elec->Get("h_signal_pt_eta_dRjet");

    // For electron case, we only use |eta| < 2.0 to calculate the efficiency
    int eta_low_bin = h_Gtt_elec_baseline_pt_eta_dRjet->GetYaxis()->FindBin(0+0.01);
    int eta_up_bin = h_Gtt_elec_baseline_pt_eta_dRjet->GetYaxis()->FindBin(2-0.01);
    
    TH1F *h_Gtt_elec_baseline_pt = (TH1F *)h_Gtt_elec_baseline_pt_eta_dRjet->ProjectionX("", eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Gtt_elec_baseline_eta = (TH1F *)h_Gtt_elec_baseline_pt_eta_dRjet->ProjectionY("")->Clone();
    TH1F *h_Gtt_elec_baseline_dRjet = (TH1F *)h_Gtt_elec_baseline_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();

    TH1F *h_Gtt_elec_signal_pt = (TH1F *)h_Gtt_elec_signal_pt_eta_dRjet->ProjectionX("", eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Gtt_elec_signal_eta = (TH1F *)h_Gtt_elec_signal_pt_eta_dRjet->ProjectionY("")->Clone();
    TH1F *h_Gtt_elec_signal_dRjet = (TH1F *)h_Gtt_elec_signal_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();

    TH1F *h_Gtt_eff_pt = (TH1F *)h_Gtt_elec_signal_pt->Clone();
    TH1F *h_Gtt_eff_eta = (TH1F *)h_Gtt_elec_signal_eta->Clone();
    TH1F *h_Gtt_eff_dRjet = (TH1F *)h_Gtt_elec_signal_dRjet->Clone();

    h_Gtt_eff_pt->Reset();
    h_Gtt_eff_eta->Reset();
    h_Gtt_eff_dRjet->Reset();

    h_Gtt_eff_pt->SetName("h_Gtt_eff_pt");
    h_Gtt_eff_eta->SetName("h_Gtt_eff_eta");
    h_Gtt_eff_dRjet->SetName("h_Gtt_eff_dRjet");

    h_Gtt_eff_pt->Divide(h_Gtt_elec_signal_pt, h_Gtt_elec_baseline_pt, 1, 1, "B");
    h_Gtt_eff_eta->Divide(h_Gtt_elec_signal_eta, h_Gtt_elec_baseline_eta, 1, 1, "B");
    h_Gtt_eff_dRjet->Divide(h_Gtt_elec_signal_dRjet, h_Gtt_elec_baseline_dRjet, 1, 1, "B");

    TH1F *h_Gtt_elec_eff_pt = (TH1F *)Gtt_elec->Get("h_eff_pt");
    TH1F *h_Gtt_elec_eff_eta = (TH1F *)Gtt_elec->Get("h_eff_eta");
    TH1F *h_Gtt_elec_eff_dRjet = (TH1F *)Gtt_elec->Get("h_eff_dRjet");

    h_Gtt_elec_eff_pt->SetName("h_Gtt_elec_eff_pt");
    h_Gtt_elec_eff_eta->SetName("h_Gtt_elec_eff_eta");
    h_Gtt_elec_eff_dRjet->SetName("h_Gtt_elec_eff_dRjet");

    TCanvas *real_efficiency = new TCanvas("real_efficiency_plot", "Real Efficiency", 1200, 400);
    real_efficiency->Divide(3, 1);

    real_efficiency->cd(1);
    gPad->SetLeftMargin(0.12);
    h_Gtt_elec_eff_pt->SetMarkerColor(kOrange);
    h_Gtt_elec_eff_pt->SetMarkerStyle(kFullCross);
    h_Gtt_elec_eff_pt->SetLineColor(kOrange);
    h_Gtt_elec_eff_pt->SetMaximum(1.1);
    h_Gtt_elec_eff_pt->SetMinimum(0);
    h_Gtt_elec_eff_pt->SetTitle("");
    h_Gtt_elec_eff_pt->SetXTitle("p_{T} [GeV]");
    h_Gtt_elec_eff_pt->SetYTitle("Real electron efficiency");
    h_Gtt_elec_eff_pt->GetYaxis()->SetTitleOffset(1.5);
    h_Gtt_elec_eff_pt->SetStats(kFALSE);
    h_Gtt_elec_eff_pt->Draw();

    h_Gtt_eff_pt->SetMarkerColor(kOrange+3);
    h_Gtt_eff_pt->SetMarkerStyle(kOpenStar);
    h_Gtt_eff_pt->SetLineColor(kOrange+3);
    h_Gtt_eff_pt->Draw("same");

    TLegend *leg1 = new TLegend(0.3, 0.2, 0.7, 0.5);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.04);
    leg1->AddEntry("h_Gtt_elec_eff_pt", "Gtt MC");
    leg1->AddEntry("h_Gtt_eff_pt", "Gtt MC (0 < |#eta| < 2.0)");
    leg1->Draw();

    real_efficiency->cd(2);
    gPad->SetLeftMargin(0.12);
    h_Gtt_elec_eff_eta->SetMarkerColor(kOrange);
    h_Gtt_elec_eff_eta->SetMarkerStyle(kFullCross);
    h_Gtt_elec_eff_eta->SetLineColor(kOrange);
    h_Gtt_elec_eff_eta->SetMaximum(1.1);
    h_Gtt_elec_eff_eta->SetMinimum(0);
    h_Gtt_elec_eff_eta->GetXaxis()->SetRangeUser(0, 2.0);
    h_Gtt_elec_eff_eta->SetTitle("");
    h_Gtt_elec_eff_eta->SetXTitle("|#eta|");
    h_Gtt_elec_eff_eta->SetYTitle("Real electron efficiency");
    h_Gtt_elec_eff_eta->GetYaxis()->SetTitleOffset(1.5);
    h_Gtt_elec_eff_eta->SetStats(kFALSE);
    h_Gtt_elec_eff_eta->Draw();

    h_Gtt_eff_eta->SetMarkerColor(kOrange+3);
    h_Gtt_eff_eta->SetMarkerStyle(kOpenStar);
    h_Gtt_eff_eta->SetLineColor(kOrange+3);
    h_Gtt_eff_eta->Draw("same");

    TLegend *leg2 = new TLegend(0.3, 0.2, 0.7, 0.5);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry("h_Gtt_elec_eff_eta", "Gtt MC");
    leg2->AddEntry("h_Gtt_eff_eta", "Gtt MC (0 < |#eta| < 2.0)");
    leg2->Draw();

    real_efficiency->cd(3);
    gPad->SetLeftMargin(0.12);
    h_Gtt_elec_eff_dRjet->SetMarkerColor(kOrange);
    h_Gtt_elec_eff_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_eff_dRjet->SetLineColor(kOrange);
    h_Gtt_elec_eff_dRjet->SetMaximum(1.1);
    h_Gtt_elec_eff_dRjet->SetMinimum(0);
    h_Gtt_elec_eff_dRjet->SetTitle("");
    h_Gtt_elec_eff_dRjet->SetXTitle("#Delta R(e, jet)");
    h_Gtt_elec_eff_dRjet->SetYTitle("Real electron efficiency");
    h_Gtt_elec_eff_dRjet->GetYaxis()->SetTitleOffset(1.5);
    h_Gtt_elec_eff_dRjet->SetStats(kFALSE);
    h_Gtt_elec_eff_dRjet->Draw();

    h_Gtt_eff_dRjet->SetMarkerColor(kOrange+3);
    h_Gtt_eff_dRjet->SetMarkerStyle(kOpenStar);
    h_Gtt_eff_dRjet->SetLineColor(kOrange+3);
    h_Gtt_eff_dRjet->Draw("same");

    TLegend *leg3 = new TLegend(0.3, 0.2, 0.7, 0.5);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetTextSize(0.04);
    leg3->AddEntry("h_Gtt_elec_eff_dRjet", "Gtt MC");
    leg3->AddEntry("h_Gtt_eff_dRjet", "Gtt MC (0 < |#eta| < 2.0)");
    leg3->Draw();

    real_efficiency->SaveAs("real_efficiency_Gtt_elec_comparison.pdf", "pdf");
}

void yt_real_efficiency_vs_dR_electron()
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zee_truth = TFile::Open(path + Zee_truth_root);
    TFile *ttbar_elec = TFile::Open(path + ttbar_elec_root);
    TFile *Gtt_elec = TFile::Open(path + Gtt_elec_root);

    TH3F *h_elec_baseline_pt_eta_dRjet = (TH3F *)data_elec->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_Zee_pt_baseline_eta_dRjet = (TH3F *)Zee_TandP->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_Zee_truth_baseline_pt_eta_dRjet = (TH3F *)Zee_truth->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_ttbar_elec_baseline_pt_eta_dRjet = (TH3F *)ttbar_elec->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_Gtt_elec_baseline_pt_eta_dRjet = (TH3F *)Gtt_elec->Get("h_baseline_pt_eta_dRjet");

    TH3F *h_elec_signal_pt_eta_dRjet = (TH3F *)data_elec->Get("h_signal_pt_eta_dRjet");
    TH3F *h_Zee_pt_signal_eta_dRjet = (TH3F *)Zee_TandP->Get("h_signal_pt_eta_dRjet");
    TH3F *h_Zee_truth_signal_pt_eta_dRjet = (TH3F *)Zee_truth->Get("h_signal_pt_eta_dRjet");
    TH3F *h_ttbar_elec_signal_pt_eta_dRjet = (TH3F *)ttbar_elec->Get("h_signal_pt_eta_dRjet");
    TH3F *h_Gtt_elec_signal_pt_eta_dRjet = (TH3F *)Gtt_elec->Get("h_signal_pt_eta_dRjet");

    // For electron case, we only use |eta| < 2.0 to calculate the efficiency
    int eta_low_bin = h_elec_baseline_pt_eta_dRjet->GetYaxis()->FindBin(0+0.01);
    int eta_up_bin = h_elec_baseline_pt_eta_dRjet->GetYaxis()->FindBin(2-0.01);

    TH1F *h_elec_baseline_dRjet = (TH1F *)h_elec_baseline_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zee_baseline_dRjet = (TH1F *)h_Zee_pt_baseline_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zee_truth_baseline_dRjet = (TH1F *)h_Zee_truth_baseline_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_ttbar_elec_baseline_dRjet = (TH1F *)h_ttbar_elec_baseline_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Gtt_elec_baseline_dRjet = (TH1F *)h_Gtt_elec_baseline_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();

    TH1F *h_elec_signal_dRjet = (TH1F *)h_elec_signal_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zee_signal_dRjet = (TH1F *)h_Zee_pt_signal_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zee_truth_signal_dRjet = (TH1F *)h_Zee_truth_signal_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_ttbar_elec_signal_dRjet = (TH1F *)h_ttbar_elec_signal_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Gtt_elec_signal_dRjet = (TH1F *)h_Gtt_elec_signal_pt_eta_dRjet->ProjectionZ("", 0, -1, eta_low_bin, eta_up_bin)->Clone();

    TH1F *h_elec_eff_dRjet = (TH1F *)h_elec_signal_dRjet->Clone();
    TH1F *h_Zee_eff_dRjet = (TH1F *)h_Zee_signal_dRjet->Clone();
    TH1F *h_Zee_truth_eff_dRjet = (TH1F *)h_Zee_truth_signal_dRjet->Clone();
    TH1F *h_ttbar_elec_eff_dRjet = (TH1F *)h_ttbar_elec_signal_dRjet->Clone();
    TH1F *h_Gtt_elec_eff_dRjet = (TH1F *)h_Gtt_elec_signal_dRjet->Clone();

    h_elec_eff_dRjet->Reset();
    h_Zee_eff_dRjet->Reset();
    h_Zee_truth_eff_dRjet->Reset();
    h_ttbar_elec_eff_dRjet->Reset();
    h_Gtt_elec_eff_dRjet->Reset();

    h_elec_eff_dRjet->SetName("h_elec_eff_dRjet");
    h_Zee_eff_dRjet->SetName("h_Zee_eff_dRjet");
    h_Zee_truth_eff_dRjet->SetName("h_Zee_truth_eff_dRjet");
    h_ttbar_elec_eff_dRjet->SetName("h_ttbar_elec_eff_dRjet");
    h_Gtt_elec_eff_dRjet->SetName("h_Gtt_elec_eff_dRjet");

    h_elec_eff_dRjet->Divide(h_elec_signal_dRjet, h_elec_baseline_dRjet, 1, 1, "B");
    h_Zee_eff_dRjet->Divide(h_Zee_signal_dRjet, h_Zee_baseline_dRjet, 1, 1, "B");
    h_Zee_truth_eff_dRjet->Divide(h_Zee_truth_signal_dRjet, h_Zee_truth_baseline_dRjet, 1, 1, "B");
    h_ttbar_elec_eff_dRjet->Divide(h_ttbar_elec_signal_dRjet, h_ttbar_elec_baseline_dRjet, 1, 1, "B");
    h_Gtt_elec_eff_dRjet->Divide(h_Gtt_elec_signal_dRjet, h_Gtt_elec_baseline_dRjet, 1, 1, "B");

    TCanvas *real_efficiency_elec = new TCanvas("real_efficiency_plot", "Real Efficiency", 600, 600);
    real_efficiency_elec->SetLeftMargin(0.12);

    h_elec_eff_dRjet->SetMarkerColor(kBlack);
    h_elec_eff_dRjet->SetMarkerStyle(kFullCircle);
    h_elec_eff_dRjet->SetLineColor(kBlack);
    h_elec_eff_dRjet->SetMaximum(1.1);
    h_elec_eff_dRjet->SetMinimum(0);
    h_elec_eff_dRjet->SetTitle("");
    h_elec_eff_dRjet->SetXTitle("#Delta R(e, jet)");
    h_elec_eff_dRjet->SetYTitle("Real electron efficiency");
    h_elec_eff_dRjet->SetTitleOffset(1.5);
    h_elec_eff_dRjet->SetStats(kFALSE);
    h_elec_eff_dRjet->Draw();

    h_Zee_eff_dRjet->SetMarkerColor(kRed);
    h_Zee_eff_dRjet->SetMarkerStyle(kFullSquare);
    h_Zee_eff_dRjet->SetLineColor(kRed);
    h_Zee_eff_dRjet->Draw("same");

    h_Zee_truth_eff_dRjet->SetMarkerColor(kBlue);
    h_Zee_truth_eff_dRjet->SetMarkerStyle(kFullTriangleUp);
    h_Zee_truth_eff_dRjet->SetLineColor(kBlue);
    h_Zee_truth_eff_dRjet->Draw("same");
    
    h_ttbar_elec_eff_dRjet->SetMarkerColor(kMagenta);
    h_ttbar_elec_eff_dRjet->SetMarkerStyle(kFullDiamond);
    h_ttbar_elec_eff_dRjet->SetLineColor(kMagenta);
    h_ttbar_elec_eff_dRjet->Draw("same");

    h_Gtt_elec_eff_dRjet->SetMarkerColor(kOrange);
    h_Gtt_elec_eff_dRjet->SetMarkerStyle(kFullCross);
    h_Gtt_elec_eff_dRjet->SetLineColor(kOrange);
    h_Gtt_elec_eff_dRjet->Draw("same");

    TLegend *leg1 = new TLegend(0.3, 0.2, 0.7, 0.5);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.04);
    leg1->AddEntry("h_elec_eff_dRjet", "Data");
    leg1->AddEntry("h_Zee_eff_dRjet", "Z #rightarrow ee T&P MC");
    leg1->AddEntry("h_Zee_truth_eff_dRjet", "Z #rightarrow ee truth matched MC");
    leg1->AddEntry("h_ttbar_elec_eff_dRjet", "ttbar MC (0 < |#eta| < 2.0)");
    leg1->AddEntry("h_Gtt_elec_eff_dRjet", "Gtt MC (0 < |#eta| < 2.0)");
    leg1->Draw();

    real_efficiency_elec->SaveAs("real_efficiency_vs_dR_elec.pdf", "pdf");
}

void yt_real_efficiency_vs_AvgMu()
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    // TFile *Zee_TandP = TFile::Open(path + Zee_root);
    // TFile *Zmumu_TandP = TFile::Open(path + Zmumu_root);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);
    TFile *Zee_TandP = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *Zmumu_TandP = TFile::Open(path + Zmumu_tag_trigger_matched_root);
    TFile *Zee_truth = TFile::Open(path + Zee_truth_root);
    TFile *Zmumu_truth = TFile::Open(path + Zmumu_truth_root);
    TFile *ttbar_elec = TFile::Open(path + ttbar_elec_root);
    TFile *ttbar_muon = TFile::Open(path + ttbar_muon_root);
    TFile *Gtt_elec = TFile::Open(path + Gtt_elec_root);
    TFile *Gtt_muon = TFile::Open(path + Gtt_muon_root);

    TH3F *h_elec_baseline_pt_eta_AvgMu = (TH3F *)data_elec->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_elec_signal_pt_eta_AvgMu = (TH3F *)data_elec->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_muon_baseline_pt_eta_AvgMu = (TH3F *)data_muon->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_muon_signal_pt_eta_AvgMu = (TH3F *)data_muon->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_Zee_TandP_baseline_pt_eta_AvgMu = (TH3F *)Zee_TandP->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_Zee_TandP_signal_pt_eta_AvgMu = (TH3F *)Zee_TandP->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_Zmumu_TandP_baseline_pt_eta_AvgMu = (TH3F *)Zmumu_TandP->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_Zmumu_TandP_signal_pt_eta_AvgMu = (TH3F *)Zmumu_TandP->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_Zee_truth_baseline_pt_eta_AvgMu = (TH3F *)Zee_truth->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_Zee_truth_signal_pt_eta_AvgMu = (TH3F *)Zee_truth->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_Zmumu_truth_baseline_pt_eta_AvgMu = (TH3F *)Zmumu_truth->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_Zmumu_truth_signal_pt_eta_AvgMu = (TH3F *)Zmumu_truth->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_ttbar_elec_baseline_pt_eta_AvgMu = (TH3F *)ttbar_elec->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_ttbar_elec_signal_pt_eta_AvgMu = (TH3F *)ttbar_elec->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_ttbar_muon_baseline_pt_eta_AvgMu = (TH3F *)ttbar_muon->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_ttbar_muon_signal_pt_eta_AvgMu = (TH3F *)ttbar_muon->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_Gtt_elec_baseline_pt_eta_AvgMu = (TH3F *)Gtt_elec->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_Gtt_elec_signal_pt_eta_AvgMu = (TH3F *)Gtt_elec->Get("h_signal_pt_eta_AvgMu");
    TH3F *h_Gtt_muon_baseline_pt_eta_AvgMu = (TH3F *)Gtt_muon->Get("h_baseline_pt_eta_AvgMu");
    TH3F *h_Gtt_muon_signal_pt_eta_AvgMu = (TH3F *)Gtt_muon->Get("h_signal_pt_eta_AvgMu");

    // Select pT region to check
    int pt_low_bin = 0;
    int pt_up_bin = -1;

    // For electron case, we only use |eta| < 2.0 to calculate the efficiency
    int eta_low_bin = h_elec_baseline_pt_eta_AvgMu->GetYaxis()->FindBin(0+0.01);
    int eta_up_bin = h_elec_baseline_pt_eta_AvgMu->GetYaxis()->FindBin(2-0.01);

    TH1F *h_elec_baseline_AvgMu = (TH1F *)h_elec_baseline_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_elec_signal_AvgMu = (TH1F *)h_elec_signal_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_muon_baseline_AvgMu = (TH1F *)h_muon_baseline_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_muon_signal_AvgMu = (TH1F *)h_muon_signal_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_Zee_TandP_baseline_AvgMu = (TH1F *)h_Zee_TandP_baseline_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zee_TandP_signal_AvgMu = (TH1F *)h_Zee_TandP_signal_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zmumu_TandP_baseline_AvgMu = (TH1F *)h_Zmumu_TandP_baseline_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_Zmumu_TandP_signal_AvgMu = (TH1F *)h_Zmumu_TandP_signal_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_Zee_truth_baseline_AvgMu = (TH1F *)h_Zee_truth_baseline_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zee_truth_signal_AvgMu = (TH1F *)h_Zee_truth_signal_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Zmumu_truth_baseline_AvgMu = (TH1F *)h_Zmumu_truth_baseline_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_Zmumu_truth_signal_AvgMu = (TH1F *)h_Zmumu_truth_signal_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_ttbar_elec_baseline_AvgMu = (TH1F *)h_ttbar_elec_baseline_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_ttbar_elec_signal_AvgMu = (TH1F *)h_ttbar_elec_signal_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_ttbar_muon_baseline_AvgMu = (TH1F *)h_ttbar_muon_baseline_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_ttbar_muon_signal_AvgMu = (TH1F *)h_ttbar_muon_signal_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_Gtt_elec_baseline_AvgMu = (TH1F *)h_Gtt_elec_baseline_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Gtt_elec_signal_AvgMu = (TH1F *)h_Gtt_elec_signal_pt_eta_AvgMu->ProjectionZ("", pt_low_bin, pt_up_bin, eta_low_bin, eta_up_bin)->Clone();
    TH1F *h_Gtt_muon_baseline_AvgMu = (TH1F *)h_Gtt_muon_baseline_pt_eta_AvgMu->ProjectionZ("")->Clone();
    TH1F *h_Gtt_muon_signal_AvgMu = (TH1F *)h_Gtt_muon_signal_pt_eta_AvgMu->ProjectionZ("")->Clone();

    h_elec_baseline_AvgMu->SetName("h_elec_baseline_AvgMu");
    h_elec_signal_AvgMu->SetName("h_elec_signal_AvgMu");
    h_muon_baseline_AvgMu->SetName("h_muon_baseline_AvgMu");
    h_muon_signal_AvgMu->SetName("h_muon_signal_AvgMu");
    h_Zee_TandP_baseline_AvgMu->SetName("h_Zee_TandP_baseline_AvgMu");
    h_Zee_TandP_signal_AvgMu->SetName("h_Zee_TandP_signal_AvgMu");
    h_Zmumu_TandP_baseline_AvgMu->SetName("h_Zmumu_TandP_baseline_AvgMu");
    h_Zmumu_TandP_signal_AvgMu->SetName("h_Zmumu_TandP_signal_AvgMu");
    h_Zee_truth_baseline_AvgMu->SetName("h_Zee_truth_baseline_AvgMu");
    h_Zee_truth_signal_AvgMu->SetName("h_Zee_truth_signal_AvgMu");
    h_Zmumu_truth_baseline_AvgMu->SetName("h_Zmumu_truth_baseline_AvgMu");
    h_Zmumu_truth_signal_AvgMu->SetName("h_Zmumu_truth_signal_AvgMu");
    h_ttbar_elec_baseline_AvgMu->SetName("h_ttbar_elec_baseline_AvgMu");
    h_ttbar_elec_signal_AvgMu->SetName("h_ttbar_elec_signal_AvgMu");
    h_ttbar_muon_baseline_AvgMu->SetName("h_ttbar_muon_baseline_AvgMu");
    h_ttbar_muon_signal_AvgMu->SetName("h_ttbar_muon_signal_AvgMu");
    h_Gtt_elec_baseline_AvgMu->SetName("h_Gtt_elec_baseline_AvgMu");
    h_Gtt_elec_signal_AvgMu->SetName("h_Gtt_elec_signal_AvgMu");
    h_Gtt_muon_baseline_AvgMu->SetName("h_Gtt_muon_baseline_AvgMu");
    h_Gtt_muon_signal_AvgMu->SetName("h_Gtt_muon_signal_AvgMu");

    TH1F *h_elec_eff_Avg = (TH1F *)h_elec_baseline_AvgMu->Clone();
    TH1F *h_muon_eff_Avg = (TH1F *)h_muon_baseline_AvgMu->Clone();
    TH1F *h_Zee_TandP_eff_Avg = (TH1F *)h_Zee_TandP_baseline_AvgMu->Clone();
    TH1F *h_Zmumu_TandP_eff_Avg = (TH1F *)h_Zmumu_TandP_baseline_AvgMu->Clone();
    TH1F *h_Zee_truth_eff_Avg = (TH1F *)h_Zee_truth_baseline_AvgMu->Clone();
    TH1F *h_Zmumu_truth_eff_Avg = (TH1F *)h_Zmumu_truth_baseline_AvgMu->Clone();
    TH1F *h_ttbar_elec_eff_Avg = (TH1F *)h_ttbar_elec_baseline_AvgMu->Clone();
    TH1F *h_ttbar_muon_eff_Avg = (TH1F *)h_ttbar_muon_baseline_AvgMu->Clone();
    TH1F *h_Gtt_elec_eff_Avg = (TH1F *)h_Gtt_elec_baseline_AvgMu->Clone();
    TH1F *h_Gtt_muon_eff_Avg = (TH1F *)h_Gtt_muon_baseline_AvgMu->Clone();

    // Reset the contents of the histograms
    h_elec_eff_Avg->Reset();
    h_muon_eff_Avg->Reset();
    h_Zee_TandP_eff_Avg->Reset();
    h_Zmumu_TandP_eff_Avg->Reset();
    h_Zee_truth_eff_Avg->Reset();
    h_Zmumu_truth_eff_Avg->Reset();
    h_ttbar_elec_eff_Avg->Reset();
    h_ttbar_muon_eff_Avg->Reset();
    h_Gtt_elec_eff_Avg->Reset();
    h_Gtt_muon_eff_Avg->Reset();

    h_elec_eff_Avg->SetName("h_elec_eff_Avg");
    h_muon_eff_Avg->SetName("h_muon_eff_Avg");
    h_Zee_TandP_eff_Avg->SetName("h_Zee_TandP_eff_Avg");
    h_Zmumu_TandP_eff_Avg->SetName("h_Zmumu_TandP_eff_Avg");
    h_Zee_truth_eff_Avg->SetName("h_Zee_truth_eff_Avg");
    h_Zmumu_truth_eff_Avg->SetName("h_Zmumu_truth_eff_Avg");
    h_ttbar_elec_eff_Avg->SetName("h_ttbar_elec_eff_Avg");
    h_ttbar_muon_eff_Avg->SetName("h_ttbar_muon_eff_Avg");
    h_Gtt_elec_eff_Avg->SetName("h_Gtt_elec_eff_Avg");
    h_Gtt_muon_eff_Avg->SetName("h_Gtt_muon_eff_Avg");

    h_elec_eff_Avg->Divide(h_elec_signal_AvgMu, h_elec_baseline_AvgMu, 1, 1, "B");
    h_muon_eff_Avg->Divide(h_muon_signal_AvgMu, h_muon_baseline_AvgMu, 1, 1, "B");
    h_Zee_TandP_eff_Avg->Divide(h_Zee_TandP_signal_AvgMu, h_Zee_TandP_baseline_AvgMu, 1, 1, "B");
    h_Zmumu_TandP_eff_Avg->Divide(h_Zmumu_TandP_signal_AvgMu, h_Zmumu_TandP_baseline_AvgMu, 1, 1, "B");
    h_Zee_truth_eff_Avg->Divide(h_Zee_truth_signal_AvgMu, h_Zee_truth_baseline_AvgMu, 1, 1, "B");
    h_Zmumu_truth_eff_Avg->Divide(h_Zmumu_truth_signal_AvgMu, h_Zmumu_truth_baseline_AvgMu, 1, 1, "B");
    h_ttbar_elec_eff_Avg->Divide(h_ttbar_elec_signal_AvgMu, h_ttbar_elec_baseline_AvgMu, 1, 1, "B");
    h_ttbar_muon_eff_Avg->Divide(h_ttbar_muon_signal_AvgMu, h_ttbar_muon_baseline_AvgMu, 1, 1, "B");
    h_Gtt_elec_eff_Avg->Divide(h_Gtt_elec_signal_AvgMu, h_Gtt_elec_baseline_AvgMu, 1, 1, "B");
    h_Gtt_muon_eff_Avg->Divide(h_Gtt_muon_signal_AvgMu, h_Gtt_muon_baseline_AvgMu, 1, 1, "B");

    TCanvas *real_efficiency_elec = new TCanvas("real_efficiency_elec", "Real Efficiency", 600, 600);
    real_efficiency_elec->SetLeftMargin(0.12);

    h_elec_eff_Avg->SetMarkerColor(kBlack);
    h_elec_eff_Avg->SetMarkerStyle(kFullCircle);
    h_elec_eff_Avg->SetLineColor(kBlack);
    h_elec_eff_Avg->SetMaximum(1.1);
    h_elec_eff_Avg->SetMinimum(0.7);
    h_elec_eff_Avg->SetTitle("");
    h_elec_eff_Avg->SetXTitle("<#mu>");
    h_elec_eff_Avg->SetYTitle("Real electron efficiency");
    h_elec_eff_Avg->GetYaxis()->SetTitleOffset(1.5);
    h_elec_eff_Avg->SetStats(kFALSE);
    h_elec_eff_Avg->Draw();

    h_Zee_TandP_eff_Avg->SetMarkerColor(kRed);
    h_Zee_TandP_eff_Avg->SetMarkerStyle(kFullSquare);
    h_Zee_TandP_eff_Avg->SetLineColor(kRed);
    h_Zee_TandP_eff_Avg->Draw("same");

    h_Zee_truth_eff_Avg->SetMarkerColor(kBlue);
    h_Zee_truth_eff_Avg->SetMarkerStyle(kFullTriangleUp);
    h_Zee_truth_eff_Avg->SetLineColor(kBlue);
    h_Zee_truth_eff_Avg->Draw("same");

    h_ttbar_elec_eff_Avg->SetMarkerColor(kMagenta);
    h_ttbar_elec_eff_Avg->SetMarkerStyle(kFullDiamond);
    h_ttbar_elec_eff_Avg->SetLineColor(kMagenta);
    h_ttbar_elec_eff_Avg->Draw("same");

    h_Gtt_elec_eff_Avg->SetMarkerColor(kOrange);
    h_Gtt_elec_eff_Avg->SetMarkerStyle(kFullCross);
    h_Gtt_elec_eff_Avg->SetLineColor(kOrange);
    h_Gtt_elec_eff_Avg->Draw("same");

    TLegend *leg1 = new TLegend(0.5, 0.65, 0.7, 0.89);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry("h_elec_eff_Avg", "Data");
    leg1->AddEntry("h_Zee_TandP_eff_Avg", "Z #rightarrow ee T&P MC");
    leg1->AddEntry("h_Zee_truth_eff_Avg", "Z #rightarrow ee truth matched MC");
    leg1->AddEntry("h_ttbar_elec_eff_Avg", "ttbar MC (0 < |#eta| < 2.0)");
    leg1->AddEntry("h_Gtt_elec_eff_Avg", "Gtt MC (0 < |#eta| < 2.0)");
    leg1->Draw();

    real_efficiency_elec->SaveAs("real_efficiency_vs_AvgMu_elec.pdf", "pdf");

    TCanvas *real_efficiency_muon = new TCanvas("real_efficiency_muon", "Real Efficiency", 600, 600);
    real_efficiency_muon->SetLeftMargin(0.12);

    h_muon_eff_Avg->SetMarkerColor(kBlack);
    h_muon_eff_Avg->SetMarkerStyle(kFullCircle);
    h_muon_eff_Avg->SetLineColor(kBlack);
    h_muon_eff_Avg->SetMaximum(1.1);
    h_muon_eff_Avg->SetMinimum(0.7);
    h_muon_eff_Avg->SetTitle("");
    h_muon_eff_Avg->SetXTitle("<#mu>");
    h_muon_eff_Avg->SetYTitle("Real muon efficiency");
    h_muon_eff_Avg->GetYaxis()->SetTitleOffset(1.5);
    h_muon_eff_Avg->SetStats(kFALSE);
    h_muon_eff_Avg->Draw();

    h_Zmumu_TandP_eff_Avg->SetMarkerColor(kRed);
    h_Zmumu_TandP_eff_Avg->SetMarkerStyle(kFullSquare);
    h_Zmumu_TandP_eff_Avg->SetLineColor(kRed);
    h_Zmumu_TandP_eff_Avg->Draw("same");

    h_Zmumu_truth_eff_Avg->SetMarkerColor(kBlue);
    h_Zmumu_truth_eff_Avg->SetMarkerStyle(kFullTriangleUp);
    h_Zmumu_truth_eff_Avg->SetLineColor(kBlue);
    h_Zmumu_truth_eff_Avg->Draw("same");

    h_ttbar_muon_eff_Avg->SetMarkerColor(kMagenta);
    h_ttbar_muon_eff_Avg->SetMarkerStyle(kFullDiamond);
    h_ttbar_muon_eff_Avg->SetLineColor(kMagenta);
    h_ttbar_muon_eff_Avg->Draw("same");

    h_Gtt_muon_eff_Avg->SetMarkerColor(kOrange);
    h_Gtt_muon_eff_Avg->SetMarkerStyle(kFullCross);
    h_Gtt_muon_eff_Avg->SetLineColor(kOrange);
    h_Gtt_muon_eff_Avg->Draw("same");

    TLegend *leg2 = new TLegend(0.5, 0.65, 0.7, 0.89);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry("h_muon_eff_Avg", "Data");
    leg2->AddEntry("h_Zmumu_TandP_eff_Avg", "Z #rightarrow #mu#mu T&P MC");
    leg2->AddEntry("h_Zmumu_truth_eff_Avg", "Z #rightarrow #mu#mu truth matched MC");
    leg2->AddEntry("h_ttbar_muon_eff_Avg", "ttbar MC");
    leg2->AddEntry("h_Gtt_muon_eff_Avg", "Gtt MC)");
    leg2->Draw();

    real_efficiency_muon->SaveAs("real_efficiency_vs_AvgMu_muon.pdf", "pdf");
}

// lepton = electron, muon
// x: pt, eta, dRjet
void yt_truth_match_TandP_comparison(TString lepton, TString x, bool TandP_truth_match = false)
{
    TFile *TandP_file;
    TFile *truth_match_file;
    TFile *TandP_truth_match_file;
    if (lepton == "electron") {
        // TandP_file = TFile::Open(path + Zee_root);
        TandP_file = TFile::Open(path + Zee_tag_trigger_matched_root);
        truth_match_file = TFile::Open(path + Zee_truth_root);
        if (TandP_truth_match)
            TandP_truth_match_file = TFile::Open(path + Zee_TP_truth_root);
    }
    else if (lepton == "muon") {
        // TandP_file = TFile::Open(path + Zmumu_root);
        TandP_file = TFile::Open(path + Zmumu_tag_trigger_matched_root);
        truth_match_file = TFile::Open(path + Zmumu_truth_root);
        if (TandP_truth_match)
            TandP_truth_match_file = TFile::Open(path + Zmumu_TP_truth_root);
    }

    TH1F *h_TandP_eff_pt = (TH1F *)TandP_file->Get("h_eff_pt");
    TH1F *h_TandP_eff_eta = (TH1F *)TandP_file->Get("h_eff_eta");
    TH1F *h_TandP_eff_dRjet = (TH1F *)TandP_file->Get("h_eff_dRjet");

    TH1F *h_truth_match_eff_pt = (TH1F *)truth_match_file->Get("h_eff_pt");
    TH1F *h_truth_match_eff_eta = (TH1F *)truth_match_file->Get("h_eff_eta");
    TH1F *h_truth_match_eff_dRjet = (TH1F *)truth_match_file->Get("h_eff_dRjet");

    if (TandP_truth_match) {
        TH1F *h_TandP_truth_match_eff_pt = (TH1F *)TandP_truth_match_file->Get("h_eff_pt");
        TH1F *h_TandP_truth_match_eff_eta = (TH1F *)TandP_truth_match_file->Get("h_eff_eta");
        TH1F *h_TandP_truth_match_eff_dRjet = (TH1F *)TandP_truth_match_file->Get("h_eff_dRjet");
    }

    TH1F *h_TandP_eff;
    TH1F *h_truth_match_eff;
    TH1F *h_TandP_truth_match_eff;
    TString x_title;
    float xlow, xup;

    if (x == "pt") {
        h_TandP_eff = (TH1F *)h_TandP_eff_pt->Clone();
        h_truth_match_eff = (TH1F *)h_truth_match_eff_pt->Clone();
        if (TandP_truth_match)
            h_TandP_truth_match_eff = (TH1F *)h_truth_match_eff_pt->Clone();
        x_title = "p_{T} [GeV]";
        xlow = 10.;
        xup = 200.;
    }
    else if (x == "eta") {
        h_TandP_eff = (TH1F *)h_TandP_eff_eta->Clone();
        h_truth_match_eff = (TH1F *)h_truth_match_eff_eta->Clone();
        if (TandP_truth_match)
            h_TandP_truth_match_eff = (TH1F *)h_truth_match_eff_eta->Clone();
        x_title = "|#eta|";
        xlow = 0.;
        if (lepton == "electron")
            xup = 2.0; // electron
        if (lepton == "muon")
            xup = 2.5; // muon
    }
    else if (x == "dRjet") {
        h_TandP_eff = (TH1F *)h_TandP_eff_dRjet->Clone();
        h_truth_match_eff = (TH1F *)h_truth_match_eff_dRjet->Clone();
        if (TandP_truth_match)
            h_TandP_truth_match_eff = (TH1F *)h_truth_match_eff_dRjet->Clone();
        if (lepton == "electron")
            x_title = "#Delta R(e, jet)";
        if (lepton == "muon")
            x_title = "#Delta R(#mu, jet)";
        xlow = 0.;
        xup = 4.;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

    //Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetRightMargin(0.08);
    pad1->SetLeftMargin(0.12);
    //pad1->SetGridy(); // grid lines
    //pad1->SetLogx();
    pad1->Draw();

    // lower plot will be in pad
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetRightMargin(0.08);
    pad2->SetLeftMargin(0.12);
    pad2->SetGridy(); // grid lines
    //pad2->SetLogx();
    pad2->Draw();

    pad1->cd(); // pad1 becomes the current pad
    //pad1->SetFrameLineWidth(2);

    // Draw curve here

    h_TandP_eff->SetMarkerColor(kRed);
    h_TandP_eff->SetMarkerStyle(kFullCircle);
    h_TandP_eff->SetLineColor(kRed);
    h_TandP_eff->SetMaximum(1.1);
    h_TandP_eff->SetMinimum(0.01);
    h_TandP_eff->SetTitle("");
    h_TandP_eff->SetXTitle(x_title);
    h_TandP_eff->SetYTitle("Real lepton efficiency");
    h_TandP_eff->GetXaxis()->SetRangeUser(xlow, xup);
    h_TandP_eff->GetYaxis()->SetTitleOffset(1.5);
    h_TandP_eff->SetStats(kFALSE);
    h_TandP_eff->Draw();

    h_truth_match_eff->SetMarkerColor(kBlue);
    h_truth_match_eff->SetMarkerStyle(kFullTriangleUp);
    h_truth_match_eff->SetLineColor(kBlue);
    h_truth_match_eff->Draw("same");

    if (TandP_truth_match) {
        h_TandP_truth_match_eff->SetMarkerColor(kGreen);
        h_TandP_truth_match_eff->SetMarkerStyle(kFullTriangleDown);
        h_TandP_truth_match_eff->SetLineColor(kGreen);
        h_TandP_truth_match_eff->Draw("same");
    }

    // Add legend
    TPaveText *text = new TPaveText(0.5, 0.5, 0.7, 0.6, "NDC");
    if (lepton == "electron")
        text->AddText("Z#rightarrow ee");
    else if (lepton == "muon")
        text->AddText("Z#rightarrow #mu#mu");
    text->SetTextSize(0.04);
    text->Draw("same");
    TLegend *leg = new TLegend(0.5, 0.2, 0.8, 0.5);
    leg->AddEntry(h_TandP_eff, "T&P MC", "lp");
    leg->AddEntry(h_truth_match_eff, "truth match MC", "lp");
    if (TandP_truth_match)
        leg->AddEntry(h_TandP_truth_match_eff, "TandP + truth match MC", "lp");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    pad2->cd(); // pad2 becomes the current pad

    TH1F *frame_left;
    if (lepton == "electron")
        frame_left = pad2->DrawFrame(xlow, 0.90, xup, 1.099); // electron
    if (lepton == "muon")
        frame_left = pad2->DrawFrame(xlow, 0.98, xup, 1.019); // muon
    frame_left->GetXaxis()->SetNdivisions(510);
    frame_left->GetYaxis()->SetNdivisions(405);
    frame_left->SetLineWidth(1);
    frame_left->SetXTitle(x_title);
    frame_left->GetXaxis()->SetTitleSize(20);
    frame_left->GetXaxis()->SetTitleFont(47);
    frame_left->GetXaxis()->SetTitleOffset(3.0);
    frame_left->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_left->GetXaxis()->SetLabelSize(12);
    frame_left->SetYTitle("ratio");
    frame_left->GetYaxis()->SetTitleSize(17);
    frame_left->GetYaxis()->SetTitleFont(43);
    frame_left->GetYaxis()->SetTitleOffset(1.5);
    frame_left->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame_left->GetYaxis()->SetLabelSize(12);
    frame_left->Draw();

    h_TandP_eff->Sumw2();
    h_truth_match_eff->Sumw2();
    if (TandP_truth_match)
        h_TandP_truth_match_eff->Sumw2();

    TH1F *ratio = (TH1F *)h_truth_match_eff->Clone();
    ratio->Divide(h_TandP_eff);
    ratio->SetMarkerColor(kBlue);
    ratio->SetMarkerStyle(kFullTriangleUp);
    ratio->SetLineColor(kBlue);
    ratio->Draw("e0, same");

    if (TandP_truth_match) {
        TH1F *ratio2 = (TH1F *)h_TandP_truth_match_eff->Clone();
        ratio2->Divide(h_TandP_eff);
        ratio2->SetMarkerColor(kGreen);
        ratio2->SetMarkerStyle(kFullTriangleDown);
        ratio2->SetLineColor(kGreen);
        ratio2->Draw("e0, same");
    }

    c1->SaveAs("Compare_TandP_truth_match_" + lepton + "_" + x + ".pdf", "pdf");
}

void yt_truth_match_TandP_comparison_table(TString lepton, int eta_bin_low, int eta_bin_up)
{
    TFile *TandP_file;
    TFile *truth_match_file;

    if (lepton == "electron") {
        // TandP_file = TFile::Open(path + Zee_root);
        TandP_file = TFile::Open(path + Zee_tag_trigger_matched_root);
        truth_match_file = TFile::Open(path + Zee_truth_root);
    }
    else if (lepton == "muon") {
        // TandP_file = TFile::Open(path + Zmumu_root);
        TandP_file = TFile::Open(path + Zmumu_tag_trigger_matched_root);
        truth_match_file = TFile::Open(path + Zmumu_truth_root);
    }

    TH3F *h_TandP_baseline_pt_eta_dRjet = (TH3F *)TandP_file->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_truth_match_baseline_pt_eta_dRjet = (TH3F *)truth_match_file->Get("h_baseline_pt_eta_dRjet");
    TH3F *h_TandP_signal_pt_eta_dRjet = (TH3F *)TandP_file->Get("h_signal_pt_eta_dRjet");
    TH3F *h_truth_match_signal_pt_eta_dRjet = (TH3F *)truth_match_file->Get("h_signal_pt_eta_dRjet");

    TH1F *h_TandP_baseline = (TH1F *)h_TandP_baseline_pt_eta_dRjet->ProjectionX("", eta_bin_low, eta_bin_up)->Clone();
    TH1F *h_truth_match_baseline = (TH1F *)h_truth_match_baseline_pt_eta_dRjet->ProjectionX("", eta_bin_low, eta_bin_up)->Clone();
    TH1F *h_TandP_signal = (TH1F *)h_TandP_signal_pt_eta_dRjet->ProjectionX("", eta_bin_low, eta_bin_up)->Clone();
    TH1F *h_truth_match_signal = (TH1F *)h_truth_match_signal_pt_eta_dRjet->ProjectionX("", eta_bin_low, eta_bin_up)->Clone();

    TH1F *h_TandP_eff = (TH1F *)h_TandP_signal->Clone();
    TH1F *h_truth_match_eff = (TH1F *)h_truth_match_signal->Clone();

    h_TandP_eff->Reset();
    h_truth_match_eff->Reset();

    h_TandP_eff->Sumw2();
    h_truth_match_eff->Sumw2();

    h_TandP_eff->Divide(h_TandP_signal, h_TandP_baseline, 1, 1, "B");
    h_truth_match_eff->Divide(h_truth_match_signal, h_truth_match_baseline, 1, 1, "B");

    // h_TandP_eff->SetLineColor(kBlue);
    // h_TandP_eff->Draw();
    // h_truth_match_eff->SetLineColor(kRed);
    // h_truth_match_eff->Draw("same")

    TH1F *ratio = (TH1F *)h_truth_match_eff->Clone();
    ratio->Divide(h_TandP_eff);

    // ratio->Draw();

    for(int i = 1; i <= ratio->GetXaxis()->GetNbins(); i++) {
        // do not include underflow and overflow bins
        cout << ratio->GetBinContent(i) << endl;
    }
}

//
// Background plots
//

void yt_background_subtraction_illustration()
{
    // TFile *data_elec = TFile::Open(path + electron);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TH3F *data_elec_baseline_pt_eta_mll = (TH3F *)data_elec->Get("h_baseline_pt_eta_mll");
    TH1D *data_elec_mll = (TH1D *)data_elec_baseline_pt_eta_mll->ProjectionZ("data_elec_mll");

    double bin_mll_low = data_elec_mll->FindBin(80. + 0.01);
    double bin_mll_up = data_elec_mll->FindBin(100. - 0.01);
    double bin_tail_low = data_elec_mll->FindBin(120. + 0.01);
    double bin_tail_up = data_elec_mll->FindBin(150. - 0.01);

    TH1D *h_mll_window = (TH1D *)data_elec_mll->Clone();
    h_mll_window->Reset();
    TH1D *h_tail = (TH1D *)data_elec_mll->Clone();
    h_tail->Reset();

    for (int i = 0; i < data_elec_mll->GetXaxis()->GetNbins() ; i++) {
        if (i >= bin_mll_low && i <= bin_mll_up)
            h_mll_window->SetBinContent(i, data_elec_mll->GetBinContent(i));
        if (i >= bin_tail_low && i <= bin_tail_up)
            h_tail->SetBinContent(i, data_elec_mll->GetBinContent(i));
    }
    TCanvas *Mee_plot = new TCanvas("Mee_plot", "Baseline M_{ee}", 500, 500);
    Mee_plot->SetLeftMargin(0.12);
    Mee_plot->SetLogy();
    data_elec_mll->SetMarkerColor(kBlack);
    data_elec_mll->SetMarkerStyle(kFullCircle);
    data_elec_mll->SetLineColor(kBlack);
    data_elec_mll->SetTitle("");
    data_elec_mll->SetYTitle("Events");
    data_elec_mll->GetYaxis()->SetTitleOffset(1.5);
    data_elec_mll->SetMaximum(data_elec_mll->GetMaximum() * 10);
    data_elec_mll->SetMinimum(0.1);
    data_elec_mll->SetStats(kFALSE);
    data_elec_mll->Draw();
    h_mll_window->SetFillColor(kYellow);
    h_mll_window->SetFillStyle(3003);
    h_mll_window->SetLineColor(kYellow);
    h_mll_window->Draw("hist,same");
    h_tail->SetFillColor(kYellow);
    h_tail->SetFillStyle(3003);
    h_tail->SetLineColor(kYellow);
    h_tail->Draw("hist,same");
}

void yt_background_template_mll_plot(int pt_bin_low = 0, int pt_bin_up = -1, int eta_bin_low = 0, int eta_bin_up = -1)
{
    // TFile *data_elec = TFile::Open(path + electron);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);

    TH1F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso;
    TH1F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight;
    TH1F *h_bkg_template_fail_CaloIso_and_TrackIso;

    if (pt_bin_low == 0 && pt_bin_up == -1 &&
        eta_bin_low == 0 && eta_bin_up == -1) {
        // 1 dim histograms
        h_bkg_template_fail_id_and_CaloIso_and_TrackIso = (TH1F *)data_elec->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso"); // baseline template
        h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight = (TH1F *)data_elec->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight"); // variation2 bkg template
        h_bkg_template_fail_CaloIso_and_TrackIso = (TH1F *)data_elec->Get("h_bkg_template_fail_CaloIso_and_TrackIso"); // variation1 bkg template
        
        cout << h_bkg_template_fail_id_and_CaloIso_and_TrackIso << endl;
        cout << h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight << endl;
        cout << h_bkg_template_fail_CaloIso_and_TrackIso << endl;
    }
    else {
        // 3-dim histograms
        TH3F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll = (TH3F *)data_elec->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll"); // baseline template
        TH3F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll = (TH3F *)data_elec->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll"); // variation2 bkg template
        TH3F *h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll = (TH3F *)data_elec->Get("h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll"); // variation1 bkg template

        // project on Z axis
        h_bkg_template_fail_id_and_CaloIso_and_TrackIso = (TH1F *)h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->ProjectionZ("baseline", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
        h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight = (TH1F *)h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll->ProjectionZ("template2", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
        h_bkg_template_fail_CaloIso_and_TrackIso = (TH1F *)h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll->ProjectionZ("template1", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
    }

    double pt_bin_low_value, eta_bin_low_value;
    double pt_bin_up_value, eta_bin_up_value;
    TH3F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll = (TH3F *)data_elec->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll");
    if (pt_bin_low != 0 && pt_bin_up != -1) {
        pt_bin_low_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetXaxis()->GetBinLowEdge(pt_bin_low);
        pt_bin_up_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetXaxis()->GetBinUpEdge(pt_bin_up);
    }
    else {
        pt_bin_low_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetXaxis()->GetBinLowEdge(1);
        int nbins = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetXaxis()->GetNbins();
        pt_bin_up_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetXaxis()->GetBinUpEdge(nbins - 1);
    }
    if (eta_bin_low_value != 0 && eta_bin_up_value != -1) {
        eta_bin_low_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetYaxis()->GetBinLowEdge(eta_bin_low);
        eta_bin_up_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetYaxis()->GetBinUpEdge(eta_bin_up);
    }
    else {
        eta_bin_low_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetYaxis()->GetBinLowEdge(1);
        int nbins = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetYaxis()->GetNbins();
        eta_bin_up_value = h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->GetYaxis()->GetBinUpEdge(nbins - 1);
    }
    cout << "(pt_bin_low_value, pt_bin_up_value, eta_bin_low_value, eta_bin_up_value) = ("
        << pt_bin_low_value << ", " << pt_bin_up_value << ", " << eta_bin_low_value << ", " << eta_bin_up_value << ")" << endl;

    // Convert double to string using stringstream
    stringstream sstream_pt_low, sstream_pt_up, sstream_eta_low, sstream_eta_up;
    sstream_pt_low << pt_bin_low_value;
    sstream_pt_up << pt_bin_up_value;
    sstream_eta_low << eta_bin_low_value * 100.;
    sstream_eta_up << eta_bin_up_value * 100;

    TCanvas *c_ee = new TCanvas("c_ee", "c_ee", 600, 600);
    c_ee->SetLeftMargin(0.12);

    // The current template 1 and 2 are exactly the same.
    h_bkg_template_fail_CaloIso_and_TrackIso->SetTitle("");
    h_bkg_template_fail_CaloIso_and_TrackIso->SetXTitle("M_{ee} [GeV]");
    h_bkg_template_fail_CaloIso_and_TrackIso->SetYTitle("Events");
    h_bkg_template_fail_CaloIso_and_TrackIso->GetYaxis()->SetTitleOffset(1.5);
    h_bkg_template_fail_CaloIso_and_TrackIso->SetLineColor(kBlue);
    h_bkg_template_fail_CaloIso_and_TrackIso->SetStats(kFALSE);
    h_bkg_template_fail_CaloIso_and_TrackIso->Draw("hist");
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso->SetTitle("");
    //h_bkg_template_fail_id_and_CaloIso_and_TrackIso->SetXTitle("M_{ee} [GeV]");
    //h_bkg_template_fail_id_and_CaloIso_and_TrackIso->SetYTitle("Events");
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso->SetStats(kFALSE);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso->SetLineColor(kBlack);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso->Draw("hist,same");
    //h_bkg_template_fail_id_and_CaloIso_and_TrackIso->Draw("hist");
    //h_bkg_template_fail_CaloIso_and_TrackIso->Draw("hist,same");
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight->SetLineColor(kRed);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight->Draw("hist,same");
    TLegend *leg = new TLegend(0.5, 0.6, 0.85, 0.8);
    leg->AddEntry(h_bkg_template_fail_id_and_CaloIso_and_TrackIso, "baseline template", "l");
    leg->AddEntry(h_bkg_template_fail_CaloIso_and_TrackIso, "variation 1 template", "l");
    leg->AddEntry(h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight, "variation 2 template", "l");
    leg->SetBorderSize(0);
    leg->Draw("same");

    string output_filename = "bkg_template_electron_pt" + sstream_pt_low.str() + sstream_pt_up.str() + "_eta" + sstream_eta_low.str() + sstream_eta_up.str() + ".pdf";
    c_ee->SaveAs(output_filename.c_str(), "pdf");
}

//
// distributions
//

void yt_kinematics_distribution()
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);

    // TFile *TandP_elec = TFile::Open(path + Zee_root);
    // TFile *TandP_muon = TFile::Open(path + Zmumu_root);
    TFile *TandP_elec = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *TandP_muon = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *truth_match_elec = TFile::Open(path + Zee_truth_root);
    TFile *truth_match_muon = TFile::Open(path + Zmumu_truth_root);

    TFile *TandP_truth_match_elec = TFile::Open(path + Zee_TP_truth_root);
    TFile *TandP_truth_match_muon = TFile::Open(path + Zmumu_TP_truth_root);
    
    TFile *ttbar_elec = TFile::Open(path + ttbar_elec_root);
    TFile *ttbar_muon = TFile::Open(path + ttbar_muon_root);

    TFile *Gtt_elec = TFile::Open(path + Gtt_elec_root);
    TFile *Gtt_muon = TFile::Open(path + Gtt_muon_root);

    TH1F *h_baseline_pt_data_elec = (TH1F *)data_elec->Get("h_baseline_pt");
    TH1F *h_baseline_eta_data_elec = (TH1F *)data_elec->Get("h_baseline_eta");

    TH1F *h_baseline_pt_data_muon = (TH1F *)data_muon->Get("h_baseline_pt");
    TH1F *h_baseline_eta_data_muon = (TH1F *)data_muon->Get("h_baseline_eta");

    TH1F *h_baseline_pt_TandP_elec = (TH1F *)TandP_elec->Get("h_baseline_pt");
    TH1F *h_baseline_eta_TandP_elec = (TH1F *)TandP_elec->Get("h_baseline_eta");

    TH1F *h_baseline_pt_TandP_muon = (TH1F *)TandP_muon->Get("h_baseline_pt");
    TH1F *h_baseline_eta_TandP_muon = (TH1F *)TandP_muon->Get("h_baseline_eta");

    TH1F *h_baseline_pt_truth_match_elec = (TH1F *)truth_match_elec->Get("h_baseline_pt");
    TH1F *h_baseline_eta_truth_match_elec = (TH1F *)truth_match_elec->Get("h_baseline_eta");

    TH1F *h_baseline_pt_truth_match_muon = (TH1F *)truth_match_muon->Get("h_baseline_pt");
    TH1F *h_baseline_eta_truth_match_muon = (TH1F *)truth_match_muon->Get("h_baseline_eta");

    TH1F *h_baseline_pt_TandP_truth_match_elec = (TH1F *)TandP_truth_match_elec->Get("h_baseline_pt");
    TH1F *h_baseline_eta_TandP_truth_match_elec = (TH1F *)TandP_truth_match_elec->Get("h_baseline_eta");

    TH1F *h_baseline_pt_TandP_truth_match_muon = (TH1F *)TandP_truth_match_muon->Get("h_baseline_pt");
    TH1F *h_baseline_eta_TandP_truth_match_muon = (TH1F *)TandP_truth_match_muon->Get("h_baseline_eta");

    TH1F *h_baseline_pt_ttbar_elec = (TH1F *)ttbar_elec->Get("h_baseline_pt");
    TH1F *h_baseline_eta_ttbar_elec = (TH1F *)ttbar_elec->Get("h_baseline_eta");

    TH1F *h_baseline_pt_ttbar_muon = (TH1F *)ttbar_muon->Get("h_baseline_pt");
    TH1F *h_baseline_eta_ttbar_muon = (TH1F *)ttbar_muon->Get("h_baseline_eta");

    TH1F *h_baseline_pt_Gtt_elec = (TH1F *)Gtt_elec->Get("h_baseline_pt");
    TH1F *h_baseline_eta_Gtt_elec = (TH1F *)Gtt_elec->Get("h_baseline_eta");

    TH1F *h_baseline_pt_Gtt_muon = (TH1F *)Gtt_muon->Get("h_baseline_pt");
    TH1F *h_baseline_eta_Gtt_muon = (TH1F *)Gtt_muon->Get("h_baseline_eta");

    h_baseline_pt_data_elec->Sumw2();
    h_baseline_eta_data_elec->Sumw2();
    h_baseline_pt_data_muon->Sumw2();
    h_baseline_eta_data_muon->Sumw2();
    h_baseline_pt_TandP_elec->Sumw2();
    h_baseline_eta_TandP_elec->Sumw2();
    h_baseline_pt_TandP_muon->Sumw2();
    h_baseline_eta_TandP_muon->Sumw2();
    h_baseline_pt_truth_match_elec->Sumw2();
    h_baseline_eta_truth_match_elec->Sumw2();
    h_baseline_pt_TandP_truth_match_muon->Sumw2();
    h_baseline_eta_TandP_truth_match_muon->Sumw2();
    h_baseline_pt_TandP_truth_match_elec->Sumw2();
    h_baseline_eta_TandP_truth_match_elec->Sumw2();
    h_baseline_pt_TandP_truth_match_muon->Sumw2();
    h_baseline_eta_TandP_truth_match_muon->Sumw2();
    h_baseline_pt_ttbar_elec->Sumw2();
    h_baseline_eta_ttbar_elec->Sumw2();
    h_baseline_pt_ttbar_muon->Sumw2();
    h_baseline_eta_ttbar_muon->Sumw2();
    h_baseline_pt_Gtt_elec->Sumw2();
    h_baseline_eta_Gtt_elec->Sumw2();
    h_baseline_pt_Gtt_muon->Sumw2();
    h_baseline_eta_Gtt_muon->Sumw2();

    h_baseline_pt_data_elec->Scale(1. / h_baseline_pt_data_elec->Integral());
    h_baseline_eta_data_elec->Scale(1. / h_baseline_eta_data_elec->Integral());
    h_baseline_pt_data_muon->Scale(1. / h_baseline_pt_data_muon->Integral());
    h_baseline_eta_data_muon->Scale(1. / h_baseline_eta_data_muon->Integral());
    h_baseline_pt_TandP_elec->Scale(1. / h_baseline_pt_TandP_elec->Integral());
    h_baseline_eta_TandP_elec->Scale(1. / h_baseline_eta_TandP_elec->Integral());
    h_baseline_pt_TandP_muon->Scale(1. / h_baseline_pt_TandP_muon->Integral());
    h_baseline_eta_TandP_muon->Scale(1. / h_baseline_eta_TandP_muon->Integral());
    h_baseline_pt_truth_match_elec->Scale(1. / h_baseline_pt_truth_match_elec->Integral());
    h_baseline_eta_truth_match_elec->Scale(1. / h_baseline_eta_truth_match_elec->Integral());
    h_baseline_pt_truth_match_muon->Scale(1. / h_baseline_pt_truth_match_muon->Integral());
    h_baseline_eta_truth_match_muon->Scale(1. / h_baseline_eta_truth_match_muon->Integral());
    h_baseline_pt_TandP_truth_match_elec->Scale(1. / h_baseline_pt_TandP_truth_match_elec->Integral());
    h_baseline_eta_TandP_truth_match_elec->Scale(1. / h_baseline_eta_TandP_truth_match_elec->Integral());
    h_baseline_pt_TandP_truth_match_muon->Scale(1. / h_baseline_pt_TandP_truth_match_muon->Integral());
    h_baseline_eta_TandP_truth_match_muon->Scale(1. / h_baseline_eta_TandP_truth_match_muon->Integral());
    h_baseline_pt_ttbar_elec->Scale(1. / h_baseline_pt_ttbar_elec->Integral());
    h_baseline_eta_ttbar_elec->Scale(1. / h_baseline_eta_ttbar_elec->Integral());
    h_baseline_pt_ttbar_muon->Scale(1. / h_baseline_pt_ttbar_muon->Integral());
    h_baseline_eta_ttbar_muon->Scale(1. / h_baseline_eta_ttbar_muon->Integral());
    h_baseline_pt_Gtt_elec->Scale(1. / h_baseline_pt_Gtt_elec->Integral());
    h_baseline_eta_Gtt_elec->Scale(1. / h_baseline_eta_Gtt_elec->Integral());
    h_baseline_pt_Gtt_muon->Scale(1. / h_baseline_pt_Gtt_muon->Integral());
    h_baseline_eta_Gtt_muon->Scale(1. / h_baseline_eta_Gtt_muon->Integral());

    TCanvas *kinematics = new TCanvas("kinematics", "kinematics", 1200, 1200);
    gStyle->SetOptStat(0);
    kinematics->Divide(2, 2);

    vector<double> vec_max;
    double max = 0.;

    kinematics->cd(1);
    gPad->SetLeftMargin(0.12);
    vec_max.push_back(h_baseline_pt_data_elec->GetMaximum());
    vec_max.push_back(h_baseline_pt_TandP_elec->GetMaximum());
    vec_max.push_back(h_baseline_pt_truth_match_elec->GetMaximum());
    vec_max.push_back(h_baseline_pt_ttbar_elec->GetMaximum());
    vec_max.push_back(h_baseline_pt_Gtt_elec->GetMaximum());
    sort(vec_max.begin(), vec_max.end());
    max = vec_max.back();
    h_baseline_pt_data_elec->SetTitle("");
    h_baseline_pt_data_elec->SetXTitle("p_{T} [GeV]");
    h_baseline_pt_data_elec->SetYTitle("Normalized number of probe electrons");
    h_baseline_pt_data_elec->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_pt_data_elec->SetMaximum(max * 1.1);
    h_baseline_pt_data_elec->SetMarkerStyle(kFullCircle);
    h_baseline_pt_data_elec->SetMarkerColor(kBlack);
    h_baseline_pt_data_elec->SetLineColor(kBlack);
    h_baseline_pt_data_elec->Draw("E1");
    h_baseline_pt_TandP_elec->SetMarkerStyle(kFullSquare);
    h_baseline_pt_TandP_elec->SetMarkerColor(kRed);
    h_baseline_pt_TandP_elec->SetLineColor(kRed);
    h_baseline_pt_TandP_elec->Draw("hist,same");
    h_baseline_pt_truth_match_elec->SetMarkerStyle(kFullTriangleUp);
    h_baseline_pt_truth_match_elec->SetMarkerColor(kBlue);
    h_baseline_pt_truth_match_elec->SetLineColor(kBlue);
    h_baseline_pt_truth_match_elec->Draw("hist,same");
    h_baseline_pt_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_baseline_pt_ttbar_elec->SetMarkerColor(kMagenta);
    h_baseline_pt_ttbar_elec->SetLineColor(kMagenta);
    h_baseline_pt_ttbar_elec->Draw("hist,same");
    h_baseline_pt_Gtt_elec->SetMarkerStyle(kFullCross);
    h_baseline_pt_Gtt_elec->SetMarkerColor(kOrange);
    h_baseline_pt_Gtt_elec->SetLineColor(kOrange);
    h_baseline_pt_Gtt_elec->Draw("hist,same");

    TLegend *leg1 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h_baseline_pt_data_elec, "Data");
    leg1->AddEntry(h_baseline_pt_TandP_elec, "Z #rightarrow ee T&P MC", "l");
    leg1->AddEntry(h_baseline_pt_truth_match_elec, "Z #rightarrow ee truth matched MC", "l");
    leg1->AddEntry(h_baseline_pt_ttbar_elec, "ttbar MC", "l");
    leg1->AddEntry(h_baseline_pt_Gtt_elec, "Gtt MC", "l");
    leg1->Draw();

    kinematics->cd(2);
    gPad->SetLeftMargin(0.12);
    vec_max.clear();
    vec_max.push_back(h_baseline_pt_data_muon->GetMaximum());
    vec_max.push_back(h_baseline_pt_TandP_muon->GetMaximum());
    vec_max.push_back(h_baseline_pt_truth_match_muon->GetMaximum());
    vec_max.push_back(h_baseline_pt_ttbar_muon->GetMaximum());
    vec_max.push_back(h_baseline_pt_Gtt_muon->GetMaximum());
    sort(vec_max.begin(), vec_max.end());
    max = vec_max.back();
    h_baseline_pt_data_muon->SetTitle("");
    h_baseline_pt_data_muon->SetXTitle("p_{T} [GeV]");
    h_baseline_pt_data_muon->SetYTitle("Normalized number of probe muons");
    h_baseline_pt_data_muon->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_pt_data_muon->SetMaximum(max * 1.1);
    h_baseline_pt_data_muon->SetMarkerStyle(kFullCircle);
    h_baseline_pt_data_muon->SetMarkerColor(kBlack);
    h_baseline_pt_data_muon->SetLineColor(kBlack);
    h_baseline_pt_data_muon->Draw("E1");
    h_baseline_pt_TandP_muon->SetMarkerStyle(kFullSquare);
    h_baseline_pt_TandP_muon->SetMarkerColor(kRed);
    h_baseline_pt_TandP_muon->SetLineColor(kRed);
    h_baseline_pt_TandP_muon->Draw("hist,same");
    h_baseline_pt_truth_match_muon->SetMarkerStyle(kFullTriangleUp);
    h_baseline_pt_truth_match_muon->SetMarkerColor(kBlue);
    h_baseline_pt_truth_match_muon->SetLineColor(kBlue);
    h_baseline_pt_truth_match_muon->Draw("hist,same");
    h_baseline_pt_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_baseline_pt_ttbar_muon->SetMarkerColor(kMagenta);
    h_baseline_pt_ttbar_muon->SetLineColor(kMagenta);
    h_baseline_pt_ttbar_muon->Draw("hist,same");
    h_baseline_pt_Gtt_muon->SetMarkerStyle(kFullCross);
    h_baseline_pt_Gtt_muon->SetMarkerColor(kOrange);
    h_baseline_pt_Gtt_muon->SetLineColor(kOrange);
    h_baseline_pt_Gtt_muon->Draw("hist,same");

    TLegend *leg2 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(h_baseline_pt_data_muon, "Data");
    leg2->AddEntry(h_baseline_pt_TandP_muon, "Z #rightarrow #mu#mu T&P MC", "l");
    leg2->AddEntry(h_baseline_pt_truth_match_muon, "Z #rightarrow #mu#mu truth matched MC", "l");
    leg2->AddEntry(h_baseline_pt_ttbar_muon, "ttbar MC", "l");
    leg2->AddEntry(h_baseline_pt_Gtt_muon, "Gtt MC", "l");
    leg2->Draw("same");

    kinematics->cd(3);
    gPad->SetLeftMargin(0.12);
    vec_max.clear();
    vec_max.push_back(h_baseline_eta_data_elec->GetMaximum());
    vec_max.push_back(h_baseline_eta_TandP_elec->GetMaximum());
    vec_max.push_back(h_baseline_eta_truth_match_elec->GetMaximum());
    vec_max.push_back(h_baseline_eta_ttbar_elec->GetMaximum());
    vec_max.push_back(h_baseline_eta_Gtt_elec->GetMaximum());
    sort(vec_max.begin(), vec_max.end());
    max = vec_max.back();
    h_baseline_eta_data_elec->SetTitle("");
    h_baseline_eta_data_elec->SetXTitle("|#eta|");
    h_baseline_eta_data_elec->SetYTitle("Normalized number of probe electrons");
    h_baseline_eta_data_elec->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_eta_data_elec->SetMaximum(max * 1.1);
    h_baseline_eta_data_elec->SetMarkerStyle(kFullCircle);
    h_baseline_eta_data_elec->SetMarkerColor(kBlack);
    h_baseline_eta_data_elec->SetLineColor(kBlack);
    h_baseline_eta_data_elec->Draw("E1");
    h_baseline_eta_TandP_elec->SetMarkerStyle(kFullSquare);
    h_baseline_eta_TandP_elec->SetMarkerColor(kRed);
    h_baseline_eta_TandP_elec->SetLineColor(kRed);
    h_baseline_eta_TandP_elec->Draw("hist,same");
    h_baseline_eta_truth_match_elec->SetMarkerStyle(kFullTriangleUp);
    h_baseline_eta_truth_match_elec->SetMarkerColor(kBlue);
    h_baseline_eta_truth_match_elec->SetLineColor(kBlue);
    h_baseline_eta_truth_match_elec->Draw("hist,same");
    h_baseline_eta_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_baseline_eta_ttbar_elec->SetMarkerColor(kMagenta);
    h_baseline_eta_ttbar_elec->SetLineColor(kMagenta);
    h_baseline_eta_ttbar_elec->Draw("hist,same");
    h_baseline_eta_Gtt_elec->SetMarkerStyle(kFullCross);
    h_baseline_eta_Gtt_elec->SetMarkerColor(kOrange);
    h_baseline_eta_Gtt_elec->SetLineColor(kOrange);
    h_baseline_eta_Gtt_elec->Draw("hist,same");

    TLegend *leg3 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.03);
    leg3->AddEntry(h_baseline_eta_data_elec, "Data");
    leg3->AddEntry(h_baseline_eta_TandP_elec, "Z #rightarrow ee T&P MC", "l");
    leg3->AddEntry(h_baseline_eta_truth_match_elec, "Z #rightarrow ee truth matched MC", "l");
    leg3->AddEntry(h_baseline_eta_ttbar_elec, "ttbar MC", "l");
    leg3->AddEntry(h_baseline_eta_Gtt_elec, "Gtt MC", "l");
    leg3->Draw("same");

    kinematics->cd(4);
    gPad->SetLeftMargin(0.12);
    vec_max.clear();
    vec_max.push_back(h_baseline_eta_data_muon->GetMaximum());
    vec_max.push_back(h_baseline_eta_TandP_muon->GetMaximum());
    vec_max.push_back(h_baseline_eta_truth_match_muon->GetMaximum());
    vec_max.push_back(h_baseline_eta_ttbar_muon->GetMaximum());
    vec_max.push_back(h_baseline_eta_Gtt_muon->GetMaximum());
    sort(vec_max.begin(), vec_max.end());
    max = vec_max.back();
    h_baseline_eta_data_muon->SetTitle("");
    h_baseline_eta_data_muon->SetXTitle("|#eta|");
    h_baseline_eta_data_muon->SetYTitle("Normalized number of probe muons");
    h_baseline_eta_data_muon->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_eta_data_muon->SetMaximum(max * 1.1);
    h_baseline_eta_data_muon->SetMarkerStyle(kFullCircle);
    h_baseline_eta_data_muon->SetMarkerColor(kBlack);
    h_baseline_eta_data_muon->SetLineColor(kBlack);
    h_baseline_eta_data_muon->Draw("E1");
    h_baseline_eta_TandP_muon->SetMarkerStyle(kFullCircle);
    h_baseline_eta_TandP_muon->SetMarkerColor(kRed);
    h_baseline_eta_TandP_muon->SetLineColor(kRed);
    h_baseline_eta_TandP_muon->Draw("hist,same");
    h_baseline_eta_truth_match_muon->SetMarkerStyle(kFullTriangleUp);
    h_baseline_eta_truth_match_muon->SetMarkerColor(kBlue);
    h_baseline_eta_truth_match_muon->SetLineColor(kBlue);
    h_baseline_eta_truth_match_muon->Draw("hist,same");
    h_baseline_eta_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_baseline_eta_ttbar_muon->SetMarkerColor(kMagenta);
    h_baseline_eta_ttbar_muon->SetLineColor(kMagenta);
    h_baseline_eta_ttbar_muon->Draw("hist,same");
    h_baseline_eta_Gtt_muon->SetMarkerStyle(kFullCross);
    h_baseline_eta_Gtt_muon->SetMarkerColor(kOrange);
    h_baseline_eta_Gtt_muon->SetLineColor(kOrange);
    h_baseline_eta_Gtt_muon->Draw("hist,same");

    TLegend *leg4 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg4->SetBorderSize(0);
    leg4->SetFillColor(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.03);
    leg4->AddEntry(h_baseline_eta_data_muon, "Data");
    leg4->AddEntry(h_baseline_eta_TandP_muon, "Z #rightarrow #mu#mu T&P MC", "l");
    leg4->AddEntry(h_baseline_eta_truth_match_muon, "Z #rightarrow #mu#mu truth matched MC", "l");
    leg4->AddEntry(h_baseline_eta_ttbar_muon, "ttbar MC", "l");
    leg4->AddEntry(h_baseline_eta_Gtt_muon, "Gtt MC", "l");
    leg4->Draw("same");

    kinematics->SaveAs("baseline_kinematics.pdf", "pdf");
}

void yt_deltaR_and_NJets_distributions()
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);

    // TFile *TandP_elec = TFile::Open(path + Zee_root);
    // TFile *TandP_muon = TFile::Open(path + Zmumu_root);
    TFile *TandP_elec = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *TandP_muon = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *truth_match_elec = TFile::Open(path + Zee_truth_root);
    TFile *truth_match_muon = TFile::Open(path + Zmumu_truth_root);

    TFile *TandP_truth_match_elec = TFile::Open(path + Zee_TP_truth_root);
    TFile *TandP_truth_match_muon = TFile::Open(path + Zmumu_TP_truth_root);
    
    TFile *ttbar_elec = TFile::Open(path + ttbar_elec_root);
    TFile *ttbar_muon = TFile::Open(path + ttbar_muon_root);

    TFile *Gtt_elec = TFile::Open(path + Gtt_elec_root);
    TFile *Gtt_muon = TFile::Open(path + Gtt_muon_root);

    // baseline
    TH1F *h_baseline_dRjet_data_elec = (TH1F *)data_elec->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_data_elec = (TH1F *)data_elec->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_data_muon = (TH1F *)data_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_data_muon = (TH1F *)data_muon->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_TandP_elec = (TH1F *)TandP_elec->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_TandP_elec = (TH1F *)TandP_elec->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_TandP_muon = (TH1F *)TandP_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_TandP_muon = (TH1F *)TandP_muon->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_truth_match_elec = (TH1F *)truth_match_elec->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_truth_match_elec = (TH1F *)truth_match_elec->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_truth_match_muon = (TH1F *)truth_match_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_truth_match_muon = (TH1F *)truth_match_muon->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_TandP_truth_match_elec = (TH1F *)TandP_truth_match_elec->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_TandP_truth_match_elec = (TH1F *)TandP_truth_match_elec->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_TandP_truth_match_muon = (TH1F *)TandP_truth_match_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_TandP_truth_match_muon = (TH1F *)TandP_truth_match_muon->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_ttbar_elec = (TH1F *)ttbar_elec->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_ttbar_elec = (TH1F *)ttbar_elec->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_ttbar_muon = (TH1F *)ttbar_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_ttbar_muon = (TH1F *)ttbar_muon->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_Gtt_elec = (TH1F *)Gtt_elec->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_Gtt_elec = (TH1F *)Gtt_elec->Get("h_baseline_nJets");

    TH1F *h_baseline_dRjet_Gtt_muon = (TH1F *)Gtt_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_nJets_Gtt_muon = (TH1F *)Gtt_muon->Get("h_baseline_nJets");

    h_baseline_dRjet_data_elec->Sumw2();
    h_baseline_nJets_data_elec->Sumw2();
    h_baseline_dRjet_data_muon->Sumw2();
    h_baseline_nJets_data_muon->Sumw2();
    h_baseline_dRjet_TandP_elec->Sumw2();
    h_baseline_nJets_TandP_elec->Sumw2();
    h_baseline_dRjet_TandP_muon->Sumw2();
    h_baseline_nJets_TandP_muon->Sumw2();
    h_baseline_dRjet_truth_match_elec->Sumw2();
    h_baseline_nJets_truth_match_elec->Sumw2();
    h_baseline_dRjet_truth_match_muon->Sumw2();
    h_baseline_nJets_truth_match_muon->Sumw2();
    h_baseline_dRjet_TandP_truth_match_elec->Sumw2();
    h_baseline_nJets_TandP_truth_match_elec->Sumw2();
    h_baseline_dRjet_TandP_truth_match_muon->Sumw2();
    h_baseline_nJets_TandP_truth_match_muon->Sumw2();
    h_baseline_dRjet_ttbar_elec->Sumw2();
    h_baseline_nJets_ttbar_elec->Sumw2();
    h_baseline_dRjet_ttbar_muon->Sumw2();
    h_baseline_nJets_ttbar_muon->Sumw2();
    h_baseline_dRjet_Gtt_elec->Sumw2();
    h_baseline_nJets_Gtt_elec->Sumw2();
    h_baseline_dRjet_Gtt_muon->Sumw2();
    h_baseline_nJets_Gtt_muon->Sumw2();

    h_baseline_dRjet_data_elec->Scale(1. / h_baseline_dRjet_data_elec->Integral());
    h_baseline_nJets_data_elec->Scale(1. / h_baseline_nJets_data_elec->Integral());
    h_baseline_dRjet_data_muon->Scale(1. / h_baseline_dRjet_data_muon->Integral());
    h_baseline_nJets_data_muon->Scale(1. / h_baseline_nJets_data_muon->Integral());
    h_baseline_dRjet_TandP_elec->Scale(1. / h_baseline_dRjet_TandP_elec->Integral());
    h_baseline_nJets_TandP_elec->Scale(1. / h_baseline_nJets_TandP_elec->Integral());
    h_baseline_dRjet_TandP_muon->Scale(1. / h_baseline_dRjet_TandP_muon->Integral());
    h_baseline_nJets_TandP_muon->Scale(1. / h_baseline_nJets_TandP_muon->Integral());
    h_baseline_dRjet_truth_match_elec->Scale(1. / h_baseline_dRjet_truth_match_elec->Integral());
    h_baseline_nJets_truth_match_elec->Scale(1. / h_baseline_nJets_truth_match_elec->Integral());
    h_baseline_dRjet_truth_match_muon->Scale(1. / h_baseline_dRjet_truth_match_muon->Integral());
    h_baseline_nJets_truth_match_muon->Scale(1. / h_baseline_nJets_truth_match_muon->Integral());
    h_baseline_dRjet_TandP_truth_match_elec->Scale(1. / h_baseline_dRjet_TandP_truth_match_elec->Integral());
    h_baseline_nJets_TandP_truth_match_elec->Scale(1. / h_baseline_nJets_TandP_truth_match_elec->Integral());
    h_baseline_dRjet_TandP_truth_match_muon->Scale(1. / h_baseline_dRjet_TandP_truth_match_muon->Integral());
    h_baseline_nJets_TandP_truth_match_muon->Scale(1. / h_baseline_nJets_TandP_truth_match_muon->Integral());
    h_baseline_dRjet_ttbar_elec->Scale(1. / h_baseline_dRjet_ttbar_elec->Integral());
    h_baseline_nJets_ttbar_elec->Scale(1. / h_baseline_nJets_ttbar_elec->Integral());
    h_baseline_dRjet_ttbar_muon->Scale(1. / h_baseline_dRjet_ttbar_muon->Integral());
    h_baseline_nJets_ttbar_muon->Scale(1. / h_baseline_nJets_ttbar_muon->Integral());
    h_baseline_dRjet_Gtt_elec->Scale(1. / h_baseline_dRjet_Gtt_elec->Integral());
    h_baseline_nJets_Gtt_elec->Scale(1. / h_baseline_nJets_Gtt_elec->Integral());
    h_baseline_dRjet_Gtt_muon->Scale(1. / h_baseline_dRjet_Gtt_muon->Integral());
    h_baseline_nJets_Gtt_muon->Scale(1. / h_baseline_nJets_Gtt_muon->Integral());

    TCanvas *distributions = new TCanvas("distributions", "distributions", 1200, 1200);
    //TCanvas *distributions = new TCanvas("distributions", "distributions", 1200, 600); // zoom in dR region
    gStyle->SetOptStat(0);
    distributions->Divide(2, 2);
    //distributions->Divide(2, 1); // zoom in dR region

    distributions->cd(1);
    gPad->SetLeftMargin(0.12);
    h_baseline_dRjet_data_elec->SetTitle("");
    h_baseline_dRjet_data_elec->SetXTitle("#Delta R(e, jet)");
    h_baseline_dRjet_data_elec->SetYTitle("Normalized number of probe electrons");
    h_baseline_dRjet_data_elec->GetYaxis()->SetTitleOffset(1.5);
    //h_baseline_dRjet_data_elec->SetMaximum(h_baseline_dRjet_Gtt_elec->GetMaximum() * 1.1);
    h_baseline_dRjet_data_elec->SetMaximum(h_baseline_dRjet_ttbar_elec->GetMaximum() * 1.1);
    h_baseline_dRjet_data_elec->GetXaxis()->SetRangeUser(0, 4);
    h_baseline_dRjet_data_elec->SetMarkerStyle(kFullCircle);
    h_baseline_dRjet_data_elec->SetMarkerColor(kBlack);
    h_baseline_dRjet_data_elec->SetLineColor(kBlack);
    h_baseline_dRjet_data_elec->Draw("E1");
    h_baseline_dRjet_TandP_elec->SetMarkerStyle(kFullSquare);
    h_baseline_dRjet_TandP_elec->SetMarkerColor(kRed);
    h_baseline_dRjet_TandP_elec->SetLineColor(kRed);
    h_baseline_dRjet_TandP_elec->Draw("hist,same");
    h_baseline_dRjet_truth_match_elec->SetMarkerStyle(kFullTriangleUp);
    h_baseline_dRjet_truth_match_elec->SetMarkerColor(kBlue);
    h_baseline_dRjet_truth_match_elec->SetLineColor(kBlue);
    h_baseline_dRjet_truth_match_elec->Draw("hist,same");
    h_baseline_dRjet_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_baseline_dRjet_ttbar_elec->SetMarkerColor(kMagenta);
    h_baseline_dRjet_ttbar_elec->SetLineColor(kMagenta);
    h_baseline_dRjet_ttbar_elec->Draw("hist,same");
    h_baseline_dRjet_Gtt_elec->SetMarkerStyle(kFullCross);
    h_baseline_dRjet_Gtt_elec->SetMarkerColor(kOrange);
    h_baseline_dRjet_Gtt_elec->SetLineColor(kOrange);
    h_baseline_dRjet_Gtt_elec->Draw("hist,same");

    TLegend *leg1 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h_baseline_dRjet_data_elec, "Data");
    leg1->AddEntry(h_baseline_dRjet_TandP_elec, "Z #rightarrow ee T&P MC", "l");
    leg1->AddEntry(h_baseline_dRjet_truth_match_elec, "Z #rightarrow ee truth matched MC", "l");
    leg1->AddEntry(h_baseline_dRjet_ttbar_elec, "ttbar MC", "l");
    leg1->AddEntry(h_baseline_dRjet_Gtt_elec, "Gtt MC", "l");
    leg1->Draw();

    distributions->cd(2);
    gPad->SetLeftMargin(0.12);
    h_baseline_dRjet_data_muon->SetTitle("");
    h_baseline_dRjet_data_muon->SetXTitle("#Delta R(#mu, jet)");
    h_baseline_dRjet_data_muon->SetYTitle("Normalized number of probe muons");
    h_baseline_dRjet_data_muon->GetYaxis()->SetTitleOffset(1.5);
    //h_baseline_dRjet_data_muon->SetMaximum(h_baseline_dRjet_Gtt_muon->GetMaximum() * 1.1);
    h_baseline_dRjet_data_muon->SetMaximum(h_baseline_dRjet_ttbar_muon->GetMaximum() * 1.1);
    h_baseline_dRjet_data_muon->GetXaxis()->SetRangeUser(0, 4);
    h_baseline_dRjet_data_muon->SetMarkerStyle(kFullCircle);
    h_baseline_dRjet_data_muon->SetMarkerColor(kBlack);
    h_baseline_dRjet_data_muon->SetLineColor(kBlack);
    h_baseline_dRjet_data_muon->Draw("E1");
    h_baseline_dRjet_TandP_muon->SetMarkerStyle(kFullSquare);
    h_baseline_dRjet_TandP_muon->SetMarkerColor(kRed);
    h_baseline_dRjet_TandP_muon->SetLineColor(kRed);
    h_baseline_dRjet_TandP_muon->Draw("hist,same");
    h_baseline_dRjet_truth_match_muon->SetMarkerStyle(kFullTriangleUp);
    h_baseline_dRjet_truth_match_muon->SetMarkerColor(kBlue);
    h_baseline_dRjet_truth_match_muon->SetLineColor(kBlue);
    h_baseline_dRjet_truth_match_muon->Draw("hist,same");
    h_baseline_dRjet_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_baseline_dRjet_ttbar_muon->SetMarkerColor(kMagenta);
    h_baseline_dRjet_ttbar_muon->SetLineColor(kMagenta);
    h_baseline_dRjet_ttbar_muon->Draw("hist,same");
    h_baseline_dRjet_Gtt_muon->SetMarkerStyle(kFullCross);
    h_baseline_dRjet_Gtt_muon->SetMarkerColor(kOrange);
    h_baseline_dRjet_Gtt_muon->SetLineColor(kOrange);
    h_baseline_dRjet_Gtt_muon->Draw("hist,same");

    TLegend *leg2 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(h_baseline_dRjet_data_muon, "Data");
    leg2->AddEntry(h_baseline_dRjet_TandP_muon, "Z #rightarrow #mu#mu T&P MC", "l");
    leg2->AddEntry(h_baseline_dRjet_truth_match_muon, "Z #rightarrow #mu#mu truth matched MC", "l");
    leg2->AddEntry(h_baseline_dRjet_ttbar_muon, "ttbar MC", "l");
    leg2->AddEntry(h_baseline_dRjet_Gtt_muon, "Gtt MC", "l");
    leg2->Draw("same");

    distributions->cd(3);
    gPad->SetLeftMargin(0.12);
    h_baseline_nJets_data_elec->SetTitle("");
    h_baseline_nJets_data_elec->SetXTitle("N_{jets}");
    h_baseline_nJets_data_elec->SetYTitle("Normalized number of probe electrons");
    h_baseline_nJets_data_elec->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_nJets_data_elec->SetMaximum(h_baseline_nJets_Gtt_elec->GetMaximum() * 1.1);
    h_baseline_nJets_data_elec->SetMarkerStyle(kFullCircle);
    h_baseline_nJets_data_elec->SetMarkerColor(kBlack);
    h_baseline_nJets_data_elec->SetLineColor(kBlack);
    h_baseline_nJets_data_elec->Draw("E1");
    h_baseline_nJets_TandP_elec->SetMarkerStyle(kFullSquare);
    h_baseline_nJets_TandP_elec->SetMarkerColor(kRed);
    h_baseline_nJets_TandP_elec->SetLineColor(kRed);
    h_baseline_nJets_TandP_elec->Draw("hist,same");
    h_baseline_nJets_truth_match_elec->SetMarkerStyle(kFullTriangleUp);
    h_baseline_nJets_truth_match_elec->SetMarkerColor(kBlue);
    h_baseline_nJets_truth_match_elec->SetLineColor(kBlue);
    h_baseline_nJets_truth_match_elec->Draw("hist,same");
    h_baseline_nJets_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_baseline_nJets_ttbar_elec->SetMarkerColor(kMagenta);
    h_baseline_nJets_ttbar_elec->SetLineColor(kMagenta);
    h_baseline_nJets_ttbar_elec->Draw("hist,same");
    h_baseline_nJets_Gtt_elec->SetMarkerStyle(kFullCross);
    h_baseline_nJets_Gtt_elec->SetMarkerColor(kOrange);
    h_baseline_nJets_Gtt_elec->SetLineColor(kOrange);
    h_baseline_nJets_Gtt_elec->Draw("hist,same");

    TLegend *leg3 = new TLegend(0.2, 0.6, 0.5, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.03);
    leg3->AddEntry(h_baseline_nJets_data_elec, "Data");
    leg3->AddEntry(h_baseline_nJets_TandP_elec, "Z #rightarrow ee T&P MC", "l");
    leg3->AddEntry(h_baseline_nJets_truth_match_elec, "Z #rightarrow ee truth matched MC", "l");
    leg3->AddEntry(h_baseline_nJets_ttbar_elec, "ttbar MC", "l");
    leg3->AddEntry(h_baseline_nJets_Gtt_elec, "Gtt MC", "l");
    leg3->Draw("same");

    distributions->cd(4);
    gPad->SetLeftMargin(0.12);
    h_baseline_nJets_data_muon->SetTitle("");
    h_baseline_nJets_data_muon->SetXTitle("N_{jets}");
    h_baseline_nJets_data_muon->SetYTitle("Normalized number of probe muons");
    h_baseline_nJets_data_muon->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_nJets_data_muon->SetMaximum(h_baseline_nJets_Gtt_muon->GetMaximum() * 1.1);
    h_baseline_nJets_data_muon->SetMarkerStyle(kFullCircle);
    h_baseline_nJets_data_muon->SetMarkerColor(kBlack);
    h_baseline_nJets_data_muon->SetLineColor(kBlack);
    h_baseline_nJets_data_muon->Draw("E1");
    h_baseline_nJets_TandP_muon->SetMarkerStyle(kFullCircle);
    h_baseline_nJets_TandP_muon->SetMarkerColor(kRed);
    h_baseline_nJets_TandP_muon->SetLineColor(kRed);
    h_baseline_nJets_TandP_muon->Draw("hist,same");
    h_baseline_nJets_truth_match_muon->SetMarkerStyle(kFullTriangleUp);
    h_baseline_nJets_truth_match_muon->SetMarkerColor(kBlue);
    h_baseline_nJets_truth_match_muon->SetLineColor(kBlue);
    h_baseline_nJets_truth_match_muon->Draw("hist,same");
    h_baseline_nJets_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_baseline_nJets_ttbar_muon->SetMarkerColor(kMagenta);
    h_baseline_nJets_ttbar_muon->SetLineColor(kMagenta);
    h_baseline_nJets_ttbar_muon->Draw("hist,same");
    h_baseline_nJets_Gtt_muon->SetMarkerStyle(kFullCross);
    h_baseline_nJets_Gtt_muon->SetMarkerColor(kOrange);
    h_baseline_nJets_Gtt_muon->SetLineColor(kOrange);
    h_baseline_nJets_Gtt_muon->Draw("hist,same");

    TLegend *leg4 = new TLegend(0.2, 0.6, 0.5, 0.9);
    leg4->SetBorderSize(0);
    leg4->SetFillColor(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.03);
    leg4->AddEntry(h_baseline_nJets_data_muon, "Data");
    leg4->AddEntry(h_baseline_nJets_TandP_muon, "Z #rightarrow #mu#mu T&P MC", "l");
    leg4->AddEntry(h_baseline_nJets_truth_match_muon, "Z #rightarrow #mu#mu truth matched MC", "l");
    leg4->AddEntry(h_baseline_nJets_ttbar_muon, "ttbar MC", "l");
    leg4->AddEntry(h_baseline_nJets_Gtt_muon, "Gtt MC", "l");
    leg4->Draw("same");

    distributions->SaveAs("baseline_deltaR_and_NJets.pdf", "pdf");
}

void yt_deltaR_baseline_and_signal_distributions()
{
    // TFile *data_elec = TFile::Open(path + electron);
    // TFile *data_muon = TFile::Open(path + muon);
    TFile *data_elec = TFile::Open(path + electron_tag_trigger_matched);
    TFile *data_muon = TFile::Open(path + muon_tag_trigger_matched);

    // TFile *TandP_elec = TFile::Open(path + Zee_root);
    // TFile *TandP_muon = TFile::Open(path + Zmumu_root);
    TFile *TandP_elec = TFile::Open(path + Zee_tag_trigger_matched_root);
    TFile *TandP_muon = TFile::Open(path + Zmumu_tag_trigger_matched_root);

    TFile *truth_match_elec = TFile::Open(path + Zee_truth_root);
    TFile *truth_match_muon = TFile::Open(path + Zmumu_truth_root);

    TFile *TandP_truth_match_elec = TFile::Open(path + Zee_TP_truth_root);
    TFile *TandP_truth_match_muon = TFile::Open(path + Zmumu_TP_truth_root);

    TFile *ttbar_elec = TFile::Open(path + ttbar_elec_root);
    TFile *ttbar_muon = TFile::Open(path + ttbar_muon_root);

    TFile *Gtt_elec = TFile::Open(path + Gtt_elec_root);
    TFile *Gtt_muon = TFile::Open(path + Gtt_muon_root);

    // For electron case, we only use |eta| < 2.0 to calculate the efficiency
    int eta_low_bin = ((TH1F *)data_elec->Get("h_baseline_pt_eta_dRjet"))->GetYaxis()->FindBin(0+0.01);
    int eta_up_bin = ((TH1F *)data_elec->Get("h_baseline_pt_eta_dRjet"))->GetYaxis()->FindBin(2-0.01);

    // baseline
    // We use baseline electron with |eta| < 2 to calculate efficiency
    TH1F *h_baseline_dRjet_data_elec = (TH1F *)((TH3F *)data_elec->Get("h_baseline_pt_eta_dRjet"))->ProjectionZ("h_baseline_dRjet_data_elec", 0, -1, eta_low_bin, eta_up_bin);
    TH1F *h_baseline_dRjet_TandP_elec = (TH1F *)((TH3F *)TandP_elec->Get("h_baseline_pt_eta_dRjet"))->ProjectionZ("h_baseline_dRjet_TandP_elec", 0, -1, eta_low_bin, eta_up_bin);
    TH1F *h_baseline_dRjet_truth_match_elec = (TH1F *)((TH3F *)truth_match_elec->Get("h_baseline_pt_eta_dRjet"))->ProjectionZ("h_baseline_dRjet_truth_match_elec", 0, -1, eta_low_bin, eta_up_bin);
    TH1F *h_baseline_dRjet_TandP_truth_match_elec = (TH1F *)((TH3F *)TandP_truth_match_elec->Get("h_baseline_pt_eta_dRjet"))->ProjectionZ("h_baseline_dRjet_TandP_truth_match_elec", 0, -1, eta_low_bin, eta_up_bin);
    TH1F *h_baseline_dRjet_ttbar_elec = (TH1F *)((TH3F *)ttbar_elec->Get("h_baseline_pt_eta_dRjet"))->ProjectionZ("h_baseline_dRjet_ttbar_elec", 0, -1, eta_low_bin, eta_up_bin);
    TH1F *h_baseline_dRjet_Gtt_elec = (TH1F *)((TH3F *)Gtt_elec->Get("h_baseline_pt_eta_dRjet"))->ProjectionZ("h_baseline_dRjet_Gtt_elec", 0, -1, eta_low_bin, eta_up_bin);
    // No |eta| requirement on muon cases
    TH1F *h_baseline_dRjet_data_muon = (TH1F *)data_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_dRjet_TandP_muon = (TH1F *)TandP_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_dRjet_truth_match_muon = (TH1F *)truth_match_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_dRjet_TandP_truth_match_muon = (TH1F *)TandP_truth_match_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_dRjet_ttbar_muon = (TH1F *)ttbar_muon->Get("h_baseline_dRjet");
    TH1F *h_baseline_dRjet_Gtt_muon = (TH1F *)Gtt_muon->Get("h_baseline_dRjet");

    h_baseline_dRjet_data_elec->Sumw2();
    h_baseline_dRjet_data_muon->Sumw2();
    h_baseline_dRjet_TandP_elec->Sumw2();
    h_baseline_dRjet_TandP_muon->Sumw2();
    h_baseline_dRjet_truth_match_elec->Sumw2();
    h_baseline_dRjet_truth_match_muon->Sumw2();
    h_baseline_dRjet_TandP_truth_match_elec->Sumw2();
    h_baseline_dRjet_TandP_truth_match_muon->Sumw2();
    h_baseline_dRjet_ttbar_elec->Sumw2();
    h_baseline_dRjet_ttbar_muon->Sumw2();
    h_baseline_dRjet_Gtt_elec->Sumw2();
    h_baseline_dRjet_Gtt_muon->Sumw2();

    h_baseline_dRjet_data_elec->Scale(1. / h_baseline_dRjet_data_elec->Integral());
    h_baseline_dRjet_data_muon->Scale(1. / h_baseline_dRjet_data_muon->Integral());
    h_baseline_dRjet_TandP_elec->Scale(1. / h_baseline_dRjet_TandP_elec->Integral());
    h_baseline_dRjet_TandP_muon->Scale(1. / h_baseline_dRjet_TandP_muon->Integral());
    h_baseline_dRjet_truth_match_elec->Scale(1. / h_baseline_dRjet_truth_match_elec->Integral());
    h_baseline_dRjet_truth_match_muon->Scale(1. / h_baseline_dRjet_truth_match_muon->Integral());
    h_baseline_dRjet_TandP_truth_match_elec->Scale(1. / h_baseline_dRjet_TandP_truth_match_elec->Integral());
    h_baseline_dRjet_TandP_truth_match_muon->Scale(1. / h_baseline_dRjet_TandP_truth_match_muon->Integral());
    h_baseline_dRjet_ttbar_elec->Scale(1. / h_baseline_dRjet_ttbar_elec->Integral());
    h_baseline_dRjet_ttbar_muon->Scale(1. / h_baseline_dRjet_ttbar_muon->Integral());
    h_baseline_dRjet_Gtt_elec->Scale(1. / h_baseline_dRjet_Gtt_elec->Integral());
    h_baseline_dRjet_Gtt_muon->Scale(1. / h_baseline_dRjet_Gtt_muon->Integral());

    // signal
    TH1F *h_signal_dRjet_data_elec = (TH1F *)data_elec->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_data_muon = (TH1F *)data_muon->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_TandP_elec = (TH1F *)TandP_elec->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_TandP_muon = (TH1F *)TandP_muon->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_truth_match_elec = (TH1F *)truth_match_elec->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_truth_match_muon = (TH1F *)truth_match_muon->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_TandP_truth_match_elec = (TH1F *)TandP_truth_match_elec->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_TandP_truth_match_muon = (TH1F *)TandP_truth_match_muon->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_ttbar_elec = (TH1F *)ttbar_elec->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_ttbar_muon = (TH1F *)ttbar_muon->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_Gtt_elec = (TH1F *)Gtt_elec->Get("h_signal_dRjet");
    TH1F *h_signal_dRjet_Gtt_muon = (TH1F *)Gtt_muon->Get("h_signal_dRjet");

    h_signal_dRjet_data_elec->Sumw2();
    h_signal_dRjet_data_muon->Sumw2();
    h_signal_dRjet_TandP_elec->Sumw2();
    h_signal_dRjet_TandP_muon->Sumw2();
    h_signal_dRjet_truth_match_elec->Sumw2();
    h_signal_dRjet_truth_match_muon->Sumw2();
    h_signal_dRjet_TandP_truth_match_elec->Sumw2();
    h_signal_dRjet_TandP_truth_match_muon->Sumw2();
    h_signal_dRjet_ttbar_elec->Sumw2();
    h_signal_dRjet_ttbar_muon->Sumw2();
    h_signal_dRjet_Gtt_elec->Sumw2();
    h_signal_dRjet_Gtt_muon->Sumw2();

    h_signal_dRjet_data_elec->Scale(1. / h_signal_dRjet_data_elec->Integral());
    h_signal_dRjet_data_muon->Scale(1. / h_signal_dRjet_data_muon->Integral());
    h_signal_dRjet_TandP_elec->Scale(1. / h_signal_dRjet_TandP_elec->Integral());
    h_signal_dRjet_TandP_muon->Scale(1. / h_signal_dRjet_TandP_muon->Integral());
    h_signal_dRjet_truth_match_elec->Scale(1. / h_signal_dRjet_truth_match_elec->Integral());
    h_signal_dRjet_truth_match_muon->Scale(1. / h_signal_dRjet_truth_match_muon->Integral());
    h_signal_dRjet_TandP_truth_match_elec->Scale(1. / h_signal_dRjet_TandP_truth_match_elec->Integral());
    h_signal_dRjet_TandP_truth_match_muon->Scale(1. / h_signal_dRjet_TandP_truth_match_muon->Integral());
    h_signal_dRjet_ttbar_elec->Scale(1. / h_signal_dRjet_ttbar_elec->Integral());
    h_signal_dRjet_ttbar_muon->Scale(1. / h_signal_dRjet_ttbar_muon->Integral());
    h_signal_dRjet_Gtt_elec->Scale(1. / h_signal_dRjet_Gtt_elec->Integral());
    h_signal_dRjet_Gtt_muon->Scale(1. / h_signal_dRjet_Gtt_muon->Integral());

    TCanvas *distributions = new TCanvas("distributions", "distributions", 1200, 1200);
    gStyle->SetOptStat(0);
    distributions->Divide(2, 2);

    distributions->cd(1);
    gPad->SetLeftMargin(0.12);
    h_baseline_dRjet_data_elec->SetTitle("");
    h_baseline_dRjet_data_elec->SetXTitle("#Delta R(e, jet)");
    h_baseline_dRjet_data_elec->SetYTitle("Normalized number of baseline probe electrons");
    h_baseline_dRjet_data_elec->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_dRjet_data_elec->SetMaximum(h_baseline_dRjet_ttbar_elec->GetMaximum() * 1.1);
    h_baseline_dRjet_data_elec->GetXaxis()->SetRangeUser(0, 1);
    h_baseline_dRjet_data_elec->SetMarkerStyle(kFullCircle);
    h_baseline_dRjet_data_elec->SetMarkerColor(kBlack);
    h_baseline_dRjet_data_elec->SetLineColor(kBlack);
    h_baseline_dRjet_data_elec->Draw("E1");
    h_baseline_dRjet_TandP_elec->SetMarkerStyle(kFullSquare);
    h_baseline_dRjet_TandP_elec->SetMarkerColor(kRed);
    h_baseline_dRjet_TandP_elec->SetLineColor(kRed);
    h_baseline_dRjet_TandP_elec->Draw("hist,same");
    h_baseline_dRjet_truth_match_elec->SetMarkerStyle(kFullTriangleUp);
    h_baseline_dRjet_truth_match_elec->SetMarkerColor(kBlue);
    h_baseline_dRjet_truth_match_elec->SetLineColor(kBlue);
    h_baseline_dRjet_truth_match_elec->Draw("hist,same");
    h_baseline_dRjet_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_baseline_dRjet_ttbar_elec->SetMarkerColor(kMagenta);
    h_baseline_dRjet_ttbar_elec->SetLineColor(kMagenta);
    h_baseline_dRjet_ttbar_elec->Draw("hist,same");
    h_baseline_dRjet_Gtt_elec->SetMarkerStyle(kFullCross);
    h_baseline_dRjet_Gtt_elec->SetMarkerColor(kOrange);
    h_baseline_dRjet_Gtt_elec->SetLineColor(kOrange);
    h_baseline_dRjet_Gtt_elec->Draw("hist,same");

    TLegend *leg1 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(h_baseline_dRjet_data_elec, "Data");
    leg1->AddEntry(h_baseline_dRjet_TandP_elec, "Z #rightarrow ee T&P MC", "l");
    leg1->AddEntry(h_baseline_dRjet_truth_match_elec, "Z #rightarrow ee truth matched MC", "l");
    leg1->AddEntry(h_baseline_dRjet_ttbar_elec, "ttbar MC", "l");
    leg1->AddEntry(h_baseline_dRjet_Gtt_elec, "Gtt MC", "l");
    leg1->Draw();

    distributions->cd(2);
    gPad->SetLeftMargin(0.12);
    h_baseline_dRjet_data_muon->SetTitle("");
    h_baseline_dRjet_data_muon->SetXTitle("#Delta R(#mu, jet)");
    h_baseline_dRjet_data_muon->SetYTitle("Normalized number of baseline probe muons");
    h_baseline_dRjet_data_muon->GetYaxis()->SetTitleOffset(1.5);
    h_baseline_dRjet_data_muon->SetMaximum(h_baseline_dRjet_ttbar_muon->GetMaximum() * 1.1);
    h_baseline_dRjet_data_muon->GetXaxis()->SetRangeUser(0, 1);
    h_baseline_dRjet_data_muon->SetMarkerStyle(kFullCircle);
    h_baseline_dRjet_data_muon->SetMarkerColor(kBlack);
    h_baseline_dRjet_data_muon->SetLineColor(kBlack);
    h_baseline_dRjet_data_muon->Draw("E1");
    h_baseline_dRjet_TandP_muon->SetMarkerStyle(kFullSquare);
    h_baseline_dRjet_TandP_muon->SetMarkerColor(kRed);
    h_baseline_dRjet_TandP_muon->SetLineColor(kRed);
    h_baseline_dRjet_TandP_muon->Draw("hist,same");
    h_baseline_dRjet_truth_match_muon->SetMarkerStyle(kFullTriangleUp);
    h_baseline_dRjet_truth_match_muon->SetMarkerColor(kBlue);
    h_baseline_dRjet_truth_match_muon->SetLineColor(kBlue);
    h_baseline_dRjet_truth_match_muon->Draw("hist,same");
    h_baseline_dRjet_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_baseline_dRjet_ttbar_muon->SetMarkerColor(kMagenta);
    h_baseline_dRjet_ttbar_muon->SetLineColor(kMagenta);
    h_baseline_dRjet_ttbar_muon->Draw("hist,same");
    h_baseline_dRjet_Gtt_muon->SetMarkerStyle(kFullCross);
    h_baseline_dRjet_Gtt_muon->SetMarkerColor(kOrange);
    h_baseline_dRjet_Gtt_muon->SetLineColor(kOrange);
    h_baseline_dRjet_Gtt_muon->Draw("hist,same");

    TLegend *leg2 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->AddEntry(h_baseline_dRjet_data_muon, "Data");
    leg2->AddEntry(h_baseline_dRjet_TandP_muon, "Z #rightarrow #mu#mu T&P MC", "l");
    leg2->AddEntry(h_baseline_dRjet_truth_match_muon, "Z #rightarrow #mu#mu truth matched MC", "l");
    leg2->AddEntry(h_baseline_dRjet_ttbar_muon, "ttbar MC", "l");
    leg2->AddEntry(h_baseline_dRjet_Gtt_muon, "Gtt MC", "l");
    leg2->Draw("same");

    distributions->cd(3);
    gPad->SetLeftMargin(0.12);
    h_signal_dRjet_data_elec->SetTitle("");
    h_signal_dRjet_data_elec->SetXTitle("#Delta R(e, jet)");
    h_signal_dRjet_data_elec->SetYTitle("Normalized number of signal probe electrons");
    h_signal_dRjet_data_elec->GetYaxis()->SetTitleOffset(1.5);
    h_signal_dRjet_data_elec->SetMaximum(h_signal_dRjet_ttbar_elec->GetMaximum() * 1.1);
    h_signal_dRjet_data_elec->GetXaxis()->SetRangeUser(0, 1);
    h_signal_dRjet_data_elec->SetMarkerStyle(kFullCircle);
    h_signal_dRjet_data_elec->SetMarkerColor(kBlack);
    h_signal_dRjet_data_elec->SetLineColor(kBlack);
    h_signal_dRjet_data_elec->Draw("E1");
    h_signal_dRjet_TandP_elec->SetMarkerStyle(kFullSquare);
    h_signal_dRjet_TandP_elec->SetMarkerColor(kRed);
    h_signal_dRjet_TandP_elec->SetLineColor(kRed);
    h_signal_dRjet_TandP_elec->Draw("hist,same");
    h_signal_dRjet_truth_match_elec->SetMarkerStyle(kFullTriangleUp);
    h_signal_dRjet_truth_match_elec->SetMarkerColor(kBlue);
    h_signal_dRjet_truth_match_elec->SetLineColor(kBlue);
    h_signal_dRjet_truth_match_elec->Draw("hist,same");
    h_signal_dRjet_ttbar_elec->SetMarkerStyle(kFullDiamond);
    h_signal_dRjet_ttbar_elec->SetMarkerColor(kMagenta);
    h_signal_dRjet_ttbar_elec->SetLineColor(kMagenta);
    h_signal_dRjet_ttbar_elec->Draw("hist,same");
    h_signal_dRjet_Gtt_elec->SetMarkerStyle(kFullCross);
    h_signal_dRjet_Gtt_elec->SetMarkerColor(kOrange);
    h_signal_dRjet_Gtt_elec->SetLineColor(kOrange);
    h_signal_dRjet_Gtt_elec->Draw("hist,same");

    TLegend *leg3 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.03);
    leg3->AddEntry(h_signal_dRjet_data_elec, "Data");
    leg3->AddEntry(h_signal_dRjet_TandP_elec, "Z #rightarrow ee T&P MC", "l");
    leg3->AddEntry(h_signal_dRjet_truth_match_elec, "Z #rightarrow ee truth matched MC", "l");
    leg3->AddEntry(h_signal_dRjet_ttbar_elec, "ttbar MC", "l");
    leg3->AddEntry(h_signal_dRjet_Gtt_elec, "Gtt MC", "l");
    leg3->Draw();

    distributions->cd(4);
    gPad->SetLeftMargin(0.12);
    h_signal_dRjet_data_muon->SetTitle("");
    h_signal_dRjet_data_muon->SetXTitle("#Delta R(#mu, jet)");
    h_signal_dRjet_data_muon->SetYTitle("Normalized number of signal probe muons");
    h_signal_dRjet_data_muon->GetYaxis()->SetTitleOffset(1.5);
    h_signal_dRjet_data_muon->SetMaximum(h_signal_dRjet_ttbar_muon->GetMaximum() * 1.1);
    h_signal_dRjet_data_muon->GetXaxis()->SetRangeUser(0, 1);
    h_signal_dRjet_data_muon->SetMarkerStyle(kFullCircle);
    h_signal_dRjet_data_muon->SetMarkerColor(kBlack);
    h_signal_dRjet_data_muon->SetLineColor(kBlack);
    h_signal_dRjet_data_muon->Draw("E1");
    h_signal_dRjet_TandP_muon->SetMarkerStyle(kFullSquare);
    h_signal_dRjet_TandP_muon->SetMarkerColor(kRed);
    h_signal_dRjet_TandP_muon->SetLineColor(kRed);
    h_signal_dRjet_TandP_muon->Draw("hist,same");
    h_signal_dRjet_truth_match_muon->SetMarkerStyle(kFullTriangleUp);
    h_signal_dRjet_truth_match_muon->SetMarkerColor(kBlue);
    h_signal_dRjet_truth_match_muon->SetLineColor(kBlue);
    h_signal_dRjet_truth_match_muon->Draw("hist,same");
    h_signal_dRjet_ttbar_muon->SetMarkerStyle(kFullDiamond);
    h_signal_dRjet_ttbar_muon->SetMarkerColor(kMagenta);
    h_signal_dRjet_ttbar_muon->SetLineColor(kMagenta);
    h_signal_dRjet_ttbar_muon->Draw("hist,same");
    h_signal_dRjet_Gtt_muon->SetMarkerStyle(kFullCross);
    h_signal_dRjet_Gtt_muon->SetMarkerColor(kOrange);
    h_signal_dRjet_Gtt_muon->SetLineColor(kOrange);
    h_signal_dRjet_Gtt_muon->Draw("hist,same");

    TLegend *leg4 = new TLegend(0.4, 0.6, 0.9, 0.9);
    leg4->SetBorderSize(0);
    leg4->SetFillColor(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.03);
    leg4->AddEntry(h_signal_dRjet_data_muon, "Data");
    leg4->AddEntry(h_signal_dRjet_TandP_muon, "Z #rightarrow #mu#mu T&P MC", "l");
    leg4->AddEntry(h_signal_dRjet_truth_match_muon, "Z #rightarrow #mu#mu truth matched MC", "l");
    leg4->AddEntry(h_signal_dRjet_ttbar_muon, "ttbar MC", "l");
    leg4->AddEntry(h_signal_dRjet_Gtt_muon, "Gtt MC", "l");
    leg4->Draw("same");

    distributions->SaveAs("dR_baseline_and_signal_distributions.pdf", "pdf");
}
