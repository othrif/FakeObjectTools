#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>

//#include "AtlasStyle/AtlasStyle.C"
//#include "AtlasStyle/AtlasLabels.C"
//#include "AtlasStyle/AtlasUtils.C"

#include <iostream>
using namespace std;

// User defined fitting
// Using a exponential function to fit: baseline: [60, 80]U[100, 120]
Double_t fit_exp_base(double *x, double *par)
{
    if (x[0] > 80. && x[0] < 100.) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]*TMath::Exp(-par[1]*(x[0]-60.));
}

// Using a exponential function to fit: range1: [60, 70]U[100, 120]
Double_t fit_exp_var1(double *x, double *par)
{
    if (x[0] > 70. && x[0] < 100.) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]*TMath::Exp(-par[1]*(x[0]-60.));
}

// Using a exponential function to fit: range2: [65, 75]U[100, 120]
Double_t fit_exp_var2(double *x, double *par)
{
    if (x[0] > 75. && x[0] < 100.) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0]*TMath::Exp(-par[1]*(x[0]-60.));
}

// Because we exclude the Z mass peak in the fitting function, we need to define a new function including Z mass peak.
Double_t fit_exp(double *x, double *par)
{
    return par[0]*TMath::Exp(-par[1]*(x[0]-60.));
}

void ytBackground_subtraction(TString trigger = "single_lepton_trigger",
                              TString template_type = "baseline", // baseline, template1, template2
                              TString fitting_range = "range_baseline", // range_baseline, range1, range2
                              bool mc_scale_to_data = false,
                              double fitting_range_low = 60.,
                              double fitting_range_up = 120.,
                              //double fitting_range_exclude_low = 80.,
                              //double fitting_range_exclude_up = 100.,
                              double mll_window_low = 80.,
                              double mll_window_up = 100.,
                              int pt_bin_low = 0,
                              int pt_bin_up = -1,
                              int eta_bin_low = 0,
                              int eta_bin_up = -1,
                              bool truth_match = false)
{
    // pt range: 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200
    // eta range: 0., 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47
    // So eta bins use (1, 3), (4,5), and (7, 8)

    bool debug = false;
    //bool debug = true;

    cout << "Notice: Remember to modify the input file path!" << endl;

    cout << "Current template = " << template_type << endl;

    TString path = "~/Desktop/20170310/";
    TFile *data_file;
    if (trigger == "single_lepton_trigger")
        data_file = TFile::Open(path + "hist-RLE-merged-data-elec.root");
    else if (trigger == "tag_trigger_matched" || trigger == "single_lepton_trigger_tag_trigger_matched")
        data_file = TFile::Open(path + "hist-RLE-merged-data-elec_tag_trigger_matched.root");
    else if (trigger == "dilepton_trigger")
        data_file = TFile::Open(path + "hist-RLE-merged-data-elec_dilepton_trigger.root");
    else if (trigger == "dilepton_trigger_tag_trigger_matched")
        data_file = TFile::Open(path + "hist-RLE-merged-data-elec_dilepton_trigger_tag_trigger_matched.root");

    TFile *mc_file;
    if (!truth_match)  // T&P
        mc_file = TFile::Open(path + "RLE_MC_Zee/hist-20170310.root");
    else // truth_match
        mc_file = TFile::Open(path + "RLE_MC_Zee_truth_match/hist-20170310.root");
        // TandP + truth_match
        //mc_file = TFile::Open(path + "RLE_MC_Zmumu_TandP_truth_match/hist-20170112.root");

    // Declare 1-dim histograms
    TH1F *data_baseline_mll;
    TH1F *data_signal_mll;
    TH1F *data_bkg_template;

    TH1F *mc_baseline_mll;
    TH1F *mc_signal_mll;
    TH1F *mc_bkg_template;

    if (pt_bin_low == 0 && pt_bin_up == -1) { // Use full range of pT, so using 1-dim histogram is easier.
        // Get 1-dim data histograms
        data_baseline_mll = (TH1F *)data_file->Get("h_baseline_mll");
        data_signal_mll = (TH1F *)data_file->Get("h_signal_mll");
        if (template_type == "baseline")
            data_bkg_template = (TH1F *)data_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso");
        else if (template_type == "template1")
            data_bkg_template = (TH1F *)data_file->Get("h_bkg_template_fail_CaloIso_and_TrackIso");
        else if (template_type == "template2")
            data_bkg_template = (TH1F *)data_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight");

        // Get 1-dim MC histograms
        mc_baseline_mll = (TH1F *)mc_file->Get("h_baseline_mll");
        mc_signal_mll = (TH1F *)mc_file->Get("h_signal_mll");
        if (template_type == "baseline")
            mc_bkg_template = (TH1F *)mc_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso");
        else if (template_type == "template1")
            mc_bkg_template = (TH1F *)mc_file->Get("h_bkg_template_fail_CaloIso_and_TrackIso");
        else if (template_type == "template2")
            mc_bkg_template = (TH1F *)mc_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight");
    }
    else {
        // Get 3-dim data histogram
        TH3F *data_baseline_mll_3dim = (TH3F *)data_file->Get("h_baseline_pt_eta_mll");
        TH3F *data_signal_mll_3dim = (TH3F *)data_file->Get("h_signal_pt_eta_mll");
        TH3F *data_bkg_template_3dim;
        if (template_type == "baseline")
            data_bkg_template_3dim = (TH3F *)data_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll");
        else if (template_type == "template1")
            data_bkg_template_3dim = (TH3F *)data_file->Get("h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll");
        else if (template_type == "template2")
            data_bkg_template_3dim = (TH3F *)data_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll");
        // project on Z axis
        data_baseline_mll = (TH1F *)data_baseline_mll_3dim->ProjectionZ("data_baseline", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
        data_signal_mll = (TH1F *)data_signal_mll_3dim->ProjectionZ("data_signal", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
        data_bkg_template = (TH1F *)data_bkg_template_3dim->ProjectionZ("data_template", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);

        cout << "pT range:" << data_signal_mll_3dim->GetXaxis()->GetBinLowEdge(pt_bin_low) << " to " << data_signal_mll_3dim->GetXaxis()->GetBinUpEdge(pt_bin_up) << endl;
        cout << "eta range: " << data_signal_mll_3dim->GetYaxis()->GetBinLowEdge(eta_bin_low) << " to " << data_signal_mll_3dim->GetYaxis()->GetBinUpEdge(eta_bin_up) << endl;

        // Get 3-dim MC histograms
        TH3F *mc_baseline_mll_3dim = (TH3F *)mc_file->Get("h_baseline_pt_eta_mll");
        TH3F *mc_signal_mll_3dim = (TH3F *)mc_file->Get("h_signal_pt_eta_mll");
        TH3F *mc_bkg_template_3dim;
        if (template_type == "baseline")
            mc_bkg_template_3dim = (TH3F *)mc_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll");
        else if (template_type == "template1")
            mc_bkg_template_3dim = (TH3F *)mc_file->Get("h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll");
        else if (template_type == "template2")
            mc_bkg_template_3dim = (TH3F *)mc_file->Get("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll");
        // project on Z axis
        mc_baseline_mll = (TH1F *)mc_baseline_mll_3dim->ProjectionZ("mc_baseline", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
        mc_signal_mll = (TH1F *)mc_signal_mll_3dim->ProjectionZ("mc_signal", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
        mc_bkg_template = (TH1F *)mc_bkg_template_3dim->ProjectionZ("mc_template", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up);
    }

    if (debug) {
        cout << data_baseline_mll << endl;
        cout << data_signal_mll << endl;
        cout << data_bkg_template << endl;

        cout << mc_baseline_mll << endl;
        cout << mc_signal_mll << endl;
        cout << mc_bkg_template << endl;

        TCanvas *cdata = new TCanvas("cdata", "", 600, 600);
        data_baseline_mll->Draw();
        data_signal_mll->SetLineColor(kRed);
        data_signal_mll->Draw("same");

        TCanvas *cmc = new TCanvas("cmc", "", 600, 600);
        mc_baseline_mll->Draw();
        mc_signal_mll->SetLineColor(kRed);
        mc_signal_mll->Draw("same");
    }

    // Tail range [120, 150]
    double tail_low = 120.;
    double tail_up = 150.;

    // Get the corresponding bin
    // Use +-0.01 to slightly deviate from bin boundary
    int fitting_range_low_bin = data_bkg_template->FindBin(fitting_range_low + 0.01);
    int fitting_range_up_bin = data_bkg_template->FindBin(fitting_range_up - 0.01);
    //int fitting_range_exclude_low_bin = data_bkg_template->FindBin(fitting_range_exclude_low - 0.01);
    //int fitting_range_exclude_up_bin = data_bkg_template->FindBin(fitting_range_exclude_up + 0.01);

    int mll_window_low_bin = data_baseline_mll->FindBin(mll_window_low + 0.01);
    int mll_window_up_bin = data_baseline_mll->FindBin(mll_window_up - 0.01);

    int tail_low_bin = data_baseline_mll->FindBin(tail_low + 0.01);
    int tail_up_bin = data_baseline_mll->FindBin(tail_up - 0.01);

    if (debug) {
        cout << "fitting_range_low_bin=" << fitting_range_low_bin << endl;
        cout << "fitting_range_up_bin=" << fitting_range_up_bin << endl;
        //cout << "fitting_range_exclude_low_bin=" << fitting_range_exclude_low_bin << endl;
        //cout << "fitting_range_exclude_up_bin=" << fitting_range_exclude_up_bin << endl;

        cout << "mll_window_low_bin=" << mll_window_low_bin << endl;
        cout << "mll_window_up_bin=" << mll_window_up_bin << endl;

        cout << "tail_low_bin=" << tail_low_bin << endl;
        cout << "tail_up_bin=" << tail_up_bin << endl;
    }

    //
    // Scale MC to data using a Gaussian fit of the Z peak (85 < mll < 95)
    //
    if (mc_scale_to_data == true) {
        // baseline
        data_baseline_mll->Fit("gaus", "0", "", 85, 95);
        TF1 *func_ee_baseline = data_baseline_mll->GetFunction("gaus");
        double data_baseline_mll_peak_area = func_ee_baseline->Integral(85., 95.);
        double mc_baseline_mll_peak_area = mc_baseline_mll->Integral(mc_baseline_mll->GetXaxis()->FindBin(85. + 0.01), mc_baseline_mll->GetXaxis()->FindBin(95. - 0.01));
        mc_baseline_mll->Scale(data_baseline_mll_peak_area / mc_baseline_mll_peak_area);
        // signal
        data_signal_mll->Fit("gaus", "0", "", 85, 95);
        TF1 *func_ee_signal = data_signal_mll->GetFunction("gaus");
        double data_signal_mll_peak_area = func_ee_signal->Integral(85., 95.);
        double mc_signal_mll_peak_area = mc_signal_mll->Integral(mc_signal_mll->GetXaxis()->FindBin(85. + 0.01), mc_signal_mll->GetXaxis()->FindBin(95. - 0.01));
        mc_signal_mll->Scale(data_signal_mll_peak_area / mc_signal_mll_peak_area);
    }

    //
    // Number of events in the Zmass peak (80 < mll < 100)
    // Data
    double data_baseline_mll_events_in_the_Z_peak = data_baseline_mll->Integral(mll_window_low_bin, mll_window_up_bin);
    double data_signal_mll_events_in_the_Z_peak = data_signal_mll->Integral(mll_window_low_bin, mll_window_up_bin);
    double data_signal_to_baseline_ratio_in_the_Z_peak = data_signal_mll_events_in_the_Z_peak / data_baseline_mll_events_in_the_Z_peak;
    double data_non_signal_events_in_the_Z_peak = data_baseline_mll_events_in_the_Z_peak - data_signal_mll_events_in_the_Z_peak;

    // Data template
    double data_bkg_template_in_the_Z_peak = data_bkg_template->Integral(mll_window_low_bin, mll_window_up_bin);

    // MC
    double mc_baseline_mll_events_in_the_Z_peak = mc_baseline_mll->Integral(mll_window_low_bin, mll_window_up_bin);
    double mc_signal_mll_events_in_the_Z_peak = mc_signal_mll->Integral(mll_window_low_bin, mll_window_up_bin);
    double mc_signal_to_baseline_ratio_in_the_Z_peak = mc_signal_mll_events_in_the_Z_peak / mc_baseline_mll_events_in_the_Z_peak;
    double mc_non_signal_events_in_the_Z_peak = mc_baseline_mll_events_in_the_Z_peak - mc_signal_mll_events_in_the_Z_peak;

    if (debug) {
        cout << "data_baseline_mll_events_in_the_Z_peak=" << data_baseline_mll_events_in_the_Z_peak << endl;
        cout << "data_signal_mll_events_in_the_Z_peak=" << data_signal_mll_events_in_the_Z_peak << endl;
        cout << "data_signal_to_baseline_ratio_in_the_Z_peak=" << data_signal_to_baseline_ratio_in_the_Z_peak << endl;
        cout << "data_non_signal_events_in_the_Z_peak=" << data_non_signal_events_in_the_Z_peak << endl;

        cout << "data_bkg_template_in_the_Z_peak=" << data_bkg_template_in_the_Z_peak << endl;

        cout << "mc_baseline_mll_events_in_the_Z_peak=" << mc_baseline_mll_events_in_the_Z_peak << endl;
        cout << "mc_signal_mll_events_in_the_Z_peak=" << mc_signal_mll_events_in_the_Z_peak << endl;
        cout << "mc_signal_to_baseline_ratio_in_the_Z_peak=" << mc_signal_to_baseline_ratio_in_the_Z_peak << endl;
        cout << "mc_non_signal_events_in_the_Z_peak=" << mc_non_signal_events_in_the_Z_peak << endl;
    }

    //
    // Number of events in the tail (120 < mll < 150)
    // Data
    double data_baseline_mll_events_in_the_tail = data_baseline_mll->Integral(tail_low_bin, tail_up_bin);
    double data_signal_mll_events_in_the_tail = data_signal_mll->Integral(tail_low_bin, tail_up_bin);
    double data_signal_to_baseline_ratio_in_the_tail = data_signal_mll_events_in_the_tail / data_baseline_mll_events_in_the_tail;
    double data_non_signal_events_in_the_tail = data_baseline_mll_events_in_the_tail - data_signal_mll_events_in_the_tail;

    // Data template
    double data_bkg_template_in_the_tail = data_bkg_template->Integral(tail_low_bin, tail_up_bin);

    // MC
    double mc_baseline_mll_events_in_the_tail = mc_baseline_mll->Integral(tail_low_bin, tail_up_bin);
    double mc_signal_mll_events_in_the_tail = mc_signal_mll->Integral(tail_low_bin, tail_up_bin);
    double mc_signal_to_baseline_ratio_in_the_tail = mc_signal_mll_events_in_the_tail / mc_baseline_mll_events_in_the_tail;
    double mc_non_signal_events_in_the_tail = mc_baseline_mll_events_in_the_tail - mc_signal_mll_events_in_the_tail;

    if (debug) {
        cout << "data_baseline_mll_events_in_the_tail=" << data_baseline_mll_events_in_the_tail << endl;
        cout << "data_signal_mll_events_in_the_tail=" << data_signal_mll_events_in_the_tail << endl;
        cout << "data_signal_to_baseline_ratio_in_the_tail=" << data_signal_to_baseline_ratio_in_the_tail << endl;
        cout << "data_non_signal_events_in_the_tail=" << data_non_signal_events_in_the_tail << endl;

        cout << "data_bkg_template_in_the_tail=" << data_bkg_template_in_the_tail << endl;

        cout << "mc_baseline_mll_events_in_the_tail=" << mc_baseline_mll_events_in_the_tail << endl;
        cout << "mc_signal_mll_events_in_the_tail=" << mc_signal_mll_events_in_the_tail << endl;
        cout << "mc_signal_to_baseline_ratio_in_the_tail=" << mc_signal_to_baseline_ratio_in_the_tail << endl;
        cout << "mc_non_signal_events_in_the_tail=" << mc_non_signal_events_in_the_tail << endl;
    }

    //
    // Define fitting function
    // We want to fit from 60 to 120 and exclude the mll window
    TF1 *fit_function;
    if (fitting_range == "range_baseline")
        fit_function = new TF1("fit_range_baseline", fit_exp_base, fitting_range_low, fitting_range_up, 2);
    else if (fitting_range == "range1")
        fit_function = new TF1("fit_range1", fit_exp_var1, fitting_range_low, fitting_range_up, 2);
    else if (fitting_range == "range2")
        fit_function = new TF1("fit_range2", fit_exp_var2, fitting_range_low, fitting_range_up, 2);
    fit_function->SetParLimits(0, 0, 100000); // force the par[0] to be positive ranging from 0 to 10^5.
    // func_Z_peak is used to make plot and calculate Z peak region.
    // Because fit_function excludes the Z peak region, we need build a new function include Z peak region.
    TF1 *func_Z_peak = new TF1("func_Z_peak", fit_exp, fitting_range_low, tail_up, 2);

    TH1F *h_fit_data_bkg_template = (TH1F *)data_bkg_template->Clone();
    if (debug)
        h_fit_data_bkg_template->Draw();
    h_fit_data_bkg_template->Sumw2();

    h_fit_data_bkg_template->Fit(fit_function, "0");
    h_fit_data_bkg_template->Fit(fit_function, "0");
    h_fit_data_bkg_template->Fit(fit_function, "0"); // fit 3 times to get better fitting results.

    // func_Z_peak should have exactly the same parameters as fitting function.
    func_Z_peak->SetParameters(fit_function->GetParameter(0), fit_function->GetParameter(1));
    if (debug) {
        func_Z_peak->SetLineColor(kYellow);
        func_Z_peak->Draw("same");
    }

    double fit_data_bkg_template_in_the_Z_peak = func_Z_peak->Integral(mll_window_low, mll_window_up);
    double fit_data_bkg_template_in_the_tail = fit_function->Integral(tail_low, tail_up);

    //
    // Statistical information
    //
    cout << "********** Statistical information **********" << endl;
    cout << "Data in the Z peak:" << endl;
    cout << "Baseline in the Z peak = " << data_baseline_mll_events_in_the_Z_peak << " +- " << sqrt(data_baseline_mll_events_in_the_Z_peak)
    << "(" << 100 * sqrt(data_baseline_mll_events_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak << "%)" << endl;
    cout << "Signal in the Z peak = " << data_signal_mll_events_in_the_Z_peak << " +- " << sqrt(data_signal_mll_events_in_the_Z_peak)
    << "(" << 100 * sqrt(data_signal_mll_events_in_the_Z_peak) / data_signal_mll_events_in_the_Z_peak << "%)" << endl;
    cout << "Signal/Baseline ratio in the Z peak = " << data_signal_to_baseline_ratio_in_the_Z_peak << " +- " << sqrt(data_signal_to_baseline_ratio_in_the_Z_peak * (1 - data_signal_to_baseline_ratio_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak) << endl;
    cout << "(Baseline - Signal) in the Z peak = " << data_non_signal_events_in_the_Z_peak << " +- " << sqrt(data_non_signal_events_in_the_Z_peak)
    << "(" << 100 * sqrt(data_non_signal_events_in_the_Z_peak) / data_non_signal_events_in_the_Z_peak << "%)" << endl;

    cout << "Template in the Z peak:" << endl;
    cout << "Template baseline in the Z peak = " << data_bkg_template_in_the_Z_peak << " +- " << sqrt(data_bkg_template_in_the_Z_peak)
    << "(" << 100 * sqrt(data_bkg_template_in_the_Z_peak) / data_bkg_template_in_the_Z_peak << "%)" << endl;

    cout << "Template in the Z peak (fitting):" << endl;
    cout << "Template baseline in the Z peak = " << fit_data_bkg_template_in_the_Z_peak << " +- " << sqrt(fit_data_bkg_template_in_the_Z_peak)
    << "(" << 100 * sqrt(fit_data_bkg_template_in_the_Z_peak) / fit_data_bkg_template_in_the_Z_peak << "%)" << endl;

    cout << "MC in the Z peak:" << endl;
    cout << "Baseline in the Z peak = " << mc_baseline_mll_events_in_the_Z_peak << " +- " << sqrt(mc_baseline_mll_events_in_the_Z_peak)
    << "(" << 100 * sqrt(mc_baseline_mll_events_in_the_Z_peak) / mc_baseline_mll_events_in_the_Z_peak << "%)" << endl;
    cout << "Signal in the Z peak = " << mc_signal_mll_events_in_the_Z_peak << " +- " << sqrt(mc_signal_mll_events_in_the_Z_peak)
    << "(" << 100 * sqrt(mc_signal_mll_events_in_the_Z_peak) / mc_signal_mll_events_in_the_Z_peak << "%)" << endl;
    cout << "Signal/Baseline ratio in the Z peak = " << mc_signal_to_baseline_ratio_in_the_Z_peak << " +- " << sqrt(mc_signal_to_baseline_ratio_in_the_Z_peak * (1 - mc_signal_to_baseline_ratio_in_the_Z_peak) / mc_baseline_mll_events_in_the_Z_peak) << endl;
    cout << "(Baseline - Signal) in the Z peak = " << mc_non_signal_events_in_the_Z_peak <<  " +- " << sqrt(mc_non_signal_events_in_the_Z_peak)
    << "(" << 100 * sqrt(mc_non_signal_events_in_the_Z_peak) / mc_non_signal_events_in_the_Z_peak << "%)" << endl;
    cout << "********************" << endl;

    cout << "Data in the tail:" << endl;
    cout << "Baseline in the tail = " << data_baseline_mll_events_in_the_tail << " +- " << sqrt(data_baseline_mll_events_in_the_tail)
    << "(" << 100 * sqrt(data_baseline_mll_events_in_the_tail) / data_baseline_mll_events_in_the_tail << "%)" << endl;
    cout << "Signal in the tail = " << data_signal_mll_events_in_the_tail << " +- " << sqrt(data_signal_mll_events_in_the_tail)
    << "(" << 100 * sqrt(data_signal_mll_events_in_the_tail) / data_signal_mll_events_in_the_tail << "%)" << endl;
    cout << "Signal/Baseline ratio in the tail = " << data_signal_to_baseline_ratio_in_the_tail << " +- " << sqrt(data_signal_to_baseline_ratio_in_the_tail * (1 - data_signal_to_baseline_ratio_in_the_tail) / data_baseline_mll_events_in_the_tail) << endl;
    cout << "(Baseline - Signal) in the tail = " << data_non_signal_events_in_the_tail << " +- " << sqrt(data_non_signal_events_in_the_tail)
    << "(" << 100 * sqrt(data_non_signal_events_in_the_tail) / data_non_signal_events_in_the_tail << "%)" << endl;
    
    cout << "Template in the tail:" << endl;
    cout << "Template baseline in the tail = " << data_bkg_template_in_the_tail << " +- " << sqrt(data_bkg_template_in_the_tail)
    << "(" << 100 * sqrt(data_bkg_template_in_the_tail) / data_bkg_template_in_the_tail << "%)" << endl;
    
    cout << "Template in the tail (fitting):" << endl;
    cout << "Template baseline in the tail = " << fit_data_bkg_template_in_the_tail << " +- " << sqrt(fit_data_bkg_template_in_the_tail)
    << "(" << 100 * sqrt(fit_data_bkg_template_in_the_tail) / fit_data_bkg_template_in_the_tail << "%)" << endl;
    
    cout << "MC in the tail:" << endl;
    cout << "Baseline in the tail = " << mc_baseline_mll_events_in_the_tail << " +- " << sqrt(mc_baseline_mll_events_in_the_tail)
    << "(" << 100 * sqrt(mc_baseline_mll_events_in_the_tail) / mc_baseline_mll_events_in_the_tail << "%)" << endl;
    cout << "Signal in the tail = " << mc_signal_mll_events_in_the_tail << " +- " << sqrt(mc_signal_mll_events_in_the_tail)
    << "(" << 100 * sqrt(mc_signal_mll_events_in_the_tail) / mc_signal_mll_events_in_the_tail << "%)" << endl;
    cout << "Signal/Baseline ratio in the tail = " << mc_signal_to_baseline_ratio_in_the_tail << " +- " << sqrt(mc_signal_to_baseline_ratio_in_the_tail * (1 - mc_signal_to_baseline_ratio_in_the_tail) / mc_baseline_mll_events_in_the_tail) << endl;
    cout << "(Baseline - Signal) in the tail = " << mc_non_signal_events_in_the_tail <<  " +- " << sqrt(mc_non_signal_events_in_the_tail)
    << "(" << 100 * sqrt(mc_non_signal_events_in_the_tail) / mc_non_signal_events_in_the_tail << "%)" << endl;
    cout << "********************" << endl;

    //
    // Estimate the signal contamination
    //
    // This part is exactly the same as line 235 to 238
    //double mc_baseline_Nevents_in_mll_window = mc_baseline_mll->Integral(mll_window_low_bin, mll_window_up_bin);
    //double mc_signal_Nevents_in_mll_window = mc_signal_mll->Integral(mll_window_low_bin, mll_window_up_bin);
    //double mc_real_eff_in_mll_window = mc_signal_Nevents_in_mll_window / mc_baseline_Nevents_in_mll_window;
    //double mc_real_eff_in_mll_window_err = sqrt(mc_real_eff_in_mll_window * (1 - mc_real_eff_in_mll_window) / mc_baseline_Nevents_in_mll_window);

    //double mc_baseline_Nevents_in_the_tail;
    //if (mc_signal_mll_events_in_the_tail != 0) {
    //    mc_baseline_Nevents_in_the_tail = mc_signal_mll_events_in_the_tail / mc_real_eff_in_mll_window;
    //    //mc_baseline_Nevents_in_the_tail = mc_signal_mll_events_in_the_tail / mc_signal_to_baseline_ratio_in_the_tail;
    //}
    //else {
        //mc_baseline_Nevents_in_the_tail = mc_baseline_mll_events_in_the_tail;
    //}
    //double data_estimated_background_Nevents_in_the_tail = data_baseline_mll_events_in_the_tail - mc_baseline_Nevents_in_the_tail;
    double data_estimated_background_Nevents_in_the_tail = data_baseline_mll_events_in_the_tail - mc_baseline_mll_events_in_the_tail;
    double data_estimated_background_Nevents_in_the_tail_err = sqrt(data_baseline_mll_events_in_the_tail + mc_baseline_mll_events_in_the_tail);
    //double signal_contimation_in_the_tail = mc_baseline_Nevents_in_the_tail / data_baseline_mll_events_in_the_tail;
    double signal_contimation_in_the_tail = mc_baseline_mll_events_in_the_tail / data_baseline_mll_events_in_the_tail;
    double signal_contimation_in_the_tail_err = signal_contimation_in_the_tail * sqrt(1. / mc_baseline_mll_events_in_the_tail + 1. / data_baseline_mll_events_in_the_tail);

    if (debug) {
        //cout << "mc_baseline_Nevents_in_mll_window = " << mc_baseline_Nevents_in_mll_window << endl;
        //cout << "mc_signal_Nevents_in_mll_window = " << mc_signal_Nevents_in_mll_window << endl;
        //cout << "mc_real_eff_in_mll_window = " << mc_real_eff_in_mll_window << " +- " << mc_real_eff_in_mll_window_err << endl;

        //cout << "mc_baseline_Nevents_in_the_tail = " << mc_baseline_Nevents_in_the_tail << endl;
        cout << "data_estimated_background_Nevents_in_the_tail = " << data_estimated_background_Nevents_in_the_tail << endl;
        cout << "data_estimated_background_Nevents_in_the_tail_err = " << data_estimated_background_Nevents_in_the_tail_err << endl;
        cout << "signal_contimation_in_the_tail = " << signal_contimation_in_the_tail << " +- " << signal_contimation_in_the_tail_err << endl;
    }

    //
    // Normalize the background template in the tail to data in the tail
    //
    double norm = data_estimated_background_Nevents_in_the_tail / fit_data_bkg_template_in_the_tail;
    double norm_err = norm * sqrt(pow(data_estimated_background_Nevents_in_the_tail_err / data_estimated_background_Nevents_in_the_tail, 2) + 1. / fit_data_bkg_template_in_the_tail);

    // When 10 GeV < pT < 15 GeV and 0 < eta < 0.8, I get a negative norm value because MC baseline in the tail
    // has more events than data estimated bkg in the tail. So I have to use 100 < mll < 120 tail region to
    // calculate the norm
    // If I use template2 and 10 GeV < pT < 15 GeV and 0.8 < eta < 1.37, the norm is very large. So I also use the
    // 100 < mll < 120 tail region to calculate the norm.
    if (norm < 0 || norm > 20) {
        int special_tail_low_bin = data_baseline_mll->FindBin(60. + 0.01);
        int special_tail_up_bin = data_baseline_mll->FindBin(70. - 0.01);
        double a = data_baseline_mll->Integral(special_tail_low_bin, special_tail_up_bin);
        double b = mc_baseline_mll->Integral(special_tail_low_bin, special_tail_up_bin);
        double c = fit_function->Integral(60., 70.);
        
        norm = (a - b) / c;
        // uncertainty of a - b: da = sqrt(a), db = sqrt(b)
        double a_minus_b = a - b;
        double d_a_minus_b = sqrt(a + b);
        norm_err = norm * sqrt(pow(d_a_minus_b / a_minus_b, 2)+ 1. / c);

        if (debug) {
            cout << "data_baseline_mll_events_in_the_60_mll_70_region=" << a << endl;
            cout << "mc_baseline_mll_events_in_the_60_mll_70_region=" << b << endl;
            cout << "fit_data_bkg_template_in_the_60_mll_70_region=" << c << endl;
        }
    }

    if (debug) {
        cout << "norm = " << norm << endl;
        cout << "norm_err = " << norm_err << endl;
    }
    data_bkg_template->Scale(norm);
    h_fit_data_bkg_template->Scale(norm);

    //
    // Calculate the background in the mll window
    //
    h_fit_data_bkg_template->Fit(fit_function, "0"); // Since we normalized the histogram, we need to do fit again.
    func_Z_peak->SetParameters(fit_function->GetParameter(0), fit_function->GetParameter(1)); // We also need to set the new parameters to func_Z_peak
    fit_data_bkg_template_in_the_Z_peak = func_Z_peak->Integral(mll_window_low, mll_window_up); // This is the bkg under Z peak estimated by the template
    double fit_data_bkg_template_in_the_Z_peak_err = sqrt(fit_data_bkg_template_in_the_Z_peak);

    double N_baseline = data_baseline_mll_events_in_the_Z_peak - fit_data_bkg_template_in_the_Z_peak;
    double N_signal = data_signal_mll_events_in_the_Z_peak;
    double real_eff = N_signal / N_baseline;

    double N_baseline_err = sqrt(data_baseline_mll_events_in_the_Z_peak + fit_data_bkg_template_in_the_Z_peak);
    double N_signal_err = sqrt(data_signal_mll_events_in_the_Z_peak);
    double real_eff_err = real_eff * sqrt(pow(N_baseline_err / N_baseline, 2) + pow(N_signal_err / N_signal, 2));

    if (debug) {
        cout << "fit_data_bkg_template_in_the_Z_peak = " << fit_data_bkg_template_in_the_Z_peak << endl;
        cout << "fit_data_bkg_template_in_the_Z_peak_err = " << fit_data_bkg_template_in_the_Z_peak_err << endl;

        cout << "N_baseline=" << N_baseline << endl;
        cout << "N_baseline_err=" << N_baseline_err << endl;
        cout << "N_signal=" << N_signal << endl;
        cout << "N_signal_err=" << N_signal_err << endl;
        cout << "real_eff=" << real_eff << endl;
        cout << "real_eff_err=" << real_eff_err << endl;
    }

    cout << "Background:" << endl;
    cout << "Background in the tail (estimated) = " << data_estimated_background_Nevents_in_the_tail << " +- " << data_estimated_background_Nevents_in_the_tail_err << endl;
    cout << "Signal contimation in the tail = " << signal_contimation_in_the_tail << " +- " << signal_contimation_in_the_tail_err << endl;
    cout << "Background in the Z peak = " << fit_data_bkg_template_in_the_Z_peak << " +- " << fit_data_bkg_template_in_the_Z_peak_err << endl;
    cout << "********************" << endl;
    cout << "Final results:" << endl;
    cout << "N(baseline) before background subtraction = " << data_baseline_mll_events_in_the_Z_peak << " +- " << sqrt(data_baseline_mll_events_in_the_Z_peak)
         << "(" << 100 * sqrt(data_baseline_mll_events_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak << "%)" << endl;
    cout << "N(baseline) after background subtraction = " << N_baseline << " +- " << N_baseline_err << "(" << 100 * N_baseline_err / N_baseline << "%)"<< endl;
    cout << "N(signal) = " << N_signal << " +- " << N_signal_err << "(" << 100 * N_signal_err / N_signal << "%)" << endl;
    cout << "real lepton efficiency (before background subtraction) = " << data_signal_to_baseline_ratio_in_the_Z_peak << " +- " << sqrt(data_signal_to_baseline_ratio_in_the_Z_peak * (1 - data_signal_to_baseline_ratio_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak) << endl;
    cout << "real lepton efficiency (after background subtraction) = " << real_eff << " +- " << real_eff_err << endl;
    cout << "********************" << endl;

    //
    // Making ratio plots
    //
    TCanvas *c1 = new TCanvas("c1", "Mll", 600, 600);

    //Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetRightMargin(0.08);
    //pad1->SetGridy(); // grid lines
    pad1->SetLogy();
    pad1->Draw();
    
    // lower plot will be in pad
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.35);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetRightMargin(0.08);
    pad2->SetGridy(); // grid lines
    //pad2->SetLogx();
    pad2->Draw();
    
    pad1->cd(); // pad1 becomes the current pad
    //pad1->SetFrameLineWidth(2);
    
    // Draw curve here
    data_baseline_mll->SetTitle("");
    data_baseline_mll->SetXTitle("M_{ll} [GeV]");
    data_baseline_mll->SetYTitle("Events");
    data_baseline_mll->SetMarkerColor(kBlack);
    data_baseline_mll->SetMarkerStyle(kFullCircle);
    data_baseline_mll->SetMarkerSize(0.5);
    data_baseline_mll->SetLineColor(kBlack);
    data_baseline_mll->SetMaximum(data_baseline_mll->GetMaximum() * 10);
    data_baseline_mll->SetMinimum(0.11);
    data_baseline_mll->SetStats(kFALSE);
    data_baseline_mll->Draw("E0");

    mc_baseline_mll->SetMarkerColor(kBlack);
    mc_baseline_mll->SetMarkerStyle(kOpenCircle);
    mc_baseline_mll->SetMarkerSize(0.5);
    mc_baseline_mll->SetLineColor(kBlack);
    mc_baseline_mll->Draw("E0,same");

    TH1F *data_after_background_subtraction = (TH1F *)data_baseline_mll->Clone();
    data_after_background_subtraction->Add(data_bkg_template, -1);
    data_after_background_subtraction->SetLineColor(kBlue);
    data_after_background_subtraction->SetMarkerColor(kBlue);
    data_after_background_subtraction->SetMarkerStyle(kFullSquare);
    data_after_background_subtraction->SetMarkerSize(0.5);
    data_after_background_subtraction->Draw("E0,same");
    
    data_bkg_template->SetMarkerColor(kRed);
    data_bkg_template->SetMarkerStyle(kFullTriangleUp);
    data_bkg_template->SetMarkerSize(0.5);
    data_bkg_template->SetLineColor(kRed);
    data_bkg_template->Draw("E0,same");
    
    func_Z_peak->SetLineColor(kGreen);
    func_Z_peak->Draw("same");
    
    // Create a dummy histogram in order to fill color between fitting_range_exclude_low < mll < fitting_range_exclude_up
    double fitting_range_exclude_low;
    double fitting_range_exclude_up;
    if (fitting_range == "range_baseline") {
        fitting_range_exclude_low = 80.;
        fitting_range_exclude_up = 100.;
    }
    else if (fitting_range == "range1") {
        fitting_range_exclude_low = 70.;
        fitting_range_exclude_up = 100.;
    }
    else if (fitting_range == "range2") {
        fitting_range_exclude_low = 75.;
        fitting_range_exclude_up = 100.;
    }
    TH1F *h_fitting_exclude_range = new TH1F("h_80_100", "h_80_100", 1, fitting_range_exclude_low, fitting_range_exclude_up);
    h_fitting_exclude_range->SetFillColor(kYellow);
    h_fitting_exclude_range->SetFillColorAlpha(kYellow, 0.35);
    h_fitting_exclude_range->SetFillStyle(1001); // solid
    h_fitting_exclude_range->SetLineColor(kYellow);
    h_fitting_exclude_range->Fill(90, data_baseline_mll->GetMaximum() * 10);
    h_fitting_exclude_range->Draw("hist,same");
    
    TLegend *leg = new TLegend(0.5, 0.6, 0.9, 0.9);
    leg->AddEntry(data_baseline_mll, "Data baseline", "lp");
    leg->AddEntry(mc_baseline_mll, "Zee MC baseline", "lp");
    leg->AddEntry(data_after_background_subtraction, "Data - bkg template", "lp");
    leg->AddEntry(data_bkg_template, "bkg templage", "lp");
    leg->AddEntry(func_Z_peak, "bkg templage fit", "l");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->Draw("same");
    
    pad2->cd(); // pad2 becomes the current pad
    
    TH1F *frame = pad2->DrawFrame(60, 0.5, 150, 1.49);
    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetNdivisions(405);
    frame->SetLineWidth(1);
    frame->SetXTitle("M_{ll} [GeV]");
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
    
    TH1F *ratio = (TH1F *)data_after_background_subtraction->Clone();
    ratio->Sumw2();
    ratio->Divide(mc_baseline_mll);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerStyle(kFullCircle);
    ratio->SetMarkerSize(0.5);
    ratio->SetLineColor(kBlack);
    ratio->Draw("same");

    TString pdf_filename = "bkg_subtraction_" + template_type;
    if (template_type == "baseline")
        pdf_filename += "_template_";
    else
        pdf_filename += "_";
    pdf_filename += fitting_range;
    if (pt_bin_low == 0 && pt_bin_up == -1 &&
        eta_bin_low == 0 && eta_bin_up == -1) {
        if (truth_match)
            pdf_filename += "_truth";
        // else
        //     pdf_filename += ".pdf";
    }
    else if (eta_bin_low == 0 && eta_bin_up == -1) {
        TH3F *data_baseline_mll_3dim = (TH3F *)data_file->Get("h_baseline_pt_eta_mll");
        TString bin_low_pt = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetXaxis()->GetBinLowEdge(pt_bin_low)));
        TString bin_up_pt = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetXaxis()->GetBinUpEdge(pt_bin_up)));
        TString bin_low_mll = TString::Format("%d", static_cast<int>(mll_window_low));
        TString bin_up_mll = TString::Format("%d", static_cast<int>(mll_window_up));
        if (truth_match)
            pdf_filename = pdf_filename + "_pt" + bin_low_pt + "_" + bin_up_pt + "_truth";
        else
            pdf_filename = pdf_filename + "_pt" + bin_low_pt + "_" + bin_up_pt;
    }
    else if (pt_bin_low == 0 && pt_bin_up == -1) {
        TH3F *data_baseline_mll_3dim = (TH3F *)data_file->Get("h_baseline_pt_eta_mll");
        TString bin_low_eta = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetYaxis()->GetBinLowEdge(eta_bin_low)*100.));
        TString bin_up_eta = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetYaxis()->GetBinUpEdge(eta_bin_up)*100.));
        TString bin_low_mll = TString::Format("%d", static_cast<int>(mll_window_low));
        TString bin_up_mll = TString::Format("%d", static_cast<int>(mll_window_up));
        if (truth_match)
            pdf_filename = pdf_filename + "_eta" + bin_low_eta + "_" + bin_up_eta + "_truth";
        else
            pdf_filename = pdf_filename + "_eta" + bin_low_eta + "_" + bin_up_eta;
    }
    else {
        TH3F *data_baseline_mll_3dim = (TH3F *)data_file->Get("h_baseline_pt_eta_mll");
        TString bin_low_pt = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetXaxis()->GetBinLowEdge(pt_bin_low)));
        TString bin_up_pt = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetXaxis()->GetBinUpEdge(pt_bin_up)));
        TString bin_low_eta = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetYaxis()->GetBinLowEdge(eta_bin_low)*100.));
        TString bin_up_eta = TString::Format("%d", static_cast<int>(data_baseline_mll_3dim->GetYaxis()->GetBinUpEdge(eta_bin_up)*100.));
        TString bin_low_mll = TString::Format("%d", static_cast<int>(mll_window_low));
        TString bin_up_mll = TString::Format("%d", static_cast<int>(mll_window_up));
        if (truth_match)
            pdf_filename = pdf_filename + "_mll" + bin_low_mll + "_" + bin_up_mll + "_pt" + bin_low_pt + "_" + bin_up_pt + "_eta" + bin_low_eta + "_" + bin_up_eta + "_truth";
        else
            pdf_filename = pdf_filename + "_mll" + bin_low_mll + "_" + bin_up_mll + "_pt" + bin_low_pt + "_" + bin_up_pt + "_eta" + bin_low_eta + "_" + bin_up_eta;
    }

    if (trigger == "single_lepton_trigger")
        pdf_filename += ".pdf";
    else if (trigger == "tag_trigger_matched" || trigger == "single_lepton_trigger_tag_trigger_matched")
        pdf_filename += "_tag_trigger_matched.pdf";
    else if (trigger == "dilepton_trigger")
        pdf_filename += "_dilepton_trigger.pdf";
    else if (trigger == "dilepton_trigger_tag_trigger_matched")
        pdf_filename += "_dilepton_trigger_tag_trigger_matched.pdf";

    c1->SaveAs(pdf_filename);
}
