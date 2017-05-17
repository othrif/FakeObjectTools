#include <TH1.h>
#include "AtlasStyle.C"

//std::vector<double> newBinsEl = {0, 10, 15, 20, 30, 50, 100, 200, 400};
//std::vector<double> newBinsMu = {0, 10, 15, 20, 30, 50, 100, 200, 400};

std::vector<double> newBinsEl = {0, 10, 15, 20, 30, 50, 100};
std::vector<double> newBinsMu = {0, 10, 15, 20, 30, 100};

int lnewBinsEl = newBinsEl.size() - 1;
int lnewBinsMu = newBinsMu.size() - 1;

// Helper functions for the histograms
TH1F* getHist(TFile* f, TString sampleName, TString lep, bool isDen);
TH1F* numHist(TFile* f, TString sampleName, TString lep);
TH1F* denHist(TFile* f, TString sampleName, TString lep);


// The main function for dividing num by den
//
// Valid options for leps are:
// {"el"}, {"mu"}, or {"el", "mu"} 
//
// scaleVVByX means if you want to scale WZ, ZZ up by 15%
// that you want scaleVVByX=1.15
void takeRatioForFakeFactor(std::vector<TString> leps, double scaleVVByX=1.0){

  gROOT->SetBatch();
  gStyle->SetOptTitle(0);
  SetAtlasStyle();

  TString suffix = "";
  if(scaleVVByX != 1.0){
    suffix.Form("_ScaleVVBy_%.2f",scaleVVByX);
  }

  TFile* fout = TFile::Open("out"+suffix+".root", "RECREATE");
  TFile* f1 = TFile::Open("all.root", "READ");

  fout->cd();

  for(TString lep : leps){

    TH1F* h_data_den = denHist(f1, "data", lep);
    TH1F* h_data_num = numHist(f1, "data", lep);

    TH1F* h_VV_den = denHist(f1, "wz", lep);
    TH1F* h_VV_num = numHist(f1, "wz", lep);

    TH1F* h_zz_den = denHist(f1, "zz", lep);
    TH1F* h_zz_num = numHist(f1, "zz", lep);

    TH1F* h_zjet_den = denHist(f1, "zjet", lep);
    TH1F* h_zjet_num = numHist(f1, "zjet", lep);

    TH1F* h_zgam_den = denHist(f1, "zgam", lep);
    TH1F* h_zgam_num = numHist(f1, "zgam", lep);

    TH1F* h_data_num_subtracted = (TH1F*) h_data_num->Clone(lep+"Num_data_with_other_bkgs_subtracted_pt");
    TH1F* h_data_den_subtracted = (TH1F*) h_data_den->Clone(lep+"Den_data_with_other_bkgs_subtracted_pt");

    h_VV_den->Add(h_zz_den);
    h_VV_num->Add(h_zz_num);

    // Scale VV up (or down) by this factor
    // (1.15 is up by 15%, 0.85 is down by 15%)
    h_VV_den->Scale(scaleVVByX);
    h_VV_num->Scale(scaleVVByX);

    h_data_den_subtracted->Add(h_VV_den, -1);
    h_data_num_subtracted->Add(h_VV_num, -1);

    h_zjet_den->Add(h_zgam_den);
    h_zjet_num->Add(h_zgam_num);

    std::vector<TString> other3LBackgrounds = {"ttbar", "tw", "ww", "tother", "ttv", "singletop", "tz", "vvv", "higgs"};

    for(TString sample : other3LBackgrounds){
      TH1F* h_other_den = denHist(f1, sample, lep);
      TH1F* h_other_num = numHist(f1, sample, lep);

      h_data_den_subtracted->Add(h_other_den, -1);
      h_data_num_subtracted->Add(h_other_num, -1);
    }

    // Take care of num and den style, plotting, etc.
    TCanvas* den = new TCanvas(lep+" denominator");
    TCanvas* num = new TCanvas(lep+" numerator");

    h_data_den->SetTitle(lep+" data den");
    h_data_den->SetLineColor(kBlack);
    h_data_den->SetMarkerColor(kBlack);
    h_data_den->SetMarkerStyle(20);
    h_data_num->SetTitle(lep+" data num");
    h_data_num->SetLineColor(kBlack);
    h_data_num->SetMarkerColor(kBlack);
    h_data_num->SetMarkerStyle(20);

    h_VV_den->SetTitle(lep+" VV MC den");
    h_VV_den->SetLineColor(kRed+2);
    h_VV_den->SetMarkerColor(kRed+2);
    h_VV_den->SetMarkerStyle(21);
    h_VV_num->SetTitle(lep+" VV MC num");
    h_VV_num->SetLineColor(kRed+2);
    h_VV_num->SetMarkerColor(kRed+2);
    h_VV_num->SetMarkerStyle(21);

    h_data_den_subtracted->SetTitle(lep+" data with other bkgs subtracted den");
    h_data_den_subtracted->SetLineColor(kRed);
    h_data_den_subtracted->SetMarkerColor(kRed);
    h_data_den_subtracted->SetMarkerStyle(25);
    h_data_num_subtracted->SetTitle(lep+" data with other bkgs subtracted num");
    h_data_num_subtracted->SetLineColor(kRed);
    h_data_num_subtracted->SetMarkerColor(kRed);
    h_data_num_subtracted->SetMarkerStyle(25);

    h_zjet_den->SetTitle(lep+" Z+jet and Z+#gamma MC den");
    h_zjet_den->SetLineColor(kBlue);
    h_zjet_den->SetMarkerColor(kBlue);
    h_zjet_den->SetMarkerStyle(24);
    h_zjet_num->SetTitle(lep+" Z+jet and Z+#gamma MC num");
    h_zjet_num->SetLineColor(kBlue);
    h_zjet_num->SetMarkerColor(kBlue);
    h_zjet_num->SetMarkerStyle(24);

    TPaveText* tDen = new TPaveText(0.3, 0.92, 0.7, 0.97, "NDC");
    tDen->SetTextFont(42);
    tDen->SetBorderSize(0);
    tDen->SetFillColor(gStyle->GetTitleFillColor());
    tDen->AddText(lep+" denominator");

    TPaveText* tNum = new TPaveText(0.3, 0.92, 0.7, 0.97, "NDC");
    tNum->SetTextFont(42);
    tNum->SetBorderSize(0);
    tNum->SetFillColor(gStyle->GetTitleFillColor());
    tNum->AddText(lep+" numerator");

    den->cd();
    h_data_den->Draw("pE1");
    h_data_den_subtracted->Draw("pE1 same");
    h_VV_den->Draw("pE1 same");
    h_zjet_den->Draw("pE1 same");
    tDen->Draw("same");

    num->cd();
    h_data_num->Draw("pE1");
    h_data_num_subtracted->Draw("pE1 same");
    h_VV_num->Draw("pE1 same");
    h_zjet_num->Draw("pE1 same");
    tNum->Draw("same");

    den->BuildLegend();
    num->BuildLegend();
    den->SaveAs(lep+"_den"+suffix+".pdf");
    den->SaveAs(lep+"_den"+suffix+".C");
    num->SaveAs(lep+"_num"+suffix+".pdf");
    num->SaveAs(lep+"_num"+suffix+".C");

    fout->cd();

    TH1F* h_data_out = (TH1F*) h_data_num->Clone();
    h_data_out->SetNameTitle("FakeFactor_"+lep+"_No_Bkg_Subtraction","Data");

    TH1F* h_data_out_subtracted = (TH1F*) h_data_num_subtracted->Clone();
    h_data_out_subtracted->SetNameTitle("FakeFactor_"+lep,"Data with other bkgs subtracted");

    TH1F* h_zjet_out = (TH1F*) h_zjet_num->Clone();
    h_zjet_out->SetNameTitle("FakeFactor_"+lep+"_Zjet+Zgam_MC","Z+jet and Z+#gamma MC");

    h_data_out->SetLineColor(kBlack);
    h_data_out->SetMarkerColor(kBlack);
    h_data_out->GetYaxis()->SetTitle("Fake Factor");

    h_data_out_subtracted->SetLineColor(kRed);
    h_data_out_subtracted->SetMarkerColor(kRed);
    h_data_out_subtracted->SetMarkerStyle(25);

    h_zjet_out->SetLineColor(kBlue);
    h_zjet_out->SetMarkerColor(kBlue);
    h_zjet_out->SetMarkerStyle(24);

    h_data_out->SetStats(0);
    h_data_out->Divide(h_data_den);

    h_data_out_subtracted->Divide(h_data_den_subtracted);

    h_zjet_out->Divide(h_zjet_den);

    TCanvas* can = new TCanvas("Fake Factor " + lep);
    h_data_out->Draw("pE1");
    h_data_out_subtracted->Draw("same pE1");
    h_zjet_out->Draw("same pE1");

    if(lep == "el"){
      h_data_out->GetYaxis()->SetRangeUser(0.0,0.20);
      TLegend* leg = can->BuildLegend(0.2, 0.70, 0.48, 0.87);
      leg->SetBorderSize(0);
    }
    else if(lep == "mu"){
      h_data_out->GetYaxis()->SetRangeUser(0.0,0.20);
      TLegend* leg = can->BuildLegend(0.6, 0.35, 0.88, 0.52);
      leg->SetBorderSize(0);
    }

    can->SaveAs("can_FakeFactor_"+lep+suffix+".pdf");
    delete can;
  }

  fout->Write();
  fout->Close();

  return;
}

TH1F* getHist(TFile* f, TString sampleName, TString lep, bool isDen){

  TString folder = "ttt";
  TString numOrDen = "Num";
  TString lepName = "";

  if(lep == "el") lepName = "Electron ";
  else if(lep == "mu") lepName = "Muon ";

  if(isDen){
    folder = "ltt";
    numOrDen = "Den";
  }

  TString name = lep+numOrDen+"_"+sampleName+"_pt";
  TH1F* hist = (TH1F*) f->Get(folder+"/"+name);
  if(!hist){
    std::cerr << "hist " << name << " not found!" << std::endl;
  }
  hist->SetNameTitle(name, name);
  hist->SetStats(0);
  
  // style can be changed later
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(kBlack);
  hist->SetLineColor(kBlack);
  hist->GetXaxis()->SetTitle(lepName+"p_{T} [GeV]");

  std::vector<double> newBins;
  int lnewBins;

  if(lep == "el"){
    newBins = newBinsEl;
    lnewBins = lnewBinsEl;
  }
  else if(lep == "mu"){
    newBins = newBinsMu;
    lnewBins = lnewBinsMu;
  }
  else if(lep == "all"){
    newBins = newBinsEl;
    lnewBins = lnewBinsEl;
  }
  else{
    std::cerr << "lep " << lep << " not defined! Exiting." << std::endl; abort();
  }

  // Rebin to use our fake factor binning. Then we'll move the 
  // underflow and overflows into the first and last bin, respectively
  hist = (TH1F*) hist->Rebin(lnewBins, hist->GetName(), &newBins[0]);

  int underflow_bin = 0;
  int overflow_bin  = hist->GetNbinsX()+1;

  int firstbin = underflow_bin + 1;
  int lastbin  = overflow_bin  - 1;

  // underflow and overflow bin values and errors
  double underflow_val = hist->GetBinContent(underflow_bin);
  double overflow_val  = hist->GetBinContent(overflow_bin);

  double underflow_err = hist->GetBinError(underflow_bin);
  double overflow_err  = hist->GetBinError(overflow_bin);

  // first and last (non-underflow and non-overflow) bins
  double firstbin_val = hist->GetBinContent(firstbin);
  double lastbin_val  = hist->GetBinContent(lastbin);

  double firstbin_err = hist->GetBinError(firstbin);
  double lastbin_err  = hist->GetBinError(lastbin);

  // add flows into first / last bin
  firstbin_val = firstbin_val + underflow_val;
  lastbin_val  = lastbin_val + overflow_val;

  // add errors in quadrature
  firstbin_err = sqrt(firstbin_err*firstbin_err + underflow_err*underflow_err);
  lastbin_err  = sqrt(lastbin_err * lastbin_err + overflow_err * overflow_err);

  // Reset content and error of the bins in the hist
  hist->SetBinContent(firstbin, firstbin_val);
  hist->SetBinError  (firstbin, firstbin_err);

  hist->SetBinContent(lastbin, lastbin_val);
  hist->SetBinError  (lastbin, lastbin_err);

  // Now clear the underflow and overflow
  hist->ClearUnderflowAndOverflow();

  return hist;
}

TH1F* numHist(TFile* f, TString sampleName, TString lep){
  bool isDen = false;
  return getHist(f, sampleName, lep, isDen);
}

TH1F* denHist(TFile* f, TString sampleName, TString lep){
  bool isDen = true;
  return getHist(f, sampleName, lep, isDen);
}
