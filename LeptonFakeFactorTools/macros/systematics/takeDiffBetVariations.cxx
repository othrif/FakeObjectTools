#include <TH1.h>

// FIXME: I think this macro is unused, and can be deleted
// The main function for taking the difference between two variations,
// and then saving it to a new histogram as the systematic uncertainty
void takeDiffBetVariations(TString fileName1, TString fileName2, TString histName1, TString histName2, TString outFileName){

  gROOT->SetBatch();
  gStyle->SetOptTitle(0);

  TFile* fout = TFile::Open(outFileName+".root", "RECREATE");

  TFile* f1 = TFile::Open(fileName1,"READ");
  TFile* f2 = TFile::Open(fileName2,"READ");

  TH1F* h1 = (TH1F*)f1->Get(histName1);
  TH1F* h2 = (TH1F*)f2->Get(histName2);

  fout->cd();
  TH1F* hout = (TH1F*)h1->Clone();
  hout->Reset();
  hout->SetName("Syst uncertainty");
  hout->SetTitle("Syst uncertainty");
  hout->GetYaxis()->SetTitle("Diff bet variations");
  hout->GetYaxis()->SetTitleOffset(1.25);

  // Use GetNbins()+2 since GetNbins() doesn't account for under/overflow
  int nBins = h1->GetXaxis()->GetNbins()+2;
  for(int bin = 0; bin < nBins; ++bin){
    double val1 = h1->GetBinContent(bin);
    double val2 = h2->GetBinContent(bin);

    double statErr1 = h1->GetBinError(bin);
    double statErr2 = h2->GetBinError(bin);

    // Syst uncertainty is the abs(diff) of the two variations
    // Stat uncertainty stored (but maybe never used?) is the average
    // of the stat errors
    hout->SetBinContent(bin, fabs(val1-val2));
    hout->SetBinError(bin, (statErr1+statErr2)/2 );
  }

  TCanvas* can = new TCanvas("syst uncertainty");
  hout->Draw("pE1");
  can->SaveAs("can_"+outFileName+".pdf");

  fout->Write();
  fout->Close();

  return;
}
