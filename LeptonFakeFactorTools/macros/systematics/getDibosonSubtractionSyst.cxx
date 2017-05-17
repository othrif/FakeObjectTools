#include <TH1.h>
#include "../AtlasStyle.C"

// The main function for taking the max difference between the nominal 
// and the variations, and then saving it to a new histogram as the systematic uncertainty
void getDibosonSubtractionSyst(TString fileNameNom, TString fileNameUp, TString fileNameDown, TString outFileName){

  gROOT->SetBatch();
  gStyle->SetOptTitle(0);

  SetAtlasStyle();

  TFile* fout = TFile::Open(outFileName+".root", "RECREATE");

  TFile* fNom = TFile::Open(fileNameNom,"READ");
  TFile* fUp = TFile::Open(fileNameUp,"READ");
  TFile* fDown = TFile::Open(fileNameDown,"READ");

  std::vector<TString> leps = {"el","mu"};
  for(TString lep : leps){

    TString histName = "FakeFactor_" + lep;

    TH1F* hNom  = (TH1F*)fNom->Get(histName);
    TH1F* hUp   = (TH1F*)fUp->Get(histName);
    TH1F* hDown = (TH1F*)fDown->Get(histName);

    fout->cd();
    TH1F* hout = (TH1F*)hNom->Clone();
    hout->Reset();
    hout->SetName(lep+" diboson subtraction syst");
    hout->SetTitle(lep+" diboson subtraction syst");
    hout->GetXaxis()->SetTitle(lep+" pT [GeV]");
    hout->GetYaxis()->SetTitle("Fake factor diboson subtraction syst");
    hout->GetYaxis()->SetTitleOffset(1.25);

    // Use GetNbins()+2 since GetNbins() doesn't account for under/overflow
    int nBins = hNom->GetXaxis()->GetNbins()+2;
    for(int bin = 0; bin < nBins; ++bin){
      double valNom  = hNom->GetBinContent(bin);
      double valUp   = hUp->GetBinContent(bin);
      double valDown = hDown->GetBinContent(bin);

      double shiftUp   = fabs(valNom - valUp);
      double shiftDown = fabs(valNom - valDown);

      double syst = std::max(shiftUp, shiftDown);

      hout->SetBinContent(bin, syst);
      hout->SetBinError(bin, syst*1e-10);
    }

    TCanvas* can = new TCanvas(lep+" diboson subtraction syst canvas");
    hout->Draw("pE1");
    can->SaveAs("can_"+outFileName+"_"+lep+".pdf");
    can->SaveAs("can_"+outFileName+"_"+lep+".C");
  }

  fout->Write();
  fout->Close();

  return;
}
