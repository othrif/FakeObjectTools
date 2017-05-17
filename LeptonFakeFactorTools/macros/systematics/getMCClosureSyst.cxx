#include <TH1.h>
#include "../AtlasStyle.C"

// The main function for taking the difference between two variations,
// and then saving it to a new histogram as the systematic uncertainty,
// added in quadrature with the stat err from the MC hists
void getMCClosureSyst(TString fileName1, TString fileName2, TString outFileName){

  gROOT->SetBatch();
  gStyle->SetOptTitle(0);

  SetAtlasStyle();

  TFile* fout = TFile::Open(outFileName+".root", "RECREATE");

  // Important! f1 is the nominal, f2 is the variation
  TFile* f1 = TFile::Open(fileName1,"READ");
  TFile* f2 = TFile::Open(fileName2,"READ");

  std::vector<TString> leps = {"el","mu"};
  for(TString lep : leps){

    TString histName = "FakeFactor_" + lep + "_Zjet+Zgam_MC";

    TH1F* h1 = (TH1F*)f1->Get(histName);
    TH1F* h2 = (TH1F*)f2->Get(histName);

    fout->cd();
    TH1F* hout = (TH1F*)h1->Clone();
    hout->Reset();
    hout->SetName(lep+" MC closure syst");
    hout->SetTitle(lep+" MC closure syst");
    hout->GetXaxis()->SetTitle(lep+" pT [GeV]");
    hout->GetYaxis()->SetTitle("Fake factor MC closure syst");
    hout->GetYaxis()->SetTitleOffset(1.25);

    // Use GetNbins()+2 since GetNbins() doesn't account for under/overflow
    int nBins = h1->GetXaxis()->GetNbins()+2;
    for(int bin = 0; bin < nBins; ++bin){
      double val1 = h1->GetBinContent(bin);
      double val2 = h2->GetBinContent(bin);

      double shift = fabs(val1 - val2);

      double statErr1 = h1->GetBinError(bin);
      double statErr2 = h2->GetBinError(bin);

      double statErr = sqrt( (statErr1*statErr1) + (statErr2*statErr2) );

      double syst = sqrt( (shift*shift) + (statErr*statErr) );

      // Syst uncertainty is the sqrt( abs(diff)**2 + statErr**2 ) of the two variations.
      hout->SetBinContent(bin, syst);
      hout->SetBinError(bin, syst*1e-10);
    }

    TCanvas* can = new TCanvas(lep+" MC closure syst canvas");
    hout->Draw("pE1");
    can->SaveAs("can_"+outFileName+"_"+lep+".pdf");
    can->SaveAs("can_"+outFileName+"_"+lep+".C");
  }

  fout->Write();
  fout->Close();

  return;
}
