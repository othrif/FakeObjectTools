#include <TH1.h>
#include "../AtlasStyle.C"

// The main function for taking the difference between two variations,
// and then saving it to a new histogram as the systematic uncertainty,
// added in quadrature with the stat err from the MC hists
void getTotalSyst(TString fileNameNom, std::vector<TString> fileNameVars, TString outFileName){

  gROOT->SetBatch();
  gStyle->SetOptTitle(0);
  SetAtlasStyle();

  // make a copy of the input file
  std::ifstream src(fileNameNom, std::ios::binary);
  std::ofstream dst(outFileName+".root", std::ios::binary);
  dst << src.rdbuf();
  src.close();
  dst.close();

  TFile* fNom = TFile::Open(fileNameNom, "READ");
  TFile* fout = TFile::Open(outFileName+".root", "UPDATE");

  std::vector<TString> leps = {"el","mu"};
  for(TString lep : leps){

    TH1F* h_FF = (TH1F*) fNom->Get("FakeFactor_"+lep);
    fout->cd();
    TH1F* h_FF_syst = (TH1F*) h_FF->Clone("FakeFactor_"+lep+"_Syst");
    h_FF_syst->SetTitle("FakeFactor_"+lep+"_Syst");
    h_FF_syst->Reset();
    h_FF_syst->GetYaxis()->SetTitle("Fake factor syst");
    h_FF_syst->GetYaxis()->SetTitleOffset(1.25);

    std::cout << "output hist: " << h_FF_syst->GetName() << std::endl;

    int nBins = h_FF_syst->GetXaxis()->GetNbins()+2;

    for(TString fileNameVar : fileNameVars){

      TFile* fVar = TFile::Open(fileNameVar,"READ");

      TIter next(fVar->GetListOfKeys());
      TKey* key = 0;

      while((key = (TKey*)next())){

        TClass* objClass = gROOT->GetClass(key->GetClassName());
        if(!objClass->InheritsFrom("TH1")) continue;

        TH1F* histVar = (TH1F*) key->ReadObj();
        TString histName = histVar->GetName();
        if( !histName.BeginsWith(lep) ) continue;

        std::cout << "input hist: " << histName << std::endl;

        for(int bin = 0; bin < nBins; ++bin){
          double valVar = histVar->GetBinContent(bin);
          double valNom = h_FF_syst->GetBinContent(bin);

          double syst = sqrt( (valNom*valNom) + (valVar*valVar) );

          //std::cout << "nom: " << valNom << " var: " << valVar << " syst: " << syst << std::endl;

          h_FF_syst->SetBinContent(bin, syst);
          h_FF_syst->SetBinError(bin, syst*1e-10);
        }
      }
    }

    TCanvas* can = new TCanvas(lep+" total syst canvas");
    h_FF_syst->Draw("pE1");
    can->SaveAs("can_"+outFileName+"_"+lep+".pdf");
    can->SaveAs("can_"+outFileName+"_"+lep+".C");
  }

  fout->Write();
  fout->Close();

  return;
}
