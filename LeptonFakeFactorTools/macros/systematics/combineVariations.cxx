#include <TH1.h>

double computeAverage(std::vector<double> vals);
double computeSumInQuadrature(std::vector<double> vals);
double computeRMS(std::vector<double> vals); // not actually used right now!

// FIXME: I think this macro is unused, and can be deleted
// The main function for taking the systematic uncertainties from
// individual variations and combining it into a total systematic
void combineVariations(std::vector<TString> fileNames, TString outFileName){

  gROOT->SetBatch();
  gStyle->SetOptTitle(0);

  TFile* fout = TFile::Open(outFileName+".root", "RECREATE");
  fout->cd();
  TH1F* hout = (TH1F*) TFile::Open(fileNames[0],"READ")->Get("Syst uncertainty")->Clone();
  hout->Reset();
  hout->SetName("Total syst uncertainty");
  hout->SetTitle("Total syst uncertainty");
  hout->GetYaxis()->SetTitle("Total syst uncertainty");
  hout->GetYaxis()->SetTitleOffset(1.25);

  std::vector<TH1F*> hists = {};

  for(TString fileName : fileNames){

    TFile* f = TFile::Open(fileName,"READ");
    TH1F* h = (TH1F*) f->Get("Syst uncertainty");

    hists.push_back(h);
  }

  // Use GetNbins()+2 since GetNbins() doesn't account for under/overflow
  int nBins = hout->GetXaxis()->GetNbins()+2;
  for(int bin = 0; bin < nBins; ++bin){

    std::vector<double> vals = {};
    std::vector<double> statErrs = {};

    for(TH1F* h : hists){
      vals.push_back(h->GetBinContent(bin));
      statErrs.push_back(h->GetBinError(bin));
    }

    double sumInQuadrature = computeSumInQuadrature(vals);
    double avgStatErr = computeAverage(statErrs);

    hout->SetBinContent(bin, sumInQuadrature);
    hout->SetBinError(bin, avgStatErr);
  }

  TCanvas* can = new TCanvas("Total syst uncertainty");
  hout->Draw("pE1");
  can->SaveAs("can_"+outFileName+".pdf");

  fout->Write();
  fout->Close();

  return;
}

double computeAverage(std::vector<double> vals){
  double sum = 0;
  for(double val : vals){
    sum += val;
  }
  
  return sum / (double)vals.size();
}
double computeSumInQuadrature(std::vector<double> vals){
  double squaredSum = 0;
  for(double val : vals){
    squaredSum += (val*val);
  }

  return sqrt(squaredSum);
}
double computeRMS(std::vector<double> vals){
  double avg = computeAverage(vals);
  double numerator = 0;

  for(double val : vals){
    numerator += ((avg-val)*(avg-val));
  }
  return sqrt(numerator / (double)vals.size());
}
