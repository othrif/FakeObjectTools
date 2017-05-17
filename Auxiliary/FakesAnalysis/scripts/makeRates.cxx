#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

using namespace std;

void makeRates(){
  TFile *f = new TFile("/terabig/cmerlass/xAOD_Analysis/LoopSUSY_fakes/Rates/rates_v44/fakeSS_inclusive.root");

  TH1F *h1_el = (TH1F*)f->Get("histoLoose_el0");
  TH1F *h2_el = (TH1F*)f->Get("histoTight_el0");

  TH1F *h3_el = (TH1F*)f->Get("histoLoose_el1");
  TH1F *h4_el = (TH1F*)f->Get("histoTight_el1");

  TH1F *h5_el = (TH1F*)f->Get("all_histoLoose_el");
  TH1F *h6_el = (TH1F*)f->Get("all_histoTight_el");

  TH1F *h1_mu = (TH1F*)f->Get("histoLoose_mu0");
  TH1F *h2_mu = (TH1F*)f->Get("histoTight_mu0");

  TH1F *h3_mu = (TH1F*)f->Get("histoLoose_mu1");
  TH1F *h4_mu = (TH1F*)f->Get("histoTight_mu1");

  TH1F *h5_mu = (TH1F*)f->Get("all_histoLoose_mu");
  TH1F *h6_mu = (TH1F*)f->Get("all_histoTight_mu");
  
 
  TFile *f2 = new TFile("/terabig/cmerlass/xAOD_Analysis/LoopSUSY_fakes/Rates/rates_v44/fakeSS_bkg_inclusive.root");

  TH1F *h1bis_el = (TH1F*)f2->Get("histoLoose_el0");
  TH1F *h2bis_el = (TH1F*)f2->Get("histoTight_el0");
  
  TH1F *h3bis_el = (TH1F*)f2->Get("histoLoose_el1");
  TH1F *h4bis_el = (TH1F*)f2->Get("histoTight_el1");

  TH1F *h5bis_el = (TH1F*)f2->Get("all_histoLoose_el");
  TH1F *h6bis_el = (TH1F*)f2->Get("all_histoTight_el");

  TH1F *h1bis_mu = (TH1F*)f2->Get("histoLoose_mu0");
  TH1F *h2bis_mu = (TH1F*)f2->Get("histoTight_mu0");
  
  TH1F *h3bis_mu = (TH1F*)f2->Get("histoLoose_mu1");
  TH1F *h4bis_mu = (TH1F*)f2->Get("histoTight_mu1");

  TH1F *h5bis_mu = (TH1F*)f2->Get("all_histoLoose_mu");
  TH1F *h6bis_mu = (TH1F*)f2->Get("all_histoTight_mu");


  //electrons
  //nominal ratio
  TH1F *h16_el = (TH1F*)h6_el->Clone("h16_el");
  h16_el->Add(h6bis_el,-1);
  TH1F *h15_el = (TH1F*)h5_el->Clone("h15_el");
  h15_el->Add(h5bis_el,-1);
  h16_el->Divide(h15_el);

  
  //stat UP&DOWN
  TH1F *h16_el_statUP = (TH1F*)h16_el->Clone("h16_el");
  TH1F *h16_el_statDOWN = (TH1F*)h16_el->Clone("h16_el");
  int nb_el = h16_el->GetNbinsX();
  for(int i=1; i<nb_el+1; i++) {
    double error = h16_el->GetBinError(i);
    double value = h16_el->GetBinContent(i);
    h16_el_statUP->SetBinContent(i, value + error);
    h16_el_statDOWN->SetBinContent(i,value - error);
  }


  //realUP
  TH1F *h16_el_realUP = (TH1F*)h6_el->Clone("h16_el");
  h16_el_realUP->Add(h6bis_el,-0.66);
  TH1F *h15_el_realUP = (TH1F*)h5_el->Clone("h15_el");
  h15_el_realUP->Add(h5bis_el,-0.66);
  h16_el_realUP->Divide(h15_el_realUP);

  //realDOWN
  TH1F *h16_el_realDOWN = (TH1F*)h6_el->Clone("h16_el");
  h16_el_realDOWN->Add(h6bis_el,-1.33);
  TH1F *h15_el_realDOWN = (TH1F*)h5_el->Clone("h15_el");
  h15_el_realDOWN->Add(h5bis_el,-1.33);
  h16_el_realDOWN->Divide(h15_el_realDOWN);

  //muons
  //nominal ratio
  TH1F *h16_mu = (TH1F*)h6_mu->Clone("h16_mu");
  h16_mu->Add(h6bis_mu,-1);
  TH1F *h15_mu = (TH1F*)h5_mu->Clone("h15_mu");
  h15_mu->Add(h5bis_mu,-1);
  h16_mu->Divide(h15_mu);


  //stat UP&DOWN
  TH1F *h16_mu_statUP = (TH1F*)h16_mu->Clone("h16_mu");
  TH1F *h16_mu_statDOWN = (TH1F*)h16_mu->Clone("h16_mu");
  int nb_mu = h16_mu->GetNbinsX();
  for(int i=1; i<nb_mu+1; i++) {
    double error = h16_mu->GetBinError(i);
    double value = h16_mu->GetBinContent(i);
    h16_mu_statUP->SetBinContent(i, value + error);
    h16_mu_statDOWN->SetBinContent(i,value - error);
  }


  //realUP
  TH1F *h16_mu_realUP = (TH1F*)h6_mu->Clone("h16_mu");
  h16_mu_realUP->Add(h6bis_mu,-0.66);
  TH1F *h15_mu_realUP = (TH1F*)h5_mu->Clone("h15_mu");
  h15_mu_realUP->Add(h5bis_mu,-0.66);
  h16_mu_realUP->Divide(h15_mu_realUP);

  //realDOWN
  TH1F *h16_mu_realDOWN = (TH1F*)h6_mu->Clone("h16_mu");
  h16_mu_realDOWN->Add(h6bis_mu,-1.33);
  TH1F *h15_mu_realDOWN = (TH1F*)h5_mu->Clone("h15_mu");
  h15_mu_realDOWN->Add(h5bis_mu,-1.33);
  h16_mu_realDOWN->Divide(h15_mu_realDOWN);
  
  //file writing
  string file_name2 = "FakeRatesTables.txt";
  ofstream output_file(file_name2.c_str());

  //electrons
  //param
  output_file << "const double Params_Fake_el_nominal[] = {"<< endl;
  for(int i=1; i<nb_el+1; i++) {
    if(i!=nb_el)
      output_file << h16_el->GetBinContent(i) << ", ";
    else output_file << h16_el->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;
  
  
  output_file << "const double Params_Fake_el_statUP[] = {"<< endl;
  for(int i=1; i<nb_el+1; i++) {
    if(i!=nb_el)
      output_file << h16_el_statUP->GetBinContent(i) << ", ";
    else output_file << h16_el_statUP->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;

  output_file << "const double Params_Fake_el_statDOWN[] = {"<< endl;
  for(int i=1; i<nb_el+1; i++) {
    if(i!=nb_el)
      output_file << h16_el_statDOWN->GetBinContent(i) << ", ";
    else output_file << h16_el_statDOWN->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;

  output_file << "const double Params_Fake_el_realUP[] = {"<< endl;
  for(int i=1; i<nb_el+1; i++) {
    if(i!=nb_el)
      output_file << h16_el_realUP->GetBinContent(i) << ", ";
    else output_file << h16_el_realUP->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;

  output_file << "const double Params_Fake_el_realDOWN[] = {"<< endl;
  for(int i=1; i<nb_el+1; i++) {
    if(i!=nb_el)
      output_file << h16_el_realDOWN->GetBinContent(i) << ", ";
    else output_file << h16_el_realDOWN->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;
 
  // bins
  const double *eta_bins_el = h3_el->GetXaxis()->GetXbins()->GetArray();
  int n_etabins_el = h3_el->GetXaxis()->GetNbins();
  
  output_file << "const float Params_Fake_el_EtaBins[] = {";
  for(int i=0; i<n_etabins_el +1; i++) {
    if(i!=n_etabins_el)
      output_file << eta_bins_el[i] <<", ";
    else  output_file << eta_bins_el[i];
  }
  output_file << "}"<< endl;

  const double *pt_bins_el = h1_el->GetXaxis()->GetXbins()->GetArray();
  int n_ptbins_el = h1_el->GetXaxis()->GetNbins();
  
  output_file << "const float Params_Fake_el_PtBins[] = {";
  for(int i=0; i<n_ptbins_el +1; i++) {
    if(i!=n_ptbins_el)
      output_file << pt_bins_el[i] <<", ";
    else  output_file << pt_bins_el[i];
  }
  output_file << "}"<< endl;
  output_file << ""<< endl;

  //muons
  //param
  output_file << "const double Params_Fake_mu_nominal[] = {"<< endl;
  for(int i=1; i<nb_mu+1; i++) {
    if(i!=nb_mu)
      output_file << h16_mu->GetBinContent(i) << ", ";
    else output_file << h16_mu->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;
  
  
  output_file << "const double Params_Fake_mu_statUP[] = {"<< endl;
  for(int i=1; i<nb_mu+1; i++) {
    if(i!=nb_mu)
      output_file << h16_mu_statUP->GetBinContent(i) << ", ";
    else output_file << h16_mu_statUP->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;

  output_file << "const double Params_Fake_mu_statDOWN[] = {"<< endl;
  for(int i=1; i<nb_mu+1; i++) {
    if(i!=nb_mu)
      output_file << h16_mu_statDOWN->GetBinContent(i) << ", ";
    else output_file << h16_mu_statDOWN->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;

  output_file << "const double Params_Fake_mu_realUP[] = {"<< endl;
  for(int i=1; i<nb_mu+1; i++) {
    if(i!=nb_mu)
      output_file << h16_mu_realUP->GetBinContent(i) << ", ";
    else output_file << h16_mu_realUP->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;

  output_file << "const double Params_Fake_mu_realDOWN[] = {"<< endl;
  for(int i=1; i<nb_mu+1; i++) {
    if(i!=nb_mu)
      output_file << h16_mu_realDOWN->GetBinContent(i) << ", ";
    else output_file << h16_mu_realDOWN->GetBinContent(i) << endl;
  }
  output_file << "};" << endl;
  output_file << " " << endl;
 
  // bins
  const double *eta_bins_mu = h3_mu->GetXaxis()->GetXbins()->GetArray();
  int n_etabins_mu = h3_mu->GetXaxis()->GetNbins();
  
  output_file << "const float Params_Fake_mu_EtaBins[] = {";
  for(int i=0; i<n_etabins_mu +1; i++) {
    if(i!=n_etabins_mu)
      output_file << eta_bins_mu[i] <<", ";
    else  output_file << eta_bins_mu[i];
  }
  output_file << "}"<< endl;

  const double *pt_bins_mu = h1_mu->GetXaxis()->GetXbins()->GetArray();
  int n_ptbins_mu = h1_mu->GetXaxis()->GetNbins();
  
  output_file << "const float Params_Fake_mu_PtBins[] = {";
  for(int i=0; i<n_ptbins_mu +1; i++) {
    if(i!=n_ptbins_mu)
      output_file << pt_bins_mu[i] <<", ";
    else  output_file << pt_bins_mu[i];
  }
  output_file << "}"<< endl;

  output_file.close();









}
