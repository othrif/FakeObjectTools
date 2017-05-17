#include <TH1.h>

// Helper functions for the histograms
TH1F* numHist(TFile* f, TString chan, TString sampleName, TString lep);
TH1F* addHist(TH1F* outHist, TH1F* inHist);
TH1F* numSystHist(TFile* f, TString chan, TString sampleName, TString lep);

void measureTtbarSREstimateStatSyst(){
  TFile* f1 = TFile::Open("all.root", "READ");

  std::vector<TString> topLike3LBackgrounds = {"ttbar", "tw", "ww"};
  std::vector<TString> leps = {"eee", "emm", "mee", "mmm"};
  TString chan = "ttt";

  // Four flavor configurations, one channel for three tight leptons
  std::vector<double> estimates = { 0, 0, 0, 0 };
  std::vector<double> stat_error = { 0, 0, 0, 0 };
  std::vector<double> syst_error = { 0, 0, 0, 0 };
  std::vector<double> tot_error = { 0, 0, 0, 0 };

  for(int ilep = 0; ilep < leps.size(); ++ilep){

    TString lep = leps.at(ilep);

    TH1F* h_ttbar_num = 0;
    TH1F* h_ttbar_numSyst = 0;

    for(TString sample : topLike3LBackgrounds){
      h_ttbar_num = addHist(h_ttbar_num, numHist(f1, chan, sample, lep));
      h_ttbar_numSyst = addHist(h_ttbar_numSyst, numSystHist(f1, chan, sample, lep));
    }

    float n_ttbar = h_ttbar_num->Integral(0, h_ttbar_num->GetNbinsX()+1);
    float stat_err_ttbar = sqrt((h_ttbar_num->GetSumw2())->GetSum());
    float syst_err_ttbar = sqrt((h_ttbar_numSyst->GetSumw2())->GetSum());

    estimates[ilep] = n_ttbar;
    stat_error[ilep] = stat_err_ttbar;
    syst_error[ilep] = syst_err_ttbar;
    tot_error[ilep] = sqrt( (stat_err_ttbar*stat_err_ttbar) + (syst_err_ttbar*syst_err_ttbar) );
  }

  double eee_sum = estimates[0];
  double eee_stat_err = stat_error[0];
  double eee_syst_err = syst_error[0];
  double eee_tot_err = tot_error[0];

  double emm_sum = estimates[1];
  double emm_stat_err = stat_error[1];
  double emm_syst_err = syst_error[1];
  double emm_tot_err = tot_error[1];

  double mee_sum = estimates[2];
  double mee_stat_err = stat_error[2];
  double mee_syst_err = syst_error[2];
  double mee_tot_err = tot_error[2];

  double mmm_sum = estimates[3];
  double mmm_stat_err = stat_error[3];
  double mmm_syst_err = syst_error[3];
  double mmm_tot_err = tot_error[3];

  double tot_sum = eee_sum + emm_sum + mee_sum + mmm_sum;
  double tot_stat_err = sqrt(pow(eee_stat_err, 2) + pow(emm_stat_err, 2) + pow(mee_stat_err, 2) + pow(mmm_stat_err, 2));
  double tot_syst_err = sqrt(pow(eee_syst_err, 2) + pow(emm_syst_err, 2) + pow(mee_syst_err, 2) + pow(mmm_syst_err, 2));
  double tot_tot_err = sqrt(pow(eee_tot_err, 2) + pow(emm_tot_err, 2) + pow(mee_tot_err, 2) + pow(mmm_tot_err, 2));

  TString sum_row;
  sum_row.Form(    "Total (stat+syst) & $%2.2f \\pm %2.2f \\pm %2.2f$	& $%2.2f \\pm %2.2f \\pm %2.2f$	& $%2.2f \\pm %2.2f \\pm %2.2f$	&$ %2.2f \\pm %2.2f \\pm %2.2f$	& $%2.2f \\pm %2.2f \\pm %2.2f$", eee_sum, eee_stat_err, eee_syst_err, emm_sum, emm_stat_err, emm_syst_err, mee_sum, mee_stat_err, mee_syst_err, mmm_sum, mmm_stat_err, mmm_syst_err, tot_sum, tot_stat_err, tot_syst_err);
  TString sum_row_tot;
  sum_row_tot.Form("Total (tot)       & $%2.2f \\pm %2.2f 	    $   & $%2.2f \\pm %2.2f 	      $ & $%2.2f \\pm %2.2f 	      $ &$ %2.2f \\pm %2.2f 	      $ & $ %2.2f \\pm %2.2f          $", eee_sum, eee_tot_err, emm_sum, emm_tot_err, mee_sum, mee_tot_err, mmm_sum, mmm_tot_err, tot_sum, tot_tot_err);

  std::cout << "                  &       eee      	        &        emm       	        &        mee       	        &        mmm       	        &        all               " << std::endl;
  std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << sum_row << std::endl;
  std::cout << sum_row_tot << std::endl;

  return;
}

TH1F* numHist(TFile* f, TString chan, TString sampleName, TString lep){
  TH1F* hist = (TH1F*) f->Get(chan+"/"+lep+"Num_"+sampleName+"_pt");

  return hist;
}

TH1F* addHist(TH1F* outHist, TH1F* inHist){
  
  if(!outHist){
    outHist = inHist;
  }
  else{
    outHist->Add(inHist);
  }
  
  return outHist;
}

TH1F* numSystHist(TFile* f, TString chan, TString sampleName, TString lep){
  TH1F* hist = (TH1F*) f->Get(chan+"/"+lep+"Num_"+sampleName+"_pt_syst");

  return hist;
}
