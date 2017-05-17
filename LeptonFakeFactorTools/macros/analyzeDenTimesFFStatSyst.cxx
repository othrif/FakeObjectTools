#include <TH1.h>

// Helper functions for histograms
TH1F* numHist(TFile* f, TString chan, TString sampleName, TString lep);
TH1F* addHist(TH1F* outHist, TH1F* inHist);
TH1F* numSystHist(TFile* f, TString chan, TString sampleName, TString lep);

// The main function for plotting the numerators measured as den*FF
void analyzeDenTimesFFStatSyst(bool verbose=false){
  TFile* f1 = TFile::Open("all.root", "READ");

  std::vector<TString> dataList = {"data"};
  std::vector<TString> topLike3LBackgrounds = {"ttbar", "tw", "ww"};
  std::vector<TString> otherBackgrounds = {"wz", "zz", "tother", "ttv", "singletop", "tz", "vvv", "higgs"};

  std::vector<TString> leps = {"eee", "emm", "mee", "mmm"};
  std::vector<TString> listOfChannels = {"ltt", "tlt", "ttl"};

  // Four flavor configurations, three channels for one loose + two tight leptons
  std::vector<std::vector<double>> estimates  = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
  std::vector<std::vector<double>> stat_error = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
  std::vector<std::vector<double>> syst_error = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
  std::vector<std::vector<double>> tot_error  = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

  for(int ilep = 0; ilep < leps.size(); ++ilep){

    TString lep = leps.at(ilep);

    for(int ichan = 0; ichan < listOfChannels.size(); ++ichan){

      TH1F* h_data_num = 0;
      TH1F* h_ttbar_num = 0;
      TH1F* h_other_num = 0;

      TH1F* h_data_numSyst = 0;
      TH1F* h_ttbar_numSyst = 0;
      TH1F* h_other_numSyst = 0;

      TString chan = listOfChannels.at(ichan);

      for(TString sample : dataList){
        h_data_num = addHist(h_data_num, numHist(f1, chan, sample, lep));
        h_data_numSyst = addHist(h_data_numSyst, numSystHist(f1, chan, sample, lep));
      }
      for(TString sample : topLike3LBackgrounds){
        h_ttbar_num = addHist(h_ttbar_num, numHist(f1, chan, sample, lep));
        h_ttbar_numSyst = addHist(h_ttbar_numSyst, numSystHist(f1, chan, sample, lep));
      }
      for(TString sample : otherBackgrounds){
        h_other_num = addHist(h_other_num, numHist(f1, chan, sample, lep));
        h_other_numSyst = addHist(h_other_numSyst, numSystHist(f1, chan, sample, lep));
      }

      float n_data  = h_data_num->Integral(0, h_data_num->GetNbinsX()+1);
      float n_other = h_other_num->Integral(0, h_other_num->GetNbinsX()+1);
      float n_ttbar = h_ttbar_num->Integral(0, h_ttbar_num->GetNbinsX()+1);

      if(verbose){
        std::cout << "ilep is " << ilep << " ichan is " << ichan << std::endl;
        std::cout << "n_data " << n_data << std::endl;
        std::cout << "n_other " << n_other << std::endl;
        std::cout << "n_ttbar " << n_ttbar << std::endl;
      }

      float N_ZSR = n_data - n_other - n_ttbar;

      float stat_err_data = sqrt((h_data_num->GetSumw2())->GetSum());
      float stat_err_other = sqrt((h_other_num->GetSumw2())->GetSum());
      float stat_err_ttbar = sqrt((h_ttbar_num->GetSumw2())->GetSum());

      float syst_err_data = sqrt((h_data_numSyst->GetSumw2())->GetSum());
      float syst_err_other = sqrt((h_other_numSyst->GetSumw2())->GetSum());
      float syst_err_ttbar = sqrt((h_ttbar_numSyst->GetSumw2())->GetSum());

      float N_ZSR_stat_error = sqrt( (stat_err_data * stat_err_data)
                              + (stat_err_other * stat_err_other)
                              + (stat_err_ttbar * stat_err_ttbar)
                              );

      float N_ZSR_syst_error = sqrt( (syst_err_data * syst_err_data)
                              + (syst_err_other * syst_err_other)
                              + (syst_err_ttbar * syst_err_ttbar)
                              );

      float N_ZSR_tot_error = sqrt( (N_ZSR_stat_error*N_ZSR_stat_error) + (N_ZSR_syst_error*N_ZSR_syst_error) );

      estimates [ilep][ichan] = N_ZSR;
      stat_error[ilep][ichan] = N_ZSR_stat_error;
      syst_error[ilep][ichan] = N_ZSR_syst_error;
      tot_error [ilep][ichan] = N_ZSR_tot_error;

    }
  }

  double ltt_sum = estimates[0][0] + estimates[1][0] + estimates[2][0] + estimates[3][0];
  double ltt_stat_err = sqrt(pow(stat_error[0][0], 2) + pow(stat_error[1][0], 2) + pow(stat_error[2][0], 2) + pow(stat_error[3][0], 2));
  double ltt_syst_err = sqrt(pow(syst_error[0][0], 2) + pow(syst_error[1][0], 2) + pow(syst_error[2][0], 2) + pow(syst_error[3][0], 2));
  double ltt_tot_err = sqrt(pow(tot_error[0][0], 2) + pow(tot_error[1][0], 2) + pow(tot_error[2][0], 2) + pow(tot_error[3][0], 2));

  double tlt_sum = estimates[0][1] + estimates[1][1] + estimates[2][1] + estimates[3][1];
  double tlt_stat_err = sqrt(pow(stat_error[0][1], 2) + pow(stat_error[1][1], 2) + pow(stat_error[2][1], 2) + pow(stat_error[3][1], 2));
  double tlt_syst_err = sqrt(pow(syst_error[0][1], 2) + pow(syst_error[1][1], 2) + pow(syst_error[2][1], 2) + pow(syst_error[3][1], 2));
  double tlt_tot_err = sqrt(pow(tot_error[0][1], 2) + pow(tot_error[1][1], 2) + pow(tot_error[2][1], 2) + pow(tot_error[3][1], 2));

  double ttl_sum = estimates[0][2] + estimates[1][2] + estimates[2][2] + estimates[3][2];
  double ttl_stat_err = sqrt(pow(stat_error[0][2], 2) + pow(stat_error[1][2], 2) + pow(stat_error[2][2], 2) + pow(stat_error[3][2], 2));
  double ttl_syst_err = sqrt(pow(syst_error[0][2], 2) + pow(syst_error[1][2], 2) + pow(syst_error[2][2], 2) + pow(syst_error[3][2], 2));
  double ttl_tot_err = sqrt(pow(tot_error[0][2], 2) + pow(tot_error[1][2], 2) + pow(tot_error[2][2], 2) + pow(tot_error[3][2], 2));

  double eee_sum = estimates[0][0] + estimates[0][1] + estimates[0][2];
  double eee_stat_err = sqrt(pow(stat_error[0][0], 2) + pow(stat_error[0][1], 2) + pow(stat_error[0][2], 2));
  double eee_syst_err = sqrt(pow(syst_error[0][0], 2) + pow(syst_error[0][1], 2) + pow(syst_error[0][2], 2));
  double eee_tot_err = sqrt(pow(tot_error[0][0], 2) + pow(tot_error[0][1], 2) + pow(tot_error[0][2], 2));

  double emm_sum = estimates[1][0] + estimates[1][1] + estimates[1][2];
  double emm_stat_err = sqrt(pow(stat_error[1][0], 2) + pow(stat_error[1][1], 2) + pow(stat_error[1][2], 2));
  double emm_syst_err = sqrt(pow(syst_error[1][0], 2) + pow(syst_error[1][1], 2) + pow(syst_error[1][2], 2));
  double emm_tot_err = sqrt(pow(tot_error[1][0], 2) + pow(tot_error[1][1], 2) + pow(tot_error[1][2], 2));

  double mee_sum = estimates[2][0] + estimates[2][1] + estimates[2][2];
  double mee_stat_err = sqrt(pow(stat_error[2][0], 2) + pow(stat_error[2][1], 2) + pow(stat_error[2][2], 2));
  double mee_syst_err = sqrt(pow(syst_error[2][0], 2) + pow(syst_error[2][1], 2) + pow(syst_error[2][2], 2));
  double mee_tot_err = sqrt(pow(tot_error[2][0], 2) + pow(tot_error[2][1], 2) + pow(tot_error[2][2], 2));

  double mmm_sum = estimates[3][0] + estimates[3][1] + estimates[3][2];
  double mmm_stat_err = sqrt(pow(stat_error[3][0], 2) + pow(stat_error[3][1], 2) + pow(stat_error[3][2], 2));
  double mmm_syst_err = sqrt(pow(syst_error[3][0], 2) + pow(syst_error[3][1], 2) + pow(syst_error[3][2], 2));
  double mmm_tot_err = sqrt(pow(tot_error[3][0], 2) + pow(tot_error[3][1], 2) + pow(tot_error[3][2], 2));

  double tot_sum = ltt_sum + tlt_sum + ttl_sum;
  double tot_stat_err = sqrt(pow(ltt_stat_err, 2) + pow(tlt_stat_err, 2) + pow(ttl_stat_err, 2));
  double tot_syst_err = sqrt(pow(ltt_syst_err, 2) + pow(tlt_syst_err, 2) + pow(ttl_syst_err, 2));
  double tot_tot_err = sqrt(pow(ltt_tot_err, 2) + pow(tlt_tot_err, 2) + pow(ttl_tot_err, 2));

  TString ltt_row;
  ltt_row.Form(    "$N_\\text{LTT} \\cdot F$ & $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$ ", estimates[0][0], stat_error[0][0], syst_error[0][0],estimates[1][0], stat_error[1][0], syst_error[1][0], estimates[2][0], stat_error[2][0], syst_error[2][0], estimates[3][0], stat_error[3][0], syst_error[3][0], ltt_sum, ltt_stat_err, ltt_syst_err);

  TString tlt_row;
  tlt_row.Form(    "$N_\\text{TLT} \\cdot F$ & $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$ ", estimates[0][1], stat_error[0][1], syst_error[0][1], estimates[1][1], stat_error[1][1], syst_error[1][1], estimates[2][1], stat_error[2][1], syst_error[2][1], estimates[3][1], stat_error[3][1], syst_error[3][1], tlt_sum, tlt_stat_err, tlt_syst_err);

  TString ttl_row;
  ttl_row.Form(    "$N_\\text{TTL} \\cdot F$ & $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$ ", estimates[0][2], stat_error[0][2], syst_error[0][2], estimates[1][2], stat_error[1][2], syst_error[1][2], estimates[2][2], stat_error[2][2], syst_error[2][2], estimates[3][2], stat_error[3][2], syst_error[3][2], ttl_sum, ttl_stat_err, ttl_syst_err);

  TString sum_row;
  sum_row.Form(    "Total (stat+syst)        & $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$  	& $ %2.2f \\pm %2.2f \\pm %2.2f$ ", eee_sum, eee_stat_err, eee_syst_err, emm_sum, emm_stat_err, emm_syst_err, mee_sum, mee_stat_err, mee_syst_err, mmm_sum, mmm_stat_err, mmm_syst_err, tot_sum, tot_stat_err, tot_syst_err);

  TString sum_row_tot;
  sum_row_tot.Form("Total (tot)              & $ %2.2f \\pm %2.2f         $  	& $ %2.2f \\pm %2.2f         $  	& $ %2.2f \\pm %2.2f         $  	& $ %2.2f \\pm %2.2f         $  	& $ %2.2f \\pm %2.2f         $ ", eee_sum, eee_tot_err, emm_sum, emm_tot_err, mee_sum, mee_tot_err, mmm_sum, mmm_tot_err, tot_sum, tot_tot_err);

  std::cout << "                  &       eee      	        &        emm       	        &        mee       	        &        mmm       	        &        all                       " << std::endl;
  std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  if(verbose){
    std::cout << ltt_row << std::endl;
    std::cout << tlt_row << std::endl;
    std::cout << ttl_row << std::endl;
    std::cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
  }
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
