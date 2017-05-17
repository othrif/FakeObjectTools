#include <TH1.h>

// Helper functions for histograms
TH1F* numHist(TFile* f, TString chan, TString sampleName, TString lep);
TH1F* addHist(TH1F* outHist, TH1F* inHist);

// The main function for plotting the numerators measured as den*FF
void mcclosure_denTimesMCFF(){
  TFile* f1 = TFile::Open("hadded/all.root", "READ");

  std::vector<TString> zjetList = {"zjet", "zgam"};


  std::vector<TString> leps = {"eee", "emm", "mee", "mmm"};
  std::vector<TString> listOfChannels = {"ltt", "tlt", "ttl"};

  // Four flavor configurations, three channels for one loose + two tight leptons
  std::vector<std::vector<double>> estimates = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
  std::vector<std::vector<double>> est_error = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

  for(int ilep = 0; ilep < leps.size(); ++ilep){

    TString lep = leps.at(ilep);

    for(int ichan = 0; ichan < listOfChannels.size(); ++ichan){

      TH1F* h_zjet_num = 0;

      TString chan = listOfChannels.at(ichan);

      for(TString sample : zjetList){
        h_zjet_num = addHist(h_zjet_num, numHist(f1, chan, sample, lep));
      }

      float n_zjet = h_zjet_num->Integral(0, h_zjet_num->GetNbinsX()+1);

      float err_zjet = sqrt((h_zjet_num->GetSumw2())->GetSum());

      estimates[ilep][ichan] = n_zjet;
      est_error[ilep][ichan] = err_zjet;

    }
  }

  double ltt_sum = estimates[0][0] + estimates[1][0] + estimates[2][0] + estimates[3][0];
  double ltt_err = sqrt(pow(est_error[0][0], 2) + pow(est_error[1][0], 2) + pow(est_error[2][0], 2) + pow(est_error[3][0], 2));

  double tlt_sum = estimates[0][1] + estimates[1][1] + estimates[2][1] + estimates[3][1];
  double tlt_err = sqrt(pow(est_error[0][1], 2) + pow(est_error[1][1], 2) + pow(est_error[2][1], 2) + pow(est_error[3][1], 2));

  double ttl_sum = estimates[0][2] + estimates[1][2] + estimates[2][2] + estimates[3][2];
  double ttl_err = sqrt(pow(est_error[0][2], 2) + pow(est_error[1][2], 2) + pow(est_error[2][2], 2) + pow(est_error[3][2], 2));

  double eee_sum = estimates[0][0] + estimates[0][1] + estimates[0][2];
  double eee_err = sqrt(pow(est_error[0][0], 2) + pow(est_error[0][1], 2) + pow(est_error[0][2], 2));

  double emm_sum = estimates[1][0] + estimates[1][1] + estimates[1][2];
  double emm_err = sqrt(pow(est_error[1][0], 2) + pow(est_error[1][1], 2) + pow(est_error[1][2], 2));

  double mee_sum = estimates[2][0] + estimates[2][1] + estimates[2][2];
  double mee_err = sqrt(pow(est_error[2][0], 2) + pow(est_error[2][1], 2) + pow(est_error[2][2], 2));

  double mmm_sum = estimates[3][0] + estimates[3][1] + estimates[3][2];
  double mmm_err = sqrt(pow(est_error[3][0], 2) + pow(est_error[3][1], 2) + pow(est_error[3][2], 2));

  double tot_sum = ltt_sum + tlt_sum + ttl_sum;
  double tot_err = sqrt(pow(ltt_err, 2) + pow(tlt_err, 2) + pow(ttl_err, 2));

  TString ltt_row;
  ltt_row.Form("N_ltt * F | %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f", estimates[0][0], est_error[0][0], estimates[1][0], est_error[1][0], estimates[2][0], est_error[2][0], estimates[3][0], est_error[3][0], ltt_sum, ltt_err);

  TString tlt_row;
  tlt_row.Form("N_tlt * F | %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f", estimates[0][1], est_error[0][1], estimates[1][1], est_error[1][1], estimates[2][1], est_error[2][1], estimates[3][1], est_error[3][1], tlt_sum, tlt_err);

  TString ttl_row;
  ttl_row.Form("N_ttl * F | %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f", estimates[0][2], est_error[0][2], estimates[1][2], est_error[1][2], estimates[2][2], est_error[2][2], estimates[3][2], est_error[3][2], ttl_sum, ttl_err);

  TString sum_row;
  sum_row.Form("Total     | %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f 	| %2.2f +/- %2.2f", eee_sum, eee_err, emm_sum, emm_err, mee_sum, mee_err, mmm_sum, mmm_err, tot_sum, tot_err);

  std::cout << "          |       eee      	|        emm       	|        mee       	|        mmm       	|        all     " << std::endl;
  std::cout << "-------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << ltt_row << std::endl;
  std::cout << tlt_row << std::endl;
  std::cout << ttl_row << std::endl;
  std::cout << "-------------------------------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << sum_row << std::endl;

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
