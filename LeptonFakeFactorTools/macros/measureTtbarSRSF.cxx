#include <TH1.h>

// Helper functions
TH1F* numHist(TFile* f, TString chan, TString sampleName, TString lep);
TH1F* addHist(TH1F* outHist, TH1F* inHist);
void printTtbarSF(TH1F* h_data, TH1F* h_ttbar, TH1F* h_other);


// The main function for dividing num by num
void measureTtbarSRSF(){
  TFile* f1 = TFile::Open("all.root", "READ");

  // Idea here is that for the ttbar selection (where DFOS requirement applied):
  // - in the emm channel, one of the muons is definitely the fake
  // - in the mee channel, one of the electrons is definitely the fake
  //
  // Also note in this notation that emm includes emm, mem, and mme for the
  // lepton indexing
  std::vector<TString> leps = {"mee", "emm"};

  // numerator level: only consider 3 tight leptons
  std::vector<TString> listOfChannels = {"ttt"};

  std::vector<TString> dataList = {"data"};
  std::vector<TString> topLike3LBackgrounds = {"ttbar", "tw", "ww"};
  std::vector<TString> other3LBackgrounds = {"wz", "zz", "zjet", "zgam", "tother", "ttv", "singletop", "tz", "vvv", "higgs"};

  for(TString lep : leps){
    std::cout << std::endl;

    TString fakeLep = "NULL";
    if(lep == "emm") fakeLep = "mu";
    else if(lep == "mee") fakeLep = "el";

    for(TString chan : listOfChannels){

      TH1F* h_data_num = 0;
      TH1F* h_ttbar_num = 0;
      TH1F* h_other_num = 0;

      for(TString sample : dataList){
        h_data_num = addHist(h_data_num, numHist(f1, chan, sample, lep));
      }

      for(TString sample : topLike3LBackgrounds){
        h_ttbar_num = addHist(h_ttbar_num, numHist(f1, chan, sample, lep));
      }

      for(TString sample : other3LBackgrounds){
        h_other_num = addHist(h_other_num, numHist(f1, chan, sample, lep));
      }

      std::cout << fakeLep << " " << chan << std::endl; printTtbarSF(h_data_num, h_ttbar_num, h_other_num);
    }

  }

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

void printTtbarSF(TH1F* h_data, TH1F* h_ttbar, TH1F* h_other){
    float n_data  = h_data->Integral(0, h_data->GetNbinsX()+1);
    float n_other = h_other->Integral(0, h_other->GetNbinsX()+1);
    float n_ttbar = h_ttbar->Integral(0, h_ttbar->GetNbinsX()+1);

    float n_num = n_data - n_other;
    float n_den = n_ttbar;

    float ttbarSF = n_num / n_den;

    float err_data = sqrt((h_data->GetSumw2())->GetSum());
    float err_other = sqrt((h_other->GetSumw2())->GetSum());
    float err_ttbar = sqrt((h_ttbar->GetSumw2())->GetSum());

    // If you wanted to print n events per category
    std::cout << "n_data " << n_data << std::endl;
    std::cout << "n_other " << n_other << std::endl;
    std::cout << "n_ttbar " << n_ttbar << std::endl;

    std::cout << "err_data " << err_data << std::endl;
    std::cout << "err_other " << err_other << std::endl;
    std::cout << "err_ttbar " << err_ttbar << std::endl;

    float err_num = sqrt( err_data*err_data + err_other*err_other );

    float totErrTermNum = err_num / n_num;
    float totErrTermDen  = err_ttbar / n_ttbar;

    float ttbarSF_error = fabs(ttbarSF) * sqrt( (totErrTermNum*totErrTermNum) + (totErrTermDen*totErrTermDen) );

    std::cout << " ttbar SF is " << ttbarSF << " \\pm " << ttbarSF_error << std::endl;
}
