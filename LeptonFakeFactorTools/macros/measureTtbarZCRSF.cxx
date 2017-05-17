#include <TH1.h>

// Helper functions
TH1F* denHist(TFile* f, TString chan, TString sampleName, TString lep);
TH1F* addHist(TH1F* outHist, TH1F* inHist);
void printTtbarSF(TH1F* h_data, TH1F* h_ttbar, TH1F* h_other);


// The main function for calculating the SF
void measureTtbarZCRSF(){
  TFile* f1 = TFile::Open("all.root", "READ");

  // Idea here is that for the ttbar selection (where DFOS requirement applied):
  // - in the emm channel, one of the muons is definitely the fake
  // - in the mee channel, one of the electrons is definitely the fake
  //
  // Also note in this notation that emm includes emm, mem, and mme for the
  // lepton indexing
  std::vector<TString> leps = {"mee", "emm"};

  // which lepton in the indexing is the loose one?
  std::vector<TString> listOfChannels = {"ltt", "tlt", "ttl"};

  std::vector<TString> dataList = {"data"};
  std::vector<TString> topLike3LBackgrounds = {"ttbar", "tw", "ww"};
  std::vector<TString> other3LBackgrounds = {"wz", "zz", "zjet", "zgam", "tother", "ttv", "singletop", "tz", "vvv", "higgs"};

  for(TString lep : leps){
    std::cout << std::endl;

    TString fakeLep = "NULL";
    if(lep == "emm") fakeLep = "mu";
    else if(lep == "mee") fakeLep = "el";

    // totals for the given lepFlavor
    TH1F* h_data_den_lep = 0;
    TH1F* h_ttbar_den_lep = 0;
    TH1F* h_other_den_lep = 0;

    for(TString chan : listOfChannels){

      // hist for the sample grouping type
      TH1F* h_data_den = 0;
      TH1F* h_ttbar_den = 0;
      TH1F* h_other_den = 0;

      for(TString sample : dataList){
        //std::cout << "sample is " << sample << std::endl;
        h_data_den = addHist(h_data_den, denHist(f1, chan, sample, lep));
      }

      for(TString sample : topLike3LBackgrounds){
        //std::cout << "sample is " << sample << std::endl;
        h_ttbar_den = addHist(h_ttbar_den, denHist(f1, chan, sample, lep));
      }

      for(TString sample : other3LBackgrounds){
        //std::cout << "sample is " << sample << std::endl;
        h_other_den = addHist(h_other_den, denHist(f1, chan, sample, lep));
      }

      // Include this line if you want the breakdown of ltt, tlt, and ttl for each lepFlavor
      //std::cout << fakeLep << " " << chan; printTtbarSF(h_data_den, h_ttbar_den, h_other_den);

      h_data_den_lep  = addHist(h_data_den_lep,  h_data_den);
      h_ttbar_den_lep = addHist(h_ttbar_den_lep, h_ttbar_den);
      h_other_den_lep = addHist(h_other_den_lep, h_other_den);

    }

    TString chan = "all";
    std::cout << fakeLep << " " << chan << std::endl; printTtbarSF(h_data_den_lep, h_ttbar_den_lep, h_other_den_lep);
  }

  return;
}

TH1F* denHist(TFile* f, TString chan, TString sampleName, TString lep){
  TH1F* hist = (TH1F*) f->Get(chan+"/"+lep+"Den_"+sampleName+"_pt");

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

    // If you wanted to print n events per category
    std::cout << "n_data " << n_data << std::endl;
    std::cout << "n_other " << n_other << std::endl;
    std::cout << "n_ttbar " << n_ttbar << std::endl;

    float n_num = n_data - n_other;
    float n_den = n_ttbar;

    float ttbarSF = n_num / n_den;

    float err_data = sqrt((h_data->GetSumw2())->GetSum());
    float err_other = sqrt((h_other->GetSumw2())->GetSum());
    float err_ttbar = sqrt((h_ttbar->GetSumw2())->GetSum());

    std::cout << "err_data " << err_data << std::endl;
    std::cout << "err_other " << err_other << std::endl;
    std::cout << "err_ttbar " << err_ttbar << std::endl;

    float err_num = sqrt( err_data*err_data + err_other*err_other );

    float totErrTermNum = err_num / n_num;
    float totErrTermDen  = err_ttbar / n_ttbar;

    float ttbarSF_error = fabs(ttbarSF) * sqrt( (totErrTermNum*totErrTermNum) + (totErrTermDen*totErrTermDen) );

    std::cout << " ttbar SF is " << ttbarSF << " \\pm " << ttbarSF_error << std::endl;
}

