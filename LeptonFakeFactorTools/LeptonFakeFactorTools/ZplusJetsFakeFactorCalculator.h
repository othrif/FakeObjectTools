/**
 * @file LeptonFakeFactorTools/ZplusJetsFakeFactorCalculator.h
 * @authors Joseph Reichert <Joseph.Reichert@cern.ch>   Kurt Brendlinger <Kurt.Brendlinger@cern.ch>
 * @date August 2016
 * @brief Fake Factor method:
 * Fake background estimation for EWK 3 lepton SUSY analysis, and originally 
 * used by 13 TeV SMWZ analysis
 *
**/

#ifndef ZplusJetsFakeFactorCalculator_H
#define ZplusJetsFakeFactorCalculator_H

#include "LeptonFakeFactorTools/BaseFakeFactorCalculator.h"

class ZplusJetsFakeFactorCalculator : public BaseFakeFactorCalculator
{

public: 
  /// Standard constructor

  ZplusJetsFakeFactorCalculator(const std::string& name = "ZjetsFakeFactorCalc");
  
  /// Standard destructor
  ~ZplusJetsFakeFactorCalculator();

  // Main methods
public:
  /// Initialize this class
  StatusCode initialize();

  /// Finalize this class; everything that should be done after the event loop should go here
  StatusCode finalize();

  /// The main accept method: the actual cuts are applied here
  bool accept(EventDef& evt);

  // Labeling is Wlep, Zlep1, Zlep2
  void set_subregion(TString s){
    m_subregion = s;
    if     (s=="ttt"){m_do_ttt = true ; m_do_ltt = false; m_do_tlt = false; m_do_ttl = false;}
    else if(s=="ltt"){m_do_ttt = false; m_do_ltt = true ; m_do_tlt = false; m_do_ttl = false;}
    else if(s=="tlt"){m_do_ttt = false; m_do_ltt = false; m_do_tlt = true ; m_do_ttl = false;}
    else if(s=="ttl"){m_do_ttt = false; m_do_ltt = false; m_do_tlt = false; m_do_ttl = true ;}
    else             {std::cerr << "invalid subregion set! Exiting." << std::endl; abort();}
  }

  // Include events which are tlt or ttl events with euu or uee configuration.
  // These are unlikely to be Z+jet, since it should be rare to have real Z+jets
  // where the jet passes the tight criteria while the lepton only passes 
  // the loose criteria (and of course Z->eu is impossible)
  void set_includeUnlikelyZjetEvents(bool b){m_includeUnlikelyZjetEvents = b;};

protected:

  void initAllFakeFactorVsPtHists(void);
  void fillAllFakeFactorVsPtHists(double pt, TString chanFlavor, bool isElectron, bool isNumPlot, double weight, double weightTotErr);
  void fillHistAndAddErrInQuadrature(double pt, TH1F& hist, TH1F& hist_syst, double weight, double weightTotErr);
  
  double addInQuadrature(double val1, double val2);
  double addInQuadrature(double val1, double val2, double val3);

  bool m_includeUnlikelyZjetEvents = false;

  TString m_subregion = "";
  bool m_do_ttt = false;
  bool m_do_ltt = false;
  bool m_do_tlt = false;
  bool m_do_ttl = false;

  TH1F* m_h_allNum_pt = 0;
  TH1F* m_h_elNum_pt = 0;
  TH1F* m_h_muNum_pt = 0;
  TH1F* m_h_eeeNum_pt = 0;
  TH1F* m_h_emmNum_pt = 0;
  TH1F* m_h_meeNum_pt = 0;
  TH1F* m_h_mmmNum_pt = 0;

  TH1F* m_h_allDen_pt = 0;
  TH1F* m_h_elDen_pt = 0;
  TH1F* m_h_muDen_pt = 0;
  TH1F* m_h_eeeDen_pt = 0;
  TH1F* m_h_emmDen_pt = 0;
  TH1F* m_h_meeDen_pt = 0;
  TH1F* m_h_mmmDen_pt = 0;

  TH1F* m_h_allNum_pt_syst = 0;
  TH1F* m_h_elNum_pt_syst = 0;
  TH1F* m_h_muNum_pt_syst = 0;
  TH1F* m_h_eeeNum_pt_syst = 0;
  TH1F* m_h_emmNum_pt_syst = 0;
  TH1F* m_h_meeNum_pt_syst = 0;
  TH1F* m_h_mmmNum_pt_syst = 0;

  TH1F* m_h_allDen_pt_syst = 0;
  TH1F* m_h_elDen_pt_syst = 0;
  TH1F* m_h_muDen_pt_syst = 0;
  TH1F* m_h_eeeDen_pt_syst = 0;
  TH1F* m_h_emmDen_pt_syst = 0;
  TH1F* m_h_meeDen_pt_syst = 0;
  TH1F* m_h_mmmDen_pt_syst = 0;

  int m_num_eee = 0;
  int m_num_eem = 0;
  int m_num_eme = 0;
  int m_num_emm = 0;
  int m_num_mee = 0;
  int m_num_mem = 0;
  int m_num_mme = 0;
  int m_num_mmm = 0;

  std::vector<int> m_cutflow;

};
//----------------------------------------------------------------------------------------
#endif
