/**
 * @file LeptonFakeFactorTools/BaseFakeFactorCalculator.h
 * @authors Joseph Reichert <Joseph.Reichert@cern.ch>   Kurt Brendlinger <Kurt.Brendlinger@cern.ch>
 * @date August 2016
 * @brief Fake Factor method:
 * Fake background estimation for EWK 3 lepton SUSY analysis, and originally 
 * used by 13 TeV SMWZ analysis
 *
**/

#ifndef BaseFakeFactorCalculator_H
#define BaseFakeFactorCalculator_H

#include "LeptonFakeFactorTools/EventDef.h"
#include "LeptonFakeFactorTools/ApplyFakeFactor.h"
#include "AsgTools/MessageCheck.h"
#include "AsgTools/AsgMessaging.h"
#include <TLorentzVector.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TEnv.h>
#include <fstream>
#include <iostream>
#include <string> // for string
#include <vector> // for vector
#include <cmath>

class BaseFakeFactorCalculator : public asg::AsgMessaging
{

public: 
  /// Standard constructor
  BaseFakeFactorCalculator(const std::string& name = "BaseFakeFactorCalculator");
  
  /// Standard destructor
  virtual ~BaseFakeFactorCalculator();

  // Main methods

  /// Initialize this class
  StatusCode initialize();

  /// Finalize this class; everything that should be done after the event loop should go here
  StatusCode finalize();

  /// The main accept method: the actual cuts are applied here
  /// Note that these pure virtual functions *must* be implemented
  /// by any class that inherits from the base class
  virtual bool accept(EventDef& evt) = 0;

  void set_sample(TString s){m_sample = s;};
  void set_configFile(TString s){m_configFile = s;};

protected:

  TString m_name = "";

  // Helper function for default fake factor vs. pt histograms
  TH1F* initFakeFactorVsPtHist(TH1F* hist, TString name);

  // Parse config file
  void configFromFile(bool& property, const std::string& propname, TEnv& rEnv){
    property = rEnv.GetValue(propname.c_str(), (int) property);
  }
  void configFromFile(double& property, const std::string& propname, TEnv& rEnv){
    property = rEnv.GetValue(propname.c_str(), property);
  }
  void configFromFile(int& property, const std::string& propname, TEnv& rEnv){
    property = rEnv.GetValue(propname.c_str(), property);
  }
  void configFromFile(std::string& property, const std::string& propname, TEnv& rEnv){
    property = rEnv.GetValue(propname.c_str(), property.c_str());
  }

  // Config file and its properties
  TString m_configFile = "";

  double m_MtCutLower  = -FLT_MAX;
  double m_MtCutUpper  =  FLT_MAX;
  double m_METCutLower = -FLT_MAX;
  double m_METCutUpper =  FLT_MAX;
  double m_Zwindow     =  FLT_MAX;
  double m_ZwindowVeto =  -1;

  bool m_pairLeptonsUsingMinMt = false;

  double m_ptlllCutUpper = FLT_MAX;

  double m_ptLep1CutLower = -FLT_MAX;
  double m_ptLep2CutLower = -FLT_MAX;
  double m_ptLep3CutLower = -FLT_MAX;

  bool m_reqZeroJets = false;
  bool m_reqGtEq1Jet = false;
  bool m_reqZeroBJets = false;
  bool m_reqGtEq1BJet = false;

  double m_ptJet1CutLower = -FLT_MAX;

  bool m_do_Zjet  = false;
  bool m_do_ttbar = false;

  bool m_do_SR  = false;
  bool m_doApplyFakeFactor = false;

  std::string m_elFakeFactor_histName = "";
  std::string m_muFakeFactor_histName = "";

  double m_ttbar_SF_el = 1;
  double m_ttbar_SF_err_el = 0;
  double m_ttbar_SF_mu = 1;
  double m_ttbar_SF_err_mu= 0;

  // Non-conf properties
  TString m_sample = "";
  TFile* m_outputFile = 0;

  // For applying the fake factors
  ApplyFakeFactor* m_applyFakeFactorToolSignalLep = 0;
  ApplyFakeFactor* m_applyFakeFactorToolOtherLep = 0;
};
//----------------------------------------------------------------------------------------
#endif
