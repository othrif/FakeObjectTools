#include "LeptonFakeFactorTools/BaseFakeFactorCalculator.h"

//=============================================================================
// Constructor
//=============================================================================
BaseFakeFactorCalculator::BaseFakeFactorCalculator(const std::string& name) :
  asg::AsgMessaging(name)
{

  m_name = name;
  return;
}




//=============================================================================
// Destructor
//=============================================================================
BaseFakeFactorCalculator::~BaseFakeFactorCalculator()
{
  return;
}


StatusCode BaseFakeFactorCalculator::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;

  TEnv rEnv;
  int success = -1;
  success = rEnv.ReadFile(m_configFile.Data(), kEnvAll);
  if( success != 0 ){
    ATH_MSG_FATAL("unable to read conf file " << m_configFile); 
    return StatusCode::FAILURE;
  }

  configFromFile(m_MtCutLower, "MtCutLower", rEnv);
  configFromFile(m_MtCutUpper, "MtCutUpper", rEnv);
  configFromFile(m_METCutLower, "METCutLower", rEnv);
  configFromFile(m_METCutUpper, "METCutUpper", rEnv);
  configFromFile(m_Zwindow, "Zwindow", rEnv);
  configFromFile(m_ZwindowVeto, "ZwindowVeto", rEnv);

  configFromFile(m_pairLeptonsUsingMinMt, "pairLeptonsUsingMinMt", rEnv);

  configFromFile(m_ptlllCutUpper, "ptlllCutUpper", rEnv);

  configFromFile(m_ptLep1CutLower, "ptLep1CutLower", rEnv); 
  configFromFile(m_ptLep2CutLower, "ptLep2CutLower", rEnv); 
  configFromFile(m_ptLep3CutLower, "ptLep3CutLower", rEnv); 

  configFromFile(m_reqZeroJets, "reqZeroJets", rEnv);
  configFromFile(m_reqGtEq1Jet, "reqGtEq1Jet", rEnv);
  configFromFile(m_reqZeroBJets, "reqZeroBJets", rEnv);
  configFromFile(m_reqGtEq1BJet, "reqGtEq1BJet", rEnv);

  configFromFile(m_ptJet1CutLower, "ptJet1CutLower", rEnv);

  configFromFile(m_do_Zjet, "do_Zjet", rEnv);
  configFromFile(m_do_ttbar, "do_ttbar", rEnv);
  configFromFile(m_do_SR, "do_SR", rEnv);
  configFromFile(m_doApplyFakeFactor, "doApplyFakeFactor", rEnv);

  configFromFile(m_elFakeFactor_histName, "elFakeFactor_histName", rEnv);
  configFromFile(m_muFakeFactor_histName, "muFakeFactor_histName", rEnv);

  configFromFile(m_ttbar_SF_el, "ttbar_SF_el", rEnv);
  configFromFile(m_ttbar_SF_err_el, "ttbar_SF_err_el", rEnv);
  configFromFile(m_ttbar_SF_mu, "ttbar_SF_mu", rEnv);
  configFromFile(m_ttbar_SF_err_mu, "ttbar_SF_err_mu", rEnv);

  if(m_sample == ""){
    ATH_MSG_FATAL("Sample is unnamed! Exiting.");
    sc = StatusCode::FAILURE;
    return sc;
  }

  m_outputFile = TFile::Open(m_sample + "_" + m_name + "_out.root", "RECREATE");

  if(m_doApplyFakeFactor){
    m_applyFakeFactorToolSignalLep = new ApplyFakeFactor("ApplyFakeFactor_SignalLep");
    m_applyFakeFactorToolOtherLep = new ApplyFakeFactor("ApplyFakeFactor_OtherLep");

    m_applyFakeFactorToolSignalLep->set_savedFakeFactorFileName("$ROOTCOREBIN/../LeptonFakeFactorTools/data/FakeFactors_April12_2017_2.4.29_IncludingJVTSFs.root");
    m_applyFakeFactorToolOtherLep->set_savedFakeFactorFileName("$ROOTCOREBIN/../LeptonFakeFactorTools/data/FakeFactors_April12_2017_2.4.29_IncludingJVTSFs.root");

    // Set histogram names to read in, for each instance of the tool
    m_applyFakeFactorToolSignalLep->set_saved_elFakeFactorName(m_elFakeFactor_histName);
    m_applyFakeFactorToolSignalLep->set_saved_muFakeFactorName(m_muFakeFactor_histName);
    m_applyFakeFactorToolOtherLep ->set_saved_elFakeFactorName(m_elFakeFactor_histName);
    m_applyFakeFactorToolOtherLep ->set_saved_muFakeFactorName(m_muFakeFactor_histName);

    // Initialize
    m_applyFakeFactorToolSignalLep->initialize().ignore();
    m_applyFakeFactorToolOtherLep->initialize().ignore();
  }

  return sc;
}

StatusCode BaseFakeFactorCalculator::finalize()
{
  // use an int as a StatusCode
  StatusCode sc = StatusCode::SUCCESS;

  if(m_applyFakeFactorToolSignalLep){
    m_applyFakeFactorToolSignalLep->finalize().ignore();
    delete m_applyFakeFactorToolSignalLep;
  }
  if(m_applyFakeFactorToolOtherLep){
    m_applyFakeFactorToolOtherLep->finalize().ignore();
    delete m_applyFakeFactorToolOtherLep;
  }

  // Note! The Close() deletes from memory any hists associated with the file as well
  m_outputFile->Write();
  m_outputFile->Close(); 

  return sc;
}

TH1F* BaseFakeFactorCalculator::initFakeFactorVsPtHist(TH1F* hist, TString name){

  hist = new TH1F(name, name, 400, 0, 400);
  hist->Sumw2();
  return hist;
}
