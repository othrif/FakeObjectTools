#ifndef FakesAnalysis_MyxAODAnalysis_H
#define FakesAnalysis_MyxAODAnalysis_H

#include <EventLoop/Algorithm.h>
#include <FakesAnalysis/parametric_histos.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include <TROOT.h>
#include <TH1.h>
#include <TTree.h>

#include "PATInterfaces/SystematicVariation.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"
#include "PATInterfaces/SystematicSet.h"

#include "SUSYTools/SUSYCrossSection.h"
#include "SUSYTools/ISUSYObjDef_xAODTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"

#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODCore/ShallowCopy.h"

#include "PileupReweighting/PileupReweightingTool.h"


struct myAnalysisCollections {
  // Containers
  xAOD::ElectronContainer* _electrons;
  xAOD::PhotonContainer* _photons;
  xAOD::MuonContainer* _muons;
  xAOD::JetContainer* _jets;
};

namespace ST {
class SUSYObjDef_xAOD;
}

using namespace ST;

// GRL
class GoodRunsListSelectionTool;

class MyxAODAnalysis : public EL::Algorithm
{

  SUSY::CrossSectionDB *my_XsecDB;  //!
  SUSYObjDef_xAOD *objTool; //!
  GoodRunsListSelectionTool *m_grl; //!
  CP::PileupReweightingTool* m_Pileup; //!

private:

  std::vector<std::string> trig_list; //!
  std::vector<std::string> trig_list_el; //!
  std::vector<std::string> trig_list_mu; //!
  std::vector<bool> trig_pass; //!
  std::vector<ST::SystInfo> systInfoList; //!

  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  bool m_debug;
  double target_lumi;
  bool doClassification;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  //electrons containers
  // xAOD::ElectronContainer* electrons_nominal; //!
  //xAOD::ShallowAuxContainer* electrons_nominal_aux; //!

  //muons containers
  //xAOD::MuonContainer* muons_nominal; //!
  //xAOD::ShallowAuxContainer* muons_nominal_aux; //!

  //jets containers
  //xAOD::JetContainer* jets_nominal; //!
  //xAOD::ShallowAuxContainer* jets_nominal_aux; //!

  //photons containers
  //xAOD::PhotonContainer* photons_nominal; //!
  // xAOD::ShallowAuxContainer* photons_nominal_aux; //!

  uint64_t nEventsProcessed; //!
  double sumOfWeights; //!
  double sumOfWeightsSquared; //!

  bool isMC; //!
  int whichYear; //!

  double m_xsec; //!
  double m_PUwei; //!
  double m_nPrimVx; //!
  double m_WeightEvents; //!
  double m_eleSF_weight; //!
  double m_muSF_weight; //!
  double m_bTag_weight; //!
  double m_JVT_weight; //!
  double el_trigSF; //!
  double mu_trigSF; //!


  //things that are going to end up in a separate class
  int n_electrons; //!
  int n_muons; //!
  int n_poslep; //!
  int n_neglep; //!
  int n_emuSS; //!
  int n_afterCleaning; //!
  int n_2lep; //!
  int n_emu; //!
  int n_jets; //!
  int n_bjets; //!
  xAOD::Muon* myMuon; //!
  xAOD::Electron* myElectron; //!


  int m_flavProbe; //!
  double m_ptProbe; //!
  double m_etaProbe; //!
  bool m_isTightProbe; //!

  parametric_histos* histo_loose_el; //!
  parametric_histos* histo_tight_el; //!
  parametric_histos* histo_loose_mu; //!
  parametric_histos* histo_tight_mu; //!
  TH1F* h_class; //!
  TH2F* h_classVSpt_ele; //!
  TH2F* h_classVSpt_mu; //!



  xAOD::TEvent *m_event;  //!

  // this is a standard constructor
  MyxAODAnalysis ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

 
  void FindProbe(myAnalysisCollections);
  void Fill_histos();
  void Fill_classification();
  void Register_parametric_histos(parametric_histos*);


  // this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);
};

#endif
