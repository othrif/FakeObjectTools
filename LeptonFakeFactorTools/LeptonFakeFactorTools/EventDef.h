#ifndef EventDef_H
#define EventDef_H

#include "AsgTools/MessageCheck.h"
#include "AsgTools/AsgMessaging.h"
#include <xAODEventInfo/EventInfo.h>
#include <xAODEgamma/Electron.h>
#include <xAODMuon/Muon.h>
#include <xAODMissingET/MissingET.h>

#include <TLorentzVector.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <string> // for string
#include <vector> // for vector
#include <cmath>

// Exists to efficiently allow you to pass arguments from
// your ntuple or xAOD into the FF calculation
// Note! All units assumed to be in GeV, not MeV
class EventDef : public asg::AsgMessaging
{

public:
  /// Standard constructor
  EventDef(const std::string& name = "EventDef");

  /// Standard destructor
  ~EventDef();

  //=============================================================================
  // Quantities to set for use as input to the tools
  //=============================================================================

  // Event Number
  unsigned long long EventNumber;

  // Lepton quantities
  int                nBaselineLeptons = 0;
  std::vector<float> lepPt = {};
  std::vector<float> lepEta = {};
  std::vector<float> lepPhi = {};
  std::vector<int>   lepFlavor = {};
  std::vector<int>   lepCharge = {};
  std::vector<int>   lepMatchesTrigger = {};
  std::vector<int>   lepIsTruth = {};

  // MET quantities
  float              met_Et = 0;
  float              met_phi = 0;

  // Jet quantities
  int                nJets = 0;
  int                nBJets = 0;
  float              leadJetPt = 0;

  // ID and isolation operating points satisfied
  std::vector<int>   lepPassSigID = {};
  std::vector<int>   lepPassOtherID = {};
  std::vector<int>   lepPassAntiID = {};

  // Did trigger pass? Need to specify single lepton vs. dilepton for
  // the trigger matching to be done properly
  bool               passSingleLeptonTrigger = false;
  bool               passDileptonTrigger = false;

  // Is the event SR-like? This includes the events which actually satisfy
  // the SR requirement as well as those which may have an antiID lepton
  // in place of a signal lepton
  bool               isSRLike = false;

  // Weights and SFs
  float              weight = 0;
  float              weight_syst_err = 0;

  // For passing the matched indices to the output
  int                lep0_index = -1;
  int                lep1_index = -1;
  int                lep2_index = -1;

  //=============================================================================
  // Helper functions for xAOD inputs
  //=============================================================================

  // Fill EventDef with xAOD::EventInfo, as well as other event-level quantities
  void fillEventInfo(const xAOD::EventInfo* ei, double totalEventWeight, bool singleLeptonTrigPassed, bool dileptonTrigPassed, bool passedSRCuts);

  // Fill EventDef with lepton information. Note that the second function is used in case no "Other" lepton criteria is used
  void fillLepton(const xAOD::IParticle* p, bool isBaselineLep, bool isSigLep, bool isOtherLep, bool isAntiIDLep, bool isTrigMatched);
  void fillLepton(const xAOD::IParticle* p, bool isBaselineLep, bool isSigLep, bool isAntiIDLep, bool isTrigMatched){
    return fillLepton(p, isBaselineLep, isSigLep, isSigLep, isAntiIDLep, isTrigMatched);
  }

  // Fill EventDef with MET information
  void fillMET(const xAOD::MissingET* MET_obj); 

  //=============================================================================
  // Helper functions for derived quantities
  //=============================================================================

  bool isElectron(int i);
  bool isMuon(int i);

  bool isSFOS(int lep1, int lep2);
  bool isDFOS(int lep1, int lep2);

  TLorentzVector getLeptonTLV(int i);
  TLorentzVector getMETTLV(void);

  double getMll(int lep1, int lep2);
  double getMT(int lep);

  const double GeV = 0.001;
};

#endif
