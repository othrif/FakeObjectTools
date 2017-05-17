#include "LeptonFakeFactorTools/EventDef.h"

//=============================================================================
// Constructor
//=============================================================================
EventDef::EventDef(const std::string& name) :
  asg::AsgMessaging(name)
{
  (void) name;
  return;
}

//=============================================================================
// Destructor
//=============================================================================
EventDef::~EventDef()
{
  return;
}

//=============================================================================
// Helper functions for xAOD inputs
//=============================================================================

void EventDef::fillEventInfo(const xAOD::EventInfo* ei, double totalEventWeight, bool singleLeptonTrigPassed, bool dileptonTrigPassed, bool passedSRCuts){

  EventNumber = ei->eventNumber();
  weight = totalEventWeight;
  passSingleLeptonTrigger = singleLeptonTrigPassed;
  passDileptonTrigger = dileptonTrigPassed;
  isSRLike = passedSRCuts;
}

void EventDef::fillLepton(const xAOD::IParticle* p, bool isBaselineLep, bool isSigLep, bool isOtherLep, bool isAntiIDLep, bool isTrigMatched){

  double pt = p->pt() * GeV;
  double eta = p->eta();
  double phi = p->phi();
  
  int flavor = 0;
  int charge = 0;

  bool isTruth = false;

  const xAOD::Electron* el = dynamic_cast<const xAOD::Electron*>(p);
  const xAOD::Muon*     mu = dynamic_cast<const xAOD::Muon*>(p);

  if(el){
    flavor = 1;
    charge = el->charge();

    // FIXME this isn't the whole story...
    int truthType = el->auxdata<int>("truthType");
    if(truthType == 2){
      isTruth = true;
    }
  }
  else if(mu){
    flavor = 2;
    charge = mu->charge();

    // FIXME this isn't the whole story...
    int truthType = mu->auxdata<int>("truthType");
    if(truthType == 6){
      isTruth = true;
    }
  }
  else{
    std::cerr << "Error! This particle " << p << " is not an electron or a muon!!" << std::endl;
    abort();
  }

  if(isBaselineLep) nBaselineLeptons += 1;

  lepPt.push_back(pt);
  lepEta.push_back(eta);
  lepPhi.push_back(phi);

  lepFlavor.push_back(flavor);
  lepCharge.push_back(charge);

  lepMatchesTrigger.push_back(isTrigMatched);
  lepIsTruth.push_back(isTruth);

  lepPassSigID.push_back(isSigLep);
  lepPassOtherID.push_back(isOtherLep);
  lepPassAntiID.push_back(isAntiIDLep);

}

void EventDef::fillMET(const xAOD::MissingET* MET_obj){
  met_Et  = MET_obj->met() * GeV;
  met_phi = MET_obj->phi();
}

//=============================================================================
// Helper functions for derived quantities
//=============================================================================

// Based on the pennSoftLepton enum:
// (but only e and mu are written to the ntuples,
// and e && mu never occurs)
// 0 = None
// 1 = e
// 2 = mu
// 3 = tau
// 4 = jet
// 5 = photon

bool EventDef::isElectron(int i){
  return lepFlavor.at(i) == 1;
}

bool EventDef::isMuon(int i){
  return lepFlavor.at(i) == 2;
}

bool EventDef::isSFOS(int lep1, int lep2){
  if(lepFlavor.at(lep1) != lepFlavor.at(lep2)) return false;
  if(lepCharge.at(lep1) == lepCharge.at(lep2)) return false;

  return true;
}

bool EventDef::isDFOS(int lep1, int lep2){
  if(lepFlavor.at(lep1) == lepFlavor.at(lep2)) return false;
  if(lepCharge.at(lep1) == lepCharge.at(lep2)) return false;

  return true;
}

TLorentzVector EventDef::getLeptonTLV(int i){

  double pt  = lepPt.at(i);
  double eta = lepEta.at(i);
  double phi = lepPhi.at(i);
  bool isEle = isElectron(i);

  double mass = 0.000511; // electron mass in GeV
  if(!isEle) mass = 0.105658; // muon mass in GeV

  TLorentzVector lep;
  lep.SetPtEtaPhiM(pt, eta, phi, mass);

  return lep;
}

TLorentzVector EventDef::getMETTLV(void){

  double MET = met_Et;
  double phi = met_phi;

  TLorentzVector tlv;
  tlv.SetPtEtaPhiM(MET, 0, phi, MET);

  return tlv;
}

double EventDef::getMll(int lep1, int lep2){
  TLorentzVector tlv_lep1 = getLeptonTLV(lep1);
  TLorentzVector tlv_lep2 = getLeptonTLV(lep2);
  return (tlv_lep1 + tlv_lep2).M();
}

double EventDef::getMT(int lep){
  TLorentzVector tlv_lep = getLeptonTLV(lep);
  TLorentzVector tlv_MET = getMETTLV();

  double pt_l = tlv_lep.Pt();
  double MET  = tlv_MET.Pt();

  double dphi = tlv_MET.DeltaPhi(tlv_lep);
  
  double MTsq = 2*pt_l*MET*( 1-cos(dphi) );
  return sqrt(MTsq);
}
