#ifndef ytRealLeptonsEfficiency_H
#define ytRealLeptonsEfficiency_H

#include <EventLoop/Algorithm.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>

#include <iostream>
#include <vector>
using namespace std;

class ytRealLeptonsEfficiency : public EL::Algorithm
{
    // put your configuration variables here as public variables.
    // that way they can be set directly from CINT and python.
public:
    // float cutValue;
    bool isMC;
    bool isData;

    TString trigger;
    TString lepton;
    TString process;

    bool tag_trigger_matched;

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
public:
    // Tree *myTree; //!
    // TH1 *myHist; //!

    TTree          *fChain;   //!pointer to the analyzed TTree or TChain

    // Declaration of leaf types
    Bool_t          HLT_e24_lhmedium_nod0_ivarloose; //!
    Bool_t          HLT_e24_lhtight_nod0_ivarloose; //!
    Bool_t          HLT_e24_lhmedium_nod0_L1EM20VH; //!
    Bool_t          HLT_e26_lhtight_iloose; //!
    Bool_t          HLT_e26_lhtight_ivarloose; //!
    Bool_t          HLT_e26_lhtight_nod0_iloose; //!
    Bool_t          HLT_e26_lhtight_nod0_ivarloose; //!
    Bool_t          HLT_e60_lhmedium; //!
    Bool_t          HLT_e60_lhmedium_nod0; //!
    Bool_t          HLT_e120_lhloose_nod0; //!
    Bool_t          HLT_e140_lhloose_nod0; //!
    Bool_t          HLT_2e17_lhvloose_nod0; //!
    Bool_t          HLT_2e15_lhvloose_nod0_L12EM13VH; //!
    Bool_t          HLT_e24_lhmedium_e9_lhmedium; //!
    Bool_t          HLT_e24_lhmedium_L1EM20VH; //!
    Bool_t          HLT_e24_lhmedium_iloose_L1EM20VH; //!
    Bool_t          HLT_e12_lhvloose_L12EM10VH; //!
    Bool_t          HLT_e17_lhloose_2e9_lhloose; //!
    Bool_t          HLT_mu24_ivarmedium; //!
    Bool_t          HLT_mu24_imedium; //!
    Bool_t          HLT_mu24_ivarloose; //!
    Bool_t          HLT_mu24_iloose; //!
    Bool_t          HLT_mu26_ivarmedium; //!
    Bool_t          HLT_mu20_ivarmedium_L1MU15; //!
    Bool_t          HLT_mu20_imedium_L1MU15; //!
    Bool_t          HLT_mu20_ivarloose_L1MU15; //!
    Bool_t          HLT_mu20_iloose_L1MU15; //!
    Bool_t          HLT_mu20_L1MU15; //!
    Bool_t          HLT_mu20_mu8noL1; //!
    Bool_t          HLT_mu22_mu8noL1; //!
    Bool_t          HLT_mu20_2mu4noL1; //!
    Bool_t          HLT_mu22_2mu4noL1; //!
    Bool_t          HLT_mu40; //!
    Bool_t          HLT_mu50; //!
    Bool_t          HLT_2mu10; //!
    Bool_t          HLT_2mu10_nomucomb; //!
    Bool_t          HLT_2mu14; //!
    Bool_t          HLT_2mu14_nomucomb; //!
    Bool_t          HLT_3mu6; //!
    Bool_t          HLT_3mu6_msonly; //!
    Bool_t          HLT_xe100_L1XE50; //!
    Bool_t          HLT_xe80_mht_L1XE50; //!
    Bool_t          HLT_xe90_mht_L1XE50; //!
    Bool_t          HLT_xe100_mht_L1XE50; //!
    Bool_t          HLT_xe100_tc_em_L1XE50; //!
    Bool_t          HLT_xe110_mht_L1XE50; //!
    Bool_t          HLT_xe80_tc_lcw_L1XE50; //!
    Bool_t          HLT_xe90_tc_lcw_L1XE50; //!
    Bool_t          HLT_xe80_tc_lcw_wEFMu_L1XE50; //!
    Bool_t          HLT_e7_lhmedium_mu24; //!
    Bool_t          HLT_e7_lhmedium_nod0_mu24; //!
    Bool_t          HLT_e17_lhloose_nod0_mu14; //!
    Bool_t          HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1; //!
    Bool_t          HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1; //!
    Bool_t          HLT_2e12_lhloose_L12EM10VH; //!
    Bool_t          HLT_e17_lhloose_mu14; //!
    Bool_t          HLT_mu18_mu8noL1; //!
    Bool_t          HLT_xe70; //!
    ULong64_t       EventNumber; //!
    Int_t           ChannelNumber; //!
    Double_t        AvgMu; //!
    Double_t        EventWeight; //!
    Double_t        PRWWeight; //!
    Int_t           PRWrandomRunNumber; //!
    Double_t        TriggerDileptonSF; //!
    Int_t           bcid; //!
    Int_t           LB; //!
    Int_t           RunNb; //!
    Int_t           DetError; //!
    Int_t           NMu; //!
    vector<double>  *Mu_eta; //!
    vector<double>  *Mu_phi; //!
    vector<double>  *Mu_pT; //!
    vector<double>  *Mu_SFw; //!
    vector<float>   *Mu_IsoSFw; //!
    vector<int>     *Mu_charge; //!
    vector<double>  *Mu_d0pvtx; //!
    vector<double>  *Mu_sigd0; //!
    vector<double>  *Mu_z0sinTheta; //!
    vector<bool>    *Mu_isBad; //!
    vector<bool>    *Mu_passOR; //!
    vector<bool>    *Mu_isTight; //!
    vector<bool>    *Mu_isCosmic; //!
    vector<int>     *Mu_type; //!
    vector<int>     *Mu_origin; //!
    vector<double>  *Mu_ptcone20; //!
    vector<double>  *Mu_ptcone30; //!
    vector<double>  *Mu_ptcone40; //!
    vector<double>  *Mu_ptvarcone20; //!
    vector<double>  *Mu_ptvarcone30; //!
    vector<double>  *Mu_ptvarcone40; //!
    vector<double>  *Mu_topoetcone20; //!
    vector<double>  *Mu_topoetcone30; //!
    vector<double>  *Mu_topoetcone40; //!
    vector<bool>    *Mu_passIsoLooseTO; //!
    vector<bool>    *Mu_passIsoLoose; //!
    vector<bool>    *Mu_passIsoTight; //!
    vector<bool>    *Mu_passIsoGrad; //!
    vector<bool>    *Mu_passIsoGradCustomTight; //!
    vector<bool>    *Mu_passIsoGradCustom; //!
    vector<bool>    *Mu_passIsoGradLoose; //!
    Float_t         MuTrigSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50; //!
    vector<bool>    *Mu_trigMatch_mu26_ivarmedium; //!
    vector<bool>    *Mu_trigMatch_mu20_iloose_L1MU15; //!
    vector<bool>    *Mu_trigMatch_mu40; //!
    vector<bool>    *Mu_trigMatch_mu50; //!
    vector<bool>    *Mu_trigMatch_mu8noL1; //!
    vector<bool>    *Mu_trigMatch_mu14; //!
    vector<bool>    *Mu_trigMatch_mu18; //!
    vector<bool>    *Mu_trigMatch_mu18_mu8noL1; //!
    vector<bool>    *Mu_trigMatch_e17_lhloose_mu14; //!
    vector<bool>    *Mu_trigMatch_e17_lhloose_nod0_mu14; //!
    vector<bool>    *Mu_trigMatch_mu20_mu8noL1; //!
    vector<bool>    *Mu_trigMatch_mu22_mu8noL1; //!
    vector<bool>    *Mu_trigMatch_mu24_iloose; //!
    vector<bool>    *Mu_trigMatch_mu24_ivarloose; //!
    vector<bool>    *Mu_trigMatch_mu24_iloose_L1MU15; //!
    vector<bool>    *Mu_trigMatch_mu24_ivarloose_L1MU15; //!
    vector<vector<bool> > *Mu_trigMatchPair_mu18_mu8noL1; //!
    vector<vector<bool> > *Mu_trigMatchPair_mu20_mu8noL1; //!
    vector<vector<bool> > *Mu_trigMatchPair_mu22_mu8noL1; //!
    Int_t           NEl; //!
    vector<double>  *El_eta; //!
    vector<double>  *El_etaclus; //!
    vector<double>  *El_phi; //!
    vector<double>  *El_pT; //!
    vector<double>  *El_E; //!
    vector<int>     *El_charge; //!
    vector<double>  *El_sigd0; //!
    vector<double>  *El_z0sinTheta; //!
    vector<double>  *El_d0pvtx; //!
    vector<bool>    *El_passOR; //!
    vector<double>  *El_SFwMediumLH; //!
    vector<float>   *El_IsoSFwMediumLH; //!
    vector<double>  *El_SFwTightLH; //!
    vector<double>  *El_SFwLooseAndBLayerLH; //!
    vector<double>  *El_SFwTrigMediumLH_e12_lhloose_L1EM10VH; //!
    vector<double>  *El_SFwTrigMediumLH_e17_lhloose; //!
    vector<double>  *El_SFwTrigMediumLH_single; //!
    vector<double>  *El_SFwTrigLooseAndBLayerLH_e12_lhloose_L1EM10VH; //!
    vector<bool>    *El_isLooseAndBLayerLH; //!
    vector<bool>    *El_isMediumLH; //!
    vector<bool>    *El_isTightLH; //!
    vector<int>     *El_nBLayerHits; //!
    vector<int>     *El_expectBLayerHit; //!
    vector<int>     *El_truthType; //!
    vector<int>     *El_truthOrigin; //!
    vector<int>     *El_bkgTruthType; //!
    vector<int>     *El_bkgTruthOrigin; //!
    vector<int>     *El_bkgMotherPdgId; //!
    vector<int>     *El_firstEgMotherTruthType; //!
    vector<int>     *El_firstEgMotherTruthOrigin; //!
    vector<int>     *El_firstEgMotherPdgId; //!
    vector<int>     *El_chFlip; //!
    vector<double>  *El_ptcone20; //!
    vector<double>  *El_ptcone30; //!
    vector<double>  *El_ptcone40; //!
    vector<double>  *El_ptvarcone20; //!
    vector<double>  *El_ptvarcone30; //!
    vector<double>  *El_ptvarcone40; //!
    vector<double>  *El_topoetcone20; //!
    vector<double>  *El_topoetcone30; //!
    vector<double>  *El_topoetcone40; //!
    vector<bool>    *El_passIsoLooseTO; //!
    vector<bool>    *El_passIsoLoose; //!
    vector<bool>    *El_passIsoTight; //!
    vector<bool>    *El_passIsoGrad; //!
    vector<bool>    *El_passIsoGradCustomTight; //!
    vector<bool>    *El_passIsoGradCustom; //!
    vector<bool>    *El_passIsoGradLoose; //!
    vector<bool>    *El_trigMatch_e12_lhloose_L1EM10VH; //!
    vector<bool>    *El_trigMatch_e17_lhloose; //!
    vector<bool>    *El_trigMatch_e60_lhmedium; //!
    vector<bool>    *El_trigMatch_e24_lhmedium_L1EM20VH; //!
    vector<bool>    *El_trigMatch_e24_lhmedium_iloose_L1EM20VH; //!
    vector<bool>    *El_trigMatch_e26_lhtight_nod0_ivarloose; //!
    vector<bool>    *El_trigMatch_2e12_lhloose_L12EM10VH; //!
    vector<bool>    *El_trigMatch_2e15_lhloose_L12EM10VH; //!
    vector<bool>    *El_trigMatch_2e15_lhvloose_L12EM13VH; //!
    vector<bool>    *El_trigMatch_2e15_lhvloose_nod0_L12EM13VH; //!
    vector<bool>    *El_trigMatch_2e17_lhvloose_nod0; //!
    vector<bool>    *El_trigMatch_e17_lhloose_mu14; //!
    vector<bool>    *El_trigMatch_e17_lhloose_nod0_mu14; //!
    vector<bool>    *El_trigMatch_e24_lhmedium_nod0_ivarloose; //!
    vector<bool>    *El_trigMatch_e24_lhtight_nod0_ivarloose; //!
    vector<bool>    *El_trigMatch_e60_lhmedium_nod0; //!
    vector<float>   *El_passChargeFlipTaggerBDTmedium; //!
    vector<float>   *El_passChargeFlipTaggerBDTloose; //!
    Int_t           NJet; //!
    vector<double>  *Jet_eta; //!
    vector<double>  *Jet_phi; //!
    vector<double>  *Jet_pT; //!
    vector<double>  *Jet_E; //!
    vector<double>  *Jet_quality; //!
    vector<double>  *Jet_JVT; //!
    vector<double>  *Jet_JVTsf; //!
    vector<double>  *Jet_MV2c20; //!
    vector<double>  *Jet_MV2c10; //!
    vector<double>  *Jet_SFw; //!
    vector<int>     *Jet_ConeTruthLabel; //!
    vector<int>     *Jet_PartonTruthLabel; //!
    vector<int>     *Jet_HadronConeExclTruthLabel; //!
    vector<double>  *Jet_deltaR; //!
    vector<int>     *Jet_nTrk; //!
    vector<bool>    *Jet_passOR; //!
    Float_t         totalJVTsf; //!
    Float_t         Etmiss_CST_Etx; //!
    Float_t         Etmiss_CST_Ety; //!
    Float_t         Etmiss_CST_Et; //!
    Float_t         Etmiss_TST_Etx; //!
    Float_t         Etmiss_TST_Ety; //!
    Float_t         Etmiss_TST_Et; //!
    Float_t         Etmiss_TSTterm_Etx; //!
    Float_t         Etmiss_TSTterm_Ety; //!
    Float_t         Etmiss_TSTterm_Et; //!
    Float_t         Etmiss_Truth_Etx; //!
    Float_t         Etmiss_Truth_Ety; //!
    Float_t         Etmiss_Truth_Et; //!
    Int_t           NTruthJet; //!
    vector<double>  *TruthJet_eta; //!
    vector<double>  *TruthJet_phi; //!
    vector<double>  *TruthJet_pT; //!
    vector<double>  *TruthJet_E; //!
    vector<int>     *TruthJet_ConeTruthLabel; //!
    vector<int>     *TruthJet_PartonTruthLabel; //!
    Int_t           NTruthL; //!
    vector<double>  *TruthL_eta; //!
    vector<double>  *TruthL_phi; //!
    vector<double>  *TruthL_pT; //!
    vector<int>     *TruthL_id; //!
    Int_t           SUSY_Spart_pdgId1; //!
    Int_t           SUSY_Spart_pdgId2; //!
    Int_t           SUSY_Gluino_decay1; //!
    Int_t           SUSY_Gluino_decay2; //!
    Float_t         GenFiltHT; //!
    Float_t         GenFiltMET; //!
    Float_t         PV_z; //!
    Int_t           Nvtx; //!
    Float_t         TruthX1; //!
    Float_t         TruthX2; //!
    Float_t         TruthQ; //!
    Float_t         SherpaNjetWeight; //!
    Int_t           TruthPDGID1; //!
    Int_t           TruthPDGID2; //!
    Double_t        normalization; //!
    Float_t         luminosity; //!
    Float_t         cross_section; //!
    Float_t         k_factor; //!
    Float_t         filter_efficiency; //!
    Float_t         cross_section_kfactor_efficiency; //!
    Double_t        event_weight_sum; //!
    Double_t        event_weight; //!
    Double_t        pileup_weight; //!
    Double_t        baseline_weight; //!
    Double_t        signal_weight; //!
    Int_t           run_number; //!
    Float_t         new_AvgMu; //!
    vector<bool>    *El_isBaseline; //!
    vector<bool>    *El_isSignal; //!
    vector<bool>    *El_isZTag; //!
    vector<bool>    *El_isZProbe; //!
    vector<int>     *El_Tag_index; //!
    vector<int>     *El_Probe_index; //!
    vector<double>  *El_ZTandP_mll; //!
    vector<bool>    *El_isZProbe_TagTriggerMatched; //!
    vector<bool>    *El_Tag_Trigger_Matched; //!
    vector<bool>    *El_Probe_Trigger_Matched; //!
    vector<double>  *El_DR_closest_Jet; //!
    vector<bool>    *Mu_isBaseline; //!
    vector<bool>    *Mu_isSignal; //!
    vector<bool>    *Mu_isZTag; //!
    vector<bool>    *Mu_isZProbe; //!
    vector<int>     *Mu_Tag_index; //!
    vector<int>     *Mu_Probe_index; //!
    vector<double>  *Mu_ZTandP_mll; //!
    vector<bool>    *Mu_isZProbe_TagTriggerMatched; //!
    vector<bool>    *Mu_Tag_Trigger_Matched; //!
    vector<bool>    *Mu_Probe_Trigger_Matched; //!
    vector<double>  *Mu_DR_closest_Jet; //!
    Double_t        baseline_mll; //!
    Double_t        signal_mll; //!

    // List of branches
    TBranch        *b_HLT_e24_lhmedium_nod0_ivarloose; //!
    TBranch        *b_HLT_e24_lhtight_nod0_ivarloose; //!
    TBranch        *b_HLT_e24_lhmedium_nod0_L1EM20VH; //!
    TBranch        *b_HLT_e26_lhtight_iloose; //!
    TBranch        *b_HLT_e26_lhtight_ivarloose; //!
    TBranch        *b_HLT_e26_lhtight_nod0_iloose; //!
    TBranch        *b_HLT_e26_lhtight_nod0_ivarloose; //!
    TBranch        *b_HLT_e60_lhmedium; //!
    TBranch        *b_HLT_e60_lhmedium_nod0; //!
    TBranch        *b_HLT_e120_lhloose_nod0; //!
    TBranch        *b_HLT_e140_lhloose_nod0; //!
    TBranch        *b_HLT_2e17_lhvloose_nod0; //!
    TBranch        *b_HLT_2e15_lhvloose_nod0_L12EM13VH; //!
    TBranch        *b_HLT_e24_lhmedium_e9_lhmedium; //!
    TBranch        *b_HLT_e24_lhmedium_L1EM20VH; //!
    TBranch        *b_HLT_e24_lhmedium_iloose_L1EM20VH; //!
    TBranch        *b_HLT_e12_lhvloose_L12EM10VH; //!
    TBranch        *b_HLT_e17_lhloose_2e9_lhloose; //!
    TBranch        *b_HLT_mu24_ivarmedium; //!
    TBranch        *b_HLT_mu24_imedium; //!
    TBranch        *b_HLT_mu24_ivarloose; //!
    TBranch        *b_HLT_mu24_iloose; //!
    TBranch        *b_HLT_mu26_ivarmedium; //!
    TBranch        *b_HLT_mu20_ivarmedium_L1MU15; //!
    TBranch        *b_HLT_mu20_imedium_L1MU15; //!
    TBranch        *b_HLT_mu20_ivarloose_L1MU15; //!
    TBranch        *b_HLT_mu20_iloose_L1MU15; //!
    TBranch        *b_HLT_mu20_L1MU15; //!
    TBranch        *b_HLT_mu20_mu8noL1; //!
    TBranch        *b_HLT_mu22_mu8noL1; //!
    TBranch        *b_HLT_mu20_2mu4noL1; //!
    TBranch        *b_HLT_mu22_2mu4noL1; //!
    TBranch        *b_HLT_mu40; //!
    TBranch        *b_HLT_mu50; //!
    TBranch        *b_HLT_2mu10; //!
    TBranch        *b_HLT_2mu10_nomucomb; //!
    TBranch        *b_HLT_2mu14; //!
    TBranch        *b_HLT_2mu14_nomucomb; //!
    TBranch        *b_HLT_3mu6; //!
    TBranch        *b_HLT_3mu6_msonly; //!
    TBranch        *b_HLT_xe100_L1XE50; //!
    TBranch        *b_HLT_xe80_mht_L1XE50; //!
    TBranch        *b_HLT_xe90_mht_L1XE50; //!
    TBranch        *b_HLT_xe100_mht_L1XE50; //!
    TBranch        *b_HLT_xe100_tc_em_L1XE50; //!
    TBranch        *b_HLT_xe110_mht_L1XE50; //!
    TBranch        *b_HLT_xe80_tc_lcw_L1XE50; //!
    TBranch        *b_HLT_xe90_tc_lcw_L1XE50; //!
    TBranch        *b_HLT_xe80_tc_lcw_wEFMu_L1XE50; //!
    TBranch        *b_HLT_e7_lhmedium_mu24; //!
    TBranch        *b_HLT_e7_lhmedium_nod0_mu24; //!
    TBranch        *b_HLT_e17_lhloose_nod0_mu14; //!
    TBranch        *b_HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1; //!
    TBranch        *b_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1; //!
    TBranch        *b_HLT_2e12_lhloose_L12EM10VH; //!
    TBranch        *b_HLT_e17_lhloose_mu14; //!
    TBranch        *b_HLT_mu18_mu8noL1; //!
    TBranch        *b_HLT_xe70; //!
    TBranch        *b_EventNumber; //!
    TBranch        *b_ChannelNumber; //!
    TBranch        *b_AvgMu; //!
    TBranch        *b_EventWeight; //!
    TBranch        *b_PRWWeight; //!
    TBranch        *b_PRWrandomRunNumber; //!
    TBranch        *b_TriggerDileptonSF; //!
    TBranch        *b_bcid; //!
    TBranch        *b_LB; //!
    TBranch        *b_RunNb; //!
    TBranch        *b_DetError; //!
    TBranch        *b_NMu; //!
    TBranch        *b_Mu_eta; //!
    TBranch        *b_Mu_phi; //!
    TBranch        *b_Mu_pT; //!
    TBranch        *b_Mu_SFw; //!
    TBranch        *b_Mu_IsoSFw; //!
    TBranch        *b_Mu_charge; //!
    TBranch        *b_Mu_d0pvtx; //!
    TBranch        *b_Mu_sigd0; //!
    TBranch        *b_Mu_z0sinTheta; //!
    TBranch        *b_Mu_isBad; //!
    TBranch        *b_Mu_passOR; //!
    TBranch        *b_Mu_isTight; //!
    TBranch        *b_Mu_isCosmic; //!
    TBranch        *b_Mu_type; //!
    TBranch        *b_Mu_origin; //!
    TBranch        *b_Mu_ptcone20; //!
    TBranch        *b_Mu_ptcone30; //!
    TBranch        *b_Mu_ptcone40; //!
    TBranch        *b_Mu_ptvarcone20; //!
    TBranch        *b_Mu_ptvarcone30; //!
    TBranch        *b_Mu_ptvarcone40; //!
    TBranch        *b_Mu_topoetcone20; //!
    TBranch        *b_Mu_topoetcone30; //!
    TBranch        *b_Mu_topoetcone40; //!
    TBranch        *b_Mu_passIsoLooseTO; //!
    TBranch        *b_Mu_passIsoLoose; //!
    TBranch        *b_Mu_passIsoTight; //!
    TBranch        *b_Mu_passIsoGrad; //!
    TBranch        *b_Mu_passIsoGradCustomTight; //!
    TBranch        *b_Mu_passIsoGradCustom; //!
    TBranch        *b_Mu_passIsoGradLoose; //!
    TBranch        *b_MuTrigSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50; //!
    TBranch        *b_Mu_trigMatch_mu26_ivarmedium; //!
    TBranch        *b_Mu_trigMatch_mu20_iloose_L1MU15; //!
    TBranch        *b_Mu_trigMatch_mu40; //!
    TBranch        *b_Mu_trigMatch_mu50; //!
    TBranch        *b_Mu_trigMatch_mu8noL1; //!
    TBranch        *b_Mu_trigMatch_mu14; //!
    TBranch        *b_Mu_trigMatch_mu18; //!
    TBranch        *b_Mu_trigMatch_mu18_mu8noL1; //!
    TBranch        *b_Mu_trigMatch_e17_lhloose_mu14; //!
    TBranch        *b_Mu_trigMatch_e17_lhloose_nod0_mu14; //!
    TBranch        *b_Mu_trigMatch_mu20_mu8noL1; //!
    TBranch        *b_Mu_trigMatch_mu22_mu8noL1; //!
    TBranch        *b_Mu_trigMatch_mu24_iloose; //!
    TBranch        *b_Mu_trigMatch_mu24_ivarloose; //!
    TBranch        *b_Mu_trigMatch_mu24_iloose_L1MU15; //!
    TBranch        *b_Mu_trigMatch_mu24_ivarloose_L1MU15; //!
    TBranch        *b_Mu_trigMatchPair_mu18_mu8noL1; //!
    TBranch        *b_Mu_trigMatchPair_mu20_mu8noL1; //!
    TBranch        *b_Mu_trigMatchPair_mu22_mu8noL1; //!
    TBranch        *b_NEl; //!
    TBranch        *b_El_eta; //!
    TBranch        *b_El_etaclus; //!
    TBranch        *b_El_phi; //!
    TBranch        *b_El_pT; //!
    TBranch        *b_El_E; //!
    TBranch        *b_El_charge; //!
    TBranch        *b_El_sigd0; //!
    TBranch        *b_El_z0sinTheta; //!
    TBranch        *b_El_d0pvtx; //!
    TBranch        *b_El_passOR; //!
    TBranch        *b_El_SFwMediumLH; //!
    TBranch        *b_El_IsoSFwMediumLH; //!
    TBranch        *b_El_SFwTightLH; //!
    TBranch        *b_El_SFwLooseAndBLayerLH; //!
    TBranch        *b_El_SFwTrigMediumLH_e12_lhloose_L1EM10VH; //!
    TBranch        *b_El_SFwTrigMediumLH_e17_lhloose; //!
    TBranch        *b_El_SFwTrigMediumLH_single; //!
    TBranch        *b_El_SFwTrigLooseAndBLayerLH_e12_lhloose_L1EM10VH; //!
    TBranch        *b_El_isLooseAndBLayerLH; //!
    TBranch        *b_El_isMediumLH; //!
    TBranch        *b_El_isTightLH; //!
    TBranch        *b_El_nBLayerHits; //!
    TBranch        *b_El_expectBLayerHit; //!
    TBranch        *b_El_truthType; //!
    TBranch        *b_El_truthOrigin; //!
    TBranch        *b_El_bkgTruthType; //!
    TBranch        *b_El_bkgTruthOrigin; //!
    TBranch        *b_El_bkgMotherPdgId; //!
    TBranch        *b_El_firstEgMotherTruthType; //!
    TBranch        *b_El_firstEgMotherTruthOrigin; //!
    TBranch        *b_El_firstEgMotherPdgId; //!
    TBranch        *b_El_chFlip; //!
    TBranch        *b_El_ptcone20; //!
    TBranch        *b_El_ptcone30; //!
    TBranch        *b_El_ptcone40; //!
    TBranch        *b_El_ptvarcone20; //!
    TBranch        *b_El_ptvarcone30; //!
    TBranch        *b_El_ptvarcone40; //!
    TBranch        *b_El_topoetcone20; //!
    TBranch        *b_El_topoetcone30; //!
    TBranch        *b_El_topoetcone40; //!
    TBranch        *b_El_passIsoLooseTO; //!
    TBranch        *b_El_passIsoLoose; //!
    TBranch        *b_El_passIsoTight; //!
    TBranch        *b_El_passIsoGrad; //!
    TBranch        *b_El_passIsoGradCustomTight; //!
    TBranch        *b_El_passIsoGradCustom; //!
    TBranch        *b_El_passIsoGradLoose; //!
    TBranch        *b_El_trigMatch_e12_lhloose_L1EM10VH; //!
    TBranch        *b_El_trigMatch_e17_lhloose; //!
    TBranch        *b_El_trigMatch_e60_lhmedium; //!
    TBranch        *b_El_trigMatch_e24_lhmedium_L1EM20VH; //!
    TBranch        *b_El_trigMatch_e24_lhmedium_iloose_L1EM20VH; //!
    TBranch        *b_El_trigMatch_e26_lhtight_nod0_ivarloose; //!
    TBranch        *b_El_trigMatch_2e12_lhloose_L12EM10VH; //!
    TBranch        *b_El_trigMatch_2e15_lhloose_L12EM10VH; //!
    TBranch        *b_El_trigMatch_2e15_lhvloose_L12EM13VH; //!
    TBranch        *b_El_trigMatch_2e15_lhvloose_nod0_L12EM13VH; //!
    TBranch        *b_El_trigMatch_2e17_lhvloose_nod0; //!
    TBranch        *b_El_trigMatch_e17_lhloose_mu14; //!
    TBranch        *b_El_trigMatch_e17_lhloose_nod0_mu14; //!
    TBranch        *b_El_trigMatch_e24_lhmedium_nod0_ivarloose; //!
    TBranch        *b_El_trigMatch_e24_lhtight_nod0_ivarloose; //!
    TBranch        *b_El_trigMatch_e60_lhmedium_nod0; //!
    TBranch        *b_El_passChargeFlipTaggerBDTmedium; //!
    TBranch        *b_El_passChargeFlipTaggerBDTloose; //!
    TBranch        *b_NJet; //!
    TBranch        *b_Jet_eta; //!
    TBranch        *b_Jet_phi; //!
    TBranch        *b_Jet_pT; //!
    TBranch        *b_Jet_E; //!
    TBranch        *b_Jet_quality; //!
    TBranch        *b_Jet_JVT; //!
    TBranch        *b_Jet_JVTsf; //!
    TBranch        *b_Jet_MV2c20; //!
    TBranch        *b_Jet_MV2c10; //!
    TBranch        *b_Jet_SFw; //!
    TBranch        *b_Jet_ConeTruthLabel; //!
    TBranch        *b_Jet_PartonTruthLabel; //!
    TBranch        *b_Jet_HadronConeExclTruthLabel; //!
    TBranch        *b_Jet_deltaR; //!
    TBranch        *b_Jet_nTrk; //!
    TBranch        *b_Jet_passOR; //!
    TBranch        *b_totalJVTsf; //!
    TBranch        *b_Etmiss_CST_Etx; //!
    TBranch        *b_Etmiss_CST_Ety; //!
    TBranch        *b_Etmiss_CST_Et; //!
    TBranch        *b_Etmiss_TST_Etx; //!
    TBranch        *b_Etmiss_TST_Ety; //!
    TBranch        *b_Etmiss_TST_Et; //!
    TBranch        *b_Etmiss_TSTterm_Etx; //!
    TBranch        *b_Etmiss_TSTterm_Ety; //!
    TBranch        *b_Etmiss_TSTterm_Et; //!
    TBranch        *b_Etmiss_Truth_Etx; //!
    TBranch        *b_Etmiss_Truth_Ety; //!
    TBranch        *b_Etmiss_Truth_Et; //!
    TBranch        *b_NTruthJet; //!
    TBranch        *b_TruthJet_eta; //!
    TBranch        *b_TruthJet_phi; //!
    TBranch        *b_TruthJet_pT; //!
    TBranch        *b_TruthJet_E; //!
    TBranch        *b_TruthJet_ConeTruthLabel; //!
    TBranch        *b_TruthJet_PartonTruthLabel; //!
    TBranch        *b_NTruthL; //!
    TBranch        *b_TruthL_eta; //!
    TBranch        *b_TruthL_phi; //!
    TBranch        *b_TruthL_pT; //!
    TBranch        *b_TruthL_id; //!
    TBranch        *b_SUSY_Spart_pdgId1; //!
    TBranch        *b_SUSY_Spart_pdgId2; //!
    TBranch        *b_SUSY_Gluino_decay1; //!
    TBranch        *b_SUSY_Gluino_decay2; //!
    TBranch        *b_GenFiltHT; //!
    TBranch        *b_GenFiltMET; //!
    TBranch        *b_PV_z; //!
    TBranch        *b_Nvtx; //!
    TBranch        *b_TruthX1; //!
    TBranch        *b_TruthX2; //!
    TBranch        *b_TruthQ; //!
    TBranch        *b_SherpaNjetWeight; //!
    TBranch        *b_TruthPDGID1; //!
    TBranch        *b_TruthPDGID2; //!
    TBranch        *b_normalization; //!
    TBranch        *b_luminosity; //!
    TBranch        *b_cross_section; //!
    TBranch        *b_k_factor; //!
    TBranch        *b_filter_efficiency; //!
    TBranch        *b_cross_section_kfactor_efficiency; //!
    TBranch        *b_event_weight_sum; //!
    TBranch        *b_event_weight; //!
    TBranch        *b_pileup_weight; //!
    TBranch        *b_baseline_weight; //!
    TBranch        *b_signal_weight; //!
    TBranch        *b_run_number; //!
    TBranch        *b_new_AvgMu; //!
    TBranch        *b_El_isBaseline; //!
    TBranch        *b_El_isSignal; //!
    TBranch        *b_El_isZTag; //!
    TBranch        *b_El_isZProbe; //!
    TBranch        *b_El_Tag_index; //!
    TBranch        *b_El_Probe_index; //!
    TBranch        *b_El_ZTandP_mll; //!
    TBranch        *b_El_isZProbe_TagTriggerMatched; //!
    TBranch        *b_El_Tag_Trigger_Matched; //!
    TBranch        *b_El_Probe_Trigger_Matched; //!
    TBranch        *b_El_DR_closest_Jet; //!
    TBranch        *b_Mu_isBaseline; //!
    TBranch        *b_Mu_isSignal; //!
    TBranch        *b_Mu_isZTag; //!
    TBranch        *b_Mu_isZProbe; //!
    TBranch        *b_Mu_Tag_index; //!
    TBranch        *b_Mu_Probe_index; //!
    TBranch        *b_Mu_ZTandP_mll; //!
    TBranch        *b_Mu_isZProbe_TagTriggerMatched; //!
    TBranch        *b_Mu_Tag_Trigger_Matched; //!
    TBranch        *b_Mu_Probe_Trigger_Matched; //!
    TBranch        *b_Mu_DR_closest_Jet; //!
    TBranch        *b_baseline_mll; //!
    TBranch        *b_signal_mll; //!

    // Histograms
    // Unweighted histograms are for data
    TH1F *h_Nvtx; //!
    TH1F *h_AvgMu; //!
    TH1F *h_AvgMu_OSee; //!
    TH1F *h_AvgMu_OSmumu; //!
    TH1F *h_zPV; //!
    //TH1F *h_NLepts; //!
    //TH1F *h_NJets; //!

    TH1F *h_baseline_OSee_mll; //!
    TH1F *h_baseline_OSmumu_mll; //!
    TH1F *h_signal_OSee_mll; //!
    TH1F *h_signal_OSmumu_mll; //!

    // Baseline histograms
    // 1-dim histograms
    TH1F *h_baseline_mll; //!
    TH1F *h_baseline_pt; //!
    TH1F *h_baseline_eta; //!
    //TH1F *h_baseline_d0pvtx; //!
    //TH1F *h_baseline_sigd0; //!
    //TH1F *h_baseline_d0err; //!
    //TH1F *h_baseline_z0sinTheta; //!
    //TH1F *h_baseline_ptvarcone20; //!
    //TH1F *h_baseline_ptvarcone30; //!
    //TH1F *h_baseline_topoetcone20; //!
    TH1F *h_baseline_nJets; //!
    TH1F *h_baseline_dRjet; //!
    //TH1F *h_baseline_Etmiss; //!
    //TH1F *h_baseline_meff; //!

    TH1F *h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt; //!
    TH1F *h_baseline_CaloIso_cut_pt; //!
    TH1F *h_baseline_TrackIso_cut_pt; //!
    TH1F *h_baseline_Iso_cut_pt; //!
    TH1F *h_baseline_z0_cut_pt; //!
    TH1F *h_baseline_sigd0_cut_pt; //!

    // 2-dim histograms
    TH2F *h_baseline_pt_eta; //! x: pt, y: eta
    TH2F *h_baseline_pt_nJets; //! x: pt, y: nJets
    TH2F *h_baseline_pt_dRjet; //! x: pt, y: dRjet
    //TH2F *h_baseline_pt_Etmiss; //! x: pt, y: Etmiss
    //TH2F *h_baseline_pt_meff; //! x: pt, y: meff

    // 3-dim histograms
    TH3F *h_baseline_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_baseline_pt_eta_dRjet; //! x: pt, y: eta, z: dRjet
    TH3F *h_baseline_pt_dRjet_mll; //! x: pt, y: dRjet, z: mll
    TH3F *h_baseline_eta_dRjet_mll; //! x: eta, y: dRjet, z: mll

    TH3F *h_baseline_pt_eta_AvgMu; //! x: pt, y: eta, z: AvgMu
    TH3F *h_baseline_pt_dRjet_AvgMu; //! x: pt, y: dRjet, z: AvgMu
    TH3F *h_baseline_eta_dRjet_AvgMu; //! x: eta, y: dRjet, z: AvgMu

    TH3F *h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_baseline_CaloIso_cut_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_baseline_TrackIso_cut_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_baseline_Iso_cut_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_baseline_z0_cut_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_baseline_sigd0_cut_pt_eta_mll; //! x: pt, y: eta, z: mll

    // Signal histograms
    // 1-dim histograms
    TH1F *h_signal_mll; //!
    TH1F *h_signal_pt; //!
    TH1F *h_signal_eta; //!
    //TH1F *h_signal_d0pvtx; //!
    //TH1F *h_signal_sigd0; //!
    //TH1F *h_signal_d0err; //!
    //TH1F *h_signal_z0sinTheta; //!
    //TH1F *h_signal_ptvarcone20; //!
    //TH1F *h_signal_ptvarcone30; //!
    //TH1F *h_signal_topoetcone20; //!
    TH1F *h_signal_nJets; //!
    TH1F *h_signal_dRjet; //!
    //TH1F *h_signal_Etmiss; //!
    //TH1F *h_signal_meff; //!

    // 2-dim histograms
    TH2F *h_signal_pt_eta; //! x: pt, y: eta
    TH2F *h_signal_pt_nJets; //! x: pt, y: nJets
    TH2F *h_signal_pt_dRjet; //! x: pt, y: dRjet
    //TH2F *h_signal_pt_Etmiss; //! x: pt, y: Etmiss
    //TH2F *h_signal_pt_meff; //! x: pt, y: meff

    // 3-dim histograms
    TH3F *h_signal_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_signal_pt_eta_dRjet; //! x: pt, y: eta, z: dRjet
    TH3F *h_signal_pt_dRjet_mll; //! x: pt, y: dRjet, z: mll
    TH3F *h_signal_eta_dRjet_mll; //! x: eta, y: dRjet, z: mll

    TH3F *h_signal_pt_eta_AvgMu; //! x: pt, y: eta, z: AvgMu
    TH3F *h_signal_pt_dRjet_AvgMu; //! x: pt, y: dRjet, z: AvgMu
    TH3F *h_signal_eta_dRjet_AvgMu; //! x: eta, y: dRjet, z: AvgMu

    // Background template histograms
    TH1F *h_bkg_template_fail_truth; //!
    TH1F *h_bkg_template_fail_id_only; //!
    TH1F *h_bkg_template_fail_CaloIso_only; //!
    TH1F *h_bkg_template_fail_TrackIso_only; //!
    TH1F *h_bkg_template_fail_id_and_CaloIso; //!
    TH1F *h_bkg_template_fail_id_and_TrackIso; //!
    TH1F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso; //!
    TH1F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight; //!
    TH1F *h_bkg_template_fail_CaloIso_and_TrackIso; //!

    // 3-dim background template histograms
    TH3F *h_bkg_template_fail_truth_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_id_only_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_CaloIso_only_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_TrackIso_only_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_id_and_CaloIso_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_id_and_TrackIso_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll; //! x: pt, y: eta, z: mll
    TH3F *h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll; //! x: pt, y: eta, z: mll

    // Cut efficiency histograms
    TH1F *h_cut_eff_LooseAndBLayerLLH_to_MediumLLH; //!
    TH1F *h_cut_eff_CaloIso; //!
    TH1F *h_cut_eff_TrackIso; //!
    TH1F *h_cut_eff_Iso; //!
    TH1F *h_cut_eff_z0; //!
    TH1F *h_cut_eff_sigd0; //!

    // Real lepton efficiency histograms
    TH1F *h_eff_pt; //!
    TH1F *h_eff_eta; //!
    TH1F *h_eff_nJets; //!
    TH1F *h_eff_dRjet; //!
    //TH1F *h_eff_Etmiss; //!
    //TH1F *h_eff_meff; //!

    TCanvas *baseline_mll_plot; //!
    TCanvas *cut_efficiency_plot; //!
    TCanvas *real_efficiency_plot; //!
    //TLine *line; //!
    //TLegend *leg; //!

    // Binning
    vector<float> m_mll_bins; //!
    vector<float> m_pt_bins; //!
    vector<float> m_eta_bins; //!
    vector<float> m_deltaR_bins; //!
    vector<float> m_NJet_bins; //!
    //vector<float> m_Etmiss_bins; //!
    //vector<float> m_meff_bins; //!
    vector<float> m_AvgMu_bins; //!

    int n_mll_bins; //!
    int n_pt_bins; //!
    int n_eta_bins; //!
    int n_deltaR_bins; //!
    int n_NJet_bins; //!
    //int n_Etmiss_bins; //!
    //int n_meff_bins; //!
    int n_AvgMu_bins; //!

    static const float m_default_mll_bins[91]; //
    static const float m_default_pt_bins[14]; //!
    static const float m_default_eta_bins[21]; //!
    static const float m_default_deltaR_bins[21]; //!
    static const float m_default_NJet_bins[11]; //!
    //static const float m_default_Etmiss_bins[9]; //!
    //static const float m_default_meff_bins[7]; //!
    static const float m_default_AvgMu_bins[51]; //!

    static const float m_muon_eta_bins[21]; //!
    static const float m_abs_eta_bins[11]; //!
    static const float m_muon_abs_eta_bins[11]; //!

    static const float m_electron_coarse_pt_bins[12]; //!
    static const float m_electron_coarse_eta_bins[6]; //!
    static const float m_muon_coarse_pt_bins[12]; //!
    static const float m_muon_coarse_eta_bins[5]; //!
    static const float m_coarse_deltaR_bins[7]; //!
    static const float m_fine_deltaR_bins[36]; //!

    // this is a standard constructor
    ytRealLeptonsEfficiency ();

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

    // Set methods
    void set_isMC(bool b) { isMC = b; cout << "isMC=" << isMC << endl; }
    void set_isData(bool b) { isData = b; cout << "isData=" << isData << endl; }
    void set_trigger(TString str) { trigger = str; cout << "trigger=" << trigger << endl; }
    void set_tag_trigger_matched(bool b) { tag_trigger_matched = b; cout << "tag_trigger_matched=" << tag_trigger_matched << endl; }
    void set_lepton(TString str) { lepton = str; cout << "lepton=" << lepton << endl; }
    void set_process(TString str) { process = str; cout << "process=" << process << endl; }

    void set_binning_default();
    void set_binning(int i, const float *bins, vector<float> &vec);
    void set_binning_pt(int i, const float *bins) { set_binning(i, bins, m_pt_bins); }
    void set_binning_eta(int i, const float *bins) { set_binning(i, bins, m_eta_bins); }
    void set_binning_mll(int i, const float *bins) { set_binning(i, bins, m_mll_bins); }
    void set_binning_deltaR(int i, const float *bins) { set_binning(i, bins, m_deltaR_bins); }
    void set_binning_AvgMu(int i, const float *bins) { set_binning(i, bins, m_AvgMu_bins); }

    // Background subtraction methods
    void do_background_subtraction() { /* Need to be implemented */ }

    // loop over leptons
    void loop_over_electrons();
    void loop_over_muons();

    // plots
    void make_baseline_mll_plot();
    void make_cut_efficiency_plot();
    void make_electron_real_efficiency_plot();
    void make_muon_real_efficiency_plot();

    // this is needed to distribute the algorithm to the workers
    ClassDef(ytRealLeptonsEfficiency, 1);
};

#endif
