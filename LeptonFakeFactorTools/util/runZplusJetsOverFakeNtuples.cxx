#include <iostream>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include <TFile.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TSystem.h>

#include "CxxUtils/make_unique.h"
#include <LeptonFakeFactorTools/ZplusJetsFakeFactorCalculator.h>
#include <LeptonFakeFactorTools/EventDef.h>

// suitable ONLY for a TTreeReader named "reader"
// branchName is the name of the branch in your ntuple, 
// if it differs from the variable name needed for the 
// fake factor tool
#ifndef READ_TREE_GENERIC
#define READ_TREE_GENERIC(typeName, name, branchName) \
  TTreeReaderValue<typeName> name(reader, #branchName)
#endif

// suitable ONLY for a TTreeReader named "reader"
#ifndef READ_TREE
#define READ_TREE(typeName, name) \
  READ_TREE_GENERIC(typeName, name, name)
#endif

int main(int argc, char* argv[]){

  bool writeHistFitterTree = false;

  if(argc < 2){
    std::cerr << "No input file specified! Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  TString fileName = argv[1];
  TFile *inFile = TFile::Open(fileName, "READ");
  if(!inFile){
    std::cerr << "Invalid root file! Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // Name of the sample: should do something a bit fancier, since a data file
  // named anything but data.root will not be handled properly
  TString sample = "";
  fileName = fileName.Remove(0, fileName.Last('/')+1); // remove the full path
  sample = fileName.ReplaceAll(".root","");

  std::cout << "Running sample: " << sample << std::endl;

  // WLep, ZLeadingLep, ZSubleadingLep
  TString subregion = "ttt";
  if(argc > 2){
    subregion = argv[2];
  }

  TString configFile = "$ROOTCOREBIN/../LeptonFakeFactorTools/data/Inclusive/Zjet_FF_calculation.conf";
  if(argc > 3){
    configFile = argv[3];
  }
  gSystem->ExpandPathName(configFile);

  TString fakeFactorType = "Nominal";
  if(argc > 4){
    fakeFactorType = argv[4];
  }

  std::cout << "subregion is " << subregion << std::endl;
  std::cout << "configFile is " << configFile << std::endl;

  ZplusJetsFakeFactorCalculator* ffCalcTool = new ZplusJetsFakeFactorCalculator();

  // Generic things
  ffCalcTool->set_sample(sample);
  ffCalcTool->set_subregion(subregion);
  ffCalcTool->set_configFile(configFile);
  ffCalcTool->msg().setLevel(MSG::INFO);
  if(!ffCalcTool->initialize().isSuccess()) return EXIT_FAILURE;

  // Use TTreeReader to read in the events. For more info see:
  // https://root.cern.ch/doc/master/classTTreeReader.html
  // In particular note that to get the Mll float value you'll
  // need to do *Mll etc.
  TTreeReader reader("tree_NoSys", inFile);
  READ_TREE(int,                   RunNumber);
  READ_TREE(float,                 weight);
  READ_TREE(float,                 metEt);
  READ_TREE(float,                 metPhi);
  READ_TREE(int,                   nJets);
  READ_TREE(int,                   nBJets);
  READ_TREE(float,                 leadJetPt);
  READ_TREE(std::vector<float>,    lepPt);
  READ_TREE(std::vector<float>,    lepEta);
  READ_TREE(std::vector<float>,    lepPhi);
  READ_TREE(std::vector<float>,    lepD0Sig);
  READ_TREE(std::vector<float>,    lepZ0SinTheta);
  READ_TREE(std::vector<int>,      lepFlavor);
  READ_TREE(std::vector<int>,      lepCharge);
  READ_TREE(std::vector<int>,      lepPassOR);
  READ_TREE(std::vector<int>,      lepTruthType);
  READ_TREE(std::vector<int>,      lepTruthOrigin);
  READ_TREE(std::vector<int>,      lepPassBL);
  READ_TREE(std::vector<int>,      lepNPix);
  READ_TREE(std::vector<int>,      lepVeryLoose);
  READ_TREE(std::vector<int>,      lepLoose);
  READ_TREE(std::vector<int>,      lepMedium);
  READ_TREE(std::vector<int>,      lepTight);
  READ_TREE(std::vector<int>,      lepIsoLoose);
  READ_TREE(std::vector<int>,      lepIsoTight);
  READ_TREE(std::vector<int>,      lepIsoGradient);
  READ_TREE(std::vector<int>,      lepIsoGradientLoose);
  READ_TREE(std::vector<int>,      lepIsoLooseTrackOnly);
  READ_TREE(std::vector<int>,      lepSignal);
  READ_TREE(std::vector<int>,      lepBaseline);
  READ_TREE(std::vector<int>,      lepMatchesTrigger);

  // For output tree
  int firstLepFlavor;
  int firstLepCharge; 
  int firstLepType;
  int firstLepOrigin; 
  float firstLepPt;
  float firstLepEta;
  float firstLepPhi;
  int secondLepFlavor;
  int secondLepCharge; 
  int secondLepType;
  int secondLepOrigin; 
  float secondLepPt;
  float secondLepEta;
  float secondLepPhi;
  int thirdLepFlavor;
  int thirdLepCharge; 
  int thirdLepType;
  int thirdLepOrigin; 
  float thirdLepPt;
  float thirdLepEta;
  float thirdLepPhi;
  float eventWeight; 
  float eventWeightSyst; 
  float eventMT; 
  float eventMll; 
  float eventMet; 
  int eventNJets;

  TFile outfile(sample+"_"+subregion+"_"+"thirdLep.root","RECREATE");
  TTree* outtree = new TTree("thirdLep","thirdLep");

  outtree->Branch("Flavor1", &firstLepFlavor);
  outtree->Branch("Charge1", &firstLepCharge);
  outtree->Branch("Type1", &firstLepType);
  outtree->Branch("Origin1", &firstLepOrigin);
  outtree->Branch("Pt1", &firstLepPt);
  outtree->Branch("Eta1", &firstLepEta);
  outtree->Branch("Phi1", &firstLepPhi);
  outtree->Branch("Flavor2", &secondLepFlavor);
  outtree->Branch("Charge2", &secondLepCharge);
  outtree->Branch("Type2", &secondLepType);
  outtree->Branch("Origin2", &secondLepOrigin);
  outtree->Branch("Pt2", &secondLepPt);
  outtree->Branch("Eta2", &secondLepEta);
  outtree->Branch("Phi2", &secondLepPhi);
  outtree->Branch("Flavor3", &thirdLepFlavor);
  outtree->Branch("Charge3", &thirdLepCharge);
  outtree->Branch("Type3", &thirdLepType);
  outtree->Branch("Origin3", &thirdLepOrigin);
  outtree->Branch("Pt3", &thirdLepPt);
  outtree->Branch("Eta3", &thirdLepEta);
  outtree->Branch("Phi3", &thirdLepPhi);
  outtree->Branch("weight", &eventWeight);
  outtree->Branch("weightSyst", &eventWeightSyst);
  outtree->Branch("MT", &eventMT);
  outtree->Branch("Mll", &eventMll);
  outtree->Branch("MET", &eventMet);
  outtree->Branch("nJets", &eventNJets);

  TFile* HF_File = NULL;
  TTree* HF_Tree = NULL;

  // Dummy vars to set default values
  float  HF_dummyFloat_n999 = -999;
  float  HF_dummyFloat_n1 = -1;
  float  HF_dummyFloat_0 = 0;
  float  HF_dummyFloat_1 = 1;
  double HF_dummyDouble = -999;
  int    HF_dummyInt_n999 = -999;
  int    HF_dummyInt_n1 = -1;
  int    HF_dummyInt_0 = 0;
  bool   HF_dummyBool = false;

  // Leptons
  float HF_lept1Pt;
  float HF_lept2Pt;
  float HF_lept3Pt;
  float HF_lept1Eta;
  float HF_lept2Eta;
  float HF_lept3Eta;
  float HF_lept1Phi;
  float HF_lept2Phi;
  float HF_lept3Phi;
  float HF_lept1q;
  float HF_lept2q;
  float HF_lept3q;
  int HF_lept1Flav;
  int HF_lept2Flav;
  int HF_lept3Flav;

  // Jets
  float HF_jet1Pt;
  int   HF_njets;
  int   HF_nbjets;
  int   HF_nCentralLightJets;

  // MET
  float HF_MET;
  float HF_Mt;

  // Others
  int    HF_nSigLep;
  double HF_eventweight; // Has to be the same for all regions, must include the scale factor to 1fb^-1
  int    HF_runNumber;

  // Trigger
  bool HF_passtrigger = true; // preselection for the ntuples
  bool HF_passtriggermatch = true; // if this failed, then the FF tool wouldn't have returned true

  // three-lepton specific
  float HF_leptPtSumVec;
  float HF_Mlll;
  float HF_bestZcandidate;
  float HF_L3Mll;
  float HF_L3Mt;

  // FF syst
  double HF_syst_FF;
  

  if(writeHistFitterTree){
    HF_File = new TFile(sample+"_"+subregion+"_HFT.root","RECREATE");
    HF_Tree = new TTree("FF_CENTRAL", "FF_CENTRAL");

    std::vector<TString> dummyFloatVars_0 = {
       "lept4Pt", "lept4Eta", "lept4Phi", "lept4q", "jet2Pt", "jet3Pt", "bJet1Pt", "bJet2Pt"
      ,"L2Mll"  ,"L2MT2", "L2dileptonpt", "L2mJJ", "L2dPhiLL", "L2dRLL"
      ,"L2dRJJ", "L2dPhimetjets", "L2dPhiWZ", "L2dPhimetW", "L2dPhimetZ"};

    std::vector<TString> dummyFloatVars_n1 = {"Meff"};

    std::vector<TString> dummyFloatVars_n999 = {
       "RJPTCM", "RJPTISR", "RJPTI", "RJRISR", "RJdphiISRI", "RJMZ", "RJMJ"
      ,"RJPTCM_VR", "RJPTISR_VR", "RJPTI_VR", "RJRISR_VR", "RJdphiISRI_VR"
      ,"RJMZ_VR", "RJMJ_VR", "RJH2PP", "RJH2PP_VR", "RJH4PP", "RJHT4PP"
      ,"RJH5PP", "RJH5PP_VR", "RJHT5PP", "RJHT5PP_VR", "RJR_minH2P_minH3P"
      ,"RJR_minH2P_minH3P_VR", "RJR_HT4PP_H4PP", "RJRPT_HT5PP"
      ,"RJRPT_HT5PP_VR", "RJdphiVP", "RJdphiVP_VR", "RJdphiPPV"};

    std::vector<TString> dummyFloatVars_1 = {
       "syst_FT_EFF_B_down", "syst_FT_EFF_B_up", "syst_FT_EFF_C_down"
      ,"syst_FT_EFF_C_up", "syst_FT_EFF_Light_down", "syst_FT_EFF_Light_up"
      ,"syst_FT_EFF_extrapolation_down", "syst_FT_EFF_extrapolation_up"
      ,"syst_FT_EFF_extrapolationFromCharm_down"
      ,"syst_FT_EFF_extrapolationFromCharm_up"
      ,"syst_EL_EFF_ID_TOTAL_UncorrUncertainty_down"
      ,"syst_EL_EFF_ID_TOTAL_UncorrUncertainty_up"
      ,"syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_down"
      ,"syst_EL_EFF_Iso_TOTAL_UncorrUncertainty_up"
      ,"syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_down"
      ,"syst_EL_EFF_Reco_TOTAL_UncorrUncertainty_up"
      ,"syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_down"
      ,"syst_EL_EFF_Trigger_TOTAL_UncorrUncertainty_up"
      ,"syst_MUON_EFF_STAT_down", "syst_MUON_EFF_STAT_up"
      ,"syst_MUON_EFF_SYS_down", "syst_MUON_EFF_SYS_up"
      ,"syst_MUON_EFF_TrigStatUncertainty_down"
      ,"syst_MUON_EFF_TrigStatUncertainty_up"
      ,"syst_MUON_EFF_TrigSystUncertainty_down"
      ,"syst_MUON_EFF_TrigSystUncertainty_up", "syst_MUON_ISO_STAT_down"
      ,"syst_MUON_ISO_STAT_up", "syst_MUON_ISO_SYS_down"
      ,"syst_MUON_ISO_SYS_up", "syst_MM_WGT_down", "syst_MM_WGT_up"
      ,"syst_MM_STAT_down", "syst_MM_STAT_up"};

    std::vector<TString> dummyDoubleVars = {"L2dPhimetLL"};

    std::vector<TString> dummyIntVars_n1 = {"lept4Flav", "eventNumber", "mcChannel"};

    std::vector<TString> dummyIntVars_0 = {
       "L2nCentralLightJets"
      ,"L2nCentralLightJets30", "L2nCentralLightJets60", "L2nCentralBJets"
      ,"L2nForwardJets", "L2finalState", "SR", "nVtx"};

    std::vector<TString> dummyIntVars_n999 = {
       "RJSR2L1H", "RJSR2L2I", "RJSR2L3C"
      ,"RJSR2L4C", "RJCR2LVVH", "RJVR2LVVH", "RJCR2LTopH", "RJVR2LTopH"
      ,"RJCR2LVVC", "RJVR2LVVC", "RJCR2LTopC", "RJVR2LTopC", "RJVR2LZ2jC"
      ,"RJVR2LZ3jC", "RJSR3L1H", "RJSR3L2I", "RJSR3L3C", "RJSR3L4C"
      ,"RJCR3LVVH", "RJVR3LVVH", "RJCR3LVVC", "RJVR3LVVC", "RJNjS", "RJNjISR"
      ,"RJNjS_VR", "RJNjISR_VR"};

    std::vector<TString> dummyBoolVars = {"L2isEMU", "L2isMUMU", "L2isEE"};

    // Branches which I don't need
    for(TString var : dummyFloatVars_n999){
      HF_Tree->Branch(var, &HF_dummyFloat_n999);
    }
    for(TString var : dummyFloatVars_n1){
      HF_Tree->Branch(var, &HF_dummyFloat_n1);
    }
    for(TString var : dummyFloatVars_0){
      HF_Tree->Branch(var, &HF_dummyFloat_0);
    }
    for(TString var : dummyFloatVars_1){
      HF_Tree->Branch(var, &HF_dummyFloat_1);
    }
    for(TString var : dummyDoubleVars){
      HF_Tree->Branch(var, &HF_dummyDouble);
    }
    for(TString var : dummyIntVars_n999){
      HF_Tree->Branch(var, &HF_dummyInt_n999);
    }
    for(TString var : dummyIntVars_n1){
      HF_Tree->Branch(var, &HF_dummyInt_n1);
    }
    for(TString var : dummyIntVars_0){
      HF_Tree->Branch(var, &HF_dummyInt_0);
    }
    for(TString var : dummyBoolVars){
      HF_Tree->Branch(var, &HF_dummyBool);
    }

    HF_Tree->Branch("lept1Pt", &HF_lept1Pt);
    HF_Tree->Branch("lept2Pt", &HF_lept2Pt);
    HF_Tree->Branch("lept3Pt", &HF_lept3Pt);
    HF_Tree->Branch("lept1Eta", &HF_lept1Eta);
    HF_Tree->Branch("lept2Eta", &HF_lept2Eta);
    HF_Tree->Branch("lept3Eta", &HF_lept3Eta);
    HF_Tree->Branch("lept1Phi", &HF_lept1Phi);
    HF_Tree->Branch("lept2Phi", &HF_lept2Phi);
    HF_Tree->Branch("lept3Phi", &HF_lept3Phi);
    HF_Tree->Branch("lept1q", &HF_lept1q);
    HF_Tree->Branch("lept2q", &HF_lept2q);
    HF_Tree->Branch("lept3q", &HF_lept3q);
    HF_Tree->Branch("lept1Flav", &HF_lept1Flav);
    HF_Tree->Branch("lept2Flav", &HF_lept2Flav);
    HF_Tree->Branch("lept3Flav", &HF_lept3Flav);
    HF_Tree->Branch("jet1Pt", &HF_jet1Pt);
    HF_Tree->Branch("njets", &HF_njets);
    HF_Tree->Branch("nbjets", &HF_nbjets);
    HF_Tree->Branch("nCentralLightJets", &HF_nCentralLightJets);
    HF_Tree->Branch("MET", &HF_MET);
    HF_Tree->Branch("Mt", &HF_Mt);
    HF_Tree->Branch("nSigLep", &HF_nSigLep);
    HF_Tree->Branch("eventweight", &HF_eventweight);
    HF_Tree->Branch("runNumber", &HF_runNumber);
    HF_Tree->Branch("passtrigger", &HF_passtrigger);
    HF_Tree->Branch("passtriggermatch", &HF_passtriggermatch);
    HF_Tree->Branch("leptPtSumVec", &HF_leptPtSumVec);
    HF_Tree->Branch("Mlll", &HF_Mlll);
    HF_Tree->Branch("bestZcandidate", &HF_bestZcandidate);
    HF_Tree->Branch("L3Mll", &HF_L3Mll);
    HF_Tree->Branch("L3Mt", &HF_L3Mt);
    HF_Tree->Branch("syst_FF", &HF_syst_FF);
  }


  // Loop over events!
  int nEventsProcessed = 0;
  int nEventsPassed = 0;

  while(reader.Next()) {
    ++nEventsProcessed;

    // create EventDef class instance whose members are all the default values.
    // (of course, EventDef* evt = new EventDef(); would also work, if .'s replaced with ->'s)
    EventDef evt;

    // Now fill derived quantities for the tool to use
    int nLeps = (*lepPt).size(); // all of the vectors should be the same size...

    // Use MCTruthClassifier here
    std::vector<int> lepIsTruth(nLeps, 0);

    // Look at ID and isolation to fill whether a given lepton passes
    // the signal ID, other ID, or antiID, as well as the appropriate
    // scale factors
    std::vector<int> lepPassSigID(nLeps, 0);
    std::vector<int> lepPassAntiID(nLeps, 0);

    int nBaselineLeptons = 0;

    for(int i = 0; i < nLeps; ++i){

      // Fill ID and SF vectors below:
      if( (*lepFlavor).at(i) == 1 ){ // electrons

        double d0sig = fabs( (*lepD0Sig).at(i) );
        double z0sinTheta = fabs( (*lepZ0SinTheta).at(i) );

        bool passesVeryLoose = (*lepVeryLoose).at(i) && (z0sinTheta < 1.0);
        bool passesAllCuts = (d0sig < 5) && (*lepIsoGradientLoose).at(i) && (*lepMedium).at(i);
        bool passesAntiID = passesVeryLoose && !passesAllCuts;

        bool passesSignal = (*lepSignal).at(i); // can add extra cuts if desired

        if((*lepBaseline).at(i) && (*lepLoose).at(i) && (*lepPassBL).at(i) && (*lepPassOR).at(i)) nBaselineLeptons++;

        // Now fill the vectors
        lepPassSigID.at(i) = passesSignal;
        lepPassAntiID.at(i) = passesAntiID && (*lepPassOR).at(i);

        lepIsTruth.at(i) = ( (*lepTruthType).at(i) == 2 ); // MCTruthClassifier

      }
      else if( (*lepFlavor).at(i) == 2){ // muons

        double d0sig = fabs( (*lepD0Sig).at(i) );
        double z0sinTheta = fabs( (*lepZ0SinTheta).at(i) );

        bool passesMedium = (*lepMedium).at(i) && (z0sinTheta < 1.0);
        bool passesAllCuts = (d0sig < 3) && (*lepIsoGradientLoose).at(i);
        bool passesAntiID = passesMedium && !passesAllCuts;

        bool passesSignal = (*lepSignal).at(i); // can add extra cuts if desired

        if((*lepBaseline).at(i) && (*lepPassOR).at(i)) nBaselineLeptons++;

        // Now fill the vectors
        lepPassSigID.at(i) = passesSignal;
        lepPassAntiID.at(i) = passesAntiID;

        lepIsTruth.at(i) = ( (*lepTruthType).at(i) == 6 ); // MCTruthClassifier
      }
    }

    // fill in the struct
    evt.nBaselineLeptons          = nBaselineLeptons;
    evt.lepFlavor                 = *lepFlavor;
    evt.lepCharge                 = *lepCharge;
    evt.lepIsTruth                = lepIsTruth;
    evt.lepMatchesTrigger         = *lepMatchesTrigger;
    evt.lepPt                     = *lepPt;
    evt.lepEta                    = *lepEta;
    evt.lepPhi                    = *lepPhi;

    evt.met_Et                    = *metEt;
    evt.met_phi                   = *metPhi;

    evt.nJets                     = *nJets;
    evt.nBJets                    = *nBJets;
    evt.leadJetPt                 = *leadJetPt;

    evt.lepPassSigID              = lepPassSigID;
    evt.lepPassOtherID            = lepPassSigID;
    evt.lepPassAntiID             = lepPassAntiID;

    // Note! This is because our ntuple preselection required
    // the dilepton trigger to be fired; we don't use the single
    // lepton trigger in our analysis, so we can set that to false here
    evt.passDileptonTrigger       = true;
    evt.passSingleLeptonTrigger   = false;

    evt.weight                    = *weight;

    if(ffCalcTool->accept(evt)){
      ++nEventsPassed;

      int fakeOrWIndex = evt.lep0_index;
      if((*lepPt).at(evt.lep1_index) < (*lepPt).at(fakeOrWIndex)) fakeOrWIndex = evt.lep1_index;
      if((*lepPt).at(evt.lep2_index) < (*lepPt).at(fakeOrWIndex)) fakeOrWIndex = evt.lep2_index;

      int lepLead_index = -1;
      int lepSubLead_index = -1;
      int lepSubSubLead_index = -1;

      std::vector<int> indices = {evt.lep0_index, evt.lep1_index, evt.lep2_index};

      for(int indexLead : indices){
        for(int indexSubLead : indices){
          if(indexLead == indexSubLead) continue;

          for(int indexSubSubLead : indices){
            if(indexLead == indexSubSubLead)    continue;
            if(indexSubLead == indexSubSubLead) continue;

            if( (*lepPt).at(indexLead) >= (*lepPt).at(indexSubLead) && (*lepPt).at(indexSubLead) >= (*lepPt).at(indexSubSubLead) ){
              lepLead_index = indexLead;
              lepSubLead_index = indexSubLead;
              lepSubSubLead_index = indexSubSubLead;
            }
          }
        }
      }

      if(lepLead_index == -1 || lepSubLead_index == -1 || lepSubSubLead_index == -1){
        std::cerr << "Somehow couldn't sort the indices! Fix this! Aborting." << std::endl;
        abort();
      }


      firstLepFlavor           = (*lepFlavor).at(lepLead_index);
      firstLepCharge           = (*lepCharge).at(lepLead_index);
      firstLepType             = (*lepTruthType).at(lepLead_index);
      firstLepOrigin           = (*lepTruthOrigin).at(lepLead_index);
      firstLepPt               = (*lepPt).at(lepLead_index);
      firstLepEta              = (*lepEta).at(lepLead_index);
      firstLepPhi              = (*lepPhi).at(lepLead_index);
      secondLepFlavor           = (*lepFlavor).at(lepSubLead_index);
      secondLepCharge           = (*lepCharge).at(lepSubLead_index);
      secondLepType             = (*lepTruthType).at(lepSubLead_index);
      secondLepOrigin           = (*lepTruthOrigin).at(lepSubLead_index);
      secondLepPt               = (*lepPt).at(lepSubLead_index);
      secondLepEta              = (*lepEta).at(lepSubLead_index);
      secondLepPhi              = (*lepPhi).at(lepSubLead_index);
      thirdLepFlavor           = (*lepFlavor).at(lepSubSubLead_index);
      thirdLepCharge           = (*lepCharge).at(lepSubSubLead_index);
      thirdLepType             = (*lepTruthType).at(lepSubSubLead_index);
      thirdLepOrigin           = (*lepTruthOrigin).at(lepSubSubLead_index);
      thirdLepPt               = (*lepPt).at(lepSubSubLead_index);
      thirdLepEta              = (*lepEta).at(lepSubSubLead_index);
      thirdLepPhi              = (*lepPhi).at(lepSubSubLead_index);
      eventWeight              = evt.weight;
      eventWeightSyst          = evt.weight_syst_err;
      eventMT                  = evt.getMT(evt.lep0_index);
      eventMll                 = evt.getMll(evt.lep1_index,evt.lep2_index);
      eventMet                 = evt.met_Et; 
      eventNJets               = evt.nJets; 

      outtree->Fill();

      if(writeHistFitterTree){

        const double GeV = 1000;

        // Leptons
        HF_lept1Pt = (*lepPt).at(lepLead_index) * GeV;
        HF_lept2Pt = (*lepPt).at(lepSubLead_index) * GeV;
        HF_lept3Pt = (*lepPt).at(lepSubSubLead_index) * GeV;
        HF_lept1Eta = (*lepEta).at(lepLead_index);
        HF_lept2Eta = (*lepEta).at(lepSubLead_index);
        HF_lept3Eta = (*lepEta).at(lepSubSubLead_index);
        HF_lept1Phi = (*lepPhi).at(lepLead_index);
        HF_lept2Phi = (*lepPhi).at(lepSubLead_index);
        HF_lept3Phi = (*lepPhi).at(lepSubSubLead_index);
        HF_lept1q = (*lepCharge).at(lepLead_index);
        HF_lept2q = (*lepCharge).at(lepSubLead_index);
        HF_lept3q = (*lepCharge).at(lepSubSubLead_index);

        // Note the -1 here because my enum for lepFlavor is 
        // 0 for unknown, 1 for el, 2 for mu, but the HF tree
        // expects 0 for el, 1 for mu.
        HF_lept1Flav = (*lepFlavor).at(lepLead_index) - 1;
        HF_lept2Flav = (*lepFlavor).at(lepSubLead_index) - 1;
        HF_lept3Flav = (*lepFlavor).at(lepSubSubLead_index) - 1;

        // Jets
        HF_jet1Pt = evt.leadJetPt * GeV;
        HF_njets = evt.nJets;
        HF_nbjets = evt.nBJets;
        HF_nCentralLightJets = evt.nJets - evt.nBJets;

        // MET
        HF_MET = evt.met_Et * GeV;

        // Others
        HF_nSigLep = ((int) (*lepSignal).at(lepLead_index)) + ((int) (*lepSignal).at(lepSubLead_index)) + ((int) (*lepSignal).at(lepSubSubLead_index));
        // weights, scaled for 1fb-1!
        HF_eventweight = evt.weight * (1/36.075);
        HF_syst_FF = evt.weight_syst_err * (1/36.075);

        // For antiID MC that we use for subtracting off the non-Zjet bkgs
        // we want to set the weight negative (so we can hadd these trees)
        // and the syst_FF to zero (since it should be tiny anyway)
        if(subregion != "ttt" && sample != "data"){
          HF_eventweight = -1 * HF_eventweight;
          HF_syst_FF = -1 * HF_syst_FF;
        }

        HF_runNumber = *RunNumber;

        // three-lepton specific
        TLorentzVector tlv_sum = evt.getLeptonTLV(lepLead_index) + evt.getLeptonTLV(lepSubLead_index) + evt.getLeptonTLV(lepSubSubLead_index);
        HF_leptPtSumVec = tlv_sum.Pt() * GeV;
        HF_Mlll = tlv_sum.M() * GeV;

        // based on min mT
        HF_L3Mll = evt.getMll(evt.lep1_index,evt.lep2_index) * GeV; // based on the W,Z lep pairing scheme
        HF_L3Mt = evt.getMT(evt.lep0_index) * GeV;

        double traditionalMll = HF_L3Mll;
        double traditionalMt  = HF_L3Mt;
        double deltaMll = fabs(traditionalMll - 91188);

        // try the first additional potential pair
        if(evt.isSFOS(evt.lep0_index,evt.lep1_index)){
          double tmpMll = evt.getMll(evt.lep0_index,evt.lep1_index) * GeV;
          if(fabs(tmpMll - 91188) < deltaMll){
            traditionalMll = tmpMll;
            deltaMll = fabs(tmpMll - 91188);
            traditionalMt = evt.getMT(evt.lep2_index) * GeV;
          }
        }

        // now the second additional potential pair
        if(evt.isSFOS(evt.lep0_index,evt.lep2_index)){
          double tmpMll = evt.getMll(evt.lep0_index,evt.lep2_index) * GeV;
          if(fabs(tmpMll - 91188) < deltaMll){
            traditionalMll = tmpMll;
            deltaMll = fabs(tmpMll - 91188);
            traditionalMt = evt.getMT(evt.lep1_index) * GeV;
          }
        }

        // alternative pairing
        HF_bestZcandidate = traditionalMll; // note the units are already right
        HF_Mt = traditionalMt; // note the units are already right

        HF_Tree->Fill();

      }
    }
  }

  std::cout << "nEventsProcessed is " << nEventsProcessed << " and nEventsPassed is " << nEventsPassed << std::endl;

  outfile.Write();
  outfile.Close();

  if(HF_File){
    HF_File->Write();
    HF_File->Close();
  }

  if(!ffCalcTool->finalize().isSuccess()) return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
