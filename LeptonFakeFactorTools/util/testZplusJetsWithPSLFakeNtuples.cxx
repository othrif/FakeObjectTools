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

  TString configFile = "$ROOTCOREBIN/../LeptonFakeFactorTools/data/FakeFactorConfigs_WID_Electron.conf";
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
  TTreeReader reader("physics", inFile);
  READ_TREE(unsigned long long, EventNumber);
  READ_TREE(int,                nBaselineLeptons);
  READ_TREE(std::vector<int>,   lepFlavor);
  READ_TREE(std::vector<int>,   lepCharge);
  READ_TREE(std::vector<int>,   lepPassOR);
  READ_TREE(std::vector<int>,   lepTruthDetailed);
  READ_TREE(std::vector<int>,   lepMatchesTrigger);
  READ_TREE(std::vector<float>, lepPt);
  READ_TREE(std::vector<float>, lepEta);
  READ_TREE(std::vector<float>, lepPhi);
  READ_TREE(std::vector<float>, lepD0);
  READ_TREE(std::vector<float>, lepD0Significance);
  READ_TREE(std::vector<float>, lepZ0SinTheta);
  READ_TREE(std::vector<float>, lepEleEtaBE);
  READ_TREE(float,              met_Et);
  READ_TREE(float,              met_phi);
  READ_TREE(std::vector<int>,   TightLH);
  READ_TREE(std::vector<int>,   MediumLH);
  READ_TREE(std::vector<int>,   LooseLH);
  READ_TREE(std::vector<int>,   LooseAndBLayerLH);
  READ_TREE(std::vector<int>,   VeryLooseLH);
  READ_TREE(std::vector<int>,   Tight);
  READ_TREE(std::vector<int>,   Medium);
  READ_TREE(std::vector<int>,   Loose);
  READ_TREE(std::vector<int>,   VeryLoose);
  READ_TREE(std::vector<int>,   Gradient);
  READ_TREE(std::vector<int>,   GradientLoose);
  READ_TREE(std::vector<int>,   LooseTrackOnly);
  READ_TREE(int,                HLT_e24_lhmedium_L1EM18VH);
  READ_TREE(int,                HLT_e24_lhmedium_L1EM20VH);
  READ_TREE(int,                HLT_e60_lhmedium);
  READ_TREE(int,                HLT_e120_lhloose);
  READ_TREE(int,                HLT_mu20_iloose_L1MU15);
  READ_TREE(int,                HLT_mu50);
  READ_TREE(float,              TotalWeightNoSF);
  READ_TREE(std::vector<float>, lepEleFullSF_Reco_TightLLH_d0z0_v8_isolGradient);
  READ_TREE(std::vector<float>, lepEleFullSF_Reco_MediumLLH_d0z0_v8_isolGradientLoose);
  READ_TREE(std::vector<float>, lepEleFullSF_Reco_LooseAndBLayerLLH_d0z0_isolLooseTrackOnly);
  READ_TREE(std::vector<float>, lepMuFullSF_Medium_d0z0_isolGradientLoose);


  // Loop over events!
  int nEventsProcessed = 0;
  int nEventsPassed= 0;

  while(reader.Next()) {
    ++nEventsProcessed;

    // create EventDef class instance whose members are all the default values.
    // (of course, EventDef* evt = new EventDef(); would also work, if .'s replaced with ->'s)
    EventDef evt;

    // take the logical OR of the single lepton triggers
    bool passTrigger = (*HLT_e24_lhmedium_L1EM18VH && sample != "data") || (*HLT_e24_lhmedium_L1EM20VH && sample == "data") || *HLT_e60_lhmedium || *HLT_e120_lhloose || *HLT_mu20_iloose_L1MU15 || *HLT_mu50;

    // Now fill derived quantities for the tool to use
    int nLeps = (*lepPt).size(); // all of the vectors should be the same size...

    // Look at truth information, when applicable
    // Enum from pennSoftLepton is:
    //  NoTruthLabel
    //  ,LeptonW
    //  ,LeptonZ
    //  ,LeptonSUSY
    //  ,LeptonTtbar
    //  ,LeptonJpsi
    //  ,LeptonUnknown
    //  ,LeptonTau
    //  (and everything after is not a lepton,
    //  so 1-7 are all truth leptons)
    std::vector<int> lepIsTruth(nLeps, 0);

    // Look at ID and isolation to fill whether a given lepton passes
    // the signal ID, other ID, or antiID, as well as the appropriate
    // scale factors
    std::vector<int> lepPassSigID(nLeps, 0);
    std::vector<int> lepPassOtherID(nLeps, 0);
    std::vector<int> lepPassAntiID(nLeps, 0);

    for(int i = 0; i < nLeps; ++i){

      int truth_info = (*lepTruthDetailed).at(i);

      // set truth = true! (based on the enum from above)
      if(truth_info > 0 && truth_info < 8){
        lepIsTruth.at(i) = 1;
      }

      // Fill ID and SF vectors below:
      // Note this is all appropriate only for the SM WZ 2015 analysis.
      // Also, if your ntuple has a branch for isSignalLepton and
      // isAntiIDLepton, then you can completely skip the next
      // ~200 lines of code...

      if( (*lepFlavor).at(i) == 1 ){ // electrons

        double pt = (*lepPt).at(i);
        double absEtaBE = fabs( (*lepEleEtaBE).at(i) );
        double d0sig = fabs( (*lepD0Significance).at(i) );
        double z0sinTheta = fabs( (*lepZ0SinTheta).at(i) );

        bool passesSigIsolation = false;   bool passesSigID = false;
        bool passesOtherIsolation = false; bool passesOtherID = false;

        double ptminSig = 999; double ptminOther = 999; double ptminAntiID = 999;

        bool passesVeryLoose = (*VeryLooseLH).at(i);
        bool passesAllCuts = (d0sig < 5) && (*GradientLoose).at(i) && (*MediumLH).at(i);
        bool passesAntiID = passesVeryLoose && !passesAllCuts;

        // If fakeFactorType is WID or ZID, then we set the cuts for measuring
        // the fake factor for those configurations. Otherwise, we'll do the nominal
        // selection, where the Wlep has the WID criteria applied and the Zleps have
        // the ZID criteria applied. Note that all of this is only necessary because 
        // different ID and pt thresholds were used for the different leptons in the 
        // SMWZ 2015 analysis!
        if(fakeFactorType == "WID_electron"){
          ptminSig = 20; ptminOther = 15; ptminAntiID = 20;

          passesSigID = (*TightLH).at(i);
          passesSigIsolation = (*Gradient).at(i);
          passesOtherID = (*MediumLH).at(i);
          passesOtherIsolation = (*GradientLoose).at(i);
        }
        else if(fakeFactorType == "ZID_electron"){
          ptminSig = 15; ptminOther = 15; ptminAntiID = 15;

          passesSigID = (*MediumLH).at(i);
          passesSigIsolation = (*GradientLoose).at(i);
          passesOtherID = (*MediumLH).at(i);
          passesOtherIsolation = (*GradientLoose).at(i);
        }
        else if(fakeFactorType == "WID_muon"){
          ptminSig = 20; ptminOther = 15; ptminAntiID = 20;

          passesSigID = (*TightLH).at(i);
          passesSigIsolation = (*Gradient).at(i);
          passesOtherID = (*LooseAndBLayerLH).at(i);
          passesOtherIsolation = (*LooseTrackOnly).at(i);
        }
        else if(fakeFactorType == "ZID_muon"){
          ptminSig = 15; ptminOther = 15; ptminAntiID = 15;

          passesSigID = (*MediumLH).at(i);
          passesSigIsolation = (*GradientLoose).at(i);
          passesOtherID = (*LooseAndBLayerLH).at(i);
          passesOtherIsolation = (*LooseTrackOnly).at(i);
        }
        else if(fakeFactorType == "Nominal"){
          ptminSig = 20; ptminOther = 15; ptminAntiID = 15;

          passesSigID = (*TightLH).at(i);
          passesSigIsolation = (*Gradient).at(i);
          passesOtherID = (*MediumLH).at(i);
          passesOtherIsolation = (*GradientLoose).at(i);
        }
        else{
          std::cerr << "fakeFactorType " << fakeFactorType << " not defined!" << std::endl;
          return EXIT_FAILURE;
        }

        // Now apply the cuts!
        if( (pt > ptminSig) && 
            (absEtaBE < 1.37 ||1.52 < absEtaBE) && 
            (d0sig < 5) &&
            (z0sinTheta < 0.5) &&
            (passesSigIsolation) &&
            (passesSigID) &&
            (*lepPassOR).at(i) ){

          lepPassSigID.at(i) = true;
        }

        if( (pt > ptminOther) && 
            (absEtaBE < 1.37 ||1.52 < absEtaBE) && 
            (d0sig < 5) &&
            (z0sinTheta < 0.5) &&
            (passesOtherIsolation) &&
            (passesOtherID) &&
            (*lepPassOR).at(i) ){

          lepPassOtherID.at(i) = true;
        }

        if( (pt > ptminAntiID) && 
            (absEtaBE < 1.37 ||1.52 < absEtaBE) && 
            (passesAntiID) &&
            (*lepPassOR).at(i) ){

          lepPassAntiID.at(i) = true;
        }

      }
      else if( (*lepFlavor).at(i) == 2){ // muons

        double pt = (*lepPt).at(i);
        double d0sig = fabs( (*lepD0Significance).at(i) );
        double z0sinTheta = fabs( (*lepZ0SinTheta).at(i) );

        bool passesID = (*Medium).at(i);
        bool passesIsolation = (*GradientLoose).at(i);

        bool passesMedium = (*Medium).at(i);
        bool passesAllCuts = (d0sig < 3) && (z0sinTheta < 0.5) && ((*GradientLoose).at(i));

        bool passesAntiID = passesMedium && !passesAllCuts;

        double ptminSig    = 20;
        double ptminOther  = 15;
        double ptminAntiID = 20;

        if(fakeFactorType == "WID_electron"){
          ptminSig = 20; ptminOther = 15; ptminAntiID = 20;
        }
        else if(fakeFactorType == "ZID_electron"){
          ptminSig = 15; ptminOther = 15; ptminAntiID = 15;
        }
        else if(fakeFactorType == "WID_muon"){
          ptminSig = 20; ptminOther = 15; ptminAntiID = 20;
        }
        else if(fakeFactorType == "ZID_muon"){
          ptminSig = 15; ptminOther = 15; ptminAntiID = 15;
        }
        else if(fakeFactorType == "Nominal"){
          ptminSig = 20; ptminOther = 15; ptminAntiID = 15;
        }

        if( (pt > ptminSig) && 
            (d0sig < 3) &&
            (z0sinTheta < 0.5) &&
            (passesIsolation) &&
            (passesID) &&
            (*lepPassOR).at(i) ){

          lepPassSigID.at(i) = true;
        }

        if( (pt > ptminOther) && 
            (d0sig < 3) &&
            (z0sinTheta < 0.5) &&
            (passesIsolation) &&
            (passesID) &&
            (*lepPassOR).at(i) ){

          lepPassOtherID.at(i) = true;
        }

        if( (pt > ptminAntiID) && 
            (passesAntiID) ){

          lepPassAntiID.at(i) = true;
        }

      }
      else{
        std::cout << "This lepton isn't an electron or a muon! Exiting." << std::endl; 
        return EXIT_FAILURE;
      }
    }

    // fill in the struct
    evt.EventNumber               = *EventNumber;
    evt.nBaselineLeptons          = *nBaselineLeptons;
    evt.lepFlavor                 = *lepFlavor;
    evt.lepCharge                 = *lepCharge;
    evt.lepIsTruth                = lepIsTruth;
    evt.lepMatchesTrigger         = *lepMatchesTrigger;
    evt.lepPt                     = *lepPt;
    evt.lepEta                    = *lepEta;
    evt.lepPhi                    = *lepPhi;

    evt.met_Et                    = *met_Et;
    evt.met_phi                   = *met_phi;

    evt.lepPassSigID              = lepPassSigID;
    evt.lepPassOtherID            = lepPassOtherID;
    evt.lepPassAntiID             = lepPassAntiID;

    // Note! This is because this analysis did not use the dilepton triggers
    evt.passSingleLeptonTrigger   = passTrigger;
    evt.passDileptonTrigger       = false;

    evt.weight                    = *TotalWeightNoSF;

    if(ffCalcTool->accept(evt)){
      ++nEventsPassed;
    }
  }

  std::cout << "nEventsProcessed is " << nEventsProcessed << " and nEventsPassed is " << nEventsPassed << std::endl;

  if(!ffCalcTool->finalize().isSuccess()) return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
