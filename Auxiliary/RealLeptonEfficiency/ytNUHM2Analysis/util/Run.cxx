// Usage:
// Run cutflow/skim/optimization isData sample=abc [PROOF/Condor]
// Run cutflow/skim/optimization isMC sample=abc [PROOF/Condor]
// Run cutflow/skim/optimization isMC 4topSM/Zee/Zmumu/ttbar/GG_ttn1 [PROOF/Condor]
//

#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/CondorDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>
#include <TH1.h>

#include "ytNUHM2Analysis/ytEventSelection.h"

#include <iostream>
#include <string>
using namespace std;

int main( int argc, char* argv[] ) {

    // Take the submit directory from the input if provided:
    //std::string submitDir = "submitDir";
    std::string submitDir;
    //if( argc > 1 ) submitDir = argv[ 1 ];

    bool isMC = false;
    bool isData = false;

    bool isCutflow = false;
    bool isSkim = false;
    bool isOptimization = false;

    string process;

    bool use_Condor = false;
    bool use_Grid = false;
    bool use_PROOF = false;

    for (int i = 1; i < argc; i++) {
        //const char *key = argv[i];
        const char* key = strtok(argv[i], "=");
        const char* val = strtok(0, " ");
        // Check MC or data.
        if (strcmp(key, "isMC") == 0)
            isMC = true;
        else if (strcmp(key, "isData") == 0)
            isData = true;
        // Check cutflow
        else if (strcmp(key, "cutflow") == 0)
            isCutflow = true;
        // Check skim
        else if (strcmp(key, "skim") == 0)
            isSkim = true;
        // Check SR optimization
        else if (strcmp(key, "optimization") == 0)
            isOptimization = true;
        //
        // Choose samples to run.
        // This is used to run over the data samples because data samples have different period.
        else if (strcmp(key, "sample") == 0)
            process = val;
        //
        //
        // We use 4topSM to compare cutflow.
        else if (strcmp(key, "4topSM") == 0)
            process = "4topSM";
        // Zee, Zmumu, ttbar, GG_ttn1 are used to study real lepton efficiency.
        else if (strcmp(key, "Zee") == 0)
            process = "Zee";
        else if (strcmp(key, "Zmumu") == 0)
            process = "Zmumu";
        else if (strcmp(key, "ttbar") == 0)
            process = "ttbar";
        else if (strcmp(key, "GG_ttn1") == 0)
            process = "GG_ttn1";
        //
        // signal NUHM2
        // strong
        else if (strcmp(key, "NUHM2_m12_300_strong") == 0)
            process = "NUHM2_m12_300_strong";
        else if (strcmp(key, "NUHM2_m12_350_strong") == 0)
            process = "NUHM2_m12_350_strong";
        else if (strcmp(key, "NUHM2_m12_400_strong") == 0)
            process = "NUHM2_m12_400_strong";
        else if (strcmp(key, "NUHM2_m12_500_strong") == 0)
            process = "NUHM2_m12_500_strong";
        else if (strcmp(key, "NUHM2_m12_600_strong") == 0)
            process = "NUHM2_m12_600_strong";
        else if (strcmp(key, "NUHM2_m12_700_strong") == 0)
            process = "NUHM2_m12_700_strong";
        else if (strcmp(key, "NUHM2_m12_800_strong") == 0)
            process = "NUHM2_m12_800_strong";
        // weak
        else if (strcmp(key, "NUHM2_m12_300_weak") == 0)
            process = "NUHM2_m12_300_weak";
        else if (strcmp(key, "NUHM2_m12_350_weak") == 0)
            process = "NUHM2_m12_350_weak";
        else if (strcmp(key, "NUHM2_m12_400_weak") == 0)
            process = "NUHM2_m12_400_weak";
        else if (strcmp(key, "NUHM2_m12_500_weak") == 0)
            process = "NUHM2_m12_500_weak";
        else if (strcmp(key, "NUHM2_m12_600_weak") == 0)
            process = "NUHM2_m12_600_weak";
        else if (strcmp(key, "NUHM2_m12_700_weak") == 0)
            process = "NUHM2_m12_700_weak";
        else if (strcmp(key, "NUHM2_m12_800_weak") == 0)
            process = "NUHM2_m12_800_weak";
        //
        //
        // bkg samples
        else if (strcmp(key, "ttW") == 0)
            process = "ttW";
        else if (strcmp(key, "ttee") == 0)
            process = "ttee";
        else if (strcmp(key, "ttmumu") == 0)
            process = "ttmumu";
        else if (strcmp(key, "tttautau") == 0)
            process = "tttautau";
        else if (strcmp(key, "llvvjj_ss_EW4") == 0)
            process = "llvvjj_ss_EW4";
        else if (strcmp(key, "llvvjj_ss_EW6") == 0)
            process = "llvvjj_ss_EW6";
        else if (strcmp(key, "lllvjj_EW6") == 0)
            process = "lllvjj_EW6";
        else if (strcmp(key, "lllljj_EW6") == 0)
            process = "lllljj_EW6";
        else if (strcmp(key, "ggllll") == 0)
            process = "ggllll";
        else if (strcmp(key, "NNPDF30NNLO_llll") == 0)
            process = "llll";
        else if (strcmp(key, "NNPDF30NNLO_lllv") == 0)
            process = "lllv";
        else if (strcmp(key, "3top_SM") == 0)
            process = "3top_SM";
        else if (strcmp(key, "WH125_inc") == 0)
            process = "WH125_inc";
        else if (strcmp(key, "ZH125_inc") == 0)
            process = "ZH125_inc";
        else if (strcmp(key, "ttH125_di") == 0)
            process = "ttH125_di";
        else if (strcmp(key, "ttH125_se") == 0)
            process = "ttH125_se";
        else if (strcmp(key, "ttH125_al") == 0)
            process = "ttH125_al";
        else if (strcmp(key, "WWW_3l3v") == 0)
            process = "WWW_3l3v";
        else if (strcmp(key, "WWZ_4l2v") == 0)
            process = "WWZ_4l2v";
        else if (strcmp(key, "WWZ_2l4v") == 0)
            process = "WWZ_2l4v";
        else if (strcmp(key, "WZZ_5l1v") == 0)
            process = "WZZ_5l1v";
        else if (strcmp(key, "WZZ_3l3v") == 0)
            process = "WZZ_3l3v";
        else if (strcmp(key, "ZZZ_6l0v") == 0)
            process = "ZZZ_6l0v";
        else if (strcmp(key, "ZZZ_4l2v") == 0)
            process = "ZZZ_4l2v";
        else if (strcmp(key, "ZZZ_2l4v") == 0)
            process = "ZZZ_2l4v";
        else if (strcmp(key, "tZ_4fl_tchan_noAllHad") == 0)
            process = "tZ_4fl_tchan_noAllHad";
        else if (strcmp(key, "ttbarWW") == 0)
            process = "ttbarWW";
        else if (strcmp(key, "tWZDR") == 0)
            process = "tWZDR";
        //
        //
        // Specify the driver to run.
        else if (strcmp(key, "Condor") == 0)
            use_Condor = true;
        else if (strcmp(key, "Grid") == 0)
            use_Grid = true;
        else if (strcmp(key, "PROOF") == 0)
            use_PROOF = true;
    }

    // Print out input arguments
    if (isCutflow)
        printf("Running cutflow for %s\n", isMC ? process.c_str() : "Data");
    else if (isSkim)
        printf("Running skim for %s\n", isMC ? process.c_str() : "Data");
    else if (isOptimization)
        printf("Running SR optimization for %s\n", isMC ? process.c_str() : "Data");
    cout << "process = " << process << endl;

    printf("isMC = %s, isData = %s\n", isMC ? "true" : "false", isData ? "true" : "false");

    if (isMC && !process.empty()) {
        if (isCutflow)
            submitDir = "cutflow_MC_" + process;
        else if (isSkim)
            submitDir = "skim_MC_" + process;
        else if (isOptimization)
            submitDir = "optimization_MC_" + process;
    }
    else if (isData && !process.empty()) {
        if (isCutflow)
            submitDir = "cutflow_Data_" + process;
        else if (isSkim)
            submitDir = "skim_Data_" + process;
        else if (isOptimization)
            submitDir = "optimization_Data_" + process;
    }
    cout << "submitDir = " << submitDir << endl;

    if (use_Condor) {
        printf("Submit jobs to CondorDriver...\n");
    }
    else if (use_Grid) {
        printf("Submit jobs to PrunDriver...\n");
    }
    else if (use_PROOF) {
        printf("Submit jobs to ProofDriver...\n");
    }
    else {
        printf("Submit jobs to DirectDriver...\n");
    }

    // Construct the samples to run on:
    SH::SampleHandler sh;

    // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
    const char* inputFilePath;
    const char* NUHM2_inputFilePath;

    if (isMC) {
        cout << "Read MC files..." << endl;
        inputFilePath = "../MC"; // no slash (/) at the end.
        NUHM2_inputFilePath = "/raid05/atlas/data/NUHM2/strongMC/p2666"; // no slash (/) at the end.

        //
        // For cutflow study
        if (process == "4topSM") {
            // SH::ScanDir().filePattern("4topSM.root").scan(sh, inputFilePath);
            // SH::ScanDir().samplePattern("user.*.410080.*4topSM*").scan(sh, inputFilePath); // Get all root files in this dataset
            SH::ScanDir().filePattern("4topSM_merged.root").scan(sh, inputFilePath);
        }
        //
        // For real lepton efficiency study
        else if (process == "Zee") {
            SH::ScanDir().filePattern("Zee_merged.root").scan(sh, inputFilePath); // Get specific root file
            // SH::ScanDir().samplePattern("user.*.361106.*Zee*").scan(sh, inputFilePath); // Get all root files in this dataset
        }
        else if (process == "Zmumu") {
            SH::ScanDir().filePattern("Zmumu_merged.root").scan(sh, inputFilePath); // Get specific root file
            // SH::ScanDir().samplePattern("user.*.361107.*Zmumu*").scan(sh, inputFilePath); // Get all root files in this dataset
        }
        else if (process == "ttbar") {
            SH::ScanDir().filePattern("ttbar_merged.root").scan(sh, inputFilePath); // Get specific root file
            // SH::ScanDir().samplePattern("user.*.410000.*ttbar*nonallhad*").scan(sh, inputFilePath); // Get all root files in this dataset
        }
        else if (process == "GG_ttn1") {
            SH::ScanDir().filePattern("GG_ttn1_merged.root").scan(sh, inputFilePath); // Get specific root file
            //SH::ScanDir().samplePattern("user.*.*GG_ttn1_*.root").scan(sh, inputFilePath); // Get all root files in this dataset
        }
        //
        // For signal optimization
        // strong
        else if (process == "NUHM2_m12_300_strong") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370600.MGPy8EG_A14N23LO_NUHM2_m12_300_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_350_strong") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370601.MGPy8EG_A14N23LO_NUHM2_m12_350_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_400_strong") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370602.MGPy8EG_A14N23LO_NUHM2_m12_400_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_500_strong") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370603.MGPy8EG_A14N23LO_NUHM2_m12_500_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_600_strong") {
            SH::ScanDir().samplePattern("user.yushen.Jul05.v39.370604.MGPy8EG_A14N23LO_NUHM2_m12_600_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_700_strong") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370605.MGPy8EG_A14N23LO_NUHM2_m12_700_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_800_strong") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370606.MGPy8EG_A14N23LO_NUHM2_m12_800_strong.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        // weak
        else if (process == "NUHM2_m12_300_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370617.MGPy8EG_A14N23LO_NUHM2_m12_300_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_350_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370618.MGPy8EG_A14N23LO_NUHM2_m12_350_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_400_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370619.MGPy8EG_A14N23LO_NUHM2_m12_400_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_500_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370620.MGPy8EG_A14N23LO_NUHM2_m12_500_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_600_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370621.MGPy8EG_A14N23LO_NUHM2_m12_600_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_700_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370622.MGPy8EG_A14N23LO_NUHM2_m12_700_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        else if (process == "NUHM2_m12_800_weak") {
            SH::ScanDir().samplePattern("user.yushen.Aug04.v44.370623.MGPy8EG_A14N23LO_NUHM2_m12_800_weak.DAOD_SUSY2.a766_r7676_output.root").scan(sh, NUHM2_inputFilePath); // Get all root files in this dataset
        }
        //
        // bkg
        else if (process == "ttW") {
            SH::ScanDir().filePattern("ttW_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ttee") {
            SH::ScanDir().filePattern("ttee_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ttmumu") {
            SH::ScanDir().filePattern("ttmumu_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "tttautau") {
            SH::ScanDir().filePattern("tttautau_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "llvvjj_ss_EW4") {
            SH::ScanDir().filePattern("llvvjj_ss_EW4_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "llvvjj_ss_EW6") {
            SH::ScanDir().filePattern("llvvjj_ss_EW6_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "lllvjj_EW6") {
            SH::ScanDir().filePattern("lllvjj_EW6_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "lllljj_EW6") {
            SH::ScanDir().filePattern("lllljj_EW6_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ggllll") {
            SH::ScanDir().filePattern("ggllll_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "llll") {
            SH::ScanDir().filePattern("llll_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "lllv") {
            SH::ScanDir().filePattern("lllv_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "3top_SM") {
            SH::ScanDir().filePattern("3top_SM_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "WH125_inc") {
            SH::ScanDir().filePattern("WH125_inc_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ZH125_inc") {
            SH::ScanDir().filePattern("ZH125_inc_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ttH125_di") {
            SH::ScanDir().filePattern("ttH125_dilep_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ttH125_se") {
            SH::ScanDir().filePattern("ttH125_semilep_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ttH125_al") {
            SH::ScanDir().filePattern("ttH125_allhad_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "WWW_3l3v") {
            SH::ScanDir().filePattern("WWW_3l3v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "WWZ_4l2v") {
            SH::ScanDir().filePattern("WWZ_4l2v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "WWZ_2l4v") {
            SH::ScanDir().filePattern("WWZ_2l4v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "WZZ_5l1v") {
            SH::ScanDir().filePattern("WZZ_5l1v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "WZZ_3l3v") {
            SH::ScanDir().filePattern("WZZ_3l3v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ZZZ_6l0v") {
            SH::ScanDir().filePattern("ZZZ_6l0v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ZZZ_4l2v") {
            SH::ScanDir().filePattern("ZZZ_4l2v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ZZZ_2l4v") {
            SH::ScanDir().filePattern("ZZZ_2l4v_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "tZ_4fl_tchan_noAllHad") {
            SH::ScanDir().filePattern("tZ_4fl_tchan_noAllHad_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "ttbarWW") {
            SH::ScanDir().filePattern("ttbarWW_merged.root").scan(sh, inputFilePath);
        }
        else if (process == "tWZDR") {
            SH::ScanDir().filePattern("tWZDR_merged.root").scan(sh, inputFilePath);
        }
        // sample
        else if (!process.empty()) {
            SH::ScanDir().samplePattern("user.*_" + process + "*_output.root").scan(sh, inputFilePath); // Get all root files in this dataset
        }
    }
    else if (isData) {
        cout << "Read Data files..." << endl;
        inputFilePath = "../Data"; // no slash (/) at the end.
        //SH::ScanDir().scan(sh, inputFilePath); // Get all datasets in inputFilePath
        SH::ScanDir().filePattern(process + ".root").scan(sh, inputFilePath); // Get specific root file
        //SH::ScanDir().filePattern("merged_all_data.root").scan(sh, inputFilePath); // Get specific root file
        //SH::ScanDir().samplePattern("user.*.physics_Main.DAOD_SUSY2.*").scan(sh, inputFilePath); // Get all root files in this dataset
    }

    // Set the name of the input TTree.
    sh.setMetaString("nc_tree", "AnaNtup");

    // Print what we found:
    sh.print(); // list all the root files in the dataset
    // cout << "sh.size()=" << sh.size() << endl; // number of dataset

    //
    // Get DerivationStat_Weights from input files
    //
    double derivation_stat_weights = 0;
    SH::Sample *sample = sh.at(sh.size() - 1); // Because we only have one dataset at here.
    // The above line will cause problem when we have several sample loaded in to the SH.

    // cout << "sample name=" << sample->name() << endl; // dataset name
    // cout << "numFiles()=" << sample->numFiles() << endl; // number of root files in dataset
    for (unsigned int i = 0; i < sample->numFiles() ; i++) {
        string fileName = sample->fileName(i); // this returns file://root file name and path
        string remove = "file://"; // need to remove "file://" part
        string::size_type find_remove_part = fileName.find(remove);
        if (find_remove_part != string::npos)
            fileName.erase(find_remove_part, remove.length()); // now contains root file name and path only
        //cout << "fileName(" << i << ")=" << fileName << endl;
        TFile *file = new TFile(fileName.c_str());
        TH1D *h1 = (TH1D *)file->Get("DerivationStat_Weights");
        derivation_stat_weights += h1->GetBinContent(1);
    }
    // cout << "DerivationStat_Weights=" << derivation_stat_weights << endl;

    // Get the dataset name
    string dataset_name = sh.at(0)->fileName(0);
    bool isFullSim = false, isAF2Sim = false;
    if (dataset_name.find("_r") != string::npos) {
        isFullSim = true;
        isAF2Sim = false;
    }
    if (dataset_name.find("_a") != string::npos ||
        dataset_name.find(".a") != string::npos) {
        isFullSim = false;
        isAF2Sim = true;
    }
    printf("isFullSim = %s, isAF2Sim = %s\n", isFullSim ? "true" : "false", isAF2Sim ? "true" : "false");

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    //job.options()->setDouble (EL::Job::optMaxEvents, 50);

    // Add our analysis to the job:
    ytEventSelection *alg = new ytEventSelection();

    alg->set_isMC(isMC);
    alg->set_isData(isData);
    alg->set_isSkim(isSkim);
    alg->set_isOptimization(isOptimization);
    alg->set_isFullSim(isFullSim);
    alg->set_isAF2Sim(isAF2Sim);
    alg->set_process(process);
    alg->set_tag_pt_threshold(25000.);
    alg->set_derivation_stat_weights(derivation_stat_weights);

    const double luminosity = 36.5; // unit: 1/fb 2015+2016
    alg->set_luminosity(luminosity);

    job.algsAdd( alg );

    if (use_Condor) {
        // Run the jobs using the Condor driver:
        EL::CondorDriver driver;
        // some commands for setting up root on the nodes
        driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase ; source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh ; eval rcSetup Base,2.5.1";
        driver.submitOnly( job, submitDir );
    }
    else if (use_Grid) {
        // Run the jobs using the Grid driver:
        EL::PrunDriver driver;
        // Specify how to name the grid output datasets
        // Note that %nickname% is populated when you do voms-proxy init, so this does not have to be replaced by hand
        driver.options()->setString("nc_outputSampleName", "user.%nickname%.%in:name[2]%.%in:name[3]%.%in:name[6]%.");
        driver.submitOnly( job, submitDir );
    }
    else if (use_PROOF) {
        EL::ProofDriver driver;
        driver.submit( job, submitDir );
    }
    else {
        // Run the job using the local/direct driver:
        EL::DirectDriver driver;
        driver.submit( job, submitDir );
    }

    return 0;
}
