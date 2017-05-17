#ifndef YT_SKIM_H
#define YT_SKIM_H

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "ytNUHM2Analysis/Leptons.h"
#include "ytNUHM2Analysis/Jet.h"
#include "ytNUHM2Analysis/yt_useful_functions.h"

#include <iostream>
#include <vector>
#include <string>
using namespace std;

struct ZTandP_pair
{
    int     event_number;
    // Tag leptons
    int     tag_index;
    int     tag_flavor;
    int     tag_charge;
    double  tag_pt;
    double  tag_eta;
    double  tag_phi;
    bool    tag_is_trigger_matched;
    double  tag_id_SFw;
    double  tag_iso_SFw;
    double  tag_trigger_SFw;
    // Probe leptons
    int     probe_index;
    int     probe_flavor;
    int     probe_charge;
    double  probe_pt;
    double  probe_eta;
    double  probe_phi;
    bool    probe_is_trigger_matched;
    double  probe_id_SFw;
    double  probe_iso_SFw;
    double  probe_trigger_SFw;
    bool    probe_is_baseline;
    bool    probe_is_signal;
    //double  probe_deltaR_closest_jet;
    // Tag and Probe pair
    double  baseline_mll;
    double  signal_mll;
};

class yt_skim
{
public:
    // flags
    bool isMC;
    bool isData;

    // user defined variables
    string process;

    //double derivation_stat_weights;
    long event_number;

    float luminosity;
    float cross_section;
    float k_factor;
    float filter_efficiency;
    float cross_section_kfactor_efficiency;

    double event_weight_sum;
    double event_weight;
    double pileup_weight;
    double baseline_weight;
    double signal_weight;
    int    run_number;
    float  new_AvgMu;
    float  tag_pt_threshold;
    //double Etmiss;

    // Declare the output
    TFile *output_file;
    TTree *output_tree;

public:

    vector<Electron>    vec_elec;
    vector<Muon>        vec_muon;
    vector<Lepton>      vec_lept;
    vector<Jet>         vec_jets;

    vector<Electron>    vec_baseline_elec;
    vector<Muon>        vec_baseline_muon;
    vector<Lepton>      vec_baseline_lept;
    vector<Jet>         vec_baseline_jets;

    vector<Electron>    vec_signal_elec;
    vector<Muon>        vec_signal_muon;
    vector<Lepton>      vec_signal_lept;
    vector<Jet>         vec_signal_jets;

    vector<Jet>         vec_signal_jets_no_eta_cut;
    vector<Jet>         vec_signal_bjet;

public:
    // Declaring the new variables
    // Physics object variables
    vector<ZTandP_pair> *EE_pair;
    vector<ZTandP_pair> *MuMu_pair;

    vector<bool>    *El_isBaseline;
    vector<bool>    *El_isSignal;
    vector<bool>    *El_isZTag;
    vector<bool>    *El_isZProbe;
    vector<int>     *El_Tag_index;
    vector<int>     *El_Probe_index;
    vector<double>  *El_ZTandP_mll;
    vector<bool>    *El_isZProbe_TagTriggerMatched;
    vector<bool>    *El_Tag_Trigger_Matched;
    vector<bool>    *El_Probe_Trigger_Matched;
    vector<double>  *El_DR_closest_Jet;

    vector<bool>    *Mu_isBaseline;
    vector<bool>    *Mu_isSignal;
    vector<bool>    *Mu_isZTag;
    vector<bool>    *Mu_isZProbe;
    vector<int>     *Mu_Tag_index;
    vector<int>     *Mu_Probe_index;
    vector<double>  *Mu_ZTandP_mll;
    vector<bool>    *Mu_isZProbe_TagTriggerMatched;
    vector<bool>    *Mu_Tag_Trigger_Matched;
    vector<bool>    *Mu_Probe_Trigger_Matched;
    vector<double>  *Mu_DR_closest_Jet;

    //vector<bool>    *Jet_isBaseline;
    //vector<bool>    *Jet_isSignal;
    //vector<bool>    *bJet_isSignal;

    // Events variables
    double normalization;
    double baseline_mll; // if probe is baseline lepton
    double signal_mll; // if probe is signal lepton
    //double jets_mll;
    //double Ht;
    //double meff;

public:
    // Counter
    int N_at_least_two_baseline_lepton_events;
    int N_exactly_two_baseline_lepton_events;
    int N_tot_baseline_electrons;
    int N_tot_baseline_muons;
    int N_tot_TandP_pair_electron;
    int N_tot_TandP_pair_muon;

public:
    yt_skim();
    virtual ~yt_skim() {}

    void set_isMC(bool b) { isMC = b; }
    void set_isData(bool b) { isData = b; }
    void set_process(string s) { process = s; }
    void set_event_number(long l) {event_number = l;}
    void set_luminosity(float f) { luminosity = f; }
    void set_cross_section(float f) { cross_section = f; }
    void set_k_factor(float f) { k_factor = f; }
    void set_filter_efficiency(float f) { filter_efficiency = f; }
    void set_cross_section_kfactor_efficiency(float f) { cross_section_kfactor_efficiency = f; }
    void set_event_weight_sum(double d) { event_weight_sum = d; }
    void set_event_weight(double d) { event_weight = d; }
    void set_pileup_weight(double d) { pileup_weight = d; }
    void set_baseline_weight(double d) { baseline_weight = d; }
    void set_signal_weight(double d) { signal_weight = d; }
    void set_run_number(int i) { run_number = i; }
    void set_new_AvgMu(float f) { new_AvgMu = f; }
    void set_tag_pt_threshold(float f) { tag_pt_threshold = f; }
    //void set_Etmiss(double d) { Etmiss = d; }

    void clear_all_object_vectors();
    void fill_all_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets,
                                 vector<Electron> baseline_elec, vector<Muon> baseline_muon, vector<Lepton> baseline_lept, vector<Jet> baseline_jets,
                                 vector<Electron> signal_elec, vector<Muon> signal_muon, vector<Lepton> signal_lept, vector<Jet> signal_jets);
    void copy_raw_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets);
    void copy_baseline_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets);
    void copy_signal_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets);
    void fill_signal_jets_no_eta_cut(vector<Jet> jets);
    void fill_signal_bjet(vector<Jet> signal_jets);

    void initialize(TTree *tree);
    void execute(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets,
                 vector<Electron> baseline_elec, vector<Muon> baseline_muon, vector<Lepton> baseline_lept, vector<Jet> baseline_jets,
                 vector<Electron> signal_elec, vector<Muon> signal_muon, vector<Lepton> signal_lept, vector<Jet> signal_jets);
    void finalize();

    void reset_vectors();
    void calculate_new_variables();
    void Z_tag_and_probe();

    double calculate_deltaR(Lepton lept, vector<Jet> jets);

    void debug_print();

    ClassDef(yt_skim, 1);
};

#endif // #define YT_SKIM_H
