#ifndef YT_USEFUL_FUNCTIONS_H
#define YT_USEFUL_FUNCTIONS_H

#include "ytNUHM2Analysis/Leptons.h"
#include "ytNUHM2Analysis/Jet.h"

#include <iomanip>

// Functions declarations
/*
void set_isMC(bool b) { isMC = b; }
void set_isData(bool b) { isData = b; }

void set_isSkim(bool b) { isSkim = b; }
void set_isOptimization(bool b) { isOptimization = b; }

void set_isFullSim(bool b) { isFullSim = b; }
void set_isAF2Sim(bool b) { isAF2Sim = b; }

void set_process(string s) { process = s; }
void set_luminosity(float f) { luminosity = f; }

void set_cross_section(float f) { cross_section = f; }
void set_k_factor(float f) { k_factor = f; }
void set_filter_efficiency(float f) { filter_efficiency = f; }
void set_cross_section_kfactor_efficiency(float f) { cross_section_kfactor_efficiency = f; }

void set_event_weight(double d) { event_weight = d; }
void set_pileup_weight(double d) { pileup_weight = d; }

void set_tag_pt_threshold(float f) { tag_pt_threshold = f; }
*/
//template<typename T> bool sort_descending_Pt(T obj1, T obj2); // NOTE: cannot use const ref for arguments.
//template<> bool sort_descending_Pt(TLorentzVector tlv1, TLorentzVector tlv2); // explicit specialization

double calculate_Ht(vector<Lepton> signal_leptons, vector<Jet> signal_jets);
double calculate_Meff(double Ht, double Etmiss);
double calculate_mll(vector<Lepton> leptons);
double calculate_mjj(vector<Jet> jets);

// Weight related methods
double ID_weight(vector<Electron> elec, bool signal = true);
double ID_weight(vector<Muon> muon);
double Iso_weight(vector<Electron> elec);
double Iso_weight(vector<Muon> muon);
double lept_weight(double elec_ID_weight, double elec_iso_weight, double muon_ID_weight, double muon_iso_weight);
double jets_weight(vector<Jet> jets);
double calculate_weight(double event_weight, double PRW_weight, double lept_weight, double jets_weight);
double calculate_weight(double event_weight, double PRW_weight, vector<Electron> elec, vector<Muon> muon, vector<Jet> jets);
//void update(int cut, bool passed, double weight);
//void debug_sum_of_weight_print();

float deltaR(TLorentzVector tlv1, TLorentzVector tlv2);

#endif // #ifndef YT_USEFUL_FUNCTIONS_H
