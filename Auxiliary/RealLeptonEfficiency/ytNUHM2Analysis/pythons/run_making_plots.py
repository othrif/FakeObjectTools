#!/usr/bin/python
import os, ROOT
os.sys.argv.append( '-b' ) #Batch mode, same as python <somescript>.py -b

#cwd = os.getcwd()
#print cwd
print "NOTE!!! Fix the input path before running."
ROOT.gROOT.LoadMacro("../scripts/ytPlots.C")
#ROOT.gROOT.LoadMacro("/Users/ytshen/Documents/Working/OU/HEP/my_codes/SS_v47/ytNUHM2Analysis/scripts/ytPlots.C")

#RLE_data_elec = "/raid05/users/shen/Ximo_ntuples/v47/Cutflow/Results/20170112/hist-RLE-merged-data-elec.root"
#RLE_data_muon = "/raid05/users/shen/Ximo_ntuples/v47/Cutflow/Results/20170112/hist-RLE-merged-data-muon.root"
RLE_data_elec = "../../Results/20170222/hist-RLE-merged-data-elec.root"
RLE_data_muon = "../../Results/20170222/hist-RLE-merged-data-muon.root"

ROOT.gROOT.ProcessLine("yt_mll_plots()")
ROOT.gROOT.ProcessLine("yt_baseline_mll_plots()")
ROOT.gROOT.ProcessLine("yt_signal_mll_plots()")
# without normalizing MC
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(1,1,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(2,2,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(3,3,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(4,4,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(5,5,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(6,6,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(7,7,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(8,8,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(9,9,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(10,10,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(11,11,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(12,12,true,false)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(13,13,true,false)")
# without normalizing MC
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(1,1,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(2,2,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(3,3,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(4,4,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(5,5,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(6,6,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(7,7,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(8,8,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(9,9,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(10,10,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(11,11,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(12,12,true,false)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(13,13,true,false)")
# with normalizing MC
ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(1,1,true,true)")
ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(2,2,true,true)")
ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(3,3,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(4,4,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(5,5,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(6,6,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(7,7,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(8,8,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(9,9,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(10,10,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(11,11,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(12,12,true,true)")
# ROOT.gROOT.ProcessLine("yt_baseline_mll_ratio_plots(13,13,true,true)")
# with normalizing MC
ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(1,1,true,true)")
ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(2,2,true,true)")
ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(3,3,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(4,4,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(5,5,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(6,6,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(7,7,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(8,8,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(9,9,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(10,10,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(11,11,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(12,12,true,true)")
# ROOT.gROOT.ProcessLine("yt_signal_mll_ratio_plots(13,13,true,true)")

#ROOT.gROOT.ProcessLine("yt_make_cut_efficiency_plot(\"RLE_data_elec\",\"electron\")")
#ROOT.gROOT.ProcessLine("yt_make_cut_efficiency_plot(\"RLE_data_muon\",\"muon\")")
ROOT.gROOT.ProcessLine("yt_make_cut_efficiency_plot_2(\"../../Results/20170222/hist-RLE-merged-data-elec.root\",\"electron\")")
ROOT.gROOT.ProcessLine("yt_make_cut_efficiency_plot_2(\"../../Results/20170222/hist-RLE-merged-data-muon.root\",\"muon\")")
#ROOT.gROOT.ProcessLine("yt_make_electron_real_efficiency_plot(\"RLE_data_elec\")")
#ROOT.gROOT.ProcessLine("yt_make_muon_real_efficiency_plot(\"RLE_data_muon\")")
ROOT.gROOT.ProcessLine("yt_make_electron_real_efficiency_plot()")
ROOT.gROOT.ProcessLine("yt_make_muon_real_efficiency_plot()")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_plots(true,true,true)")
ROOT.gROOT.ProcessLine("yt_make_real_efficiency_plots(true,true,false)")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_ratio_plots(\"electron\", \"pt\",true,true,true)")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_ratio_plots(\"electron\", \"eta\",true,true,true)")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_ratio_plots(\"electron\", \"dRjet\",true,true,true)")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_ratio_plots(\"muon\", \"pt\",true,true,true)")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_ratio_plots(\"muon\", \"eta\",true,true,true)")
# ROOT.gROOT.ProcessLine("yt_make_real_efficiency_ratio_plots(\"muon\", \"dRjet\",true,true,true)")

ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"electron\", 1, 3)") # 0 < |eta| < 0.8
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"electron\", 4, 5)") # 0.8 < |eta| < 1.37
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"electron\", 7, 8)") # 1.52 < |eta| < 2
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"electron\", 1, 8)") # 0 < |eta| < 2

ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"muon\", 1, 1)") # 0 < |eta| < 0.6
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"muon\", 2, 2)") # 0.6 < |eta| < 1.2
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"muon\", 3, 3)") # 1.2 < |eta| < 1.8
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"muon\", 4, 4)") # 1.8 < |eta| < 2.5
ROOT.gROOT.ProcessLine("yt_real_efficiency_with_trigger_eta_range(true, \"muon\", 1, 4)") # 0 < |eta| < 2.5

ROOT.gROOT.ProcessLine("yt_real_efficiency_function_of_eta_with_different_pt()")

# ROOT.gROOT.ProcessLine("yt_real_efficiency_vs_AvgMu()")
# ROOT.gROOT.ProcessLine("yt_real_efficiency_vs_dR_electron()")

ROOT.gROOT.ProcessLine("yt_truth_match_TandP_comparison(\"electron\", \"pt\",true)")
ROOT.gROOT.ProcessLine("yt_truth_match_TandP_comparison(\"electron\", \"eta\",true)")
ROOT.gROOT.ProcessLine("yt_truth_match_TandP_comparison(\"electron\", \"dRjet\",true)")
ROOT.gROOT.ProcessLine("yt_truth_match_TandP_comparison(\"muon\", \"pt\",true)")
ROOT.gROOT.ProcessLine("yt_truth_match_TandP_comparison(\"muon\", \"eta\",true)")
ROOT.gROOT.ProcessLine("yt_truth_match_TandP_comparison(\"muon\", \"dRjet\",true)")

# ROOT.gROOT.ProcessLine("yt_Gtt_electron_real_efficiency_individual_cut_study()")
# ROOT.gROOT.ProcessLine("yt_Gtt_electron_real_efficiency_combine_cuts_study()")
# ROOT.gROOT.ProcessLine("yt_Gtt_electron_real_efficiency_comparison()")

ROOT.gROOT.ProcessLine("yt_background_subtraction_illustration()")
ROOT.gROOT.ProcessLine("yt_background_template_mll_plot(1,1,1,8)")
ROOT.gROOT.ProcessLine("yt_background_template_mll_plot(2,2,1,8)")
ROOT.gROOT.ProcessLine("yt_kinematics_distribution()")
ROOT.gROOT.ProcessLine("yt_deltaR_and_NJets_distributions()")
ROOT.gROOT.ProcessLine("yt_deltaR_baseline_and_signal_distributions()")
