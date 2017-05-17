#!/usr/bin/python
import os
os.system(". ../scripts/run_all_skim.sh")
os.system("./merge_data.py -b --skimmed")
os.system("./make_AvgMu_and_mll.py -b --skimmed")
os.system(". ../scripts/run_all_RLE.sh")
os.system("./merge_data.py -b --RLE")
os.system("../make_AvgMu_and_mll.py -b --RLE")
os.system("./merge_data.py -b --RLE-trigger tag_trigger_matched")
os.system("./merge_data.py -b --RLE-trigger dilepton_trigger")
os.system("./merge_data.py -b --RLE-trigger dilepton_trigger_tag_trigger_matched")
os.system("./run_making_plots.py")
#os.system("./run_background_subtraction.py -b --RLE")
os.system("./run_background_subtraction.py -b --RLE-trigger tag_trigger_matched")
# os.system("./run_background_subtraction.py -b --RLE-trigger tag_trigger_matched")
# os.system("./run_background_subtraction.py -b --RLE-trigger dilepton_trigger")
# os.system("./run_background_subtraction.py -b --RLE-trigger dilepton_trigger_tag_trigger_matched")
os.system("./run_elec_systematics.py")
os.system("./run_muon_systematics.py")
os.system("./run_trigger_systematics.py")
os.system("./run_trigger_bkg_subtraction.py")
os.system("./run_relative_differences_of_efficiency.py")
os.system("./make_matrix_method_input.py")
os.system("./make_final_RLE_plots.py")