#!/usr/bin/python
import os

def main():
    do_background_subtraction()
    get_results_after_background_subtraction()

def do_background_subtraction():
    commands = []
    # 10 < pT < 15: 0 < |eta| < 0.8
    command_a01 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,1,1,1,3)\' > baseline_template_range_baseline_pT10_15_eta008_tag_trigger_matched.txt'
    command_a02 = 'root -q -b \'../scripts/ytBackground_subtraction.C("single_lepton_trigger","baseline","range_baseline",true,60,120,80,100,1,1,1,3)\' > baseline_template_range_baseline_pT10_15_eta008_single_lepton_trigger.txt'
    command_a03 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger","baseline","range_baseline",true,60,120,80,100,1,1,1,3)\' > baseline_template_range_baseline_pT10_15_eta008_dilepton_trigger.txt'
    command_a04 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger_tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,1,1,1,3)\' > baseline_template_range_baseline_pT10_15_eta008_dilepton_trigger_tag_trigger_matched.txt'
    # 10 < pT < 15: 0.8 < |eta| < 1.37
    command_a05 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,1,1,4,5)\' > baseline_template_range_baseline_pT10_15_eta08137_tag_trigger_matched.txt'
    command_a06 = 'root -q -b \'../scripts/ytBackground_subtraction.C("single_lepton_trigger","baseline","range_baseline",true,60,120,80,100,1,1,4,5)\' > baseline_template_range_baseline_pT10_15_eta08137_single_lepton_trigger.txt'
    command_a07 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger","baseline","range_baseline",true,60,120,80,100,1,1,4,5)\' > baseline_template_range_baseline_pT10_15_eta08137_dilepton_trigger.txt'
    command_a08 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger_tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,1,1,4,5)\' > baseline_template_range_baseline_pT10_15_eta08137_dilepton_trigger_tag_trigger_matched.txt'
    # 10 < pT < 15: 1.52 < |eta| < 2.01
    command_a09 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,1,1,7,8)\' > baseline_template_range_baseline_pT10_15_eta152201_tag_trigger_matched.txt'
    command_a10 = 'root -q -b \'../scripts/ytBackground_subtraction.C("single_lepton_trigger","baseline","range_baseline",true,60,120,80,100,1,1,7,8)\' > baseline_template_range_baseline_pT10_15_eta152201_single_lepton_trigger.txt'
    command_a11 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger","baseline","range_baseline",true,60,120,80,100,1,1,7,8)\' > baseline_template_range_baseline_pT10_15_eta152201_dilepton_trigger.txt'
    command_a12 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger_tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,1,1,7,8)\' > baseline_template_range_baseline_pT10_15_eta152201_dilepton_trigger_tag_trigger_matched.txt'

    # 15 < pT < 20: 0 < |eta| < 0.8
    command_b01 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,2,2,1,3)\' > baseline_template_range_baseline_pT15_20_eta008_tag_trigger_matched.txt'
    command_b02 = 'root -q -b \'../scripts/ytBackground_subtraction.C("single_lepton_trigger","baseline","range_baseline",true,60,120,80,100,2,2,1,3)\' > baseline_template_range_baseline_pT15_20_eta008_single_lepton_trigger.txt'
    command_b03 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger","baseline","range_baseline",true,60,120,80,100,2,2,1,3)\' > baseline_template_range_baseline_pT15_20_eta008_dilepton_trigger.txt'
    command_b04 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger_tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,2,2,1,3)\' > baseline_template_range_baseline_pT15_20_eta008_dilepton_trigger_tag_trigger_matched.txt'
    # 15 < pT < 20: 0.8 < |eta| < 1.37
    command_b05 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,2,2,4,5)\' > baseline_template_range_baseline_pT15_20_eta08137_tag_trigger_matched.txt'
    command_b06 = 'root -q -b \'../scripts/ytBackground_subtraction.C("single_lepton_trigger","baseline","range_baseline",true,60,120,80,100,2,2,4,5)\' > baseline_template_range_baseline_pT15_20_eta08137_single_lepton_trigger.txt'
    command_b07 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger","baseline","range_baseline",true,60,120,80,100,2,2,4,5)\' > baseline_template_range_baseline_pT15_20_eta08137_dilepton_trigger.txt'
    command_b08 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger_tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,2,2,4,5)\' > baseline_template_range_baseline_pT15_20_eta08137_dilepton_trigger_tag_trigger_matched.txt'
    # 15 < pT < 20: 1.52 < |eta| < 2.01
    command_b09 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,2,2,7,8)\' > baseline_template_range_baseline_pT15_20_eta152201_tag_trigger_matched.txt'
    command_b10 = 'root -q -b \'../scripts/ytBackground_subtraction.C("single_lepton_trigger","baseline","range_baseline",true,60,120,80,100,2,2,7,8)\' > baseline_template_range_baseline_pT15_20_eta152201_single_lepton_trigger.txt'
    command_b11 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger","baseline","range_baseline",true,60,120,80,100,2,2,7,8)\' > baseline_template_range_baseline_pT15_20_eta152201_dilepton_trigger.txt'
    command_b12 = 'root -q -b \'../scripts/ytBackground_subtraction.C("dilepton_trigger_tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,2,2,7,8)\' > baseline_template_range_baseline_pT15_20_eta152201_dilepton_trigger_tag_trigger_matched.txt'

    commands.append(command_a01)
    commands.append(command_a02)
    commands.append(command_a03)
    commands.append(command_a04)
    commands.append(command_a05)
    commands.append(command_a06)
    commands.append(command_a07)
    commands.append(command_a08)
    commands.append(command_a09)
    commands.append(command_a10)
    commands.append(command_a11)
    commands.append(command_a12)

    commands.append(command_b01)
    commands.append(command_b02)
    commands.append(command_b03)
    commands.append(command_b04)
    commands.append(command_b05)
    commands.append(command_b06)
    commands.append(command_b07)
    commands.append(command_b08)
    commands.append(command_b09)
    commands.append(command_b10)
    commands.append(command_b11)
    commands.append(command_b12)

    for command in commands:
        os.system(command)

def get_results_after_background_subtraction():
    print "10 < pT < 15, 0 < |eta| < 0.8"
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta008_tag_trigger_matched.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta008_single_lepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta008_dilepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta008_dilepton_trigger_tag_trigger_matched.txt')
    print "10 < pT < 15, 0.8 < |eta| < 1.37"
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta08137_tag_trigger_matched.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta08137_single_lepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta08137_dilepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta08137_dilepton_trigger_tag_trigger_matched.txt')
    print "10 < pT < 15, 1.52 < |eta| < 2.01"
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta152201_tag_trigger_matched.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta152201_single_lepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta152201_dilepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT10_15_eta152201_dilepton_trigger_tag_trigger_matched.txt')
 
    print "15 < pT < 20, 0 < |eta| < 0.8"
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta008_tag_trigger_matched.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta008_single_lepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta008_dilepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta008_dilepton_trigger_tag_trigger_matched.txt')
    print "15 < pT < 20, 0.8 < |eta| < 1.37"
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta08137_tag_trigger_matched.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta08137_single_lepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta08137_dilepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta08137_dilepton_trigger_tag_trigger_matched.txt')
    print "15 < pT < 20, 1.37 < |eta| < 2.01"
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta152201_tag_trigger_matched.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta152201_single_lepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta152201_dilepton_trigger.txt')
    print extract_background_subtraction_results('baseline_template_range_baseline_pT15_20_eta152201_dilepton_trigger_tag_trigger_matched.txt')


def extract_background_subtraction_results(input_file):
    RLE_before_bkg = 0
    RLE_after_bkg = 0

    with open(input_file) as fopen:
        for line in fopen.readlines():
            if "real lepton efficiency (before background subtraction)" in line:
                RLE_before_bkg = line[line.find("=")+2:line.find("+-")-1]
            elif "real lepton efficiency (after background subtraction)" in line: 
                RLE_after_bkg = line[line.find("=")+2:line.find("+-")-1]
    fopen.close()

    return RLE_before_bkg, RLE_after_bkg

if __name__ == "__main__":
	main()
