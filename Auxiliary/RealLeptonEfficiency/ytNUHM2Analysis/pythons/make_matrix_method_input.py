#!/usr/bin/python
import os, math

def main():
    # remember to modify the input file path
    get_electron_real_efficiency()
    get_muon_real_efficiency()

    read_files = ["realElecEfficiencyTables.h",
                  "realMuEfficiencyTables.h"]

    # merge twe files
    with open("realEfficiencyTables.h", "wb") as outfile:
        for fname in read_files:
            with open(fname, "rb") as infile:
                outfile.write(infile.read())



def get_electron_real_efficiency():
    variation_files = ["/Users/ytshen/Desktop/20170310/logs/el_var_templates.txt",
                       "/Users/ytshen/Desktop/20170310/logs/el_var_fitting_ranges.txt",
                       "/Users/ytshen/Desktop/20170310/logs/el_var_mll_windows.txt"]

    pt_ranges = []

    central_value_eta_bin1 = []
    stat_uncertainty_eta_bin1 = []
    syst_uncertainty_eta_bin1 = []
    syst_uncertainty_eta_bin1_var_templates = []
    syst_uncertainty_eta_bin1_var_fitting_ranges = []
    syst_uncertainty_eta_bin1_var_mll_windows = []

    central_value_eta_bin2 = []
    stat_uncertainty_eta_bin2 = []
    syst_uncertainty_eta_bin2 = []
    syst_uncertainty_eta_bin2_var_templates = []
    syst_uncertainty_eta_bin2_var_fitting_ranges = []
    syst_uncertainty_eta_bin2_var_mll_windows = []

    central_value_eta_bin3 = []
    stat_uncertainty_eta_bin3 = []
    syst_uncertainty_eta_bin3 = []
    syst_uncertainty_eta_bin3_var_templates = []
    syst_uncertainty_eta_bin3_var_fitting_ranges = []
    syst_uncertainty_eta_bin3_var_mll_windows = []

    # Read the central value, statistic uncertainty, and systematic uncertainties caussing by variating templates
    with open(variation_files[0]) as fopen:
        lines = fopen.readlines()
        for line in lines:
            if "0.0 < eta < 0.8" in line:
                pt_ranges.append(str(line[0: line.find("GeV,")+3]))
                central_value_eta_bin1.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                stat_uncertainty_eta_bin1.append(float(line[line.find("$pm")+5 : line.find("(stat)")-1]))
                syst_uncertainty_eta_bin1_var_templates.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "0.8 < eta < 1.37" in line:
                central_value_eta_bin2.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                stat_uncertainty_eta_bin2.append(float(line[line.find("$pm")+5 : line.find("(stat)")-1]))
                syst_uncertainty_eta_bin2_var_templates.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "1.52 < eta < 2.01" in line:
                central_value_eta_bin3.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                stat_uncertainty_eta_bin3.append(float(line[line.find("$pm")+5 : line.find("(stat)")-1]))
                syst_uncertainty_eta_bin3_var_templates.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))

    # Now read the systematic uncertainties caussing by variating fitting ranges
    with open(variation_files[1]) as fopen:
        lines = fopen.readlines()
        for line in lines:
            if "0.0 < eta < 0.8" in line:
                syst_uncertainty_eta_bin1_var_fitting_ranges.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "0.8 < eta < 1.37" in line:
                syst_uncertainty_eta_bin2_var_fitting_ranges.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "1.52 < eta < 2.01" in line:
                syst_uncertainty_eta_bin3_var_fitting_ranges.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))

    # Now read the systematic uncertainties caussing by variating mll windows
    with open(variation_files[2]) as fopen:
        lines = fopen.readlines()
        for line in lines:
            if "0.0 < eta < 0.8" in line:
                syst_uncertainty_eta_bin1_var_mll_windows.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "0.8 < eta < 1.37" in line:
                syst_uncertainty_eta_bin2_var_mll_windows.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "1.52 < eta < 2.01" in line:
                syst_uncertainty_eta_bin3_var_mll_windows.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))

    # Calculate total systematic uncertainties using the quadratic sum
    syst_uncertainty_eta_bin1 = [math.sqrt(syst_uncertainty_eta_bin1_var_templates[i]**2 +
                                           syst_uncertainty_eta_bin1_var_fitting_ranges[i]**2 +
                                           syst_uncertainty_eta_bin1_var_mll_windows[i]**2) 
                                for i in range(0, len(syst_uncertainty_eta_bin1_var_templates))]

    syst_uncertainty_eta_bin2 = [math.sqrt(syst_uncertainty_eta_bin2_var_templates[i]**2 +
                                           syst_uncertainty_eta_bin2_var_fitting_ranges[i]**2 +
                                           syst_uncertainty_eta_bin2_var_mll_windows[i]**2) 
                                for i in range(0, len(syst_uncertainty_eta_bin2_var_templates))]

    syst_uncertainty_eta_bin3 = [math.sqrt(syst_uncertainty_eta_bin3_var_templates[i]**2 +
                                           syst_uncertainty_eta_bin3_var_fitting_ranges[i]**2 +
                                           syst_uncertainty_eta_bin3_var_mll_windows[i]**2) 
                                for i in range(0, len(syst_uncertainty_eta_bin1_var_templates))]

    # Calculate relative systemtic uncertainties
    relative_syst_uncert_bin1 = [syst_uncertainty_eta_bin1[i] / central_value_eta_bin1[i]
                                 for i in range(0, len(central_value_eta_bin1))]
    relative_syst_uncert_bin2 = [syst_uncertainty_eta_bin2[i] / central_value_eta_bin2[i]
                                 for i in range(0, len(central_value_eta_bin2))]
    relative_syst_uncert_bin3 = [syst_uncertainty_eta_bin3[i] / central_value_eta_bin3[i]
                                 for i in range(0, len(central_value_eta_bin3))]

    fo = open('realElecEfficiencyTables.h','w')
    fo.write("//\n// Electron real efficiency\n//\n")

    fo.write("const double Params_RealEffEl_Eff[] = {\n")
    fo.write("//eta=[0,0.8], [0.8,1.37], [1.52,2.01]\n")
    for i in range(0, len(central_value_eta_bin1)):
        if i != len(central_value_eta_bin1) - 1:
            fo.write(str(format(round(central_value_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin3[i], 4), '.4f')) + ", " +\
                     "//" + pt_ranges[i] + "\n")
        else:
            fo.write(str(format(round(central_value_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin3[i], 4), '.4f')) + "  "\
                     "//" + pt_ranges[i] + "\n")
    fo.write("};\n\n")

    # Use the absolute uncertainty for the statistical uncertainties
    # because matrix method need to use the absolute values
    fo.write("const double Params_RealEffEl_Stat[]={ // Use absolute uncertainties\n")
    fo.write("//eta=[0,0.8], [0.8,1.37], [1.52,2.01]\n")
    for i in range(0, len(stat_uncertainty_eta_bin1)):
        if i != len(stat_uncertainty_eta_bin1) - 1:
            fo.write(str(format(round(stat_uncertainty_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin3[i], 4), '.4f')) + ", " +\
                     "//" + pt_ranges[i] + "\n")
        else:
            fo.write(str(format(round(stat_uncertainty_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin3[i], 4), '.4f')) + "  "\
                     "//" + pt_ranges[i] + "\n")
    fo.write("};\n\n")

    # Use the relative uncertainties for the systematical uncertainties
    # because matrix method uses the relative one for systematics
    fo.write("const double Params_RealEffEl_SystCorr_relative[]={ // Use relative uncertainties\n")
    fo.write("//eta=[0,0.8], [0.8,1.37], [1.52,2.01]\n")
    for i in range(0, len(syst_uncertainty_eta_bin1)):
        if i != len(syst_uncertainty_eta_bin1) - 1:
            fo.write(str(format(round(relative_syst_uncert_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin3[i], 4), '.4f')) + ", " +\
                     "//" + pt_ranges[i] + "\n")
        else:
            fo.write(str(format(round(relative_syst_uncert_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin3[i], 4), '.4f')) + "  "\
                     "//" + pt_ranges[i] + "\n")

    fo.write("};\n\n")

    #
    # Read the uncertainty in busy environment (pT, dR)
    #
    # Use relative uncertainties for busy environment
    # Because I already calculate the relative uncertainties in ytRelative_difference_of_efficiency.C
    # , I just need to read the values from the input file and write them into the output file.
    busy_files = ["/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt1020_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt2030_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt3040_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt4050_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt5060_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt6080_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt80150_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt150200_eta020.txt"]

    busy_uncertainties = []

    for busy_file in busy_files:
        uncertainty = []
        with open(busy_file) as fopen:
            lines = fopen.readlines()
            read_flag = False
            for line in lines:
                if "Relative difference of efficiency GG_ttn1_elec:" in line:
                    read_flag = True
                if read_flag:
                    if "bin=1:" in line:
                        uncertainty.append(float(line[line.find("bin=1:")+7 : line.find("$pm$")-1]))
                    elif "bin=2:" in line:
                        uncertainty.append(float(line[line.find("bin=2:")+7 : line.find("$pm$")-1]))
                    elif "bin=3:" in line:
                        uncertainty.append(float(line[line.find("bin=3:")+7 : line.find("$pm$")-1]))
                    elif "bin=4:" in line:
                        uncertainty.append(float(line[line.find("bin=4:")+7 : line.find("$pm$")-1]))
                    elif "bin=5:" in line:
                        uncertainty.append(float(line[line.find("bin=5:")+7 : line.find("$pm$")-1]))
                    elif "bin=6:" in line:
                        uncertainty.append(float(line[line.find("bin=6:")+7 : line.find("$pm$")-1]))
                    elif "bin=7:" in line:
                        uncertainty.append(float(line[line.find("bin=7:")+7 : line.find("$pm$")-1]))
                    elif "bin=8:" in line:
                        uncertainty.append(float(line[line.find("bin=8:")+7 : line.find("$pm$")-1]))
                    elif "bin=9:" in line:
                        uncertainty.append(float(line[line.find("bin=9:")+7 : line.find("$pm$")-1]))
                    elif "bin=10:" in line:
                        uncertainty.append(float(line[line.find("bin=10:")+8 : line.find("$pm$")-1]))
                    elif "bin=11:" in line:
                        read_flag = False
        busy_uncertainties.append(uncertainty)
    
    fo.write("const double Params_RealEffEl_SystBusy_relative[] = { // Use relative uncertainties\n")
    fo.write("//dR=[0.,0.1], [0.1,0.15], [0.15,0.2], [0.2,0.3], [0.3,0.35], [0.35,0.4], [0.4,0.6], [0.6,1.], [1.,2.], [2.,4.]\n")
    for i in range(0, len(busy_uncertainties)):
        if i != len(busy_uncertainties) - 1:
            fo.write(str(format(round(busy_uncertainties[i][0], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][1], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][2], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][3], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][4], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][5], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][6], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][7], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][8], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][9], 4), '.4f')) + ", \n")
        else:
            fo.write(str(format(round(busy_uncertainties[i][0], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][1], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][2], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][3], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][4], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][5], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][6], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][7], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][8], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][9], 4), '.4f')) + "\n")
    fo.write("};\n\n")

    #fo.write("const float Params_RealEffEl_EtaBins[] = {0, 0.8, 1.37, 1.52, 2.01};\n" +\
    fo.write("const float Params_RealEffEl_EtaBins[] = {0, 0.8, 1.4, 2.01};\n" +\
             "const float Params_RealEffEl_PtBins[] = {10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 150, 200};\n" +\
             "const float Params_RealEffEl_dRBins[] = {0., 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.6, 1., 2., 4.};\n" +\
             "const unsigned int Params_RealEffEl_nEtaBins = sizeof(Params_RealEffEl_EtaBins)/sizeof(float)-1;\n" +\
             "const unsigned int Params_RealEffEl_nPtBins = sizeof(Params_RealEffEl_PtBins)/sizeof(float)-1;\n" +\
             "const unsigned int Params_RealEffEl_ndRBins = sizeof(Params_RealEffEl_dRBins)/sizeof(float)-1;\n")


def get_muon_real_efficiency():
    variation_file = "/Users/ytshen/Desktop/20170310/logs/mu_efficiency_systematics.txt"

    pt_ranges = []

    central_value_eta_bin1 = []
    total_uncertainty_eta_bin1 = []
    stat_uncertainty_eta_bin1 = []
    syst_uncertainty_eta_bin1 = []

    central_value_eta_bin2 = []
    total_uncertainty_eta_bin2 = []
    stat_uncertainty_eta_bin2 = []
    syst_uncertainty_eta_bin2 = []

    central_value_eta_bin3 = []
    total_uncertainty_eta_bin3 = []
    stat_uncertainty_eta_bin3 = []
    syst_uncertainty_eta_bin3 = []

    central_value_eta_bin4 = []
    total_uncertainty_eta_bin4 = []
    stat_uncertainty_eta_bin4 = []
    syst_uncertainty_eta_bin4 = []

    with open(variation_file) as fopen:
        lines = fopen.readlines()
        for line in lines:
            if "0 < eta < 0.6" in line:
                pt_ranges.append(str(line[0 : line.find("GeV,")+3]))
                central_value_eta_bin1.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin1.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
                stat_uncertainty_eta_bin1.append(float(line[line.find("(tot)")+8: line.find("(stat)")-1]))
                syst_uncertainty_eta_bin1.append(float(line[line.find("(stat)")+12: line.find("(syst)")-1]))
            elif "0.6 < eta < 1.2" in line:
                central_value_eta_bin2.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin2.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
                stat_uncertainty_eta_bin2.append(float(line[line.find("(tot)")+8: line.find("(stat)")-1]))
                syst_uncertainty_eta_bin2.append(float(line[line.find("(stat)")+12: line.find("(syst)")-1]))
            elif "1.2 < eta < 1.8" in line:
                central_value_eta_bin3.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin3.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
                stat_uncertainty_eta_bin3.append(float(line[line.find("(tot)")+8: line.find("(stat)")-1]))
                syst_uncertainty_eta_bin3.append(float(line[line.find("(stat)")+12: line.find("(syst)")-1]))
            elif "1.8 < eta < 2.5" in line:
                central_value_eta_bin4.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin4.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
                stat_uncertainty_eta_bin4.append(float(line[line.find("(tot)")+8: line.find("(stat)")-1]))
                syst_uncertainty_eta_bin4.append(float(line[line.find("(stat)")+12: line.find("(syst)")-1]))

    fo = open('realMuEfficiencyTables.h','w')
    fo.write("//\n// Mun real efficiency\n//\n")

    fo.write("const double Params_RealEffMu_Eff[] = {\n")
    fo.write("//eta=[0,0.6], [0.6,1.2], [1.2,1.8], [1.8,2.5]\n")
    for i in range(0, len(central_value_eta_bin1)):
        if i != len(central_value_eta_bin1) - 1:
            fo.write(str(format(round(central_value_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin3[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin4[i], 4), '.4f')) + ", " +\
                     "//" + pt_ranges[i] + "\n")
        else:
            fo.write(str(format(round(central_value_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin3[i], 4), '.4f')) + ", " +\
                     str(format(round(central_value_eta_bin4[i], 4), '.4f')) + "  "\
                     "//" + pt_ranges[i] + "\n")
    fo.write("};\n\n")

    # Use the absolute uncertainty for the statistical uncertainties
    # because matrix method need to use the absolute values
    fo.write("const double Params_RealEffMu_Stat[]={ // Use absolute uncertainties\n")
    fo.write("//eta=[0,0.6], [0.6,1.2], [1.2,1.8], [1.8,2.5]\n")
    for i in range(0, len(stat_uncertainty_eta_bin1)):
        if i != len(stat_uncertainty_eta_bin1) - 1:
            fo.write(str(format(round(stat_uncertainty_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin3[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin4[i], 4), '.4f')) + ", " +\
                     "//" + pt_ranges[i] + "\n")
        else:
            fo.write(str(format(round(stat_uncertainty_eta_bin1[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin2[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin3[i], 4), '.4f')) + ", " +\
                     str(format(round(stat_uncertainty_eta_bin4[i], 4), '.4f')) + "  "\
                     "//" + pt_ranges[i] + "\n")
    fo.write("};\n\n")

    # Use the relative uncertainties for the systematical uncertainties
    # because matrix method uses the relative one for systematics
    fo.write("const double Params_RealEffMu_SystCorr_relative[]={ // Use relative uncertainties\n")
    fo.write("//eta=[0,0.6], [0.6,1.2], [1.2,1.8], [1.8,2.5]\n")
    for i in range(0, len(syst_uncertainty_eta_bin1)):
        if i != len(syst_uncertainty_eta_bin1) - 1:
            relative_syst_uncert_bin1 = syst_uncertainty_eta_bin1[i] / central_value_eta_bin1[i]
            relative_syst_uncert_bin2 = syst_uncertainty_eta_bin2[i] / central_value_eta_bin2[i]
            relative_syst_uncert_bin3 = syst_uncertainty_eta_bin3[i] / central_value_eta_bin3[i]
            relative_syst_uncert_bin4 = syst_uncertainty_eta_bin4[i] / central_value_eta_bin4[i]
            fo.write(str(format(round(relative_syst_uncert_bin1, 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin2, 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin3, 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin4, 4), '.4f')) + ", " +\
                     "//" + pt_ranges[i] + "\n")
        else:
            fo.write(str(format(round(relative_syst_uncert_bin1, 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin2, 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin3, 4), '.4f')) + ", " +\
                     str(format(round(relative_syst_uncert_bin4, 4), '.4f')) + "  "\
                     "//" + pt_ranges[i] + "\n")

    fo.write("};\n\n")

    #
    # Read the uncertainty in busy environment (pT, dR)
    #
    # Use relative uncertainties for busy environment
    # Because I already calculate the relative uncertainties in ytRelative_difference_of_efficiency.C
    # , I just need to read the values from the input file and write them into the output file.
    busy_files = ["/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt1020_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt2030_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt3040_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt4050_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt5060_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt6080_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt80150_eta020.txt",
                  "/Users/ytshen/Desktop/20170310/logs/relative_difference_of_efficiency_pt150200_eta020.txt"]

    busy_uncertainties = []

    for busy_file in busy_files:
        uncertainty = []
        with open(busy_file) as fopen:
            lines = fopen.readlines()
            read_flag = False
            for line in lines:
                if "Relative difference of efficiency GG_ttn1_muon:" in line:
                    read_flag = True
                if read_flag:
                    if "bin=1:" in line:
                        uncertainty.append(float(line[line.find("bin=1:")+7 : line.find("$pm$")-1]))
                    elif "bin=2:" in line:
                        uncertainty.append(float(line[line.find("bin=2:")+7 : line.find("$pm$")-1]))
                    elif "bin=3:" in line:
                        uncertainty.append(float(line[line.find("bin=3:")+7 : line.find("$pm$")-1]))
                    elif "bin=4:" in line:
                        uncertainty.append(float(line[line.find("bin=4:")+7 : line.find("$pm$")-1]))
                    elif "bin=5:" in line:
                        uncertainty.append(float(line[line.find("bin=5:")+7 : line.find("$pm$")-1]))
                    elif "bin=6:" in line:
                        uncertainty.append(float(line[line.find("bin=6:")+7 : line.find("$pm$")-1]))
                    elif "bin=7:" in line:
                        uncertainty.append(float(line[line.find("bin=7:")+7 : line.find("$pm$")-1]))
                    elif "bin=8:" in line:
                        uncertainty.append(float(line[line.find("bin=8:")+7 : line.find("$pm$")-1]))
                    elif "bin=9:" in line:
                        uncertainty.append(float(line[line.find("bin=9:")+7 : line.find("$pm$")-1]))
                    elif "bin=10:" in line:
                        uncertainty.append(float(line[line.find("bin=10:")+8 : line.find("$pm$")-1]))
                    elif "bin=11:" in line:
                        read_flag = False
        busy_uncertainties.append(uncertainty)

    fo.write("const double Params_RealEffMu_SystBusy_relative[] = { // Use relative uncertainties\n")
    fo.write("//dR=[0.,0.1], [0.1,0.15], [0.15,0.2], [0.2,0.3], [0.3,0.35], [0.35,0.4], [0.4,0.6], [0.6,1.], [1.,2.], [2.,4.]\n")
    for i in range(0, len(busy_uncertainties)):
        if i != len(busy_uncertainties) - 1:
            fo.write(str(format(round(busy_uncertainties[i][0], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][1], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][2], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][3], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][4], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][5], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][6], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][7], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][8], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][9], 4), '.4f')) + ", \n")
        else:
            fo.write(str(format(round(busy_uncertainties[i][0], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][1], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][2], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][3], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][4], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][5], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][6], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][7], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][8], 4), '.4f')) + ", " +\
                     str(format(round(busy_uncertainties[i][9], 4), '.4f')) + "\n")
    fo.write("};\n\n")

    fo.write("const float Params_RealEffMu_EtaBins[] = {0, 0.6, 1.2, 1.8, 2.5};\n" +\
             "const float Params_RealEffMu_PtBins[] = {10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 150, 200};\n" +\
             "const float Params_RealEffMu_dRBins[] = {0., 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.6, 1., 2., 4.};\n" +\
             "const unsigned int Params_RealEffMu_nEtaBins = sizeof(Params_RealEffMu_EtaBins)/sizeof(float)-1;\n" +\
             "const unsigned int Params_RealEffMu_nPtBins = sizeof(Params_RealEffMu_PtBins)/sizeof(float)-1;\n" +\
             "const unsigned int Params_RealEffMu_ndRBins = sizeof(Params_RealEffMu_dRBins)/sizeof(float)-1;\n")

if __name__ == "__main__":
    main()
