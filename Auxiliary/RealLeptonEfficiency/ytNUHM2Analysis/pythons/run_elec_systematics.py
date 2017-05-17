#!/usr/bin/python
import os

def main():
    variate_template()
    variate_fitting_range()
    variate_mll_window()

templates = ['baseline_template', 'template1', 'template2']
fitting_ranges = ['range_baseline', 'range1', 'range2']
mll_windows = [[80, 100], [75, 105], [85, 95]]
#pt_ranges = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200]
pt_ranges = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 150, 200]
el_eta_ranges = [0., 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47]
#mu_eta_ranges = [0, 0.6, 1.2, 1.8, 2.5]
eta_ranges = [[1,3], [4,5], [7,8]]

# Get the output_filename
def mll_str(mll_window_low, mll_window_up):
    "Generate a string which contains the lower and upper bounds of mll."
    return "_mll" + str(mll_window_low) + str(mll_window_up)


def pt_str(pt_low, pt_up):
    "Generate a string which contains the lower and upper bounds of pt."
    return "_pT" + str(pt_low) + str(pt_up)


def eta_str(eta_low, eta_up):
    "Generate a string which contains the lower and upper bounds of eta."
    return "_eta" + str(int(eta_low * 100)) + str(int(eta_up * 100))


def calculate_real_electron_efficiencies_with_systematics(filename1, filename2, filename3):
    "Calculate the real lepton efficiencies from files and make a table for calculating systematics."
    central_values = []
    uncertainties = []
    with open(filename1) as fopen1:
        lines1 = fopen1.readlines()
        for line in lines1:
            if "pT1015" in filename1 or "pT1520" in filename1:
                if "real lepton efficiency (after background subtraction)" in line:
                    #print line
                    central_values.append(line[line.find("=")+2: line.find("+-")-1])
                    uncertainties.append(line[line.find("+-")+3: len(line)-1])
            else:
                if "real lepton efficiency (before background subtraction)" in line:
                    #print line
                    central_values.append(line[line.find("=")+1: line.find("+-")-1])
                    uncertainties.append(line[line.find("+-")+3: len(line)-1])

    with open(filename2) as fopen2:
        lines2 = fopen2.readlines()
        for line in lines2:
            if "pT1015" in filename2 or "pT1520" in filename2:
                if "real lepton efficiency (after background subtraction)" in line:
                    #print line
                    central_values.append(line[line.find("=")+2: line.find("+-")-1])
                    uncertainties.append(line[line.find("+-")+3: len(line)-1])
            else:
                if "real lepton efficiency (before background subtraction)" in line:
                    #print line
                    central_values.append(line[line.find("=")+2: line.find("+-")-1])
                    uncertainties.append(line[line.find("+-")+3: len(line)-1])

    with open(filename3) as fopen3:
        lines3 = fopen3.readlines()
        for line in lines3:
            if "pT1015" in filename3 or "pT1520" in filename3:
                if "real lepton efficiency (after background subtraction)" in line:
                    #print line
                    central_values.append(line[line.find("=")+2: line.find("+-")-1])
                    uncertainties.append(line[line.find("+-")+3: len(line)-1])
            else:
                if "real lepton efficiency (before background subtraction)" in line:
                    #print line
                    central_values.append(line[line.find("=")+2: line.find("+-")-1])
                    uncertainties.append(line[line.find("+-")+3: len(line)-1])

    #print str(central_values[0]) + " +- " + str(uncertainties[0]) + ", " + str(central_values[1]) + " +- " + str(uncertainties[1]) + ", " + str(central_values[2]) + " +- " + str(uncertainties[2])
    #systematics = (float(central_values[1]) - float(central_values[0]))**2 + (float(central_values[2]) - float(central_values[0]))**2
    systematics1 = float(central_values[1]) - float(central_values[0])
    systematics2 = float(central_values[2]) - float(central_values[0])
    systematics = max(abs(systematics1), abs(systematics2))
    #print "Real lepton efficiency = " + central_values[0] + " $pm$ " + uncertainties[0] + " (stat) $pm$ " + str(systematics) + " (syst)"
    return "Real lepton efficiency = " + central_values[0] + " $pm$ " + uncertainties[0] + " (stat) $pm$ " + str(systematics) + " (syst)"


def variate_template():
    "variate template and fix fitting range and mll."
    print "Now variate templates: baseline_template, template1, template2."
    for template in templates:
        for pt_range in pt_ranges:
            pt_range_index = pt_ranges.index(pt_range) + 1
            if pt_range_index < len(pt_ranges):
                pT = pt_str(pt_range, pt_ranges[pt_range_index])
                for eta_range in eta_ranges:
                    eta_range_low = eta_range[0]
                    eta_range_up = eta_range[1]
                    eta_low = el_eta_ranges[eta_range_low - 1]
                    eta_up = el_eta_ranges[eta_range_up]
                    eta = eta_str(eta_low, eta_up)
                    output_filename = template + "_range_baseline_mll80100" + pT + eta + ".txt"
                    if template is 'baseline_template':
                        if pt_range_index is 11:
                            command1 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,' + str(pt_range_index) + ',' + str(13) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                        else:
                            command1 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,80,100,' + str(pt_range_index) + ',' + str(pt_range_index) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                    else:
                        if pt_range_index is 11:
                            command1 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","' + template + '","range_baseline",true,60,120,80,100,' + str(pt_range_index) + ',' + str(11) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                        else:
                            command1 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","' + template + '","range_baseline",true,60,120,80,100,' + str(pt_range_index) + ',' + str(pt_range_index) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                    #print command1
                    os.system(command1)

    fo = open("el_var_templates.txt", "w")

    #print "Variate templates:"
    fo.write("Variate templates:\n")
    for pt_range in pt_ranges:
        pt_range_index = pt_ranges.index(pt_range) + 1
        if pt_range_index < len(pt_ranges):
            pT = pt_str(pt_range, pt_ranges[pt_range_index])
            for eta_range in eta_ranges:
                eta_range_low = eta_range[0]
                eta_range_up = eta_range[1]
                eta_low = el_eta_ranges[eta_range_low - 1]
                eta_up = el_eta_ranges[eta_range_up]
                eta = eta_str(eta_low, eta_up)
                file1 = "baseline_template_range_baseline_mll80100" + pT + eta + ".txt"
                file2 = "template1_range_baseline_mll80100" + pT + eta + ".txt"
                file3 = "template2_range_baseline_mll80100" + pT + eta + ".txt"
                #print str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " + str(eta_low) + " < eta < " + str(eta_up),
                #print calculate_real_electron_efficiencies_with_systematics(file1, file2, file3)
                write_to_file = str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " + str(eta_low) + " < eta < " + str(eta_up) + " " + \
                                calculate_real_electron_efficiencies_with_systematics(file1, file2, file3) + "\n"
                fo.write(write_to_file)

    fo.close()


def variate_fitting_range():
    "variate fitting range and fix template and mll."
    print "Now variate the fitting ranges: range_baseline, range1, and range2."
    for fitting_range in fitting_ranges:
        for pt_range in pt_ranges:
            pt_range_index = pt_ranges.index(pt_range) + 1
            if pt_range_index < len(pt_ranges):
                pT = pt_str(pt_range, pt_ranges[pt_range_index])
                for eta_range in eta_ranges:
                    eta_range_low = eta_range[0]
                    eta_range_up = eta_range[1]
                    eta_low = el_eta_ranges[eta_range_low - 1]
                    eta_up = el_eta_ranges[eta_range_up]
                    eta = eta_str(eta_low, eta_up)
                    output_filename = "baseline_template_" + fitting_range + "_mll80100" + pT + eta + ".txt"
                    if pt_range_index is 11:
                        command2 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","' + fitting_range + '",true,60,120,80,100,' + str(pt_range_index) + ',' + str(13) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                    else:
                        command2 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","' + fitting_range + '",true,60,120,80,100,' + str(pt_range_index) + ',' + str(pt_range_index) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                    #print command2
                    os.system(command2)

    fo = open("el_var_fitting_ranges.txt", "w")

    #print "Variate fitting ranges:"
    fo.write("Variate fitting ranges:\n")
    for pt_range in pt_ranges:
        pt_range_index = pt_ranges.index(pt_range) + 1
        if pt_range_index < len(pt_ranges):
            pT = pt_str(pt_range, pt_ranges[pt_range_index])
            for eta_range in eta_ranges:
                eta_range_low = eta_range[0]
                eta_range_up = eta_range[1]
                eta_low = el_eta_ranges[eta_range_low - 1]
                eta_up = el_eta_ranges[eta_range_up]
                eta = eta_str(eta_low, eta_up)
                file1 = "baseline_template_range_baseline_mll80100" + pT + eta + ".txt"
                file2 = "baseline_template_range1_mll80100" + pT + eta + ".txt"
                file3 = "baseline_template_range2_mll80100" + pT + eta + ".txt"
                #print str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " + str(eta_low) + " < eta < " + str(eta_up),
                #calculate_real_electron_efficiencies_with_systematics(file1, file2, file3)
                write_to_file = str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " + str(eta_low) + " < eta < " + str(eta_up) + " " + \
                                calculate_real_electron_efficiencies_with_systematics(file1, file2, file3) + "\n"
                fo.write(write_to_file)

    fo.close()


def variate_mll_window():
    "variate mll window and fix template and fitting range."
    print "Now variate the mll windows: [80, 100], [75, 105], and [85, 95]."
    for mll_window in mll_windows:
        mll = mll_str(mll_window[0], mll_window[1])
        for pt_range in pt_ranges:
            pt_range_index = pt_ranges.index(pt_range) + 1
            if pt_range_index < len(pt_ranges):
                pT = pt_str(pt_range, pt_ranges[pt_range_index])
                for eta_range in eta_ranges:
                    eta_range_low = eta_range[0]
                    eta_range_up = eta_range[1]
                    eta_low = el_eta_ranges[eta_range_low - 1]
                    eta_up = el_eta_ranges[eta_range_up]
                    eta = eta_str(eta_low, eta_up)
                    output_filename = "baseline_template_range_baseline" + mll + pT + eta + ".txt"
                    if pt_range_index is 11:
                        command3 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,' + str(mll_window[0]) + ',' + str(mll_window[1]) + ',' + str(pt_range_index) + ',' + str(13) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                    else:
                        command3 = 'root -q -b \'../scripts/ytBackground_subtraction.C("tag_trigger_matched","baseline","range_baseline",true,60,120,' + str(mll_window[0]) + ',' + str(mll_window[1]) + ',' + str(pt_range_index) + ',' + str(pt_range_index) + ',' + str(eta_range_low) + ',' + str(eta_range_up) + ')\' > ' + output_filename
                    #print command3
                    os.system(command3)

    fo = open("el_var_mll_windows.txt", "w")

    #print "Variate mll windows:"
    fo.write("Variate mll windows:\n")
    for pt_range in pt_ranges:
        pt_range_index = pt_ranges.index(pt_range) + 1
        if pt_range_index < len(pt_ranges):
            pT = pt_str(pt_range, pt_ranges[pt_range_index])
            for eta_range in eta_ranges:
                eta_range_low = eta_range[0]
                eta_range_up = eta_range[1]
                eta_low = el_eta_ranges[eta_range_low - 1]
                eta_up = el_eta_ranges[eta_range_up]
                eta = eta_str(eta_low, eta_up)
                file1 = "baseline_template_range_baseline_mll80100" + pT + eta + ".txt"
                file2 = "baseline_template_range_baseline_mll75105" + pT + eta + ".txt"
                file3 = "baseline_template_range_baseline_mll8595" + pT + eta + ".txt"
                #print str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " + str(eta_low) + " < eta < " + str(eta_up),
                #calculate_real_electron_efficiencies_with_systematics(file1, file2, file3)
                write_to_file = str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " + str(eta_low) + " < eta < " + str(eta_up) + " " + \
                                calculate_real_electron_efficiencies_with_systematics(file1, file2, file3) + "\n"
                fo.write(write_to_file)

    fo.close()


if __name__ == "__main__":
    main()
