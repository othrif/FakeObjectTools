#!/usr/bin/python
import os, ROOT, math

def main():
    variate_mll_windows()

def calculate_real_muon_efficiency(mll_window_low = 80.,
                                   mll_window_up = 100.,
                                   pt_bin_low = 0,
                                   pt_bin_up = -1,
                                   eta_bin_low = 0,
                                   eta_bin_up = -1):

    "Calculate the real muon efficiency and systematics"

    debug = False
    verbose = False

    #print "Notice: Remember to modify the input file path!"
    path = "~/Desktop/20170310/"
    data_file = ROOT.TFile.Open(path + "hist-RLE-merged-data-muon_tag_trigger_matched.root")

    if debug is True:
        print hex(id(data_file))

    # Get 3-dim data histograms
    data_baseline_mll_3dim = data_file.Get("h_baseline_pt_eta_mll")
    data_signal_mll_3dim = data_file.Get("h_signal_pt_eta_mll")
    # project on Z axis
    data_baseline_mll = data_baseline_mll_3dim.ProjectionZ("data_baseline", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)
    data_signal_mll = data_signal_mll_3dim.ProjectionZ("data_signal", pt_bin_low, pt_bin_up, eta_bin_low, eta_bin_up)

    #print "pT range:" + str(int(data_signal_mll_3dim.GetXaxis().GetBinLowEdge(pt_bin_low))) + " to " \
    #                  + str(int(data_signal_mll_3dim.GetXaxis().GetBinUpEdge(pt_bin_up)))
    #print "eta range: " + str(data_signal_mll_3dim.GetYaxis().GetBinLowEdge(eta_bin_low)) + " to " \
    #                    + str(data_signal_mll_3dim.GetYaxis().GetBinUpEdge(eta_bin_up))
    
    if debug is True:
        print hex(id(data_baseline_mll)), hex(id(data_signal_mll))

    mll_window_low_bin = data_baseline_mll.FindBin(mll_window_low + 0.01);
    mll_window_up_bin = data_baseline_mll.FindBin(mll_window_up - 0.01);

    if debug is True:
        print mll_window_low_bin, mll_window_up_bin

    #
    # Number of events in the Zmass peak (80 < mll < 100)
    # Data
    data_baseline_mll_events_in_the_Z_peak = data_baseline_mll.Integral(mll_window_low_bin, mll_window_up_bin);
    data_signal_mll_events_in_the_Z_peak = data_signal_mll.Integral(mll_window_low_bin, mll_window_up_bin);
    data_signal_to_baseline_ratio_in_the_Z_peak = data_signal_mll_events_in_the_Z_peak / data_baseline_mll_events_in_the_Z_peak;
    data_non_signal_events_in_the_Z_peak = data_baseline_mll_events_in_the_Z_peak - data_signal_mll_events_in_the_Z_peak;

    if debug is True:
        print "data_baseline_mll_events_in_the_Z_peak=" + str(data_baseline_mll_events_in_the_Z_peak)
        print "data_signal_mll_events_in_the_Z_peak=" + str(data_signal_mll_events_in_the_Z_peak)
        print "data_signal_to_baseline_ratio_in_the_Z_peak=" + str(data_signal_to_baseline_ratio_in_the_Z_peak)
        print "data_non_signal_events_in_the_Z_peak=" + str(data_non_signal_events_in_the_Z_peak)

        #
        # Statistical information
        #
        print "********** Statistical information **********"
        print "Data in the Z peak:"
        print "Baseline in the Z peak = " + str(data_baseline_mll_events_in_the_Z_peak) + " +- " + str(math.sqrt(data_baseline_mll_events_in_the_Z_peak)),
        print "(" + str(100 * math.sqrt(data_baseline_mll_events_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak) + "%)"

        print "Signal in the Z peak = " + str(data_signal_mll_events_in_the_Z_peak) + " +- " + str(math.sqrt(data_signal_mll_events_in_the_Z_peak)),
        print "(" + str(100 * math.sqrt(data_signal_mll_events_in_the_Z_peak) / data_signal_mll_events_in_the_Z_peak) + "%)"

        print "Signal/Baseline ratio in the Z peak = " + str(data_signal_to_baseline_ratio_in_the_Z_peak) + " +- ",
        print str(math.sqrt(data_signal_to_baseline_ratio_in_the_Z_peak * (1 - data_signal_to_baseline_ratio_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak))
    
        print "(Baseline - Signal) in the Z peak = " + str(data_non_signal_events_in_the_Z_peak) + " +- " + str(math.sqrt(data_non_signal_events_in_the_Z_peak)),
        print "(" + str(100 * math.sqrt(data_non_signal_events_in_the_Z_peak) / data_non_signal_events_in_the_Z_peak) + "%)"
        print "********************"

    if verbose is True:
        print "Final results:"
        print "N(baseline) before background subtraction = " + str(data_baseline_mll_events_in_the_Z_peak) + " +- " + str(math.sqrt(data_baseline_mll_events_in_the_Z_peak)),
        print "(" + str(100 * math.sqrt(data_baseline_mll_events_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak) + "%)"
        print "N(signal) = " + str(data_signal_mll_events_in_the_Z_peak) + " +- " + str(math.sqrt(data_signal_mll_events_in_the_Z_peak)),
        print "(" + str(100 * math.sqrt(data_signal_mll_events_in_the_Z_peak) / data_signal_mll_events_in_the_Z_peak) + "%)"
        print "real lepton efficiency (before background subtraction) = " + str(data_signal_to_baseline_ratio_in_the_Z_peak) + " +- ",
        print math.sqrt(data_signal_to_baseline_ratio_in_the_Z_peak * (1 - data_signal_to_baseline_ratio_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak)
        print "********************"

    central_value = data_signal_to_baseline_ratio_in_the_Z_peak
    uncertainty = math.sqrt(data_signal_to_baseline_ratio_in_the_Z_peak * (1 - data_signal_to_baseline_ratio_in_the_Z_peak) / data_baseline_mll_events_in_the_Z_peak)
    return central_value, uncertainty

def mll_str(mll_window_low, mll_window_up):
    "Generate a string which contains the lower and upper bounds of mll."
    return "_mll" + str(mll_window_low) + str(mll_window_up)


def pt_str(pt_low, pt_up):
    "Generate a string which contains the lower and upper bounds of pt."
    return "_pT" + str(pt_low) + str(pt_up)


def eta_str(eta_low, eta_up):
    "Generate a string which contains the lower and upper bounds of eta."
    return "_eta" + str(int(eta_low * 100)) + str(int(eta_up * 100))


mll_windows = [[80, 100], [75, 105], [85, 95]]
pt_ranges = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 150, 200]
mu_eta_ranges = [0, 0.6, 1.2, 1.8, 2.5]


def variate_mll_windows():
    "variate mll window to calculate the systematics"

    fo = open("mu_efficiency_systematics.txt", "w")

    for pt_range in pt_ranges:
        pt_range_index = pt_ranges.index(pt_range) + 1
        if pt_range_index < len(pt_ranges):
            pT = pt_str(pt_range, pt_ranges[pt_range_index])
            for eta_range in mu_eta_ranges:
                eta_range_index = mu_eta_ranges.index(eta_range)
                if eta_range_index < len(mu_eta_ranges) - 1:
                    eta = eta_str(eta_range, mu_eta_ranges[mu_eta_ranges.index(eta_range) + 1])
                    #print "muon_variation" + mll_str(80, 100) + pT + eta
                    central_value_0, uncertainty_0 = calculate_real_muon_efficiency(80, 100, pt_range_index, pt_range_index, eta_range_index + 1, eta_range_index + 1)
                    #print central_value_0, uncertainty_0
                    
                    #print "muon_variation" + mll_str(75, 85) + pT + eta
                    central_value_1, uncertainty_1 = calculate_real_muon_efficiency(75, 105, pt_range_index, pt_range_index, eta_range_index + 1, eta_range_index + 1)
                    #print central_value_1, uncertainty_1
                    
                    #print "muon_variation" + mll_str(85, 95) + pT + eta
                    central_value_2, uncertainty_2 = calculate_real_muon_efficiency(85, 95, pt_range_index, pt_range_index, eta_range_index + 1, eta_range_index + 1)
                    #print central_value_2, uncertainty_2

                    systematics = math.sqrt((central_value_1 - central_value_0)**2 + (central_value_2 - central_value_0)**2)
                    total_uncertainty = math.sqrt(uncertainty_0**2 + systematics**2)
                    # print to screen
                    print str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, ",
                    print str(eta_range) + " < eta < " + str(mu_eta_ranges[eta_range_index + 1]) + " ", 
                    print "Real lepton efficiency = " + str(central_value_0) + " $pm$ " + str(total_uncertainty) + " (tot) : " + str(uncertainty_0) + " (stat) $pm$ " + str(systematics) + " (syst)"
                    # write to a file
                    fo.write(str(pt_range) + " GeV < pT < " + str(pt_ranges[pt_range_index]) + " GeV, " +\
                             str(eta_range) + " < eta < " + str(mu_eta_ranges[eta_range_index + 1]) + " " +\
                             "Real lepton efficiency = " + str(central_value_0) + " $pm$ " + str(total_uncertainty) + " (tot) : " + str(uncertainty_0) + " (stat) $pm$ " + str(systematics) + " (syst)\n")

    fo.close()

if __name__ == "__main__":
    main()
