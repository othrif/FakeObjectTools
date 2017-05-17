#!/usr/bin/python
import os, ROOT, math

path = "~/Desktop/20170310/"

electron = "hist-RLE-merged-data-elec.root";
muon = "hist-RLE-merged-data-muon.root";

electron_tag_trigger_matched = "hist-RLE-merged-data-elec_tag_trigger_matched.root";
muon_tag_trigger_matched = "hist-RLE-merged-data-muon_tag_trigger_matched.root";

electron_dilepton_trigger = "hist-RLE-merged-data-elec_dilepton_trigger.root";
muon_dilepton_trigger = "hist-RLE-merged-data-muon_dilepton_trigger.root";

electron_dilepton_trigger_tag_trigger_matched = "hist-RLE-merged-data-elec_dilepton_trigger_tag_trigger_matched.root";
muon_dilepton_trigger_tag_trigger_matched = "hist-RLE-merged-data-muon_dilepton_trigger_tag_trigger_matched.root";

pt_ranges = [10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 150, 200]

eta_ranges = {
    "elec":[0, 0.8, 1.37, 1.52, 2.01],
    "muon":[0., 0.6, 1.2, 1.8, 2.5]
}

def main():
    trigger_systematic_uncertainty("electron")
    trigger_systematic_uncertainty("muon")

def trigger_systematic_uncertainty(lepton):
    # Naming:
    # sin: single_lepton_trigger
    # sin_tag: single_lepton_trigger_tag_trigger_matched
    # di: dilepton_trigge
    # di_tag: dilepton_trigger_tag_trigger_matched

    file_sin = None
    file_sin_tag = None
    file_di = None
    file_di_tag = None

    eta_range = None

    fo_name = ""

    if lepton == "electron":
        file_sin = ROOT.TFile(path + electron)
        file_sin_tag = ROOT.TFile(path + electron_tag_trigger_matched)
        file_di = ROOT.TFile(path + electron_dilepton_trigger)
        file_di_tag = ROOT.TFile(path + electron_dilepton_trigger_tag_trigger_matched)

        h_sin_baseline_pt_eta_dRjet = file_sin.Get("h_baseline_pt_eta_dRjet")
        h_sin_signal_pt_eta_dRjet = file_sin.Get("h_signal_pt_eta_dRjet")

        h_sin_tag_baseline_pt_eta_dRjet = file_sin_tag.Get("h_baseline_pt_eta_dRjet")
        h_sin_tag_signal_pt_eta_dRjet = file_sin_tag.Get("h_signal_pt_eta_dRjet")

        h_di_baseline_pt_eta_dRjet = file_di.Get("h_baseline_pt_eta_dRjet")
        h_di_signal_pt_eta_dRjet = file_di.Get("h_signal_pt_eta_dRjet")

        h_di_tag_baseline_pt_eta_dRjet = file_di_tag.Get("h_baseline_pt_eta_dRjet")
        h_di_tag_signal_pt_eta_dRjet = file_di_tag.Get("h_signal_pt_eta_dRjet")

        eta_range = eta_ranges["elec"]

        fo_name = "el_trigger_systematic.txt"
    elif lepton == "muon":
        file_sin = ROOT.TFile(path + muon)
        file_sin_tag = ROOT.TFile(path + muon_tag_trigger_matched)
        file_di = ROOT.TFile(path + muon_dilepton_trigger)
        file_di_tag = ROOT.TFile(path + muon_dilepton_trigger_tag_trigger_matched)

        h_sin_baseline_pt_eta_dRjet = file_sin.Get("h_baseline_pt_eta_dRjet")
        h_sin_signal_pt_eta_dRjet = file_sin.Get("h_signal_pt_eta_dRjet")

        h_sin_tag_baseline_pt_eta_dRjet = file_sin_tag.Get("h_baseline_pt_eta_dRjet")
        h_sin_tag_signal_pt_eta_dRjet = file_sin_tag.Get("h_signal_pt_eta_dRjet")

        h_di_baseline_pt_eta_dRjet = file_di.Get("h_baseline_pt_eta_dRjet")
        h_di_signal_pt_eta_dRjet = file_di.Get("h_signal_pt_eta_dRjet")

        h_di_tag_baseline_pt_eta_dRjet = file_di_tag.Get("h_baseline_pt_eta_dRjet")
        h_di_tag_signal_pt_eta_dRjet = file_di_tag.Get("h_signal_pt_eta_dRjet")    

        eta_range = eta_ranges["muon"]
    
        fo_name = "mu_trigger_systematic.txt"

    fo = open(fo_name, "w")
    fo.write("Variate triggers:\n")

    for i in range(1, h_sin_signal_pt_eta_dRjet.GetXaxis().GetNbins() + 1):
        for eta in eta_range:
            if lepton == "electron" and eta == 1.37:
                continue
            eta_index = eta_range.index(eta)
            if eta_index < len(eta_range) - 1:
                eta_next = eta_range[eta_index + 1]

                eta_bin_low = h_sin_signal_pt_eta_dRjet.GetYaxis().FindBin(eta + 0.01)
                eta_bin_up = h_sin_signal_pt_eta_dRjet.GetYaxis().FindBin(eta_next - 0.01)

                h_sin_baseline = h_sin_baseline_pt_eta_dRjet.ProjectionX("h_sin_baseline", eta_bin_low, eta_bin_up)
                h_sin_signal = h_sin_signal_pt_eta_dRjet.ProjectionX("h_sin_signal", eta_bin_low, eta_bin_up)

                h_sin_tag_baseline = h_sin_tag_baseline_pt_eta_dRjet.ProjectionX("h_sin_tag_baseline", eta_bin_low, eta_bin_up)
                h_sin_tag_signal = h_sin_tag_signal_pt_eta_dRjet.ProjectionX("h_sin_tag_signal", eta_bin_low, eta_bin_up)

                h_di_baseline = h_di_baseline_pt_eta_dRjet.ProjectionX("h_di_baseline", eta_bin_low, eta_bin_up)
                h_di_signal = h_di_signal_pt_eta_dRjet.ProjectionX("h_di_signal", eta_bin_low, eta_bin_up)

                h_di_tag_baseline = h_di_tag_baseline_pt_eta_dRjet.ProjectionX("h_di_tag_baseline", eta_bin_low, eta_bin_up)
                h_di_tag_signal = h_di_tag_signal_pt_eta_dRjet.ProjectionX("h_di_tag_signal", eta_bin_low, eta_bin_up)

                h_sin_eff = h_sin_signal.Clone()
                h_sin_tag_eff = h_sin_tag_signal.Clone()
                h_di_eff = h_di_signal.Clone()
                h_di_tag_eff = h_di_tag_signal.Clone()

                h_sin_eff.Sumw2()
                h_sin_tag_eff.Sumw2()
                h_di_eff.Sumw2()
                h_di_tag_eff.Sumw2()

                h_sin_eff.SetName("h_sin_eff")
                h_sin_tag_eff.SetName("h_sin_tag_eff")
                h_di_eff.SetName("h_di_eff")
                h_di_tag_eff.SetName("h_di_tag_eff")

                h_sin_eff.Divide(h_sin_baseline)
                h_sin_tag_eff.Divide(h_sin_tag_baseline)
                h_di_eff.Divide(h_di_baseline)
                h_di_tag_eff.Divide(h_di_tag_baseline)

            #for i in range(1, h_sin_eff.GetXaxis().GetNbins() + 1):
                # central_value = h_sin_eff.GetBinContent(i)
                # var1 = h_sin_tag_eff.GetBinContent(i) - central_value
                central_value = h_sin_tag_eff.GetBinContent(i)
                var1 = h_sin_eff.GetBinContent(i) - central_value
                var2 = h_di_eff.GetBinContent(i) - central_value
                var3 = h_di_tag_eff.GetBinContent(i) - central_value
                uncertainty = math.sqrt(var1**2 + var2**2 + var3**2)
                uncertainty_in_percentage = uncertainty / central_value * 100
                pt_low = int(h_sin_eff.GetXaxis().GetBinLowEdge(i))
                pt_up = int(h_sin_eff.GetXaxis().GetBinUpEdge(i))

                write_to_file = str(pt_low) + " GeV < pT < " + str(pt_up) + " GeV, " + str(eta) + " < |eta| < " + str(eta_next) + \
                                " trigger systematic uncertainty = " + str(uncertainty) + " (" + str(uncertainty_in_percentage) + "%)\n"
                fo.write(write_to_file)

    fo.close()

if __name__ == "__main__":
	main()
