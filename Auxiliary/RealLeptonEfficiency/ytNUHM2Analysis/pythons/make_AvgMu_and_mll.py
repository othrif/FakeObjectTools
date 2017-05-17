#!/usr/bin/python
import os, ROOT, getopt
os.sys.argv.append( '-b' ) #Batch mode, same as python <somescript>.py -b

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "b", ["skimmed", "RLE", "RLE-trigger="])
    except getopt.GetoptError:
        print 'make_AvgMu_and_mll.py -b --skimmed/--RLE' # -b: batch mode
        os.sys.exit(2)

    # Empty histograms
    h_data_AvgMu_OSee = ""
    h_data_AvgMu_OSmumu = ""
    h_Zee_AvgMu_OSee = ""
    h_Zmumu_AvgMu_OSmumu = ""

    h_data_baseline_OSee_mll = ""
    h_data_baseline_OSmumu_mll = ""
    h_data_signal_OSee_mll = ""
    h_data_signal_OSmumu_mll = ""

    h_Zee_baseline_OSee_mll_weighted = ""
    h_Zmumu_baseline_OSmumu_mll_weighted = ""
    h_Zee_signal_OSee_mll_weighted = ""
    h_Zmumu_signal_OSmumu_mll_weighted = ""

    prefix = ""
    suffix = ""

    for opt, arg in opts:
        if opt == "--skimmed":
            # For skimmed
            f_skimmed_data = ROOT.TFile("hist-merged-Data.root")
            f_skimmed_Zee = ROOT.TFile("skim_MC_Zee/hist-MC.root")
            f_skimmed_Zmumu = ROOT.TFile("skim_MC_Zmumu/hist-MC.root")

            h_skimmed_data_AvgMu_OSee = f_skimmed_data.Get("h_AvgMu_OSee")
            h_skimmed_data_AvgMu_OSmumu = f_skimmed_data.Get("h_AvgMu_OSmumu")
            h_skimmed_Zee_AvgMu_OSee = f_skimmed_Zee.Get("h_AvgMu_OSee")
            h_skimmed_Zmumu_AvgMu_OSmumu = f_skimmed_Zmumu.Get("h_AvgMu_OSmumu")

            h_skimmed_data_baseline_OSee_mll = f_skimmed_data.Get("h_baseline_OSee_mll")
            h_skimmed_data_baseline_OSmumu_mll = f_skimmed_data.Get("h_baseline_OSmumu_mll")
            h_skimmed_data_signal_OSee_mll = f_skimmed_data.Get("h_signal_OSee_mll")
            h_skimmed_data_signal_OSmumu_mll = f_skimmed_data.Get("h_signal_OSmumu_mll")

            h_skimmed_Zee_baseline_OSee_mll_weighted = f_skimmed_Zee.Get("h_baseline_OSee_mll")
            h_skimmed_Zmumu_baseline_OSmumu_mll_weighted = f_skimmed_Zmumu.Get("h_baseline_OSmumu_mll")
            h_skimmed_Zee_signal_OSee_mll_weighted = f_skimmed_Zee.Get("h_signal_OSee_mll")
            h_skimmed_Zmumu_signal_OSmumu_mll_weighted = f_skimmed_Zmumu.Get("h_signal_OSmumu_mll")

            h_data_AvgMu_OSee = h_skimmed_data_AvgMu_OSee
            h_data_AvgMu_OSmumu = h_skimmed_data_AvgMu_OSmumu
            h_Zee_AvgMu_OSee = h_skimmed_Zee_AvgMu_OSee
            h_Zmumu_AvgMu_OSmumu = h_skimmed_Zmumu_AvgMu_OSmumu

            h_data_baseline_OSee_mll = h_skimmed_data_baseline_OSee_mll
            h_data_baseline_OSmumu_mll = h_skimmed_data_baseline_OSmumu_mll
            h_data_signal_OSee_mll = h_skimmed_data_signal_OSee_mll
            h_data_signal_OSmumu_mll = h_skimmed_data_signal_OSmumu_mll

            h_Zee_baseline_OSee_mll_weighted = h_skimmed_Zee_baseline_OSee_mll_weighted
            h_Zmumu_baseline_OSmumu_mll_weighted = h_skimmed_Zmumu_baseline_OSmumu_mll_weighted
            h_Zee_signal_OSee_mll_weighted = h_skimmed_Zee_signal_OSee_mll_weighted
            h_Zmumu_signal_OSmumu_mll_weighted = h_skimmed_Zmumu_signal_OSmumu_mll_weighted

            prefix = "skimmed_"
        elif opt == "--RLE" or opt == "--RLE-trigger":
            # For RLE
            # Create NULL pointer
            f_RLE_data_el = None
            f_RLE_data_mu = None
            if opt == "--RLE":
                f_RLE_data_el = ROOT.TFile("hist-RLE-merged-data-elec.root")
                f_RLE_data_mu = ROOT.TFile("hist-RLE-merged-data-muon.root")
            elif opt == "--RLE-trigger":
                trigger = arg
                f_RLE_data_el = ROOT.TFile("hist-RLE-merged-data-elec_" + trigger + ".root")
                f_RLE_data_mu = ROOT.TFile("hist-RLE-merged-data-muon_" + trigger + ".root")
                suffix = "_" + trigger
            f_RLE_Zee = ROOT.TFile("RLE_MC_Zee/hist-20170222.root")
            f_RLE_Zmumu = ROOT.TFile("RLE_MC_Zmumu/hist-20170222.root")

            h_RLE_data_AvgMu_OSee = f_RLE_data_el.Get("h_AvgMu_OSee")
            h_RLE_data_AvgMu_OSmumu = f_RLE_data_mu.Get("h_AvgMu_OSmumu")
            h_RLE_Zee_AvgMu_OSee = f_RLE_Zee.Get("h_AvgMu_OSee")
            h_RLE_Zmumu_AvgMu_OSmumu = f_RLE_Zmumu.Get("h_AvgMu_OSmumu")

            h_RLE_data_baseline_OSee_mll = f_RLE_data_el.Get("h_baseline_OSee_mll")
            h_RLE_data_baseline_OSmumu_mll = f_RLE_data_mu.Get("h_baseline_OSmumu_mll")
            h_RLE_data_signal_OSee_mll = f_RLE_data_el.Get("h_signal_OSee_mll")
            h_RLE_data_signal_OSmumu_mll = f_RLE_data_mu.Get("h_signal_OSmumu_mll")

            h_RLE_Zee_baseline_OSee_mll_weighted = f_RLE_Zee.Get("h_baseline_OSee_mll")
            h_RLE_Zmumu_baseline_OSmumu_mll_weighted = f_RLE_Zmumu.Get("h_baseline_OSmumu_mll")
            h_RLE_Zee_signal_OSee_mll_weighted = f_RLE_Zee.Get("h_signal_OSee_mll")
            h_RLE_Zmumu_signal_OSmumu_mll_weighted = f_RLE_Zmumu.Get("h_signal_OSmumu_mll")

            h_data_AvgMu_OSee = h_RLE_data_AvgMu_OSee
            h_data_AvgMu_OSmumu = h_RLE_data_AvgMu_OSmumu
            h_Zee_AvgMu_OSee = h_RLE_Zee_AvgMu_OSee
            h_Zmumu_AvgMu_OSmumu = h_RLE_Zmumu_AvgMu_OSmumu

            h_data_baseline_OSee_mll = h_RLE_data_baseline_OSee_mll
            h_data_baseline_OSmumu_mll = h_RLE_data_baseline_OSmumu_mll
            h_data_signal_OSee_mll = h_RLE_data_signal_OSee_mll
            h_data_signal_OSmumu_mll = h_RLE_data_signal_OSmumu_mll

            h_Zee_baseline_OSee_mll_weighted = h_RLE_Zee_baseline_OSee_mll_weighted
            h_Zmumu_baseline_OSmumu_mll_weighted = h_RLE_Zmumu_baseline_OSmumu_mll_weighted
            h_Zee_signal_OSee_mll_weighted = h_RLE_Zee_signal_OSee_mll_weighted
            h_Zmumu_signal_OSmumu_mll_weighted = h_RLE_Zmumu_signal_OSmumu_mll_weighted

            prefix = "RLE_"

    c1 = ROOT.TCanvas("AvgMu_OSee")
    c1.SetLogy()
    h_data_AvgMu_OSee.SetLineColor(1) #kBlack
    h_data_AvgMu_OSee.SetMarkerColor(1)
    h_data_AvgMu_OSee.SetMarkerSize(1)
    h_data_AvgMu_OSee.SetMarkerStyle(20)
    h_data_AvgMu_OSee.Draw("E1")
    h_Zee_AvgMu_OSee.SetLineColor(2) #kRed
    h_Zee_AvgMu_OSee.Draw("hist,same")
    c1.SaveAs(prefix + "AvgMu_OSee" + suffix + ".pdf")

    c2 = ROOT.TCanvas("AvgMu_OSmumu")
    c2.SetLogy()
    h_data_AvgMu_OSmumu.SetLineColor(1)
    h_data_AvgMu_OSmumu.SetMarkerColor(1)
    h_data_AvgMu_OSmumu.SetMarkerSize(1)
    h_data_AvgMu_OSmumu.SetMarkerStyle(20)
    h_data_AvgMu_OSmumu.Draw("E1")
    h_Zmumu_AvgMu_OSmumu.SetLineColor(4) #kBlue
    h_Zmumu_AvgMu_OSmumu.Draw("hist,same")
    c2.SaveAs(prefix + "AvgMu_OSmumu" + suffix + ".pdf")

    c3 = ROOT.TCanvas("baseline_Mee")
    h_data_baseline_OSee_mll.SetLineColor(1)
    h_data_baseline_OSee_mll.SetMarkerColor(1)
    h_data_baseline_OSee_mll.SetMarkerSize(1)
    h_data_baseline_OSee_mll.SetMarkerStyle(20)
    h_data_baseline_OSee_mll.Draw("E1")
    h_Zee_baseline_OSee_mll_weighted.SetLineColor(2) #kRed
    h_Zee_baseline_OSee_mll_weighted.Draw("hist,same")
    c3.SaveAs(prefix + "baseline_Mee" + suffix + ".pdf")

    c4 = ROOT.TCanvas("baseline_Mmumu")
    h_data_baseline_OSmumu_mll.SetLineColor(1)
    h_data_baseline_OSmumu_mll.SetMarkerColor(1)
    h_data_baseline_OSmumu_mll.SetMarkerSize(1)
    h_data_baseline_OSmumu_mll.SetMarkerStyle(20)
    h_data_baseline_OSmumu_mll.Draw("E1")
    h_Zmumu_baseline_OSmumu_mll_weighted.SetLineColor(4) #kBlue
    h_Zmumu_baseline_OSmumu_mll_weighted.Draw("hist,same")
    c4.SaveAs(prefix + "baseline_Mmumu" + suffix + ".pdf")

    c5 = ROOT.TCanvas("signal_Mee")
    h_data_signal_OSee_mll.SetLineColor(1)
    h_data_signal_OSee_mll.SetMarkerColor(1)
    h_data_signal_OSee_mll.SetMarkerSize(1)
    h_data_signal_OSee_mll.SetMarkerStyle(20)
    h_data_signal_OSee_mll.Draw("E1")
    h_Zee_signal_OSee_mll_weighted.SetLineColor(2) #kRed
    h_Zee_signal_OSee_mll_weighted.Draw("hist,same")
    c5.SaveAs(prefix + "signal_Mee" + suffix + ".pdf")

    c6 = ROOT.TCanvas("signal_Mmumu")
    h_data_signal_OSmumu_mll.SetLineColor(1)
    h_data_signal_OSmumu_mll.SetMarkerColor(1)
    h_data_signal_OSmumu_mll.SetMarkerSize(1)
    h_data_signal_OSmumu_mll.SetMarkerStyle(20)
    h_data_signal_OSmumu_mll.Draw("E1")
    h_Zmumu_signal_OSmumu_mll_weighted.SetLineColor(4) #kBlue
    h_Zmumu_signal_OSmumu_mll_weighted.Draw("hist,same")
    c6.SaveAs(prefix + "signal_Mmumu" + suffix + ".pdf")

if __name__ == "__main__":
    main(os.sys.argv[1:])
