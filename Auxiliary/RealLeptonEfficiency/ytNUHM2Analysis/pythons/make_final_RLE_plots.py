#!/usr/bin/python
import ROOT, array, math

def make_real_electron_efficiency_total_systematic_plot():
    central_value_eta_bin1 = []
    stat_err_eta_bin1 = []
    template_syst_err_eta_bin1 = []
    fitting_range_syst_err_eta_bin1 = []
    mll_syst_err_eta_bin1 = []

    central_value_eta_bin2 = []
    stat_err_eta_bin2 = []
    template_syst_err_eta_bin2 = []
    fitting_range_syst_err_eta_bin2 = []
    mll_syst_err_eta_bin2 = []

    central_value_eta_bin3 = []
    stat_err_eta_bin3 = []
    template_syst_err_eta_bin3 = []
    fitting_range_syst_err_eta_bin3 = []
    mll_syst_err_eta_bin3 = []

    #variation_files = ["../logs/elec_vary_templates.txt", "../logs/elec_vary_fitting_ranges.txt", "../logs/elec_vary_mll_windows.txt"]
    variation_files = ["/Users/ytshen/Desktop/20170310/logs/el_var_templates.txt",
                       "/Users/ytshen/Desktop/20170310/logs/el_var_fitting_ranges.txt",
                       "/Users/ytshen/Desktop/20170310/logs/el_var_mll_windows.txt"]

    with open(variation_files[0]) as fopen1:
        lines = fopen1.readlines()
        for line in lines:
            if "0.0 < eta < 0.8" in line:
                central_value_eta_bin1.append(float(line[line.find("=")+2 : line.find("$pm$")-1]))
                stat_err_eta_bin1.append(float(line[line.find("$pm$")+5 : line.find("(stat)")-1]))
                template_syst_err_eta_bin1.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "0.8 < eta < 1.37" in line:
                central_value_eta_bin2.append(float(line[line.find("=")+2 : line.find("$pm$")-1]))
                stat_err_eta_bin2.append(float(line[line.find("$pm$")+5 : line.find("(stat)")-1]))
                template_syst_err_eta_bin2.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "1.52 < eta < 2.01" in line:
                central_value_eta_bin3.append(float(line[line.find("=")+2 : line.find("$pm$")-1]))
                stat_err_eta_bin3.append(float(line[line.find("$pm$")+5 : line.find("(stat)")-1]))
                template_syst_err_eta_bin3.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))

    with open(variation_files[1]) as fopen2:
        lines = fopen2.readlines()
        for line in lines:
            if "0.0 < eta < 0.8" in line:
                fitting_range_syst_err_eta_bin1.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "0.8 < eta < 1.37" in line:
                fitting_range_syst_err_eta_bin2.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "1.52 < eta < 2.01" in line:
                fitting_range_syst_err_eta_bin3.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))

    with open(variation_files[2]) as fopen3:
        lines = fopen3.readlines()
        for line in lines:
            if "0.0 < eta < 0.8" in line:
                mll_syst_err_eta_bin1.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "0.8 < eta < 1.37" in line:
                mll_syst_err_eta_bin2.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))
            elif "1.52 < eta < 2.01" in line:
                mll_syst_err_eta_bin3.append(float(line[line.find("(stat)")+12 : line.find("(syst)")-1]))

    total_uncertainty_eta_bin1 = []
    total_uncertainty_eta_bin2 = []
    total_uncertainty_eta_bin3 = []

    for i in range(0, len(template_syst_err_eta_bin1)):
        total_uncertainty_eta_bin1.append(math.sqrt(template_syst_err_eta_bin1[i]**2 + fitting_range_syst_err_eta_bin1[i]**2 + mll_syst_err_eta_bin1[i]**2))
        total_uncertainty_eta_bin2.append(math.sqrt(template_syst_err_eta_bin2[i]**2 + fitting_range_syst_err_eta_bin2[i]**2 + mll_syst_err_eta_bin2[i]**2))
        total_uncertainty_eta_bin3.append(math.sqrt(template_syst_err_eta_bin3[i]**2 + fitting_range_syst_err_eta_bin3[i]**2 + mll_syst_err_eta_bin3[i]**2))

    pt_bin_center = array.array("d", [(10.+15.)/2., (15.+20.)/2., (20.+25.)/2, (25.+30.)/2, (30.+35.)/2, (35.+40.)/2, (40.+50.)/2, (50.+60.)/2, (60.+70.)/2, (70.+80.)/2, (80.+200.)/2])
    pt_bin_err = array.array("d", [(15.-10.)/2., (20.-15.)/2, (25.-20.)/2, (30.-25.)/2, (35.-30.)/2, (40.-35.)/2, (50.-40.)/2, (60.-50.)/2, (70.-60.)/2, (80.-70.)/2, (200.-80.)/2])

    eta_bin_1 = array.array("d", central_value_eta_bin1)
    eta_bin_1_err = array.array("d", total_uncertainty_eta_bin1)

    eta_bin_2 = array.array("d", central_value_eta_bin2)
    eta_bin_2_err = array.array("d", total_uncertainty_eta_bin2)

    eta_bin_3 = array.array("d", central_value_eta_bin3)
    eta_bin_3_err = array.array("d", total_uncertainty_eta_bin3)

    canvas = ROOT.TCanvas("c_ee", "Real electron efficiency", 600, 600)
    canvas.SetLeftMargin(0.12)
    #canvas.SetGrid()
    mg = ROOT.TMultiGraph()
    mg.SetTitle("")

    gr1 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_1, pt_bin_err, eta_bin_1_err)
    gr1.SetName("gr1")
    gr1.SetMarkerColor(1) #kBlack
    gr1.SetMarkerStyle(20) #kFullCircle
    gr1.SetLineColor(1)

    gr2 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_2, pt_bin_err, eta_bin_2_err)
    gr2.SetName("gr2")
    gr2.SetMarkerColor(2) #kRed
    gr2.SetMarkerStyle(21) #kFullSquare
    gr2.SetLineColor(2)

    gr3 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_3, pt_bin_err, eta_bin_3_err)
    gr3.SetName("gr3")
    gr3.SetMarkerColor(4) #kBlue
    gr3.SetMarkerStyle(22) #kFullTriangleUp
    gr3.SetLineColor(4)

    mg.Add(gr1)
    mg.Add(gr2)
    mg.Add(gr3)
    mg.Draw("AP")
    mg.GetXaxis().SetTitle("p_{T} [GeV]")
    mg.GetYaxis().SetTitle("Electron Real Efficiency")
    mg.GetYaxis().SetTitleOffset(1.5)
    mg.GetXaxis().SetLimits(10, 200);
    mg.SetMinimum(0.5)
    mg.SetMaximum(1.1)

    leg = ROOT.TLegend(0.5, 0.2, 0.9, 0.5)
    leg.SetHeader("Data 2015 + 2016: 36.5 fb^{-1}", "C")
    leg.AddEntry("gr1", "0. < |#eta| < 0.8", "LP")
    leg.AddEntry("gr2", "0.8 < |#eta| < 1.37", "LP")
    leg.AddEntry("gr3", "1.52 < |#eta| < 2.01", "LP")

    leg.SetBorderSize(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.03);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.Draw("same")

    line = ROOT.TLine(10., 1., 200., 1.)
    line.SetLineStyle(2)
    line.Draw()

    canvas.SaveAs("real_electron_efficiency_total_systematics.pdf", "pdf")



def make_real_muon_efficiency_total_systematic_plot():
    central_value_eta_bin1 = []
    total_uncertainty_eta_bin1 = []
    central_value_eta_bin2 = []
    total_uncertainty_eta_bin2 = []
    central_value_eta_bin3 = []
    total_uncertainty_eta_bin3 = []
    central_value_eta_bin4 = []
    total_uncertainty_eta_bin4 = []
    #with open("../logs/muon_variations.txt") as fopen:
    with open("/Users/ytshen/Desktop/20170310/logs/mu_efficiency_systematics.txt") as fopen:
        lines = fopen.readlines()
        for line in lines:
            if "0 < eta < 0.6" in line:
                central_value_eta_bin1.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin1.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
            elif "0.6 < eta < 1.2" in line:
                central_value_eta_bin2.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin2.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
            elif "1.2 < eta < 1.8" in line:
                central_value_eta_bin3.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin3.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))
            elif "1.8 < eta < 2.5" in line:
                central_value_eta_bin4.append(float(line[line.find("=")+2 : line.find("$pm")-1]))
                total_uncertainty_eta_bin4.append(float(line[line.find("$pm")+5 : line.find("(tot)")-1]))		

    pt_bin_center = array.array("d", [(10.+15.)/2., (15.+20.)/2., (20.+25.)/2, (25.+30.)/2, (30.+35.)/2, (35.+40.)/2, (40.+50.)/2, (50.+60.)/2, (60.+70.)/2, (70.+80.)/2, (80.+200.)/2])
    pt_bin_err = array.array("d", [(15.-10.)/2., (20.-15.)/2, (25.-20.)/2, (30.-25.)/2, (35.-30.)/2, (40.-35.)/2, (50.-40.)/2, (60.-50.)/2, (70.-60.)/2, (80.-70.)/2, (200.-80.)/2])

    eta_bin_1 = array.array("d", central_value_eta_bin1)
    eta_bin_1_err = array.array("d", total_uncertainty_eta_bin1)

    eta_bin_2 = array.array("d", central_value_eta_bin2)
    eta_bin_2_err = array.array("d", total_uncertainty_eta_bin2)

    eta_bin_3 = array.array("d", central_value_eta_bin3)
    eta_bin_3_err = array.array("d", total_uncertainty_eta_bin3)

    eta_bin_4 = array.array("d", central_value_eta_bin4)
    eta_bin_4_err = array.array("d", total_uncertainty_eta_bin4)

    canvas = ROOT.TCanvas("c_mumu", "Real muon efficiency", 600, 600)
    canvas.SetLeftMargin(0.12)
    #canvas.SetGrid()
    mg = ROOT.TMultiGraph()
    mg.SetTitle("")

    gr1 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_1, pt_bin_err, eta_bin_1_err)
    gr1.SetName("gr1")
    gr1.SetMarkerColor(1) #kBlack
    gr1.SetMarkerStyle(20) #kFullCircle
    gr1.SetLineColor(1)

    gr2 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_2, pt_bin_err, eta_bin_2_err)
    gr2.SetName("gr2")
    gr2.SetMarkerColor(2) #kRed
    gr2.SetMarkerStyle(21) #kFullSquare
    gr2.SetLineColor(2)

    gr3 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_3, pt_bin_err, eta_bin_3_err)
    gr3.SetName("gr3")
    gr3.SetMarkerColor(4) #kBlue
    gr3.SetMarkerStyle(22) #kFullTriangleUp
    gr3.SetLineColor(4)

    gr4 = ROOT.TGraphErrors(len(pt_bin_center), pt_bin_center, eta_bin_4, pt_bin_err, eta_bin_4_err)
    gr4.SetName("gr4")
    gr4.SetMarkerColor(3) #kGreen
    gr4.SetMarkerStyle(23) #kFullTriangleDown
    gr4.SetLineColor(3)

    mg.Add(gr1)
    mg.Add(gr2)
    mg.Add(gr3)
    mg.Add(gr4)
    mg.Draw("AP")
    mg.GetXaxis().SetTitle("p_{T} [GeV]")
    mg.GetYaxis().SetTitle("Muon Real Efficiency")
    mg.GetYaxis().SetTitleOffset(1.5)
    mg.GetXaxis().SetLimits(10, 200);
    mg.SetMinimum(0.5)
    mg.SetMaximum(1.1)

    leg = ROOT.TLegend(0.5, 0.2, 0.9, 0.5)
    leg.SetHeader("Data 2015 + 2016: 36.5 fb^{-1}", "C")
    leg.AddEntry("gr1", "0. < |#eta| < 0.6", "LP")
    leg.AddEntry("gr2", "0.6 < |#eta| < 1.2", "LP")
    leg.AddEntry("gr3", "1.2 < |#eta| < 1.8", "LP")
    leg.AddEntry("gr4", "1.8 < |#eta| < 2.5", "LP")

    leg.SetBorderSize(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.03);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.Draw("same")

    line = ROOT.TLine(10., 1., 200., 1.)
    line.SetLineStyle(2)
    line.Draw()

    canvas.SaveAs("real_muon_efficiency_total_systematics.pdf", "pdf")

make_real_electron_efficiency_total_systematic_plot()
make_real_muon_efficiency_total_systematic_plot()
