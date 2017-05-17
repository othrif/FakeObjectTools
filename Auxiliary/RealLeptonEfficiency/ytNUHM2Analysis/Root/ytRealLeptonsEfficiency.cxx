#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ytNUHM2Analysis/ytRealLeptonsEfficiency.h>

// this is needed to distribute the algorithm to the workers
ClassImp(ytRealLeptonsEfficiency)

// Define the binnings
const float ytRealLeptonsEfficiency::m_default_mll_bins[91] = {
    60., 61., 62., 63., 64., 65., 66., 67., 68., 69.,
    70., 71., 72., 73., 74., 75., 76., 77., 78., 79.,
    80., 81., 82., 83., 84., 85., 86., 87., 88., 89.,
    90., 91., 92., 93., 94., 95., 96., 97., 98., 99.,
    100., 101., 102., 103., 104., 105., 106., 107., 108., 109.,
    110., 111., 112., 113., 114., 115., 116., 117., 118., 119.,
    120., 121., 122., 123., 124., 125., 126., 127., 128., 129.,
    130., 131., 132., 133., 134., 135., 136., 137., 138., 139.,
    140., 141., 142., 143., 144., 145., 146., 147., 148., 149.,
    150.
};

const float ytRealLeptonsEfficiency::m_default_pt_bins[14] = {
    10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 120, 150, 200
};

const float ytRealLeptonsEfficiency::m_default_eta_bins[21]  = {
    -2.47, -2.37, -2.01, -1.81, -1.52, -1.37, -1.15, -0.8, -0.6, -0.1, 0., 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47
};

const float ytRealLeptonsEfficiency::m_default_deltaR_bins[21] = {
    0, 0.2, 0.4, 0.6, 0.8,
    1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.2, 2.4, 2.6, 2.8,
    3.0, 3.2, 3.4, 3.6, 3.8,
    4.0
};

const float ytRealLeptonsEfficiency::m_default_NJet_bins[11] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
};
/*
const float ytRealLeptonsEfficiency::m_default_Etmiss_bins[9] = {
    0., 50. , 100., 150., 200. , 250., 300., 350., 400.
};
const float ytRealLeptonsEfficiency::m_default_meff_bins[7] = {
    0., 250., 500., 750., 1000., 1250., 1500.
};
*/
const float ytRealLeptonsEfficiency::m_default_AvgMu_bins[51] = {
    0., 1., 2., 3., 4., 5., 6., 7., 8., 9.,
    10., 11., 12., 13., 14., 15., 16., 17., 18., 19.,
    20., 21., 22., 23., 24., 25., 26., 27., 28., 29.,
    30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
    40., 41., 42., 43., 44., 45., 46., 47., 48., 49.,
    50.
};

const float ytRealLeptonsEfficiency::m_muon_eta_bins[21] = {
    -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5
};

const float ytRealLeptonsEfficiency::m_abs_eta_bins[11] = {
    0., 0.1, 0.6, 0.8, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47
};

const float ytRealLeptonsEfficiency::m_muon_abs_eta_bins[11] = {
    0, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5
};

const float ytRealLeptonsEfficiency::m_electron_coarse_pt_bins[12] = {
    10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 150, 200
};

const float ytRealLeptonsEfficiency::m_electron_coarse_eta_bins[6] = {
    0, 0.8, 1.37, 1.52, 2.01, 2.47
};

const float ytRealLeptonsEfficiency::m_muon_coarse_pt_bins[12] = {
    10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 150, 200
};

const float ytRealLeptonsEfficiency::m_muon_coarse_eta_bins[5] = {
    //0, 0.5, 1.25 , 1.75 , 2.5
    0, 0.6, 1.2, 1.8, 2.5
};
/*
const float ytRealLeptonsEfficiency::m_coarse_deltaR_bins[7] = {
    0, 0.2, 0.4, 0.6, 1.0, 2.0, 4.0
};
*/
const float ytRealLeptonsEfficiency::m_fine_deltaR_bins[36] = {
    0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.2, 2.4, 2.6, 2.8,
    3.0, 3.2, 3.4, 3.6, 3.8,
    4.0
};

ytRealLeptonsEfficiency::ytRealLeptonsEfficiency ()
{
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

    //const char * function_name = "ytRealLeptonsEfficiency()";
    //Info(function_name, "Function calls");

    this->set_binning_default();

    n_mll_bins = m_mll_bins.size() - 1;
    n_pt_bins = m_pt_bins.size() - 1;
    n_eta_bins = m_eta_bins.size() - 1;
    n_deltaR_bins = m_deltaR_bins.size() - 1;
    n_NJet_bins = m_NJet_bins.size() - 1;
    //n_Etmiss_bins = m_Etmiss_bins.size() - 1;
    //n_meff_bins = m_meff_bins.size() - 1;
    n_AvgMu_bins = m_AvgMu_bins.size() - 1;
/*
    // For debug
    cout << "n_mll_bins=" << n_mll_bins << endl;
    cout << "n_pt_bins=" << n_pt_bins << endl;
    cout << "n_eta_bins=" << n_eta_bins << endl;
    cout << "n_deltaR_bins=" << n_deltaR_bins << endl;
    cout << "n_NJet_bins=" << n_NJet_bins << endl;
    cout << "n_Etmiss_bins=" << n_Etmiss_bins << endl;
    cout << "n_meff_bins=" << n_meff_bins << endl;
*/
}



EL::StatusCode ytRealLeptonsEfficiency::setupJob (EL::Job& job)
{
    // Here you put code that sets up the job on the submission object
    // so that it is ready to work with your algorithm, e.g. you can
    // request the D3PDReader service or add output files.  Any code you
    // put here could instead also go into the submission script.  The
    // sole advantage of putting it here is that it gets automatically
    // activated/deactivated when you add/remove the algorithm from your
    // job, which may or may not be of value to you.

    //const char * function_name = "setupJob()";
    //Info(function_name, "Function calls");

    return EL::StatusCode::SUCCESS;
}




EL::StatusCode ytRealLeptonsEfficiency::histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

    //const char * function_name = "histInitialize()";
    //Info(function_name, "Function calls");

    // Use fine deltaR binning
    this->set_binning_deltaR(sizeof(m_fine_deltaR_bins) / sizeof(m_fine_deltaR_bins[0]), m_fine_deltaR_bins);
    n_deltaR_bins = m_deltaR_bins.size() - 1;

    if (lepton == "electron") {
        this->set_binning_eta(sizeof(m_abs_eta_bins) / sizeof(m_abs_eta_bins[0]), m_abs_eta_bins);
    }
    if (lepton == "muon") {
        this->set_binning_eta(sizeof(m_muon_coarse_eta_bins) / sizeof(m_muon_coarse_eta_bins[0]), m_muon_coarse_eta_bins);
    }
    n_eta_bins = m_eta_bins.size() - 1;
/*
    // For debug
    cout << "m_mll_bins are" << endl;
    for (unsigned int i = 0; i < m_mll_bins.size(); i++) {
        cout << m_mll_bins[i] << endl;
    }
    cout << "n_mll_bins=" << n_mll_bins << endl;
    cout << "m_pt_bins are" << endl;
    for (unsigned int i = 0; i < m_pt_bins.size(); i++) {
        cout << m_pt_bins[i] << endl;
    }
    cout << "n_pt_bins=" << n_pt_bins << endl;
    cout << "m_eta_bins are" << endl;
    for (unsigned int i = 0; i < m_eta_bins.size(); i++) {
        cout << m_eta_bins[i] << endl;
    }
    cout << "n_eta_bins=" << n_eta_bins << endl;
    cout << "m_deltaR_bins are" << endl;
    for (unsigned int i = 0; i < m_deltaR_bins.size(); i++) {
        cout << m_deltaR_bins[i] << endl;
    }
    cout << "n_deltaR_bins=" << n_deltaR_bins << endl;
*/
    h_Nvtx          = new TH1F("h_Nvtx", "Nvtx;Nvtx;Events", 50, 0 , 50);
    h_AvgMu         = new TH1F("h_AvgMu", "AvgMu;<#mu>;Events", 50, 0 , 50);
    h_AvgMu_OSee    = new TH1F("h_AvgMu_OSee", "AvgMu;<#mu>;Events", 50, 0 , 50);
    h_AvgMu_OSmumu  = new TH1F("h_AvgMu_OSmumu", "AvgMu;<#mu>;Events", 50, 0 , 50);
    h_zPV           = new TH1F("h_zPV", "zPV", 200, -200, 200);
    //h_NLepts    = new TH1F("h_NLepts", "Number of leptons;N_{leptons};Events",10, 0, 10);
    //h_NJets     = new TH1F("h_NJets", "Number of jets;N_{jets};Events", 40, 0, 40);

    h_Nvtx->Sumw2();
    h_AvgMu->Sumw2();
    h_AvgMu_OSee->Sumw2();
    h_AvgMu_OSmumu->Sumw2();
    h_zPV->Sumw2();
    //h_NLepts->Sumw2();
    //h_NJets->Sumw2();

    wk()->addOutput(h_Nvtx);
    wk()->addOutput(h_AvgMu);
    wk()->addOutput(h_AvgMu_OSee);
    wk()->addOutput(h_AvgMu_OSmumu);
    wk()->addOutput(h_zPV);
    //wk()->addOutput(h_NLepts);
    //wk()->addOutput(h_NJets);

    h_baseline_OSee_mll     = new TH1F("h_baseline_OSee_mll", "mll;M_{ll} [GeV];Events", 90, 60, 150);
    h_baseline_OSmumu_mll   = new TH1F("h_baseline_OSmumu_mll", "mll;M_{ll} [GeV];Events", 90, 60, 150);
    h_signal_OSee_mll       = new TH1F("h_signal_OSee_mll", "mll;M_{ll} [GeV];Events", 90, 60, 150);
    h_signal_OSmumu_mll     = new TH1F("h_signal_OSmumu_mll", "mll;M_{ll} [GeV];Events", 90, 60, 150);

    h_baseline_OSee_mll->Sumw2();
    h_baseline_OSmumu_mll->Sumw2();
    h_signal_OSee_mll->Sumw2();
    h_signal_OSmumu_mll->Sumw2();

    wk()->addOutput(h_baseline_OSee_mll);
    wk()->addOutput(h_baseline_OSmumu_mll);
    wk()->addOutput(h_signal_OSee_mll);
    wk()->addOutput(h_signal_OSmumu_mll);

    //h_Nvtx_weighted         = new TH1F("h_Nvtx_weighted", "Nvtx (weighted);Nvtx;Events", 50, 0 , 50);
    //h_AvgMu_weighted        = new TH1F("h_AvgMu_weighted", "AvgMu (weighted);<#mu>;Events", 50, 0 , 50);
    //h_AvgMu_OSee_weighted   = new TH1F("h_AvgMu_OSee_weighted", "AvgMu (weighted);<#mu>;Events", 50, 0 , 50);
    //h_AvgMu_OSmumu_weighted = new TH1F("h_AvgMu_OSmumu_weighted", "AvgMu (weighted);<#mu>;Events", 50, 0 , 50);
    //h_zPV_weighted          = new TH1F("h_zPV_weighted", "zPV (weighted)", 200, -200, 200);
    //h_NLepts_weighted    = new TH1F("h_NLepts_weighted", "Number of leptons (weighted);N_{leptons};Events",10, 0, 10);
    //h_NJets_weighted     = new TH1F("h_NJets_weighted", "Number of jets (weighted);N_{jets};Events", 40, 0, 40);

    //h_Nvtx_weighted->Sumw2();
    //h_AvgMu_weighted->Sumw2();
    //h_AvgMu_OSee_weighted->Sumw2();
    //h_AvgMu_OSmumu_weighted->Sumw2();
    //h_zPV_weighted->Sumw2();
    //h_NLepts_weighted->Sumw2();
    //h_NJets_weighted->Sumw2();

    //wk()->addOutput(h_Nvtx_weighted);
    //wk()->addOutput(h_AvgMu_weighted);
    //wk()->addOutput(h_AvgMu_OSee_weighted);
    //wk()->addOutput(h_AvgMu_OSmumu_weighted);
    //wk()->addOutput(h_zPV_weighted);
    //wk()->addOutput(h_NLepts_weighted);
    //wk()->addOutput(h_NJets_weighted);

    //h_baseline_OSee_mll_weighted    = new TH1F("h_baseline_OSee_mll_weighted", "mll;M_{ll} [GeV];Events", 90, 60, 150);
    //h_baseline_OSmumu_mll_weighted  = new TH1F("h_baseline_OSmumu_mll_weighted", "mll;M_{ll} [GeV];Events", 90, 60, 150);
    //h_signal_OSee_mll_weighted      = new TH1F("h_signal_OSee_mll_weighted", "mll;M_{ll} [GeV];Events", 90, 60, 150);
    //h_signal_OSmumu_mll_weighted    = new TH1F("h_signal_OSmumu_mll_weighted", "mll;M_{ll} [GeV];Events", 90, 60, 150);

    //h_baseline_OSee_mll_weighted->Sumw2();
    //h_baseline_OSmumu_mll_weighted->Sumw2();
    //h_signal_OSee_mll_weighted->Sumw2();
    //h_signal_OSmumu_mll_weighted->Sumw2();

    //wk()->addOutput(h_baseline_OSee_mll_weighted);
    //wk()->addOutput(h_baseline_OSmumu_mll_weighted);
    //wk()->addOutput(h_signal_OSee_mll_weighted);
    //wk()->addOutput(h_signal_OSmumu_mll_weighted);

    // Baseline histograms
    // 1-dim histograms
    h_baseline_mll          = new TH1F("h_baseline_mll", "mll;M_{ll} [GeV];Events", n_mll_bins, &m_mll_bins[0]);
    h_baseline_pt           = new TH1F("h_baseline_pt", "pt;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);
    h_baseline_eta          = new TH1F("h_baseline_eta", "eta;#eta;Events", n_eta_bins, &m_eta_bins[0]);
    //h_baseline_d0pvtx       = new TH1F("h_baseline_d0pvtx", "d0pvtx", 100, -0.15, 0.15);
    //h_baseline_sigd0        = new TH1F("h_baseline_sigd0", "sigd0", 100, 0, 10);
    //h_baseline_d0err        = new TH1F("h_baseline_d0err", "d0err", 100, -0.1, 0.1);
    //h_baseline_z0sinTheta   = new TH1F("h_baseline_z0SinTheta", "fabs(z0sinTheta)", 100, 0, 0.5);
    //h_baseline_ptvarcone20  = new TH1F("h_baseline_ptvarcone20", "ptvarcone20", 100, 0, 0.2);
    //h_baseline_ptvarcone30  = new TH1F("h_baseline_ptvarcone30", "ptvarcone30", 100, 0, 0.2);
    //h_baseline_topoetcone20 = new TH1F("h_baseline_topoetcone20", "topoetcone20", 100, 0, 0.2);
    h_baseline_nJets        = new TH1F("h_baseline_nJets", "NJet", n_NJet_bins, &m_NJet_bins[0]);
    h_baseline_dRjet        = new TH1F("h_baseline_dRjet", "dRjet", n_deltaR_bins, &m_deltaR_bins[0]);
    //h_baseline_Etmiss       = new TH1F("h_baseline_Etmiss", "Etmiss", n_Etmiss_bins, &m_Etmiss_bins[0]);
    //h_baseline_meff         = new TH1F("h_baseline_meff", "meff", n_meff_bins, &m_meff_bins[0]);

    h_baseline_mll->Sumw2();
    h_baseline_pt->Sumw2();
    h_baseline_eta->Sumw2();
    //h_baseline_d0pvtx->Sumw2();
    //h_baseline_sigd0->Sumw2();
    //h_baseline_d0err->Sumw2();
    //h_baseline_z0sinTheta->Sumw2();
    //h_baseline_ptvarcone20->Sumw2();
    //h_baseline_ptvarcone30->Sumw2();
    //h_baseline_topoetcone20->Sumw2();
    h_baseline_nJets->Sumw2();
    h_baseline_dRjet->Sumw2();
    //h_baseline_Etmiss->Sumw2();
    //h_baseline_meff->Sumw2();

    wk()->addOutput(h_baseline_mll);
    wk()->addOutput(h_baseline_pt);
    wk()->addOutput(h_baseline_eta);
    //wk()->addOutput(h_baseline_d0pvtx);
    //wk()->addOutput(h_baseline_sigd0);
    //wk()->addOutput(h_baseline_d0err);
    //wk()->addOutput(h_baseline_z0sinTheta);
    //wk()->addOutput(h_baseline_ptvarcone20);
    //wk()->addOutput(h_baseline_ptvarcone30);
    //wk()->addOutput(h_baseline_topoetcone20);
    wk()->addOutput(h_baseline_nJets);
    wk()->addOutput(h_baseline_dRjet);
    //wk()->addOutput(h_baseline_Etmiss);
    //wk()->addOutput(h_baseline_meff);

    h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt = new TH1F("h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt", "LooseAndBLayerLLH_to_MediumLLH_cut;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);
    h_baseline_CaloIso_cut_pt  = new TH1F("h_baseline_CaloIso_cut_pt", "CaloIso_cut;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);
    h_baseline_TrackIso_cut_pt = new TH1F("h_baseline_TrackIso_cut_pt", "TrackIso_cut;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);
    h_baseline_Iso_cut_pt      = new TH1F("h_baseline_Iso_cut_pt", "Iso_cut;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);
    h_baseline_z0_cut_pt       = new TH1F("h_baseline_z0_cut_pt", "z0_cut;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);
    h_baseline_sigd0_cut_pt    = new TH1F("h_baseline_sigd0_cut_pt", "sigd0_cut;p_{T} [GeV];Events", n_pt_bins, &m_pt_bins[0]);

    h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt->Sumw2();
    h_baseline_CaloIso_cut_pt->Sumw2();
    h_baseline_TrackIso_cut_pt->Sumw2();
    h_baseline_Iso_cut_pt->Sumw2();
    h_baseline_z0_cut_pt->Sumw2();
    h_baseline_sigd0_cut_pt->Sumw2();

    wk()->addOutput(h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt);
    wk()->addOutput(h_baseline_CaloIso_cut_pt);
    wk()->addOutput(h_baseline_TrackIso_cut_pt);
    wk()->addOutput(h_baseline_Iso_cut_pt);
    wk()->addOutput(h_baseline_z0_cut_pt);
    wk()->addOutput(h_baseline_sigd0_cut_pt);

    // 2-dim histograms
    h_baseline_pt_eta       = new TH2F("h_baseline_pt_eta", "p_{T} vs #eta;p_{T} [GeV];#eta", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0]);
    h_baseline_pt_nJets     = new TH2F("h_baseline_pt_nJets", "p_{T} vs NJet;p_{T} [GeV];NJet", n_pt_bins, &m_pt_bins[0], n_NJet_bins, &m_NJet_bins[0]);
    h_baseline_pt_dRjet     = new TH2F("h_baseline_pt_dRjet", "p_{T} vs dRjet;p_{T} [GeV];dRjet", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0]);
    //h_baseline_pt_Etmiss    = new TH2F("h_baseline_pt_Etmiss", "p_{T} vs Etmiss;p_{T} [GeV];Etmiss", n_pt_bins, &m_pt_bins[0], n_Etmiss_bins, &m_Etmiss_bins[0]);
    //h_baseline_pt_meff      = new TH2F("h_baseline_pt_meff", "p_{T} vs meff;p_{T} [GeV];meff", n_pt_bins, &m_pt_bins[0], n_meff_bins, &m_meff_bins[0]);

    h_baseline_pt_eta->Sumw2();
    h_baseline_pt_nJets->Sumw2();
    h_baseline_pt_dRjet->Sumw2();
    //h_baseline_pt_Etmiss->Sumw2();
    //h_baseline_pt_meff->Sumw2();

    wk()->addOutput(h_baseline_pt_eta);
    wk()->addOutput(h_baseline_pt_nJets);
    wk()->addOutput(h_baseline_pt_dRjet);
    //wk()->addOutput(h_baseline_pt_Etmiss);
    //wk()->addOutput(h_baseline_pt_meff);

    // 3-dim histograms
    h_baseline_pt_eta_mll   = new TH3F("h_baseline_pt_eta_mll", "Baseline 3-dim;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_pt_eta_dRjet = new TH3F("h_baseline_pt_eta_dRjet", "Baseline 3-dim;p_{T} [GeV];#eta;#Delta R(#ell, jet)", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_deltaR_bins, &m_deltaR_bins[0]);
    h_baseline_pt_dRjet_mll = new TH3F("h_baseline_pt_dRjet_mll", "Baseline 3-dim;p_{T} [GeV];#Delta R(#ell, jet);M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_eta_dRjet_mll = new TH3F("h_baseline_eta_dRjet_mll", "Baseline 3-dim;#eta;#Delta R(#ell, jet);M_{ll} [GeV]", n_eta_bins, &m_eta_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_mll_bins, &m_mll_bins[0]);

    h_baseline_pt_eta_mll->Sumw2();
    h_baseline_pt_eta_dRjet->Sumw2();
    h_baseline_pt_dRjet_mll->Sumw2();
    h_baseline_eta_dRjet_mll->Sumw2();

    wk()->addOutput(h_baseline_pt_eta_mll);
    wk()->addOutput(h_baseline_pt_eta_dRjet);
    wk()->addOutput(h_baseline_pt_dRjet_mll);
    wk()->addOutput(h_baseline_eta_dRjet_mll);

    h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll = new TH3F("h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll", "LooseAndBLayerLLH_to_MediumLLH_cut;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_CaloIso_cut_pt_eta_mll  = new TH3F("h_baseline_CaloIso_cut_pt_eta_mll", "CaloIso_cut;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_TrackIso_cut_pt_eta_mll = new TH3F("h_baseline_TrackIso_cut_pt_eta_mll", "TrackIso_cut;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_Iso_cut_pt_eta_mll      = new TH3F("h_baseline_Iso_cut_pt_eta_mll", "Iso_cut;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_z0_cut_pt_eta_mll       = new TH3F("h_baseline_z0_cut_pt_eta_mll", "z0_cut;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_baseline_sigd0_cut_pt_eta_mll    = new TH3F("h_baseline_sigd0_cut_pt_eta_mll", "sigd0_cut;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);

    h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll->Sumw2();
    h_baseline_CaloIso_cut_pt_eta_mll->Sumw2();
    h_baseline_TrackIso_cut_pt_eta_mll->Sumw2();
    h_baseline_Iso_cut_pt_eta_mll->Sumw2();
    h_baseline_z0_cut_pt_eta_mll->Sumw2();
    h_baseline_sigd0_cut_pt_eta_mll->Sumw2();

    wk()->addOutput(h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll);
    wk()->addOutput(h_baseline_CaloIso_cut_pt_eta_mll);
    wk()->addOutput(h_baseline_TrackIso_cut_pt_eta_mll);
    wk()->addOutput(h_baseline_Iso_cut_pt_eta_mll);
    wk()->addOutput(h_baseline_z0_cut_pt_eta_mll);
    wk()->addOutput(h_baseline_sigd0_cut_pt_eta_mll);

    h_baseline_pt_eta_AvgMu = new TH3F("h_baseline_pt_eta_AvgMu", "AvgMu;p_{T} [GeV];#eta;<#mu>", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_AvgMu_bins, &m_AvgMu_bins[0]);
    h_baseline_pt_dRjet_AvgMu = new TH3F("h_baseline_pt_dRjet_AvgMu", "AvgMu;p_{T} [GeV];#Delta R(#ell, jet);<#mu>", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_AvgMu_bins, &m_AvgMu_bins[0]);
    h_baseline_eta_dRjet_AvgMu = new TH3F("h_baseline_eta_dRjet_AvgMu", "AvgMu;#eta;#Delta R(#ell, jet);<#mu>", n_eta_bins, &m_eta_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_AvgMu_bins, &m_AvgMu_bins[0]);

    h_baseline_pt_eta_AvgMu->Sumw2();
    h_baseline_pt_dRjet_AvgMu->Sumw2();
    h_baseline_eta_dRjet_AvgMu->Sumw2();

    wk()->addOutput(h_baseline_pt_eta_AvgMu);
    wk()->addOutput(h_baseline_pt_dRjet_AvgMu);
    wk()->addOutput(h_baseline_eta_dRjet_AvgMu);

    // Signal histograms
    // 1-dim histograms
    h_signal_mll            = new TH1F("h_signal_mll", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_signal_pt             = new TH1F("h_signal_pt", "pt;p_{T} [GeV];Events;", n_pt_bins, &m_pt_bins[0]);
    h_signal_eta            = new TH1F("h_signal_eta", "eta;#eta;Events", n_eta_bins, &m_eta_bins[0]);
    //h_signal_d0pvtx         = new TH1F("h_signal_d0pvtx", "d0pvtx", 100, -0.15, 0.15);
    //h_signal_sigd0          = new TH1F("h_signal_sigd0", "sigd0", 100, 0, 10);
    //h_signal_d0err          = new TH1F("h_signal_d0err", "d0err", 100, -0.1, 0.1);
    //h_signal_z0sinTheta     = new TH1F("h_signal_z0SinTheta", "fabs(z0sinTheta)", 100, 0, 0.5);
    //h_signal_ptvarcone20    = new TH1F("h_signal_ptvarcone20", "ptvarcone20", 100, 0, 0.2);
    //h_signal_ptvarcone30    = new TH1F("h_signal_ptvarcone30", "ptvarcone30", 100, 0, 0.2);
    //h_signal_topoetcone20   = new TH1F("h_signal_topoetcone20", "topoetcone20", 100, 0, 0.2);
    h_signal_nJets          = new TH1F("h_signal_nJets", "NJet", n_NJet_bins, &m_NJet_bins[0]);
    h_signal_dRjet          = new TH1F("h_signal_dRjet", "dRjet", n_deltaR_bins, &m_deltaR_bins[0]);
    //h_signal_Etmiss         = new TH1F("h_signal_Etmiss", "Etmiss", n_Etmiss_bins, &m_Etmiss_bins[0]);
    //h_signal_meff           = new TH1F("h_signal_meff", "meff", n_meff_bins, &m_meff_bins[0]);

    h_signal_mll->Sumw2();
    h_signal_pt->Sumw2();
    h_signal_eta->Sumw2();
    //h_signal_d0pvtx->Sumw2();
    //h_signal_sigd0->Sumw2();
    //h_signal_d0err->Sumw2();
    //h_signal_z0sinTheta->Sumw2();
    //h_signal_ptvarcone20->Sumw2();
    //h_signal_ptvarcone30->Sumw2();
    //h_signal_topoetcone20->Sumw2();
    h_signal_nJets->Sumw2();
    h_signal_dRjet->Sumw2();
    //h_signal_Etmiss->Sumw2();
    //h_signal_meff->Sumw2();

    wk()->addOutput(h_signal_mll);
    wk()->addOutput(h_signal_pt);
    wk()->addOutput(h_signal_eta);
    //wk()->addOutput(h_signal_d0pvtx);
    //wk()->addOutput(h_signal_sigd0);
    //wk()->addOutput(h_signal_d0err);
    //wk()->addOutput(h_signal_z0sinTheta);
    //wk()->addOutput(h_signal_ptvarcone20);
    //wk()->addOutput(h_signal_ptvarcone30);
    //wk()->addOutput(h_signal_topoetcone20);
    wk()->addOutput(h_signal_nJets);
    wk()->addOutput(h_signal_dRjet);
    //wk()->addOutput(h_signal_Etmiss);
    //wk()->addOutput(h_signal_meff);
   
    // 2-dim histograms
    h_signal_pt_eta     = new TH2F("h_signal_pt_eta", "p_{T} vs #eta;p_{T} [GeV];#eta", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0]);
    h_signal_pt_nJets   = new TH2F("h_signal_pt_nJets", "p_{T} vs NJet;p_{T} [GeV];NJet", n_pt_bins, &m_pt_bins[0], n_NJet_bins, &m_NJet_bins[0]);
    h_signal_pt_dRjet   = new TH2F("h_signal_pt_dRjet", "p_{T} vs dRjet;p_{T} [GeV];dRjet", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0]);
    //h_signal_pt_Etmiss  = new TH2F("h_signal_pt_Etmiss", "p_{T} vs Etmiss;p_{T} [GeV];Etmiss", n_pt_bins, &m_pt_bins[0], n_Etmiss_bins, &m_Etmiss_bins[0]);
    //h_signal_pt_meff    = new TH2F("h_signal_pt_meff", "p_{T} vs meff;p_{T} [GeV];meff", n_pt_bins, &m_pt_bins[0], n_meff_bins, &m_meff_bins[0]);

    h_signal_pt_eta->Sumw2();
    h_signal_pt_nJets->Sumw2();
    h_signal_pt_dRjet->Sumw2();
    //h_signal_pt_Etmiss->Sumw2();
    //h_signal_pt_meff->Sumw2();

    wk()->addOutput(h_signal_pt_eta);
    wk()->addOutput(h_signal_pt_nJets);
    wk()->addOutput(h_signal_pt_dRjet);
    //wk()->addOutput(h_signal_pt_Etmiss);
    //wk()->addOutput(h_signal_pt_meff);

    // 3-dim histograms
    h_signal_pt_eta_mll = new TH3F("h_signal_pt_eta_mll", "Signal 3-dim;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_signal_pt_eta_dRjet = new TH3F("h_signal_pt_eta_dRjet", "Signal 3-dim;p_{T} [GeV];#eta;#Delta R(#ell, jet)", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_deltaR_bins, &m_deltaR_bins[0]);
    h_signal_pt_dRjet_mll = new TH3F("h_signal_pt_dRjet_mll", "Signal 3-dim;p_{T} [GeV];#Delta R(#ell, jet);M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_signal_eta_dRjet_mll = new TH3F("h_signal_eta_dRjet_mll", "Signal 3-dim;#eta;#Delta R(#ell, jet);M_{ll} [GeV]", n_eta_bins, &m_eta_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_mll_bins, &m_mll_bins[0]);

    h_signal_pt_eta_mll->Sumw2();
    h_signal_pt_eta_dRjet->Sumw2();
    h_signal_pt_dRjet_mll->Sumw2();
    h_signal_eta_dRjet_mll->Sumw2();

    wk()->addOutput(h_signal_pt_eta_mll);
    wk()->addOutput(h_signal_pt_eta_dRjet);
    wk()->addOutput(h_signal_pt_dRjet_mll);
    wk()->addOutput(h_signal_eta_dRjet_mll);

    // Cut efficiency histograms
    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH = new TH1F("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH", "LooseAndBLayerLLH_to_MediumLLH;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    h_cut_eff_CaloIso   = new TH1F("h_cut_eff_CaloIso", "CaloIso;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    h_cut_eff_TrackIso  = new TH1F("h_cut_eff_TrackIso", "TrackIso;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    h_cut_eff_Iso       = new TH1F("h_cut_eff_Iso", "Iso;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    h_cut_eff_z0        = new TH1F("h_cut_eff_z0", "z0;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    h_cut_eff_sigd0     = new TH1F("h_cut_eff_sigd0", "sigd0;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);

    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Sumw2();
    h_cut_eff_CaloIso->Sumw2();
    h_cut_eff_TrackIso->Sumw2();
    h_cut_eff_Iso->Sumw2();
    h_cut_eff_z0->Sumw2();
    h_cut_eff_sigd0->Sumw2();

    wk()->addOutput(h_cut_eff_LooseAndBLayerLLH_to_MediumLLH);
    wk()->addOutput(h_cut_eff_CaloIso);
    wk()->addOutput(h_cut_eff_TrackIso);
    wk()->addOutput(h_cut_eff_Iso);
    wk()->addOutput(h_cut_eff_z0);
    wk()->addOutput(h_cut_eff_sigd0);

    h_signal_pt_eta_AvgMu = new TH3F("h_signal_pt_eta_AvgMu", "AvgMu;p_{T} [GeV];#eta;<#mu>", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_AvgMu_bins, &m_AvgMu_bins[0]);
    h_signal_pt_dRjet_AvgMu = new TH3F("h_signal_pt_dRjet_AvgMu", "AvgMu;p_{T} [GeV];#Delta R(#ell, jet);<#mu>", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_AvgMu_bins, &m_AvgMu_bins[0]);
    h_signal_eta_dRjet_AvgMu = new TH3F("h_signal_eta_dRjet_AvgMu", "AvgMu;#eta;#Delta R(#ell, jet);<#mu>", n_eta_bins, &m_eta_bins[0], n_deltaR_bins, &m_deltaR_bins[0], n_AvgMu_bins, &m_AvgMu_bins[0]);

    h_signal_pt_eta_AvgMu->Sumw2();
    h_signal_pt_dRjet_AvgMu->Sumw2();
    h_signal_eta_dRjet_AvgMu->Sumw2();

    wk()->addOutput(h_signal_pt_eta_AvgMu);
    wk()->addOutput(h_signal_pt_dRjet_AvgMu);
    wk()->addOutput(h_signal_eta_dRjet_AvgMu);

    // Background template histograms
    h_bkg_template_fail_truth           = new TH1F("h_bkg_template_fail_truth", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_only         = new TH1F("h_bkg_template_fail_id_only", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_CaloIso_only    = new TH1F("h_bkg_template_fail_CaloIso_only", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_TrackIso_only   = new TH1F("h_bkg_template_fail_TrackIso_only", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_CaloIso  = new TH1F("h_bkg_template_fail_id_and_CaloIso", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_TrackIso = new TH1F("h_bkg_template_fail_id_and_TrackIso", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso         = new TH1F("h_bkg_template_fail_id_and_CaloIso_and_TrackIso", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight   = new TH1F("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_CaloIso_and_TrackIso                = new TH1F("h_bkg_template_fail_CaloIso_and_TrackIso", "mll;M_{ll} [GeV];Events;", n_mll_bins, &m_mll_bins[0]);

    h_bkg_template_fail_truth->Sumw2();
    h_bkg_template_fail_id_only->Sumw2();
    h_bkg_template_fail_CaloIso_only->Sumw2();
    h_bkg_template_fail_TrackIso_only->Sumw2();
    h_bkg_template_fail_id_and_CaloIso->Sumw2();
    h_bkg_template_fail_id_and_TrackIso->Sumw2();
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso->Sumw2();
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight->Sumw2();
    h_bkg_template_fail_CaloIso_and_TrackIso->Sumw2();

    wk()->addOutput(h_bkg_template_fail_truth);
    wk()->addOutput(h_bkg_template_fail_id_only);
    wk()->addOutput(h_bkg_template_fail_CaloIso_only);
    wk()->addOutput(h_bkg_template_fail_TrackIso_only);
    wk()->addOutput(h_bkg_template_fail_id_and_CaloIso);
    wk()->addOutput(h_bkg_template_fail_id_and_TrackIso);
    wk()->addOutput(h_bkg_template_fail_id_and_CaloIso_and_TrackIso);
    wk()->addOutput(h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight);
    wk()->addOutput(h_bkg_template_fail_CaloIso_and_TrackIso);

    h_bkg_template_fail_truth_pt_eta_mll = new TH3F("h_bkg_template_fail_truth_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_only_pt_eta_mll = new TH3F("h_bkg_template_fail_id_only_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_CaloIso_only_pt_eta_mll = new TH3F("h_bkg_template_fail_CaloIso_only_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_TrackIso_only_pt_eta_mll = new TH3F("h_bkg_template_fail_TrackIso_only_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_CaloIso_pt_eta_mll = new TH3F("h_bkg_template_fail_id_and_CaloIso_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_TrackIso_pt_eta_mll = new TH3F("h_bkg_template_fail_id_and_TrackIso_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll = new TH3F("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll = new TH3F("h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);
    h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll = new TH3F("h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll", "mll;p_{T} [GeV];#eta;M_{ll} [GeV]", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0], n_mll_bins, &m_mll_bins[0]);

    h_bkg_template_fail_truth_pt_eta_mll->Sumw2();
    h_bkg_template_fail_id_only_pt_eta_mll->Sumw2();
    h_bkg_template_fail_CaloIso_only_pt_eta_mll->Sumw2();
    h_bkg_template_fail_TrackIso_only_pt_eta_mll->Sumw2();
    h_bkg_template_fail_id_and_CaloIso_pt_eta_mll->Sumw2();
    h_bkg_template_fail_id_and_TrackIso_pt_eta_mll->Sumw2();
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->Sumw2();
    h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll->Sumw2();
    h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll->Sumw2();

    wk()->addOutput(h_bkg_template_fail_truth_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_id_only_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_CaloIso_only_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_TrackIso_only_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_id_and_CaloIso_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_id_and_TrackIso_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll);
    wk()->addOutput(h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll);

    // Real lepton efficiency histograms
    // 1-dim histograms
    h_eff_pt        = new TH1F("h_eff_pt", "eff_pt;p_{T} [GeV];Efficiency", n_pt_bins, &m_pt_bins[0]);
    h_eff_eta       = new TH1F("h_eff_eta", "eff_eta;|#eta|;Efficiency", n_eta_bins, &m_eta_bins[0]);
    h_eff_nJets     = new TH1F("h_eff_nJets", "eff_nJet;NJet;Efficiency", n_NJet_bins, &m_NJet_bins[0]);
    h_eff_dRjet     = new TH1F("h_eff_dRjet", "eff_dRjet;dRjet;Efficiency", n_deltaR_bins, &m_deltaR_bins[0]);
    //h_eff_Etmiss    = new TH1F("h_eff_Etmiss", "eff_Etmiss;Etmiss;Efficiency", n_Etmiss_bins, &m_Etmiss_bins[0]);
    //h_eff_meff      = new TH1F("h_eff_meff", "eff_meff;meff;Efficiency", n_meff_bins, &m_meff_bins[0]);

    h_eff_pt->Sumw2();
    h_eff_eta->Sumw2();
    h_eff_nJets->Sumw2();
    h_eff_dRjet->Sumw2();

    wk()->addOutput(h_eff_pt);
    wk()->addOutput(h_eff_eta);
    wk()->addOutput(h_eff_nJets);
    wk()->addOutput(h_eff_dRjet);
    //wk()->addOutput(h_eff_Etmiss);
    //wk()->addOutput(h_eff_meff);

    // 2-dim histograms
    //h_2d_eff_pt_eta     = new TH2F("h_2d_eff_pt_eta", "eff_pt_eta;p_{T} [GeV];#eta;Efficiency", n_pt_bins, &m_pt_bins[0], n_eta_bins, &m_eta_bins[0]);
    //h_2d_eff_pt_nJets   = new TH2F("h_2d_eff_pt_nJets", "eff_pt_nJets;p_{T} [GeV];NJet;Efficiency", n_pt_bins, &m_pt_bins[0], n_NJet_bins, &m_NJet_bins[0]);
    //h_2d_eff_pt_dRjet   = new TH2F("h_2d_eff_pt_dRjet", "eff_pt_dRjet;p_{T} [GeV];dRjet;Efficiency", n_pt_bins, &m_pt_bins[0], n_deltaR_bins, &m_deltaR_bins[0]);
    //h_2d_eff_pt_Etmiss  = new TH2F("h_2d_eff_pt_Etmiss", "eff_pt_Etmiss;p_{T} [GeV];Etmiss;Efficiency", n_pt_bins, &m_pt_bins[0], n_Etmiss_bins, &m_Etmiss_bins[0]);
    //h_2d_eff_pt_meff    = new TH2F("h_2d_eff_pt_meff", "eff_pt_meff;p_{T} [GeV];meff;Efficiency", n_pt_bins, &m_pt_bins[0], n_meff_bins, &m_meff_bins[0]);

    //h_2d_eff_pt_eta->Sumw2();
    //h_2d_eff_pt_nJets->Sumw2();

    //wk()->addOutput(h_2d_eff_pt_eta);
    //wk()->addOutput(h_2d_eff_pt_nJets);
    //wk()->addOutput(h_2d_eff_pt_dRjet);
    //wk()->addOutput(h_2d_eff_pt_Etmiss);
    //wk()->addOutput(h_2d_eff_pt_meff);

    return EL::StatusCode::SUCCESS;
}




EL::StatusCode ytRealLeptonsEfficiency::fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed

    //const char * function_name = "fileExecute()";
    //Info(function_name, "Function calls");

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode ytRealLeptonsEfficiency::changeInput (bool firstFile)
{
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.

    //const char * function_name = "changeInput()";
    //Info(function_name, "Function calls");

    //Init (wk()->tree());
    //wk()->tree();
    return EL::StatusCode::SUCCESS;
}




EL::StatusCode ytRealLeptonsEfficiency::initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.

    //const char * function_name = "initialize()";
    //Info(function_name, "Function calls");

    // Set object pointer
    Mu_eta = 0;
    Mu_phi = 0;
    Mu_pT = 0;
    Mu_SFw = 0;
    Mu_IsoSFw = 0;
    Mu_charge = 0;
    Mu_d0pvtx = 0;
    Mu_sigd0 = 0;
    Mu_z0sinTheta = 0;
    Mu_isBad = 0;
    Mu_passOR = 0;
    Mu_isTight = 0;
    Mu_isCosmic = 0;
    Mu_type = 0;
    Mu_origin = 0;
    Mu_ptcone20 = 0;
    Mu_ptcone30 = 0;
    Mu_ptcone40 = 0;
    Mu_ptvarcone20 = 0;
    Mu_ptvarcone30 = 0;
    Mu_ptvarcone40 = 0;
    Mu_topoetcone20 = 0;
    Mu_topoetcone30 = 0;
    Mu_topoetcone40 = 0;
    Mu_passIsoLooseTO = 0;
    Mu_passIsoLoose = 0;
    Mu_passIsoTight = 0;
    Mu_passIsoGrad = 0;
    Mu_passIsoGradCustomTight = 0;
    Mu_passIsoGradCustom = 0;
    Mu_passIsoGradLoose = 0;
    Mu_trigMatch_mu26_ivarmedium = 0;
    Mu_trigMatch_mu20_iloose_L1MU15 = 0;
    Mu_trigMatch_mu40 = 0;
    Mu_trigMatch_mu50 = 0;
    Mu_trigMatch_mu8noL1 = 0;
    Mu_trigMatch_mu14 = 0;
    Mu_trigMatch_mu18 = 0;
    Mu_trigMatch_mu18_mu8noL1 = 0;
    Mu_trigMatch_e17_lhloose_mu14 = 0;
    Mu_trigMatch_e17_lhloose_nod0_mu14 = 0;
    Mu_trigMatch_mu20_mu8noL1 = 0;
    Mu_trigMatch_mu22_mu8noL1 = 0;
    Mu_trigMatch_mu24_iloose = 0;
    Mu_trigMatch_mu24_ivarloose = 0;
    Mu_trigMatch_mu24_iloose_L1MU15 = 0;
    Mu_trigMatch_mu24_ivarloose_L1MU15 = 0;
    Mu_trigMatchPair_mu18_mu8noL1 = 0;
    Mu_trigMatchPair_mu20_mu8noL1 = 0;
    Mu_trigMatchPair_mu22_mu8noL1 = 0;
    El_eta = 0;
    El_etaclus = 0;
    El_phi = 0;
    El_pT = 0;
    El_E = 0;
    El_charge = 0;
    El_sigd0 = 0;
    El_z0sinTheta = 0;
    El_d0pvtx = 0;
    El_passOR = 0;
    El_SFwMediumLH = 0;
    El_IsoSFwMediumLH = 0;
    El_SFwTightLH = 0;
    El_SFwLooseAndBLayerLH = 0;
    El_SFwTrigMediumLH_e12_lhloose_L1EM10VH = 0;
    El_SFwTrigMediumLH_e17_lhloose = 0;
    El_SFwTrigMediumLH_single = 0;
    El_SFwTrigLooseAndBLayerLH_e12_lhloose_L1EM10VH = 0;
    El_isLooseAndBLayerLH = 0;
    El_isMediumLH = 0;
    El_isTightLH = 0;
    El_nBLayerHits = 0;
    El_expectBLayerHit = 0;
    El_truthType = 0;
    El_truthOrigin = 0;
    El_bkgTruthType = 0;
    El_bkgTruthOrigin = 0;
    El_bkgMotherPdgId = 0;
    El_firstEgMotherTruthType = 0;
    El_firstEgMotherTruthOrigin = 0;
    El_firstEgMotherPdgId = 0;
    El_chFlip = 0;
    El_ptcone20 = 0;
    El_ptcone30 = 0;
    El_ptcone40 = 0;
    El_ptvarcone20 = 0;
    El_ptvarcone30 = 0;
    El_ptvarcone40 = 0;
    El_topoetcone20 = 0;
    El_topoetcone30 = 0;
    El_topoetcone40 = 0;
    El_passIsoLooseTO = 0;
    El_passIsoLoose = 0;
    El_passIsoTight = 0;
    El_passIsoGrad = 0;
    El_passIsoGradCustomTight = 0;
    El_passIsoGradCustom = 0;
    El_passIsoGradLoose = 0;
    El_trigMatch_e12_lhloose_L1EM10VH = 0;
    El_trigMatch_e17_lhloose = 0;
    El_trigMatch_e60_lhmedium = 0;
    El_trigMatch_e24_lhmedium_L1EM20VH = 0;
    El_trigMatch_e24_lhmedium_iloose_L1EM20VH = 0;
    El_trigMatch_e26_lhtight_nod0_ivarloose = 0;
    El_trigMatch_2e12_lhloose_L12EM10VH = 0;
    El_trigMatch_2e15_lhloose_L12EM10VH = 0;
    El_trigMatch_2e15_lhvloose_L12EM13VH = 0;
    El_trigMatch_2e15_lhvloose_nod0_L12EM13VH = 0;
    El_trigMatch_2e17_lhvloose_nod0 = 0;
    El_trigMatch_e17_lhloose_mu14 = 0;
    El_trigMatch_e17_lhloose_nod0_mu14 = 0;
    El_trigMatch_e24_lhmedium_nod0_ivarloose = 0;
    El_trigMatch_e24_lhtight_nod0_ivarloose = 0;
    El_trigMatch_e60_lhmedium_nod0 = 0;
    El_passChargeFlipTaggerBDTmedium = 0;
    El_passChargeFlipTaggerBDTloose = 0;
    Jet_eta = 0;
    Jet_phi = 0;
    Jet_pT = 0;
    Jet_E = 0;
    Jet_quality = 0;
    Jet_JVT = 0;
    Jet_JVTsf = 0;
    Jet_MV2c20 = 0;
    Jet_MV2c10 = 0;
    Jet_SFw = 0;
    Jet_ConeTruthLabel = 0;
    Jet_PartonTruthLabel = 0;
    Jet_HadronConeExclTruthLabel = 0;
    Jet_deltaR = 0;
    Jet_nTrk = 0;
    Jet_passOR = 0;
    TruthJet_eta = 0;
    TruthJet_phi = 0;
    TruthJet_pT = 0;
    TruthJet_E = 0;
    TruthJet_ConeTruthLabel = 0;
    TruthJet_PartonTruthLabel = 0;
    TruthL_eta = 0;
    TruthL_phi = 0;
    TruthL_pT = 0;
    TruthL_id = 0;
    El_isBaseline = 0;
    El_isSignal = 0;
    El_isZTag = 0;
    El_isZProbe = 0;
    El_Tag_index = 0;
    El_Probe_index = 0;
    El_ZTandP_mll = 0;
    El_isZProbe_TagTriggerMatched = 0;
    El_Tag_Trigger_Matched = 0;
    El_Probe_Trigger_Matched = 0;
    El_DR_closest_Jet = 0;
    Mu_isBaseline = 0;
    Mu_isSignal = 0;
    Mu_isZTag = 0;
    Mu_isZProbe = 0;
    Mu_Tag_index = 0;
    Mu_Probe_index = 0;
    Mu_ZTandP_mll = 0;
    Mu_isZProbe_TagTriggerMatched = 0;
    Mu_Tag_Trigger_Matched = 0;
    Mu_Probe_Trigger_Matched = 0;
    Mu_DR_closest_Jet = 0;

    // Set branch addresses and branch pointers
    fChain = wk()->tree();
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("HLT_e24_lhmedium_nod0_ivarloose", &HLT_e24_lhmedium_nod0_ivarloose, &b_HLT_e24_lhmedium_nod0_ivarloose);
    fChain->SetBranchAddress("HLT_e24_lhtight_nod0_ivarloose", &HLT_e24_lhtight_nod0_ivarloose, &b_HLT_e24_lhtight_nod0_ivarloose);
    fChain->SetBranchAddress("HLT_e24_lhmedium_nod0_L1EM20VH", &HLT_e24_lhmedium_nod0_L1EM20VH, &b_HLT_e24_lhmedium_nod0_L1EM20VH);
    fChain->SetBranchAddress("HLT_e26_lhtight_iloose", &HLT_e26_lhtight_iloose, &b_HLT_e26_lhtight_iloose);
    fChain->SetBranchAddress("HLT_e26_lhtight_ivarloose", &HLT_e26_lhtight_ivarloose, &b_HLT_e26_lhtight_ivarloose);
    fChain->SetBranchAddress("HLT_e26_lhtight_nod0_iloose", &HLT_e26_lhtight_nod0_iloose, &b_HLT_e26_lhtight_nod0_iloose);
    fChain->SetBranchAddress("HLT_e26_lhtight_nod0_ivarloose", &HLT_e26_lhtight_nod0_ivarloose, &b_HLT_e26_lhtight_nod0_ivarloose);
    fChain->SetBranchAddress("HLT_e60_lhmedium", &HLT_e60_lhmedium, &b_HLT_e60_lhmedium);
    fChain->SetBranchAddress("HLT_e60_lhmedium_nod0", &HLT_e60_lhmedium_nod0, &b_HLT_e60_lhmedium_nod0);
    fChain->SetBranchAddress("HLT_e120_lhloose_nod0", &HLT_e120_lhloose_nod0, &b_HLT_e120_lhloose_nod0);
    fChain->SetBranchAddress("HLT_e140_lhloose_nod0", &HLT_e140_lhloose_nod0, &b_HLT_e140_lhloose_nod0);
    fChain->SetBranchAddress("HLT_2e17_lhvloose_nod0", &HLT_2e17_lhvloose_nod0, &b_HLT_2e17_lhvloose_nod0);
    fChain->SetBranchAddress("HLT_2e15_lhvloose_nod0_L12EM13VH", &HLT_2e15_lhvloose_nod0_L12EM13VH, &b_HLT_2e15_lhvloose_nod0_L12EM13VH);
    fChain->SetBranchAddress("HLT_e24_lhmedium_e9_lhmedium", &HLT_e24_lhmedium_e9_lhmedium, &b_HLT_e24_lhmedium_e9_lhmedium);
    fChain->SetBranchAddress("HLT_e24_lhmedium_L1EM20VH", &HLT_e24_lhmedium_L1EM20VH, &b_HLT_e24_lhmedium_L1EM20VH);
    fChain->SetBranchAddress("HLT_e24_lhmedium_iloose_L1EM20VH", &HLT_e24_lhmedium_iloose_L1EM20VH, &b_HLT_e24_lhmedium_iloose_L1EM20VH);
    fChain->SetBranchAddress("HLT_e12_lhvloose_L12EM10VH", &HLT_e12_lhvloose_L12EM10VH, &b_HLT_e12_lhvloose_L12EM10VH);
    fChain->SetBranchAddress("HLT_e17_lhloose_2e9_lhloose", &HLT_e17_lhloose_2e9_lhloose, &b_HLT_e17_lhloose_2e9_lhloose);
    fChain->SetBranchAddress("HLT_mu24_ivarmedium", &HLT_mu24_ivarmedium, &b_HLT_mu24_ivarmedium);
    fChain->SetBranchAddress("HLT_mu24_imedium", &HLT_mu24_imedium, &b_HLT_mu24_imedium);
    fChain->SetBranchAddress("HLT_mu24_ivarloose", &HLT_mu24_ivarloose, &b_HLT_mu24_ivarloose);
    fChain->SetBranchAddress("HLT_mu24_iloose", &HLT_mu24_iloose, &b_HLT_mu24_iloose);
    fChain->SetBranchAddress("HLT_mu26_ivarmedium", &HLT_mu26_ivarmedium, &b_HLT_mu26_ivarmedium);
    fChain->SetBranchAddress("HLT_mu20_ivarmedium_L1MU15", &HLT_mu20_ivarmedium_L1MU15, &b_HLT_mu20_ivarmedium_L1MU15);
    fChain->SetBranchAddress("HLT_mu20_imedium_L1MU15", &HLT_mu20_imedium_L1MU15, &b_HLT_mu20_imedium_L1MU15);
    fChain->SetBranchAddress("HLT_mu20_ivarloose_L1MU15", &HLT_mu20_ivarloose_L1MU15, &b_HLT_mu20_ivarloose_L1MU15);
    fChain->SetBranchAddress("HLT_mu20_iloose_L1MU15", &HLT_mu20_iloose_L1MU15, &b_HLT_mu20_iloose_L1MU15);
    fChain->SetBranchAddress("HLT_mu20_L1MU15", &HLT_mu20_L1MU15, &b_HLT_mu20_L1MU15);
    fChain->SetBranchAddress("HLT_mu20_mu8noL1", &HLT_mu20_mu8noL1, &b_HLT_mu20_mu8noL1);
    fChain->SetBranchAddress("HLT_mu22_mu8noL1", &HLT_mu22_mu8noL1, &b_HLT_mu22_mu8noL1);
    fChain->SetBranchAddress("HLT_mu20_2mu4noL1", &HLT_mu20_2mu4noL1, &b_HLT_mu20_2mu4noL1);
    fChain->SetBranchAddress("HLT_mu22_2mu4noL1", &HLT_mu22_2mu4noL1, &b_HLT_mu22_2mu4noL1);
    fChain->SetBranchAddress("HLT_mu40", &HLT_mu40, &b_HLT_mu40);
    fChain->SetBranchAddress("HLT_mu50", &HLT_mu50, &b_HLT_mu50);
    fChain->SetBranchAddress("HLT_2mu10", &HLT_2mu10, &b_HLT_2mu10);
    fChain->SetBranchAddress("HLT_2mu10_nomucomb", &HLT_2mu10_nomucomb, &b_HLT_2mu10_nomucomb);
    fChain->SetBranchAddress("HLT_2mu14", &HLT_2mu14, &b_HLT_2mu14);
    fChain->SetBranchAddress("HLT_2mu14_nomucomb", &HLT_2mu14_nomucomb, &b_HLT_2mu14_nomucomb);
    fChain->SetBranchAddress("HLT_3mu6", &HLT_3mu6, &b_HLT_3mu6);
    fChain->SetBranchAddress("HLT_3mu6_msonly", &HLT_3mu6_msonly, &b_HLT_3mu6_msonly);
    fChain->SetBranchAddress("HLT_xe100_L1XE50", &HLT_xe100_L1XE50, &b_HLT_xe100_L1XE50);
    fChain->SetBranchAddress("HLT_xe80_mht_L1XE50", &HLT_xe80_mht_L1XE50, &b_HLT_xe80_mht_L1XE50);
    fChain->SetBranchAddress("HLT_xe90_mht_L1XE50", &HLT_xe90_mht_L1XE50, &b_HLT_xe90_mht_L1XE50);
    fChain->SetBranchAddress("HLT_xe100_mht_L1XE50", &HLT_xe100_mht_L1XE50, &b_HLT_xe100_mht_L1XE50);
    fChain->SetBranchAddress("HLT_xe100_tc_em_L1XE50", &HLT_xe100_tc_em_L1XE50, &b_HLT_xe100_tc_em_L1XE50);
    fChain->SetBranchAddress("HLT_xe110_mht_L1XE50", &HLT_xe110_mht_L1XE50, &b_HLT_xe110_mht_L1XE50);
    fChain->SetBranchAddress("HLT_xe80_tc_lcw_L1XE50", &HLT_xe80_tc_lcw_L1XE50, &b_HLT_xe80_tc_lcw_L1XE50);
    fChain->SetBranchAddress("HLT_xe90_tc_lcw_L1XE50", &HLT_xe90_tc_lcw_L1XE50, &b_HLT_xe90_tc_lcw_L1XE50);
    fChain->SetBranchAddress("HLT_xe80_tc_lcw_wEFMu_L1XE50", &HLT_xe80_tc_lcw_wEFMu_L1XE50, &b_HLT_xe80_tc_lcw_wEFMu_L1XE50);
    fChain->SetBranchAddress("HLT_e7_lhmedium_mu24", &HLT_e7_lhmedium_mu24, &b_HLT_e7_lhmedium_mu24);
    fChain->SetBranchAddress("HLT_e7_lhmedium_nod0_mu24", &HLT_e7_lhmedium_nod0_mu24, &b_HLT_e7_lhmedium_nod0_mu24);
    fChain->SetBranchAddress("HLT_e17_lhloose_nod0_mu14", &HLT_e17_lhloose_nod0_mu14, &b_HLT_e17_lhloose_nod0_mu14);
    fChain->SetBranchAddress("HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1", &HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1, &b_HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1);
    fChain->SetBranchAddress("HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1", &HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1, &b_HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1);
    fChain->SetBranchAddress("HLT_2e12_lhloose_L12EM10VH", &HLT_2e12_lhloose_L12EM10VH, &b_HLT_2e12_lhloose_L12EM10VH);
    fChain->SetBranchAddress("HLT_e17_lhloose_mu14", &HLT_e17_lhloose_mu14, &b_HLT_e17_lhloose_mu14);
    fChain->SetBranchAddress("HLT_mu18_mu8noL1", &HLT_mu18_mu8noL1, &b_HLT_mu18_mu8noL1);
    fChain->SetBranchAddress("HLT_xe70", &HLT_xe70, &b_HLT_xe70);
    fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
    fChain->SetBranchAddress("ChannelNumber", &ChannelNumber, &b_ChannelNumber);
    fChain->SetBranchAddress("AvgMu", &AvgMu, &b_AvgMu);
    fChain->SetBranchAddress("EventWeight", &EventWeight, &b_EventWeight);
    fChain->SetBranchAddress("PRWWeight", &PRWWeight, &b_PRWWeight);
    fChain->SetBranchAddress("PRWrandomRunNumber", &PRWrandomRunNumber, &b_PRWrandomRunNumber);
    fChain->SetBranchAddress("TriggerDileptonSF", &TriggerDileptonSF, &b_TriggerDileptonSF);
    fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
    fChain->SetBranchAddress("LB", &LB, &b_LB);
    fChain->SetBranchAddress("RunNb", &RunNb, &b_RunNb);
    fChain->SetBranchAddress("DetError", &DetError, &b_DetError);
    fChain->SetBranchAddress("NMu", &NMu, &b_NMu);
    fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
    fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
    fChain->SetBranchAddress("Mu_pT", &Mu_pT, &b_Mu_pT);
    fChain->SetBranchAddress("Mu_SFw", &Mu_SFw, &b_Mu_SFw);
    fChain->SetBranchAddress("Mu_IsoSFw", &Mu_IsoSFw, &b_Mu_IsoSFw);
    fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
    fChain->SetBranchAddress("Mu_d0pvtx", &Mu_d0pvtx, &b_Mu_d0pvtx);
    fChain->SetBranchAddress("Mu_sigd0", &Mu_sigd0, &b_Mu_sigd0);
    fChain->SetBranchAddress("Mu_z0sinTheta", &Mu_z0sinTheta, &b_Mu_z0sinTheta);
    fChain->SetBranchAddress("Mu_isBad", &Mu_isBad, &b_Mu_isBad);
    fChain->SetBranchAddress("Mu_passOR", &Mu_passOR, &b_Mu_passOR);
    fChain->SetBranchAddress("Mu_isTight", &Mu_isTight, &b_Mu_isTight);
    fChain->SetBranchAddress("Mu_isCosmic", &Mu_isCosmic, &b_Mu_isCosmic);
    fChain->SetBranchAddress("Mu_type", &Mu_type, &b_Mu_type);
    fChain->SetBranchAddress("Mu_origin", &Mu_origin, &b_Mu_origin);
    fChain->SetBranchAddress("Mu_ptcone20", &Mu_ptcone20, &b_Mu_ptcone20);
    fChain->SetBranchAddress("Mu_ptcone30", &Mu_ptcone30, &b_Mu_ptcone30);
    fChain->SetBranchAddress("Mu_ptcone40", &Mu_ptcone40, &b_Mu_ptcone40);
    fChain->SetBranchAddress("Mu_ptvarcone20", &Mu_ptvarcone20, &b_Mu_ptvarcone20);
    fChain->SetBranchAddress("Mu_ptvarcone30", &Mu_ptvarcone30, &b_Mu_ptvarcone30);
    fChain->SetBranchAddress("Mu_ptvarcone40", &Mu_ptvarcone40, &b_Mu_ptvarcone40);
    fChain->SetBranchAddress("Mu_topoetcone20", &Mu_topoetcone20, &b_Mu_topoetcone20);
    fChain->SetBranchAddress("Mu_topoetcone30", &Mu_topoetcone30, &b_Mu_topoetcone30);
    fChain->SetBranchAddress("Mu_topoetcone40", &Mu_topoetcone40, &b_Mu_topoetcone40);
    fChain->SetBranchAddress("Mu_passIsoLooseTO", &Mu_passIsoLooseTO, &b_Mu_passIsoLooseTO);
    fChain->SetBranchAddress("Mu_passIsoLoose", &Mu_passIsoLoose, &b_Mu_passIsoLoose);
    fChain->SetBranchAddress("Mu_passIsoTight", &Mu_passIsoTight, &b_Mu_passIsoTight);
    fChain->SetBranchAddress("Mu_passIsoGrad", &Mu_passIsoGrad, &b_Mu_passIsoGrad);
    fChain->SetBranchAddress("Mu_passIsoGradCustomTight", &Mu_passIsoGradCustomTight, &b_Mu_passIsoGradCustomTight);
    fChain->SetBranchAddress("Mu_passIsoGradCustom", &Mu_passIsoGradCustom, &b_Mu_passIsoGradCustom);
    fChain->SetBranchAddress("Mu_passIsoGradLoose", &Mu_passIsoGradLoose, &b_Mu_passIsoGradLoose);
    fChain->SetBranchAddress("MuTrigSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50", &MuTrigSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50, &b_MuTrigSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50);
    fChain->SetBranchAddress("Mu_trigMatch_mu26_ivarmedium", &Mu_trigMatch_mu26_ivarmedium, &b_Mu_trigMatch_mu26_ivarmedium);
    fChain->SetBranchAddress("Mu_trigMatch_mu20_iloose_L1MU15", &Mu_trigMatch_mu20_iloose_L1MU15, &b_Mu_trigMatch_mu20_iloose_L1MU15);
    fChain->SetBranchAddress("Mu_trigMatch_mu40", &Mu_trigMatch_mu40, &b_Mu_trigMatch_mu40);
    fChain->SetBranchAddress("Mu_trigMatch_mu50", &Mu_trigMatch_mu50, &b_Mu_trigMatch_mu50);
    fChain->SetBranchAddress("Mu_trigMatch_mu8noL1", &Mu_trigMatch_mu8noL1, &b_Mu_trigMatch_mu8noL1);
    fChain->SetBranchAddress("Mu_trigMatch_mu14", &Mu_trigMatch_mu14, &b_Mu_trigMatch_mu14);
    fChain->SetBranchAddress("Mu_trigMatch_mu18", &Mu_trigMatch_mu18, &b_Mu_trigMatch_mu18);
    fChain->SetBranchAddress("Mu_trigMatch_mu18_mu8noL1", &Mu_trigMatch_mu18_mu8noL1, &b_Mu_trigMatch_mu18_mu8noL1);
    fChain->SetBranchAddress("Mu_trigMatch_e17_lhloose_mu14", &Mu_trigMatch_e17_lhloose_mu14, &b_Mu_trigMatch_e17_lhloose_mu14);
    fChain->SetBranchAddress("Mu_trigMatch_e17_lhloose_nod0_mu14", &Mu_trigMatch_e17_lhloose_nod0_mu14, &b_Mu_trigMatch_e17_lhloose_nod0_mu14);
    fChain->SetBranchAddress("Mu_trigMatch_mu20_mu8noL1", &Mu_trigMatch_mu20_mu8noL1, &b_Mu_trigMatch_mu20_mu8noL1);
    fChain->SetBranchAddress("Mu_trigMatch_mu22_mu8noL1", &Mu_trigMatch_mu22_mu8noL1, &b_Mu_trigMatch_mu22_mu8noL1);
    fChain->SetBranchAddress("Mu_trigMatch_mu24_iloose", &Mu_trigMatch_mu24_iloose, &b_Mu_trigMatch_mu24_iloose);
    fChain->SetBranchAddress("Mu_trigMatch_mu24_ivarloose", &Mu_trigMatch_mu24_ivarloose, &b_Mu_trigMatch_mu24_ivarloose);
    fChain->SetBranchAddress("Mu_trigMatch_mu24_iloose_L1MU15", &Mu_trigMatch_mu24_iloose_L1MU15, &b_Mu_trigMatch_mu24_iloose_L1MU15);
    fChain->SetBranchAddress("Mu_trigMatch_mu24_ivarloose_L1MU15", &Mu_trigMatch_mu24_ivarloose_L1MU15, &b_Mu_trigMatch_mu24_ivarloose_L1MU15);
    fChain->SetBranchAddress("Mu_trigMatchPair_mu18_mu8noL1", &Mu_trigMatchPair_mu18_mu8noL1, &b_Mu_trigMatchPair_mu18_mu8noL1);
    fChain->SetBranchAddress("Mu_trigMatchPair_mu20_mu8noL1", &Mu_trigMatchPair_mu20_mu8noL1, &b_Mu_trigMatchPair_mu20_mu8noL1);
    fChain->SetBranchAddress("Mu_trigMatchPair_mu22_mu8noL1", &Mu_trigMatchPair_mu22_mu8noL1, &b_Mu_trigMatchPair_mu22_mu8noL1);
    fChain->SetBranchAddress("NEl", &NEl, &b_NEl);
    fChain->SetBranchAddress("El_eta", &El_eta, &b_El_eta);
    fChain->SetBranchAddress("El_etaclus", &El_etaclus, &b_El_etaclus);
    fChain->SetBranchAddress("El_phi", &El_phi, &b_El_phi);
    fChain->SetBranchAddress("El_pT", &El_pT, &b_El_pT);
    fChain->SetBranchAddress("El_E", &El_E, &b_El_E);
    fChain->SetBranchAddress("El_charge", &El_charge, &b_El_charge);
    fChain->SetBranchAddress("El_sigd0", &El_sigd0, &b_El_sigd0);
    fChain->SetBranchAddress("El_z0sinTheta", &El_z0sinTheta, &b_El_z0sinTheta);
    fChain->SetBranchAddress("El_d0pvtx", &El_d0pvtx, &b_El_d0pvtx);
    fChain->SetBranchAddress("El_passOR", &El_passOR, &b_El_passOR);
    fChain->SetBranchAddress("El_SFwMediumLH", &El_SFwMediumLH, &b_El_SFwMediumLH);
    fChain->SetBranchAddress("El_IsoSFwMediumLH", &El_IsoSFwMediumLH, &b_El_IsoSFwMediumLH);
    fChain->SetBranchAddress("El_SFwTightLH", &El_SFwTightLH, &b_El_SFwTightLH);
    fChain->SetBranchAddress("El_SFwLooseAndBLayerLH", &El_SFwLooseAndBLayerLH, &b_El_SFwLooseAndBLayerLH);
    fChain->SetBranchAddress("El_SFwTrigMediumLH_e12_lhloose_L1EM10VH", &El_SFwTrigMediumLH_e12_lhloose_L1EM10VH, &b_El_SFwTrigMediumLH_e12_lhloose_L1EM10VH);
    fChain->SetBranchAddress("El_SFwTrigMediumLH_e17_lhloose", &El_SFwTrigMediumLH_e17_lhloose, &b_El_SFwTrigMediumLH_e17_lhloose);
    fChain->SetBranchAddress("El_SFwTrigMediumLH_single", &El_SFwTrigMediumLH_single, &b_El_SFwTrigMediumLH_single);
    fChain->SetBranchAddress("El_SFwTrigLooseAndBLayerLH_e12_lhloose_L1EM10VH", &El_SFwTrigLooseAndBLayerLH_e12_lhloose_L1EM10VH, &b_El_SFwTrigLooseAndBLayerLH_e12_lhloose_L1EM10VH);
    fChain->SetBranchAddress("El_isLooseAndBLayerLH", &El_isLooseAndBLayerLH, &b_El_isLooseAndBLayerLH);
    fChain->SetBranchAddress("El_isMediumLH", &El_isMediumLH, &b_El_isMediumLH);
    fChain->SetBranchAddress("El_isTightLH", &El_isTightLH, &b_El_isTightLH);
    fChain->SetBranchAddress("El_nBLayerHits", &El_nBLayerHits, &b_El_nBLayerHits);
    fChain->SetBranchAddress("El_expectBLayerHit", &El_expectBLayerHit, &b_El_expectBLayerHit);
    fChain->SetBranchAddress("El_truthType", &El_truthType, &b_El_truthType);
    fChain->SetBranchAddress("El_truthOrigin", &El_truthOrigin, &b_El_truthOrigin);
    fChain->SetBranchAddress("El_bkgTruthType", &El_bkgTruthType, &b_El_bkgTruthType);
    fChain->SetBranchAddress("El_bkgTruthOrigin", &El_bkgTruthOrigin, &b_El_bkgTruthOrigin);
    fChain->SetBranchAddress("El_bkgMotherPdgId", &El_bkgMotherPdgId, &b_El_bkgMotherPdgId);
    fChain->SetBranchAddress("El_firstEgMotherTruthType", &El_firstEgMotherTruthType, &b_El_firstEgMotherTruthType);
    fChain->SetBranchAddress("El_firstEgMotherTruthOrigin", &El_firstEgMotherTruthOrigin, &b_El_firstEgMotherTruthOrigin);
    fChain->SetBranchAddress("El_firstEgMotherPdgId", &El_firstEgMotherPdgId, &b_El_firstEgMotherPdgId);
    fChain->SetBranchAddress("El_chFlip", &El_chFlip, &b_El_chFlip);
    fChain->SetBranchAddress("El_ptcone20", &El_ptcone20, &b_El_ptcone20);
    fChain->SetBranchAddress("El_ptcone30", &El_ptcone30, &b_El_ptcone30);
    fChain->SetBranchAddress("El_ptcone40", &El_ptcone40, &b_El_ptcone40);
    fChain->SetBranchAddress("El_ptvarcone20", &El_ptvarcone20, &b_El_ptvarcone20);
    fChain->SetBranchAddress("El_ptvarcone30", &El_ptvarcone30, &b_El_ptvarcone30);
    fChain->SetBranchAddress("El_ptvarcone40", &El_ptvarcone40, &b_El_ptvarcone40);
    fChain->SetBranchAddress("El_topoetcone20", &El_topoetcone20, &b_El_topoetcone20);
    fChain->SetBranchAddress("El_topoetcone30", &El_topoetcone30, &b_El_topoetcone30);
    fChain->SetBranchAddress("El_topoetcone40", &El_topoetcone40, &b_El_topoetcone40);
    fChain->SetBranchAddress("El_passIsoLooseTO", &El_passIsoLooseTO, &b_El_passIsoLooseTO);
    fChain->SetBranchAddress("El_passIsoLoose", &El_passIsoLoose, &b_El_passIsoLoose);
    fChain->SetBranchAddress("El_passIsoTight", &El_passIsoTight, &b_El_passIsoTight);
    fChain->SetBranchAddress("El_passIsoGrad", &El_passIsoGrad, &b_El_passIsoGrad);
    fChain->SetBranchAddress("El_passIsoGradCustomTight", &El_passIsoGradCustomTight, &b_El_passIsoGradCustomTight);
    fChain->SetBranchAddress("El_passIsoGradCustom", &El_passIsoGradCustom, &b_El_passIsoGradCustom);
    fChain->SetBranchAddress("El_passIsoGradLoose", &El_passIsoGradLoose, &b_El_passIsoGradLoose);
    fChain->SetBranchAddress("El_trigMatch_e12_lhloose_L1EM10VH", &El_trigMatch_e12_lhloose_L1EM10VH, &b_El_trigMatch_e12_lhloose_L1EM10VH);
    fChain->SetBranchAddress("El_trigMatch_e17_lhloose", &El_trigMatch_e17_lhloose, &b_El_trigMatch_e17_lhloose);
    fChain->SetBranchAddress("El_trigMatch_e60_lhmedium", &El_trigMatch_e60_lhmedium, &b_El_trigMatch_e60_lhmedium);
    fChain->SetBranchAddress("El_trigMatch_e24_lhmedium_L1EM20VH", &El_trigMatch_e24_lhmedium_L1EM20VH, &b_El_trigMatch_e24_lhmedium_L1EM20VH);
    fChain->SetBranchAddress("El_trigMatch_e24_lhmedium_iloose_L1EM20VH", &El_trigMatch_e24_lhmedium_iloose_L1EM20VH, &b_El_trigMatch_e24_lhmedium_iloose_L1EM20VH);
    fChain->SetBranchAddress("El_trigMatch_e26_lhtight_nod0_ivarloose", &El_trigMatch_e26_lhtight_nod0_ivarloose, &b_El_trigMatch_e26_lhtight_nod0_ivarloose);
    fChain->SetBranchAddress("El_trigMatch_2e12_lhloose_L12EM10VH", &El_trigMatch_2e12_lhloose_L12EM10VH, &b_El_trigMatch_2e12_lhloose_L12EM10VH);
    fChain->SetBranchAddress("El_trigMatch_2e15_lhloose_L12EM10VH", &El_trigMatch_2e15_lhloose_L12EM10VH, &b_El_trigMatch_2e15_lhloose_L12EM10VH);
    fChain->SetBranchAddress("El_trigMatch_2e15_lhvloose_L12EM13VH", &El_trigMatch_2e15_lhvloose_L12EM13VH, &b_El_trigMatch_2e15_lhvloose_L12EM13VH);
    fChain->SetBranchAddress("El_trigMatch_2e15_lhvloose_nod0_L12EM13VH", &El_trigMatch_2e15_lhvloose_nod0_L12EM13VH, &b_El_trigMatch_2e15_lhvloose_nod0_L12EM13VH);
    fChain->SetBranchAddress("El_trigMatch_2e17_lhvloose_nod0", &El_trigMatch_2e17_lhvloose_nod0, &b_El_trigMatch_2e17_lhvloose_nod0);
    fChain->SetBranchAddress("El_trigMatch_e17_lhloose_mu14", &El_trigMatch_e17_lhloose_mu14, &b_El_trigMatch_e17_lhloose_mu14);
    fChain->SetBranchAddress("El_trigMatch_e17_lhloose_nod0_mu14", &El_trigMatch_e17_lhloose_nod0_mu14, &b_El_trigMatch_e17_lhloose_nod0_mu14);
    fChain->SetBranchAddress("El_trigMatch_e24_lhmedium_nod0_ivarloose", &El_trigMatch_e24_lhmedium_nod0_ivarloose, &b_El_trigMatch_e24_lhmedium_nod0_ivarloose);
    fChain->SetBranchAddress("El_trigMatch_e24_lhtight_nod0_ivarloose", &El_trigMatch_e24_lhtight_nod0_ivarloose, &b_El_trigMatch_e24_lhtight_nod0_ivarloose);
    fChain->SetBranchAddress("El_trigMatch_e60_lhmedium_nod0", &El_trigMatch_e60_lhmedium_nod0, &b_El_trigMatch_e60_lhmedium_nod0);
    fChain->SetBranchAddress("El_passChargeFlipTaggerBDTmedium", &El_passChargeFlipTaggerBDTmedium, &b_El_passChargeFlipTaggerBDTmedium);
    fChain->SetBranchAddress("El_passChargeFlipTaggerBDTloose", &El_passChargeFlipTaggerBDTloose, &b_El_passChargeFlipTaggerBDTloose);
    fChain->SetBranchAddress("NJet", &NJet, &b_NJet);
    fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
    fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
    fChain->SetBranchAddress("Jet_pT", &Jet_pT, &b_Jet_pT);
    fChain->SetBranchAddress("Jet_E", &Jet_E, &b_Jet_E);
    fChain->SetBranchAddress("Jet_quality", &Jet_quality, &b_Jet_quality);
    fChain->SetBranchAddress("Jet_JVT", &Jet_JVT, &b_Jet_JVT);
    fChain->SetBranchAddress("Jet_JVTsf", &Jet_JVTsf, &b_Jet_JVTsf);
    fChain->SetBranchAddress("Jet_MV2c20", &Jet_MV2c20, &b_Jet_MV2c20);
    fChain->SetBranchAddress("Jet_MV2c10", &Jet_MV2c10, &b_Jet_MV2c10);
    fChain->SetBranchAddress("Jet_SFw", &Jet_SFw, &b_Jet_SFw);
    fChain->SetBranchAddress("Jet_ConeTruthLabel", &Jet_ConeTruthLabel, &b_Jet_ConeTruthLabel);
    fChain->SetBranchAddress("Jet_PartonTruthLabel", &Jet_PartonTruthLabel, &b_Jet_PartonTruthLabel);
    fChain->SetBranchAddress("Jet_HadronConeExclTruthLabel", &Jet_HadronConeExclTruthLabel, &b_Jet_HadronConeExclTruthLabel);
    fChain->SetBranchAddress("Jet_deltaR", &Jet_deltaR, &b_Jet_deltaR);
    fChain->SetBranchAddress("Jet_nTrk", &Jet_nTrk, &b_Jet_nTrk);
    fChain->SetBranchAddress("Jet_passOR", &Jet_passOR, &b_Jet_passOR);
    fChain->SetBranchAddress("totalJVTsf", &totalJVTsf, &b_totalJVTsf);
    fChain->SetBranchAddress("Etmiss_CST_Etx", &Etmiss_CST_Etx, &b_Etmiss_CST_Etx);
    fChain->SetBranchAddress("Etmiss_CST_Ety", &Etmiss_CST_Ety, &b_Etmiss_CST_Ety);
    fChain->SetBranchAddress("Etmiss_CST_Et", &Etmiss_CST_Et, &b_Etmiss_CST_Et);
    fChain->SetBranchAddress("Etmiss_TST_Etx", &Etmiss_TST_Etx, &b_Etmiss_TST_Etx);
    fChain->SetBranchAddress("Etmiss_TST_Ety", &Etmiss_TST_Ety, &b_Etmiss_TST_Ety);
    fChain->SetBranchAddress("Etmiss_TST_Et", &Etmiss_TST_Et, &b_Etmiss_TST_Et);
    fChain->SetBranchAddress("Etmiss_TSTterm_Etx", &Etmiss_TSTterm_Etx, &b_Etmiss_TSTterm_Etx);
    fChain->SetBranchAddress("Etmiss_TSTterm_Ety", &Etmiss_TSTterm_Ety, &b_Etmiss_TSTterm_Ety);
    fChain->SetBranchAddress("Etmiss_TSTterm_Et", &Etmiss_TSTterm_Et, &b_Etmiss_TSTterm_Et);
    fChain->SetBranchAddress("Etmiss_Truth_Etx", &Etmiss_Truth_Etx, &b_Etmiss_Truth_Etx);
    fChain->SetBranchAddress("Etmiss_Truth_Ety", &Etmiss_Truth_Ety, &b_Etmiss_Truth_Ety);
    fChain->SetBranchAddress("Etmiss_Truth_Et", &Etmiss_Truth_Et, &b_Etmiss_Truth_Et);
    fChain->SetBranchAddress("NTruthJet", &NTruthJet, &b_NTruthJet);
    fChain->SetBranchAddress("TruthJet_eta", &TruthJet_eta, &b_TruthJet_eta);
    fChain->SetBranchAddress("TruthJet_phi", &TruthJet_phi, &b_TruthJet_phi);
    fChain->SetBranchAddress("TruthJet_pT", &TruthJet_pT, &b_TruthJet_pT);
    fChain->SetBranchAddress("TruthJet_E", &TruthJet_E, &b_TruthJet_E);
    fChain->SetBranchAddress("TruthJet_ConeTruthLabel", &TruthJet_ConeTruthLabel, &b_TruthJet_ConeTruthLabel);
    fChain->SetBranchAddress("TruthJet_PartonTruthLabel", &TruthJet_PartonTruthLabel, &b_TruthJet_PartonTruthLabel);
    fChain->SetBranchAddress("NTruthL", &NTruthL, &b_NTruthL);
    fChain->SetBranchAddress("TruthL_eta", &TruthL_eta, &b_TruthL_eta);
    fChain->SetBranchAddress("TruthL_phi", &TruthL_phi, &b_TruthL_phi);
    fChain->SetBranchAddress("TruthL_pT", &TruthL_pT, &b_TruthL_pT);
    fChain->SetBranchAddress("TruthL_id", &TruthL_id, &b_TruthL_id);
    fChain->SetBranchAddress("SUSY_Spart_pdgId1", &SUSY_Spart_pdgId1, &b_SUSY_Spart_pdgId1);
    fChain->SetBranchAddress("SUSY_Spart_pdgId2", &SUSY_Spart_pdgId2, &b_SUSY_Spart_pdgId2);
    fChain->SetBranchAddress("SUSY_Gluino_decay1", &SUSY_Gluino_decay1, &b_SUSY_Gluino_decay1);
    fChain->SetBranchAddress("SUSY_Gluino_decay2", &SUSY_Gluino_decay2, &b_SUSY_Gluino_decay2);
    fChain->SetBranchAddress("GenFiltHT", &GenFiltHT, &b_GenFiltHT);
    fChain->SetBranchAddress("GenFiltMET", &GenFiltMET, &b_GenFiltMET);
    fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
    fChain->SetBranchAddress("Nvtx", &Nvtx, &b_Nvtx);
    fChain->SetBranchAddress("TruthX1", &TruthX1, &b_TruthX1);
    fChain->SetBranchAddress("TruthX2", &TruthX2, &b_TruthX2);
    fChain->SetBranchAddress("TruthQ", &TruthQ, &b_TruthQ);
    fChain->SetBranchAddress("SherpaNjetWeight", &SherpaNjetWeight, &b_SherpaNjetWeight);
    fChain->SetBranchAddress("TruthPDGID1", &TruthPDGID1, &b_TruthPDGID1);
    fChain->SetBranchAddress("TruthPDGID2", &TruthPDGID2, &b_TruthPDGID2);
    fChain->SetBranchAddress("normalization", &normalization, &b_normalization);
    fChain->SetBranchAddress("luminosity", &luminosity, &b_luminosity);
    fChain->SetBranchAddress("cross_section", &cross_section, &b_cross_section);
    fChain->SetBranchAddress("k_factor", &k_factor, &b_k_factor);
    fChain->SetBranchAddress("filter_efficiency", &filter_efficiency, &b_filter_efficiency);
    fChain->SetBranchAddress("cross_section_kfactor_efficiency", &cross_section_kfactor_efficiency, &b_cross_section_kfactor_efficiency);
    fChain->SetBranchAddress("event_weight_sum", &event_weight_sum, &b_event_weight_sum);
    fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
    fChain->SetBranchAddress("pileup_weight", &pileup_weight, &b_pileup_weight);
    fChain->SetBranchAddress("baseline_weight", &baseline_weight, &b_baseline_weight);
    fChain->SetBranchAddress("signal_weight", &signal_weight, &b_signal_weight);
    fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
    fChain->SetBranchAddress("new_AvgMu", &new_AvgMu, &b_new_AvgMu);
    fChain->SetBranchAddress("El_isBaseline", &El_isBaseline, &b_El_isBaseline);
    fChain->SetBranchAddress("El_isSignal", &El_isSignal, &b_El_isSignal);
    fChain->SetBranchAddress("El_isZTag", &El_isZTag, &b_El_isZTag);
    fChain->SetBranchAddress("El_isZProbe", &El_isZProbe, &b_El_isZProbe);
    fChain->SetBranchAddress("El_Tag_index", &El_Tag_index, &b_El_Tag_index);
    fChain->SetBranchAddress("El_Probe_index", &El_Probe_index, &b_El_Probe_index);
    fChain->SetBranchAddress("El_ZTandP_mll", &El_ZTandP_mll, &b_El_ZTandP_mll);
    fChain->SetBranchAddress("El_isZProbe_TagTriggerMatched", &El_isZProbe_TagTriggerMatched, &b_El_isZProbe_TagTriggerMatched);
    fChain->SetBranchAddress("El_Tag_Trigger_Matched", &El_Tag_Trigger_Matched, &b_El_Tag_Trigger_Matched);
    fChain->SetBranchAddress("El_Probe_Trigger_Matched", &El_Probe_Trigger_Matched, &b_El_Probe_Trigger_Matched);
    fChain->SetBranchAddress("El_DR_closest_Jet", &El_DR_closest_Jet, &b_El_DR_closest_Jet);
    fChain->SetBranchAddress("Mu_isBaseline", &Mu_isBaseline, &b_Mu_isBaseline);
    fChain->SetBranchAddress("Mu_isSignal", &Mu_isSignal, &b_Mu_isSignal);
    fChain->SetBranchAddress("Mu_isZTag", &Mu_isZTag, &b_Mu_isZTag);
    fChain->SetBranchAddress("Mu_isZProbe", &Mu_isZProbe, &b_Mu_isZProbe);
    fChain->SetBranchAddress("Mu_Tag_index", &Mu_Tag_index, &b_Mu_Tag_index);
    fChain->SetBranchAddress("Mu_Probe_index", &Mu_Probe_index, &b_Mu_Probe_index);
    fChain->SetBranchAddress("Mu_ZTandP_mll", &Mu_ZTandP_mll, &b_Mu_ZTandP_mll);
    fChain->SetBranchAddress("Mu_isZProbe_TagTriggerMatched", &Mu_isZProbe_TagTriggerMatched, &b_Mu_isZProbe_TagTriggerMatched);
    fChain->SetBranchAddress("Mu_Tag_Trigger_Matched", &Mu_Tag_Trigger_Matched, &b_Mu_Tag_Trigger_Matched);
    fChain->SetBranchAddress("Mu_Probe_Trigger_Matched", &Mu_Probe_Trigger_Matched, &b_Mu_Probe_Trigger_Matched);
    fChain->SetBranchAddress("Mu_DR_closest_Jet", &Mu_DR_closest_Jet, &b_Mu_DR_closest_Jet);
    fChain->SetBranchAddress("baseline_mll", &baseline_mll, &b_baseline_mll);
    fChain->SetBranchAddress("signal_mll", &signal_mll, &b_signal_mll);

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode ytRealLeptonsEfficiency::execute ()
{
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    //const char * function_name = "execute()";
    //Info(function_name, "Function calls");

    wk()->tree()->GetEntry(wk()->treeEntry());

    // Dump histograms
    double weight = 1.;
    if (isData) {
        weight = 1.;
    }
    else if (isMC) {
        weight = normalization;
    }

    h_Nvtx->Fill(Nvtx, weight);
    if (isData) {  // for Data: weight = 1.;
        h_AvgMu->Fill(new_AvgMu * 1.0 / 1.09, weight);
    }
    else if (isMC) {
        h_AvgMu->Fill(new_AvgMu, weight);
    }
    h_zPV->Fill(PV_z, weight);

    TLorentzVector tlv_lept0, tlv_lept1;
    if (NEl == 2 &&
        El_charge->at(0) * El_charge->at(1) == -1) {
        // OSee
        if (isData) {  // for Data: weight = 1.;
            h_AvgMu_OSee->Fill(new_AvgMu * 1.0 / 1.09, weight);
        }
        else if (isMC) {
            h_AvgMu_OSee->Fill(new_AvgMu, weight);
        }

        tlv_lept0.SetPtEtaPhiE(El_pT->at(0), El_eta->at(0), El_phi->at(0), El_E->at(0));
        tlv_lept1.SetPtEtaPhiE(El_pT->at(1), El_eta->at(1), El_phi->at(1), El_E->at(1));
        double mll = (tlv_lept0 + tlv_lept1).M();
        if (El_isBaseline->at(0) == 1 && El_isBaseline->at(1) == 1)
            h_baseline_OSee_mll->Fill(mll / 1000., weight);
        if (El_isSignal->at(0) == 1 && El_isSignal->at(1) == 1)
            h_signal_OSee_mll->Fill(mll / 1000., weight);
    }
    if (NMu == 2 &&
        Mu_charge->at(0) * Mu_charge->at(1) == -1) {
        // OSmumu
        if (isData) { // for Data: weight = 1.;
            h_AvgMu_OSmumu->Fill(new_AvgMu * 1.0 / 1.09, weight);
        }
        else if (isMC) {
            h_AvgMu_OSmumu->Fill(new_AvgMu, weight);
        }

        float Mu_Mass = 105.6583715;
        tlv_lept0.SetPtEtaPhiM(Mu_pT->at(0), Mu_eta->at(0), Mu_phi->at(0), Mu_Mass);
        tlv_lept1.SetPtEtaPhiM(Mu_pT->at(1), Mu_eta->at(1), Mu_phi->at(1), Mu_Mass);
        double mll = (tlv_lept0 + tlv_lept1).M();
        if (Mu_isBaseline->at(0) == 1 && Mu_isBaseline->at(1) == 1)
            h_baseline_OSmumu_mll->Fill(mll / 1000., weight);
        if (Mu_isSignal->at(0) == 1 && Mu_isSignal->at(1) == 1)
            h_signal_OSmumu_mll->Fill(mll / 1000., weight);
    }
/*
    if (vec_lept.size() > 0)
        h_NLepts->Fill(vec_lept.at(0).get_number(), weight);
    else
        h_NLepts->Fill(0., weight);

    if (vec_jets.size() > 0)
        h_NJets->Fill(vec_jets.at(0).get_number(), weight);
    else
        h_NJets->Fill(0., weight);
*/

/*
    if (isData) {
        h_Nvtx->Fill(Nvtx);
        h_AvgMu->Fill(new_AvgMu * 1.0 / 1.09);
        h_zPV->Fill(PV_z);

        TLorentzVector tlv_lept0, tlv_lept1;
        if (NEl == 2 &&
            El_charge->at(0) * El_charge->at(1) == -1) {
            // OSee
            h_AvgMu_OSee->Fill(new_AvgMu * 1.0 / 1.09);

            tlv_lept0.SetPtEtaPhiE(El_pT->at(0), El_eta->at(0), El_phi->at(0), El_E->at(0));
            tlv_lept1.SetPtEtaPhiE(El_pT->at(1), El_eta->at(1), El_phi->at(1), El_E->at(1));
            double mll = (tlv_lept0 + tlv_lept1).M();
            if (El_isBaseline->at(0) == 1 && El_isBaseline->at(1) == 1)
                h_baseline_OSee_mll->Fill(mll / 1000.);
            if (El_isSignal->at(0) == 1 && El_isSignal->at(1) == 1)
                h_signal_OSee_mll->Fill(mll / 1000.);
        }
        if (NMu == 2 &&
            Mu_charge->at(0) * Mu_charge->at(1) == -1) {
            // OSmumu
            h_AvgMu_OSmumu->Fill(new_AvgMu * 1.0 / 1.09);

            float Mu_Mass = 105.6583715;
            tlv_lept0.SetPtEtaPhiM(Mu_pT->at(0), Mu_eta->at(0), Mu_phi->at(0), Mu_Mass);
            tlv_lept1.SetPtEtaPhiM(Mu_pT->at(1), Mu_eta->at(1), Mu_phi->at(1), Mu_Mass);
            double mll = (tlv_lept0 + tlv_lept1).M();
            if (Mu_isBaseline->at(0) == 1 && Mu_isBaseline->at(1) == 1)
                h_baseline_OSmumu_mll->Fill(mll / 1000.);
            if (Mu_isSignal->at(0) == 1 && Mu_isSignal->at(1) == 1)
                h_signal_OSmumu_mll->Fill(mll / 1000.);
        }
//
        if (vec_lept.size() > 0)
            h_NLepts->Fill(vec_lept.at(0).get_number());
        else
            h_NLepts->Fill(0.);

        if (vec_jets.size() > 0)
            h_NJets->Fill(vec_jets.at(0).get_number());
        else
            h_NJets->Fill(0.);
//
    }
    else if (isMC) {
        //double weight = luminosity * cross_section_kfactor_efficiency * 1000. * EventWeight * pileup_weight / derivation_stat_weights;
        double weight = normalization;
        h_Nvtx_weighted->Fill(Nvtx, weight);
        h_AvgMu_weighted->Fill(new_AvgMu, weight);
        h_zPV_weighted->Fill(PV_z, weight);

        TLorentzVector tlv_lept0, tlv_lept1;
        if (NEl == 2 &&
            El_charge->at(0) * El_charge->at(1) == -1) {
            // OSee
            h_AvgMu_OSee_weighted->Fill(new_AvgMu, weight);

            tlv_lept0.SetPtEtaPhiE(El_pT->at(0), El_eta->at(0), El_phi->at(0), El_E->at(0));
            tlv_lept1.SetPtEtaPhiE(El_pT->at(1), El_eta->at(1), El_phi->at(1), El_E->at(1));
            double mll = (tlv_lept0 + tlv_lept1).M();
            if (El_isBaseline->at(0) == 1 && El_isBaseline->at(1) == 1)
                h_baseline_OSee_mll_weighted->Fill(mll / 1000., weight);
            if (El_isSignal->at(0) == 1 && El_isSignal->at(1) == 1)
                h_signal_OSee_mll_weighted->Fill(mll / 1000., weight);
        }
        if (NMu == 2 &&
            Mu_charge->at(0) * Mu_charge->at(1) == -1) {
            // OSmumu
            h_AvgMu_OSmumu_weighted->Fill(new_AvgMu, weight);

            float Mu_Mass = 105.6583715;
            tlv_lept0.SetPtEtaPhiM(Mu_pT->at(0), Mu_eta->at(0), Mu_phi->at(0), Mu_Mass);
            tlv_lept1.SetPtEtaPhiM(Mu_pT->at(1), Mu_eta->at(1), Mu_phi->at(1), Mu_Mass);
            double mll = (tlv_lept0 + tlv_lept1).M();
            if (Mu_isBaseline->at(0) == 1 && Mu_isBaseline->at(1) == 1)
                h_baseline_OSmumu_mll_weighted->Fill(mll / 1000., weight);
            if (Mu_isSignal->at(0) == 1 && Mu_isSignal->at(1) == 1)
                h_signal_OSmumu_mll_weighted->Fill(mll / 1000., weight);
        }
//
        if (vec_lept.size() > 0)
            h_NLepts_weighted->Fill(vec_lept.at(0).get_number(), weight);
        else
            h_NLepts_weighted->Fill(0., weight);

        if (vec_jets.size() > 0)
            h_NJets_weighted->Fill(vec_jets.at(0).get_number(), weight);
        else
            h_NJets_weighted->Fill(0., weight);
//
    }
*/

    // Loop over all leptons
    if (lepton == "electron") {
        // tag leptons must satisfy single lepton trigger requirement.
        // single lepton triggers:
        if (trigger == "single_lepton_trigger") {
            // Electron: e24_lhmedium_L1EM20VH (2015)
            //           e26_lhtight_nod0_ivarloose (2016)
            if (run_number < 290000) { // 2015 data
                if (!HLT_e24_lhmedium_iloose_L1EM20VH) // use HLT_e24_lhmedium_iloose_L1EM20VH instead of using HLT_e24_lhmedium_L1EM20VH
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
            else if (run_number > 290000) { // 2016 data
                if (!HLT_e26_lhtight_nod0_ivarloose)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
        }
        else if (trigger == "dilepton_trigger") {
            // Electron: HLT_2e12_lhloose_L12EM10VH (2015)
            //           HLT_2e17_lhvloose_nod0 (2016)
            if (run_number < 290000) { // 2015 data)
                if (!HLT_2e12_lhloose_L12EM10VH)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
            else if (run_number > 290000) { // 2016 data
                if (!HLT_2e17_lhvloose_nod0)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
        }
        else if (trigger == "dilepton_trigger_but_fail_single_lepton_trigger") {
            if (run_number < 290000) { // 2015 data)
                if ( !(HLT_2e12_lhloose_L12EM10VH && !HLT_e24_lhmedium_iloose_L1EM20VH) ) // use HLT_e24_lhmedium_iloose_L1EM20VH instead of using HLT_e24_lhmedium_L1EM20VH
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
            else if (run_number > 290000) { // 2016 data
                if ( !(HLT_2e17_lhvloose_nod0 && !HLT_e26_lhtight_nod0_ivarloose) )
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
        }
        else if (trigger == "other_trigger") {
            // do something here.
            cout << "Warning: Neet to implement this." << endl;
        }

        this->loop_over_electrons();
    }
    else if (lepton == "muon") {
        // tag leptons must satisfy single lepton trigger requirement.
        // single lepton triggers:
        if (trigger == "single_lepton_trigger") {
            // Muon: mu20_iloose_L1MU15 (2015)
            //       mu26_ivarmedium (2016)
            if (run_number < 290000) { // 2015 data
                if (!HLT_mu20_iloose_L1MU15)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
            else if (run_number > 290000) { // 2016 data
                if (!HLT_mu26_ivarmedium)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
        }
        else if (trigger == "dilepton_trigger") {
            // Muon: HLT_mu18_mu8noL1 (2015)
            //       HLT_mu22_mu8noL1 (2016)
            if (run_number < 290000) { // 2015 data
                if (!HLT_mu18_mu8noL1)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
            else if (run_number > 290000) { // 2016 data
                if (!HLT_mu22_mu8noL1)
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
        }
        else if (trigger == "dilepton_trigger_but_fail_single_lepton_trigger") {
            if (run_number < 290000) { // 2015 data
                if ( !(HLT_mu18_mu8noL1 && !HLT_mu20_iloose_L1MU15) )
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
            else if (run_number > 290000) { // 2016 data
                if ( !(HLT_mu22_mu8noL1 && !HLT_mu26_ivarmedium) )
                    return EL::StatusCode::SUCCESS; // Go to next event
            }
        }
        else if(trigger == "other_trigger") {
            // do something here.
            cout << "Warning: Neet to implement this." << endl;
        }

        this->loop_over_muons();
    }

    // Calculate the cut efficiency
    // option="B": Binomial errors are computed
    h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Divide(h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt, h_baseline_pt, 1, 1, "B");
    h_cut_eff_CaloIso->Divide(h_baseline_CaloIso_cut_pt, h_baseline_pt, 1, 1, "B");
    h_cut_eff_TrackIso->Divide(h_baseline_TrackIso_cut_pt, h_baseline_pt, 1, 1, "B");
    h_cut_eff_Iso->Divide(h_baseline_Iso_cut_pt, h_baseline_pt, 1, 1, "B");
    h_cut_eff_z0->Divide(h_baseline_z0_cut_pt, h_baseline_pt, 1, 1, "B");
    h_cut_eff_sigd0->Divide(h_baseline_sigd0_cut_pt, h_baseline_pt, 1, 1, "B");

    // Calculate the real lepton efficiency
    h_eff_pt->Divide(h_signal_pt, h_baseline_pt, 1, 1, "B");
    h_eff_eta->Divide(h_signal_eta, h_baseline_eta, 1, 1, "B");
    h_eff_nJets->Divide(h_signal_nJets, h_baseline_nJets, 1, 1, "B");
    h_eff_dRjet->Divide(h_signal_dRjet, h_baseline_dRjet, 1, 1, "B");
    //h_eff_Etmiss->Divide(h_signal_Etmiss, h_baseline_Etmiss, 1, 1, "B");
    //h_eff_meff->Divide(h_signal_meff, h_baseline_meff, 1, 1, "B");

    //h_2d_eff_pt_eta->Divide(h_signal_pt_eta, h_baseline_pt_eta,1, 1, "B");
    //h_2d_eff_pt_nJets->Divide(h_signal_pt_nJets, h_baseline_pt_nJets,1, 1, "B");
    //h_2d_eff_pt_dRjet->Divide(h_signal_pt_dRjet, h_baseline_pt_dRjet,1, 1, "B");
    //h_2d_eff_pt_Etmiss->Divide(h_signal_pt_Etmiss, h_baseline_pt_Etmiss,1, 1, "B");
    //_2d_eff_pt_meff->Divide(h_signal_pt_meff, h_baseline_pt_meff,1, 1, "B");

    return EL::StatusCode::SUCCESS;
}




EL::StatusCode ytRealLeptonsEfficiency::postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.

    //const char * function_name = "postExecute()";
    //Info(function_name, "Function calls");

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode ytRealLeptonsEfficiency::finalize ()
{
    // This method is the mirror image of initialize(), meaning it gets
    // called after the last event has been processed on the worker node
    // and allows you to finish up any objects you created in
    // initialize() before they are written to disk.  This is actually
    // fairly rare, since this happens separately for each worker node.
    // Most of the time you want to do your post-processing on the
    // submission node after all your histogram outputs have been
    // merged.  This is different from histFinalize() in that it only
    // gets called on worker nodes that processed input events.

    //const char * function_name = "finalize()";
    //Info(function_name, "Function calls");

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode ytRealLeptonsEfficiency::histFinalize ()
{
    // This method is the mirror image of histInitialize(), meaning it
    // gets called after the last event has been processed on the worker
    // node and allows you to finish up any objects you created in
    // histInitialize() before they are written to disk.  This is
    // actually fairly rare, since this happens separately for each
    // worker node.  Most of the time you want to do your
    // post-processing on the submission node after all your histogram
    // outputs have been merged.  This is different from finalize() in
    // that it gets called on all worker nodes regardless of whether
    // they processed input events.

    //const char * function_name = "histFinalize()";
    //Info(function_name, "Function calls");

    this->make_baseline_mll_plot();
    this->make_cut_efficiency_plot();

    if (lepton == "electron") {
        this->make_electron_real_efficiency_plot();
    }
    if (lepton == "muon") {
        this->make_muon_real_efficiency_plot();
    }

    return EL::StatusCode::SUCCESS;
}



void ytRealLeptonsEfficiency::set_binning_default()
{
    m_mll_bins      = vector<float> (m_default_mll_bins, m_default_mll_bins + sizeof(m_default_mll_bins) / sizeof(float));
    m_pt_bins       = vector<float> (m_default_pt_bins, m_default_pt_bins + sizeof(m_default_pt_bins) / sizeof(float));
    m_eta_bins      = vector<float> (m_default_eta_bins, m_default_eta_bins + sizeof(m_default_eta_bins) / sizeof(float));
    m_deltaR_bins   = vector<float> (m_default_deltaR_bins, m_default_deltaR_bins + sizeof(m_default_deltaR_bins) / sizeof(float));
    m_NJet_bins     = vector<float> (m_default_NJet_bins, m_default_NJet_bins + sizeof(m_default_NJet_bins) / sizeof(float));
    //m_Etmiss_bins   = vector<float> (m_default_Etmiss_bins, m_default_Etmiss_bins + sizeof(m_default_Etmiss_bins) / sizeof(float));
    //m_meff_bins     = vector<float> (m_default_meff_bins, m_default_meff_bins + sizeof(m_default_meff_bins) / sizeof(float));
    m_AvgMu_bins    = vector<float> (m_default_AvgMu_bins, m_default_AvgMu_bins + sizeof(m_default_AvgMu_bins) / sizeof(float));
}



void ytRealLeptonsEfficiency::set_binning(int n, const float *bins, vector<float> &vec)
{
    vec.clear();
    for (int i = 0; i < n; i++) {
        vec.push_back(bins[i]);
    }
}



void ytRealLeptonsEfficiency::loop_over_electrons()
{
    for (int n_el = 0; n_el < NEl; n_el++) {
        if (!El_isBaseline->at(n_el)) // reject non baseline electrons
            continue;
        // All the electrons here are baseline electrons

        // Truth match (only apply on the non Z MC, such as Gtt and ttbar).
        // We want to keep isolated Electron (2) and these isolated electrons must come from top (10), W bosons (12), and Z bosons (13).
        bool truth_match = false;
        if (isMC) {
            if (process == "GG_ttn1") {
                if (El_truthType->at(n_el) == 2 ||
                    (El_truthType->at(n_el) == 4 && (El_bkgTruthOrigin->at(n_el) == 10 || El_bkgTruthOrigin->at(n_el) == 12 || El_bkgTruthOrigin->at(n_el) == 13 || El_bkgTruthOrigin->at(n_el) == 40)))
                    // the bkgOriginal is related charge filp. Ximo said just add the requirements at here.
                    truth_match = true;
            }
            else { // process != GG_ttn1
                if ((El_truthType->at(n_el) == 2 && (El_truthOrigin->at(n_el) == 10 || El_truthOrigin->at(n_el) == 12 || El_truthOrigin->at(n_el) == 13)) ||
                    (El_truthType->at(n_el) == 4 && (El_bkgTruthOrigin->at(n_el) == 10 || El_bkgTruthOrigin->at(n_el) == 12 || El_bkgTruthOrigin->at(n_el) == 13 || El_bkgTruthOrigin->at(n_el) == 40)))
                    // the bkgOriginal is related charge filp. Ximo said just add the requirements at here.
                    truth_match = true;
            }
        }

        // // Tag trigger matching
        // if (tag_trigger_matched) {
        //     if (!El_isZProbe_TriggerMatched->at(n_el))
        //         continue;
        // }

        // Background template
        double calo_isolation = El_topoetcone20->at(n_el) / El_pT->at(n_el);
        double track_isolation = El_ptvarcone20->at(n_el) / El_pT->at(n_el);

        float weight = 1.; // Data weight = 1.
        if (isMC) {
            weight = normalization * baseline_weight;
        }

        if (isMC && !truth_match) {
            h_bkg_template_fail_truth->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_bkg_template_fail_truth_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        if (!El_isMediumLH->at(n_el)) { // Background template requires the electron object fails the MediumLLH
            h_bkg_template_fail_id_only->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_bkg_template_fail_id_only_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
            // reverse calo isolation
            if (calo_isolation > 0.10) {
                h_bkg_template_fail_id_and_CaloIso->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
                h_bkg_template_fail_id_and_CaloIso_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
            }
            // reverse track isolation
            if (track_isolation > 0.06) {
                h_bkg_template_fail_id_and_TrackIso->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
                h_bkg_template_fail_id_and_TrackIso_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
            }
            // reverse calo and track isolation (baseline bkg template)
            if (calo_isolation > 0.15 &&
                track_isolation > 0.08) {
                h_bkg_template_fail_id_and_CaloIso_and_TrackIso->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
                h_bkg_template_fail_id_and_CaloIso_and_TrackIso_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
            }
            // reverse calo and track isolation with tighter requirements (variation2 bkg template)
            if (calo_isolation > 0.20 &&
                track_isolation > 0.15) {
                h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
                h_bkg_template_fail_id_and_CaloIso_and_TrackIso_tight_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
            }
        }
        // reverse calo isolation
        if (calo_isolation > 0.20) {
            h_bkg_template_fail_CaloIso_only->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_bkg_template_fail_CaloIso_only_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        // reverse track isolation
        if (track_isolation > 0.15) {
            h_bkg_template_fail_TrackIso_only->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_bkg_template_fail_TrackIso_only_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        // reverse calo and track isolation only (variation1 bkg template)
        if (calo_isolation > 0.06 &&
            track_isolation > 0.06) {
            h_bkg_template_fail_CaloIso_and_TrackIso->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_bkg_template_fail_CaloIso_and_TrackIso_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }

        // 3-dim Cut efficiency
        // Apply MediumLLH
        // if (El_isMediumLH->at(n_el)) {
        //     h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        // }
        // Apply topoetcone20/pt < 0.06
        // if ((El_topoetcone20->at(n_el) /  El_pT->at(n_el)) < 0.06) {
        //     h_baseline_CaloIso_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        // }
        // Apply ptvarcone20/pt < 0.06
        // if ((El_ptvarcone20->at(n_el) / El_pT->at(n_el)) < 0.06) {
        //     h_baseline_TrackIso_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        // }
        // Apply topoetcone20/pt < 0.06 && ptvarcone20/pt < 0.06 at the same time
        // if ((El_topoetcone20->at(n_el) /  El_pT->at(n_el)) < 0.06 &&
        //     (El_ptvarcone20->at(n_el) / El_pT->at(n_el)) < 0.06) {
        //     h_baseline_Iso_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        // }
        // Apply |z0*sin(theta)| < 0.5
        // if (fabs(El_z0sinTheta->at(n_el)) < 0.5) {
        //     h_baseline_z0_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        // }

        // mll window: 60 < mll < 150
        h_baseline_mll->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
        h_baseline_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        // Reset weight for signal leptons
        if (isMC) {
            weight = normalization * signal_weight;
        }
        if (El_isSignal->at(n_el)) {
            h_signal_mll->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_signal_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }

        if (isData) {
            // Use probe electrons only
            if (!El_isZProbe->at(n_el))
                continue;
        }
        else if (isMC) {
            if (process == "Zee") {
                // Use probe electron only
                if (!El_isZProbe->at(n_el))
                    continue;
            }
            else if (process == "Zee_TandP_truth_matched") {
                // Use truth match on top of probe electron
                if (!El_isZProbe->at(n_el))
                    continue;
                if (!truth_match)
                    continue;
            }
            else { // process == "GG_ttn1" || process == "ttbar" || process == "Zee_truth_match" || process == "Zmumu_truth_match"
                // Use truth match
                //cout << "process is " << process << endl;
                if (!truth_match)
                    continue;
            }
        }

        // Tag trigger matching
        if (tag_trigger_matched) {
            if (!El_isZProbe_TagTriggerMatched->at(n_el))
                continue;
        }
/*
        if (isData) {
            h_baseline_mll->Fill(El_ZTandP_mll->at(n_el) / 1000.);
            h_baseline_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000.);

            if (El_isSignal->at(n_el)) {
                h_signal_mll->Fill(El_ZTandP_mll->at(n_el) / 1000.);
                h_signal_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000.);
            }

            // Use probe electrons only
            if (!El_isZProbe->at(n_el))
                continue;
        }
        else if (isMC) {
            h_baseline_mll->Fill(El_ZTandP_mll->at(n_el) / 1000., weight);
            h_baseline_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);

            if (El_isSignal->at(n_el)) {
                h_signal_mll->Fill(El_ZTandP_mll->at(n_el) / 1000., normalization * signal_weight);
                h_signal_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., normalization * signal_weight);
            }

            if (process == "Zee") {
                // Use probe electron only
                if (!El_isZProbe->at(n_el))
                    continue;
            }
            else if (process == "Zee_TandP_truth_matched") {
                // Use truth match on top of probe electron
                if (!El_isZProbe->at(n_el))
                    continue;
                if (!truth_match)
                    continue;
            }
            else { // process == "GG_ttn1" || process == "ttbar" || process == "Zee_truth_match" || process == "Zmumu_truth_match"
                // Use truth match
                //cout << "process is " << process << endl;
                if (!truth_match)
                    continue;
            }
        }
*/
        // Reset weight for signal leptons
        if (isMC) {
            weight = normalization * baseline_weight;
        }

        // We use |eta| < 2 for electrons
        if (fabs(El_eta->at(n_el)) > 2.)
            continue;

        // 1-dim histograms
        //h_baseline_mll->Fill(El_ZTandP_mll->at(n_el) / 1000.); // only dump 80 < mll < 100
        h_baseline_pt->Fill(El_pT->at(n_el) / 1000., weight);
        h_baseline_eta->Fill(fabs(El_eta->at(n_el)), weight);
        //h_baseline_d0pvtx->Fill(El_d0pvtx->at(n_el), weight);
        //h_baseline_sigd0->Fill(El_sigd0->at(n_el), weight);
        //h_baseline_d0err->Fill(El_d0pvtx->at(n_el) / El_sigd0->at(n_el), weight);
        //h_baseline_z0sinTheta->Fill(El_z0sinTheta->at(n_el), weight);
        //h_baseline_ptvarcone20->Fill(El_ptvarcone20->at(n_el) / 1000., weight);
        //h_baseline_ptvarcone30->Fill(El_ptvarcone30->at(n_el) / 1000., weight);
        //h_baseline_topoetcone20->Fill(El_topoetcone20->at(n_el) / 1000., weight);
        h_baseline_nJets->Fill(NJet, weight);
        h_baseline_dRjet->Fill(El_DR_closest_Jet->at(n_el), weight);
        //h_baseline_Etmiss->Fill(Etmiss_TST_Et / 1000., weight);
        //h_baseline_meff->Fill(meff / 1000., weight);
        // 2-dim histograms
        h_baseline_pt_eta->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), weight);
        h_baseline_pt_nJets->Fill(El_pT->at(n_el) / 1000., NJet, weight);
        h_baseline_pt_dRjet->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), weight);
        //h_baseline_pt_Etmiss->Fill(El_pT->at(n_el) / 1000., Etmiss_TST_Et / 1000., weight);
        //h_baseline_pt_meff->Fill(El_pT->at(n_el) / 1000., meff / 1000., weight);
        // 3-dim histograms
        //h_baseline_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);  // only dump 80 < mll < 100
        h_baseline_pt_eta_dRjet->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), weight);
        h_baseline_pt_dRjet_mll->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), El_ZTandP_mll->at(n_el) / 1000., weight);
        h_baseline_eta_dRjet_mll->Fill(fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), El_ZTandP_mll->at(n_el) / 1000., weight);
        if (isData) {
            h_baseline_pt_eta_AvgMu->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), new_AvgMu * 1.0 / 1.09, weight);
            h_baseline_pt_dRjet_AvgMu->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), new_AvgMu * 1.0 / 1.09, weight);
            h_baseline_eta_dRjet_AvgMu->Fill(fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), new_AvgMu * 1.0 / 1.09, weight);
        }
        else if (isMC) {
            h_baseline_pt_eta_AvgMu->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), new_AvgMu, weight);
            h_baseline_pt_dRjet_AvgMu->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), new_AvgMu, weight);
            h_baseline_eta_dRjet_AvgMu->Fill(fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), new_AvgMu, weight);
        }

        // Cut efficiency
        // We use |eta| < 2 for electrons
        // Apply MediumLLH
        if (El_isMediumLH->at(n_el) &&
            fabs(El_eta->at(n_el)) < 2.) {
            h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt->Fill(El_pT->at(n_el) / 1000., weight);
        }
        // Apply topoetcone20/pt < 0.06
        if ((El_topoetcone20->at(n_el) /  El_pT->at(n_el)) < 0.06 &&
            fabs(El_eta->at(n_el)) < 2.) {
            h_baseline_CaloIso_cut_pt->Fill(El_pT->at(n_el) / 1000., weight);
        }
        // Apply ptvarcone20/pt < 0.06
        if ((El_ptvarcone20->at(n_el) / El_pT->at(n_el)) < 0.06 &&
            fabs(El_eta->at(n_el)) < 2.) {
            h_baseline_TrackIso_cut_pt->Fill(El_pT->at(n_el) / 1000., weight);
        }
        // Apply topoetcone20/pt < 0.06 && ptvarcone20/pt < 0.06 at the same time
        if ((El_topoetcone20->at(n_el) /  El_pT->at(n_el)) < 0.06 &&
            (El_ptvarcone20->at(n_el) / El_pT->at(n_el)) < 0.06 &&
            fabs(El_eta->at(n_el)) < 2.) {
            h_baseline_Iso_cut_pt->Fill(El_pT->at(n_el) / 1000., weight);
        }
        // Apply |z0*sin(theta)| < 0.5
        if (fabs(El_z0sinTheta->at(n_el)) < 0.5 &&
            fabs(El_eta->at(n_el)) < 2.) {
            h_baseline_z0_cut_pt->Fill(El_pT->at(n_el) / 1000., weight);
        }

        // 3-dim Cut efficiency
        // Apply MediumLLH
        if (El_isMediumLH->at(n_el)) {
            h_baseline_LooseAndBLayerLLH_to_MediumLLH_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        // Apply topoetcone20/pt < 0.06
        if ((El_topoetcone20->at(n_el) /  El_pT->at(n_el)) < 0.06) {
            h_baseline_CaloIso_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        // Apply ptvarcone20/pt < 0.06
        if ((El_ptvarcone20->at(n_el) / El_pT->at(n_el)) < 0.06) {
            h_baseline_TrackIso_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        // Apply topoetcone20/pt < 0.06 && ptvarcone20/pt < 0.06 at the same time
        if ((El_topoetcone20->at(n_el) /  El_pT->at(n_el)) < 0.06 &&
            (El_ptvarcone20->at(n_el) / El_pT->at(n_el)) < 0.06) {
            h_baseline_Iso_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }
        // Apply |z0*sin(theta)| < 0.5
        if (fabs(El_z0sinTheta->at(n_el)) < 0.5) {
            h_baseline_z0_cut_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight);
        }

        // Reset weight for signal leptons
        if (isMC) {
            weight = normalization * signal_weight;
        }
        // All the electrons here are signal electrons
        if (El_isSignal->at(n_el)) {
            // 1-dim histograms
            //h_signal_mll->Fill(El_ZTandP_mll->at(n_el) / 1000., weight); // only dump 80 < mll < 100
            h_signal_pt->Fill(El_pT->at(n_el) / 1000., weight);
            h_signal_eta->Fill(fabs(El_eta->at(n_el)), weight);
            //h_signal_d0pvtx->Fill(El_d0pvtx->at(n_el), weight);
            //h_signal_sigd0->Fill(El_sigd0->at(n_el), weight);
            //h_signal_d0err->Fill(El_d0pvtx->at(n_el) / El_sigd0->at(n_el), weight);
            //h_signal_z0sinTheta->Fill(El_z0sinTheta->at(n_el), weight);
            //h_signal_ptvarcone20->Fill(El_ptvarcone20->at(n_el) / 1000., weight);
            //h_signal_ptvarcone30->Fill(El_ptvarcone30->at(n_el) / 1000., weight);
            //h_signal_topoetcone20->Fill(El_topoetcone20->at(n_el) / 1000., weight);
            h_signal_nJets->Fill(NJet, weight);
            h_signal_dRjet->Fill(El_DR_closest_Jet->at(n_el), weight);
            //h_signal_Etmiss->Fill(Etmiss_TST_Et / 1000., weight);
            //h_signal_meff->Fill(meff / 1000., weight);
            // 2-dim histograms
            h_signal_pt_eta->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), weight);
            h_signal_pt_nJets->Fill(El_pT->at(n_el) / 1000., NJet, weight);
            h_signal_pt_dRjet->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), weight);
            //h_signal_pt_Etmiss->Fill(El_pT->at(n_el) / 1000.,Etmiss_TST_Et / 1000., weight);
            //h_signal_pt_meff->Fill(El_pT->at(n_el) / 1000., meff / 1000., weight);
            // 3-dim histograms
            //h_signal_pt_eta_mll->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_ZTandP_mll->at(n_el) / 1000., weight); // only dump 80 < mll < 100
            h_signal_pt_eta_dRjet->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), weight);
            h_signal_pt_dRjet_mll->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), El_ZTandP_mll->at(n_el) / 1000., weight);
            h_signal_eta_dRjet_mll->Fill(fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), El_ZTandP_mll->at(n_el) / 1000., weight);
            if (isData) {
                h_signal_pt_eta_AvgMu->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), new_AvgMu * 1.0 / 1.09, weight);
                h_signal_pt_dRjet_AvgMu->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), new_AvgMu * 1.0 / 1.09, weight);
                h_signal_eta_dRjet_AvgMu->Fill(fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), new_AvgMu * 1.0 / 1.09, weight);
            }
            else if (isMC) {
                h_signal_pt_eta_AvgMu->Fill(El_pT->at(n_el) / 1000., fabs(El_eta->at(n_el)), new_AvgMu, weight);
                h_signal_pt_dRjet_AvgMu->Fill(El_pT->at(n_el) / 1000., El_DR_closest_Jet->at(n_el), new_AvgMu, weight);
                h_signal_eta_dRjet_AvgMu->Fill(fabs(El_eta->at(n_el)), El_DR_closest_Jet->at(n_el), new_AvgMu, weight);
            }
        }
    }
}



void ytRealLeptonsEfficiency::loop_over_muons()
{
    for (int n_mu = 0; n_mu < NMu; n_mu++) {
        if (!Mu_isBaseline->at(n_mu)) // reject non baseline muons
            continue;
        // All the muons here are baseline muons

        // Truth match (only apply on the non Z MC, such as Gtt and ttbar).
        // We want to keep isolated Muons (6) and these isolated muons must come from top (10), W bosons (12), and Z bosons (13).
        bool truth_match = false;
        if (isMC) {
            if (process == "GG_ttn1") {
                if (Mu_type->at(n_mu) == 6)
                    truth_match = true;
            }
            else { // process != GG_ttn1
                if (Mu_type->at(n_mu) == 6 &&
                    (Mu_origin->at(n_mu) == 10 || Mu_origin->at(n_mu) == 12 || Mu_origin->at(n_mu) == 13))
                    truth_match = true;
            }
        }

        float weight = 1.; // Data weight = 1.
        if (isMC) {
            weight = normalization * baseline_weight;
        }

        // mll window: 60 < mll < 150
        h_baseline_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        h_baseline_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        // Reset weight for signal leptons
        if (isMC) {
            weight = normalization * signal_weight;
        }
        if (Mu_isSignal->at(n_mu)) {
            h_signal_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000., weight);
            h_signal_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        }

        if (isData) {
            // Use probe muons only
            if (!Mu_isZProbe->at(n_mu))
                continue;
        }
        else if (isMC) {
            if (process == "Zmumu") {
                // Use probe muons only
                if (!Mu_isZProbe->at(n_mu))
                    continue;
            }
            else if (process == "Zmumu_TandP_truth_matched") {
                // Use truth match on top of probe electron
                if (!Mu_isZProbe->at(n_mu))
                    continue;
                if (!truth_match)
                    continue;
            }
            else { // process == "GG_ttn1" || process == "ttbar" || process == "Zee_truth_match" || process == "Zmumu_truth_match"
                // Use truth match
                //cout << "process is " << process << endl;
                if (!truth_match)
                    continue;
            }
        }

        // Tag trigger matching
        if (tag_trigger_matched) {
            if (!Mu_isZProbe_TagTriggerMatched->at(n_mu))
                continue;
        }
/*
        if (isData) {
            h_baseline_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000.);
            h_baseline_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000.);

            if (Mu_isSignal->at(n_mu)) {
                h_signal_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000.);
                h_signal_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000.);
            }

            // Use probe muons only
            if (!Mu_isZProbe->at(n_mu))
                continue;
        }
        else if (isMC) {
            h_baseline_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000., weight);
            h_baseline_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight);

            if (Mu_isSignal->at(n_mu)) {
                h_signal_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000., normalization * signal_weight);
                h_signal_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., normalization * signal_weight);
            }

            if (process == "Zmumu") {
                // Use probe muons only
                if (!Mu_isZProbe->at(n_mu))
                    continue;
            }
            else if (process == "Zmumu_TandP_truth_matched") {
                // Use truth match on top of probe electron
                if (!Mu_isZProbe->at(n_mu))
                    continue;
                if (!truth_match)
                    continue;
            }
            else { // process == "GG_ttn1" || process == "ttbar" || process == "Zee_truth_match" || process == "Zmumu_truth_match"
                // Use truth match
                //cout << "process is " << process << endl;
                if (!truth_match)
                    continue;
            }
        }
*/
        // Reset weight for signal leptons
        if (isMC) {
            weight = normalization * baseline_weight;
        }

        // 1-dim histograms
        //h_baseline_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000., weight); // only dump 80 < mll < 100
        h_baseline_pt->Fill(Mu_pT->at(n_mu) / 1000., weight);
        h_baseline_eta->Fill(fabs(Mu_eta->at(n_mu)), weight);
        //h_baseline_d0pvtx->Fill(Mu_d0pvtx->at(n_mu), weight);
        //h_baseline_sigd0->Fill(Mu_sigd0->at(n_mu), weight);
        //h_baseline_d0err->Fill(Mu_d0pvtx->at(n_mu) / Mu_sigd0->at(n_mu), weight);
        //h_baseline_z0sinTheta->Fill(Mu_z0sinTheta->at(n_mu), weight);
        //h_baseline_ptvarcone20->Fill(Mu_ptvarcone20->at(n_mu) / 1000., weight);
        //h_baseline_ptvarcone30->Fill(Mu_ptvarcone30->at(n_mu) / 1000., weight);
        //h_baseline_topoetcone20->Fill(Mu_topoetcone20->at(n_mu) / 1000., weight);
        h_baseline_nJets->Fill(NJet, weight);
        h_baseline_dRjet->Fill(Mu_DR_closest_Jet->at(n_mu), weight);
        //h_baseline_Etmiss->Fill(Etmiss_TST_Et / 1000., weight);
        //h_baseline_meff->Fill(meff / 1000., weight);
        // 2-dim histograms
        h_baseline_pt_eta->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), weight);
        h_baseline_pt_nJets->Fill(Mu_pT->at(n_mu) / 1000., NJet, weight);
        h_baseline_pt_dRjet->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), weight);
        //h_baseline_pt_Etmiss->Fill(Mu_pT->at(n_mu) / 1000.,Etmiss_TST_Et / 1000., weight);
        //h_baseline_pt_meff->Fill(Mu_pT->at(n_mu) / 1000., meff / 1000., weight);
        // 3-dim histograms
        //h_baseline_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight); // only dump 80 < mll < 100
        h_baseline_pt_eta_dRjet->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), weight);
        h_baseline_pt_dRjet_mll->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        h_baseline_eta_dRjet_mll->Fill(fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        if (isData) {
            h_baseline_pt_eta_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), new_AvgMu * 1.0 / 1.09, weight);
            h_baseline_pt_dRjet_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), new_AvgMu * 1.0 / 1.09, weight);
            h_baseline_eta_dRjet_AvgMu->Fill(fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), new_AvgMu * 1.0 / 1.09, weight);
        }
        else if (isMC) {
            h_baseline_pt_eta_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), new_AvgMu, weight);
            h_baseline_pt_dRjet_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), new_AvgMu, weight);
            h_baseline_eta_dRjet_AvgMu->Fill(fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), new_AvgMu, weight);
        }

        // cut efficiency
        // Apply ptvarcone30/pt < 0.06
        if ((Mu_ptvarcone30->at(n_mu) / Mu_pT->at(n_mu)) < 0.06) {
            h_baseline_TrackIso_cut_pt->Fill(Mu_pT->at(n_mu) / 1000., weight);
        }
        // Apply |z0*sin(theta)| < 0.5
        if (fabs(Mu_z0sinTheta->at(n_mu)) < 0.5) {
            h_baseline_z0_cut_pt->Fill(Mu_pT->at(n_mu) / 1000., weight);
        }
        // Apply |d0sig| < 3
        if (fabs(Mu_sigd0->at(n_mu)) < 3) {
            h_baseline_sigd0_cut_pt->Fill(Mu_pT->at(n_mu) / 1000., weight);
        }

        // 3-dim cut efficiency
        // Apply ptvarcone30/pt < 0.06
        if ((Mu_ptvarcone30->at(n_mu) / Mu_pT->at(n_mu)) < 0.06) {
            h_baseline_TrackIso_cut_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        }
        // Apply |z0*sin(theta)| < 0.5
        if (fabs(Mu_z0sinTheta->at(n_mu)) < 0.5) {
            h_baseline_z0_cut_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        }
        // Apply |d0sig| < 3
        if (fabs(Mu_sigd0->at(n_mu)) < 3) {
            h_baseline_sigd0_cut_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
        }

        // Reset weight for signal leptons
        if (isMC) {
            weight = normalization * signal_weight;
        }
        // All the muons here are signal muons
        if (Mu_isSignal->at(n_mu)) {
            // 1-dim histograms
            //h_signal_mll->Fill(Mu_ZTandP_mll->at(n_mu) / 1000., weight); // only dump 80 < mll < 100
            h_signal_pt->Fill(Mu_pT->at(n_mu) / 1000., weight);
            h_signal_eta->Fill(fabs(Mu_eta->at(n_mu)), weight);
            //h_signal_d0pvtx->Fill(Mu_d0pvtx->at(n_mu), weight);
            //h_signal_sigd0->Fill(Mu_sigd0->at(n_mu), weight);
            //h_signal_d0err->Fill(Mu_d0pvtx->at(n_mu) / Mu_sigd0->at(n_mu), weight);
            //h_signal_z0sinTheta->Fill(Mu_z0sinTheta->at(n_mu), weight);
            //h_signal_ptvarcone20->Fill(Mu_ptvarcone20->at(n_mu) / 1000., weight);
            //h_signal_ptvarcone30->Fill(Mu_ptvarcone30->at(n_mu) / 1000., weight);
            //h_signal_topoetcone20->Fill(Mu_topoetcone20->at(n_mu) / 1000., weight);
            h_signal_nJets->Fill(NJet, weight);
            h_signal_dRjet->Fill(Mu_DR_closest_Jet->at(n_mu), weight);
            //h_signal_Etmiss->Fill(Etmiss_TST_Et / 1000., weight);
            //h_signal_meff->Fill(meff / 1000., weight);
            // 2-dim histograms
            h_signal_pt_eta->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), weight);
            h_signal_pt_nJets->Fill(Mu_pT->at(n_mu) / 1000., NJet, weight);
            h_signal_pt_dRjet->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), weight);
            //h_signal_pt_Etmiss->Fill(Mu_pT->at(n_mu) / 1000.,Etmiss_TST_Et / 1000., weight);
            //h_signal_pt_meff->Fill(Mu_pT->at(n_mu) / 1000., meff / 1000., weight);
            // 3-dim histograms
            //h_signal_pt_eta_mll->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_ZTandP_mll->at(n_mu) / 1000., weight); // only dump 80 < mll < 100
            h_signal_pt_eta_dRjet->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), weight);
            h_signal_pt_dRjet_mll->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
            h_signal_eta_dRjet_mll->Fill(fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), Mu_ZTandP_mll->at(n_mu) / 1000., weight);
            if (isData) {
                h_signal_pt_eta_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), new_AvgMu * 1.0 / 1.09, weight);
                h_signal_pt_dRjet_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), new_AvgMu * 1.0 / 1.09, weight);
                h_signal_eta_dRjet_AvgMu->Fill(fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), new_AvgMu * 1.0 / 1.09, weight);
            }
            else if (isMC) {
                h_signal_pt_eta_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., fabs(Mu_eta->at(n_mu)), new_AvgMu, weight);
                h_signal_pt_dRjet_AvgMu->Fill(Mu_pT->at(n_mu) / 1000., Mu_DR_closest_Jet->at(n_mu), new_AvgMu, weight);
                h_signal_eta_dRjet_AvgMu->Fill(fabs(Mu_eta->at(n_mu)), Mu_DR_closest_Jet->at(n_mu), new_AvgMu, weight);
            }
        }
    }
}



void ytRealLeptonsEfficiency::make_baseline_mll_plot()
{
    baseline_mll_plot = new TCanvas("baseline_mll", "Baseline M_{ll}", 500, 500);
    baseline_mll_plot->SetLeftMargin(0.12);
    baseline_mll_plot->cd();

    TH1F *h_mll = (TH1F *)h_baseline_pt_eta_mll->ProjectionZ("");
    h_mll->SetLineColor(kBlack);
    h_mll->SetMarkerColor(kBlack);
    h_mll->SetMarkerStyle(kFullCircle);
    h_mll->SetStats(kFALSE);
    h_mll->GetXaxis()->SetRangeUser(60, 150);
    h_mll->SetTitle("");
    h_mll->SetXTitle("M_{ll}");
    h_mll->SetYTitle("Number of probe leptos");
    h_mll->GetYaxis()->SetTitleOffset(1.5);
    h_mll->Draw();

    double elec_max = h_mll->GetMaximum();

    // Create a dummy histogram in order to fill color between 80 < mll < 100
    TH1F *h_80_100 = new TH1F("h_80_100", "h_80_100", 1, 80, 100);
    h_80_100->SetFillColor(kYellow);
    h_80_100->SetFillColorAlpha(kYellow, 0.35);
    h_80_100->SetFillStyle(1001); // solid
    h_80_100->SetLineColor(kYellow);
    h_80_100->Fill(90, elec_max * 10);
    h_80_100->Draw("hist,same");

    wk()->addOutput(baseline_mll_plot);
}



void ytRealLeptonsEfficiency::make_cut_efficiency_plot()
{
    cut_efficiency_plot = new TCanvas("cut_efficiency", "Cut efficiency comparison", 500, 500);
    cut_efficiency_plot->SetLeftMargin(0.12);
    cut_efficiency_plot->cd();
    //gStyle->SetOptStat(0);

    h_eff_pt->SetLineColor(kBlack);
    h_eff_pt->SetMarkerColor(kBlack);
    h_eff_pt->SetMarkerStyle(kFullCircle);
    h_eff_pt->SetStats(kFALSE);
    h_eff_pt->SetMaximum(1.1);
    h_eff_pt->SetMinimum(0.);
    h_eff_pt->SetTitle("");
    h_eff_pt->GetYaxis()->SetTitleOffset(1.5);
    h_eff_pt->Draw();

    if (lepton == "electron") {
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetLineColor(kRed);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetMarkerColor(kRed);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->SetMarkerStyle(kFullSquare);
        h_cut_eff_LooseAndBLayerLLH_to_MediumLLH->Draw("same");

        h_cut_eff_CaloIso->SetLineColor(kGreen);
        h_cut_eff_CaloIso->SetMarkerColor(kGreen);
        h_cut_eff_CaloIso->SetMarkerStyle(kFullTriangleUp);
        h_cut_eff_CaloIso->Draw("same");
    }

    h_cut_eff_TrackIso->SetLineColor(kBlue);
    h_cut_eff_TrackIso->SetMarkerColor(kBlue);
    h_cut_eff_TrackIso->SetMarkerStyle(kFullTriangleDown);
    h_cut_eff_TrackIso->Draw("same");

    h_cut_eff_z0->SetLineColor(kMagenta);
    h_cut_eff_z0->SetMarkerColor(kMagenta);
    h_cut_eff_z0->SetMarkerStyle(kFullDiamond);
    h_cut_eff_z0->Draw("same");

    if (lepton == "muon") {
        h_cut_eff_sigd0->SetLineColor(kCyan);
        h_cut_eff_sigd0->SetMarkerColor(kCyan);
        h_cut_eff_sigd0->SetMarkerStyle(kFullCross);
        h_cut_eff_sigd0->Draw("same");
    }

    float X_max = h_eff_pt->GetXaxis()->GetXmax();
    float X_min = h_eff_pt->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.5, 0.2, 0.8, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry("h_eff_pt", "Real Efficiency");
    if (lepton == "electron") {
        leg->AddEntry("h_cut_eff_LooseAndBLayerLLH_to_MediumLLH", "LooseLH to MediumLH");
        leg->AddEntry("h_cut_eff_CaloIso", "Calo isolation");
    }
    leg->AddEntry("h_cut_eff_TrackIso", "Track isolation");
    leg->AddEntry("h_cut_eff_z0", "|z_{0} sin(#theta)| < 0.5");
    if (lepton == "muon") {
        leg->AddEntry("h_cut_eff_sigd0", "|d_{0}/#sigma_{d0}| < 3");
    }
    leg->Draw();

    wk()->addOutput(cut_efficiency_plot);
}



void ytRealLeptonsEfficiency::make_electron_real_efficiency_plot()
{
    // For electrons, there are 3 eta ranges: [0, 0.8], [0.8, 1.37], [1.52, 2]
    // Corresponding to eta bins: [1, 3], [4, 5], and [7, 8]
    // In order to get the correct bin, we have to use lower_edge + 0.01 and upper_edge - 0.01 to find the bin.
    // This is because FindBin() sometimes assign the lower_edge (upper_edge) to the correct_bin-1 (correct_bin+1).
    TH1F *h_baseline_eta_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_1", h_baseline_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(0.8 - 0.01));
    TH1F *h_baseline_eta_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_2", h_baseline_pt_eta->GetYaxis()->FindBin(0.8 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(1.37 - 0.01));
    TH1F *h_baseline_eta_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_3", h_baseline_pt_eta->GetYaxis()->FindBin(1.52 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(2.0 - 0.01));

    TH1F *h_signal_eta_range_1 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_1", h_signal_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(0.8 - 0.01));
    TH1F *h_signal_eta_range_2 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_2", h_signal_pt_eta->GetYaxis()->FindBin(0.8 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(1.37 - 0.01));
    TH1F *h_signal_eta_range_3 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_3", h_signal_pt_eta->GetYaxis()->FindBin(1.52 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(2.0 - 0.01));

    TH1F *ratio_1 = (TH1F *)h_signal_eta_range_1->Clone();
    TH1F *ratio_2 = (TH1F *)h_signal_eta_range_2->Clone();
    TH1F *ratio_3 = (TH1F *)h_signal_eta_range_3->Clone();

    ratio_1->SetName("ratio_1");
    ratio_2->SetName("ratio_2");
    ratio_3->SetName("ratio_3");

    ratio_1->Divide(h_baseline_eta_range_1);
    ratio_2->Divide(h_baseline_eta_range_2);
    ratio_3->Divide(h_baseline_eta_range_3);

    real_efficiency_plot = new TCanvas("real_efficiency_plot", "Electron Real Efficiency", 500, 500);
    real_efficiency_plot->SetLeftMargin(0.12);
    real_efficiency_plot->cd();

    ratio_1->SetLineColor(kBlack);
    ratio_1->SetMarkerColor(kBlack);
    ratio_1->SetMarkerStyle(kFullCircle);
    ratio_1->SetXTitle("p_{T} [GeV]");
    ratio_1->SetYTitle("Electron Real Efficiency");
    ratio_1->GetYaxis()->SetTitleOffset(1.5);
    ratio_1->SetStats(kFALSE);
    ratio_1->SetMaximum(1.1);
    ratio_1->SetMinimum(0.);
    ratio_1->SetTitle("");
    ratio_1->Draw();

    ratio_2->SetLineColor(kRed);
    ratio_2->SetMarkerColor(kRed);
    ratio_2->SetMarkerStyle(kFullSquare);
    ratio_2->Draw("same");

    ratio_3->SetLineColor(kBlue);
    ratio_3->SetMarkerColor(kBlue);
    ratio_3->SetMarkerStyle(kFullTriangleUp);
    ratio_3->Draw("same");

    float X_max = ratio_1->GetXaxis()->GetXmax();
    float X_min = ratio_1->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.5, 0.2, 0.7, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry("ratio_1", "0 < |#eta| < 0.8");
    leg->AddEntry("ratio_2", "0.8 < |#eta| < 1.37");
    leg->AddEntry("ratio_3", "1.52 < |#eta| < 2.0");
    leg->Draw();

    wk()->addOutput(real_efficiency_plot);
}



void ytRealLeptonsEfficiency::make_muon_real_efficiency_plot()
{
    // For muons, there are 4 eta ranges: [0, 0.6], [0.6, 1.2], [1.2, 1.8], [1.8, 2.5]
    // Corresponding to eta bins: [1, 1], [2, 2], [3, 3], and [4, 4]
    // In order to get the correct bin, we have to use lower_edge + 0.01 and upper_edge - 0.01 to find the bin.
    // This is because FindBin() sometimes assign the lower_edge (upper_edge) to the correct_bin-1 (correct_bin+1).
    TH1F *h_baseline_eta_range_1 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_1", h_baseline_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(0.6 - 0.01));
    TH1F *h_baseline_eta_range_2 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_2", h_baseline_pt_eta->GetYaxis()->FindBin(0.6 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(1.2 - 0.01));
    TH1F *h_baseline_eta_range_3 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_3", h_baseline_pt_eta->GetYaxis()->FindBin(1.2 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(1.8 - 0.01));
    TH1F *h_baseline_eta_range_4 = (TH1F *)h_baseline_pt_eta->ProjectionX("h_eta_range_4", h_baseline_pt_eta->GetYaxis()->FindBin(1.8 + 0.01), h_baseline_pt_eta->GetYaxis()->FindBin(2.5 - 0.01));

    TH1F *h_signal_eta_range_1 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_1", h_signal_pt_eta->GetYaxis()->FindBin(0. + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(0.6 - 0.01));
    TH1F *h_signal_eta_range_2 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_2", h_signal_pt_eta->GetYaxis()->FindBin(0.6 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(1.2 - 0.01));
    TH1F *h_signal_eta_range_3 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_3", h_signal_pt_eta->GetYaxis()->FindBin(1.2 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(1.8 - 0.01));
    TH1F *h_signal_eta_range_4 = (TH1F *)h_signal_pt_eta->ProjectionX("h_signal_eta_range_4", h_signal_pt_eta->GetYaxis()->FindBin(1.8 + 0.01), h_signal_pt_eta->GetYaxis()->FindBin(2.5 - 0.01));

    TH1F *ratio_1 = (TH1F *)h_signal_eta_range_1->Clone();
    TH1F *ratio_2 = (TH1F *)h_signal_eta_range_2->Clone();
    TH1F *ratio_3 = (TH1F *)h_signal_eta_range_3->Clone();
    TH1F *ratio_4 = (TH1F *)h_signal_eta_range_4->Clone();

    ratio_1->SetName("ratio_1");
    ratio_2->SetName("ratio_2");
    ratio_3->SetName("ratio_3");
    ratio_4->SetName("ratio_4");

    ratio_1->Divide(h_baseline_eta_range_1);
    ratio_2->Divide(h_baseline_eta_range_2);
    ratio_3->Divide(h_baseline_eta_range_3);
    ratio_4->Divide(h_baseline_eta_range_4);

    real_efficiency_plot = new TCanvas("real_efficiency_plot", "Muon Real Efficiency", 500, 500);
    real_efficiency_plot->SetLeftMargin(0.12);
    real_efficiency_plot->cd();

    ratio_1->SetLineColor(kBlack);
    ratio_1->SetMarkerColor(kBlack);
    ratio_1->SetMarkerStyle(kFullCircle);
    ratio_1->SetXTitle("p_{T} [GeV]");
    ratio_1->SetYTitle("Muon Real Efficiency");
    ratio_1->GetYaxis()->SetTitleOffset(1.5);
    ratio_1->SetStats(kFALSE);
    ratio_1->SetMaximum(1.1);
    ratio_1->SetMinimum(0.);
    ratio_1->SetTitle("");
    ratio_1->Draw();

    ratio_2->SetLineColor(kRed);
    ratio_2->SetMarkerColor(kRed);
    ratio_2->SetMarkerStyle(kFullSquare);
    ratio_2->Draw("same");

    ratio_3->SetLineColor(kBlue);
    ratio_3->SetMarkerColor(kBlue);
    ratio_3->SetMarkerStyle(kFullTriangleUp);
    ratio_3->Draw("same");

    ratio_4->SetLineColor(kGreen);
    ratio_4->SetMarkerColor(kGreen);
    ratio_4->SetMarkerStyle(kFullTriangleDown);
    ratio_4->Draw("same");

    float X_max = ratio_1->GetXaxis()->GetXmax();
    float X_min = ratio_1->GetXaxis()->GetXmin();

    TLine *line = new TLine(X_min, 1., X_max,1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    TLegend *leg = new TLegend(0.5, 0.2, 0.7, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry("ratio_1", "0 < |#eta| < 0.6");
    leg->AddEntry("ratio_2", "0.6 < |#eta| < 1.2");
    leg->AddEntry("ratio_3", "1.2 < |#eta| < 1.8");
    leg->AddEntry("ratio_4", "1.8 < |#eta| < 2.5");
    leg->Draw();

    wk()->addOutput(real_efficiency_plot);
}
