#define YT_SKIM_CXX
#include "ytNUHM2Analysis/yt_skim.h"

ClassImp(yt_skim)

yt_skim::yt_skim()
{
    // Physics object variables
    EE_pair     = new vector<ZTandP_pair>();
    MuMu_pair   = new vector<ZTandP_pair>();

    El_isBaseline   = new vector<bool>();
    El_isSignal     = new vector<bool>();
    El_isZTag       = new vector<bool>();
    El_isZProbe     = new vector<bool>();
    El_Tag_index    = new vector<int>();
    El_Probe_index  = new vector<int>();
    El_ZTandP_mll   = new vector<double>();
    El_isZProbe_TagTriggerMatched = new vector<bool>();
    El_Tag_Trigger_Matched = new vector<bool>();
    El_Probe_Trigger_Matched = new vector<bool>();
    El_DR_closest_Jet = new vector<double>();

    Mu_isBaseline   = new vector<bool>();
    Mu_isSignal     = new vector<bool>();
    Mu_isZTag       = new vector<bool>();
    Mu_isZProbe     = new vector<bool>();
    Mu_Tag_index    = new vector<int>();
    Mu_Probe_index  = new vector<int>();
    Mu_ZTandP_mll   = new vector<double>();
    Mu_isZProbe_TagTriggerMatched = new vector<bool>();
    Mu_Tag_Trigger_Matched = new vector<bool>();
    Mu_Probe_Trigger_Matched = new vector<bool>();
    Mu_DR_closest_Jet = new vector<double>();

    //Jet_isBaseline  = new vector<bool>();
    //Jet_isSignal    = new vector<bool>();
    //bJet_isSignal   = new vector<bool>();

    // Events variables
    normalization = 1.;
    baseline_mll = 0.;
    signal_mll = 0.;
    //jets_mll = 0.;
    //Ht = 0.;
    //meff = 0.;

    // Counter
    N_at_least_two_baseline_lepton_events = 0;
    N_exactly_two_baseline_lepton_events = 0;
    N_tot_baseline_electrons = 0;
    N_tot_baseline_muons = 0;
    N_tot_TandP_pair_electron = 0;
    N_tot_TandP_pair_muon = 0;
}

void yt_skim::clear_all_object_vectors()
{
    // clear all vectors

    vec_elec.clear();
    vec_muon.clear();
    vec_lept.clear();
    vec_jets.clear();

    vec_baseline_elec.clear();
    vec_baseline_muon.clear();
    vec_baseline_lept.clear();
    vec_baseline_jets.clear();

    vec_signal_elec.clear();
    vec_signal_muon.clear();
    vec_signal_lept.clear();
    vec_signal_jets.clear();

    vec_signal_jets_no_eta_cut.clear();
    vec_signal_bjet.clear();
}

void yt_skim::fill_all_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets,
                                      vector<Electron> baseline_elec, vector<Muon> baseline_muon, vector<Lepton> baseline_lept, vector<Jet> baseline_jets,
                                      vector<Electron> signal_elec, vector<Muon> signal_muon, vector<Lepton> signal_lept, vector<Jet> signal_jets)
{
    // in event loop:
    // copy vectors from yt_selector.
    this->copy_raw_object_vectors(elec, muon, lept, jets);
    this->copy_baseline_object_vectors(baseline_elec, baseline_muon, baseline_lept, baseline_jets);
    this->copy_signal_object_vectors(signal_elec, signal_muon, signal_lept, signal_jets);
    this->fill_signal_jets_no_eta_cut(jets);
    this->fill_signal_bjet(signal_jets);
}

void yt_skim::copy_raw_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets)
{
    vec_elec = elec;
    vec_muon = muon;
    vec_lept = lept;
    vec_jets = jets;
}

void yt_skim::copy_baseline_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets)
{
    vec_baseline_elec = elec;
    vec_baseline_muon = muon;
    vec_baseline_lept = lept;
    vec_baseline_jets = jets;
}

void yt_skim::copy_signal_object_vectors(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets)
{
    vec_signal_elec = elec;
    vec_signal_muon = muon;
    vec_signal_lept = lept;
    vec_signal_jets = jets;
}

void yt_skim::fill_signal_jets_no_eta_cut(vector<Jet> jets)
{
    for (auto & jet_itr : jets) {
        if (jet_itr.get_pt() > 20000. &&
            jet_itr.get_passOR() == true) {
            if (jet_itr.get_pt() < 60000. &&
                fabs(jet_itr.get_eta()) < 2.4 &&
                jet_itr.get_JVT() < 0.59) {
                continue;
            }
            else {
                if (jet_itr.get_quality() < 0.5)
                    vec_signal_jets_no_eta_cut.push_back(jet_itr);
            }
        }
    }
}

void yt_skim::fill_signal_bjet(vector<Jet> signal_jets)
{
    for (auto & signal_jets_itr : signal_jets) {
        if (signal_jets_itr.get_isBjet() == true)
            vec_signal_bjet.push_back(signal_jets_itr);
    }
}

void yt_skim::initialize(TTree *tree)
{
    // Declare the output
    TString output_path("/raid05/users/shen/Ximo_ntuples/v50/Skimmed/");
    if (isMC) {
        output_file = TFile::Open(output_path + "/MC_probes_" + process + ".root", "RECREATE");
    }
    else if (isData) {
        output_file = TFile::Open(output_path + "/data_probes_" + process + ".root", "RECREATE");
    }
    output_tree = (TTree*)tree->CloneTree(0);

    // Building the new branches
    // Normalization
    output_tree->Branch("normalization", &normalization);
    output_tree->Branch("luminosity", &luminosity);
    output_tree->Branch("cross_section", &cross_section);
    output_tree->Branch("k_factor", &k_factor);
    output_tree->Branch("filter_efficiency", &filter_efficiency);
    output_tree->Branch("cross_section_kfactor_efficiency", &cross_section_kfactor_efficiency);
    output_tree->Branch("event_weight_sum", &event_weight_sum);
    output_tree->Branch("event_weight", &event_weight);
    output_tree->Branch("pileup_weight", &pileup_weight);
    output_tree->Branch("baseline_weight", &baseline_weight);
    output_tree->Branch("signal_weight", &signal_weight);
    output_tree->Branch("run_number", &run_number);
    output_tree->Branch("new_AvgMu", &new_AvgMu);
    // Physics object variables
    output_tree->Branch("EE_pair", &EE_pair);
    output_tree->Branch("MuMu_pair", &MuMu_pair);

    output_tree->Branch("El_isBaseline",                &El_isBaseline);
    output_tree->Branch("El_isSignal",                  &El_isSignal);
    output_tree->Branch("El_isZTag",                    &El_isZTag);
    output_tree->Branch("El_isZProbe",                  &El_isZProbe);
    output_tree->Branch("El_Tag_index",                 &El_Tag_index);
    output_tree->Branch("El_Probe_index",               &El_Probe_index);
    output_tree->Branch("El_ZTandP_mll",                &El_ZTandP_mll);
    output_tree->Branch("El_isZProbe_TagTriggerMatched",&El_isZProbe_TagTriggerMatched);
    output_tree->Branch("El_Tag_Trigger_Matched",       &El_Tag_Trigger_Matched);
    output_tree->Branch("El_Probe_Trigger_Matched",     &El_Probe_Trigger_Matched);
    output_tree->Branch("El_DR_closest_Jet",            &El_DR_closest_Jet);

    output_tree->Branch("Mu_isBaseline",                &Mu_isBaseline);
    output_tree->Branch("Mu_isSignal",                  &Mu_isSignal);
    output_tree->Branch("Mu_isZTag",                    &Mu_isZProbe);
    output_tree->Branch("Mu_isZProbe",                  &Mu_isZProbe);
    output_tree->Branch("Mu_Tag_index",                 &Mu_Tag_index);
    output_tree->Branch("Mu_Probe_index",               &Mu_Probe_index);
    output_tree->Branch("Mu_ZTandP_mll",                &Mu_ZTandP_mll);
    output_tree->Branch("Mu_isZProbe_TagTriggerMatched",&Mu_isZProbe_TagTriggerMatched);
    output_tree->Branch("Mu_Tag_Trigger_Matched",       &Mu_Tag_Trigger_Matched);
    output_tree->Branch("Mu_Probe_Trigger_Matched",     &Mu_Probe_Trigger_Matched);
    output_tree->Branch("Mu_DR_closest_Jet",            &Mu_DR_closest_Jet);

    //output_tree->Branch("Jet_isBaseline",     &Jet_isBaseline);
    //output_tree->Branch("Jet_isSignal",   &Jet_isSignal);
    //output_tree->Branch("bJet_isSignal",  &bJet_isSignal);

    // Events numbers variables
    output_tree->Branch("baseline_mll", &baseline_mll);
    output_tree->Branch("signal_mll",   &signal_mll);
    //output_tree->Branch("jets_mll",   &jets_mll);
    //output_tree->Branch("Ht",         &Ht);
    //output_tree->Branch("meff",       &meff);

}

void yt_skim::execute(vector<Electron> elec, vector<Muon> muon, vector<Lepton> lept, vector<Jet> jets,
                      vector<Electron> baseline_elec, vector<Muon> baseline_muon, vector<Lepton> baseline_lept, vector<Jet> baseline_jets,
                      vector<Electron> signal_elec, vector<Muon> signal_muon, vector<Lepton> signal_lept, vector<Jet> signal_jets)
{
    // clear all the vector members
    this->clear_all_object_vectors();
    // fill all the vector members
    this->fill_all_object_vectors(elec, muon, lept, jets,
                                  baseline_elec, baseline_muon, baseline_lept, baseline_jets,
                                  signal_elec, signal_muon, signal_lept, signal_jets);
    // calculate normalization for this event
    if (isMC)
        normalization = luminosity * (cross_section_kfactor_efficiency * 1000.) * (event_weight / event_weight_sum) * pileup_weight;
    else if (isData)
        normalization = 1.;

    if (baseline_lept.size() >= 2) {
        N_at_least_two_baseline_lepton_events++;
        if (baseline_lept.size() == 2)
            N_exactly_two_baseline_lepton_events++;

        N_tot_baseline_electrons += vec_baseline_elec.size();
        N_tot_baseline_muons += vec_baseline_muon.size();

        this->reset_vectors();
        this->calculate_new_variables();
        this->Z_tag_and_probe();

        // fill all new branches
        output_tree->Fill();
    }
}

void yt_skim::finalize()
{
    cout << "Write into file" << endl;
    output_file->cd();
    output_tree->Write();
    output_file->Close();
}

void yt_skim::reset_vectors()
{
/*
    // Create dummy instance of ZTandP_pair
    // Use dummy to file EE_pair and MuMu_pair
    ZTandP_pair dummy;
    dummy.tag_index = -1;
    dummy.tag_flavor = 0;
    dummy.tag_charge = 0;
    dummy.tag_pt = 0.;
    dummy.tag_eta = 10.;
    dummy.tag_phi = 10.;
    dummy.tag_is_trigger_matched = false;
    dummy.tag_id_SFw = 1.;
    dummy.tag_iso_SFw = 1.;
    dummy.tag_trigger_SFw = 1.;
    //
    dummy.probe_index = -1;
    dummy.probe_flavor = 0;
    dummy.probe_charge = 0;
    dummy.probe_pt = 0.;
    dummy.probe_eta = 10.;
    dummy.probe_phi = 10.;
    dummy.probe_is_trigger_matched = false;
    dummy.probe_id_SFw = 1.;
    dummy.probe_iso_SFw = 1.;
    dummy.probe_trigger_SFw = 1.;
    dummy.probe_is_baseline = false;
    dummy.probe_is_signal = false;
    //dummy.probe_deltaR_closest_jet = 100;
    dummy.baseline_mll = 0.;
    dummy.signal_mll = 0.;
*/
    // reset vectors in event loop
    EE_pair->clear();
    MuMu_pair->clear();
/*
    for (unsigned int i = 0; i < vec_elec.size(); i++) {
        EE_pair->push_back(dummy);
    }

    for (unsigned int i = 0; i < vec_muon.size(); i++) {
        MuMu_pair->push_back(dummy);
    }
*/
    El_isBaseline->clear();
    El_isSignal->clear();
    El_isZTag->clear();
    El_isZProbe->clear();
    El_Tag_index->clear();
    El_Probe_index->clear();
    El_ZTandP_mll->clear();
    El_isZProbe_TagTriggerMatched->clear();
    El_Tag_Trigger_Matched->clear();
    El_Probe_Trigger_Matched->clear();
    El_DR_closest_Jet->clear();

    Mu_isBaseline->clear();
    Mu_isSignal->clear();
    Mu_isZTag->clear();
    Mu_isZProbe->clear();
    Mu_Tag_index->clear();
    Mu_Probe_index->clear();
    Mu_ZTandP_mll->clear();
    Mu_isZProbe_TagTriggerMatched->clear();
    Mu_Tag_Trigger_Matched->clear();
    Mu_Probe_Trigger_Matched->clear();
    Mu_DR_closest_Jet->clear();

    //Jet_isBaseline->clear();
    //Jet_isSignal->clear();
    //bJet_isSignal->clear();

    for (unsigned int i = 0; i < vec_elec.size(); i++) {
        El_isBaseline->push_back(false);
        El_isSignal->push_back(false);
        El_isZTag->push_back(false);
        El_isZProbe->push_back(false);
        El_Tag_index->push_back(100);
        El_Probe_index->push_back(100);
        El_ZTandP_mll->push_back(0);
        El_isZProbe_TagTriggerMatched->push_back(false);
        El_Tag_Trigger_Matched->push_back(false);
        El_Probe_Trigger_Matched->push_back(false);
        El_DR_closest_Jet->push_back(100);
    }


    for (unsigned int i = 0; i < vec_muon.size(); i++) {
        Mu_isBaseline->push_back(false);
        Mu_isSignal->push_back(false);
        Mu_isZTag->push_back(false);
        Mu_isZProbe->push_back(false);
        Mu_Tag_index->push_back(100);
        Mu_Probe_index->push_back(100);
        Mu_ZTandP_mll->push_back(0);
        Mu_isZProbe_TagTriggerMatched->push_back(false);
        Mu_Tag_Trigger_Matched->push_back(false);
        Mu_Probe_Trigger_Matched->push_back(false);
        Mu_DR_closest_Jet->push_back(100);
    }
/*
    for (unsigned int i = 0; i < vec_jets.size(); i++) {
        Jet_isBaseline->push_back(false);
        Jet_isSignal->push_back(false);
        bJet_isSignal->push_back(false);
    } 
*/
}

void yt_skim::calculate_new_variables()
{
    baseline_mll = calculate_mll(vec_baseline_lept);
    signal_mll = calculate_mll(vec_signal_lept);
    //jets_mll = calculate_mjj(vec_signal_jets);
    //Ht = calculate_Ht(vec_signal_lept, vec_signal_jets);
    //meff = calculate_Meff(Ht, Etmiss);

    // Filling the baseline / signal objects indicators / DRjets
    for (auto & elec_itr : vec_baseline_elec) {
        int current_index = elec_itr.get_index();
        El_isBaseline->at(current_index) = true;
        El_DR_closest_Jet->at(current_index) = calculate_deltaR(elec_itr, vec_signal_jets_no_eta_cut);
    }

    for (auto & muon_itr : vec_baseline_muon) {
        int current_index = muon_itr.get_index();
        Mu_isBaseline->at(current_index) = true;
        Mu_DR_closest_Jet->at(current_index) = calculate_deltaR(muon_itr, vec_signal_jets_no_eta_cut);
    }

    for (auto & elec_itr : vec_signal_elec) {
        El_isSignal->at(elec_itr.get_index()) = true;
    }

    for (auto & muon_itr : vec_signal_muon) {
        Mu_isSignal->at(muon_itr.get_index()) = true;
    }
/*
    for (auto & jet_itr : vec_baseline_jets) {
        Jet_isBaseline->at(jet_itr.get_index()) = true;
    }

    for (auto & jet_itr : vec_signal_jets) {
        Jet_isSignal->at(jet_itr.get_index()) = true;
    }

    for (auto & bjet_itr : vec_signal_bjet) {
        bJet_isSignal->at(bjet_itr.get_index()) = true;
    }
*/
}

/*
double yt_skim::calculate_mll(vector<Lepton> leptons)
{
    // calculate the invariant mass of the 2 leading leptons
    if (leptons.size() < 2)
        return -999.;

    TLorentzVector tlv_lept_0, tlv_lept_1;
    if (leptons.at(0).get_flavor() == 11)
        tlv_lept_0.SetPtEtaPhiE(leptons.at(0).get_pt(), leptons.at(0).get_eta(), leptons.at(0).get_phi(), leptons.at(0).get_E());
    else if (leptons.at(0).get_flavor() == 13)
        tlv_lept_0.SetPtEtaPhiM(leptons.at(0).get_pt(), leptons.at(0).get_eta(), leptons.at(0).get_phi(), Mu_Mass);

    if (leptons.at(0).get_flavor() == 11)
        tlv_lept_1.SetPtEtaPhiE(leptons.at(1).get_pt(), leptons.at(1).get_eta(), leptons.at(1).get_phi(), leptons.at(1).get_E());
    else if (leptons.at(0).get_flavor() == 13)
        tlv_lept_1.SetPtEtaPhiM(leptons.at(1).get_pt(), leptons.at(1).get_eta(), leptons.at(1).get_phi(), Mu_Mass);

    double mll = (tlv_lept_0 + tlv_lept_1).M();
    return mll;
}

double yt_skim::calculate_mjj(vector<Jet> jets)
{
    // calculate the invariant mass of the 2 leading jets
    if (jets.size() < 2)
        return -999.;

    TLorentzVector tlv_jets_0, tlv_jets_1;
    tlv_jets_0.SetPtEtaPhiE(jets.at(0).get_pt(), jets.at(0).get_eta(), jets.at(0).get_phi(), jets.at(0).get_E());
    tlv_jets_1.SetPtEtaPhiE(jets.at(1).get_pt(), jets.at(1).get_eta(), jets.at(1).get_phi(), jets.at(1).get_E());
    double mjj = (tlv_jets_0 + tlv_jets_1).M();
    return mjj;
}

double yt_skim::calculate_Ht(vector<Lepton> leptons, vector<Jet> jets)
{
    double Ht = 0;
    for (auto & lep : leptons) {
        Ht = Ht + lep.get_pt();
    }
    for (auto & jet : jets) {
        Ht = Ht + jet.get_pt();
    }
    return Ht;
}

double yt_skim::calculate_Meff(double Ht, double Etmiss)
{
    // the scalar sum of the signal leptons pT, all signal jets pT, and ETmiss
    return  Ht + Etmiss;
}
*/
double yt_skim::calculate_deltaR(Lepton lept, vector<Jet> jets)
{
    // Get the deltaR between the considered lepton and the closest jet

    double deltaR = 100.;
    int flavor = lept.get_flavor();
    for (auto & jet_itr : jets) {
        // computing deltaR
        TLorentzVector lep_tlv;
        TLorentzVector jet_tlv;
        if (flavor == 11)
            lep_tlv.SetPtEtaPhiE(lept.get_pt(), lept.get_eta(), lept.get_phi(), lept.get_E());
        else if (flavor == 13)
            lep_tlv.SetPtEtaPhiM(lept.get_pt(), lept.get_eta(), lept.get_phi(), lept.get_M());
        jet_tlv.SetPtEtaPhiE(jet_itr.get_pt(), jet_itr.get_eta(), jet_itr.get_phi(), jet_itr.get_E());
        double lep_rapidity = lep_tlv.Rapidity();
        double jet_rapidity = jet_tlv.Rapidity();
        double Deta = fabs(lep_rapidity - jet_rapidity);
        double Dphi = fabs(lept.get_phi() - jet_itr.get_phi());
        if (Dphi > TMath::Pi()) Dphi = 2 * TMath::Pi() - Dphi;
        double DR = sqrt(pow(Deta, 2) + pow(Dphi, 2));

        // DR cut
        if ( DR < deltaR ) {
            deltaR = DR;
        }
    }

    return deltaR;
}

void yt_skim::Z_tag_and_probe()
{
    // Z T&P method

    // Electrons:
    for (auto & probe_elec_itr : vec_baseline_elec) {
        double best_match_mll = 1000000.;
        ZTandP_pair candidate;
        for (auto & tag_elec_itr : vec_signal_elec) {
            if (tag_elec_itr.get_index() == probe_elec_itr.get_index())
                continue;
            if (tag_elec_itr.get_pt() < tag_pt_threshold)
                continue;
            if (fabs(tag_elec_itr.get_eta()) > 2.)
                continue;
            // Opposite Charge requirement
            if (tag_elec_itr.get_charge() == probe_elec_itr.get_charge())
                continue;
            
            candidate.event_number = event_number;

            candidate.tag_index = tag_elec_itr.get_index();
            candidate.tag_flavor = 11;
            candidate.tag_charge = tag_elec_itr.get_charge();
            candidate.tag_pt = tag_elec_itr.get_pt();
            candidate.tag_eta = tag_elec_itr.get_eta();
            candidate.tag_phi = tag_elec_itr.get_phi();
            candidate.tag_id_SFw = tag_elec_itr.get_SFwMediumLH();
            candidate.tag_iso_SFw = tag_elec_itr.get_IsoSFwMediumLH();
            candidate.tag_trigger_SFw = 0.; // no information

            candidate.probe_index = probe_elec_itr.get_index();
            candidate.probe_flavor = 11;
            candidate.probe_charge = probe_elec_itr.get_charge();
            candidate.probe_pt = probe_elec_itr.get_pt();
            candidate.probe_eta = probe_elec_itr.get_eta();
            candidate.probe_phi = probe_elec_itr.get_phi();
            candidate.probe_id_SFw = probe_elec_itr.get_SFwLooseAndBLayerLH(); // baseline require LooseAndBLayerLH
            candidate.probe_iso_SFw = probe_elec_itr.get_IsoSFwMediumLH();
            candidate.probe_trigger_SFw = 0.; // no information
            candidate.probe_is_baseline = probe_elec_itr.get_baseline();
            candidate.probe_is_signal = probe_elec_itr.get_isSignal();

            // Check does the tag trigger matched
            // electrons:
            // 2015: HLT_e24_lhmedium_L1EM20VH
            // 2016: HLT_e26_lhtight_nod0_ivarloose
            bool is_tag_trigger_matched = false;
            bool is_probe_trigger_matched = false;
            if (run_number < 290000) { // 2015 data
                if (tag_elec_itr.get_trigMatch_e24_lhmedium_iloose_L1EM20VH()) { // use get_trigMatch_e24_lhmedium_iloose_L1EM20VH() because get_trigMatch_e24_lhmedium_L1EM20VH() doesn't work correctly.
                    is_tag_trigger_matched = true;
                    candidate.tag_is_trigger_matched = true;
                }
                if (probe_elec_itr.get_trigMatch_e24_lhmedium_iloose_L1EM20VH()) {
                    is_probe_trigger_matched = true;
                    candidate.probe_is_trigger_matched = true;
                }
            }
            else if (run_number > 290000) {
                if (tag_elec_itr.get_trigMatch_e26_lhtight_nod0_ivarloose()) {
                    is_tag_trigger_matched = true;
                    candidate.tag_is_trigger_matched = true;
                }
                if (probe_elec_itr.get_trigMatch_e26_lhtight_nod0_ivarloose()) {
                    is_probe_trigger_matched = true;
                    candidate.probe_is_trigger_matched = true;
                }
            }

            TLorentzVector tlv_tag;
            TLorentzVector tlv_probe;
            tlv_tag.SetPtEtaPhiM(tag_elec_itr.get_pt(), tag_elec_itr.get_eta(), tag_elec_itr.get_phi(), tag_elec_itr.get_M());
            tlv_probe.SetPtEtaPhiM(probe_elec_itr.get_pt(), probe_elec_itr.get_eta(), probe_elec_itr.get_phi(), probe_elec_itr.get_M());
            double current_mll = (tlv_tag + tlv_probe).M();
            if (fabs(current_mll - 91.1876 * 1000.) < fabs(best_match_mll - 91.1876 * 1000.))
                best_match_mll = current_mll;
            // mll window cut
            if (current_mll > 80000. && current_mll < 100000.) {
                if (probe_elec_itr.get_baseline()) {
                    candidate.baseline_mll = current_mll;
                }
                if (probe_elec_itr.get_isSignal()) {
                    candidate.signal_mll = current_mll;
                }

                El_isZTag->at(tag_elec_itr.get_index()) = true;
                El_isZProbe->at(probe_elec_itr.get_index()) = true;
                El_Tag_index->at(tag_elec_itr.get_index()) = tag_elec_itr.get_index();
                El_Probe_index->at(probe_elec_itr.get_index()) = probe_elec_itr.get_index();

                if (is_tag_trigger_matched) {
                    El_isZProbe_TagTriggerMatched->at(probe_elec_itr.get_index()) = true;
                    El_Tag_Trigger_Matched->at(tag_elec_itr.get_index()) = true;
                }
                if (is_probe_trigger_matched) {
                    El_Probe_Trigger_Matched->at(probe_elec_itr.get_index()) = true;
                }
                N_tot_TandP_pair_electron++;
            }
        }
        El_ZTandP_mll->at(probe_elec_itr.get_index()) = best_match_mll;
        EE_pair->push_back(candidate);
    }

    // Muons
    for (auto & probe_muon_itr : vec_baseline_muon) {
        double best_match_mll = 1000000.;
        ZTandP_pair candidate;
        for (auto & tag_muon_itr : vec_signal_muon) {
            if (tag_muon_itr.get_index() == probe_muon_itr.get_index())
                continue;
            if (tag_muon_itr.get_pt() < tag_pt_threshold)
                continue;
            // Opposite Charge requirement
            if (tag_muon_itr.get_charge() == probe_muon_itr.get_charge())
                continue;

            candidate.event_number = event_number;

            candidate.tag_index = tag_muon_itr.get_index();
            candidate.tag_flavor = 11;
            candidate.tag_charge = tag_muon_itr.get_charge();
            candidate.tag_pt = tag_muon_itr.get_pt();
            candidate.tag_eta = tag_muon_itr.get_eta();
            candidate.tag_phi = tag_muon_itr.get_phi();
            candidate.tag_id_SFw = tag_muon_itr.get_SFw();
            candidate.tag_iso_SFw = tag_muon_itr.get_IsoSFw();
            candidate.tag_trigger_SFw = 0.; // no information

            candidate.probe_index = probe_muon_itr.get_index();
            candidate.probe_flavor = 11;
            candidate.probe_charge = probe_muon_itr.get_charge();
            candidate.probe_pt = probe_muon_itr.get_pt();
            candidate.probe_eta = probe_muon_itr.get_eta();
            candidate.probe_phi = probe_muon_itr.get_phi();
            candidate.probe_id_SFw = probe_muon_itr.get_SFw();
            candidate.probe_iso_SFw = probe_muon_itr.get_IsoSFw();
            candidate.probe_trigger_SFw = 0.; // no information
            candidate.probe_is_baseline = probe_muon_itr.get_baseline();
            candidate.probe_is_signal = probe_muon_itr.get_isSignal();

            // Check does the tag tragger matched
            // muons:
            // 2015: HLT_mu20_iloose_L1MU15
            // 2016: HLT_mu26_ivarmedium
            bool is_tag_trigger_matched = false;
            bool is_probe_trigger_matched = false;
            if (run_number < 290000) { // 2015 data
                if (tag_muon_itr.get_trigMatch_mu20_iloose_L1MU15()) {
                    is_tag_trigger_matched = true;
                    candidate.tag_is_trigger_matched = true;
                }
                if (probe_muon_itr.get_trigMatch_mu20_iloose_L1MU15()) {
                    is_probe_trigger_matched = true;
                    candidate.probe_is_trigger_matched = true;
                }
            }
            else if (run_number > 290000) { // 2016 data
                if (tag_muon_itr.get_trigMatch_mu26_ivarmedium()) {
                    is_tag_trigger_matched = true;
                    candidate.tag_is_trigger_matched = true;
                }
                if (probe_muon_itr.get_trigMatch_mu26_ivarmedium()) {
                    is_probe_trigger_matched = true;
                    candidate.probe_is_trigger_matched = true;
                }
            }

            TLorentzVector tlv_tag;
            TLorentzVector tlv_probe;
            tlv_tag.SetPtEtaPhiM(tag_muon_itr.get_pt(), tag_muon_itr.get_eta(), tag_muon_itr.get_phi(), tag_muon_itr.get_M());
            tlv_probe.SetPtEtaPhiM(probe_muon_itr.get_pt(), probe_muon_itr.get_eta(), probe_muon_itr.get_phi(), probe_muon_itr.get_M());
            double current_mll = (tlv_tag + tlv_probe).M();
            if (fabs(current_mll - 91.1876 * 1000.) < fabs(best_match_mll - 91.1876 * 1000.))
                best_match_mll = current_mll;
            // mll window cut
            if (current_mll > 80000. && current_mll < 100000.) {
                if (probe_muon_itr.get_baseline()) {
                    candidate.baseline_mll = current_mll;
                }
                if (probe_muon_itr.get_isSignal()) {
                    candidate.signal_mll = current_mll;
                }

                Mu_isZTag->at(tag_muon_itr.get_index()) = true;
                Mu_isZProbe->at(probe_muon_itr.get_index()) = true;
                Mu_Tag_index->at(tag_muon_itr.get_index()) = tag_muon_itr.get_index();
                Mu_Probe_index->at(probe_muon_itr.get_index()) = probe_muon_itr.get_index();

                if (is_tag_trigger_matched) {
                    Mu_isZProbe_TagTriggerMatched->at(probe_muon_itr.get_index()) = true;
                    Mu_Tag_Trigger_Matched->at(tag_muon_itr.get_index()) = true;
                }
                if (is_probe_trigger_matched) {
                    Mu_Probe_Trigger_Matched->at(probe_muon_itr.get_index()) = true;
                }
                N_tot_TandP_pair_muon++;
            }
        }
        Mu_ZTandP_mll->at(probe_muon_itr.get_index()) = best_match_mll;
        MuMu_pair->push_back(candidate);
    }
}

void yt_skim::debug_print()
{
    cout << "isMC=" << isMC << endl;
    cout << "isData=" << isData << endl;
    cout << "process=" << process << endl;
    cout << "luminosity=" << luminosity << endl;
    cout << "cross_section=" << cross_section << endl;
    cout << "k_factor=" << k_factor << endl;
    cout << "filter_efficiency=" << filter_efficiency << endl;
    cout << "derivation_stat_weights=" << event_weight_sum << endl;

    cout << "N_at_least_two_baseline_lepton_events=" << N_at_least_two_baseline_lepton_events << endl;
    cout << "N_exactly_two_baseline_lepton_events=" << N_exactly_two_baseline_lepton_events << endl;
    cout << "N_tot_baseline_electrons=" << N_tot_baseline_electrons << endl; 
    cout << "N_tot_baseline_muons=" << N_tot_baseline_muons << endl;
    cout << "N_tot_TandP_pair_electron=" << N_tot_TandP_pair_electron << endl;
    cout << "N_tot_TandP_pair_muon=" << N_tot_TandP_pair_muon << endl;
}
