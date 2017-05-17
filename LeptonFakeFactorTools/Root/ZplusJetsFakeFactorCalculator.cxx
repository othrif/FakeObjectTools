#include "LeptonFakeFactorTools/ZplusJetsFakeFactorCalculator.h"

//=============================================================================
// Constructor
//=============================================================================
ZplusJetsFakeFactorCalculator::ZplusJetsFakeFactorCalculator(const std::string& name) :
  BaseFakeFactorCalculator(name)
{
  return;
}

//=============================================================================
// Destructor
//=============================================================================
ZplusJetsFakeFactorCalculator::~ZplusJetsFakeFactorCalculator()
{
  return;
}

//=============================================================================
// Initialize
//=============================================================================
StatusCode ZplusJetsFakeFactorCalculator::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  sc = BaseFakeFactorCalculator::initialize();
  if(!sc){
    ATH_MSG_FATAL("BaseFakeFactorCalculator initialize failed!");
    return sc;
  }

  // Mostly for debugging purposes; eventually I'll either remove this or improve this
  m_cutflow = std::vector<int>(16,0);

  // Create output file and histograms
  m_outputFile->mkdir(m_subregion)->cd();
  initAllFakeFactorVsPtHists();

  if(!m_do_ttt && !m_do_ltt && !m_do_tlt && !m_do_ttl){
    ATH_MSG_FATAL("Must specify a fake factor subregion! Options are ttt, ltt, tlt, or ttl");
    sc = StatusCode::FAILURE;
    return sc;
  }

  return sc;
}

//=============================================================================
// Finalize
//=============================================================================
StatusCode ZplusJetsFakeFactorCalculator::finalize()
{
  // use an int as a StatusCode
  StatusCode sc = BaseFakeFactorCalculator::finalize();
  if(!sc){
    ATH_MSG_FATAL("BaseFakeFactorCalculator finalize failed!");
    return sc;
  }

  ATH_MSG_INFO("-------------------");
  ATH_MSG_INFO("eee | " << m_num_eee);
  ATH_MSG_INFO("eem | " << m_num_eem); 
  ATH_MSG_INFO("eme | " << m_num_eme); 
  ATH_MSG_INFO("emm | " << m_num_emm); 
  ATH_MSG_INFO("mee | " << m_num_mee); 
  ATH_MSG_INFO("mem | " << m_num_mem); 
  ATH_MSG_INFO("mme | " << m_num_mme); 
  ATH_MSG_INFO("mmm | " << m_num_mmm); 
  ATH_MSG_INFO("-------------------"); 

  for(int i = 0; i < (int)m_cutflow.size(); ++i){
    ATH_MSG_INFO("cut number " << i << ": " << m_cutflow[i]);
  }

  return sc;
}

//=============================================================================
// The main accept() function
//=============================================================================
bool ZplusJetsFakeFactorCalculator::accept(EventDef& evt){

  // Keep numerators blind for data
  if(m_do_SR && m_do_ttt && m_sample == "data") return false;

  m_cutflow[0]++;

  // Trigger requirement
  if( !(evt.passSingleLeptonTrigger || evt.passDileptonTrigger) ) return false;

  m_cutflow[1]++;

  // For Z+jets, ignore events with more than 3 baseline leptons. 
  // Note that one could use an antiID selection that is not a 
  // subset of the baseline, in which case there may be 
  // fewer than 3 baseline leptons in a valid Z+jet event
  if(evt.nBaselineLeptons > 3) return false;

  m_cutflow[2]++;

  // For nJet binning
  if(m_reqZeroJets && evt.nJets != 0) return false;
  if(m_reqGtEq1Jet && evt.nJets  < 1) return false;

  // Bjet veto
  if(m_reqZeroBJets && evt.nBJets != 0) return false;
  if(m_reqGtEq1BJet && evt.nBJets  < 1) return false;

  // For lead jet pt cut
  if(evt.nJets >= 0){
    if(evt.leadJetPt < m_ptJet1CutLower) return false;
  }

  m_cutflow[3]++;

  if(evt.met_Et > m_METCutUpper) return false;
  if(evt.met_Et < m_METCutLower) return false;

  m_cutflow[4]++;

  // Z-lead, Z-sublead, W
  // Find the signal leptons closest to the Z-peak
  // Indices correspond to the same indices as the ntuples
  int SFOS_lep1index = -1;
  int SFOS_lep2index = -1;
  int index_Wcand = -1;
  int index_fake = -1;
  int index_numobj = -1;

  // Leptons saved in our EventDef
  int nlep = evt.lepPt.size();
  std::vector<bool> lepPassWSel(nlep,false);
  std::vector<bool> lepPassZSel(nlep,false);
  std::vector<bool> lepPassAntiID(nlep,false);
  std::vector<bool> lepRemove(nlep,false);

  int n_interestingLeps = 0;
  int n_Wsignal = 0;
  int n_Zsignal = 0;
  int n_antiID = 0;

  // loop over leptons in the event
  for(int i = 0; i < nlep; ++i){
    bool passesWCriteria      = evt.lepPassSigID.at(i);
    bool passesZCriteria      = evt.lepPassOtherID.at(i);
    bool passesAntiIDCriteria = evt.lepPassAntiID.at(i);
    bool interestingLepton = false;

    // Did it pass the signal lepton criteria for leptons
    // matching to the W boson?
    if(passesWCriteria){
      interestingLepton = true;
      lepPassWSel[i] = true;
      ++n_Wsignal;
    }

    // Did it pass the signal lepton criteria for leptons
    // matching to the Z boson?
    if(passesZCriteria){
      interestingLepton = true;
      lepPassZSel[i] = true;
      ++n_Zsignal;
    }

    // Did it pass the antiID criteria?
    if(passesAntiIDCriteria){
      interestingLepton = true;
      lepPassAntiID[i] = true;
      ++n_antiID;
    }

    if(interestingLepton){
      ++n_interestingLeps; // n leptons that passed W, Z, or antiID criteria
    }
    else{
      lepRemove[i] = true; // ignore these leptons forever into the future!
    }
  }

  ATH_MSG_VERBOSE("here nlep " << nlep << " n_interestingLeps " << n_interestingLeps << " nWlep " << n_Wsignal << " nAntiID " << n_antiID << " nZlep " << n_Zsignal);

  // Targeting Z+jets, so there better be at least 2 leptons passing the Z ID criteria!
  if(n_Zsignal < 2) return false;
  m_cutflow[5]++;

  // At least one W (numerator) or one antiID (denominator) lepton
  if(n_Wsignal == 0 && n_antiID == 0) return false;
  m_cutflow[6]++;

  // Exactly three leptons passing the W, Z, antiID criteria
  if(n_interestingLeps != 3) return false;

  double MT = FLT_MAX;
  double min_delta = FLT_MAX;

  // i = leading Z lepton
  // j = subleading Z lepton
  // k = W lepton (or fake)
  for(int i = 0; i < nlep; ++i){
    for(int j = 0; j < nlep; ++j){
      if(i==j) continue;

      for(int k = 0; k < nlep; ++k){
        if(i==k) continue;
        if(j==k) continue;

        // skip the removed leptons
        if(lepRemove[i] || lepRemove[j] || lepRemove[k]) continue;

        // enforce leading / subleading order for the Z leptons
        if(evt.lepPt.at(i) < evt.lepPt.at(j)) continue;

        // SFOS Z-lepton pair
        if(m_do_Zjet || (m_do_ttbar && m_do_SR)){
          if(!evt.isSFOS(i,j)) continue;
        }

        // For ttbar selection, require i,j to be a DFOS pair, and exclude all SFOS pairs
        if(m_do_ttbar && !m_do_SR){

          // Impose different flavor, opposite sign for i,j
          if(!evt.isDFOS(i,j)) continue;

          // Exclude possibility of i,k or j,k as SFOS pairs
          if(evt.isSFOS(i,k)) continue;
          if(evt.isSFOS(j,k)) continue;
        }

        // mll and mT
        double tmp_mll = evt.getMll(i,j);
        double tmp_MT  = evt.getMT(k);

        double tmp_delta = fabs(tmp_mll - 91.188);

        // if m_pairLeptonsUsingMinMt: minimize mT for assigning the leptons
        // else: minimize mass diff w.r.t. true Z mass to assign the leptons
        // (note in either case, i and j already had to pass the SFOS cut)
        if(m_pairLeptonsUsingMinMt){
          if(tmp_MT > MT) continue;

          min_delta = tmp_delta;
          MT = tmp_MT;

          SFOS_lep1index = i;
          SFOS_lep2index = j;
          index_Wcand = k;

        }
        else{
          if(tmp_delta > min_delta) continue;

          min_delta = tmp_delta;
          MT = tmp_MT;

          SFOS_lep1index = i;
          SFOS_lep2index = j;
          index_Wcand = k;
        }
      } // k
    } // j
  } // i

  m_cutflow[7]++;

  // mll cut
  if(min_delta > m_Zwindow/2.) return false;

  // Here, we make sure our mll and m3L are outside of the Zwindow
  // Due to low mass resonances, we also require mll and m3L above 20 GeV
  if(m_ZwindowVeto > 0){
    if(min_delta < m_ZwindowVeto/2.) return false;
    if(evt.getMll(SFOS_lep1index, SFOS_lep2index) < 20) return false;

    double m3L = (evt.getLeptonTLV(SFOS_lep1index) + evt.getLeptonTLV(SFOS_lep2index) + evt.getLeptonTLV(index_Wcand)).M();
    if( fabs(m3L - 91.188) < m_ZwindowVeto/2. ) return false;
    if( m3L < 20 ) return false;
  }

  // MT cuts
  if(MT < m_MtCutLower) return false;
  if(MT > m_MtCutUpper) return false;

  // Valid indices for each lepton
  if(SFOS_lep1index < 0) return false;
  if(SFOS_lep2index < 0) return false;
  if(index_Wcand < 0) return false;

  // pT requirements
  std::vector<double> vecOfLepPt = {evt.lepPt.at(index_Wcand),evt.lepPt.at(SFOS_lep1index),evt.lepPt.at(SFOS_lep2index)};
  std::sort(vecOfLepPt.begin(),vecOfLepPt.end(),std::greater<double>()); // sort in descending order
  double ptLep1 = vecOfLepPt.at(0);
  double ptLep2 = vecOfLepPt.at(1);
  double ptLep3 = vecOfLepPt.at(2);

  if(ptLep1 < m_ptLep1CutLower) return false;
  if(ptLep2 < m_ptLep2CutLower) return false;
  if(ptLep3 < m_ptLep3CutLower) return false;

  double ptlll = (evt.getLeptonTLV(SFOS_lep1index) + evt.getLeptonTLV(SFOS_lep2index) + evt.getLeptonTLV(index_Wcand)).Pt();
  if(ptlll > m_ptlllCutUpper) return false;

  m_cutflow[8]++;

  // for bookkeeping purposes
  TString chanFlavor = "";
  if     (evt.isElectron(index_Wcand))    chanFlavor += "e";
  else if(evt.isMuon(index_Wcand))        chanFlavor += "m";
  if     (evt.isElectron(SFOS_lep1index)) chanFlavor += "e";
  else if(evt.isMuon(SFOS_lep1index))     chanFlavor += "m";
  if     (evt.isElectron(SFOS_lep2index)) chanFlavor += "e";
  else if(evt.isMuon(SFOS_lep2index))     chanFlavor += "m";

  // What is the truth information for each lepton?
  // If e.g. truth_leadz, then the leading Z lepton is a real lepton.
  bool truth_leadz    = evt.lepIsTruth.at(SFOS_lep1index);
  bool truth_subleadz = evt.lepIsTruth.at(SFOS_lep2index);
  bool truth_W        = evt.lepIsTruth.at(index_Wcand);
  int nTruthLeps = (int)truth_leadz + (int)truth_subleadz + (int)truth_W;

  int index_truthFake = -1;
  int type_fake = -1;

  if( nTruthLeps == 3 ){
    // ERROR! Three real leptons!
    index_truthFake = index_Wcand;
  }
  else if( nTruthLeps == 1 ){
    // ERROR! Only one real lepton!
    index_truthFake = index_Wcand;
  }
  else if( nTruthLeps == 0 ){
    // ERROR! No real leptons!
    index_truthFake = index_Wcand;
  }
  else{
    if(!truth_leadz){
      // Leading Z lepton is fake
      index_truthFake = SFOS_lep1index;
      type_fake = evt.lepFlavor.at(SFOS_lep1index);
    }
    else if(!truth_subleadz){
      // Subleading Z lepton is fake
      index_truthFake = SFOS_lep2index;
      type_fake = evt.lepFlavor.at(SFOS_lep2index);
    }
    else if(!truth_W){
      // W lepton is fake
      index_truthFake = index_Wcand;
      type_fake = evt.lepFlavor.at(index_Wcand);
    }
  }

  // If you only want to consider events with one fake and two real
  // leptons, in some particular configuration, uncomment one of the
  // below lines.
  //if(index_truthFake != index_Wcand) return false; // FRR
  //if(index_truthFake != SFOS_lep1index) return false; // RFR
  //if(index_truthFake != SFOS_lep2index) return false; // RRF


  // Categorize using truth information when looking at the MC.
  // Currently only used for the ttbar estimate, to apply the
  // ttbar SF in the SR, since the ttbar ttt scenario
  // is the only case where we are completely uncertain of which 
  // lepton might actually be the fake
  if(m_do_ttt && (m_sample == "ttbar" || m_sample == "tw" || m_sample == "ww")){

    // if the three reco leptons did not truth match to exactly three truth leptons
    if(type_fake < 0){

      // guaranteed to be an e or mu fake in eee or mmm
      if(chanFlavor == "eee") type_fake = 1;
      else if(chanFlavor == "mmm") type_fake = 2;
      else{

        int lowestPtIndex = index_Wcand;

        if(!evt.lepIsTruth.at(SFOS_lep1index) && evt.lepPt.at(SFOS_lep1index) < evt.lepPt.at(lowestPtIndex)) lowestPtIndex = SFOS_lep1index;
        if(!evt.lepIsTruth.at(SFOS_lep2index) && evt.lepPt.at(SFOS_lep2index) < evt.lepPt.at(lowestPtIndex)) lowestPtIndex = SFOS_lep2index;


        type_fake = evt.lepFlavor.at(lowestPtIndex);

        ATH_MSG_WARNING("Ambiguity in which lepton is fake! Taking the flavor of the lowest pt fake lepton. " 
            << " pt: "
            << evt.lepPt.at(SFOS_lep1index) << " "
            << evt.lepPt.at(SFOS_lep2index) << " "
            << evt.lepPt.at(index_Wcand) << " "
            << " truth: "
            << evt.lepIsTruth.at(SFOS_lep1index) << " "
            << evt.lepIsTruth.at(SFOS_lep2index) << " "
            << evt.lepIsTruth.at(index_Wcand) << " "
            << " flavor: "
            << evt.lepFlavor.at(SFOS_lep1index) << " "
            << evt.lepFlavor.at(SFOS_lep2index) << " "
            << evt.lepFlavor.at(index_Wcand) << " "
            << " charge: "
            << evt.lepCharge.at(SFOS_lep1index) << " "
            << evt.lepCharge.at(SFOS_lep2index) << " "
            << evt.lepCharge.at(index_Wcand) << " "
            << " weight: "
            << evt.weight << " "
            << " lowest pt flavor "
            << type_fake
        );
      }   
    }

    // Be careful that you've set the SFs for the ttbar SR appropriately!!!
    // Also note that the order of multiplying SF_err*weight and SF*weight matters!
    if(type_fake == 1){
      evt.weight_syst_err = addInQuadrature(evt.weight_syst_err, m_ttbar_SF_err_el*evt.weight);
      evt.weight *= m_ttbar_SF_el;
    }
    if(type_fake == 2){
      evt.weight_syst_err = addInQuadrature(evt.weight_syst_err, m_ttbar_SF_err_mu*evt.weight);
      evt.weight *= m_ttbar_SF_mu;
    }
  }

  m_cutflow[9]++;

  // Now we can check the lepton quality, after assigning them to the W or Z.
  // Start with trig matching:
  bool trigMatch_leadZ = evt.lepMatchesTrigger.at(SFOS_lep1index);
  bool trigMatch_subleadZ = evt.lepMatchesTrigger.at(SFOS_lep2index);
  bool trigMatch_W = evt.lepMatchesTrigger.at(index_Wcand);

  // is it a numerator event?
  if(m_do_ttt){

    // Check if any of the leptons fired the appropriate trigger.
    // Skip W candidate for the fake factor region's trigger requirement,
    // since that is done for the denominator as well.
    if(!m_do_SR && m_do_Zjet){
      if(evt.passSingleLeptonTrigger){
        if(!trigMatch_leadZ && !trigMatch_subleadZ) return false;
      }
      else if(evt.passDileptonTrigger){
        if( !(trigMatch_leadZ && trigMatch_subleadZ) ) return false;
      }
    }
    else{
      if(evt.passSingleLeptonTrigger){
        if(!trigMatch_leadZ && !trigMatch_subleadZ && !trigMatch_W) return false;
      }
      else if(evt.passDileptonTrigger){
        if( !((trigMatch_leadZ && trigMatch_subleadZ) || (trigMatch_leadZ && trigMatch_W) || (trigMatch_subleadZ && trigMatch_W)) ) return false;
      }
    }

    // Require appropriate criteria for each lep
    if(!lepPassWSel[index_Wcand]) return false;
    if(!lepPassZSel[SFOS_lep1index]) return false;
    if(!lepPassZSel[SFOS_lep2index]) return false;

    // numerator object!
    index_numobj = index_Wcand;
  }
  else if(m_do_ltt){

    // Trig match the tight leptons
    if(evt.passSingleLeptonTrigger){
      if(!trigMatch_leadZ && !trigMatch_subleadZ) return false;
    }
    else if(evt.passDileptonTrigger){
      if( !(trigMatch_leadZ && trigMatch_subleadZ) ) return false;
    }

    // Require appropriate criteria for each lep
    if(!lepPassAntiID[index_Wcand]) return false;
    if(!lepPassZSel[SFOS_lep1index]) return false;
    if(!lepPassZSel[SFOS_lep2index]) return false;

    // fake object!
    index_fake = index_Wcand;
  }
  else if(m_do_tlt){

    // Trig match the tight leptons
    if(evt.passSingleLeptonTrigger){
      if(!trigMatch_W && !trigMatch_subleadZ) return false;
    }
    else if(evt.passDileptonTrigger){
      if( !(trigMatch_W && trigMatch_subleadZ) ) return false;
    }

    // Require appropriate criteria for each lep
    if(!lepPassWSel[index_Wcand]) return false;
    if(!lepPassAntiID[SFOS_lep1index]) return false;
    if(!lepPassZSel[SFOS_lep2index]) return false;

    // fake object!
    index_fake = SFOS_lep1index;

    // in some cases, require an additional SFOS pair consistent with the Z,
    // since tlt or ttl events with euu or uee are unlikely to be Z+jet
    if(m_do_Zjet && !m_includeUnlikelyZjetEvents){
      if(!evt.isSFOS(SFOS_lep2index,index_Wcand)) return false;
    }
  }
  else if(m_do_ttl){

    // Trig match the tight leptons
    if(evt.passSingleLeptonTrigger){
      if(!trigMatch_W && !trigMatch_leadZ) return false;
    }
    else if(evt.passDileptonTrigger){
      if( !(trigMatch_W && trigMatch_leadZ) ) return false;
    }

    // Require appropriate criteria for each lep
    if(!lepPassWSel[index_Wcand]) return false;
    if(!lepPassZSel[SFOS_lep1index]) return false;
    if(!lepPassAntiID[SFOS_lep2index]) return false;

    // fake object!
    index_fake = SFOS_lep2index;

    // in some cases, require an additional SFOS pair consistent with the Z,
    // since tlt or ttl events with euu or uee are unlikely to be Z+jet
    if(m_do_Zjet && !m_includeUnlikelyZjetEvents){
      if(!evt.isSFOS(SFOS_lep1index,index_Wcand)) return false;
    }
  }

  // Apply the ttbar data-to-mc ZCR SFs derived previously, but only to the denominator
  if(!m_do_ttt && (m_sample == "ttbar" || m_sample == "tw" || m_sample == "ww") ){
    if(evt.isElectron(index_fake)){
      evt.weight_syst_err = addInQuadrature(evt.weight_syst_err, m_ttbar_SF_err_el*evt.weight);
      evt.weight *= m_ttbar_SF_el;
    }
    else if(evt.isMuon(index_fake)){
      evt.weight_syst_err = addInQuadrature(evt.weight_syst_err, m_ttbar_SF_err_mu*evt.weight);
      evt.weight *= m_ttbar_SF_mu;
    }
  }

  m_cutflow[10]++;

  m_cutflow[11]++;

  m_cutflow[12]++;
  m_cutflow[13]++;

  m_cutflow[14]++;


  // FILL IN FAKE PLOTS
  if(index_fake >= 0){

    double pt = evt.lepPt.at(index_fake);
    bool isElectron = evt.isElectron(index_fake);

    // Denominator plots
    bool isNumPlot = false; // i.e. isDenPlot
    fillAllFakeFactorVsPtHists(pt, chanFlavor, isElectron, isNumPlot, evt.weight, evt.weight_syst_err);

    // Numerator via ff*denominator
    if(m_doApplyFakeFactor){

      // Modify weight by ff, but retain nominal weight for later plots
      double ff = 1;
      double ffStatErr = 0;
      double ffSystErr = 0;
      LepEnum::LepType typeOfLep = LepEnum::Electron;
      if(!isElectron)  typeOfLep = LepEnum::Muon;

      if(m_do_ltt){
        ff = m_applyFakeFactorToolSignalLep->apply(pt, typeOfLep);
        ffStatErr = m_applyFakeFactorToolSignalLep->getStatErr(pt, typeOfLep);
        ffSystErr = m_applyFakeFactorToolSignalLep->getSystErr(pt, typeOfLep);
      }
      else if(m_do_tlt || m_do_ttl){
        ff = m_applyFakeFactorToolOtherLep->apply(pt, typeOfLep);
        ffStatErr = m_applyFakeFactorToolOtherLep->getStatErr(pt, typeOfLep);
        ffSystErr = m_applyFakeFactorToolOtherLep->getSystErr(pt, typeOfLep);
      }

      ATH_MSG_VERBOSE(typeOfLep << " ff is " << ff);

      // Scale the weight by the FF
      double ffWeight = ff*evt.weight;

      // Note here that evt.weight_syst_err is effectively just the ttbar SF error * evt.weight,
      // so the total error is the sum in quadrature of the properly weighted ttbar error,
      // the FF stat error, and the FF syst error
      double ffStatErrWeight = ffStatErr*evt.weight;
      double ffSystErrWeight = ffSystErr*evt.weight;
      double ffTotErrWeight = addInQuadrature(ff*evt.weight_syst_err, ffStatErrWeight, ffSystErrWeight);

      // Now point the weight and weigt_syst_err at these values
      evt.weight = ffWeight;
      evt.weight_syst_err = ffTotErrWeight;

      ATH_MSG_VERBOSE("ffWeight is " << ffWeight << " ffStatErrWeight is " << ffStatErrWeight << " ffSystErrWeight is " << ffSystErrWeight << " ffTotErrWeight is " << ffTotErrWeight);

      isNumPlot = true;
      fillAllFakeFactorVsPtHists(pt, chanFlavor, isElectron, isNumPlot, ffWeight, ffTotErrWeight);
    }
  }

  // Numerator plots
  if(index_numobj >= 0 && !m_doApplyFakeFactor){

    double pt = evt.lepPt.at(index_numobj);
    bool isElectron = evt.isElectron(index_numobj);
    bool isNumPlot = true;

    fillAllFakeFactorVsPtHists(pt, chanFlavor, isElectron, isNumPlot, evt.weight, evt.weight_syst_err);
  }

  if     (chanFlavor == "eee") ++m_num_eee;
  else if(chanFlavor == "eem") ++m_num_eem;
  else if(chanFlavor == "eme") ++m_num_eme;
  else if(chanFlavor == "emm") ++m_num_emm;
  else if(chanFlavor == "mee") ++m_num_mee;
  else if(chanFlavor == "mem") ++m_num_mem;
  else if(chanFlavor == "mme") ++m_num_mme;
  else if(chanFlavor == "mmm") ++m_num_mmm;

  ATH_MSG_VERBOSE("Hmmm.... " 
      << " pt: "
      << evt.lepPt.at(index_Wcand) << " "
      << evt.lepPt.at(SFOS_lep1index) << " "
      << evt.lepPt.at(SFOS_lep2index) << " "
      << " truth: "
      << evt.lepIsTruth.at(index_Wcand) << " "
      << evt.lepIsTruth.at(SFOS_lep1index) << " "
      << evt.lepIsTruth.at(SFOS_lep2index) << " "
      << " flavor: "
      << evt.lepFlavor.at(index_Wcand) << " "
      << evt.lepFlavor.at(SFOS_lep1index) << " "
      << evt.lepFlavor.at(SFOS_lep2index) << " "
      << " charge: "
      << evt.lepCharge.at(index_Wcand) << " "
      << evt.lepCharge.at(SFOS_lep1index) << " "
      << evt.lepCharge.at(SFOS_lep2index) << " "
      << " weight: "
      << evt.weight << " "
      << " lowest pt flavor "
      << type_fake
  );

  evt.lep0_index = index_Wcand;
  evt.lep1_index = SFOS_lep1index;
  evt.lep2_index = SFOS_lep2index;

  m_cutflow[15]++;

  ATH_MSG_VERBOSE("made it " << chanFlavor);
  return true;
}



//=============================================================================
// Helper functions
//=============================================================================

void ZplusJetsFakeFactorCalculator::initAllFakeFactorVsPtHists(void){
  m_h_allNum_pt = initFakeFactorVsPtHist(m_h_allNum_pt, "allNum_"+m_sample+"_pt");
  m_h_allDen_pt = initFakeFactorVsPtHist(m_h_allDen_pt, "allDen_"+m_sample+"_pt");
  m_h_allNum_pt_syst = initFakeFactorVsPtHist(m_h_allNum_pt_syst, "allNum_"+m_sample+"_pt_syst");
  m_h_allDen_pt_syst = initFakeFactorVsPtHist(m_h_allDen_pt_syst, "allDen_"+m_sample+"_pt_syst");

  if(m_do_Zjet || m_do_ttbar){
    m_h_elNum_pt = initFakeFactorVsPtHist(m_h_elNum_pt, "elNum_"+m_sample+"_pt");
    m_h_elDen_pt = initFakeFactorVsPtHist(m_h_elDen_pt, "elDen_"+m_sample+"_pt");
    m_h_muNum_pt = initFakeFactorVsPtHist(m_h_muNum_pt, "muNum_"+m_sample+"_pt");
    m_h_muDen_pt = initFakeFactorVsPtHist(m_h_muDen_pt, "muDen_"+m_sample+"_pt");

    m_h_eeeNum_pt = initFakeFactorVsPtHist(m_h_eeeNum_pt, "eeeNum_"+m_sample+"_pt");
    m_h_eeeDen_pt = initFakeFactorVsPtHist(m_h_eeeDen_pt, "eeeDen_"+m_sample+"_pt");
    m_h_emmNum_pt = initFakeFactorVsPtHist(m_h_emmNum_pt, "emmNum_"+m_sample+"_pt");
    m_h_emmDen_pt = initFakeFactorVsPtHist(m_h_emmDen_pt, "emmDen_"+m_sample+"_pt");
    m_h_meeNum_pt = initFakeFactorVsPtHist(m_h_meeNum_pt, "meeNum_"+m_sample+"_pt");
    m_h_meeDen_pt = initFakeFactorVsPtHist(m_h_meeDen_pt, "meeDen_"+m_sample+"_pt");
    m_h_mmmNum_pt = initFakeFactorVsPtHist(m_h_mmmNum_pt, "mmmNum_"+m_sample+"_pt");
    m_h_mmmDen_pt = initFakeFactorVsPtHist(m_h_mmmDen_pt, "mmmDen_"+m_sample+"_pt");

    m_h_elNum_pt_syst = initFakeFactorVsPtHist(m_h_elNum_pt_syst, "elNum_"+m_sample+"_pt_syst");
    m_h_elDen_pt_syst = initFakeFactorVsPtHist(m_h_elDen_pt_syst, "elDen_"+m_sample+"_pt_syst");
    m_h_muNum_pt_syst = initFakeFactorVsPtHist(m_h_muNum_pt_syst, "muNum_"+m_sample+"_pt_syst");
    m_h_muDen_pt_syst = initFakeFactorVsPtHist(m_h_muDen_pt_syst, "muDen_"+m_sample+"_pt_syst");

    m_h_eeeNum_pt_syst = initFakeFactorVsPtHist(m_h_eeeNum_pt_syst, "eeeNum_"+m_sample+"_pt_syst");
    m_h_eeeDen_pt_syst = initFakeFactorVsPtHist(m_h_eeeDen_pt_syst, "eeeDen_"+m_sample+"_pt_syst");
    m_h_emmNum_pt_syst = initFakeFactorVsPtHist(m_h_emmNum_pt_syst, "emmNum_"+m_sample+"_pt_syst");
    m_h_emmDen_pt_syst = initFakeFactorVsPtHist(m_h_emmDen_pt_syst, "emmDen_"+m_sample+"_pt_syst");
    m_h_meeNum_pt_syst = initFakeFactorVsPtHist(m_h_meeNum_pt_syst, "meeNum_"+m_sample+"_pt_syst");
    m_h_meeDen_pt_syst = initFakeFactorVsPtHist(m_h_meeDen_pt_syst, "meeDen_"+m_sample+"_pt_syst");
    m_h_mmmNum_pt_syst = initFakeFactorVsPtHist(m_h_mmmNum_pt_syst, "mmmNum_"+m_sample+"_pt_syst");
    m_h_mmmDen_pt_syst = initFakeFactorVsPtHist(m_h_mmmDen_pt_syst, "mmmDen_"+m_sample+"_pt_syst");
  }
}

void ZplusJetsFakeFactorCalculator::fillAllFakeFactorVsPtHists(double pt, TString chanFlavor, bool isElectron, bool isNumPlot, double weight, double weight_syst_err){

  // Fill plots for:
  // 1) combined num (or den)
  // 2) electron OR muon num (or den)
  // 3) eee channel num (or den)
  // OR emm/mem/mme channel num (or den)
  // OR mee/eme/eem channel num (or den)
  // OR mmm channel num (or den)

  // point to the num hists to start
  TH1F* h_all_pt = m_h_allNum_pt;
  TH1F* h_el_pt  = m_h_elNum_pt;
  TH1F* h_mu_pt  = m_h_muNum_pt;
  TH1F* h_eee_pt = m_h_eeeNum_pt;
  TH1F* h_emm_pt = m_h_emmNum_pt;
  TH1F* h_mee_pt = m_h_meeNum_pt;
  TH1F* h_mmm_pt = m_h_mmmNum_pt;

  TH1F* h_all_pt_syst = m_h_allNum_pt_syst;
  TH1F* h_el_pt_syst  = m_h_elNum_pt_syst;
  TH1F* h_mu_pt_syst  = m_h_muNum_pt_syst;
  TH1F* h_eee_pt_syst = m_h_eeeNum_pt_syst;
  TH1F* h_emm_pt_syst = m_h_emmNum_pt_syst;
  TH1F* h_mee_pt_syst = m_h_meeNum_pt_syst;
  TH1F* h_mmm_pt_syst = m_h_mmmNum_pt_syst;

  // point instead to den hists if appropriate
  if(!isNumPlot){
    h_all_pt = m_h_allDen_pt;
    h_el_pt  = m_h_elDen_pt;
    h_mu_pt  = m_h_muDen_pt;
    h_eee_pt = m_h_eeeDen_pt;
    h_emm_pt = m_h_emmDen_pt;
    h_mee_pt = m_h_meeDen_pt;
    h_mmm_pt = m_h_mmmDen_pt;

    h_all_pt_syst = m_h_allDen_pt_syst;
    h_el_pt_syst  = m_h_elDen_pt_syst;
    h_mu_pt_syst  = m_h_muDen_pt_syst;
    h_eee_pt_syst = m_h_eeeDen_pt_syst;
    h_emm_pt_syst = m_h_emmDen_pt_syst;
    h_mee_pt_syst = m_h_meeDen_pt_syst;
    h_mmm_pt_syst = m_h_mmmDen_pt_syst;
  }

  fillHistAndAddErrInQuadrature(pt, *h_all_pt, *h_all_pt_syst, weight, weight_syst_err);

  if(isElectron){
    fillHistAndAddErrInQuadrature(pt, *h_el_pt, *h_el_pt_syst, weight, weight_syst_err);
  }
  else{
    fillHistAndAddErrInQuadrature(pt, *h_mu_pt, *h_mu_pt_syst, weight, weight_syst_err);
  }

  if(chanFlavor == "eee"){
    fillHistAndAddErrInQuadrature(pt, *h_eee_pt, *h_eee_pt_syst, weight, weight_syst_err);
  }
  else if(chanFlavor == "emm" || chanFlavor == "mem" || chanFlavor == "mme"){
    fillHistAndAddErrInQuadrature(pt, *h_emm_pt, *h_emm_pt_syst, weight, weight_syst_err);
  }
  else if(chanFlavor == "mee" || chanFlavor == "eme" || chanFlavor == "eem"){
    fillHistAndAddErrInQuadrature(pt, *h_mee_pt, *h_mee_pt_syst, weight, weight_syst_err);
  }
  else if(chanFlavor == "mmm"){
    fillHistAndAddErrInQuadrature(pt, *h_mmm_pt, *h_mmm_pt_syst, weight, weight_syst_err);
  }

  return;
}

void ZplusJetsFakeFactorCalculator::fillHistAndAddErrInQuadrature(double pt, TH1F& hist, TH1F& hist_syst, double weight, double weight_syst_err){

  // No matter what, we want to fill the histogram
  hist.Fill(pt, weight);

  /*
  // Now we add in the weight from the FF uncertainty, whenever applicable
  int bin = hist.FindBin(pt);
  double currentBinErr = hist.GetBinError(bin);
  double totBinErr = addInQuadrature(currentBinErr, weight_syst_err);
  hist.SetBinError(bin, totBinErr);
  */

  // Then we can use JUST the err from this hist to assess syst separately from stat
  hist_syst.Fill(pt, weight_syst_err);

  return;
}

double ZplusJetsFakeFactorCalculator::addInQuadrature(double val1, double val2){
  return sqrt( (val1*val1) + (val2*val2) );
}

double ZplusJetsFakeFactorCalculator::addInQuadrature(double val1, double val2, double val3){
  return addInQuadrature( addInQuadrature(val1, val2), val3 );
}
