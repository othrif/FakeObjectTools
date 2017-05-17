#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <FakesAnalysis/MyxAODAnalysis.h>
#include <FakesAnalysis/parametric_histos.h>
#include "EventLoop/OutputStream.h"
#include "FourMomUtils/xAODP4Helpers.h"
#include <TTreeFormula.h>
#include <memory.h>
#include <TSystem.h>
#include <TRandom3.h>

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"
#include "xAODBTagging/BTagging.h"
#include "JetMomentTools/JetVertexTaggerTool.h"

#include <AsgTools/MessageCheck.h>
#include "SUSYTools/SUSYObjDef_xAOD.h"

// Amg include
#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "xAODEgamma/EgammaxAODHelpers.h"

//Trigger Include
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"

//GRL include
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#endif // ROOTCORE

using namespace xAOD;

static SG::AuxElement::Accessor<int>  acc_Origin("truthOrigin");
static SG::AuxElement::Accessor<int>  acc_Type("truthType");

// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

MyxAODAnalysis :: MyxAODAnalysis (): my_XsecDB(0), objTool(0), m_Pileup(0), m_grl(0)
{
  trig_list = {
    "HLT_e24_lhmedium_L1EM20VH",
    "HLT_e26_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium",
    "HLT_e120_lhloose", 
    "HLT_mu20_iloose_L1MU15",
    "HLT_mu24_ivarloose_L1MU15",
    "HLT_mu60_0eta105_msonly"
    "HLT_e24_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium_nod0",
    "HLT_e60_medium",
    "HLT_e140_lhloose_nod0",
    "HLT_e300_etcut",
    "HLT_mu24_iloose",
    "HLT_mu24_iloose_L1MU15", 
    "HLT_mu24_ivarloose",
    "HLT_mu24_ivarloose_L1MU15", 
    "HLT_mu40",
    "HLT_mu50",
    "HLT_mu24_ivarmedium",
    "HLT_mu24_imedium",
    "HLT_mu24_ivarmedium/imedium",
    "HLT_mu26_ivarmedium/imedium",
    "HLT_mu26_ivarmedium",
  };

  trig_list_el = {
    "HLT_e24_lhmedium_L1EM20VH",
    "HLT_e26_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium",
    "HLT_e120_lhloose", 
    "HLT_e24_lhtight_nod0_ivarloose",
    "HLT_e60_lhmedium_nod0",
    "HLT_e60_medium",
    "HLT_e140_lhloose_nod0",
    "HLT_e300_etcut",
  };

  trig_list_mu = {
    "HLT_mu20_iloose_L1MU15",
    "HLT_mu24_ivarloose_L1MU15",
    "HLT_mu60_0eta105_msonly"
    "HLT_mu24_iloose",
    "HLT_mu24_iloose_L1MU15", 
    "HLT_mu24_ivarloose",
    "HLT_mu24_ivarloose_L1MU15", 
    "HLT_mu40",
    "HLT_mu50",
    "HLT_mu24_ivarmedium",
    "HLT_mu24_imedium",
    "HLT_mu24_ivarmedium/imedium",
    "HLT_mu26_ivarmedium/imedium",
    "HLT_mu26_ivarmedium",
  };

  
  


}

EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job &job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD ();
  
  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file
  
  // tell EventLoop about our output:
  // EL::OutputStream out ("output");
  //job.outputAdd (out);
  
  return EL::StatusCode::SUCCESS;  


}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

  //here the parametrisation binning is decided

  //pt parametrisation
  std::vector<double> xbins_el = {7.,12.,20.,35.,50.,1000};
  std::vector<double> xbins_mu = {5.,12.,20.,35.,50.,1000};
 

  //eta parametrisation
  std::vector<double> ybins_el = {0.,0.7,1.37,1.52,2.,2.5};
  std::vector<double> ybins_mu = {0.,0.7,1.37,1.52,2.,2.5};

  std::vector<std::vector<double> > bins_el;
  std::vector<std::vector<double> > bins_mu;

  bins_el.push_back(xbins_el);
  bins_el.push_back(ybins_el);
  bins_mu.push_back(xbins_mu);
  bins_mu.push_back(ybins_mu);

  TString name1("Loose_el");
  TString name2("Tight_el");
  TString name3("Loose_mu");
  TString name4("Tight_mu");

  //creation of the histograms
  histo_loose_el = new parametric_histos(bins_el,name1);
  histo_tight_el = new parametric_histos(bins_el,name2);
  histo_loose_mu = new parametric_histos(bins_mu,name3);
  histo_tight_mu = new parametric_histos(bins_mu,name4);

  h_class = new TH1F("h_class","h_class",6,0.5,6.5);
  wk()->addOutput(h_class);

  h_classVSpt_ele = new TH2F("h_classVSpt_ele","h_classVSpt_ele",3,0.5,3.5,3,&(xbins_el[0]));
  wk()->addOutput(h_classVSpt_ele);
  h_classVSpt_mu = new TH2F("h_classVSpt_mu","h_classVSpt_mu",3,0.5,3.5,3,&(xbins_mu[0]));
  wk()->addOutput(h_classVSpt_mu);




  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{

    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.

    const char *APP_NAME = "changeInput()";

    m_event = wk()->xaodEvent();
    Info( APP_NAME, "Number of events in this file = %lli", m_event->getEntries() ); // print long long int

    //Read the CutBookkeeper container
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if (!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()) {
        Error( APP_NAME, "Failed to retrieve CutBookkeepers from MetaData! Exiting.");
    }

    const xAOD::CutBookkeeper* allEventsCBK = 0;
    int maxcycle = -1;

    //let's find the right CBK (latest with StreamAOD input before derivations)
    for ( auto cbk : *completeCBC ) {
      if ( cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD" && cbk->cycle() > maxcycle) {
	maxcycle = cbk->cycle();
	allEventsCBK = cbk;
      }
    }

    nEventsProcessed  = 0;
    sumOfWeights  = 0.;
    sumOfWeightsSquared = 0.;

    if (allEventsCBK) {
      nEventsProcessed  = allEventsCBK->nAcceptedEvents();
      sumOfWeights        = allEventsCBK->sumOfEventWeights();
      sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
      Info( APP_NAME, "CutBookkeepers Accepted %lu SumWei %f sumWei2 %f ", nEventsProcessed, sumOfWeights, sumOfWeightsSquared);
    } else {
      Info( APP_NAME, "No relevent CutBookKeepers found" );
    }



    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.
    const char *APP_NAME = "MyxAODAnalysis::initialize()";

    bool isAtlfast = false;
    bool isData = false;
    if (wk()->metaData()->castDouble("isatlfast") == 1) isAtlfast = true;
    if (wk()->metaData()->castDouble("isdata") == 1) isData = true;

    m_event = wk()->xaodEvent();

    ST::ISUSYObjDef_xAODTool::DataSource datasource = (isData ? ST::ISUSYObjDef_xAODTool::Data : (isAtlfast ? ST::ISUSYObjDef_xAODTool::AtlfastII : ST::ISUSYObjDef_xAODTool::FullSim));

    objTool = new ST::SUSYObjDef_xAOD( "SUSYObjDef_xAOD" );
    objTool->msg().setLevel( MSG::FATAL);
   
    //set SUSYTools config file
    std::string MyConfigFile = "FakesAnalysis/SUSYTools_fakes.conf";
    

    // Configure the SUSYObjDef instance
    ANA_CHECK( objTool->setProperty("ConfigFile", MyConfigFile) );
    ANA_CHECK(objTool->setProperty("DataSource", datasource) ) ;
   
    // PU reweight tool
    std::vector<std::string> prw_conf;
    
    prw_conf.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root");
 
    std::vector<std::string> prw_lumicalc;
    prw_lumicalc.push_back("FakesAnalysis/ilumicalc_histograms_None_276262-284484.root");
    prw_lumicalc.push_back("FakesAnalysis/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-005.root");

    m_Pileup  = new CP::PileupReweightingTool("Pileup");
    m_Pileup->msg().setLevel( MSG::FATAL );

    Info( APP_NAME, "Initializing PRW");

    ANA_CHECK( m_Pileup->setProperty( "ConfigFiles", prw_conf) );
    ANA_CHECK( m_Pileup->setProperty( "LumiCalcFiles", prw_lumicalc) );
    ANA_CHECK( m_Pileup->setProperty( "DataScaleFactor", 1.0 / 1.09) );
    ANA_CHECK( m_Pileup->setProperty( "DataScaleFactorUP", 1.0 / 1.0) );
    ANA_CHECK( m_Pileup->setProperty( "DataScaleFactorDOWN", 1.0 / 1.18) );
    ANA_CHECK( m_Pileup->initialize());

    ANA_CHECK( objTool->setProperty("PRWConfigFiles", prw_conf) );
    ANA_CHECK( objTool->setProperty("PRWLumiCalcFiles", prw_lumicalc) );

  
    //SUSYTools initialization
    if (m_debug) ANA_CHECK( objTool->setProperty("OutputLevel", MSG::VERBOSE));
    if ( objTool->initialize() != StatusCode::SUCCESS) {
        Error(APP_NAME, "Cannot intialize SUSYObjDef_xAOD..." );
        Error(APP_NAME, "Exiting... " );
        return EL::StatusCode::FAILURE;
    } else {
      Info(APP_NAME, "SUSYObjDef_xAOD initialized... " );
    }
    if (m_debug) Info(APP_NAME, "SUSYTools initialized" );

    //read xsec from database
    my_XsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));
    if (m_debug) Info(APP_NAME, "xsec DB initialized" );

    //GRL initialization
    Info( APP_NAME, "Initializing GRL");
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/FakesAnalysis/data15_13TeV.periodAllYear_DetStatus-v79-repro20-02_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml"));
    vecStringGRL.push_back(gSystem->ExpandPathName("$ROOTCOREBIN/data/FakesAnalysis/data16_13TeV.periodAllYear_DetStatus-v88-pro20-21_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml"));
    ANA_CHECK( m_grl->setProperty( "GoodRunsListVec", vecStringGRL) );
    ANA_CHECK( m_grl->setProperty( "PassThrough", false) ); // if true (default) will ignore result of GRL and will just pass all events
    if (!m_grl->initialize().isSuccess()) { // check this isSuccess
      Error(APP_NAME, "Failed to properly initialize the GRL. Exiting." );
      return EL::StatusCode::FAILURE;
    }


    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: execute ()
{

  
    const char *APP_NAME = "execute()";

    if (m_debug) Info(APP_NAME, "New event" );
   

     
    //----------------------------
    // Event information
    //---------------------------
    if (m_debug) Info(APP_NAME, "Retrieving event info collection." );

   

    const xAOD::EventInfo *eventInfo = 0;
    if ( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ) {
      Error(APP_NAME, "Failed to retrieve event info collection. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    // check if the event is data or MC
    isMC = eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION );
    ANA_CHECK( objTool->ApplyPRWTool());
    whichYear = objTool->treatAsYear();

    // MC weights
    if (isMC) {
     m_xsec = my_XsecDB->xsectTimesEff(eventInfo->mcChannelNumber());
     float sherpa22_wei;
     if (eventInfo->mcChannelNumber() >= 363102 && eventInfo->mcChannelNumber() <= 363435) {
       sherpa22_wei = objTool->getSherpaVjetsNjetsWeight();
     }
     else sherpa22_wei = 1;
     m_WeightEvents = (eventInfo->mcEventWeights()).at(0) * sherpa22_wei;
   
     ANA_CHECK( m_Pileup->apply(*eventInfo, false) );
     m_PUwei= eventInfo->auxdata< float >( "PileupWeight" );

    }
    else {
      m_xsec = 1.;
      m_WeightEvents = 1;
      m_PUwei = 1.;
    }

    // apply GRL on data
    if (!isMC && !m_grl->passRunLB(*eventInfo)) {
      return EL::StatusCode::SUCCESS; // go to next event
    }

    const xAOD::VertexContainer* vertices = 0;
    if ( m_event->retrieve( vertices, "PrimaryVertices").isSuccess() ) {
      m_nPrimVx = vertices->size();
    }
    else {
        Error(APP_NAME, "Failed to retrieve Vertex container. Exiting." );
        return EL::StatusCode::FAILURE;
    }


    //------------
    // CLEANING CUTS HERE
    //------------
 

    if (!isMC &&
            ((eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error ) ||
             (eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error ) ||
             (eventInfo->errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error) ||
             (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ))) {
        return EL::StatusCode::SUCCESS; // go to next event
    }


    //Get the event Vertex
    if (m_debug) Info(APP_NAME, "Primary Vertex");
    const xAOD::Vertex* PrimVx = 0;
    PrimVx = objTool->GetPrimVtx();
    if (!PrimVx) return EL::StatusCode::SUCCESS; 


    //Filling Trigger decision
    bool at_least_one_trig = false;
    for (UInt_t it = 0; it < trig_list.size(); it++) {
        bool passed = objTool->IsTrigPassed(trig_list.at(it));
        trig_pass.push_back(passed);
        if (passed) at_least_one_trig = true;
        eventInfo->auxdecor< std::vector<bool> >("TrigPass") = trig_pass;
    }
    //If no trigger has fired, don't bother going ahead
    if (!at_least_one_trig) {
        return EL::StatusCode::SUCCESS; // go to next event
    }

    //------------
    // CALIBRATION AND OR
    //------------

    // Get the nominal object containers from the event
    // Electrons
    if (m_debug) Info(APP_NAME, "Electrons");
    xAOD::ElectronContainer* electrons_nominal(0);
    xAOD::ShallowAuxContainer* electrons_nominal_aux(0);
    ANA_CHECK( objTool->GetElectrons(electrons_nominal, electrons_nominal_aux) );

    // Muons
    if (m_debug) Info(APP_NAME, "Muons");
    xAOD::MuonContainer* muons_nominal(0);
    xAOD::ShallowAuxContainer* muons_nominal_aux(0);
    ANA_CHECK( objTool->GetMuons(muons_nominal, muons_nominal_aux) );

    // Jets
    if (m_debug) Info(APP_NAME, "Jets");
    xAOD::JetContainer* jets_nominal(0);
    xAOD::ShallowAuxContainer* jets_nominal_aux(0);
    ANA_CHECK(objTool->GetJets(jets_nominal, jets_nominal_aux));

    // Photons (this is here just for the OR)
    if (m_debug) Info(APP_NAME, "Photons");
    xAOD::PhotonContainer* photons_nominal(0);
    xAOD::ShallowAuxContainer* photons_nominal_aux(0); 
    ANA_CHECK( objTool->GetPhotons(photons_nominal, photons_nominal_aux) );

  
    if (m_debug) Info(APP_NAME, "Overlap Removal");
    ANA_CHECK( objTool->OverlapRemoval(electrons_nominal, muons_nominal, jets_nominal, photons_nominal) );


    myAnalysisCollections myCollections;
    myCollections._electrons = electrons_nominal;
    myCollections._photons = photons_nominal;
    myCollections._muons = muons_nominal;
    myCollections._jets = jets_nominal;


    //--------
    //SCALE FACTORS (trigger SF needs to be calculated in FindProbe(), as you need to run on containers!)
    //--------

    if (isMC) {
      m_eleSF_weight = objTool->GetTotalElectronSF(*electrons_nominal, true, true, false, true);
      m_muSF_weight = objTool->GetTotalMuonSF(*muons_nominal, true, true, "");
      m_bTag_weight = objTool->BtagSF(jets_nominal);
      m_JVT_weight = objTool->JVT_SF(jets_nominal);
    }
    else {
      m_eleSF_weight = 1.;
      m_muSF_weight = 1.;
      m_bTag_weight = 1.;
      m_JVT_weight = 1.;
    }

    

    //--------
    //SELECTION
    //--------


    FindProbe(myCollections);
    if(doClassification) Fill_classification();
    Fill_histos();
    

    //------
    //FINAL CLEANING
    //-----
   

 
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
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


  histo_loose_el->Execute();
  Register_parametric_histos(histo_loose_el);
  histo_tight_el->Execute();
  Register_parametric_histos(histo_tight_el);
  histo_loose_mu->Execute();
  Register_parametric_histos(histo_loose_mu);
  histo_tight_mu->Execute();
  Register_parametric_histos(histo_tight_mu);

    const char *APP_NAME = "finalize()";
    Info(APP_NAME, "All done. Now cleaning up...");

    if ( my_XsecDB ) {
        delete my_XsecDB;
        my_XsecDB = 0;
    }

    if ( objTool ) {
        delete objTool;
        objTool= 0;
    }

    if ( m_grl ) {
        delete m_grl;
        m_grl= 0;
    }

    if ( m_Pileup) {
    delete m_Pileup;
    m_Pileup = 0;
  }


    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
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
    return EL::StatusCode::SUCCESS;
}



void  MyxAODAnalysis :: FindProbe (myAnalysisCollections myCollections) {

  if (m_debug) std::cout<<"Enter filter"<<std::endl;
  
  n_electrons = 0;
  n_muons = 0;
  n_poslep = 0;
  n_neglep = 0;
  n_jets = 0;
  n_bjets = 0;

  el_trigSF = 1;
  mu_trigSF = 1;

  bool cosmic = false;
  bool badjet = false;

  
  double trigSF = objTool->GetTotalElectronSF(*myCollections._electrons, false, false, true, false, objTool->TrigSingleLep());

  for (const auto& el : *myCollections._electrons) {
    if (m_debug) std::cout<<"Enter electron container"<<std::endl;
     
    if (el->auxdata< char >("baseline") == 1 &&  el->auxdata< char >("passOR") == 1) {
      //std::cout<<"el eta " <<el->eta() <<std::endl;
      objTool->TrigMatch(el, trig_list_el);
      n_electrons++;
      if (el->trackParticle()->charge() > 0) n_poslep++;
      else n_neglep++; 
      myElectron = el;
    }
  }

  for (const auto& mu : *myCollections._muons) {

     if (mu->auxdata< char >("baseline") == 1 &&  mu->auxdata< char >("passOR") == 1) {
       //std::cout<<"mu eta " <<mu->eta() <<std::endl;
       objTool->TrigMatch(mu, trig_list_mu);
       n_muons++;
       if (mu->primaryTrackParticle()->charge() > 0) n_poslep++;
       else n_neglep++;
       myMuon = mu;
       if ( mu->auxdata<char>("cosmic") == 1) {
	 cosmic = true;
       }
     }
  }

  if (isMC) {
    if(whichYear == 2015)
      mu_trigSF = objTool->GetTotalMuonTriggerSF( *myCollections._muons, "HLT_mu20_iloose_L1MU15_OR_HLT_mu50"); 
    else mu_trigSF = objTool->GetTotalMuonTriggerSF( *myCollections._muons, "HLT_mu26_ivarmedium_OR_HLT_mu50");         
  } 
  if (mu_trigSF==0) mu_trigSF =1;
    

  for (const auto& jet : *myCollections._jets) {
    if ( (int)jet->auxdata< char >("bad"))  badjet = true;
    if ( jet->auxdata< char >("signal") == 1  && jet->pt() > 25000. && ( fabs( jet->eta()) < 2.5) && jet->auxdata< char >("passOR") == 1 ) n_jets++;
    objTool->IsBJet(*jet);
    if ( jet->auxdata< char >("signal") == 1  && jet->pt() > 25000. && ( fabs( jet->eta()) < 2.5) && jet->auxdata< char >("passOR") == 1 && jet->auxdata< char >("bjet") == 1) n_bjets++;
  }
  
  //std::cout << "n bjets " << n_bjets << std::endl;
  
  if (m_debug) std::cout<<"n_electrons "<<n_electrons<<" n_muons "<<n_muons<<" n_poslep "<<n_poslep<< " n_neglep "<<n_neglep<<std::endl;
  

  m_flavProbe = -1;
  m_ptProbe = -1;
  m_etaProbe = -1;
  m_isTightProbe = false;

  //std::cout << "pt probe: " << m_ptProbe << std::endl;

  if(cosmic) return;
  if(badjet) return;
  //if(n_jets <2) return;
  if(n_bjets <1) return;
  if( n_electrons != 1 || n_muons != 1 || (n_poslep != 2 && n_neglep !=2)) return;
  //if( n_electrons != 1 || n_muons != 1) return;



  double iso_topo_el = myElectron->isolation( xAOD::Iso::topoetcone20);
  double iso_var_el = myElectron->isolation( xAOD::Iso::ptvarcone20);

  double iso_topo_mu = myMuon->isolation( xAOD::Iso::topoetcone20);
  double iso_var_mu = myMuon->isolation( xAOD::Iso::ptvarcone30);

  //std::cout << "mu topo " << iso_topo_mu << " mu va " << iso_var_mu << std::endl;
    


  if(isMC) {
    int muon_origin;
    int electron_origin;

    electron_origin = acc_Origin(*myElectron);
    const xAOD::TrackParticle* trackParticle = myMuon->primaryTrackParticle();
    if (trackParticle) {
      if (acc_Origin.isAvailable(*trackParticle)) muon_origin = acc_Origin(*trackParticle);
      else muon_origin = 0;
    }
    else muon_origin = 0;

    if(!doClassification) {
      if ( muon_origin != 10 &&  muon_origin != 12 &&  muon_origin != 13 &&  muon_origin != 43) {    
        return; //avoid fakes double counting
      }
      else if ( electron_origin != 10 &&  electron_origin != 12 &&  electron_origin != 13 &&  electron_origin != 43) {    
        return; //avoid fakes double counting
      }
    }
  }


  if( myMuon->pt() > myElectron->pt() ) {
    bool is_tag_mu = false;
    if ( myMuon->pt() > 30000 && myMuon->auxdata< char >("signal") && iso_var_mu/(myMuon->pt()) <0.005) is_tag_mu = true;

    bool is_match_mu = false;
    int n_trig_matched = 0;
    for (UInt_t itrig = 0; itrig < trig_list_mu.size(); itrig++) {
        if(myMuon->auxdata<int>(trig_list_mu.at(itrig))) n_trig_matched++ ;
      }
    if(n_trig_matched>0) is_match_mu = true;

    if (is_tag_mu && is_match_mu ) {
      m_ptProbe = myElectron->pt();
      m_etaProbe = myElectron->eta();
      m_flavProbe = 11 * myElectron->trackParticle()->charge();
      if(myElectron->auxdata< char >("signal")) m_isTightProbe = true;
    }
  }
  else {

    bool is_tag_ele = false;
    if (myElectron->pt() > 30000 && myElectron->auxdata< char >("signal") && iso_topo_el/(myElectron->pt()) <0.005  && iso_var_el/(myElectron->pt()) <0.005) is_tag_ele = true;

    bool is_match_ele =false;
    int n_trig_matched = 0;
    for (UInt_t itrig = 0; itrig < trig_list_el.size(); itrig++) {
        if(myElectron->auxdata<int>(trig_list_el.at(itrig))) n_trig_matched++ ;
      }
    if(n_trig_matched>0) is_match_ele = true;

    if (is_tag_ele && is_match_ele ) {
      m_ptProbe = myMuon->pt();
      m_etaProbe = myMuon->eta();
      m_flavProbe = 13 * myMuon->primaryTrackParticle()->charge();
      if(myMuon->auxdata< char >("signal")) m_isTightProbe = true;
    }
  }
  //if(m_ptProbe != -1) std::cout<< "eta probe " << m_etaProbe <<std::endl;
  
  return;

}



void  MyxAODAnalysis :: Fill_classification() {
  //check that I've selected a probe
  if (m_ptProbe == -1)return;

  double lumi_weight = 1;
  if(isMC) lumi_weight = target_lumi/sumOfWeights;
  double weight = m_xsec * lumi_weight * m_PUwei * m_WeightEvents * m_eleSF_weight * m_muSF_weight * m_bTag_weight * m_JVT_weight * el_trigSF * mu_trigSF;


  int myLepType = 0;
  // 1 -> prompt electron
  // 2 -> HF electron 
  // 3 -> LF electron
  // 4 -> prompt muon
  // 5 -> HF muon
  // 6 -> LF muon

  //electrons classification
  if(fabs(m_flavProbe) == 11) {
    
     int origin = acc_Origin(*myElectron);
     int type = acc_Type(*myElectron);
     int bkgMotherPdgId = myElectron->auxdata<int>("bkgMotherPdgId");
     
     if(type==2) myLepType = 1;
     else if(origin==5 && fabs(bkgMotherPdgId)==11) myLepType = 1;
     else if(origin==26 || origin==33 || origin==25 || origin==32) myLepType = 2;
     else myLepType = 3;

     //if(myLepType==3) std::cout << "electron type " << type << " origin " << origin << std::endl;
     h_classVSpt_ele->Fill(myLepType,m_ptProbe/1000,weight);
    }
  
  //muon classification

  if(fabs(m_flavProbe) == 13) {

    int origin = 0;
    int type = 0;

    const xAOD::TrackParticle* trackParticle = myMuon->primaryTrackParticle();
    if (trackParticle) {
     origin = acc_Origin(*trackParticle);
     type = acc_Type(*trackParticle);
    }
     if(type==6) myLepType = 4;
     else if(origin==26 || origin==33 || origin==25 || origin==32) myLepType = 5;
     else myLepType = 6;

     //std::cout << "muon type " << type << " origin " << origin << std::endl;
     h_classVSpt_mu->Fill(myLepType-3,m_ptProbe/1000,weight);

  }

  h_class->Fill(myLepType,weight);


}




void  MyxAODAnalysis :: Fill_histos() {
  
  double lumi_weight = 1;
  if(isMC) lumi_weight = target_lumi/sumOfWeights;
  double weight = m_xsec * lumi_weight * m_PUwei * m_WeightEvents * m_eleSF_weight * m_muSF_weight * m_bTag_weight * m_JVT_weight * el_trigSF * mu_trigSF;

 // std::cout <<  m_xsec << " " << lumi_weight << " " << m_PUwei << " " << m_WeightEvents << " " << m_eleSF_weight << " " << m_muSF_weight << " " << m_bTag_weight << " " << m_JVT_weight << " " << el_trigSF << " " << mu_trigSF << std::endl;
 
 
 // std::cout << "pt probe: " << m_ptProbe << std::endl;

  if (m_ptProbe == -1)return;
 
  //std::cout << "pt probe: " << m_ptProbe << std::endl;

  if(fabs(m_flavProbe) == 11) {
    histo_loose_el->Fill_values(0, m_ptProbe/1000);
    histo_loose_el->Fill_values(1,std::abs(m_etaProbe));
    histo_loose_el->Fill_weights(weight);

    if(m_isTightProbe) {
      histo_tight_el->Fill_values(0, m_ptProbe/1000);
      histo_tight_el->Fill_values(1,std::abs(m_etaProbe));
      histo_tight_el->Fill_weights(weight);
    }
  }

  if(fabs(m_flavProbe) == 13){
    histo_loose_mu->Fill_values(0, m_ptProbe/1000);
    histo_loose_mu->Fill_values(1,std::abs(m_etaProbe));
    histo_loose_mu->Fill_weights(weight);

    if(m_isTightProbe){
      histo_tight_mu->Fill_values(0, m_ptProbe/1000);
      histo_tight_mu->Fill_values(1,std::abs(m_etaProbe));
      histo_tight_mu->Fill_weights(weight);
    }
  }

  return;

}

void MyxAODAnalysis :: Register_parametric_histos(parametric_histos *histo)
{

  for (int i = 0; i < histo->GetNhistos(); i++) {  //THIS NEEDS TO BE IMPLEMENTED
    wk()->addOutput(histo->Get1Dhisto(i));
  }
  
  TH1F *h_all = histo->GetGeneralHisto();
  wk()->addOutput(h_all);
  
  return;
}



