/**
 * @file OJEfficiencyLoader.cxx
 * @author Thomas Gillam <thomas.gillam@cern.ch>
 * @date May, 2013
 * @brief Load fake rates computed by Julien & Otilia
 *
 **/ 

#include "GeneralisedMatrixMethod/OJEfficiencyLoader.h"

namespace GeneralisedMatrixMethod {
#include "TwoLepSSFakeBkg_Tables.hxx"
}

#include <cmath>


namespace GeneralisedMatrixMethod {
  /*************************************************************************
          Accessors of numerical results stored in TwoLepSSFakeBkg_Tables.hxx
   *************************************************************************/
  OJEfficiencyLoader::Efficiency OJEfficiencyLoader::realEfficiencyEl(float eta, float ptInGeV, float drjet)
  {
    Efficiency eff;
    eta=fabs(eta);
    unsigned int i_eta=0,i_pt=0, i_dr=0;  // JO: added i_dr
    for(;i_eta<Params_RealEffEl_nEtaBins-1;i_eta++)
      if(eta<=Params_RealEffEl_EtaBins[i_eta+1]) break;
    for(;i_pt<Params_RealEffEl_nPtBins-1;i_pt++)
      if(ptInGeV<=Params_RealEffEl_PtBins[i_pt+1]) break;
    for(;i_dr<Params_RealEffEl_ndRBins-1;i_dr++) // JO: added
      if(drjet<=Params_RealEffEl_dRBins[i_dr+1]) break; // JO: added
    unsigned int bin=i_pt*Params_RealEffEl_nEtaBins+i_eta;
	unsigned int bin_dr = i_pt*Params_RealEffEl_ndRBins + i_dr; // JO: added
    eff.e = Params_RealEffEl_Eff[bin];
    eff.eStat = Params_RealEffEl_Stat[bin];
    eff.eSystUncorr = Params_RealEffEl_SystUncorr[bin];
	double ec = Params_RealEffEl_Eff[bin]*Params_RealEffEl_SystCorr_relative[bin]; // JO: new
	double eb = Params_RealEffEl_Eff[bin]*Params_RealEffEl_SystBusy_relative[bin_dr]; // JO: new
    eff.eSystCorr = sqrt(ec*ec + eb*eb); // JO: new
    eff.bin_index = -1; // JO: unused
    return eff;		
  }

  OJEfficiencyLoader::Efficiency OJEfficiencyLoader::realEfficiencyMu(float eta, float ptInGeV, float drjet)
  {
    Efficiency eff;
    eta=fabs(eta);
    unsigned int i_eta=0,i_pt=0, i_dr=0; // JO: added i_dr
    for(;i_eta<Params_RealEffMu_nEtaBins-1;i_eta++)
      if(eta<=Params_RealEffMu_EtaBins[i_eta+1]) break;
    for(;i_pt<Params_RealEffMu_nPtBins-1;i_pt++)
      if(ptInGeV<=Params_RealEffMu_PtBins[i_pt+1]) break;
    for(;i_dr<Params_RealEffMu_ndRBins-1;i_dr++) // JO: added
      if(drjet<=Params_RealEffMu_dRBins[i_dr+1]) break; // JO: added
    unsigned int bin=i_pt*Params_RealEffMu_nEtaBins+i_eta;
	unsigned int bin_dr = i_pt*Params_RealEffMu_ndRBins + i_dr; // JO: added
    eff.e = Params_RealEffMu_Eff[bin];
    eff.eStat = Params_RealEffMu_Stat[bin];
    eff.eSystUncorr = Params_RealEffMu_SystUncorr[bin];
	double ec = Params_RealEffMu_Eff[bin]*Params_RealEffMu_SystCorr_relative[bin]; // JO: new
	double eb = Params_RealEffMu_Eff[bin]*Params_RealEffMu_SystBusy_relative[bin_dr]; // JO: new
    eff.eSystCorr = sqrt(ec*ec + eb*eb); // JO: new
    eff.bin_index = -1; // JO: unused
    return eff;
  }
   
  OJEfficiencyLoader::Efficiency OJEfficiencyLoader::fakeRateEl(float eta, float ptInGeV, float, int)
  {
    Efficiency eff;
    eta=fabs(eta);
    bool isExtrapolated=false;
    unsigned int i=0;
    for(;;i++)
    {
      const Bin2D&b=Params_FakeRateEl_Bins[i];
      if(ptInGeV>=b.minPt && ptInGeV<b.maxPt && eta>=b.minEta && eta<b.maxEta) break;
      if(i<(Params_FakeRateEl_nBins-1)) continue;
      isExtrapolated=true;
      break;
    }
	// JO: return same result regardless of bjet multiplicity
    // if(hasBJet)
    // {
      // eff.e = Params_FakeRateEl_hasB_Eff[i];
      // eff.bin_index = i;
      // eff.eStat = Params_FakeRateEl_hasB_Stat[i];
      // eff.eSystUncorr = Params_FakeRateEl_hasB_SystUncorr[i];
      // //if(isExtrapolated) eff.ESystExtrap = Params_FakeRateEl_hasB_SystOverflow;
      // //else eff.ESystExtrap = 0;
      // eff.eSystCorr = Params_FakeRateEl_hasB_SystCorr[i];
    // }
    // else
    // {
      // eff.e = Params_FakeRateEl_noB_Eff[i];
      // eff.bin_index = i;
      // eff.eStat = Params_FakeRateEl_noB_Stat[i];
      // eff.eSystUncorr = Params_FakeRateEl_noB_SystUncorr[i];
      // //if(isExtrapolated) eff.ESystExtrap = Params_FakeRateEl_noB_SystOverflow;
      // //else eff.ESystExtrap = 0;
      // eff.eSystCorr = Params_FakeRateEl_noB_SystCorr[i];	
    // }
    eff.e = Params_FakeRateEl_Eff[i];
    eff.bin_index = i;
    eff.eStat = Params_FakeRateEl_Stat[i];
    eff.eSystUncorr = Params_FakeRateEl_SystUncorr[i];
    eff.eSystCorr = Params_FakeRateEl_SystCorr[i];
    //if(isExtrapolated) eff.ESystExtrap = Params_FakeRateEl_hasB_SystOverflow;
    //else eff.ESystExtrap = 0;
    return eff;
  }

  OJEfficiencyLoader::Efficiency OJEfficiencyLoader::fakeRateMu(float eta, float ptInGeV,  float, int)
  {
    Efficiency eff;
    eta=fabs(eta);
    bool isExtrapolated=false;
    unsigned int i=0;
    for(;;i++)
    {
      const Bin2D&b=Params_FakeRateMu_Bins[i];
      if(ptInGeV>=b.minPt && ptInGeV<b.maxPt && eta>=b.minEta && eta<b.maxEta) break;
      if(i<(Params_FakeRateMu_nBins-1)) continue;
      isExtrapolated=true;
      break;
    }
    eff.e = Params_FakeRateMu_Eff[i];
    eff.bin_index = i;
    eff.eStat = Params_FakeRateMu_Stat[i];
    eff.eSystUncorr = Params_FakeRateMu_SystUncorr[i];
    //if(isExtrapolated) eff.ESystExtrap = Params_FakeRateMu_SystOverflow;
    //else eff.ESystExtrap = 0;
    eff.eSystCorr = Params_FakeRateMu_SystCorr[i];
    return eff;
  }

  unsigned int OJEfficiencyLoader::numFakeBinsEl() const
  {
    return Params_FakeRateEl_nBins;
  }

  unsigned int OJEfficiencyLoader::numFakeBinsMu() const
  {
    return Params_FakeRateMu_nBins;
  }
}
