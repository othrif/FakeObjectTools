/**
 * @file FakeLeptBkg/Root/TwoLepSSFakeBkg.cxx
 * @author Julien Maurer <jmaurer@cern.ch>  Otilia Ducu <oducu@cern.ch>
 * @date May, 2012
 * @brief Fake background estimation for same-sign di-lepton analysis: matrix method and charge mis-ID  
 **/

#include <iostream>
using std::cout;
using std::endl;

#include "FakeLeptBkg/TwoLepSSFakeBkg.h"
#include "TwoLepSSFakeBkg_Tables.hxx"

#include <cmath>

// Functions returning efficiencies
struct EFF
{
	double E,EStat,ESystCorr,ESystUncorr,ESystExtrap;
	int bin_index;
};

EFF GetRealEff_El(float eta,float pt);
EFF GetRealEff_Mu(float eta,float pt);
EFF GetChargeMisID(bool isTight,float eta,float pt);
EFF GetFakeRate_El(float eta,float pt,bool hasBJet);
EFF GetFakeRate_Mu(float eta,float pt,bool hasBJet);


// Matrix method computation
double __mxm_GetWeight(int channel,double lep1Eta,double lep1Pt,bool lep1IsTight,double lep2Eta,double lep2Pt,
		bool lep2IsTight,bool use3DMxM,bool hasBJet,double*wgt_ErrBySource,bool splitErrBySource);
double __mxm4_GetNfake(double npp,double npf,double nfp,double nff,double f1,double f2,double r1,double r2);
void __mxm4_GetErrNFake(double*ErrBySource,double npp,double npf,double nfp,double nff,const EFF&f1,const EFF&f2,const EFF&r1,const EFF&r2,int channel,bool hasBJet);
double __mxm3_GetNfake(double npp,double npf,double nfp,double f1,double f2,double r1,double r2);
double __mxm3_GetErrNFake(double npp,double npf,double nfp,double f1,double f2,double r1,double r2,
			  double err_f1,double err_f2,double err_r1,double err_r2,bool UseCorrelErr);

// Constants 

#define CHANNEL_ELEL 0
#define CHANNEL_MUMU 1
#define CHANNEL_ELMU 2

// Error sources
#define ERRSOURCE_MXM_EL_UNCORR (0)
#define ERRSOURCE_MXM_EL_CORR (ERRSOURCE_MXM_EL_UNCORR+2*Params_FakeRateEl_nBins)
#define ERRSOURCE_MXM_MU_UNCORR (ERRSOURCE_MXM_EL_CORR+1)
//#define ERRSOURCE_MXM_MU_UNCORR_LAST (ERRSOURCE_MXM_MU_UNCORR+Params_FakeRateMu_nBins)
#define ERRSOURCE_MXM_MU_CORR (ERRSOURCE_MXM_MU_UNCORR+Params_FakeRateMu_nBins)
#define ERRSOURCE_MXM_NERRORS (ERRSOURCE_MXM_MU_CORR+1)
#define ERRSOURCE_MXM_EL_NERRORS (ERRSOURCE_MXM_EL_CORR+1)
#define ERRSOURCE_MXM_MU_NERRORS (ERRSOURCE_MXM_NERRORS-ERRSOURCE_MXM_EL_NERRORS)


#define ERRSOURCE_MISID_NERRORS (2*Params_ChargeMisID_Tight_nEtaBins*Params_ChargeMisID_Tight_nPtBins + 1)

/******************************************
          Class member methods
 ******************************************/

const unsigned int TwoLepSSFakeBkg::nErrorSources_MxM(ERRSOURCE_MXM_NERRORS);
const unsigned int TwoLepSSFakeBkg::nErrorSources_MxM_El(ERRSOURCE_MXM_EL_NERRORS);
const unsigned int TwoLepSSFakeBkg::nErrorSources_MxM_Mu(ERRSOURCE_MXM_MU_NERRORS);
const unsigned int TwoLepSSFakeBkg::nErrorSources_OStoSS(ERRSOURCE_MISID_NERRORS);
 
TwoLepSSFakeBkg::TwoLepSSFakeBkg(void)
{
	m_SplitErrBySource_MxM = false;
	m_SplitErrBySource_OStoSS = false;
	m_UseMxM3D = false;
}

TwoLepSSFakeBkg::~TwoLepSSFakeBkg(void)
{
}

void TwoLepSSFakeBkg::Use3DMatrixMethod(bool useMxM3D)
{
	m_UseMxM3D = useMxM3D;
}

void TwoLepSSFakeBkg::UseSplitErrBySource_MxM(bool split)
{
	m_SplitErrBySource_MxM = split;
}

void TwoLepSSFakeBkg::UseSplitErrBySource_OStoSS(bool split)
{
	m_SplitErrBySource_OStoSS = split;
}

double TwoLepSSFakeBkg::GetWeight_OStoSS_ElMu(double elEta,double elPt,bool elIsTight,double *wgt_err) const
{
	const EFF w = GetChargeMisID(elIsTight,elEta,elPt);
	if(wgt_err)
	{
		if(m_SplitErrBySource_OStoSS)
		{
			wgt_err[0] = sqrt(w.ESystCorr*w.ESystCorr + w.ESystExtrap*w.ESystExtrap);
			for(int i=1;i<ERRSOURCE_MISID_NERRORS;++i)
			{
				if((i-1)!=w.bin_index) wgt_err[i] = 0;
				else wgt_err[i] = sqrt(w.EStat*w.EStat + w.ESystUncorr*w.ESystUncorr);
			}
		}
		else
		{
			if(wgt_err) *wgt_err = sqrt(w.EStat*w.EStat + w.ESystUncorr*w.ESystUncorr + w.ESystCorr*w.ESystCorr + w.ESystExtrap*w.ESystExtrap);			
		}
	}
	return w.E;
}

double TwoLepSSFakeBkg::GetWeight_OStoSS_ElEl(double el1Eta,double el1Pt,bool el1IsTight,double el2Eta,double el2Pt,bool el2IsTight,double *wgt_err) const
{
	const EFF w1 = GetChargeMisID(el1IsTight,el1Eta,el1Pt);
	const EFF w2 = GetChargeMisID(el2IsTight,el2Eta,el2Pt);
	if(wgt_err)
	{
		if(m_SplitErrBySource_OStoSS)
		{
			const double K1 = (1-2*w2.E), K2 = (1-2*w1.E);
			double w1_errSq = w1.EStat*w1.EStat + w1.ESystUncorr*w1.ESystUncorr + w1.ESystCorr*w1.ESystCorr + w1.ESystExtrap*w1.ESystExtrap;
			if(w1.bin_index!=w2.bin_index)
			{
				wgt_err[0] = sqrt(K1*K1*(w1.ESystCorr*w1.ESystCorr+w1.ESystExtrap*w1.ESystExtrap) + K2*K2*(w2.ESystCorr*w2.ESystCorr+w2.ESystExtrap*w2.ESystExtrap));
				for(int i=1;i<ERRSOURCE_MISID_NERRORS;++i)
				{
					if((i-1)==w1.bin_index) wgt_err[i] = K1*sqrt(w1.EStat*w1.EStat + w1.ESystUncorr*w1.ESystUncorr);
					else if((i-1)==w2.bin_index) wgt_err[i] = K2*sqrt(w2.EStat*w2.EStat + w2.ESystUncorr*w2.ESystUncorr);
					else wgt_err[i] = 0.;
				}
			}
			else
			{
				wgt_err[0] = K1*sqrt(w1.ESystCorr*w1.ESystCorr+w1.ESystExtrap*w1.ESystExtrap) 
					+ K2*sqrt(w2.ESystCorr*w2.ESystCorr+w2.ESystExtrap*w2.ESystExtrap);
				for(int i=1;i<ERRSOURCE_MISID_NERRORS;++i)
				{
					if((i-1)==w1.bin_index) wgt_err[i] = K1*sqrt(w1.EStat*w1.EStat + w1.ESystUncorr*w1.ESystUncorr) 
							+ K2*sqrt(w2.EStat*w2.EStat + w2.ESystUncorr*w2.ESystUncorr);
					else wgt_err[i] = 0.;
				}
			}
		}
		else
		{
			double w1_errSq = w1.EStat*w1.EStat + w1.ESystUncorr*w1.ESystUncorr + w1.ESystCorr*w1.ESystCorr + w1.ESystExtrap*w1.ESystExtrap;
			if(w1.bin_index!=w2.bin_index)
			{
				double w2_errSq = w2.EStat*w2.EStat + w2.ESystUncorr*w2.ESystUncorr + w2.ESystCorr*w2.ESystCorr + w2.ESystExtrap*w2.ESystExtrap;
				*wgt_err = sqrt(pow(1-2*w2.E,2)*w1_errSq + pow(1-2*w1.E,2)*w2_errSq);
			}
			else *wgt_err = fabs(2.-4.*w1.E)*sqrt(w1_errSq);
		}
	}
	return w1.E+w2.E-2*w1.E*w2.E;  
}

double TwoLepSSFakeBkg::GetWeight_MxM_ElEl(double el1Eta,double el1Pt,bool el1IsTight,double el2Eta,double el2Pt,bool el2IsTight,bool hasBJet,double *wgt_err) const
{
	return __mxm_GetWeight(CHANNEL_ELEL,el1Eta,el1Pt,el1IsTight,el2Eta,el2Pt,el2IsTight,m_UseMxM3D,hasBJet,wgt_err,m_SplitErrBySource_MxM);
}

double TwoLepSSFakeBkg::GetWeight_MxM_ElMu(double elEta,double elPt,bool elIsTight,double muEta,double muPt,bool muIsTight,bool hasBJet,double *wgt_err) const
{
	return __mxm_GetWeight(CHANNEL_ELMU,elEta,elPt,elIsTight,muEta,muPt,muIsTight,m_UseMxM3D,hasBJet,wgt_err,m_SplitErrBySource_MxM);
}

double TwoLepSSFakeBkg::GetWeight_MxM_MuMu(double mu1Eta,double mu1Pt,bool mu1IsTight,double mu2Eta,double mu2Pt,bool mu2IsTight,bool hasBJet,double *wgt_err) const
{
	return __mxm_GetWeight(CHANNEL_MUMU,mu1Eta,mu1Pt,mu1IsTight,mu2Eta,mu2Pt,mu2IsTight,m_UseMxM3D,hasBJet,wgt_err,m_SplitErrBySource_MxM);
}

/******************************************
          Internal helper functions
 ******************************************/

// Apply matrix method
double __mxm_GetWeight(int channel,double lep1Eta,double lep1Pt,bool lep1IsTight,double lep2Eta,double lep2Pt,
		bool lep2IsTight,bool use3DMxM,bool hasBJet,double*wgt_ErrBySource,bool splitErrBySource)
{
	EFF f1,f2,r1,r2;
	switch(channel)
	{
	case CHANNEL_ELEL:
		r1 = GetRealEff_El(lep1Eta,lep1Pt);
		r2 = GetRealEff_El(lep2Eta,lep2Pt);
		f1 = GetFakeRate_El(lep1Eta,lep1Pt,hasBJet);
		f2 = GetFakeRate_El(lep2Eta,lep2Pt,hasBJet);
		break;
	case CHANNEL_ELMU:
		r1 = GetRealEff_El(lep1Eta,lep1Pt);
		r2 = GetRealEff_Mu(lep2Eta,lep2Pt);
		f1 = GetFakeRate_El(lep1Eta,lep1Pt,hasBJet);
		f2 = GetFakeRate_Mu(lep2Eta,lep2Pt,hasBJet);
		break;
	case CHANNEL_MUMU:
		r1 = GetRealEff_Mu(lep1Eta,lep1Pt);
		r2 = GetRealEff_Mu(lep2Eta,lep2Pt);
		f1 = GetFakeRate_Mu(lep1Eta,lep1Pt,hasBJet);
		f2 = GetFakeRate_Mu(lep2Eta,lep2Pt,hasBJet);
		break;
	}
	double npp=0,npf=0,nfp=0,nff=0;
	if(lep1IsTight && lep2IsTight) npp=1;
	else if(lep1IsTight) npf=1;
	else if(lep2IsTight) nfp=1;
	else nff=1;
	double wgt=0;
	if(use3DMxM)
	{
		if(nff>0) return 0;
		wgt = __mxm3_GetNfake(npp,npf,nfp,f1.E,f2.E,r1.E,r2.E);
		double r1_err = sqrt(pow(r1.EStat,2)+pow(r1.ESystUncorr,2)+pow(r1.ESystExtrap,2)+pow(r1.ESystUncorr,2));
		double r2_err = sqrt(pow(r2.EStat,2)+pow(r2.ESystUncorr,2)+pow(r2.ESystExtrap,2)+pow(r2.ESystUncorr,2));
		double f1_err = sqrt(pow(f1.EStat,2)+pow(f1.ESystUncorr,2)+pow(f1.ESystExtrap,2)+pow(f1.ESystUncorr,2));
		double f2_err = sqrt(pow(f2.EStat,2)+pow(f2.ESystUncorr,2)+pow(f2.ESystExtrap,2)+pow(f2.ESystUncorr,2));
		if(wgt_ErrBySource) wgt_ErrBySource[0] = __mxm3_GetErrNFake(npp,npf,nfp,f1.E,f2.E,r1.E,r2.E,f1_err,f2_err,r1_err,r2_err,channel!=CHANNEL_ELMU);
	}
	else
	{
		wgt = __mxm4_GetNfake(npp,npf,nfp,nff,f1.E,f2.E,r1.E,r2.E);
		if(wgt_ErrBySource)
		{
			if(splitErrBySource)
			{
				for(unsigned int i=0;i<ERRSOURCE_MXM_NERRORS;i++) wgt_ErrBySource[i]=0;			
				__mxm4_GetErrNFake(wgt_ErrBySource,npp,npf,nfp,nff,f1,f2,r1,r2,channel,hasBJet);
			}
			else
			{
				double ErrBySource[ERRSOURCE_MXM_NERRORS];
				for(unsigned int i=0;i<ERRSOURCE_MXM_NERRORS;i++) ErrBySource[i]=0;
				__mxm4_GetErrNFake(ErrBySource,npp,npf,nfp,nff,f1,f2,r1,r2,channel,hasBJet);
				double errSq=0;
				for(unsigned int i=0;i<ERRSOURCE_MXM_NERRORS;i++) errSq += ErrBySource[i]*ErrBySource[i];
				wgt_ErrBySource[0] = sqrt(errSq);
			}
		}
	}
	return wgt;
}

// Sum all contributions with possibly one fake-lepton
double __mxm4_GetNfake(double npp,double npf,double nfp,double nff,double f1,double f2,double r1,double r2) {
  double nrealfake = r1*f2/(r1-f1)/(r2-f2)*((f1-1.)*(1.-r2)*npp+(1.-f1)*r2*npf+f1*(1.-r2)*nfp-f1*r2*nff);
  double nfakereal = f1*r2/(r1-f1)/(r2-f2)*((r1-1.)*(1.-f2)*npp+(1.-r1)*f2*npf+r1*(1.-f2)*nfp-r1*f2*nff);
  double nfakefake = f1*f2/(r1-f1)/(r2-f2)*((1.-r1)*(1.-r2)*npp+(r1-1.)*r2*npf+r1*(r2-1.)*nfp+r1*r2*nff);
  return nrealfake + nfakereal + nfakefake;
}

void __mxm4_GetErrNFake(double*ErrBySource,double npp,double npf,double nfp,double nff,const EFF&f1,const EFF&f2,const EFF&r1,const EFF&r2,int channel,bool hasBJet)
{
	if(!ErrBySource) return;
	const double dNFr1 = -1*(f1.E*(npp-f2.E*(npf+npp)+f1.E*(-nfp-npp+f2.E*(nff+nfp+npf+npp)))*r2.E)/((f1.E-r1.E)*(f1.E-r1.E)*(f2.E-r2.E));
	const double dNFr2 = -1*(f2.E*(npp-f2.E*(npf+npp)+f1.E*(-nfp-npp+f2.E*(nff+nfp+npf+npp)))*r1.E)/((f1.E-r1.E)*(f2.E-r2.E)*(f2.E-r2.E));
	const double dNFf1 = (r1.E*(npp-nfp*r1.E-npp*r1.E+f2.E*(npf*(-1+r1.E)+npp*(-1+r1.E)+(nff+nfp)*r1.E))*r2.E)/((f1.E-r1.E)*(f1.E-r1.E)*(f2.E-r2.E));
	const double dNFf2 = (r1.E*r2.E*(npp-npf*r2.E-npp*r2.E+f1.E*(nfp*(-1+r2.E)+npp*(-1+r2.E)+(nff+npf)*r2.E)))/((f1.E-r1.E)*(f2.E-r2.E)*(f2.E-r2.E));
	const double err_r1 = sqrt(r1.EStat*r1.EStat + r1.ESystUncorr*r1.ESystUncorr + r1.ESystCorr*r1.ESystCorr);
	const double err_r2 = sqrt(r2.EStat*r2.EStat + r2.ESystUncorr*r2.ESystUncorr + r2.ESystCorr*r2.ESystCorr);
	int bjoff = hasBJet?1:0;
	switch(channel)
	{
		case CHANNEL_ELEL:
			if(f1.bin_index!=f2.bin_index)
			{
				ErrBySource[ERRSOURCE_MXM_EL_UNCORR+2*f1.bin_index+bjoff] = fabs(dNFf1)*sqrt(pow(f1.EStat,2)+pow(f1.ESystUncorr,2));
				ErrBySource[ERRSOURCE_MXM_EL_UNCORR+2*f2.bin_index+bjoff] = fabs(dNFf2)*sqrt(pow(f2.EStat,2)+pow(f2.ESystUncorr,2));
			}
			else ErrBySource[ERRSOURCE_MXM_EL_UNCORR+2*f1.bin_index+bjoff] = fabs(dNFf1+dNFf2)*sqrt(pow(f1.EStat,2)+pow(f1.ESystUncorr,2));		
			//ErrBySource[ERRSOURCE_MXM_EL_UNCORR_LAST+bjoff] = fabs(dNFf1*f1.ESystExtrap + dNFf2*f2.ESystExtrap);	// not implemented for the moment
			ErrBySource[ERRSOURCE_MXM_EL_CORR] = sqrt(pow(dNFf1*f1.ESystCorr+dNFf2*f2.ESystCorr,2) + pow(fabs(dNFr1)*err_r1+fabs(dNFr2)*err_r2,2));
			break;
		case CHANNEL_MUMU:
			if(f1.bin_index!=f2.bin_index)
			{
				ErrBySource[ERRSOURCE_MXM_MU_UNCORR+f1.bin_index] = fabs(dNFf1)*sqrt(pow(f1.EStat,2)+pow(f1.ESystUncorr,2));
				ErrBySource[ERRSOURCE_MXM_MU_UNCORR+f2.bin_index] = fabs(dNFf2)*sqrt(pow(f2.EStat,2)+pow(f2.ESystUncorr,2));
			}
			else ErrBySource[ERRSOURCE_MXM_MU_UNCORR+f1.bin_index] = fabs(dNFf1+dNFf2)*sqrt(pow(f1.EStat,2)+pow(f1.ESystUncorr,2));				
			//ErrBySource[ERRSOURCE_MXM_MU_UNCORR_LAST] = fabs(dNFf1*f1.ESystExtrap + dNFf2*f2.ESystExtrap);
			ErrBySource[ERRSOURCE_MXM_MU_CORR] = sqrt(pow(dNFf1*f1.ESystCorr+dNFf2*f2.ESystCorr,2) + pow(fabs(dNFr1)*err_r1+fabs(dNFr2)*err_r2,2));
			break;
		case CHANNEL_ELMU:
			ErrBySource[ERRSOURCE_MXM_EL_UNCORR+2*f1.bin_index+bjoff] = fabs(dNFf1)*sqrt(pow(f1.EStat,2)+pow(f1.ESystUncorr,2));
			ErrBySource[ERRSOURCE_MXM_MU_UNCORR+f2.bin_index] = fabs(dNFf2)*sqrt(pow(f2.EStat,2)+pow(f2.ESystUncorr,2));
			//ErrBySource[ERRSOURCE_MXM_EL_UNCORR_LAST+bjoff] = fabs(dNFf1)*f1.ESystExtrap;				
			//ErrBySource[ERRSOURCE_MXM_MU_UNCORR_LAST] = fabs(dNFf2)*f2.ESystExtrap;
			ErrBySource[ERRSOURCE_MXM_EL_CORR] = sqrt(pow(dNFf1*f1.ESystCorr,2) + pow(dNFr1*err_r1,2));
			ErrBySource[ERRSOURCE_MXM_MU_CORR] = sqrt(pow(dNFf2*f2.ESystCorr,2) + pow(dNFr1*err_r2,2));
			break;
	}
}

// Assumes there are no event with 2 fakes, always one true lepton
double __mxm3_GetNfake(double npp,double npf,double nfp,double f1,double f2,double r1,double r2) {
  double nrealfake = f2*((1-r2)*npp-r2*npf)/(f2-r2);
  double nfakereal = f1*((1-r1)*npp-r1*nfp)/(f1-r1);
  return nrealfake + nfakereal;
}

double __mxm3_GetErrNFake(double npp,double npf,double nfp,double f1,double f2,double r1,double r2,
			  double err_f1,double err_f2,double err_r1,double err_r2,bool UseCorrelErr)
{
  double dNFr1 = -((f1*(-npp+f1*(nfp+npp)))/pow(f1-r1,2));
  double dNFr2 = -((f2*(-npp+f2*(npf+npp)))/pow(f2-r2,2));
  double dNFf1 = (r1*(npp*(-1+r1)+nfp*r1))/pow(f1-r1,2);
  double dNFf2 = (r2*(npp*(-1+r2)+npf*r2))/pow(f2-r2,2);
  if (UseCorrelErr)
    return sqrt(pow(dNFr1*err_r1 + dNFr2*err_r2,2)+pow(dNFf1*err_f1+dNFf2*err_f2,2));
  else
    return sqrt(pow(dNFr1*err_r1,2)+pow(dNFr2*err_r2,2)+pow(dNFf1*err_f1,2)+pow(dNFf2*err_f2,2)); 	
}

/*************************************************************************
        Accessors of numerical results stored in TwoLepSSFakeBkg_Tables.hxx
 *************************************************************************/
EFF GetRealEff_El(float eta,float pt)
{
	EFF eff = {0,0,0,0,0,0};
	eta=fabs(eta);
	unsigned int i_eta=0,i_pt=0;
	for(;i_eta<Params_RealEffEl_nEtaBins-1;i_eta++)
		if(eta<=Params_RealEffEl_EtaBins[i_eta+1]) break;
	for(;i_pt<Params_RealEffEl_nPtBins-1;i_pt++)
		if(pt<=Params_RealEffEl_PtBins[i_pt+1]) break;
	unsigned int bin=i_pt*Params_RealEffEl_nEtaBins+i_eta;
	eff.E = Params_RealEffEl_Eff[bin];
	eff.EStat = Params_RealEffEl_Stat[bin];
	eff.ESystUncorr = Params_RealEffEl_SystUncorr[bin];
	eff.ESystCorr = 0;
	eff.bin_index = bin;
	return eff;		
}

EFF GetRealEff_Mu(float eta,float pt)
{
	EFF eff = {0,0,0,0,0,0};
	eta=fabs(eta);
	unsigned int i_eta=0,i_pt=0;
	for(;i_eta<Params_RealEffMu_nEtaBins-1;i_eta++)
		if(eta<=Params_RealEffMu_EtaBins[i_eta+1]) break;
	for(;i_pt<Params_RealEffMu_nPtBins-1;i_pt++)
		if(pt<=Params_RealEffMu_PtBins[i_pt+1]) break;
	unsigned int bin=i_pt*Params_RealEffMu_nEtaBins+i_eta;
	eff.E = Params_RealEffMu_Eff[bin];
	eff.EStat = Params_RealEffMu_Stat[bin];
	eff.ESystUncorr = Params_RealEffMu_SystUncorr[bin];
	eff.ESystCorr = 0;
	eff.bin_index = bin;
	return eff;
}
 
EFF GetFakeRate_El(float eta,float pt,bool hasBJet)
{
	EFF eff = {0,0,0,0,0};
	eta=fabs(eta);
	bool isExtrapolated=false;
	unsigned int i=0;
	for(;;i++)
	{
		const Bin2D&b=Params_FakeRateEl_Bins[i];
		if(pt>=b.minPt && pt<b.maxPt && eta>=b.minEta && eta<b.maxEta) break;
		if(i<(Params_FakeRateEl_nBins-1)) continue;
		isExtrapolated=true;
		break;
	}
	/*
	if(hasBJet)
	{
		eff.E = Params_FakeRateEl_hasB_Eff[i];
		eff.bin_index = i;
		eff.EStat = Params_FakeRateEl_hasB_Stat[i];
		eff.ESystUncorr = Params_FakeRateEl_hasB_SystUncorr[i];
		//if(isExtrapolated) eff.ESystExtrap = Params_FakeRateEl_hasB_SystOverflow;
		//else eff.ESystExtrap = 0;
		eff.ESystCorr = Params_FakeRateEl_hasB_SystCorr[i];
	}
	else
	{
		eff.E = Params_FakeRateEl_noB_Eff[i];
		eff.bin_index = i;
		eff.EStat = Params_FakeRateEl_noB_Stat[i];
		eff.ESystUncorr = Params_FakeRateEl_noB_SystUncorr[i];
		//if(isExtrapolated) eff.ESystExtrap = Params_FakeRateEl_noB_SystOverflow;
		//else eff.ESystExtrap = 0;
		eff.ESystCorr = Params_FakeRateEl_noB_SystCorr[i];	
	}*/
	eff.E = Params_FakeRateEl_Eff[i];
	eff.bin_index = i;
	eff.EStat = Params_FakeRateEl_Stat[i];
	eff.ESystUncorr = Params_FakeRateEl_SystUncorr[i];
	//if(isExtrapolated) eff.ESystExtrap = Params_FakeRateEl_noB_SystOverflow;
	//else eff.ESystExtrap = 0;
	eff.ESystCorr = Params_FakeRateEl_SystCorr[i];	
	return eff;
}

EFF GetFakeRate_Mu(float eta,float pt,bool)
{
	EFF eff = {0,0,0,0,0};
	eta=fabs(eta);
	bool isExtrapolated=false;
	unsigned int i=0;
	for(;;i++)
	{
		const Bin2D&b=Params_FakeRateMu_Bins[i];
		if(pt>=b.minPt && pt<b.maxPt && eta>=b.minEta && eta<b.maxEta) break;
		if(i<(Params_FakeRateMu_nBins-1)) continue;
		isExtrapolated=true;
		break;
	}
	eff.E = Params_FakeRateMu_Eff[i];
	eff.bin_index = i;
	eff.EStat = Params_FakeRateMu_Stat[i];
	eff.ESystUncorr = Params_FakeRateMu_SystUncorr[i];
	//if(isExtrapolated) eff.ESystExtrap = Params_FakeRateMu_SystOverflow;
	//else eff.ESystExtrap = 0;
	eff.ESystCorr = Params_FakeRateMu_SystCorr[i];
	return eff;
}

EFF GetChargeMisID(bool isTight,float eta,float pt)
{	
	EFF eff = {0,0,0,0,0};
	unsigned int i_eta=0,i_pt=0;
	if(isTight)
	{
		if(Params_ChargeMisID_Tight_EtaBins[0]>=0.) eta=fabs(eta);
		for(;i_eta<Params_ChargeMisID_Tight_nEtaBins-1;i_eta++)
				if(eta<=Params_ChargeMisID_Tight_EtaBins[i_eta+1]) break;
		for(;i_pt<Params_ChargeMisID_Tight_nPtBins-1;i_pt++)
				if(pt<=Params_ChargeMisID_Tight_PtBins[i_pt+1]) break;
		unsigned int bin=i_pt*Params_ChargeMisID_Tight_nEtaBins+i_eta;	
		eff.E = Params_ChargeMisID_Tight_Rate[bin];
		eff.EStat = Params_ChargeMisID_Tight_Stat[bin];
		eff.ESystUncorr = Params_ChargeMisID_Tight_SystUncorr[bin];
		eff.ESystCorr = Params_ChargeMisID_Tight_SystCorr[bin];
		eff.bin_index = bin;
	}
	else
	{
		if(Params_ChargeMisID_Loose_EtaBins[0]>=0.) eta=fabs(eta);
		for(;i_eta<Params_ChargeMisID_Loose_nEtaBins-1;i_eta++)
				if(eta<=Params_ChargeMisID_Loose_EtaBins[i_eta+1]) break;
		for(;i_pt<Params_ChargeMisID_Loose_nPtBins-1;i_pt++)
				if(pt<=Params_ChargeMisID_Loose_PtBins[i_pt+1]) break;
		unsigned int bin=i_pt*Params_ChargeMisID_Loose_nEtaBins+i_eta;	
		eff.E = Params_ChargeMisID_Loose_Rate[bin];
		eff.EStat = Params_ChargeMisID_Loose_Stat[bin];
		eff.ESystUncorr = Params_ChargeMisID_Loose_SystUncorr[bin];
		eff.ESystCorr = Params_ChargeMisID_Loose_SystCorr[bin];
		eff.bin_index = bin + Params_ChargeMisID_Tight_nEtaBins*Params_ChargeMisID_Tight_nPtBins;		
	}
	return eff;
}

/// Static class members: description of error sources

const char* TwoLepSSFakeBkg::GetErrorSourceDescription_MxM(unsigned int index)
{
	if(ERRSOURCE_MXM_NERRORS!=15) return "No idea, implementation is not up-to-date.";
	if(index>=ERRSOURCE_MXM_EL_UNCORR && index<ERRSOURCE_MXM_EL_CORR)
		return (index&2) ? 
			"Bin-dependent: Electron fake rate (with B) stat+uncorrelated syst" :
			"Bin-dependent: Electron fake rate (no B) stat+uncorrelated syst";
	else if(index==ERRSOURCE_MXM_EL_CORR)
		return "Overall: Electron fake rate correlated syst (MC subtr) + Electron real efficiency";
	else if(index>=ERRSOURCE_MXM_MU_UNCORR && index<ERRSOURCE_MXM_MU_CORR)
		return "Bin-dependent: Muon fake rate stat+uncorrelated syst";
	//else if(index==ERRSOURCE_MXM_MU_UNCORR_LAST)
	//	return "Last pT bin: Muon fake rate extrapolation to higher pT";
	else if(index==ERRSOURCE_MXM_MU_CORR)
		return "Overall: Muon fake rate correlated syst (MC subtr) + Muon real efficiency";
	else return "ERROR: not recognized";
}

const char* TwoLepSSFakeBkg::GetErrorSourceDescription_OStoSS(unsigned int index)
{
	if(ERRSOURCE_MISID_NERRORS!=1) return "No idea, implementation is not up-to-date.";
	if(index==0)
		return "All combined";
	else return "ERROR: not recognized";
}