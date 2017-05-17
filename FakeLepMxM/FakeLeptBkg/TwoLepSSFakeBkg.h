/**
 * @file FakeLeptBkg/FakeLeptBkg/TwoLepSSFakeBkg.h
 * @authors Julien Maurer <jmaurer@cern.ch>   Otilia Ducu <oducu@cern.ch>
 * @date May, 2012
 * @brief Fake background estimation for same-sign di-lepton analysis: matrix method and charge mis-ID
 * 
 * 
 * IMPORTANT: All energies to be set in GeV!! 
 * 
 *  ------
 *  Usage:
 *  ------
 *   
 *  -- Background from charge mis-ID
 * 
 *      -- rescaling of opposite sign events with 2 true leptons
 *         from MC (Z, ttbar, dibosons...) or data by electron charge mis-ID rate
 *      
 *      -- use functions GetWeight_OStoSS_El*() to get event weight and syst error
 *   
 *      -- keep only events for which el_type==2 or 4 (for both electrons in ee channel)
 *
 *
 *  -- Background from fakes
 *
 *      -- matrix method with 2 leptons (4D)
 *
 *      -- use functions GetWeight_MxM_*() to get event weight and syst error
 * 
 *      -- needed arguments; lepton pT/eta/pass tight, presence of bjets in event
 * 
 *      -- however, background from charge mis-ID has to be subtracted to not interfere
 * 
 *      -- recipe to do this:
 * 
 *          --  matrix method linear wrt number of events
 * 
 *          --  apply also matrix method on MC events used for charge mis-ID 
 *
 *          --  add events to data MxM estimation with a weight (-1 * weight_misID * weight_mxm)
 *
 *  -- Splitted uncertainties
 * 
 *      --  Possible to retrieve different sources of uncorrelated uncertainties splitted, or everything combined
 *
 *      --  Toggled by UseSplitErrBySource_MxM(true/false, default=false)  -- currently, splitting only supported for matrix method, not charge flip
 *
 *      --  If active, pass an array of size TwoLepSSFakeBkg::nErrorSources_MxM to function GetWeight_MxM_*()
 *
 *      --  Meaning of different sources can be retrieved by a call to GetErrorSourceDescription_MxM
 *
 *  -- Additional comments:
 * 
 *     --  Not supporting d0/sig < 5 (electrons) anymore, removed argument from the functions
 *
 *     -- Assuming that background events have always at least one true lepton (W/Z+jet, ttbar...),
 *        matrix method can be simplified to use only 3 samples: N[pass,pass], N[pass,fail], N[fail,pass]
 *
 *         -- allows to use single lepton trigger with tighter cuts than the ones defining "loose" leptons
 *
 *         -- can be set by Use3DMatrixMethod(true)
 *
 *         -- no support of splitted systematics for 3D matrix though in the current version!
 */


#ifndef FAKELEPTBKG_TWOLEPSSFAKEBKG_H
#define FAKELEPTBKG_TWOLEPSSFAKEBKG_H


class TwoLepSSFakeBkg
{
protected:
	bool m_SplitErrBySource_MxM, m_SplitErrBySource_OStoSS, m_UseMxM3D;
public:
	static const unsigned int nErrorSources_MxM;
	static const unsigned int nErrorSources_MxM_El;
	static const unsigned int nErrorSources_MxM_Mu;
	static const unsigned int nErrorSources_OStoSS;
	
	TwoLepSSFakeBkg(void);
	~TwoLepSSFakeBkg();

	// to set options
	void UseSplitErrBySource_MxM(bool split);	
	void UseSplitErrBySource_OStoSS(bool split);		
	void Use3DMatrixMethod(bool useMxM3D);

	// to rescale OS events to SS (electron charge flip)
	double GetWeight_OStoSS_ElMu(double elEta,double elPt,bool elIsTight,double *wgt_err) const;
	double GetWeight_OStoSS_ElEl(double el1Eta,double el1Pt,bool el1IsTight,double el2Eta,double el2Pt,bool el2IsTight,double *wgt_err) const;

	// to apply matrix method
	double GetWeight_MxM_ElEl(double el1Eta,double el1Pt,bool el1IsTight,double el2Eta,double el2Pt,bool el2IsTight,bool hasBJets,double *wgt_err) const;
	double GetWeight_MxM_ElMu(double elEta,double elPt,bool elIsTight,double muEta,double muPt,bool muIsTight,bool hasBJets,double *wgt_err) const;
	double GetWeight_MxM_MuMu(double mu1Eta,double mu1Pt,bool mu1IsTight,double mu2Eta,double mu2Pt,bool mu2IsTight,bool hasBJets,double *wgt_err) const;

	static const char* GetErrorSourceDescription_MxM(unsigned int index);
	static const char* GetErrorSourceDescription_OStoSS(unsigned int index);  
};

#endif
