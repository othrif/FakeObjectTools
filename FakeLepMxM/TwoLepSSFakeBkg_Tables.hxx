/**
 * @file TwoLepSSFakeBkg_Tables.hxx
 * @author Julien Maurer <jmaurer@cern.ch>  Otilia Ducu <oducu@cern.ch>
 * @date October 2016
 * @brief Fake background estimation for same-sign di-lepton analysis: tables of efficiencies used by the tool
 * @Objects used to perform the measurements : see section 4 of https://cds.cern.ch/record/2151944/files/ATL-COM-PHYS-2016-495.pdf
 **/

#pragma once

struct Bin2D {float minPt,maxPt,minEta,maxEta;}; 

const double Params_RealEffEl_Eff[] = {
0.638774, 0.757052, 0.68449, 
0.710553, 0.779129, 0.753873, 
0.775979, 0.826167, 0.810887, 
0.832765, 0.873852, 0.859228, 
0.872636, 0.897416, 0.878121, 
0.903054, 0.913886, 0.899835, 
0.933192, 0.937673, 0.927504, 
0.952533, 0.956723, 0.943655, 
0.968721, 0.972381, 0.959945, 
0.981162, 0.979692, 0.970789
};
const double Params_RealEffEl_Stat[] = {
0.0053, 0.0043, 0.0047, 
0.0028, 0.0027, 0.003, 
0.0016, 0.0017, 0.002, 
0.00097, 0.001, 0.0013, 
0.00064, 0.00076, 0.001, 
0.00045, 0.00058, 0.00082, 
0.00027, 0.00035, 0.00053, 
0.00047, 0.0006, 0.00097, 
0.00057, 0.00072, 0.0013, 
0.00078, 0.0011, 0.0019
};
const double Params_RealEffEl_SystUncorr[] = {
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.,
0., 0., 0.
};
const double Params_RealEffEl_SystCorr_relative[] = {
0.052583, 0.0522614, 0.0665238, 
0.0303517, 0.0353324, 0.0451503, 
0.0470247, 0.0484841, 0.0523718, 
0.042914, 0.043442, 0.0468909, 
0.0220235, 0.0240769, 0.025946, 
0.0206512, 0.0216172, 0.0225401, 
0.0200614, 0.0201566, 0.0203817, 
0.0101897, 0.0100856, 0.010287, 
0.00507717, 0.00507976, 0.00634429, 
0.0052182, 0.00525583, 0.00579173
};
const double Params_RealEffEl_SystBusy_relative[] = { // JO: parametrization is different (pT,dR) instead of (pT,eta)! hard-coded, for the moment
	0.08, 0.04, // 10-15 GeV
	0.08, 0.04, // 15-20 GeV
	0.08, 0.04, // 20-25 GeV
	0.08, 0.04, // 25-30 GeV
	0.08, 0.04, // 30-35 GeV
	0.08, 0.04, // 35-40 GeV
	0.08, 0.04, // 40-50 GeV
	0.08,0.04, // 50-60 GeV
	0.05,0.05, // 60-80 GeV
	0.05, 0.05 // >80 GeV
};

//const float Params_RealEffEl_EtaBins[] = {0, 0.8, 1.37, 1.52, 2};
const float Params_RealEffEl_EtaBins[] = {0, 0.8, 1.45, 2};
const float Params_RealEffEl_PtBins[] = {10.,15.,20.,25.,30.,35.,40.,50.,60.,80.,7000.};
const float Params_RealEffEl_dRBins[] = {0.,0.6,1e12};
const unsigned int Params_RealEffEl_nEtaBins = sizeof(Params_RealEffEl_EtaBins)/sizeof(float)-1;
const unsigned int Params_RealEffEl_nPtBins = sizeof(Params_RealEffEl_PtBins)/sizeof(float)-1;
const unsigned int Params_RealEffEl_ndRBins = sizeof(Params_RealEffEl_dRBins)/sizeof(float)-1;

const double Params_RealEffMu_Eff[] = {
0.828302, 0.830669, 0.825705, 0.816492, 
0.849776, 0.849229, 0.844155, 0.840655, 
0.894633, 0.895476, 0.894483, 0.891838, 
0.927425, 0.926303, 0.927246, 0.924489, 
0.945173, 0.942838, 0.942075, 0.939322, 
0.964994, 0.964842, 0.962678, 0.958919, 
0.980801, 0.981705, 0.981806, 0.97885, 
0.987251, 0.987747, 0.988328, 0.986063, 
0.989665, 0.990581, 0.990799, 0.988181, 
0.991588, 0.992368, 0.990901, 0.989205, 
0.991853, 0.992388, 0.991963, 0.990263
};
const double Params_RealEffMu_Stat[]={
0.0034, 0.0028, 0.0029, 0.0029, 
0.0018, 0.0019, 0.0021, 0.0021, 
0.001, 0.0011, 0.0012, 0.0012, 
0.00059, 0.00068, 0.00073, 0.00075, 
0.00041, 0.00046, 0.00053, 0.00055, 
0.00029, 0.0003, 0.00035, 0.00039, 
0.00016, 0.00016, 0.00016, 0.00019, 
0.00028, 0.00027, 0.00026, 0.00032, 
0.00045, 0.00042, 0.00043, 0.00054, 
0.00063, 0.00059, 0.00066, 0.00081, 
0.00059, 0.00057, 0.00061, 0.00078
};
const double Params_RealEffMu_SystUncorr[] = {
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.,
0., 0., 0., 0.
};
const double Params_RealEffMu_SystCorr_relative[]={
0.01, 0.01, 0.01, 0.01, 
0.005, 0.005, 0.005, 0.005, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001, 
0.001, 0.001, 0.001, 0.001
};
const double Params_RealEffMu_SystBusy_relative[] = { // JO: parametrization is different (pT,dR) instead of (pT,eta)! hard-coded, for the moment
	0.4, 0.1, // 10-15 GeV
	0.2, 0.07, // 15-20 GeV
	0.2, 0.07, // 20-25 GeV
	0.2, 0.07, // 25-30 GeV
	0.2, 0.07, // 30-35 GeV
	0.1, 0.05, // 35-40 GeV
	0.1, 0.05, // 40-50 GeV
	0.05,0.03, // 50-60 GeV
	0.05,0.03, // 60-70 GeV
	0.05,0.03, // 70-80 GeV
	0.01, 0.01 // >80 GeV
};

const float Params_RealEffMu_EtaBins[] = {0, 0.6, 1.2, 1.8, 2.5};
const float Params_RealEffMu_PtBins[] ={10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,80.,7000.}; 
const float Params_RealEffMu_dRBins[] ={0.,0.6,1e12}; 
const unsigned int Params_RealEffMu_nEtaBins = sizeof(Params_RealEffMu_EtaBins)/sizeof(float)-1;
const unsigned int Params_RealEffMu_nPtBins = sizeof(Params_RealEffMu_PtBins)/sizeof(float)-1;
const unsigned int Params_RealEffMu_ndRBins = sizeof(Params_RealEffMu_dRBins)/sizeof(float)-1;


const double Params_FakeRateEl_Eff[] = {0.1100, 0.0829, 0.1087, 0.2943, 0.2505}; 
const double Params_FakeRateEl_Stat[] = {0.0118, 0.0176, 0.0255, 0.0574, 0.0850}; 
const double Params_FakeRateEl_SystUncorr[] = {0, 0, 0, 0, 0};
const double Params_FakeRateEl_SystCorr[] = {0.0550, 0.0416, 0.0560, 0.1487, 0.1366}; 
const Bin2D Params_FakeRateEl_Bins[] = {{10.,15.,0.,2.5},{15.,20.,0.,2.5},{20.,30.,0.,2.5},{30.,40.,0.,2.5},{40.,7000.,0.,2.5}};
const unsigned int Params_FakeRateEl_nBins=(sizeof(Params_FakeRateEl_Bins)/sizeof(Bin2D));

const double Params_FakeRateMu_Eff[] = {0.1487, 0.1120, 0.1040, 0.1535, 0.1715, 0.1755}; 
const double Params_FakeRateMu_Stat[] = {0.0169, 0.0225, 0.0247, 0.0545, 0.0750, 0.0504}; 
const double Params_FakeRateMu_SystUncorr[] = {0., 0., 0., 0., 0., 0.};  
const double Params_FakeRateMu_SystCorr[] = {0.0744, 0.0561, 0.0524, 0.0790, 0.0907, 0.0890}; 

const Bin2D Params_FakeRateMu_Bins[] = {{10.,15.,0.,2.5},{15.,20.,0.,2.5},{20.,30.,0.,2.5},{30.,40.,0.,2.5},{40.,51.,0.,2.5},{51.,7000.,0.,2.5}};
const unsigned int Params_FakeRateMu_nBins=(sizeof(Params_FakeRateMu_Bins)/sizeof(Bin2D));

double Params_ChargeMisID_Tight_Rate[] = {
/// Data 9.5 ifb
0.00121, 0.00095, 0.00356, 
0.00090, 0.00107, 0.00533, 
0.00072, 0.00194, 0.00807, 
0.00046, 0.00164, 0.00788, 
0.00044, 0.00189, 0.01013, 
0.00072, 0.00350, 0.01451, 
0.00133, 0.00387, 0.02124, 
0.00106, 0.00380, 0.01882, 
0.00210, 0.00415, 0.02632, 
0.00243, 0.00708, 0.03641
};
double Params_ChargeMisID_Tight_Stat[] = {
/// Data 9.5 ifb
0.00059, 0.00035, 0.00050, 
0.00012, 0.00009, 0.00021, 
0.00005, 0.00008, 0.00018, 
0.00002, 0.00006, 0.00017, 
0.00007, 0.00015, 0.00041, 
0.00014, 0.00031, 0.00082, 
0.00023, 0.00047, 0.00139, 
0.00030, 0.00064, 0.00183, 
0.00048, 0.00088, 0.00277, 
0.00034, 0.00074, 0.00232
};
double Params_ChargeMisID_Tight_SystUncorr[] = {
/// Hubert : data 5.8 ifb
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.
};
double Params_ChargeMisID_Tight_SystCorr[] = {
/// Data 9.5 ifb
0.00121*0.67215, 0.00095*0.50298, 0.00356*0.48377, 
0.00090*0.20954, 0.00107*0.11182, 0.00533*0.09796, 
0.00072*0.06456, 0.00194*0.03373, 0.00807*0.02385, 
0.00046*0.00848, 0.00164*0.01358, 0.00788*0.01413, 
0.00044*0.04779, 0.00189*0.02931, 0.01013*0.05647, 
0.00072*0.21242, 0.00350*0.02504, 0.01451*0.02730, 
0.00133*0.12567, 0.00387*0.11202, 0.02124*0.02284, 
0.00106*0.14416, 0.00380*0.04946, 0.01882*0.02919, 
0.00210*0.04218, 0.00415*0.02926, 0.02632*0.07026, 
0.00243*0.02441, 0.00708*0.04875, 0.03641*0.02202
};

/// Hubert:
const float Params_ChargeMisID_Tight_EtaBins[] = {0.0, 0.8, 1.37, 2.0};
const float Params_ChargeMisID_Tight_PtBins[] ={10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 7000};

const unsigned int Params_ChargeMisID_Tight_nEtaBins = sizeof(Params_ChargeMisID_Tight_EtaBins)/sizeof(float)-1;
const unsigned int Params_ChargeMisID_Tight_nPtBins = sizeof(Params_ChargeMisID_Tight_PtBins)/sizeof(float)-1;


double Params_ChargeMisID_Loose_Rate[] = {
/// Data 5.8 ifb
0.00439, 0.00787, 0.01895, 
0.00498, 0.00804, 0.02651, 
0.00633, 0.01434, 0.04757, 
0.00821, 0.01712, 0.05790, 
0.01223, 0.02309, 0.08458, 
0.02050, 0.03326, 0.09072, 
0.01890, 0.04367, 0.10730, 
0.02370, 0.06395, 0.12821, 
0.04791, 0.06695, 0.11017, 
0.02703, 0.06276, 0.10322
};
double Params_ChargeMisID_Loose_Stat[] = {
/// Data 5.8 ifb
0.00112, 0.00121, 0.00183, 
0.00033, 0.00045, 0.00080, 
0.00025, 0.00046, 0.00092, 
0.00036, 0.00056, 0.00135, 
0.00091, 0.00148, 0.00398, 
0.00224, 0.00406, 0.00655, 
0.00378, 0.00757, 0.01098, 
0.00642, 0.01290, 0.02441, 
0.01292, 0.02005, 0.02003, 
0.00960, 0.01938, 0.02294
};
double Params_ChargeMisID_Loose_SystUncorr[] = {
/// Hubert : data 5.8 ifb
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.,
	0, 0., 0.
};
double Params_ChargeMisID_Loose_SystCorr[] = {
/// Data 5.8 ifb
0.00439*2.11808, 0.00787*1.07742, 0.01895*0.55493, 
0.00498*0.38913, 0.00804*0.27304, 0.02651*0.15252, 
0.00633*0.06633, 0.01434*0.05300, 0.04757*0.03716, 
0.00821*0.01472, 0.01712*0.01541, 0.05790*0.01773, 
0.01223*0.03143, 0.02309*0.03724, 0.08458*0.08316, 
0.02050*0.05942, 0.03326*0.03620, 0.09072*0.07947, 
0.01890*0.13074, 0.04367*0.03113, 0.10730*0.04686, 
0.02370*0.08522, 0.06395*0.12885, 0.12821*0.05011, 
0.04791*0.07514, 0.06695*0.08125, 0.11017*0.08856, 
0.02703*0.16673, 0.06276*0.04936, 0.10322*0.02315
};

/// Hubert:
const float Params_ChargeMisID_Loose_EtaBins[] = {0.0, 0.8, 1.37, 2.0};
const float Params_ChargeMisID_Loose_PtBins[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 7000};

const unsigned int Params_ChargeMisID_Loose_nEtaBins = sizeof(Params_ChargeMisID_Loose_EtaBins)/sizeof(float)-1;
const unsigned int Params_ChargeMisID_Loose_nPtBins = sizeof(Params_ChargeMisID_Loose_PtBins)/sizeof(float)-1;
