/*****************************************************************************/
/*                                                                           */
/* File Name        : FakeLepMCTemplate.h                                    */
/* Author           : Othmane Rifki			                                 */
/* Email            : othmane.rifki@cern.ch			                         */
/* Description      : Header file for the FakeLepMCTemplate class            */
/*                                                                           */
/***** C 2016 ****************************************************************/

#ifndef FakeLepMCTemplate_h
#define FakeLepMCTemplate_h

#include <iostream> 
// These are standard C++ header files.
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <list>
#include <vector>
#include <map>
#include <sys/types.h>
#include <sys/stat.h>

// These are ROOT header files.
#include "TMath.h"
#include "Rtypes.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "THStack.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TRandom.h"

using namespace std;

class FakeLepMCTemplate
{
 public:   

  //caw
  vector< TH1D *> data;
  map< int, vector< TH1D *> > bkg;  // < ifile, CR vector >
  vector< TH1D *> hsum;
  int nfiles;
  int ncr;
  int nsf;

  FakeLepMCTemplate();               
  int DoFit(double *corr, double *err);
  int Initialize(const std::string path, const std::string *files, const std::string *CRnames, const int NumFiles, const int NumCRs, const int NumSF);

 private:
  double myFunction(double par[]);
  static void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg);

  Double_t CalLikelihood(TH1D* MC, TH1D* data);

};

#endif
