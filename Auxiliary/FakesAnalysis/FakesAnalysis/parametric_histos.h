#ifndef parametric_histos_h
#define parametric_histos_h

#include <iostream>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include <TFile.h>

class parametric_histos {

 public:
  //costruttore
  parametric_histos();
  parametric_histos(float* xbins, int nxbins, std::vector<double>& x);
  parametric_histos(float* xbins, int nxbins, std::vector<double>& x, std::vector<double>& weight);
  parametric_histos(float* xbins, float* ybins, int nxbins, int nybins, std::vector< std::vector< double> >& val);
  parametric_histos(std::vector< std::vector<double> > bins, TString  histos_name );
  //distruttore
  ~parametric_histos();
  //metodi
   void Write();
   void Write2(const char * name_file);
   int GetxBinNumber(double x);
   void AddBinContent(int bin);
   void AddBinContent(int bin, double weight);
   std::vector<TH1F*> Get1Dhistos();
   TH1F *Get1Dhisto(int i){ return histos[i]; }
   TH1F* GetGeneralHisto(){ return h_all; }
   void Fill_values(int par_number, double value);
   void Fill_weights(double weight);
   void Execute();
   int GetNhistos(){return m_npar;}
 private:

  std::vector<double> m_x;
  std::vector<double> m_W;
  std::vector<double> m_y;

  std::vector<std::vector<double> > m_val;
  std::vector<std::vector<double> > m_bins;
  std::vector<std::vector<double> > m_binPos;
  std::vector<int> m_nbins;
  std::vector<TH1F*> histos;
  
  float *m_xbins;
  float *m_ybins;
  int m_nentries;
  int m_nxbins;
  int m_nybins;
  int m_npar;
  TH1F  *h_x;
  TH1F  *h_y;
  TH1F  *h_xy;
  TH2F *h_xy2D;
  TH1F *h_all;
  TString m_histosName;


};

#endif // __parametric_histos_h__
