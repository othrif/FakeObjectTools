#define parametric_histos_cxx
#include "FakesAnalysis/parametric_histos.h"
#include <sstream>
#include<string>
#include "TString.h"
//costruttori:
parametric_histos::parametric_histos(){} //default

parametric_histos::parametric_histos(float* xbins, int nxbins,  std::vector<double>& x){ //1 parameter


  m_x = x;
  m_xbins = new float[nxbins];
  m_xbins = xbins;
  m_nentries = x.size();
  m_nxbins = nxbins;
  std::cout<< m_nxbins << std::endl;
  h_x = new TH1F("x_hist", "",m_nxbins ,xbins);
  for( int i=0; i< m_nentries; i++){
    h_x->Fill(m_x[i]/1000);
  }

}

parametric_histos::parametric_histos(float* xbins, int nxbins,  std::vector<double>& x, std::vector<double>& weight ){ //1parameter + weight

  m_x = x;
  m_W = weight;
  m_xbins = new float[nxbins];
  m_xbins = xbins;
  m_nentries = x.size();
  m_nxbins = nxbins;
  std::cout<< m_nxbins << std::endl;
  h_x = new TH1F("x_hist", "",m_nxbins ,xbins);
  for( int i=0; i< m_nentries; i++){
    h_x->Fill(m_x[i]/1000, m_W[i]);
  }

}
parametric_histos::parametric_histos(float* xbins, float* ybins,  int nxbins,  int nybins, std::vector< std::vector<double> >& val){ //2 parameters

  m_x = val[0];
  m_y = val[1];
  m_xbins = new float[nxbins];
  m_ybins = new float[nybins];
  m_xbins = xbins;
  m_ybins = ybins;
  m_nentries = m_x.size(); // missing a check for y dimension
  m_nxbins = nxbins;
  std::cout << "nxbin "<<m_nxbins << std::endl;
  m_nybins = nybins;
  std::cout << "nybin "<<m_nybins << std::endl;
  h_x = new TH1F("x_hist", "",m_nxbins ,xbins);
  h_y = new TH1F("y_hist", "",m_nybins ,ybins);
  h_xy = new TH1F("xy_hist", "",m_nybins*m_nxbins, -0.5, m_nybins*m_nxbins - 0.5);
  h_xy2D = new TH2F("xy2D_hist","",m_nxbins ,xbins,m_nybins ,ybins);
  TAxis* xAxis = h_x->GetXaxis();
  TAxis* yAxis = h_y->GetXaxis();
  for( int i=0; i< m_nentries; i++){
    h_x->Fill(m_x[i]);
    h_y->Fill(m_y[i]);
    h_xy2D->Fill(m_x[i],m_y[i]);
     int ibinX = xAxis->FindBin(m_x[i]) - 1;
    // std::cout << "pt "<< m_x[i] << std::endl;
    // std::cout << "x bin "<< ibinX << std::endl;
    int ibinY = yAxis->FindBin(m_y[i]) - 1;
    std::cout << "eta "<< m_y[i] << std::endl;
    std::cout << "y bin " << ibinY<< std::endl;
    int ibin = ibinX + ibinY*m_nxbins;
    std::cout << ibin << std::endl;
    h_xy->Fill(ibin);
  
  }
}
parametric_histos::parametric_histos(std::vector< std::vector<double> > bins, TString histo_name){ //general

  m_npar = bins.size(); 
  std::cout<<"n par " << m_npar <<std::endl;
  m_bins = bins;
  m_histosName = histo_name;
  for(int i =0; i < m_npar; i++) {
    std::vector <double> par_vec;
    m_val.push_back(par_vec);
    m_nbins.push_back(bins[i].size()-1);
  }

}


parametric_histos::~parametric_histos(){}



void parametric_histos::Fill_values(int par_number, double value){ 
  m_val[par_number].push_back(value);
}

void parametric_histos::Fill_weights(double weight){ 
  m_W.push_back(weight);
}

void parametric_histos::Execute(){
  // std::cout << "enter execute" << std::endl;
  const unsigned int vectorsize = m_val[0].size(); //here you need to make a controll
 
  for(unsigned int i=0; i<m_val.size(); ++i) {
    if(m_val[i].size() != vectorsize){
      std::cout << "vectors with different size!" << std::endl;
      return;
    }
  }
  m_nentries = m_val[0].size(); 
  // std::cout << "decide size" << std::endl;
  for( int i=0; i<m_npar; i++){
    TString name = "histo"+m_histosName;
    name +=i;
    // std::cout << name << std::endl;
    //std::cout << m_nbins[i] << std::endl;
    double *bins_vec = new double[m_nbins[i]+1];
    for(int j=0; j < m_nbins[i]+1; j++){ 
      bins_vec[j] = m_bins[i][j];
      //std::cout<< bins_vec[j] << std::endl;
    }
    TH1F *histo = new TH1F(name,"", m_nbins[i], bins_vec);
    histo->Sumw2();
    TAxis* axis = histo->GetXaxis();
    std::vector<double> bin_position;
    //std::cout<<"histo creato"<< std::endl;
    for (int j = 0; j < m_nentries; j++){
      //	std::cout << m_val[i][j] << std::endl;
      //	std::cout << m_W[j] << std::endl;
      histo->Fill(m_val[i][j], m_W[j]);
      int ibin = axis->FindBin(m_val[i][j]) - 1;
      bin_position.push_back(ibin);
    }
   
    histos.push_back(histo);
    m_binPos.push_back(bin_position);
   
      
  }
   
  int totNbin = 1;
  for(int i = 0; i < m_npar; i++){
    totNbin = totNbin * (m_nbins[i]);
  }
  
  h_all = new TH1F("all_histo"+m_histosName, "",totNbin, -0.5, totNbin - 0.5);
  h_all->Sumw2();
 
  for (int i = 0; i < m_nentries; i++){
    int ibin = 0;
    for( int j = 0; j < m_npar; j++){
      int jbin=m_binPos[j][i];
      for( int k = 0; k< j; k++){
	jbin = jbin*(m_nbins[k]);
      }
      ibin = ibin +jbin;
    }
    h_all->Fill(ibin, m_W[i]);
  }

}


void parametric_histos::Write(){ 
 TFile* output = new TFile("prova.root","recreate");
 h_x->Write();
 h_y->Write();
 h_xy->Write();
 h_xy2D->Write();
 output->Close();
}

void parametric_histos::Write2(const char * name_file){ 
 TFile* output = new TFile(name_file,"UPDATE");
 for (int i =0; i < m_npar; i++) histos[i]->Write();
 h_all->Write();
 output->Close();
}

int parametric_histos::GetxBinNumber(double x){
  int i=0;
  int N =1;
  while(x/1000>m_xbins[i]){
    N = i+1;
    i++;
  }
  return N;
}

void parametric_histos::AddBinContent(int bin){
  h_x->AddBinContent(bin);
}

void parametric_histos::AddBinContent(int bin, double weight){
  h_x->AddBinContent(bin, weight);
}

std::vector<TH1F*> parametric_histos::Get1Dhistos(){
  return histos;
}
