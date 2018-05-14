/*****************************************************************************/
/*                                                                           */
/* File Name        : FakeLepMCTemplate.cxx                                  */
/* Author           : Othmane Rifki			                                 */
/* Email            : othmane.rifki@cern.ch			                         */
/* Description      : Source file for the FakeLepMCTemplate class            */
/*                                                                           */
/***** C 2016 ****************************************************************/

#include "FakeLepMCTemplate/FakeLepMCTemplate.h"

static FakeLepMCTemplate *LKF_obj;

int FakeLepMCTemplate::Initialize(const std::string path, const std::string *files, const std::string *CRnames, const int NumFiles, const int NumCRs, const int NumSFs){

  nfiles = NumFiles;
  ncr = NumCRs;
  nsf = NumSFs;

  // pass the histograms
  TFile *f_data;
  char temp_file[512];

  for( int i = 0; i < NumFiles ; i++ ){
	std::sprintf(temp_file, "%s/%s.root", path.c_str(),files[i].c_str());
    f_data = new TFile( temp_file, "READ" );
    if( !f_data->IsOpen() )throw "input file not opened";

    for( long j = 0; j < NumCRs; j++ ){
      TH1D * t = (TH1D*) f_data->Get( CRnames[j].data() )->Clone();
      if( t == 0 )throw "missing CR histogram";
      t->SetDirectory( 0 );
      if( i == 0 ){ //first file is data
        data.push_back( t );
        TH1D * ts = (TH1D *) t->Clone();
        ts->SetDirectory( 0 );
        hsum.push_back( ts );
      }
      else{
        bkg[i].push_back( t );
      }
    }
    f_data->Close();
  }

  return 0;
}

int FakeLepMCTemplate::DoFit(double *corr, double *err){

  // Initializing the variables
  const int nfit = nsf;

  TFitter* minimizer = new TFitter(nfit);
  minimizer->SetFCN(FakeLepMCTemplate::minuitFunction);
  double arglist[10];

  //Setting up minuit parameters
  Double_t *vstart = new Double_t[nfit];
  Double_t *step = new Double_t[nfit];
  TString *name = new TString[nfit];

  for(int n=0; n<nfit; n++){
	vstart[n] = 1.;
	step[n] = 0.01;
	name[n] = Form("template%d",n+1);
  }

  for (int j = 0; j < nfit; j++) {
	minimizer->SetParameter(j, name[j], vstart[j], step[j], 0., 10.);}

  minimizer->ExecuteCommand("CAL", arglist, 1);
  minimizer->ExecuteCommand("SET NOG", arglist, 0);
  minimizer->ExecuteCommand("MINIMIZE", arglist, 0);
  minimizer->ExecuteCommand("MIGRAD", arglist, 0);
  minimizer->ExecuteCommand("HESSE", arglist, 0);
  minimizer->ExecuteCommand("CAL", arglist, 1);


  for (int j = 0; j < nfit; j++) {
	corr[j] = minimizer->GetParameter(j);
	err[j] = minimizer->GetParError(j);
  }

  delete[] vstart;
  delete[] step;
  delete[] name;

  return 0;
}

double FakeLepMCTemplate::myFunction(double par[]){

  double sf;

  // Initialize the histograms with data and reset to 0
  for( auto &h : hsum ){
    h->Reset();
  }
  // Reweight distributions with scale factors
  //-------------------------------------------------------------------------------------------------

  for( int i = 1; i < nfiles; i++ ) {
	sf = 1.0;

	for(int n = 0; n < nsf; n++){
	  if (i == 2+n)
		sf = fabs(par[n]);
	}

	for( int j = 0; j < ncr; j++ ){
	  hsum.at( j )->Add( bkg[i][j], sf);
	}

  }

  //-------------------------------------------------------------------------------------------------

  // Now we calculate the probability
  double result = 0; // chi2

  for( int i = 0; i < (int)data.size(); i++ ){
    result += CalLikelihood( hsum[i], data[i] ); // CalLikelihood(MC, data)
  }

  return result;
}

void FakeLepMCTemplate::minuitFunction(int& nDim, double* gout, double& result, double par[], int flg){
  result = LKF_obj->myFunction(par);
  if(nDim || gout || flg){} // to avoid warning
}

FakeLepMCTemplate::FakeLepMCTemplate()
{
  LKF_obj = this;
}

//-----------------------------------------------------------------------------
// Likelihood calculation
//-----------------------------------------------------------------------------
Double_t FakeLepMCTemplate::CalLikelihood(TH1D* MC, TH1D* data) {

  double f = 0;

  for (int i = 1; i <= MC->GetNbinsX(); i++) {
	if (MC->GetBinContent(i) > 0) {
	  //Uncertainty on the predicted number of events
	  double x;
	  double loc_prob = 0;
	  double gaus;
	  double gaus_norm = 0;

	  //This is integral is within 2 sigma.
	  double max, min;
	  max = MC->GetBinContent(i) + 2. * MC->GetBinError(i);
	  min = MC->GetBinContent(i) - 2. * MC->GetBinError(i);
	  if (min < 0)
		min = 0;

	  //we do 30 iterrations for the integral
	  for (int j = 0; j < 30; j++) {
		x = min + (max - min) * ((double) j + 0.5) / 30.;
		gaus = TMath::Gaus(x, MC->GetBinContent(i), MC->GetBinError(i), kFALSE);
		gaus_norm += gaus;
		loc_prob += TMath::Poisson(data->GetBinContent(i), x) * gaus;
	  }

	  loc_prob = loc_prob / gaus_norm;

	  f += -2*log(loc_prob);
	}
  }

  return f;
}
