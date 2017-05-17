/**
 *
 * @file FakeLepMxM/example_dMxM.cxx
 * @authors Julien Maurer <jmaurer@cern.ch>
 * @date September, 2015
 * @ You can find below the pseudo-code for an example use of the dynamic (generalised) matrix method as implemented in the GeneralisedMatrixMethod package
 * 
 **/

/*
/*  There are three problematics to have in mind:
/*
/*  1) The dynamic (generalized) matrix method, as opposed to the usual fixed-size matrix. While the formalism is a bit different, 
/*     as it takes advantage of the particular structure of the matrix instead of doing a generic matrix inversion, 
/*     the results are rigorously identical between the two approaches. 
/*
/*
/*  2) The fact that the tool often returns several weights. This is not a specificity of the dynamic matrix method, 
/*     and should be done similarly for a fixed-size matrix - although in several cases, it can simplif away. 
/*     Take the example of a signal region requiring at least one signal lepton, and cutting on the transverse mass mt(computed with the leading signal lepton). 
/*     Events in the signal region might have exactly one signal lepton, or one signal lepton + one baseline-not-signal lepton, or two signal leptons, or...
/*     The fake lepton background estimate will therefore be: 
/*          N_fakes = N_fakes^{1 signal lepton} + N_fakes^{1 signal lepton + 1 baseline-not-signal lepton} + N_fakes^{2 signal leptons} + ...
/*     Let's now assume an event in input of the matrix method, containing one electron and one muon. 
/*     One can apply the matrix method to estimate either the term N_fakes^{1 signal electron + 1 baseline-not-signal muon}, 
/*     or the term N_fakes^{1 signal muon + 1 baseline-not-signal electron}, or the term N_fakes{1 signal electron + 1 signal muon}. 
/*     Obviously, all three terms are needed, which is why the tool will return three different weights, corresponding to each of these estimates; 
/*     the Result structure also encodes, for each weight, which leptons are to be considered as signal and which are not. 
/*     One could argue that this splitting of the weights could be hidden inside the tool, and make it return only one weight that would be the sum of the three.
/*     Unfortunately this is not possible in the general case, because the rest of the SR selection might depend on which leptons are signal or not. 
/*     For example, here, mT would be computed respectively with the electron for the first term, then the muon for the second term, 
/*     and whichever lepton has largest pT for the 3rd term; for that reason, the tool can't know in advance which combinations will lead 
/*     to events passing or not the SR requirements; this has to be done in the analysis code, as is illustrated below. 
/*     Other non-trivial examples include the computation of the effective mass, requirements on the charges of the leptons (for a same-sign selection...)
/*
/*
/*  3) The systematic uncertainties, which originate from the propagation of uncertainties on the measured fake rates and prompt lepton efficiencies. 
/*     The tool handles correlations between the leptons (kinematics, flavours) when they apply (e.g. both leptons using same fake-rate) and returns 
/*     a list of  uncertainties (see below, in Counter::Fill). These uncertainties are signed, so as to benefit from event-by-event cancellations 
/*     due to anti-correlations between the weights for events with tight/loose leptons. However, those cancellations can happen only if the 
/*     different nuisance parameters are kept separate (otherwise one needs to add them in quadrature at the event level, and the sign of the correlation is lost). 
/*     
/*/

#include "GeneralisedMatrixMethod/MatrixMethod.h"
#include <algorithm>
#include <cmath>
#include <iostream>

class Counter
{
	double yield;
	double sumw2; // sum of squared weights, for statistical uncertainty
	double* sumsys; // (signed) sums of one_sigma_up for each nuisance parameter
	int nsys; // number of nuisance parameters
public:	
	Counter() : yield(0.), sumw2(0.), sumsys(nullptr), nsys(0) {}
	~Counter() { delete[] sumsys; }
	void Fill(const GeneralisedMatrixMethod::MatrixMethod::Result& result);
	void Print() const;
};

void event_loop(Counter& counter);

int main(int argc,const char* argv[])
{
	Counter counter;
	event_loop(counter);
	counter.Print();
}

void event_loop(Counter& counter)
{
	GeneralisedMatrixMethod::MatrixMethod dyn_mxm;
	for(auto& event : events_pool)
	{
		// 1. List the baseline leptons in the event (and some of their kinematics) to be provided to the tool
		std::vector<GeneralisedMatrixMethod::MatrixMethod::Lepton MXMLepton> mxm_leptons;
		int n_baseline_leptons = event.leptons.size();
		for(auto& lepton : event.leptons)
		{
			mxm_leptons.emplace_back(
				lepton.pt_in_GeV(),
				lepton.eta(),
				helpers.deltaRmin(lepton,event.jets),
				lepton.is_electron(),
				lepton.is_signal()
			);
		}
		
		// 2. Call the tool and retrieve the list of weights and "promoted to tight" leptons
		std::vector<GeneralisedMatrixMethod::MatrixMethod::Result> mxm_results = dyn_mxm.weightsForInput(mxm_leptons, event.bjets.size()>0);
		
		// 3. Loop on all lepton combinations considered by the tool and fill SR yield counters. 
		//    Some of the event global variables (meff, mt) depend on which leptons are signal,
		//    so they need to be recomputed for each combination
		for(auto& result : mxm_results)
		{
			// 3.1 count the number of signal leptons in the current combination
			int n_signal_leptons = std::count(result.leptonTightList.begin(), result.leptonTightList.end(), true);
			if(!n_signal_leptons) continue;
			
			// 3.2 retrieve the leading (signal) lepton, use it to compute mt
			int lead_lepton_index = std::find(result.leptonTightList.begin(), result.leptonTightList.end(), true) - result.leptonTightList.begin();
			float mt = helpers.mt(event.lepton[lead_lepton_index], event.met);
			
			// 3.3 recompute meff with the signal leptons corresponding to this combination
			float meff = event.ht_jets + event.met;
			for(int i=0;i<event.leptons.size();++i)
			{
				if(result.leptonTightList[i]) meff += event.leptons[i].pt();
			}
			
			// 3.4 fill signal region "one lepton inclusive"
			if(n_signal_leptons>=1 && meff>1000. && mt>200.) counter.Fill(result);
		}
	}
}

void Counter::Fill(const GeneralisedMatrixMethod::MatrixMethod::Result& result)
{
	// 4. Initialize the first time the function is called
	if(!this->sumsys)
	{
		this->nsys = 3 + result.uncUncorr.size();
		this->sumsys = new double[this->nsys];
		std::fill_n(this->sumsys, this->nsys, 0.);
	}
	
	// 5. Increment counters (and uncertainties)
	this->yield += result.weight;
	this->sum_w2 += pow(result.weight,2); // sum of squared weights for the statistical uncertainty, 
										  // neglecting correlations between different combinations for a single event (effect shown to be negligible in the SS/3L)
	this->sumsys[0] += result.uncElCorr; // propagated global systematic uncertainty on the electron fake rates (for ex. extrapolation CR->SR, subtraction of prompt bkg)
	this->sumsys[1] += result.uncMuCorr; // propagated global systematic uncertainty on the muon fake rates (for ex. extrapolation CR->SR, subtraction of prompt bkg)
	this->sumsys[2] += result.uncEff; // propagated uncertainty on prompt lepton efficiencies
	for(int i=3;i<this->nsys;++i)
	{
		this->sumsys[i] += result.uncUncorr[i-3]; // propagated uncertainties on the lepton fake rates uncorrelated across bins (typically, the statistical uncertainties)
	}
	// if you don't want to track separately each nuisance parameter, you can sum them in quadrature at this stage, 
	// but this will lead to an overestimate of the uncertainties
}

void Counter.Print() const
{
	double sys2 = 0;
	for(int i=0;i<this->nsys;++i) sys2 += pow(sumsys[i],2); // quadratic sum of all systematic uncertainties, only done at the end of the event loop
	std::cout<<this->yield<<" +/- "<<sqrt(this->sumw2)<<" (stat) +/- "<<sqrt(sys2)<<" (syst)\n";
}