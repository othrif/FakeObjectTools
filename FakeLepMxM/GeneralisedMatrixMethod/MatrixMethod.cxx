/**
 * @file MatrixMethod.cxx
 * @author Thomas Gillam <thomas.gillam@cern.ch>
 * @date May, 2013
 * @brief Estimate fake components for events for arbitrary number of leptons
 **/ 

#include "GeneralisedMatrixMethod/MatrixMethod.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// extern bool debug_me;

namespace GeneralisedMatrixMethod {

  std::vector<MatrixMethod::Result> MatrixMethod::weightsForInput(const std::vector<Lepton> &leptons, int nBJetsInEvent)
  {
    this->leptons = &leptons;
    setNumberOfLeptons();

    if (!haveEnoughLeptons()) {
      return std::vector<Result>();
    }

    allocateMemory();
    fillLeptonInformation(nBJetsInEvent);

    std::vector<Result> results;

    zeroConfiguration(configTLOut);
    while (!isLastConfiguration(configTLOut)) {
      if (!isWantedTLOutConfiguration()) {
        incrementConfiguration(configTLOut);
        continue;
      }

      results.push_back(resultForCurrentTLOutConfiguration());

      incrementConfiguration(configTLOut);
    }

    clearMemory();
    
    return results;
  }

  void MatrixMethod::setNumberOfLeptons()
  {
    numLeptons = leptons->size();

    // Truncate the problem to save runtime
    if (numLeptons > 12) numLeptons = 12;
  }

  bool MatrixMethod::haveEnoughLeptons() const
  {
    return (numLeptons >= 2);
  }

  void MatrixMethod::allocateMemory()
  {
    isLepElectron = new bool [numLeptons];
    effReal = new float [numLeptons];
    effFake = new float [numLeptons];
    uEffReal = new float [numLeptons];
    uEffFakeStat = new float [numLeptons];
    uEffFakeUncorr = new float [numLeptons];
    uEffFakeCorr = new float [numLeptons];
    uEffFakeBinIndex = new unsigned int [numLeptons];
    configTLIn = new char [numLeptons];
    configTLOut = new char [numLeptons];
    configRF = new char [numLeptons];
    derivativesReal = new float [numLeptons];
    derivativesFake = new float [numLeptons];
  }

  void MatrixMethod::fillLeptonInformation(int nBJetsInEvent)
  {
    for (unsigned int i = 0; i < numLeptons; ++i) {
      float eta = (*leptons)[i].eta;
      float pt = (*leptons)[i].pt;
      float drjet = (*leptons)[i].drjet;
      bool isElectron = (*leptons)[i].isElectron;
      bool isTight = (*leptons)[i].isTight;

      OJEfficiencyLoader::Efficiency real, fake;
      
      //unsigned int binFlavourOffset = 0; // JO: changing bin indexing scheme
      //unsigned int binBJetOffset = (bJetsInEvent ? 1 : 0); // JO: changing bin indexing scheme
      //unsigned int rawBinScaleFactor = 1; // JO: changing bin indexing scheme
      if (isElectron) {
        isLepElectron[i] = true;
        //rawBinScaleFactor = 2; // JO: changing bin indexing scheme
        real = ojEfficiencyLoader.realEfficiencyEl(eta, pt, drjet);
        fake = ojEfficiencyLoader.fakeRateEl(eta, pt, drjet, nBJetsInEvent);
      } else {
        isLepElectron[i] = false;
        //binFlavourOffset = 2*ojEfficiencyLoader.numFakeBinsEl(); // JO: changing bin indexing scheme
        //binBJetOffset = 0; // JO: changing bin indexing scheme
        real = ojEfficiencyLoader.realEfficiencyMu(eta, pt, drjet);
        fake = ojEfficiencyLoader.fakeRateMu(eta, pt, drjet, nBJetsInEvent);
      }

      effReal[i] = real.e;
      uEffReal[i] = sqrt(real.eStat*real.eStat + real.eSystUncorr*real.eSystUncorr + real.eSystCorr*real.eSystCorr);

      effFake[i] = fake.e;
      uEffFakeStat[i] = fake.eStat;
      uEffFakeUncorr[i] = fake.eSystUncorr;
      uEffFakeCorr[i] = fake.eSystCorr;
      //uEffFakeBinIndex[i] = binFlavourOffset + binBJetOffset + rawBinScaleFactor * fake.bin_index;  // JO: changing bin indexing scheme
      uEffFakeBinIndex[i] = fake.bin_index;  // JO: changing bin indexing scheme

      if (isTight) configTLIn[i] = 1;
      else configTLIn[i] = 0;
	  
	  // if(debug_me)
	  // {
		  // std::cout<<"using real eff = "<<effReal[i]<<" and fake eff = "<< effFake[i]<<std::endl;
		  // std::cout<<"isTight = " << isTight << ", isElectron = " << isElectron << ", eta = " << eta << ", pt = " << pt << std::endl;
	  // }
	  
      //
      //// TODO refactor
      //bool requireFirstLeptonReal = false;
      //if (numLeptons > 2) requireFirstLeptonReal = true;
      //if (requireFirstLeptonReal && (i == 0)) {
      //  effReal[0] = 1.f;
      //  uEffReal[0] = 0.f;
      //  effFake[0] = 0.f;
      //  uEffFakeStat[0] = 0.f;
      //  uEffFakeUncorr[0] = 0.f;
      //  uEffFakeCorr[0] = 0.f;
      //  uEffFakeBinIndex[0] = 0;
      //  continue;
      //}
    }
  }

  void MatrixMethod::zeroConfiguration(char *config)
  {
    std::fill(config, config+numLeptons, 0);
  }

  MatrixMethod::Result MatrixMethod::resultForCurrentTLOutConfiguration()
  {
    float weight = 0.f;

    zeroConfiguration(configRF);
    zeroFloatArray(derivativesReal);
    zeroFloatArray(derivativesFake);
    while (!isLastConfiguration(configRF)) {
      if (!isWantedRFConfiguration()) {
        incrementConfiguration(configRF);
        continue;
      }

      float matrixElement = getMatrixElement();
      float inverseFactor = getInverseFactor();

      float thisWeight = matrixElement * inverseFactor;
      weight += thisWeight;
      incrementDerivatives(thisWeight);

      incrementConfiguration(configRF);
    }

    Result result(numLeptons, configTLOut, weight);

    addVarianceInformation(result);
    return result;
  }

  void MatrixMethod::zeroFloatArray(float *array)
  {
    std::fill(array, array+numLeptons, 0.f);
  }

  bool MatrixMethod::isLastConfiguration(const char *config) const
  {
    for (unsigned int i = 0; i < numLeptons; ++i ) {
      if (config[i] > 1) return true;
    }
    return false;
  }

  void MatrixMethod::incrementConfiguration(char *config)
  {
    for (unsigned int i = 0; i < numLeptons; ++i ) {
      if (config[i] == 0) {
        config[i] = 1;
        return;
      } else {
        config[i] = 0;
      }
    }
    
    // Flag as finished
    if (numLeptons > 0) config[0] = 2;
  }

  bool MatrixMethod::isWantedTLOutConfiguration() const
  {
    unsigned int nTight = numTightOrRealInConfig(configTLOut);
    return (nTight > 1);
  }

  bool MatrixMethod::isWantedRFConfiguration() const
  {
    //// TODO refactor
    //bool requireFirstLeptonReal = false;
    //if (numLeptons > 2) requireFirstLeptonReal = true;
    //if (requireFirstLeptonReal && !configRF[0]) return false;

    unsigned int nTightOutputLeptons = numTightOrRealInConfig(configTLOut);
    unsigned int nReal = numTightOrRealInConfig(configRF);
    return (nReal < nTightOutputLeptons);
    //return (!areConfigsEquivalent(configRF, configTLOut, numLeptons));
  }

  bool MatrixMethod::areConfigsEquivalent(const char *config1, const char *config2) const
  {
    for (unsigned int i = 0; i < numLeptons; ++i) {
      if (config1[i] != config2[i]) return false;
    }
    return true;
  }

  unsigned int MatrixMethod::numTightOrRealInConfig(const char *config) const
  {
    unsigned int nTight = 0;
    for (unsigned int i = 0; i < numLeptons; ++i ) {
      if (config[i] == 1) ++nTight;
    }
    return nTight;
  }

  float MatrixMethod::getMatrixElement() const
  {
    float element = 1.f;
    for (unsigned int i = 0; i < numLeptons; ++i ) {
      if (configTLOut[i] == 1 && configRF[i] == 1) {
        element *= effReal[i];
      }
      else if (configTLOut[i] == 1 && configRF[i] == 0) {
        element *= effFake[i];
      }
      else if (configTLOut[i] == 0 && configRF[i] == 1) {
        element *= 1.f - effReal[i];
      }
      else if (configTLOut[i] == 0 && configRF[i] == 0) {
        element *= 1.f - effFake[i];
      }
      else {
        std::cerr << "Error: unknown TLout, RF config" << std::endl;
        std::cerr << configToStr(configTLOut) << " " << configToStr(configRF) << std::endl;
        throw;
      }
    }
    return element;
  }

  float MatrixMethod::getInverseFactor() const
  {
    float element = 1.f;
    for (unsigned int i = 0; i < numLeptons; ++i ) {
      // Pick out cofactor
      if (configTLIn[i] == 1 && configRF[i] == 1) {
        element *= 1.f - effFake[i];
      }
      else if (configTLIn[i] == 0 && configRF[i] == 1) {
        element *= - effFake[i];
      }
      else if (configTLIn[i] == 1 && configRF[i] == 0) {
        element *= - (1.f - effReal[i]);
      }
      else if (configTLIn[i] == 0 && configRF[i] == 0) {
        element *= effReal[i];
      }
      else {
        std::cerr << "Error: unknown RF, TLin config" << std::endl;
        std::cerr << configToStr(configRF) << " " << configToStr(configTLIn) << std::endl;
        throw;
      }
      
      // Divide by determinant
      element /= effReal[i] - effFake[i];
    }
    return element;
  }

  void MatrixMethod::incrementDerivatives(float weight)
  {
    for (unsigned int i = 0; i < numLeptons; ++i) {
      float dTotaldReal = weight * (dLogPhidReal(i) + dLogInvPhidReal(i));
      float dTotaldFake = weight * (dLogPhidFake(i) + dLogInvPhidFake(i));
      
      derivativesReal[i] += dTotaldReal;
      derivativesFake[i] += dTotaldFake;
    }
  }

  float MatrixMethod::dLogPhidReal(unsigned int i) const
  {
    if (configTLOut[i] == 1 && configRF[i] == 1) {
      if (effReal[i] == 0.f) return 0.f;
      return (1.f / effReal[i]);
    }
    else if (configTLOut[i] == 1 && configRF[i] == 0) {
      return 0.f;
    }
    else if (configTLOut[i] == 0 && configRF[i] == 1) {
      if (effReal[i] == 1.f) return 0.f;
      return (1.f / (effReal[i] - 1.f));
    }
    else if (configTLOut[i] == 0 && configRF[i] == 0) {
      return 0.f;
    }
    else {
      std::cerr << "Error: unknown TLout, RF config" << std::endl;
      std::cerr << configToStr(configTLOut) << " " << configToStr(configRF) << std::endl;
      throw;
    }
  }

  float MatrixMethod::dLogPhidFake(unsigned int i) const
  {
    if (configTLOut[i] == 1 && configRF[i] == 1) {
      return 0.f;
    }
    else if (configTLOut[i] == 1 && configRF[i] == 0) {
      if (effFake[i] == 0.f) return 0.f;
      return (1.f / effFake[i]);
    }
    else if (configTLOut[i] == 0 && configRF[i] == 1) {
      return 0.f;
    }
    else if (configTLOut[i] == 0 && configRF[i] == 0) {
      if (effFake[i] == 1.f) return 0.f;
      return (1.f / (effFake[i] - 1.f));
    }
    else {
      std::cerr << "Error: unknown TLOut, RF config" << std::endl;
      std::cerr << configToStr(configTLOut) << " " << configToStr(configRF) << std::endl;
      throw;
    }
  }

  float MatrixMethod::dLogInvPhidReal(unsigned int i) const
  {
    float scale = 1.f / (effFake[i] - effReal[i]);
    if (configRF[i] == 1 && configTLIn[i] == 1) {
      return (scale * 1.f);
    }
    else if (configRF[i] == 1 && configTLIn[i] == 0) {
      return (scale * 1.f);
    }
    else if (configRF[i] == 0 && configTLIn[i] == 1) {
      if (effReal[i] == 1.f) return 0.f;
      return (scale * (effFake[i] - 1.f)/(effReal[i] - 1.f));
    }
    else if (configRF[i] == 0 && configTLIn[i] == 0) {
      if (effReal[i] == 0.f) return 0.f;
      return (scale * effFake[i] / effReal[i]);
    }
    else {
      std::cerr << "Error: unknown RF, TLIn config" << std::endl;
      std::cerr << configToStr(configRF) << " " << configToStr(configTLIn) << std::endl;
      throw;
    }
  }

  float MatrixMethod::dLogInvPhidFake(unsigned int i) const
  {
    float scale = 1.f / (effReal[i] - effFake[i]);
    if (configRF[i] == 1 && configTLIn[i] == 1) {
      if (effFake[i] == 1.f) return 0.f;
      return (scale * (effReal[i] - 1.f)/(effFake[i] - 1.f));
    }
    else if (configRF[i] == 1 && configTLIn[i] == 0) {
      if (effFake[i] == 0.f) return 0.f;
      return (scale * effReal[i] / effFake[i]);
    }
    else if (configRF[i] == 0 && configTLIn[i] == 1) {
      return (scale * 1.f);
    }
    else if (configRF[i] == 0 && configTLIn[i] == 0) {
      return (scale * 1.f);
    }
    else {
      std::cerr << "Error: unknown RF, TLin config" << std::endl;
      std::cerr << configToStr(configRF) << " " << configToStr(configTLIn) << std::endl;
      throw;
    }
  }

  void MatrixMethod::addVarianceInformation(MatrixMethod::Result &result) const
  {
    //unsigned int numBins = 2 * ojEfficiencyLoader.numFakeBinsEl() + ojEfficiencyLoader.numFakeBinsMu(); // JO: changing bin indexing scheme
    unsigned int numBins = ojEfficiencyLoader.numFakeBinsEl() + ojEfficiencyLoader.numFakeBinsMu(); // JO: changing bin indexing scheme
	// JO: changed names in the latter
    result.uncUncorr.clear(); 
    result.uncUncorr.resize(numBins, 0.f);
    result.uncElCorr = 0.f;
    result.uncMuCorr = 0.f;
    result.uncEff = 0.f;

    for (unsigned int bin = 0; bin < numBins; ++bin) {
      float sumOfDerivatives = 0.f;
      for (unsigned int lep = 0; lep < numLeptons; ++lep) {
        if (uEffFakeBinIndex[lep] != bin) continue;
        result.uncUncorr[bin] = sqrt(uEffFakeStat[lep]*uEffFakeStat[lep] + uEffFakeUncorr[lep]*uEffFakeUncorr[lep]); // JO: add sqrt()
        sumOfDerivatives += derivativesFake[lep];
      }
      result.uncUncorr[bin] *= sumOfDerivatives; // JO: variance -> signed uncertainty
    }


    float sumElFakeCorr = 0.f;
    float sumElReal = 0.f;
    float sumMuFakeCorr = 0.f;
    float sumMuReal = 0.f;
    for (unsigned int lep = 0; lep < numLeptons; ++lep) {
      if (isLepElectron[lep]) {
        sumElFakeCorr += derivativesFake[lep] * uEffFakeCorr[lep];
        sumElReal += fabs(derivativesReal[lep]) * uEffReal[lep];
      } else {
        sumMuFakeCorr += derivativesFake[lep] * uEffFakeCorr[lep];
        sumMuReal += fabs(derivativesReal[lep]) * uEffReal[lep];
      }
    }
    result.uncElCorr = sumElFakeCorr; //  JO: variance -> signed uncertainty, and remove real efficiency uncertainty
    result.uncMuCorr = sumMuFakeCorr; //  JO: variance -> signed uncertainty, and remove real efficiency uncertainty
	result.uncEff = sumElReal + sumMuReal; //  JO: treat real efficiency uncertainty separately
  }

  void MatrixMethod::clearMemory()
  {
    delete [] isLepElectron;
    delete [] effReal;
    delete [] effFake;
    delete [] uEffReal;
    delete [] uEffFakeStat;
    delete [] uEffFakeUncorr;
    delete [] uEffFakeCorr;
    delete [] uEffFakeBinIndex;
    delete [] configTLIn;
    delete [] configTLOut;
    delete [] configRF;
    delete [] derivativesReal;
    delete [] derivativesFake;
  }

  std::string MatrixMethod::configToStr(const char *config) const
  {
    std::string result;
    for (unsigned int i = 0; i < numLeptons; ++i) {
      char lep;
      if (config[i] == 0) {
        lep = '0';
      } else if (config[i] == 1) {
        lep = '1';
      } else {
        lep = config[i];
      }
      result += lep;
    }
    return result;
  }
}

