/**
 * @file MatrixMethod.h
 * @author Thomas Gillam <thomas.gillam@cern.ch>
 * @date May, 2013
 * @brief Estimate fake components for events for arbitrary number of leptons
 *
 **/ 

#ifndef GENERALISEDMATRIXMETHOD_MATRIX_METHOD_H
#define GENERALISEDMATRIXMETHOD_MATRIX_METHOD_H

#include "OJEfficiencyLoader.h"
#include <vector>
#include <string>

namespace GeneralisedMatrixMethod {
  class MatrixMethod {
    public:
      // NB That pt's should be in GeV to be compatible!
      struct Lepton {
        Lepton(float pt, float eta, float drjet, bool isElectron, bool isTight)
          : pt(pt), eta(eta), drjet(drjet), isElectron(isElectron), isTight(isTight) { }

        float pt;
        float eta;
		float drjet;
        bool isElectron;
        bool isTight;
      };

      struct Result {
        Result(unsigned int nLeptons, char *configTLOut, float weight) : weight(weight)
        {
          for (unsigned int i = 0; i < nLeptons; ++i) {
            bool flag = (bool)configTLOut[i];
            leptonTightList.push_back(flag);
          }
        }

        std::vector<bool> leptonTightList;
        float weight;
        float uncElCorr;
        float uncMuCorr;
		float uncEff;
        std::vector<float> uncUncorr;
      };

      std::vector<Result> weightsForInput(const std::vector<Lepton> &leptons, int nBJetsInEvent);

    private:
      void setNumberOfLeptons();
      bool haveEnoughLeptons() const;
      void allocateMemory();
      void fillLeptonInformation(int nBJetsInEvent);
      void zeroConfiguration(char *config);
      Result resultForCurrentTLOutConfiguration();
      void zeroFloatArray(float *array);

      void clearMemory();

      bool isLastConfiguration(const char *config) const;
      void incrementConfiguration(char *config);
      bool isWantedTLOutConfiguration() const;
      bool isWantedRFConfiguration() const;
      bool areConfigsEquivalent(const char *config1, const char *config2) const;
      unsigned int numTightOrRealInConfig(const char *config) const;
      float getMatrixElement() const;
      float getInverseFactor() const;
      void incrementDerivatives(float weight);
      float dLogPhidReal(unsigned int i) const;
      float dLogPhidFake(unsigned int i) const;
      float dLogInvPhidReal(unsigned int i) const;
      float dLogInvPhidFake(unsigned int i) const;
      void addVarianceInformation(Result &result) const;

      std::string configToStr(const char *config) const;

    private:
      OJEfficiencyLoader ojEfficiencyLoader;

      const std::vector<Lepton> *leptons;
      unsigned int numLeptons;

      bool *isLepElectron;
      float *effReal;
      float *effFake;
      float *uEffReal;
      float *uEffFakeStat;
      float *uEffFakeUncorr;
      float *uEffFakeCorr;
      unsigned int *uEffFakeBinIndex;

      char *configTLIn;
      char *configTLOut;
      char *configRF;

      float *derivativesReal;
      float *derivativesFake;
  };
}

#endif // GENERALISEDMATRIXMETHOD_MATRIX_METHOD_H
