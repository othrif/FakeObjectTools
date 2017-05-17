/**
 * @file OJEfficiencyLoader.h
 * @author Thomas Gillam <thomas.gillam@cern.ch>
 * @date May, 2013
 * @brief Load fake rates computed by Julien & Otilia
 *
 **/ 

#ifndef GENERALISEDMATRIXMETHOD_OJ_EFFICENCY_LOADER_H
#define GENERALISEDMATRIXMETHOD_OJ_EFFICENCY_LOADER_H

namespace GeneralisedMatrixMethod {
  class OJEfficiencyLoader {
    public:
      struct Efficiency {
        Efficiency() : e(0.), eStat(0.), eSystCorr(0.), eSystUncorr(0.), eSystExtrap(0.), bin_index(0) { }
        double e, eStat, eSystCorr, eSystUncorr, eSystExtrap;
        int bin_index;
      };
      
      Efficiency realEfficiencyEl(float eta, float ptInGeV, float drjet);
      Efficiency realEfficiencyMu(float eta, float ptInGeV, float drjet);
      Efficiency fakeRateEl(float eta, float ptInGeV, float drjet, int nBJetsInEvent);
      Efficiency fakeRateMu(float eta, float ptInGeV, float drjet, int nBJetsInEvent);

      unsigned int numFakeBinsEl() const;
      unsigned int numFakeBinsMu() const;
  };
}

#endif // GENERALISEDMATRIXMETHOD_OJ_EFFICENCY_LOADER_H
