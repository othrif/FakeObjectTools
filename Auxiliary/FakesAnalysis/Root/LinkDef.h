#ifndef FAKESANALYSIS_LINKDEF_H
#define FAKESANALYSIS_LINKDEF_H

#include <vector>
#include <string>
#include <map>
#include <utility>

#include <FakesAnalysis/MyxAODAnalysis.h>
#include <FakesAnalysis/parametric_histos.h>

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#endif

#ifdef __CINT__
#pragma link C++ class MyxAODAnalysis+;
#pragma link C++ class pair<string,TTree >+;
#pragma link C++ class LoopSUSY_fakes+;
#endif

#endif
