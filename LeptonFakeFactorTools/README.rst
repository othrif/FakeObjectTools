
=====================
LeptonFakeFactorTools
=====================

**Author:** Joseph Reichert

**Email:** Joseph.Reichert@cern.ch

**Description:** Package containing tools for SUSY fake lepton background estimation 
via the fake factor method.

-------------------------------

.. contents:: Table of contents


-------------------------------

-----
Setup
-----

To checkout the package:: 

    svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/SUSYPhys/FakeObjectTools/trunk/LeptonFakeFactorTools;

The package is structured in a similar manner to those in the AnalysisBase, so if 
you're using RootCore, you can rerun:: 

    rc find_packages; rc compile;

At this time, the requirements file for an Athena release does not exist, but if you 
are interested in using the package within Athena, let me know, because it should be 
simple enough to add.

-------------------------------

------------------------
Structure of the Package
------------------------

For a detailed discussion on the structure of the package, see 
`this talk <https://indico.cern.ch/event/569420/contributions/2305678/attachments/1337968/2013523/2016-09-16_Lepton_Fake_Factor_Tool.pdf>`_.
Note that some things in the package have been updated since the talk was presented,
but all of the big picture details remain the same (and I'll update the talk in the
near future).

In a nutshell, the following classes and macros exist:

EventDef.h
==========

Acts as the interface between your ntuple (or xAOD objects) and the rest of the 
fake factor classes. 

- To see how the inputs are defined according to the xAOD
  definitions, you can look at the fillEventInto(...), fillLepton(...), and
  fillMET(...) functions in Root/EventDef.cxx.

- The signal and antiID definitions are up to your analysis to define and calculate
  before filling the EventDef object. Note that an "other" ID also exists in the
  EventDef; this is in case your analysis were to require different two ID criteria 
  on the leptons (e.g. a TightLH leading electron and a MediumLH subleading electron).
  If you do not intend to use an "other" definition, you can set it equal to your
  nominal signal lepton ID.

BaseFakeFactorCalculator.h
==========================

The base class for the various fake factor selections to inherit from.

- Reads in conf files a la SUSYTools, which are located in the data/ directory.
  These conf files are where various cuts are defined (e.g. MT, MET, Mll), as well
  as whether or not you consider the given region to be your SR, the values of 
  any sample-dependent additional weights to apply, etc.

- Of course, if you want to add some variable that doesn't yet exist in the conf file, 
  you'll need to add it into the header (to declare the var) and the source code 
  (to properly read it from the conf file) of the base class, in addition 
  to making use of it inside of the derived fake factor calculator.

ApplyFakeFactor.h
=================

The tool for actually applying the fake factor.

- This class loads the previously derived fake factor histograms into memory, and has
  an apply() method to look up and return the fake factor for a given pt value.

- Note that currently, the pt must be specified in GeV!

- This tool has the technical ability to handle (pt, eta) binning, as well as handling
  tau leptons, but neither have been validated as of yet.

- This tool is implemented inside of the BaseFakeFactorCalculator, so one can use their
  derived fake factor calculator, and do everything necessary to obtain their background
  estimate in the SR. However, if it is simpler for you, the tool can also be called
  from within your analysis framework, and return the fake factor values there.

ZplusJetsFakeFactorCalculator.h
===============================

The class designed to perform a Z+jets event selection, and then extract the 
Z+jets fake factor.

- Derived from the base class, and takes the EventDef class as input. The outputs are
  numerator-level and denominator-level histograms, used for calculating the fake factor.
  Once a fake factor has already been obtained, this class can also be used to apply
  the fake factor to the appropriate CR, and extract the background estimate in the SR.

- Note that this class is used for the 3L Z+jets fake factor, as well as the 3L
  ttbar estimate (via a data-to-MC SF measured in DFOS events, which is then applied
  to the MC ttbar in the SR).

Other Code in the Package
=========================

- util/testZplusJetsWithPSLFakeNtuples.cxx is a script that serves as an example for
  how one can create an EventDef object and pass it into the fake factor calculator.

- Various macros exist in the macros/ directory for doing things such as dividing
  fake factor numerators by denominators to obtain the fake factor histogram root files,
  determining ttbar SFs in the 3L scenario, printing tables of fake background
  estimates, etc.

-------------------------------

----------------------------
Example: SM WZ 2015 Analysis
----------------------------

To go through the steps used in the SM WZ 2015 analysis for obtaining
background estimates, you can follow the steps described in test/README.
The ntuples used as input are located on AFS, and the testing scripts rely on
util/testZplusJetsWithPSLFakeNtuples.cxx

Note that the steps needed in the SM WZ 2015 analysis were actually somewhat 
more complicated than may be needed in your analysis! At the time, the low 
statistics meant the IDs used for measuring some fake factors had to be adjusted, 
not to mention that the W leptons and Z leptons used different ID requirements, 
which other analyses may not deem necessary. The antiID selection used here was
also looser than the baseline definition, which is why separate denominator and
numerator ntuples exist for the various steps! A simpler example will be added once 
it is available.


