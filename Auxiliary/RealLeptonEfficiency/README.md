This is the framework for real lepton efficiency (RLE) calculation

v04. New version of RLE package.
1. The new version of RLE packages become a plug-in module.
   Users can use their own framework and insert the skim module and 
   the RLE module in their framework.
2. The old version of RLE package is removed.

v03. Checkin CutflowsAndSkim and RLE packages 
Usage: Go into the directory and source the setup file.
Don't need to use checkout.sh

v02. Correct checkout.sh
Usage: ./checkout.sh
To checkout the necessary files.
You will get two directories CutflowsAndSkim and RLE.
Changing directory to one of them and source the setup file to setup
the rcSetp version.

v01. Initial commit
1. The current version of RLE framework is separated into two parts.
   The first part is cutflow & skim using TSelector to loop over events.
   And the second part is real lepton efficiency calculation which uses
   EventLoop (EL) algorithm to loop over events.
   The EL version of cutflow & skim is under developing. The benefit of
   using EL is that we can submit jobs to Condor or Grid to increase the
   processing performance.
2. Add checkout.sh
