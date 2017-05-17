import time, os, math, sys
import ROOT

ROOT.gROOT.SetBatch(True)

def main(options,args) :

    print 'Starting retrieve.py to retrieve file on Grid or Condor.'

    print 'Loading Root...'
    ROOT.gROOT.ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C")

    submitDir = options.dir

    if os.path.exists(submitDir) :
        ROOT.EL.Driver.wait(submitDir)

    return

if __name__ == "__main__" :
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('--dir', type='string', default='test', dest='dir', help='Output dir')

    (options,args) = p.parse_args()

    main(options,args)
    