import ROOT
import logging
import shutil
import os
import re

#Example of usage
# python Run.py -L -s -b -t --inputDS 387918 --driver grid

logging.basicConfig(level=logging.INFO)
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--submitDir", help="dir to store the output", default="submit_dir")
parser.add_option("-L", "--load", action='store_true', default=False, help="load previously saved Sample Handler")
parser.add_option("--inputDS", help="input DS from DQ2", default="none")
parser.add_option("--driver", help="select where to run", choices=("direct", "prooflite", "grid"), default="direct")
parser.add_option("--nevents", type=int, help="number of events to process for all the datasets")
parser.add_option("--skip-events", type=int, help="skip the first n events")
parser.add_option("-w", "--overwrite", action='store_false', default=True, help="overwrite previous submitDir")
parser.add_option("-d", "--debug", action='store_true', default=False, help="activate DEBUG mode")
parser.add_option("--processID", type=int, help="process ID flag (int)", default=-1)
parser.add_option("--doClass", action='store_true', default=False, help="fill probe lepton classification histo")



(options, args) = parser.parse_args()

import atexit
@atexit.register
def quite_exit():
	ROOT.gSystem.Exit(0)

logging.info("loading packages")
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

if options.overwrite:
	shutil.rmtree(options.submitDir, True)

#Set up the job for xAOD access:
ROOT.xAOD.Init().ignore();

# create a new sample handler to describe the data files we use
logging.info("creating new sample handler")
sh_all = ROOT.SH.SampleHandler()

if options.load:
	sh_all.load("SHandlers/");
else:
	if options.inputDS != "none":
		ROOT.SH.scanDQ2 (sh_all, options.inputDS);
	else:
		search_directories = []
		search_directories = ("/terabig/cmerlass/FakeRate_area/FakesAnalysis/test",)  
		# scan for datasets in the given directories
		for directory in search_directories:
			ROOT.SH.scanDir(sh_all, directory)

logging.info("%d different datasets found scanning all directories", len(sh_all))

#create the Sample Handler that will contain the actual datasets to be run over
sh_selected = ROOT.SH.SampleHandler()

for sample in sh_all:

	if (options.inputDS in sample.name()) or options.inputDS == "none":
		print "Found %s" % sample.name()		
		isatlfast = 0;
		isData = 0;

		if "data" in sample.name():
			isData = 1;
		
		prod = "mc15c"	
		print prod

		current=-1  
		first=-1
		last=-1
		for z in range(1,6):
			current = sample.name().find(".",current+1)
			if z==1:
				first = current+1
			if z==3:
				last = current
 	 
		short_name = sample.name()[first:last]

		kurrent = sample.name().find("_p2")
		etag = sample.name().find(".e")
		if sample.name()[kurrent-10]=="a":
			isatlfast = 1

                if isData == 0:
			if isatlfast == 1 :
				short_label=short_name+"."+sample.name()[etag+1:etag+7]+sample.name()[kurrent-10:kurrent-5]+sample.name()[kurrent+1:kurrent+6]
			else:
				short_label=short_name+"."+sample.name()[etag+1:etag+7]+sample.name()[kurrent-11:kurrent-5]+sample.name()[kurrent+1:kurrent+6]
		
                if isData == 1:
			short_label=short_name+"."+sample.name()[kurrent-11:kurrent-5]+sample.name()[kurrent+7:kurrent+12]

		sample.setMetaString("MC_production", prod)
		sample.setMetaString("short_name", short_name)
		sample.setMetaDouble("isatlfast", isatlfast)	
		sample.setMetaDouble("isdata", isData)	
		sh_selected.add(sample)


# print out the samples we found
logging.info("%d selected datasets", len(sh_selected))

# set the name of the tree in our files
sh_all.setMetaString("nc_tree", "CollectionTree")
  
# this is the basic description of our job
logging.info("creating new job")
job = ROOT.EL.Job()
job.sampleHandler(sh_selected)
job.options().setString(ROOT.EL.Job.optXaodAccessMode,ROOT.EL.Job.optXaodAccessMode_class);  


if options.nevents:
	logging.info("processing only %d events", options.nevents)
	job.options().setDouble(ROOT.EL.Job.optMaxEvents, options.nevents)

if options.skip_events:
	logging.info("skipping first %d events", options.skip_events)
	job.options().setDouble(ROOT.EL.Job.optSkipEvents, options.skip_events)

# add our algorithm to the job
logging.info("creating algorithms")
alg = ROOT.MyxAODAnalysis()
alg.m_debug = options.debug
alg.target_lumi = 36100
alg.doClassification = options.doClass

logging.info("adding algorithms")
job.algsAdd(alg)

# make the driver we want to use
# this one works by running the algorithm directly
logging.info("creating driver")
driver = None
if (options.driver == "direct"):
	logging.info("running on direct")
	driver = ROOT.EL.DirectDriver()
	logging.info("submit job")
	driver.submit(job, options.submitDir)
elif (options.driver == "prooflite"):
	logging.info("running on prooflite")
	driver = ROOT.EL.ProofDriver()
	logging.info("submit job")
	driver.numWorkers = 10;
	driver.submit(job, options.submitDir)
elif (options.driver == "grid"):
	logging.info("running on Grid") 
	driver = ROOT.EL.PrunDriver() 
	outname= "user."+os.environ["RUCIO_ACCOUNT"]+"."+short_label+".v06_gradientloose"
	driver.options().setString("nc_outputSampleName", outname)
	driver.options().setDouble("nc_disableAutoRetry", 0)
#	driver.options().setDouble("nc_maxNFilesPerJob", 1)
#	driver.options().setDouble("nc_nGBPerJob",6)
#	driver.options().setString("nc_optGridDestSE", "UNIBE-LHEP_LOCALGROUPDISK")
	
	logging.info("submit job")
	driver.submitOnly(job, options.submitDir)
	shutil.rmtree(options.submitDir, True)

logging.info("Done!")
