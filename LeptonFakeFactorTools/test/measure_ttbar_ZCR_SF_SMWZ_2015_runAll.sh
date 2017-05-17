listOfIterations="
Nominal
"

listOfChannels="
ltt
tlt
ttl
"

confFile="$ROOTCOREBIN/../LeptonFakeFactorTools/data/measure_ttbar_ZCR_SF/FakeFactorConfigs_Nominal.conf"
DenominatorDir="/afs/cern.ch/work/j/jreicher/public/FakeFactor_SMWZ_2015_Ntuples/Denominator/"

if [[ $HOSTNAME == *"upenn"* ]]; then
  DenominatorDir="/disk/userdata00/atlas_data2/kurb/2015/XRun2Susy_AB2345/FakesSelection_mc15b/fetch/data-SKIM/"
fi

echo $DenominatorDir

for arg in $listOfIterations;
do
  echo "starting $arg"
  mkdir "$arg";
  cd "$arg";

  for chan in $listOfChannels;
  do

    mkdir "$chan"_denominator;
    cd "$chan"_denominator;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/data.root $chan $confFile $arg; 
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/qqww.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/singletop.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/tother.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/ttbar.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/ttv.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/tw.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/tz.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/vvv.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/wz.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/zgam.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/zjet.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $DenominatorDir/zz.root $chan $confFile $arg;

    cd ..;
    hadd "$chan"_denominator/all.root "$chan"_denominator/*root;

  done

  cd ..; # get to top level dir again

done
