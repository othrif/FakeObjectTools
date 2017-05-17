listOfIterations="
Nominal
"

listOfChannels="
ttt
"
confFile="$ROOTCOREBIN/../LeptonFakeFactorTools/data/measure_ttbar_SR_estimate/FakeFactorConfigs_Nominal.conf"
NumeratorDir="/afs/cern.ch/work/j/jreicher/public/FakeFactor_SMWZ_2015_Ntuples/Numerator/"

if [[ $HOSTNAME == *"upenn"* ]]; then
  NumeratorDir="/disk/userdata00/atlas_data2/kurb/2015/XRun2Susy_AB2345/NumeratorSelection_mc15b_LooseEle/fetch/data-SKIM/"
fi

echo $NumeratorDir

for arg in $listOfIterations;
do
  echo "starting $arg"
  mkdir "$arg";
  cd "$arg";

  for chan in $listOfChannels;
  do

    mkdir "$chan"_numerator;
    cd "$chan"_numerator;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/data.root $chan $confFile $arg; 
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/qqww.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/singletop.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/tother.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/ttbar.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/ttv.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/tw.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/tz.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/vvv.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/wz.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/zgam.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/zjet.root $chan $confFile $arg;
    testZplusJetsWithPSLFakeNtuples $NumeratorDir/zz.root $chan $confFile $arg;

    cd ..;
    hadd "$chan"_numerator/all.root "$chan"_numerator/*root;

  done

  cd ..; # get to top level dir again

done
