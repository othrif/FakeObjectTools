listOfIterations="
WID_electron
ZID_electron
WID_muon
ZID_muon
"

confFile="$ROOTCOREBIN/../LeptonFakeFactorTools/data/Zjet_FF_calculation/FakeFactorConfigs_Zjet_FF_calculation.conf"
DenominatorDir="/afs/cern.ch/work/j/jreicher/public/FakeFactor_SMWZ_2015_Ntuples/Denominator/"
NumeratorDir="/afs/cern.ch/work/j/jreicher/public/FakeFactor_SMWZ_2015_Ntuples/Numerator/"

if [[ $HOSTNAME == *"upenn"* ]]; then
  DenominatorDir="/disk/userdata00/atlas_data2/kurb/2015/XRun2Susy_AB2345/FakesSelection_mc15b/fetch/data-SKIM/"
  NumeratorDir="/disk/userdata00/atlas_data2/kurb/2015/XRun2Susy_AB2345/NumeratorSelection_mc15b_LooseEle/fetch/data-SKIM/"
fi

echo $DenominatorDir
echo $NumeratorDir

for arg in $listOfIterations;
do
  echo "starting $arg"
  mkdir "$arg";
  cd "$arg";

  mkdir ttt_numerator;
  cd ttt_numerator;

  # Data
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/data.root ttt $confFile $arg; 

  # ttbar-like
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/qqww.root ttt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/ttbar.root ttt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/tw.root ttt $confFile $arg;

  # Zjet / Zgam
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/zgam.root ttt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/zjet.root ttt $confFile $arg;

  # Diboson
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/wz.root ttt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $NumeratorDir/zz.root ttt $confFile $arg;

  cd ..;
  mkdir ltt_denominator;
  cd ltt_denominator;

  # Data
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/data.root ltt $confFile $arg; 

  # ttbar-like
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/qqww.root ltt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/ttbar.root ltt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/tw.root ltt $confFile $arg;

  # Zjet / Zgam
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/zgam.root ltt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/zjet.root ltt $confFile $arg;

  # Diboson
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/wz.root ltt $confFile $arg;
  testZplusJetsWithPSLFakeNtuples $DenominatorDir/zz.root ltt $confFile $arg;

  cd ..;
  hadd ltt_denominator/all.root ltt_denominator/*root; hadd ttt_numerator/all.root ttt_numerator/*root;
  cd ..;
done
