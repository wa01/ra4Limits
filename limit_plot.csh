#! /bin/tcsh
#
if ( $#argv != 1 ) then
  echo "Missing or additional arguments"
  exit 1
endif
if ( $1 =~ */ ) then
  echo "Remove trailing /"
  exit 1
endif
if ( !( -d $1) ) then
  echo "No such directory $1"
  exit 1
endif
if ( -e $1/$1.pkl ) then
  echo "Summary pickle file exists : $1/$1.pkl"
  exit 1
endif

python combineLimits.py $1 $1/$1.pkl
python limit_scan.py -b $1/$1.pkl $1/$1-scan
set limDir = "$PWD/$1"

pushd $CMSSW_BASE/src/CMS-SUS-XPAG/PlotsSMS
set limFile = "$limDir/$1-scan.root"
cat <<EOF >!  config/T5qqqqWW_SUS15006_tmp.cfg
# AVAILABLE COLORS:
# kMagenta
# kBlue
# kOrange
# kRed
#####################################################
#FORMAT: input root histo-name line-color area-color  
#####################################################
HISTOGRAM $limFile hXsec_exp_corr
EXPECTED $limFile graph_smoothed_Exp graph_smoothed_ExpP graph_smoothed_ExpM kRed kOrange
OBSERVED $limFile graph_smoothed_Obs graph_smoothed_ObsP graph_smoothed_ObsM kBlack kGray
# Preliminary Simulation or leave empty
PRELIMINARY Preliminary
# Lumi in fb 
LUMI 2.3
# Beam energy in TeV
ENERGY 13

EOF

python python/makeSMSplots.py  config/T5qqqqWW_SUS15006_tmp.cfg SUS15006-T5qqqqWW-
ls SUS15006-T5qqqqWW-*
mv SUS15006-T5qqqqWW-* $limDir/

popd

exit
