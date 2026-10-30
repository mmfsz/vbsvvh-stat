# Source, do not execute: sets up the one environment vbsvvh-stat runs in.
#   export VBSVVH_STAT_CMSSW=/path/to/CMSSW_16_0_0   # the area with Combine (README)
#   source setup.sh
if [ -z "${VBSVVH_STAT_CMSSW:-}" ]; then
    echo "setup.sh: set VBSVVH_STAT_CMSSW to your CMSSW area with Combine and CombineHarvester (see README)" >&2
    return 1
fi
source /cvmfs/cms.cern.ch/cmsset_default.sh
pushd "$VBSVVH_STAT_CMSSW/src" > /dev/null && eval "$(scramv1 runtime -sh)" && popd > /dev/null
VBSVVH_STAT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PATH="$VBSVVH_STAT_ROOT/bin:$PATH"
export PYTHONPATH="$VBSVVH_STAT_ROOT:${PYTHONPATH:-}"
