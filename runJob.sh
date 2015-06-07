#!/bin/bash

#
# Run the job
#

Script=$0
Executable=$1
jobName=$2
nEventsPerJob=$3
job=$4
pileupLevel=$5

if [ $# -ne 5 ];then
    echo "Usage: ./runJob.sh <generateData> <jobName> <nEventsPerJob> <JobNumber> <pileupLevel>";
    exit 0;
fi

jobn=$jobName-`printf "%4.4d" $job`
ufsrootfile=`echo $jobn | awk '{printf("%s.root", $1)}'`
cardn=`echo $jobn | awk '{printf("%s.cards", $1)}'`

echo "################################################################################"
echo "EXECUTING THE FOLLOWING COMMAND...................."
echo "$Executable $cardn $nEventsPerJob $job $pileupLevel"
echo "################################################################################"

$Executable $cardn $nEventsPerJob $job $pileupLevel

#
# Check status
#

runStatus=$?
if [ "$runStatus" != "0" ]; then
    echo "$1 failed -- file not stored in HDFS"
    exit $runStatus
fi

#
# Transfer root file(s) to HDFS using lcg-cp
#

SRM_SERVER="srm://cmssrm1.hep.wisc.edu:8443/srm/v2/server?SFN="
FILE_DIR="/hdfs/store/user/$USER/UFS/$jobName"
echo $ufsrootfile*

for i in `\ls -1 $ufsrootfile*`
do

  lcg-cp -vv $i -bD srmv2 "$SRM_SERVER/$FILE_DIR/$i";
  exitstat=$?;
  if [ "$exitstat" -ne 0 ]; then
      echo "lcg-cp for $i to $FILE_DIR failed";
  else
      rm -rf "$i";
  fi

done
