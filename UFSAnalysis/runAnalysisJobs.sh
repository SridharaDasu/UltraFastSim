#!/bin/bash

cp /dev/null Config.txt
echo $OUTPUT >> Config.txt
#if [ `echo $OUTPUT | grep ZHmmbb` ]; then 
#    echo 17.44 >> Config.txt
#elif [ `echo $OUTPUT | grep ZZmmbb` ]; then 
#    echo 396.7 >> Config.txt 
#elif [ `echo $OUTPUT | grep Zmm` ]; then 
#    echo 1691000 >> Config.txt
#else
#    echo 1 >> Config.txt
#fi
#echo 100.0 >> Config.txt
cat $INPUT >> Config.txt
echo EOF >> Config.txt
#$HOME/UserCode/dasu/UltraFastSim/ZHOfflineAnalysis
$HOME/UserCode/dasu/UltraFastSim/ZHAnalysis
