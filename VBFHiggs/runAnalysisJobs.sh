#!/bin/bash
cp /dev/null Config.txt
echo $OUTPUT >> Config.txt
cat $INPUT >> Config.txt
$HOME/CMSSW_4_2_6/VBFHiggs/VBFHiggs Config.txt
