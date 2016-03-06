#!/bin/bash

cat $INPUT > input.txt
$HOME/CMSSW_3_11_2/src/UserCode/dasu/UltraFastSim/unl_zh input.txt $OUTPUT
