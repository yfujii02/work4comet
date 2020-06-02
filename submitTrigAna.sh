#!/bin/bash

NumRun=$1
for(( set=0; set<4; set++))
do
   while : ; do
     NJobs=$(ps aux |grep "python3 runTrig.py" |wc -l)
     if [ $NJobs -le 1 ];then
       break
     else
       sleep 300
     fi
   done
   for thr in 150 250 350 450 550 650 750 850
   do
      #python3 runTrig.py  0 84 $set $thr > out${set}_thr${thr}keV_deltaSegIsOne.txt &
      python3 runTrig.py  0 $NumRun $set $thr > out${set}_thr${thr}keV_deltaSegIsOne.txt &
   done
done
