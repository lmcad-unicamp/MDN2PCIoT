#!/bin/bash
#set -x

results=604
init_part=0
targetFile=2p-stm32f469xx
sourceGraph=LeNet-0604vertices_directed
lockedInput=0
klpMode=MDN2PCIoTcomm
cores=1
threads=1
same=50
samee=50
numberOfSubgraphs=1
experiment=$klpMode-$targetFile-free-$sourceGraph

mkdir $experiment
mkdir $experiment/results-$(hostname)/
mkdir $experiment/results-$(hostname)/results-$results

for i in $(seq 1);
        do { time ./MDN2PCIoT source-graphs/$sourceGraph.grf target-graphs/$targetFile.tgt $experiment/results-$(hostname)/results-$results/$i $init_part $klpMode verb $lockedInput $threads $same $samee $numberOfSubgraphs; } &> $experiment/results-$(hostname)/results-$results/result-$(hostname)-$i &
done; 
