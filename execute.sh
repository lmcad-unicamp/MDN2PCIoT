#!/bin/bash
#set -x

experiment=comm-2p-locked-after-per-layers
results=604
init_part=initial-partitionings/lenet-604v-2p-49664-per-layers
targetFile=2p-stm32f469xx
sourceGraph=LeNet-0604vertices
lockedInput=256
klpMode=DN2PCIoTcomm
cores=1
threads=1

mkdir $experiment
mkdir $experiment/results-$(hostname)/
mkdir $experiment/results-$(hostname)/results-$results

for i in $(seq 1);
        do { time ./DN2PCIoT source-graphs/$sourceGraph.grf target-graphs/$targetFile.tgt $experiment/results-$(hostname)/results-$results/$i $init_part $klpMode verb $lockedInput $threads; } &> $experiment/results-$(hostname)/results-$results/result-$(hostname)-$i &
done; 
