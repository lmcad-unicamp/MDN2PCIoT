#!/bin/bash
#set -x

# results: folder name
results=604

# init_part: 0 for random initial partitioning or the file path for a defined partitioning
init_part=0

# targetFile: file for the target graph
targetFile=2p-stm32f469xx

# sourceGraph: file for the source graph
sourceGraph=LeNet-0604vertices_directed

# lockedInput: 0 for free-input-layer experiments or the input layer number of vertices for locked-input-layer experiments (256 for LeNet 2:1 and 1024 for LeNet 1:1)
lockedInput=0

# klpMode: KLPcomm, DN2PCIoTcomm or MDN2PCIoTcomm for communication reduction or KLPiR, DN2PCIoTiR or MDN2PCIoTiR for inference rate maximization, according to the chosen algorithm
klpMode=MDN2PCIoTiR

# cores: number of executions
cores=1

# threads: number of threads for OpenMP parallelization
threads=1

# same: swap stabilization (after same times that cost has the same result, the findBestNodeExchange() function finishes and current best vertex is chosen to swap; this number is fixed for smaller graphs)
same=50

# samee: step stabilization (after samee times that cost has the same result, the step finishes and the best operation found so far is chosen)
samee=50

# numberOfSubgraphs: number of subgraphs produced in the coarsening phase minus one (numberOfSubgraphs=1 equals to the original source graph; no subgraphs are produced). This number is fixed for smaller graphs.
numberOfSubgraphs=1

# algorithm: respective algorithm according to the commit (KLP, DN2PCIoT, or MDN2PCIoT)
algorithm=MDN2PCIoT

# experiment: folder name
experiment=$klpMode-$targetFile-free-$sourceGraph

mkdir $experiment
mkdir $experiment/results-$(hostname)/
mkdir $experiment/results-$(hostname)/results-$results

for i in $(seq 1);
        do { time ./$algorithm source-graphs/$sourceGraph.grf target-graphs/$targetFile.tgt $experiment/results-$(hostname)/results-$results/$i $init_part $klpMode verb $lockedInput $threads $same $samee $numberOfSubgraphs; } &> $experiment/results-$(hostname)/results-$results/result-$(hostname)-$i &
done; 

