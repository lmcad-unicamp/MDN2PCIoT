/***************************************************************************
 *   Copyright (C) 2020 by Fabíola Martins Campos de Oliveira 		   	   *
 *   fabiola.bass@gmail.com			                           			   *
 *                       						   						   *
 *   This file is part of KLP, DN²PCIoT, and MDN²PCIoT.  				   *
 *                                      		   		   				   *
 *   KLP, DN²PCIoT, and MDN²PCIoT is free software: you can redistribute   *
 *   it and/or modify it under the terms of the GNU General Public License *
 *   as published by the Free Software Foundation, either version 3 of the * 
 *   License, or (at your option) any later version.				       *
 *									   									   *
 *   KLP, DN²PCIoT, and MDN²PCIoT is distributed in the hope that it will  *
 *   be useful,	but WITHOUT ANY WARRANTY; without even the implied 		   *	
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See * 
 *   the GNU General Public License for more details.			   		   *
 *									   									   *
 *   You should have received a copy of the GNU General Public License     *
 *   long with KLP, DN²PCIoT, and MDN²PCIoT.  If not, see 				   *
 *   <http://www.gnu.org/licenses/>.     								   *
 ***************************************************************************/

#include <cstddef>
#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include "ReadInputFiles.h"
#include "GenericPartitioningAlgorithm.h"
#include "Comm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

Comm::Comm(SourceGraph *srcG, int nPartitions, bool movingNodes, bool exchangingNodes, TargetGraph *tgtG, int nVertices, int forcedInput, bool initPart, char *verbose, char *initPartFile, int numberOfThreads) 
	: GenericPartitioningAlgorithm(nVertices, nPartitions, movingNodes, exchangingNodes, numberOfThreads),
	target(tgtG), sourceG(srcG) 
{
	validArray = (int *) calloc(nPartitions, sizeof(int));
	if(validArray == NULL) {
		cout << "\n validArray could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}

	validIndexes = (int *) calloc(nPartitions, sizeof(int));
	if(validIndexes == NULL) {
		cout << "\n validIndexes could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}

	validIndexesPerThread = (int **) malloc(getNumberOfThreads() * sizeof(int *));
	if (validIndexesPerThread == NULL) {
		cout << "\n validIndexesPerThread could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < getNumberOfThreads(); i++) {
		validIndexesPerThread[i] = (int *) malloc(nPartitions * sizeof(int));
		if (validIndexesPerThread[i] == NULL) {
			cout << "\n validIndexesPerThread[i] could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}		
	}

	targetIndexes = (int *) calloc(nPartitions, sizeof(int));
	if(targetIndexes == NULL) {
		cout << "\n targetIndexes could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < nPartitions; i++) {
		validIndexes[i] = i;
		targetIndexes[i] = i;
	}

	sharedMemoryPerPartition = (bool **) malloc(getNumberOfPartitions() * sizeof(bool *));
	if (sharedMemoryPerPartition == NULL) {
		cout << "\n sharedMemoryPerPartition could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < getNumberOfPartitions(); i++) {
		sharedMemoryPerPartition[i] = (bool *) malloc(sourceG->getNumberOfLayers() * sizeof(bool));
		if (sharedMemoryPerPartition[i] == NULL) {
			cout << "\n sharedMemoryPerPartition[i] could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}		
	}
	for (int i = 0; i < getNumberOfPartitions(); i++) {
		for (int j = 0; j < sourceG->getNumberOfLayers(); j++) {
			sharedMemoryPerPartition[i][j] = false;
		}
	}

	setInitialPartitionings(forcedInput, initPart, initPartFile);

	// to use in validPartitioning()
	sortedTargetMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
	if(sortedTargetMem == NULL) {
		cout << "\n sortedTargetMem could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < getNumberOfPartitions(); i++) {
		sortedTargetMem[i] = target->getMemory(i);
	}
	mergeSortWithIndexes(sortedTargetMem, 1, getNumberOfPartitions(), targetIndexes);
	// this makes the sorting not stable; to have a stable sort comment this block and use Source_graph->V/2 - 1 - a (or b) in FindMaxReductionInCost()
	for (int a = 0; a < getNumberOfPartitions()/2; a++) {
		int b = getNumberOfPartitions() - a - 1;
		int sortaux;

		sortaux = sortedTargetMem[a];
		sortedTargetMem[a] = sortedTargetMem[b];
		sortedTargetMem[b] = sortaux;

		sortaux = targetIndexes[a];
		targetIndexes[a] = targetIndexes[b];
		targetIndexes[b] = sortaux;
	}

	partitionGraphPerThread = (SourceGraph **) malloc(getNumberOfThreads() * sizeof(SourceGraph *));
	if (partitionGraphPerThread == NULL) {
		cout << "\n partitionGraphPerThread could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < getNumberOfThreads(); i++)
		partitionGraphPerThread[i] = new SourceGraph(getNumberOfPartitions(), getNumberOfPartitions());

	// build partitionGraph based on SourceG (input graph); eliminates redundant edges	
	partitionGraph = new SourceGraph(getNumberOfPartitions(), getNumberOfPartitions());

	double currentCost = 0;
	bool inserted;
	int cost, costRed, costMETIS = 0;
	for (int i = 0; i < getNumberOfVertices(); i++) {
		int iLayer = 0;
		for (int j = 0; j < sourceG->getNumberOfLayers(); j++) {
			if (iLayer == sourceG->getNumberOfLayers() - 1)
				break;
			if (i < sourceG->getLayerInitialPos(j + 1)) {
				break;
			} else {
				iLayer++;
			}
		}
		for (links nodeAdj = sourceG->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
			if (getCurrentPartitioning(i) != getCurrentPartitioning(nodeAdj->w) && nodeAdj->w > i) {
				inserted = partitionGraph->insertArcSrc(getCurrentPartitioning(i), getCurrentPartitioning(nodeAdj->w), nodeAdj->edgeWeight, nodeAdj->redundantNeurons, &cost, sourceG->getOriginalDepth(iLayer), iLayer, sourceG->getLayerInitialPos(iLayer), sourceG->getLayerWidth(iLayer), sourceG->getLayerHeight(iLayer), i);
				costMETIS += cost;
				if (inserted) {
					//currentCost += nodeAdj->edgeWeight;
					currentCost += cost;
					if ((strcmp(verbose, "verb")) == 0) {
						//cout << "\n i: " << i << " nodeAdj: " << nodeAdj->w << " eW: " << nodeAdj->edgeWeight << " redN: " << nodeAdj->redundantNeurons << " cost: " << currentCost;
					}
				}
			}
		}
	}
	if ((strcmp(verbose, "verb")) == 0) {
		//partitionGraph->printGraphSrc();	
	}

	cout << "\n costMETIS: " << costMETIS;
	currentCost = computeCost(getInitialPartitioning(), false);
	//cout << "\n cost: " << currentCost;
}

/*SourceGraph *Comm::getSourceGraph() {
	return sourceG;
}*/

void Comm::setSourceGraph(SourceGraph *srcG) {
	sourceG = srcG;
}

// Seed the random number generator.
void Comm::seed_rng(void)
{
    int fp = open("/dev/random", O_RDONLY);
    if (fp == -1) abort();
    unsigned seed;
    unsigned pos = 0;
    while (pos < sizeof(seed)) {
        int amt = read(fp, (char *) &seed + pos, sizeof(seed) - pos);
        if (amt <= 0) abort();
        pos += amt;
    }
    srand(seed);
    close(fp);
}

void Comm::setInitialPartitionings(int forcedInput, bool initPart, char *initPartFile) {
	int i, *j=NULL, k, r, l;

	setForcedInputSize(forcedInput);
	if (initPart == false) {
		// random, memory valid but different amount of vertices per partition
		// force input layer (image) to be in the same device/partition
		int p = 0;
		i = 0;
		while (i < getForcedInputSize()) {
			if (validArray[p] < target->getMemory(p)) {
				setVertexInitialPartitionings(i, p);
				validArray[p] += sourceG->getMemory(i);
			} else {
				p++;	
			}
			i++;
		}	
		seed_rng();
		//int shared;
		for (i = getForcedInputSize(); i < getNumberOfVertices(); i++) {
			r = rand();

			// finds the layer i belongs
			int layer = 0, shared = 0;
			for (int j = 0; j < sourceG->getNumberOfLayers(); j++) {
				if (layer == sourceG->getNumberOfLayers() - 1)
					break;
				if (i < sourceG->getLayerInitialPos(j + 1)) {
					break;
				} else {
					layer++;
				}
			}
			//cout << "\n layer: " << layer;
			if (sharedMemoryPerPartition[r % getNumberOfPartitions()][layer] == false) {
				//sharedMemoryPerPartition[r % getNumberOfPartitions()][layer] = true;
				shared = sourceG->getSharedParam(layer);
				//cout << " sMPP: false"; 
			}
			//else
				//cout << " sMPP: true";
			//cout << " shared: " << shared;
			int setp;
			if (validArray[r % getNumberOfPartitions()] + sourceG->getMemory(i) + shared < target->getMemory(r % getNumberOfPartitions())) {
				if (sharedMemoryPerPartition[r % getNumberOfPartitions()][layer] == false) {
					sharedMemoryPerPartition[r % getNumberOfPartitions()][layer] = true;
				}
				setVertexInitialPartitionings(i, r % getNumberOfPartitions());
				validArray[r % getNumberOfPartitions()] += sourceG->getMemory(i) + shared;
				//cout << "\n i: "<< i << " v[" << r % getNumberOfPartitions() << "]: " << validArray[r % getNumberOfPartitions()];
			} else {
				bool set = false;
				if (r % getNumberOfPartitions() == getNumberOfPartitions() - 1) {
					setp = - r % getNumberOfPartitions();
				} else {
					setp = 1;
				}
				int count = 0;
				while(!set) {
					if (r % getNumberOfPartitions() + setp == getNumberOfPartitions()) {
						setp = - r % getNumberOfPartitions();
						count++;
						if (count == getNumberOfPartitions()) {
							break;
						}
					}
					if (validArray[r % getNumberOfPartitions() + setp] + sourceG->getMemory(i) + shared < target->getMemory(r % getNumberOfPartitions() + setp)) {
						setVertexInitialPartitionings(i, r % getNumberOfPartitions() + setp);
						validArray[r % getNumberOfPartitions() + setp] += sourceG->getMemory(i) + shared;
						set = true;
					} else {
							setp++;
					}
				}
				if (count == getNumberOfPartitions()) {
					cout << "\n No valid initial partitioning was found! \n\n";
					exit(1);
				}
				//cout << "\n i: "<< i << " v[" << r % getNumberOfPartitions() + setp << "]: " << validArray[r % getNumberOfPartitions() + setp];
			}
		}
	}

	else if (initPart) {
		// input a specific initial partitioning
		FILE *fp;
		fp = fopen(initPartFile, "r");
		if(fp == NULL){
			cout << "initPart.txt file could not be open! \n";
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < getNumberOfVertices(); i++) {
			l = fscanf(fp, "%d ", &k);
			setVertexInitialPartitionings(i, k);
		}
		fclose(fp);
		fp=NULL;
	}
}

/* Returns the cost of partition the graph with partitioning p. */
/* For Comm, the cost is the amount of (external) communication 
	among partitions */
//double Comm::ComputeCost(Partitioning_t& p) {
double Comm::computeCost(const int *partitioning, bool openmp) {
	double currentCost = 0;
	bool inserted;

	if (openmp) {
		delete partitionGraphPerThread[omp_get_thread_num()];
		// build partitionGraph based on SourceG (input graph)
		partitionGraphPerThread[omp_get_thread_num()] = new SourceGraph(getNumberOfPartitions(), getNumberOfPartitions());

		//currentCost = 0;
		int cost = 0;
		for (int i = 0; i < getNumberOfVertices(); i++) {
			int iLayer = 0;
			for (int j = 0; j < sourceG->getNumberOfLayers(); j++) {
				if (iLayer == sourceG->getNumberOfLayers() - 1)
					break;
				if (i < sourceG->getLayerInitialPos(j + 1)) {
					break;
				} else {
					iLayer++;
				}
			}
			for (links nodeAdj = sourceG->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
				if (partitioning[i] != partitioning[nodeAdj->w] && nodeAdj->w > i) {
					inserted = partitionGraphPerThread[omp_get_thread_num()]->insertArcSrc(partitioning[i], partitioning[nodeAdj->w], nodeAdj->edgeWeight, nodeAdj->redundantNeurons, &cost, sourceG->getOriginalDepth(iLayer), iLayer, sourceG->getLayerInitialPos(iLayer), sourceG->getLayerWidth(iLayer), sourceG->getLayerHeight(iLayer), i);
					if (inserted) {
						//currentCost += nodeAdj->edgeWeight;
						currentCost += cost;
						//cout << "\n i: " << i << " nodeAdj: " << nodeAdj->w << " eW: " << nodeAdj->edgeWeight << " redN: " << nodeAdj->redundantNeurons << " cost: " << currentCost;
					}
				}
			}
		}

	} else {

		delete partitionGraph;
		// build partitionGraph based on SourceG (input graph)
		partitionGraph = new SourceGraph(getNumberOfPartitions(), getNumberOfPartitions());

		//currentCost = 0;
		int cost = 0;
		for (int i = 0; i < getNumberOfVertices(); i++) {
			int iLayer = 0;
			for (int j = 0; j < sourceG->getNumberOfLayers(); j++) {
				if (iLayer == sourceG->getNumberOfLayers() - 1)
					break;
				if (i < sourceG->getLayerInitialPos(j + 1)) {
					break;
				} else {
					iLayer++;
				}
			}
			for (links nodeAdj = sourceG->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
				if (partitioning[i] != partitioning[nodeAdj->w] && nodeAdj->w > i) {
					inserted = partitionGraph->insertArcSrc(partitioning[i], partitioning[nodeAdj->w], nodeAdj->edgeWeight, nodeAdj->redundantNeurons, &cost, sourceG->getOriginalDepth(iLayer), iLayer, sourceG->getLayerInitialPos(iLayer), sourceG->getLayerWidth(iLayer), sourceG->getLayerHeight(iLayer), i);
					if (inserted) {
						//currentCost += nodeAdj->edgeWeight;
						currentCost += cost;
						//cout << "\n i: " << i << " nodeAdj: " << nodeAdj->w << " eW: " << nodeAdj->edgeWeight << " redN: " << nodeAdj->redundantNeurons << " cost: " << currentCost;
					}
				}
			}
		}
	}

	return currentCost;

	// considers all edges to compute cost, even if it is a redundant edge, i.e., an edge that comes from the same vertex and go to the same partition
	/*unsigned a, i;
	//double cost;
	links nodeAdj;
	//SourceGraph *srcG = getSourceGraph();

	cost = 0;
	for (a = 0; a < sourceG->getNumberOfVertices(); a++) {
		//if (partitioning[a] == 0) {
			for (nodeAdj = sourceG->getAdjOfVertex(a); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
				if (partitioning[nodeAdj->w] != partitioning[a]) {
					cout << "\n a: " << a << " edge: " << nodeAdj->w << " cost: " << cost << " eW: " << nodeAdj->edgeWeight;
					cost += nodeAdj->edgeWeight; // * buildNetworkTopology(16, 4, partitioning[a], partitioning[nodeAdj->w]); João
				}
			}
		//}
	}
	//cout << "\ncost: " << cost;
	return cost/2;*/
}

void Comm::computeRedundantMemoryPerPartition(const int *p) {
	// Calculate amount of memory needed for each partition	
	for (int i = 0; i < getNumberOfPartitions(); i++) {
		validArray[i] = 0;
	}
	for (int i = 0; i < sourceG->getNumberOfVertices(); i++) {
		validArray[p[i]] += sourceG->getMemory(i);
	}
}

void Comm::computeNonredundantMemoryPerPartition(const int *p) {
	int i, j;

#	pragma omp parallel num_threads(getNumberOfThreads())
	{
#	pragma omp for
	for (i = 0; i < getNumberOfPartitions(); i++) {
		validArray[i] = 0;
		for (j = 0; j < sourceG->getNumberOfLayers(); j++) {
			sharedMemoryPerPartition[i][j] = false;
		}
	}

	// Calculate amount of memory needed for each partition
	j = 0; // layer
#	pragma omp for private(j)
	for (i = 0; i < sourceG->getNumberOfVertices(); i++) {
#		pragma omp critical
		validArray[p[i]] += sourceG->getMemory(i);

		j = 0;
		while (j < sourceG->getNumberOfLayers() - 1) {
			if (i < sourceG->getLayerInitialPos(j + 1)) {
				break;
			} else {
				j++;
			}
		}

#		pragma omp critical
		sharedMemoryPerPartition[p[i]][j] = true;
	}
#	pragma omp for private(j)
	for (i = 0; i < getNumberOfPartitions(); i++) {
		//cout << "\nThread[" << omp_get_thread_num() << "]";
		for (j = 0; j < sourceG->getNumberOfLayers(); j++) {
			if (sharedMemoryPerPartition[i][j] == true) {
				validArray[i] += sourceG->getSharedParam(j);
			}
		}
	}
	} // pragma
}

bool Comm::partitionsFitDevices(const int *p, int *sourceMem, bool openmp) {
	int a, b;

	mergeSortWithIndexes(sourceMem, 1, getNumberOfPartitions(), validIndexes);

	if (openmp) {
		// nonnaive implementation: sort validArray and targetMem in descending order
		// this makes the sorting not stable; to have a stable sort comment this block and use Source_graph->V/2 - 1 - a (or b) in FindMaxReductionInCost()
	//#	pragma omp parallel for num_threads(getNumberOfThreads()) private(b)
		for (a = 0; a < getNumberOfPartitions()/2; a++) {
			b = getNumberOfPartitions() - a - 1;

			int sortaux = sourceMem[a];
			sourceMem[a] = sourceMem[b];
			sourceMem[b] = sortaux;

			sortaux = validIndexesPerThread[omp_get_thread_num()][a];
			validIndexesPerThread[omp_get_thread_num()][a] = validIndexesPerThread[omp_get_thread_num()][b];
			validIndexesPerThread[omp_get_thread_num()][b] = sortaux;
		}
	} else {
		// nonnaive implementation: sort validArray and targetMem in descending order
		// this makes the sorting not stable; to have a stable sort comment this block and use Source_graph->V/2 - 1 - a (or b) in FindMaxReductionInCost()
	//#	pragma omp parallel for num_threads(getNumberOfThreads()) private(b)
		for (a = 0; a < getNumberOfPartitions()/2; a++) {
			b = getNumberOfPartitions() - a - 1;

			int sortaux = sourceMem[a];
			sourceMem[a] = sourceMem[b];
			sourceMem[b] = sortaux;

			sortaux = validIndexes[a];
			validIndexes[a] = validIndexes[b];
			validIndexes[b] = sortaux;
		}
	}

	for (a = 0; a < getNumberOfPartitions(); a++) {
		if (sourceMem[a] > sortedTargetMem[a]) {
			free(sourceMem);
			return false;
		}
	}

	free(sourceMem);
	return true;
}

/* Returns true if the partitioning is valid, false otherwise. */
bool Comm::validPartitioning(const int *p) {

	if (sourceG->getEnableMemory() == 0) {
		return true;
	} else if (sourceG->getEnableMemory() == 1) {
		// Calculate amount of memory needed for each partition
		computeNonredundantMemoryPerPartition(p);

		//static int *sourceMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
		int *sourceMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
		if(sourceMem == NULL) {
			cout << "\n sourceMem could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}
#		pragma omp parallel for num_threads(getNumberOfThreads())
		for (int a = 0; a < getNumberOfPartitions(); a++) {
			sourceMem[a] = validArray[a];
		}

		// Check if every partition fits to devices
		return partitionsFitDevices(p, sourceMem, false);
	}
}

void Comm::computeDiffNodeMemoryPerPartition(int *partitioning, 
								   				   unsigned nodeOrK,
												   unsigned originalNodeOrKPartition,
												   int *sourceMem) {
	int i, layer = 0;
	bool sumSharedMem = true, subtractSharedMem = true;

	// vertice (node or k) memory
	sourceMem[partitioning[nodeOrK]] += sourceG->getMemory(nodeOrK);
	sourceMem[originalNodeOrKPartition] -= sourceG->getMemory(nodeOrK);

	// finds the layer node or k belongs
	for (i = 0; i < sourceG->getNumberOfLayers(); i++) {
		if (layer == sourceG->getNumberOfLayers() - 1)
			break;
		if (nodeOrK < sourceG->getLayerInitialPos(i + 1)) {
			break;
		} else {
			layer++;
		}
	}

	int last;
	if (layer >= sourceG->getNumberOfLayers() - 1) {
		last = sourceG->getNumberOfVertices();
	} else {
		last = sourceG->getLayerInitialPos(layer + 1);
	}
	// searches if new node or k partition must sum shared memory
	for (i = sourceG->getLayerInitialPos(layer); i < last; i++) {
		if (i == nodeOrK) continue;
		if (partitioning[i] == partitioning[nodeOrK]) {
			sumSharedMem = false;
			break;
		}
	}
	if (sumSharedMem) {
		sourceMem[partitioning[nodeOrK]] += sourceG->getSharedParam(layer);
	}

	// searches if old node or k partition must subtract shared memory
	for (i = sourceG->getLayerInitialPos(layer); i < last; i++) {
		if (i == nodeOrK) continue;
		if (partitioning[i] == originalNodeOrKPartition) {
			subtractSharedMem = false;
			break;
		}
	}
	if (subtractSharedMem) {
		sourceMem[originalNodeOrKPartition] -= sourceG->getSharedParam(layer);
	}
}

void Comm::computeDiffMemoryPerPartition(int *partitioning, 
								   				 unsigned node, 
								   				 unsigned k, 
								   				 unsigned originalNodePartition, 
								   				 unsigned original_k_Partition,
												 int *sourceMem,
												 bool singleMoveOrSwap) {

	// compute for node
	computeDiffNodeMemoryPerPartition(partitioning, node, originalNodePartition, sourceMem);

	// compute for k in case of swaps	
	if (singleMoveOrSwap) {
		computeDiffNodeMemoryPerPartition(partitioning, k, original_k_Partition, sourceMem);
	}
}

/* Returns true if the partitioning is valid, false otherwise. */
bool Comm::diffValidPartitioning(int *partitioning,
								   	   unsigned node, 
								   	   unsigned k, 
								   	   unsigned originalNodePartition, 
								   	   unsigned original_k_Partition,
					bool singleOrSwap) {
	//SourceGraph *srcG = getSourceGraph();

	if (sourceG->getEnableMemory() == 0) {
		return true;
	} else if (sourceG->getEnableMemory() == 1) {

		//static int *sourceMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
		int *sourceMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
		if(sourceMem == NULL) {
			cout << "\n sourceMem could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}
		for (int a = 0; a < getNumberOfPartitions(); a++) {
			sourceMem[a] = validArray[a];
		}

		// Calculate amount of memory needed for each partition
		computeDiffMemoryPerPartition(partitioning, node, k, originalNodePartition, original_k_Partition, sourceMem, singleOrSwap);

		// Check if every partition fits to devices
		return partitionsFitDevices(partitioning, sourceMem, true);
		//return false;
	}
}

int Comm::getValidArray(int i) const {
	return validArray[i];
}

int Comm::getValidIndex(int i) const {
	return validIndexes[i];
}

void Comm::mergeSortWithIndexes(int *A, int p, int r, int *indexes){
	int j;

	//static int **Ai = (int **) calloc(r, sizeof(int *));
	int **Ai = (int **) calloc(r, sizeof(int *));
	if (Ai == NULL) {
		printf("Ai could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
/**#	pragma omp parallel num_threads(getNumberOfThreads())
	{
#	pragma omp for*/
	for (j = 0; j < r; j++) {
		Ai[j] = (int *) calloc(2, sizeof(int));
		if (Ai[j] == NULL) {
			printf("Ai[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}

//#	pragma omp for
	for (j = 0; j < r; j++) {
		Ai[j][0] = j;
		Ai[j][1] = A[j];
	}
	//} // pragma

	/*for (j = p; j <= r; j++) {
		printf("\n%d %d", Ai[j - 1][0], Ai[j - 1][1]);
	}
	printf("\n");*/

	mergeSort(Ai, p, r);

/*#	pragma omp parallel num_threads(getNumberOfThreads())
	{
#	pragma omp for*/
	for (j = p; j <= r; j++) {
		A[j - 1] = Ai[j - 1][1];
		indexes[j - 1] = Ai[j - 1][0];
	}

	/*for (j = p; j <= r; j++) {
		printf("\n%d %d", Ai[j - 1][0], Ai[j - 1][1]);
	}
	printf("\n");*/

//#	pragma omp for
	for (j = 0; j < r; j++) {
		free(Ai[j]);
		Ai[j] = NULL;
	}	
	//} // pragma

	free(Ai);
	Ai=NULL;
}

void Comm::mergeSort(int **A, int p, int r) {
	int j, q;

	if (p < r) {
		//q = floor((p+r)/2);
		q = (p + r)/2;
		mergeSort(A, p, q);
		mergeSort(A, q + 1, r);
		merge(A, p, q, r);
	}
}

void Comm::merge(int **A, int p, int q, int r) {
	int n1, n2, **L=NULL, **R=NULL, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	L = (int **) calloc(n1 + 1, sizeof(int *));
	if (L == NULL) {
		printf("L could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < n1 + 1; j++) {
		L[j] = (int *) calloc(2, sizeof(int));
		if (L[j] == NULL) {
			printf("L[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}
	R = (int **) calloc(n2 + 1, sizeof(int *));
	if (R == NULL) {
		printf("R could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < n2 + 1; j++) {
		R[j] = (int *) calloc(2, sizeof(int));
		if (R[j] == NULL) {
			printf("R[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}

	for (i = 1; i < n1 + 1; i++) {
		L[i - 1][1] = A[p + i - 2][1];
		L[i - 1][0] = A[p + i - 2][0];
	}
	for (j = 1; j < n2 + 1; j++) {
		R[j - 1][1] = A[q + j - 1][1];
		R[j - 1][0] = A[q + j - 1][0];
	}

	L[n1][1] = 2147483647;
	L[n1][0] = -1;
	R[n2][1] = 2147483647;
	R[n2][0] = -1;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++) {
		if (L[i - 1][1] <= R[j - 1][1]) {
			A[k - 1][1] = L[i - 1][1];
			A[k - 1][0] = L[i - 1][0];
			i++;
		} else {
			A[k - 1][1] = R[j - 1][1];
			A[k - 1][0] = R[j - 1][0];
			j++;
		}
	}

	for (j = 0; j < n1 + 1; j++) {
		free(L[j]);
		L[j] = NULL;
	}	
	free(L);
	L=NULL;

	for (j = 0; j < n2 + 1; j++) {
		free(R[j]);
		R[j] = NULL;
	}	
	free(R);
	R=NULL;
}

Comm::~Comm() {
	if (validArray != NULL) {
		free(validArray);
		validArray = NULL;
	}
	if (validIndexes != NULL) {
		free(validIndexes);
		validIndexes = NULL;
	}
	if (targetIndexes != NULL) {
		free(targetIndexes);
		targetIndexes = NULL;
	}
	if (sortedTargetMem != NULL) {
		free(sortedTargetMem);
		sortedTargetMem = NULL;
	}
	if (sharedMemoryPerPartition != NULL) {
		for (int i = 0; i < getNumberOfPartitions(); i++) {
			if (sharedMemoryPerPartition[i] != NULL) {
				free(sharedMemoryPerPartition[i]);
			}
		}
		free(sharedMemoryPerPartition);
		sharedMemoryPerPartition = NULL;
	}
	delete partitionGraph;
	if (partitionGraphPerThread != NULL) {
		for (int i = 0; i < getNumberOfThreads(); i++) {
			if (partitionGraphPerThread[i] != NULL) {
				free(partitionGraphPerThread[i]);
			}
		}
		free(partitionGraphPerThread);
		partitionGraphPerThread = NULL;
	}
	if (validIndexesPerThread != NULL) {
		for (int i = 0; i < getNumberOfThreads(); i++) {
			if (validIndexesPerThread[i] != NULL) {
				free(validIndexesPerThread[i]);
			}
		}
		free(validIndexesPerThread);
		validIndexesPerThread = NULL;
	}
}

