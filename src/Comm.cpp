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
#include <sys/time.h>
#include "ReadInputFiles.h"
#include "GenericPartitioningAlgorithm.h"
#include "Comm.h"

using namespace std;

Comm::Comm(SourceGraph *srcG, int nPartitions, bool movingNodes, bool exchangingNodes, TargetGraph *tgtG, int nVertices, int forcedInput, bool initPart) : GenericPartitioningAlgorithm(nVertices, nPartitions, movingNodes, exchangingNodes) {
	target = tgtG;
	setSourceGraph(srcG);

	validArray = (int *) calloc(nPartitions, sizeof(int));
	if(validArray == NULL) {
		cout << "\n validArray could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}

	setInitialPartitionings(forcedInput, initPart);

	validIndexes = (int *) calloc(nPartitions, sizeof(int));
	if(validIndexes == NULL) {
		cout << "\n validIndexes could not be allocated! \n";
		exit(EXIT_FAILURE);	
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

	// build partitionGraph based on SourceG (input graph)	
	partitionGraph = new SourceGraph(getNumberOfPartitions(), getNumberOfPartitions());

	currentCost = 0;
	bool inserted;
	for (int i = 0; i < getNumberOfVertices(); i++) {
		for (link nodeAdj = getSourceGraph()->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
			if (getCurrentPartitioning(i) != getCurrentPartitioning(nodeAdj->w) && nodeAdj->w > i) {
				//cout << "\n i: " << i << " nodeAdj: " << nodeAdj->w;
				inserted = partitionGraph->insertArcSrc(getCurrentPartitioning(i), getCurrentPartitioning(nodeAdj->w), nodeAdj->edgeWeight, i);
				if (inserted) {
					currentCost += nodeAdj->edgeWeight;
				}
			}
		}
	}
	partitionGraph->printGraphSrc();

	//cout << "\n cost: " << currentCost;
	currentCost = computeCost(getInitialPartitioning());
	//cout << "\n cost: " << currentCost;
}

SourceGraph *Comm::getSourceGraph() {
	return sourceG;
}

void Comm::setSourceGraph(SourceGraph *srcG) {
	sourceG = srcG;
}

void Comm::setInitialPartitionings(int forcedInput, bool initPart) {
	int i, *j=NULL, k, r, l;

	// random equally sized
	/*srand(time(NULL));
	j = (int *) calloc(getNumberOfPartitions(), sizeof(int));
	for (i = 0; i < getNumberOfVertices(); i++) {
		r = rand();
		//cout << "\n i: " << i << " r: " << r;
		if (j[r % getNumberOfPartitions()] < getNumberOfVertices()/getNumberOfPartitions()) {
			setVertexInitialPartitionings(i, r % getNumberOfPartitions());
			//initialPartitioning[i] = r % getNumberOfPartitions();
			j[r % getNumberOfPartitions()]++;
			//cout << "\n i: "<< i << " j[" << r % getNumberOfPartitions() << "]: " << j[r % getNumberOfPartitions()];
		} else {
			bool set = false;
			int setp;
			if (r % getNumberOfPartitions() == getNumberOfPartitions() - 1) {
				setp = - r % getNumberOfPartitions();
			} else {
				setp = 1;
			}
			while(!set) {
				if (j[r % getNumberOfPartitions() + setp] < getNumberOfVertices()/getNumberOfPartitions() && (r % getNumberOfPartitions() + setp) < getNumberOfPartitions()) {
					setVertexInitialPartitionings(i, r % getNumberOfPartitions() + setp);
					//initialPartitioning[i] = r % getNumberOfPartitions() + setp;
					j[r % getNumberOfPartitions() + setp]++;
					set = true;
				} else {
					if (r % getNumberOfPartitions() + setp == getNumberOfPartitions()) {
						setp = - r % getNumberOfPartitions();
					} else {
						setp++;
					}
				}
			}
		}
	}*/

	if (initPart == false) {
		// random, memory valid but different amount of vertices per partition
		// force input layer (image) to be in the same device/partition
		int p = 0;
		i = 0;
		setForcedInputSize(forcedInput);
		while (i < getForcedInputSize()) {
			if (validArray[p] < target->getMemory(p)) {
				setVertexInitialPartitionings(i, p);
				validArray[p] += sourceG->getMemory(i);
			} else {
				p++;	
			}
			i++;
		}	
		/*struct timeval tim;
		gettimeofday(&tim, NULL);
		double initial = tim.tv_sec+(tim.tv_usec/1000000.0);*/
		srand(time(NULL));
		// if input will not be forced, use:
		//for (i = 0; i < getNumberOfVertices(); i++) {
		for (i = getForcedInputSize(); i < getNumberOfVertices(); i++) {
			r = rand();
			if (validArray[r % getNumberOfPartitions()] + sourceG->getMemory(i) < target->getMemory(r % getNumberOfPartitions())) {
				setVertexInitialPartitionings(i, r % getNumberOfPartitions());
				validArray[r % getNumberOfPartitions()] += sourceG->getMemory(i);
			} else {
				bool set = false;
				int setp;
				if (r % getNumberOfPartitions() == getNumberOfPartitions() - 1) {
					setp = - r % getNumberOfPartitions();
				} else {
					setp = 1;
				}
				while(!set) {
					if (validArray[r % getNumberOfPartitions() + setp] + sourceG->getMemory(i) < target->getMemory(r % getNumberOfPartitions() + setp)) {
						setVertexInitialPartitionings(i, r % getNumberOfPartitions() + setp);
						validArray[r % getNumberOfPartitions() + setp] += sourceG->getMemory(i);
						set = true;
					} else {
						if (r % getNumberOfPartitions() + setp == getNumberOfPartitions()) {
							setp = - r % getNumberOfPartitions();
						} else {
							setp++;
						}
					}
				}
			}
		}
	}

	// round-Robin for 2 partitions
    /*for (i = 0; i < sourceG->getNumberOfVertices(); i++) {
        if (i % 2 == 0) {
            initialPartitioning[i] = 0;
        } else {
            initialPartitioning[i] = 1;
        }
    }*/

	// input a specific initial partitioning (2 partitionings, in order to use same initial partitioning as KL.c
	/*FILE *fp;
	fp = fopen("initPart.txt", "r");
	if(fp == NULL){
		cout << "initPart.txt file could not be open! \n";
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < getNumberOfVertices(); i++) {
		initialPartitioning[i] = 1;
	}
	for (i = 0; i < getNumberOfVertices()/2; i++) {
		j = fscanf(fp, "%d ", &k);
		initialPartitioning[k] = 0;
	}
	fclose(fp);
	fp=NULL;*/

	else if (initPart) {
		// input a specific initial partitioning
		FILE *fp;
		fp = fopen("../initPart.txt", "r");
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
double Comm::computeCost(int *partitioning) {
	double cost = 0;
	bool inserted;

	delete partitionGraph;
	// build partitionGraph based on SourceG (input graph)
	partitionGraph = new SourceGraph(getNumberOfPartitions(), getNumberOfPartitions());

	for (int i = 0; i < getNumberOfVertices(); i++) {
		for (link nodeAdj = getSourceGraph()->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
			if (partitioning[i] != partitioning[nodeAdj->w] && nodeAdj->w > i) {
				//cout << "\n i: " << i << " nodeAdj: " << nodeAdj->w;
				inserted = partitionGraph->insertArcSrc(partitioning[i], partitioning[nodeAdj->w], nodeAdj->edgeWeight, i);
				if (inserted) {
					cost += nodeAdj->edgeWeight;
				}
			}
		}
	}
	return cost;

	// considers all edges to compute cost, even if it is a redundant edge, i.e., an edge that comes from the same vertex and go to the same partition
	/*unsigned a, i;
	//double cost;
	link nodeAdj;
	SourceGraph *srcG = getSourceGraph();

	cost = 0;
	for (a = 0; a < srcG->getNumberOfVertices(); a++) {
		//if (partitioning[a] == 0) {
			for (nodeAdj = srcG->getAdjOfVertex(a); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
				if (partitioning[nodeAdj->w] != partitioning[a]) {
					//cout << "\n a: " << a << " edge: " << nodeAdj->w;
					cost += nodeAdj->edgeWeight * buildNetworkTopology(16, 4, partitioning[a], partitioning[nodeAdj->w]);
				}
			}
		//}
	}
	//cout << "\ncost: " << cost;
	return cost/2;*/
}

// João
int Comm::buildNetworkTopology(int n, int k, int i, int j){
    int dist;
               if(i == j)
                dist = 0;
            else
                if(2*i/k == 2*j/k)
                    dist = 1;
                else
                    if((2*i/k != 2*j/k) && (4*i/(k*k) == 4*j/(k*k)))
                        dist = 3;
                    else
                        if(4*i/(k*k) != 4*j/(k*k))
                            dist  = 5;

            return dist;    
}

/* For Comm, the cost is the amount of (external) communication 
	among partitions */
//double Comm::ComputeCost(Partitioning_t& p) {
double Comm::computeDiffCost(int *partitioning, 
								   unsigned node, 
								   unsigned k, 
								   unsigned originalNodePartition, 
								   unsigned original_k_Partition) {
	link nodeAdj, nodeAdjSrcG;
	SourceGraph *srcG = getSourceGraph();
	double testCost = currentCost;
	
	// TODO

	return testCost;
}

void Comm::reCost(int *partitioning) {
	currentCost = computeCost(partitioning);
}

/* Returns true if the partitioning is valid, false otherwise. */
bool Comm::validPartitioning(int *p) {
	SourceGraph *srcG = getSourceGraph();

	if (srcG->getEnableMemory() == 0) {
		return true;
	} else if (srcG->getEnableMemory() == 1) {
		int i, a, b, *targetMem, *sourceMem;

		// Calculate amount of memory needed for each partition	
		for (i = 0; i < getNumberOfPartitions(); i++) {
			validArray[i] = 0;
		}
		for (i = 0; i < srcG->getNumberOfVertices(); i++) {
			validArray[p[i]] += srcG->getMemory(i);
		}
		// end

		sourceMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
		if(sourceMem == NULL) {
			cout << "\n sourceMem could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}
		for (i = 0; i < getNumberOfPartitions(); i++) {
			sourceMem[i] = validArray[i];
		}
		targetMem = (int *) malloc(getNumberOfPartitions() * sizeof(int));
		if(targetMem == NULL) {
			cout << "\n targetMem could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}
		for (i = 0; i < getNumberOfPartitions(); i++) {
			targetMem[i] = target->getMemory(i);
		}

		mergeSortWithIndexes(sourceMem, 1, getNumberOfPartitions(), validIndexes);
		mergeSortWithIndexes(targetMem, 1, getNumberOfPartitions(), targetIndexes);

		// nonnaive implementation: sort validArray and targetMem in descending order
		// this makes the sorting not stable; to have a stable sort comment this block and use Source_graph->V/2 - 1 - a (or b) in FindMaxReductionInCost()
		for (a = 0; a < getNumberOfPartitions()/2; a++) {
			b = getNumberOfPartitions() - a - 1;

			int sortaux = sourceMem[a];
			sourceMem[a] = sourceMem[b];
			sourceMem[b] = sortaux;

			sortaux = validIndexes[a];
			validIndexes[a] = validIndexes[b];
			validIndexes[b] = sortaux;

			sortaux = targetMem[a];
			targetMem[a] = targetMem[b];
			targetMem[b] = sortaux;

			sortaux = targetIndexes[a];
			targetIndexes[a] = targetIndexes[b];
			targetIndexes[b] = sortaux;
		}

		for (i = 0; i < getNumberOfPartitions(); i++) {
			if (sourceMem[i] > targetMem[i]) {
				free(sourceMem);
				free(targetMem);
				return false;
			}
		}

		free(sourceMem);
		free(targetMem);
		return true;
	}
}

int Comm::getValidArray(int i) {
	return validArray[i];
}

int Comm::getValidIndex(int i) {
	return validIndexes[i];
}

void Comm::mergeSortWithIndexes(int *A, int p, int r, int *indexes){
	int j, **Ai=NULL;

	Ai = (int **) calloc(r, sizeof(int *));
	if (Ai == NULL) {
		printf("Ai could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < r; j++) {
		Ai[j] = (int *) calloc(2, sizeof(int));
		if (Ai[j] == NULL) {
			printf("Ai[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}

	for (j = 0; j < r; j++) {
		Ai[j][0] = j;
		Ai[j][1] = A[j];
	}

	/*for (j = p; j <= r; j++) {
		printf("\n%d %d", Ai[j - 1][0], Ai[j - 1][1]);
	}
	printf("\n");*/

	mergeSort(Ai, p, r);

	for (j = p; j <= r; j++) {
		A[j - 1] = Ai[j - 1][1];
		indexes[j - 1] = Ai[j - 1][0];
	}

	/*for (j = p; j <= r; j++) {
		printf("\n%d %d", Ai[j - 1][0], Ai[j - 1][1]);
	}
	printf("\n");*/

	for (j = 0; j < r; j++) {
		free(Ai[j]);
		Ai[j] = NULL;
	}	
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
}
