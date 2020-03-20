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

#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "GenericPartitioningAlgorithm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

GenericPartitioningAlgorithm::GenericPartitioningAlgorithm(int nVertices, int nPartitions, bool movingNodes, bool exchangingNodes, int nThreads) : numberOfPartitions(nPartitions), considerMovingNodes(movingNodes), considerExchangingNodes(exchangingNodes), numberOfVertices(nVertices), numberOfThreads(nThreads) {

	// if set is not used
	releasedVertices = (bool *) malloc(getNumberOfVertices() * sizeof(bool));
	releaseAllNodes();

	initialPartitioning = (int *) malloc(getNumberOfVertices() * sizeof(int));
	if(initialPartitioning == NULL) {
		cout << "\n initialPartitioning could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}

	currentPartitioning = (int *) malloc(getNumberOfVertices() * sizeof(int));
	if(currentPartitioning == NULL) {
		cout << "\n currentPartitioning could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}

	bestPartitioning = (int *) malloc(getNumberOfVertices() * sizeof(int));
	if(bestPartitioning == NULL) {
		cout << "\n bestPartitioning could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}

	currentPartPerThread = (int **) malloc(getNumberOfThreads() * sizeof(int *));
	if(currentPartPerThread == NULL) {
		cout << "\n currentPartPerThread could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < getNumberOfThreads(); i++) {
		currentPartPerThread[i] = (int *) malloc(getNumberOfVertices() * sizeof(int));
		if (currentPartPerThread[i] == NULL) {
			cout << "\n currentPartPerThread[i] could not be allocated! \n";
			exit(EXIT_FAILURE);	
		}		
	}
}

void GenericPartitioningAlgorithm::setForcedInputSize(int size) {
	forcedInputSize = size;
}

int GenericPartitioningAlgorithm::getCurrentPartitioning(int node) const {
	return currentPartitioning[node];	
}

void GenericPartitioningAlgorithm::setVertexInitialPartitionings(int v, int partition) {
	initialPartitioning[v] = partition;
	currentPartitioning[v] = partition;
	bestPartitioning[v] = partition;
}

/* Check what is the best move move for node i (node) considering the current Partitioning */
bool GenericPartitioningAlgorithm::findBestSingleNodeMove(NodeIndex_t node, 
														  Move_t &m, 
													      double &bestC) {
	bool foundMove = false;
	PartitionIndex_t originalPartitioning = currentPartPerThread[omp_get_thread_num()][node];
	//PartitionIndex_t originalPartitioning = currentPartitioning[node];

	for (PartitionIndex_t p = 0; p < numberOfPartitions; p++) {
		if (p == originalPartitioning) continue;

		currentPartPerThread[omp_get_thread_num()][node] = p;
		//currentPartitioning[node] = p;

		//if (!validPartitioning(currentPartitioning)) continue;
		if (!diffValidPartitioning(currentPartPerThread[omp_get_thread_num()], node, 0, originalPartitioning, 0, 0)) continue;

		//double c = computeDiffCost(currentPartitioning, node, 0);
		double c = computeCost(currentPartPerThread[omp_get_thread_num()], true);
		//double c = computeCost(currentPartitioning);
		//cout << "\nThread[" << omp_get_thread_num() << "]: " << c;
		if (!foundMove || (c < bestC)) {
			foundMove = true;
			bestC = c;
			m.move_node(node, p);
		}
	}
	/* Restore current Partitioning. */
	currentPartPerThread[omp_get_thread_num()][node] = originalPartitioning;
	//currentPartitioning[node] = originalPartitioning;
	/*if (omp_get_thread_num() == 2)
		cout << "\nfindBestSingleNodeMove: m.a: " << m.a << " cost: " << bestC << " foundMove: " << foundMove;*/
	return foundMove;
}

/* Check what is the best exchange for the node i (node) considering the current Partitioning */
bool GenericPartitioningAlgorithm::findBestNodeExchange(NodeIndex_t node, 
														Move_t &m, 
														double& bestC) {
	bool foundExchange = false;
	//PartitionIndex_t originalNodePartition = currentPartitioning[node];
	PartitionIndex_t originalNodePartition = currentPartPerThread[omp_get_thread_num()][node];
	//int i;

	//cout << "\nbestC: " << bestC;

	std::set<NodeIndex_t>::iterator it;
	//for (i = 0; i < sourceG->getNumberOfVertices(); i++) {
		//if (releasedNodes[i]) {
			//int k = i;
	for (it = releasedNodes.begin(); it != releasedNodes.end(); ++it) {
		NodeIndex_t k = *it;

		/* Skip in case i and k are the same. */
		if (node == k) continue;

		// force input to stay in the same partition, if getInputSize() > 0
		if (k < getForcedInputSize()) continue;

		/* Skip in case i and k are in the same partition. */
		if (currentPartPerThread[omp_get_thread_num()][k] == originalNodePartition) continue;

		/* Check whether exchanging k and i is profitable. */
		PartitionIndex_t original_k_Partition = currentPartPerThread[omp_get_thread_num()][k];
		/* Pretend we perform the exchange */ 
		currentPartPerThread[omp_get_thread_num()][node] = original_k_Partition;
		currentPartPerThread[omp_get_thread_num()][k] = originalNodePartition;

		//if (validPartitioning(currentPartitioning)) {
		if (diffValidPartitioning(currentPartPerThread[omp_get_thread_num()], node, k, originalNodePartition, original_k_Partition, 1)) {
			/* check the cost. */
			//double c = computeDiffCost(currentPartitioning, node, k, originalNodePartition, original_k_Partition);
			double c = computeCost(currentPartPerThread[omp_get_thread_num()], true);
			//cout << "\nnode: " << node << " k: " << k << " cost: " << c;
			if (!foundExchange || (c <= bestC)) {
				foundExchange = true;
				bestC = c;
				m.exchange_node(node, originalNodePartition, k, original_k_Partition);
			 }
		}

		/* Recover the original partitioning. */
		currentPartPerThread[omp_get_thread_num()][k] = original_k_Partition;
		currentPartPerThread[omp_get_thread_num()][node] = originalNodePartition;
		//}
	}

	/* Restore current Partitioning. */
	currentPartPerThread[omp_get_thread_num()][node] = originalNodePartition;
	//cout << "\nfindBestNodeExchange: m.a: " << m.a << " m.b: " << m.b << " cost: " << bestC << " foundExchange: " << foundExchange;
	return foundExchange;
}

bool GenericPartitioningAlgorithm::findBestMove(Move_t &move, double &cost, int &thread) {
	bool found = false;
	std::set<NodeIndex_t>::iterator it;
	NodeIndex_t i;

#	pragma omp parallel for num_threads(getNumberOfThreads()) reduction(||: found) schedule(static, 1)
	for (i = 0; i < getNumberOfVertices(); i++) {
		if (!releasedVertices[i]) continue;
	/*for (it = releasedNodes.begin(); it != releasedNodes.end(); ++it) {
		NodeIndex_t i = *it;*/

		// force input to stay in the same partition, if getInputSize() > 0
		if (i < getForcedInputSize()) continue;

		Move_t m = move; 
		double c = cost;
		int t = omp_get_thread_num();

		if (considerMovingNodes) {
			if (findBestSingleNodeMove(i, m, c)) {
#				pragma omp critical
				if (!found || (c < cost)) {
					found = true;
					//bestMove = m;
					cost = c;
					move = m;
					thread = t;
				}
			}
		}

		if (considerExchangingNodes){
			if (findBestNodeExchange(i, m, c)) {
#				pragma omp critical
				if (!found || (c <= cost)) {
					found = true;
					//bestMove = m;
					cost = c;
					move = m;
					thread = t;
				}
			}
		}
	}
	//cout << "\n findBestMove: m.a: " << move.a << " m.b: " << move.b << " cost: " << cost << " found: " << found;
	return found;
}

const int *GenericPartitioningAlgorithm::getInitialPartitioning() const {
	return initialPartitioning;
}

const int *GenericPartitioningAlgorithm::getBestPartitioning() const {
	return bestPartitioning;
}

void GenericPartitioningAlgorithm::printPartitioning(const int *partitioning) const {
	int i;
	
	cout << "\nPartitioning:\n";
	for (i = 0; i < getNumberOfVertices(); i++) {
		cout << partitioning[i] << " ";
	}
	cout << endl;
}

void GenericPartitioningAlgorithm::copyCurrentPartForThreads(){
#	pragma omp parallel for num_threads(getNumberOfThreads())
	for (int j = 0; j < getNumberOfThreads(); j++)
		for (int i = 0; i < getNumberOfVertices(); i++)
			currentPartPerThread[j][i] = currentPartitioning[i];
}

void GenericPartitioningAlgorithm::copyBackCurrentPartForThreads(int thread){
#	pragma omp parallel for num_threads(getNumberOfThreads())
	for (int i = 0; i < getNumberOfVertices(); i++)
		currentPartitioning[i] = currentPartPerThread[thread][i];
}

bool GenericPartitioningAlgorithm::run(char *instance) {
	bool bestPartitionModified = true;
	double bestC = computeCost(bestPartitioning, false);
	int i;
	FILE *wcost = NULL, *currentPart = NULL, *epochPart = NULL;
	char currentPartFileName[1000], epochPartFileName[1000];

	wcost = fopen(instance, "w");

	strcpy(currentPartFileName, instance);
	strcat(currentPartFileName, "-currentPart");

	strcpy(epochPartFileName, instance);
	strcat(epochPartFileName, "-epochPart");

	struct timeval tim;
	double initial, epoch_time;

	if (bestC >=0) {
		fprintf(wcost, " %f", bestC);
		cout << "\n" << bestC;
	} else {
		fprintf(wcost, " %f", -bestC);
		cout << "\n" << -bestC;
	}

	/* Cycle through epochs. */
	while(bestPartitionModified) { 
		gettimeofday(&tim, NULL);
		initial = tim.tv_sec+(tim.tv_usec/1000000.0);

		/* Epoch */
		bestPartitionModified = false;

		// updates currentCost 	
		//reCost(bestPartitioning);

		// updates validArray
		validPartitioning(bestPartitioning);

#		pragma omp parallel for num_threads(getNumberOfThreads())
		for (i = 0; i < getNumberOfVertices(); i++) {
			currentPartitioning[i] = bestPartitioning[i];
		}

		/* Cycle through iterations */
		releaseAllNodes();
		//while (!empty()) {
		while (!releasedNodes.empty()) {
			/* Iteration */
			Move_t m;
			double cost = bestC;
			int thread;

			copyCurrentPartForThreads();

			if (!findBestMove(m, cost, thread)) 
				break; /* Could not find a valid move -- epoch done. */

			//cout << "\nt=" << thread << " c=" << bestC << " m=" << m.a << " c=" << cost; 

			//copyBackCurrentPartForThreads(thread);

			/* Lock nodes. */
			lockNodes(m);

			/* Perform move */
			performMove(currentPartitioning, m);

			//printPartitioning(currentPartitioning);

			// updates currentCost 	
			//reCost(currentPartitioning);

			// updates validArray
			validPartitioning(currentPartitioning);

			/* Save best */
			if (cost <= bestC) {
				if (cost < bestC) {
					bestPartitionModified = true;
				}
				bestC = cost;
				currentPart = fopen(currentPartFileName, "w");
				//bestPartitioning = currentPartitioning;
				for (i = 0; i < getNumberOfVertices(); i++) {
					bestPartitioning[i] = currentPartitioning[i];
					fprintf(currentPart, "%d ", bestPartitioning[i]);
				}
				fclose(currentPart);
				currentPart = NULL;
			}
			//cout << "\n  run: m.a: " << m.a << " m.b: " << m.b << " bestC: " << bestC << "\n";
		}

		epochPart = fopen(epochPartFileName, "w");
		for (i = 0; i < getNumberOfVertices(); i++) {
			fprintf(epochPart, "%d ", bestPartitioning[i]);
		}
		fclose(epochPart);
		epochPart = NULL;

		gettimeofday(&tim, NULL);
		epoch_time = tim.tv_sec+(tim.tv_usec/1000000.0);
		//cout << "\n  epoch: bestC: " << bestC << ", epoch time: " << (epoch_time - initial) / 3600 << "\n";
		if (bestC >=0) {
			fprintf(wcost, " %f", bestC);
			cout << "\n" << bestC << ", epoch time: " << (epoch_time - initial) / 3600 << endl;
		} else {
			fprintf(wcost, " %f", -bestC);
			cout << "\n" << -bestC << ", epoch time: " << (epoch_time - initial) / 3600 << endl;
		}
	}
	fprintf(wcost, "\n");
	fclose(wcost);
	wcost = NULL;
	return bestPartitionModified;
}

GenericPartitioningAlgorithm::~GenericPartitioningAlgorithm() {
	if (initialPartitioning != NULL) {
		free(initialPartitioning);
		initialPartitioning = NULL;
	}
	if (bestPartitioning != NULL) {
		free(bestPartitioning);
		bestPartitioning = NULL;
	}
	if (currentPartitioning != NULL) {
		free(currentPartitioning);
		currentPartitioning = NULL;
	}

	if (currentPartPerThread != NULL) {
		for (int i = 0; i < getNumberOfThreads(); i++) {
			if (currentPartPerThread[i] != NULL) {
				free(currentPartPerThread[i]);
			}
		}
		free(currentPartPerThread);
		currentPartPerThread = NULL;
	}
}
