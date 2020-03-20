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
#include <sys/time.h>
#include "GenericPartitioningAlgorithm.h"

using namespace std;

GenericPartitioningAlgorithm::GenericPartitioningAlgorithm(int nVertices, int nPartitions, bool movingNodes, bool exchangingNodes) {

	setNumberOfVertices(nVertices);
	setNumberOfPartitions(nPartitions);
	setConsiderMovingNodes(movingNodes);
	setConsiderExchangingNodes(exchangingNodes);

	// if set is not used
	//releasedNodes = (NodeIndex_t *) malloc(getNumberOfVertices() * sizeof(NodeIndex_t));
	//releaseAllNodes();

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
}

void GenericPartitioningAlgorithm::setForcedInputSize(int size) {
	forcedInputSize = size;
}

int GenericPartitioningAlgorithm::getForcedInputSize() {
	return forcedInputSize;
}

void GenericPartitioningAlgorithm::setNumberOfPartitions(int nPartitions) {
	numberOfPartitions = nPartitions;
}

void GenericPartitioningAlgorithm::setNumberOfVertices(int nVertices) {
	numberOfVertices = nVertices;
}

int GenericPartitioningAlgorithm::getNumberOfVertices() {
	return numberOfVertices;
}

void GenericPartitioningAlgorithm::setConsiderMovingNodes(bool movingNodes) {
	considerMovingNodes = movingNodes;
}

void GenericPartitioningAlgorithm::setConsiderExchangingNodes(bool exchangingNodes) {
	considerExchangingNodes = exchangingNodes;
}

bool GenericPartitioningAlgorithm::getConsiderExchangingNodes() {
	return considerExchangingNodes;
}

bool GenericPartitioningAlgorithm::getConsiderMovingNodes() {
	return considerMovingNodes;
}

int GenericPartitioningAlgorithm::getCurrentPartitioning(int node) {
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
	PartitionIndex_t originalPartitioning = currentPartitioning[node];

	for (PartitionIndex_t p = 0; p < numberOfPartitions; p++) {
		if (p == originalPartitioning) continue;
		currentPartitioning[node] = p;

		if (!validPartitioning(currentPartitioning)) continue;

		//double c = computeDiffCost(currentPartitioning, node, 0);
		double c = computeCost(currentPartitioning);
		if (!foundMove || (c < bestC)) {
			foundMove = true;
			bestC = c;
			m.move_node(node, p);
		}
	}
	/* Restore current Partitioning. */
	currentPartitioning[node] = originalPartitioning;
	//cout << "\nfindBestSingleNodeMove: m.a: " << m.a << " cost: " << bestC << " foundMove: " << foundMove;
	return foundMove;
}

/* Check what is the best exchange for the node i (node) considering the current Partitioning */
bool GenericPartitioningAlgorithm::findBestNodeExchange(NodeIndex_t node, 
														Move_t &m, 
														double& bestC) {
	bool foundExchange = false;
	PartitionIndex_t originalNodePartition = currentPartitioning[node];
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
		if (currentPartitioning[k] == originalNodePartition) continue;

		/* Check whether exchanging k and i is profitable. */
		PartitionIndex_t original_k_Partition = currentPartitioning[k];
		/* Pretend we perform the exchange */ 
		currentPartitioning[node] = original_k_Partition;
		currentPartitioning[k] = originalNodePartition;

		if (validPartitioning(currentPartitioning)) {
			/* check the cost. */
			//double c = computeDiffCost(currentPartitioning, node, k, originalNodePartition, original_k_Partition);
			double c = computeCost(currentPartitioning);
			//cout << "\nnode: " << node << " k: " << k << " cost: " << c;
			if (!foundExchange || (c <= bestC)) {
				foundExchange = true;
				bestC = c;
				m.exchange_node(node, originalNodePartition, k, original_k_Partition);
			 }
		}

		/* Recover the original partitioning. */
		currentPartitioning[k] = original_k_Partition;
		currentPartitioning[node] = originalNodePartition;
		//}
	}

	/* Restore current Partitioning. */
	currentPartitioning[node] = originalNodePartition;
	//cout << "\nfindBestNodeExchange: m.a: " << m.a << " m.b: " << m.b << " cost: " << bestC << " foundExchange: " << foundExchange;
	return foundExchange;
}

bool GenericPartitioningAlgorithm::findBestMove(Move_t &move, double &cost) {

	bool found = false;
	std::set<NodeIndex_t>::iterator it;
	//int i;

	//for (i = 0; i < sourceG->getNumberOfVertices(); i++) {
		//if (releasedNodes[i]) {
	for (it = releasedNodes.begin(); it != releasedNodes.end(); ++it) {
		NodeIndex_t i = *it;

		// force input to stay in the same partition, if getInputSize() > 0
		if (i < getForcedInputSize()) continue;

		Move_t m = move; 
		double c = cost;

		if (considerMovingNodes) {
			if (findBestSingleNodeMove(i, m, c)) {
				if (!found || (c < cost)) {
					found = true;
					//bestMove = m;
					cost = c;
					move = m;
				}
			}
		}

		if (considerExchangingNodes){
			if (findBestNodeExchange(i, m, c)) {
				if (!found || (c <= cost)) {
					found = true;
					//bestMove = m;
					cost = c;
					move = m;
				}
			}
		}
		//}
	}
	//cout << "\n findBestMove: m.a: " << move.a << " m.b: " << move.b << " cost: " << cost << " found: " << found;
	return found;
}


void GenericPartitioningAlgorithm::performMove(int *partitioning, Move_t &m) {
	if (m.move_a) partitioning[m.a] = m.a_to;
	if (m.move_b) partitioning[m.b] = m.b_to;
}

void GenericPartitioningAlgorithm::lockNodes(Move_t &m) {
	if (m.move_a)
		//releasedNodes[m.a] = false;
		releasedNodes.erase(m.a);
	if (m.move_b)	
		//releasedNodes[m.b] = false;
		releasedNodes.erase(m.b);
}

void GenericPartitioningAlgorithm::releaseAllNodes() {
	/*int i;

	for (i = 0; i < getNumberOfVertices(); i++) {
		releasedNodes[i] = true;
	}*/

    /* Add all node indices to released_nodes. */ 
	for (NodeIndex_t i = 0; i < getNumberOfVertices(); i++) 
    	releasedNodes.insert(i);
}

/*bool GenericPartitioningAlgorithm::empty() {
	int i;

	for (i = 0; i < getNumberOfVertices(); i++) {
		if (releasedNodes[i] == true) {
			return false;
		}
	}
	return true;
}*/

int *GenericPartitioningAlgorithm::getInitialPartitioning() {
	return initialPartitioning;
}

int *GenericPartitioningAlgorithm::getBestPartitioning() {
	return bestPartitioning;
}

int GenericPartitioningAlgorithm::getNumberOfPartitions() {
	return numberOfPartitions;
}

void GenericPartitioningAlgorithm::printPartitioning(int *partitioning) {
	int i;
	
	cout << "\nPartitioning:\n";
	for (i = 0; i < getNumberOfVertices(); i++) {
		cout << partitioning[i] << " ";
	}
	cout << endl;
}

bool GenericPartitioningAlgorithm::run() {
	bool bestPartitionModified = true;
	double bestC = computeCost(bestPartitioning);
	int i;

	struct timeval tim;
	double initial, epoch_time;

	/* Cycle through epochs. */
	while(bestPartitionModified) { 
		gettimeofday(&tim, NULL);
		initial = tim.tv_sec+(tim.tv_usec/1000000.0);

		/* Epoch */
		bestPartitionModified = false;

		// updates currentCost 	
		reCost(bestPartitioning);

		for (i = 0; i < getNumberOfVertices(); i++) {
			currentPartitioning[i] = bestPartitioning[i];
		}

		/* Cycle through iterations */
		releaseAllNodes();
		//while (!empty()) {
		while (!releasedNodes.empty()) {
			/* Iteration */
			Move_t m;
			double cost;

			if (!findBestMove(m, cost)) 
				break; /* Could not find a valid move -- epoch done. */

			/* Lock nodes. */
			lockNodes(m);

			/* Perform move */
			performMove(currentPartitioning, m);

			//printPartitioning(currentPartitioning);

			// updates currentCost 	
			reCost(currentPartitioning);

			/* Save best */
			if (cost <= bestC) {
				if (cost < bestC) {
					bestPartitionModified = true;
				}
				bestC = cost;
				//bestPartitioning = currentPartitioning;
				for (i = 0; i < getNumberOfVertices(); i++) {
					bestPartitioning[i] = currentPartitioning[i];
				}
			}
			//cout << "\n  run: m.a: " << m.a << " m.b: " << m.b << " bestC: " << bestC << "\n";
		}
		gettimeofday(&tim, NULL);
		epoch_time = tim.tv_sec+(tim.tv_usec/1000000.0);
		cout << "\n  epoch: bestC: " << bestC << ", epoch time: " << epoch_time - initial << "\n";
	}
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
}
