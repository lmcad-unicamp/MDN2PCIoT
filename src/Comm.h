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

#ifndef COMM_H
#define COMM_H

#include "SourceGraph.h"
#include "TargetGraph.h"
#include "GenericPartitioningAlgorithm.h" // generic KL class definition

class Comm : public GenericPartitioningAlgorithm {

public:
	Comm(SourceGraph *srcG, int nPartitions, bool movingNodes, bool exchangingNodes, TargetGraph *tgtG, int nVertices, int forcedInput, bool initPart);

	/* Returns the cost of partition the graph with partitioning p. */
	virtual double computeCost(int *partitioning);
	virtual double computeDiffCost(int *partitioning, unsigned node, unsigned k, unsigned originalNodePartition, unsigned original_k_Partition);
	virtual void reCost(int *partitioning);

	/* Returns true if the partitioning is valid, false otherwise. */
	virtual bool validPartitioning(int *partitioning);

	int getValidArray(int i);
	int getValidIndex(int i);

	int buildNetworkTopology(int n, int k, int i, int j);

	SourceGraph *getSourceGraph();
	~Comm();

private:
	// sort functions (for validPartitioning)
	void mergeSortWithIndexes(int *A, int p, int r, int *indexes);
	void mergeSort(int **A, int p, int r);
	// auxiliary function that combines two sorted halves to be sorted
	void merge(int **A, int p, int q, int r);

	void setSourceGraph(SourceGraph *srcG);
	void setInitialPartitionings(int forcedInput, bool initPart);

	SourceGraph *sourceG;
	TargetGraph *target;

	// validArray contains memory needed for each partition
	int *validArray;
	// in order to calculate if some partitioning is memory valid
	int *validIndexes;
	int *targetIndexes;

	// in order to just update cost, thus calculating it faster	
	double currentCost;

	SourceGraph *partitionGraph;
};

#endif // COMM_H
