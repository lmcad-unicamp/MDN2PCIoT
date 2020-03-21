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

#ifndef INFERENCERATE_H
#define INFERENCERATE_H

#include "SourceGraph.h"
#include "TargetGraph.h"
#include "GenericPartitioningAlgorithm.h" // generic KL class definition

class InferenceRate : public GenericPartitioningAlgorithm {

public:
	InferenceRate(SourceGraph *srcG, int nPartitions, bool movingNodes, bool exchangingNodes, TargetGraph *tgtG, int nVertices, int forcedInput, bool initPart, char *verbose, char *initPartFile, int numberOfThreads, bool multilevel, int *multilevelPart, int higherLevelNumberOfVertices);

	/* Returns the cost of partition the graph with partitioning p. */
	virtual double computeCost(const int *partitioning, bool openmp, bool nonredundantEdges);

	/* Returns true if the partitioning is valid, false otherwise. */
	virtual bool validPartitioning(const int *partitioning);
	virtual bool diffValidPartitioning(int *partitioning, unsigned node, unsigned k, unsigned originalNodePartition, unsigned original_k_Partition, bool singleOrSwap);

	void computeRedundantMemoryPerPartition(const int *p);
	void computeNonredundantMemoryPerPartition(const int *p);
	bool partitionsFitDevices(const int *p, int *sourceMem, bool openmp);
	void computeDiffNodeMemoryPerPartition(int *partitioning, unsigned vertex, unsigned originalVertexPartition, int *sourceMem);
	void computeDiffMemoryPerPartition(int *partitioning, unsigned node, unsigned k, unsigned originalNodePartition, unsigned original_k_Partition, int *sourceMem, bool singleOrSwap);

	int getValidArray(int i) const;
	int getValidIndex(int i) const;
	
	~InferenceRate();

private:
	// sort functions (for validPartitioning)
	void mergeSortWithIndexes(int *A, int p, int r, int *indexes);
	void mergeSort(int **A, int p, int r);
	// auxiliary function that combines two sorted halves to be sorted
	void merge(int **A, int p, int q, int r);

	void setSourceGraph(SourceGraph *srcG);
	void setInitialPartitionings(int forcedInput, bool initPart, char *initPartFile, bool multilevel, int *multilevelPart, int higherLevelNumberOfVertices);
	void seed_rng(void);

	const SourceGraph *sourceG;	// could be const
	TargetGraph *target;	// const

	// validArray contains memory needed for each partition
	int *validArray;
	// in order to calculate if some partitioning is memory valid
	int *validIndexes;
	int **validIndexesPerThread; // OpenMP
	int *targetIndexes;

	// diffValidPartitioning
	//int *sourceMem;
	int *sortedTargetMem;

	// eliminates redundant memory by accounting for shared memory only once per partition
	bool **sharedMemoryPerPartition;

	// in order to just update cost, thus calculating it faster	
	double currentCost;

	// eliminates redundant edges (edges that come from the same vertex and go to the same partition)
	SourceGraph *partitionGraph;
	SourceGraph **partitionGraphPerThread; // OpenMP

	// stores the computational power needed by each partition
	double *partitionComputationalWeight;
	// stores the inference rate of each connection between each device
	double **connectionMatrixInferenceRate;
	// stores the inference rate of each device
	double *deviceInferenceRate;

	// stores the computational power needed by each partition
	double **partitionComputationalWeightPerThread; // OpenMP
	// stores the inference rate of each connection between each device
	double ***connectionMatrixInferenceRatePerThread; // OpenMP
	// stores the inference rate of each device
	double **deviceInferenceRatePerThread; // OpenMP
};

#endif // INFERENCERATE_H
