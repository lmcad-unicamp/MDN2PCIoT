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

#ifndef GPA
#define GPA

#include <set>
#include <vector>
#include "SourceGraph.h"

typedef unsigned NodeIndex_t;

typedef unsigned PartitionIndex_t;

typedef void graph_t;

using namespace std;

/* Partitioning is an array of Partition indices. 
   The i'th array position holds the partition index for the i'th node. */
typedef std::vector<PartitionIndex_t> Partitioning_t;

struct Move_t {

  /* Default constructor. */
  Move_t() : move_a(false), move_b(false) {}

  NodeIndex_t a;	
  PartitionIndex_t a_to;
  bool move_a;

  NodeIndex_t b;
  PartitionIndex_t b_to;
  bool move_b;

  void move_node(NodeIndex_t n, PartitionIndex_t p) {
    move_a = true;
    move_b = false;
    a = n;
    a_to = p;
  }

  void exchange_node(NodeIndex_t n1, PartitionIndex_t p1, NodeIndex_t n2, PartitionIndex_t p2) {
    move_a = true;
    move_b = true;
    a = n1;
    a_to = p2;
    b = n2;
    b_to = p1;
  }

};

class GenericPartitioningAlgorithm {

public:
	GenericPartitioningAlgorithm(int nVertices, int nPartitions, bool movingNodes, bool exchangingNodes, int numberOfThreads);
	
	void setInitialPartitioning();
	void printPartitioning(const int *partitioning) const;

	bool run(char *instance, bool multilevel, int samee, int sameee);

	/* Returns the cost of partition the graph with partitioning p. */
	virtual double computeCost(const int *partitioning, bool openmp, bool nonredundantEdges) = 0;
	//virtual double computeDiffCost(int *partitioning, unsigned node, unsigned k, unsigned originalNodePartition, unsigned original_k_Partition, bool openmp, bool nonredundantEdges, bool exchange) = 0;
	// updates currentCost
	//virtual void reCost(int *partitioning) = 0;

	/* Returns true if the partitioning is valid, false otherwise. */
	virtual bool validPartitioning(const int *partitioning) = 0;
	virtual bool diffValidPartitioning(int *partitioning, unsigned node, unsigned k, unsigned originalNodePartition, unsigned original_k_Partition, bool singleOrSwap) = 0;
	virtual int getValidArray(int i) const = 0;

	const int *getInitialPartitioning() const;
	const int *getBestPartitioning() const; 
	int getCurrentPartitioning(int node) const;

	int getVertexOfBestPartitioning(int u) const;

	bool getConsiderExchangingNodes() const {
		return considerExchangingNodes;
	}
	bool getConsiderMovingNodes() const {
		return considerMovingNodes;
	}
	int getNumberOfVertices() const {
		return numberOfVertices;
	}
	int getNumberOfPartitions() const {
		return numberOfPartitions;
	}

	void setVertexInitialPartitionings(int v, int partition);

	// to force input to be in the same partition
	void setForcedInputSize(int size);	
	int getForcedInputSize() const {
		return forcedInputSize;
	}

	// to force input to be in the same partition
	int getNumberOfThreads() const {
		return numberOfThreads;
	}
	int getThreadsForPartitions() const {
		return threadsForPartitions;
	}

	void setLockedArraySize(int lsize);
	void setLockedArray(int position, int node);
	int getLockedArray(int position) const {
		return lockedArray[position];
	}

	~GenericPartitioningAlgorithm();

private:
	bool findBestSingleNodeMove(NodeIndex_t node, Move_t &m, double &bestC);
	bool findBestNodeExchange(NodeIndex_t node, Move_t &m, double &bestC, int samee);
	bool findBestMove(Move_t &m, double &cost, int &thread, int samee, int sameee);
	
	// If a member function is defined in the body of a class definition, the member function is implicitly declared inline.
	void performMove(int *partitioning, Move_t &m) {
		if (m.move_a) partitioning[m.a] = m.a_to;
		if (m.move_b) partitioning[m.b] = m.b_to;
	}

	void lockNodes(Move_t &m)  {
		if (m.move_a) {
			releasedVertices[m.a] = false;
			releasedNodes.erase(m.a);
		}
		if (m.move_b) {
			releasedVertices[m.b] = false;
			releasedNodes.erase(m.b);
		}
	}	

	//bool empty(); (if set is not used)
	void releaseAllNodes() {

#		pragma omp parallel for num_threads(getNumberOfThreads())
		for (int i = 0; i < getNumberOfVertices(); i++) {
			releasedVertices[i] = true;
		}

	    /* Add all node indices to released_nodes. */
		for (NodeIndex_t i = 0; i < getNumberOfVertices(); i++) 
    			releasedNodes.insert(i);
	} 

	void lockInnerNode(int v) {
		releasedVertices[v] = false;
		releasedNodes.erase(v);
	}  

	void copyCurrentPartForThreads();
	void copyBackCurrentPartForThreads(int thread);

	// Multilevel refinement
	void lockInnerNodes();

	// data members:
	int *initialPartitioning;
	int *currentPartitioning;
	int *bestPartitioning;

	int **currentPartPerThread;

	const int numberOfPartitions;		// const

	std::set<NodeIndex_t> releasedNodes;
	bool *releasedVertices; // same as releasedNodes, for OpenMP

	const bool considerMovingNodes;	// const
	const bool considerExchangingNodes;	// const

	double bestCost;

	//graph_t* sourceGraph;
	const int numberOfVertices;	// const

	// to force input to be in the same partition
	int forcedInputSize; // const

	// OpenMP
	int numberOfThreads;
	int threadsForPartitions;

	// Multilevel refinement
	int *lockedArray;
	int lockedArraySize;

	// 0: comm; 1: iR
	bool objFunction;
};

#endif // GPA
