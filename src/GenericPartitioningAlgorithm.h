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
	GenericPartitioningAlgorithm(int nVertices, int nPartitions, bool movingNodes, bool exchangingNodes);
	
	void setInitialPartitioning();
	void printPartitioning(int *partitioning);

	bool run();

	/* Returns the cost of partition the graph with partitioning p. */
	virtual double computeCost(int *partitioning) = 0;
	virtual double computeDiffCost(int *partitioning, unsigned node, unsigned k, unsigned originalNodePartition, unsigned original_k_Partition) = 0;
	// updates currentCost
	virtual void reCost(int *partitioning) = 0;

	/* Returns true if the partitioning is valid, false otherwise. */
	virtual bool validPartitioning(int *partitioning) = 0;

	int *getInitialPartitioning();
	int *getBestPartitioning();
	int getCurrentPartitioning(int node);

	int getNumberOfPartitions();
	int getNumberOfVertices();
	bool getConsiderExchangingNodes();
	bool getConsiderMovingNodes();

	void setVertexInitialPartitionings(int v, int partition);

	// to force input to be in the same partition
	int getForcedInputSize();
	void setForcedInputSize(int size);

	~GenericPartitioningAlgorithm();

private:
	bool findBestSingleNodeMove(NodeIndex_t node, Move_t &m, double &bestC);
	bool findBestNodeExchange(NodeIndex_t node, Move_t &m, double &bestC);
	bool findBestMove(Move_t &m, double &cost);
	void performMove(int *partitioning, Move_t &m);

	void lockNodes(Move_t &m);
	void releaseAllNodes();    //bool empty(); (if set is not used)

	void setNumberOfPartitions(int nPartitions);
	void setConsiderMovingNodes(bool movingNodes);
	void setConsiderExchangingNodes(bool exchangingNodes);
	void setNumberOfVertices(int nVertices);

	// data members:
	int *initialPartitioning;
	int *currentPartitioning;
	int *bestPartitioning;

	int numberOfPartitions;

	std::set<NodeIndex_t> releasedNodes;

	bool considerMovingNodes;
	bool considerExchangingNodes;

	double bestCost;

	//graph_t* sourceGraph;
	int numberOfVertices;

	// to force input to be in the same partition
	int forcedInputSize;
};

#endif // GPA
