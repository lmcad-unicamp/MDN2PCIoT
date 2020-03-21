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

#ifndef COARSENGRAPH_H
#define COARSENGRAPH_H

#include "SourceGraph.h"
#include "TargetGraph.h"

class CoarsenGraph {

public:
	CoarsenGraph(SourceGraph *sourceG, TargetGraph *tgtG, bool verb, const int *partitioning, int definedNumberOfCoarsenedGraphs);

	SourceGraph* getCoarsenedGraph(int level) const {
		return sourceCoarsenedGraph[level];
	}

	int getNumberOfCoarsenedGraphs() const {
		return numberOfCoarsenedGraphs;
	}

	int *getDefinedPartitioning() const {
		return p;
	}

	void seed_rng();

	~CoarsenGraph();

private:
	// Multilevel
	void coarsenSourceGraph();
	void sortDegree(int c);

	// sort functions
	void mergeSortWithIndexes(int *A, int p, int r, int *indexes);
	void mergeSort(int **A, int p, int r);
	// auxiliary function that combines two sorted halves to be sorted
	void merge(int **A, int p, int q, int r);

	//SourceGraph *sourceG1;		// Source graph level 1
	SourceGraph **sourceCoarsenedGraph;
	TargetGraph *target;
	int numberOfCoarsenedGraphs;
	double granularityBalanceFactor;
	int *sortedDegree;
	int *sortedVerticesByDegree;
	bool verb;

	// defined partitioning
	int *p;
};

#endif // COARSENGRAPH_H
