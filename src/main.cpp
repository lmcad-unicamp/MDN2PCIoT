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
#include <string.h>
#include <stdlib.h>
#include "ReadInputFiles.h"
#include "InferenceRate.h"
#include "Comm.h"
#include "CoarsenGraph.h"

using namespace std;

int main(int argc, char *argv[]) {

	bool singleMove = true;
	bool swaps = true;
	bool initPart = false;
	int forcedInput = 0;

	if (strcmp(argv[7], "0") == 0) {
		forcedInput = 0; // 0 to do not force input
	} else if (strcmp(argv[7], "256") == 0) {
		forcedInput = 256; // 256 to force input (SBAC-PAD-2018, 2x2)
	} else if (strcmp(argv[7], "1024") == 0) {
		forcedInput = 1024; // 1024 to LeNet 1x1
	}

	if ((strcmp(argv[4], "0")) == 0) {
		initPart = false;  //false: random initial partitioning
							// true: read initial partitioning from initPart.txt
	} else {
		initPart = true;
	}

	// OpenMP
	int numberOfThreads = strtol(argv[8], NULL, 10);

	ReadInputFiles r(argv[1], argv[2], numberOfThreads); // automatic object

	SourceGraph g(*r.getSourceGraph());
	if ((strcmp(argv[6], "verb")) == 0) {
		g.printGraph();
	} else {
		g.printGraphHeader();
	}

	// Adjacency matrix
	/*for (int i = 0; i < g.getNumberOfVertices(); i++) {
		for (int j = 0; j < g.getNumberOfVertices(); j++)
			cout << g.getAdjMatrixEdgeWeight(i, j) << " ";
		cout << "\n";
	}*/

	TargetGraph t(*r.getTargetGraph());

	// multilevel graph coarsening
	bool verb = false;
	if ((strcmp(argv[6], "verb")) == 0) {
		verb = true;
	}

	int partitioning[g.getNumberOfVertices()];

	Comm initialP(&g, r.getTargetNumberOfVertices(), singleMove, swaps, &t, g.getNumberOfVertices(), forcedInput, initPart, argv[6], argv[4], numberOfThreads, false, partitioning, g.getNumberOfVertices(), false, partitioning, 0);

	//initialP.printPartitioning(initialP.getInitialPartitioning());
	
	CoarsenGraph coarsenGraph(&g, &t, verb, initialP.getInitialPartitioning(), strtol(argv[11], NULL, 10));

	if ((strcmp(argv[6], "verb")) == 0) {
		t.printGraphHeader();
		//t.printGraph();
	}

	cout << "\nsourceGraphFile: " << argv[1] << "\ntargetGraphFile: " << argv[2];
	cout << "\nSingleMoves: " << singleMove << ", swaps: " << swaps << ", locked input: " << forcedInput << ", initPart: " << initPart << ", number of threads: " << numberOfThreads;
	cout << "\ninitPartFileName: " << argv[4];

	int same, samee;
	if (coarsenGraph.getCoarsenedGraph(0)->getNumberOfVertices() < 700) {
		same = 50;
		samee = strtol(argv[10], NULL, 10);
	} else if (coarsenGraph.getCoarsenedGraph(0)->getNumberOfVertices() < 2500) {
		same = 50;
		samee = strtol(argv[10], NULL, 10);
	} else {
		same = strtol(argv[9], NULL, 10); // Before: 50;
		samee = strtol(argv[10], NULL, 10);
	}

	cout << "\nSame: " << same << ", Samee: " << samee;

	for (int c = coarsenGraph.getNumberOfCoarsenedGraphs() - 1; c >= 0; c--) {

		int nVertices = g.getNumberOfVertices();
		// random initial partitioning
		bool multilevel = false;
		// turn off pmatch: defPart = false
		bool defPart = false;
		if (c < coarsenGraph.getNumberOfCoarsenedGraphs() - 1) {
			nVertices = coarsenGraph.getCoarsenedGraph(c + 1)->getNumberOfVertices();
			multilevel = true;
			defPart = false;
		}

		// Skip larger graphs
		/*if (g.getNumberOfVertices() < 1200) {
			if (r.getTargetNumberOfVertices() > 50 && coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() > 400)// || coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() > 1600)
				continue;
		} else if (g.getNumberOfVertices() < 2500) {
			if (r.getTargetNumberOfVertices() > 50 && coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() > 800)
				continue;
		}*/ /*else {
			if (coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() > 4800)
				//continue;
		}*/

		// ----------------------//
		// ----------------------//
		// ----- Comm -----//
		if ((strcmp(argv[5], "MDN2PCIoTcomm")) == 0) {
			Comm Comm(coarsenGraph.getCoarsenedGraph(c), r.getTargetNumberOfVertices(), singleMove, swaps, &t, coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices(), forcedInput, initPart, argv[6], argv[4], numberOfThreads, multilevel, partitioning, nVertices, defPart, coarsenGraph.getDefinedPartitioning(), 0);

			//cout << "\n nVertices: " << nVertices;

			cout << "\nMDN2PCIoTcomm c: " << c - 1 << " numberOfVertices: " << coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() << " numberOfEdges: " << coarsenGraph.getCoarsenedGraph(c)->getNumberOfEdges();
			//cout << "\ninitPart: " << ;

			cout << "\n\nInitial partitioning: \n";
			Comm.printPartitioning(Comm.getInitialPartitioning());
			//cout << "initialCost (nonredundantEdges: " << Comm.computeCost(fullGraph.getInitialPartitioning(), false, true) << "\n";
			cout << "initialCost: " << Comm.computeCost(Comm.getInitialPartitioning(), false, true) << "\n";
			if (c == 0)
				cout << "initialCost (Nonredundant edges): " << Comm.computeCost(Comm.getInitialPartitioning(), false, true) << "\n";
			int memoryValid = Comm.validPartitioning(Comm.getInitialPartitioning());
			for (int i = 0; i < Comm.getNumberOfPartitions(); i++) {
				cout << "mem[" << i << "]: " << Comm.getValidArrayPerThread(0, i) << " ";
			}
			cout << "\nmemoryValid: " << memoryValid << "\n\n" << endl;

			// !KLP1
			if ((r.getTargetNumberOfVertices() < 50 && coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() < 700) || 
			 r.getTargetNumberOfVertices() < 12 ) //&& coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() < 2500)
				multilevel = false;

			Comm.run(argv[3], multilevel, same, samee);

			cout << "\nBest partitioning found:";
			Comm.printPartitioning(Comm.getBestPartitioning());
			cout << "bestCost: " << Comm.computeCost(Comm.getBestPartitioning(), false, true) << "\n";
			if (c == 0)
				cout << "bestCost (Nonredundant edges): " << Comm.computeCost(Comm.getBestPartitioning(), false, true) << "\n";
			memoryValid = Comm.validPartitioning(Comm.getBestPartitioning());
			for (int i = 0; i < Comm.getNumberOfPartitions(); i++) {
				cout << "mem[" << i << "]: " << Comm.getValidArrayPerThread(0, i) << " ";
			}
			cout << "\nmemoryValid: " << memoryValid << "\n\n";

			// for next iteration of multilevel partitioning
			for (int i = 0; i < coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices(); i++)
				partitioning[i] = Comm.getVertexOfBestPartitioning(i);

		// ----------------------------//
		// ----------------------------//
		// ----- InferenceRate -----//
		} else if ((strcmp(argv[5], "MDN2PCIoTiR")) == 0) {
			InferenceRate InferenceRate(coarsenGraph.getCoarsenedGraph(c), r.getTargetNumberOfVertices(), singleMove, swaps, &t, coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices(), forcedInput, initPart, argv[6], argv[4], numberOfThreads, multilevel, partitioning, nVertices, 1);

			cout << "\nMDN2PCIoTiR c: " << c - 1 << " numberOfVertices: " << coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() << " numberOfEdges: " << coarsenGraph.getCoarsenedGraph(c)->getNumberOfEdges();

			cout << "\n\nInitial partitioning:";

			InferenceRate.printPartitioning(InferenceRate.getInitialPartitioning());
			cout << "initialCost: " << InferenceRate.computeCost(InferenceRate.getInitialPartitioning(), false, true) << "\n";
			if (c == 0)
				cout << "initialCost (Nonredundant edges): " << InferenceRate.computeCost(InferenceRate.getInitialPartitioning(), false, true) << "\n";
			int memoryValid = InferenceRate.validPartitioning(InferenceRate.getInitialPartitioning());
			for (int i = 0; i < InferenceRate.getNumberOfPartitions(); i++) {
				cout << "mem[" << i << "]: " << InferenceRate.getValidArrayPerThread(0, i) << " ";
			}
			cout << "\nmemoryValid: " << memoryValid << "\n\n" << endl;
	
			// !KLP1
			if ((r.getTargetNumberOfVertices() < 50 && coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices() < 700) || 
			 r.getTargetNumberOfVertices() < 12)
				multilevel = false;

			InferenceRate.run(argv[3], multilevel, same, samee);

			cout << "\nBest partitioning found:";
			InferenceRate.printPartitioning(InferenceRate.getBestPartitioning());
			cout << "bestCost: " << InferenceRate.computeCost(InferenceRate.getBestPartitioning(), false, true) << "\n";
			if (c == 0)
				cout << "bestCost (Nonredundant edges): " << InferenceRate.computeCost(InferenceRate.getBestPartitioning(), false, true) << "\n";
			memoryValid = InferenceRate.validPartitioning(InferenceRate.getBestPartitioning());
			for (int i = 0; i < InferenceRate.getNumberOfPartitions(); i++) {
				cout << "mem[" << i << "]: " << InferenceRate.getValidArrayPerThread(0, i) << " ";
			}
			cout << "\nmemoryValid: " << memoryValid << "\n\n";

			// for next iteration of multilevel partitioning
			for (int i = 0; i < coarsenGraph.getCoarsenedGraph(c)->getNumberOfVertices(); i++)
				partitioning[i] = InferenceRate.getVertexOfBestPartitioning(i);
		}
	}

	return 0;
}

