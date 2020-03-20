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

using namespace std;

int main(int argc, char *argv[]) {
	//char fileName[50] = "source-graphs/LeNet-0765@vertices.grf";
	//char targetName[27] = "target-graphs/machines.tgt";

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
		//g.printGraph();
	}

	// Adjacency matrix
	/*for (int i = 0; i < g.getNumberOfVertices(); i++) {
		for (int j = 0; j < g.getNumberOfVertices(); j++)
			cout << g.getAdjMatrixEdgeWeight(i, j) << " ";
		cout << "\n";
	}*/

	TargetGraph t(*r.getTargetGraph());

	// ----------------------//
	// ----------------------//
	// ----- Comm -----//
	if ((strcmp(argv[5], "DN2PCIoTcomm")) == 0) {
		Comm Comm2(&g, r.getTargetNumberOfVertices(), singleMove, swaps, &t, g.getNumberOfVertices(), forcedInput, initPart, argv[6], argv[4], numberOfThreads);

		if ((strcmp(argv[6], "verb")) == 0) {
			t.printGraph();
		}
		cout << "\nsourceGraphFile: " << argv[1] << "\ntargetGraphFile: " << argv[2];
		cout << "\nSingleMoves: " << singleMove << ", swaps: " << swaps << ", locked input: " << forcedInput << ", initPart: " << initPart << ", DN2PCIoTcomm, number of threads: " << numberOfThreads;
		cout << "\ninitPartFileName: " << argv[4];

		cout << "\n\nInitial partitioning:";
		Comm2.printPartitioning(Comm2.getInitialPartitioning());
		cout << "initialCost: " << Comm2.computeCost(Comm2.getInitialPartitioning(), false) << "\n";
		int memoryValid = Comm2.validPartitioning(Comm2.getInitialPartitioning());
		for (int i = 0; i < Comm2.getNumberOfPartitions(); i++) {
			cout << "mem[" << i << "]: " << Comm2.getValidArray(i) << " ";
		}
		cout << "\nmemoryValid: " << memoryValid << "\n\n" << endl;
	
		Comm2.run(argv[3]);

		cout << "\nBest partitioning found:";
		Comm2.printPartitioning(Comm2.getBestPartitioning());
		cout << "bestCost: " << Comm2.computeCost(Comm2.getBestPartitioning(), false) << "\n";
		memoryValid = Comm2.validPartitioning(Comm2.getBestPartitioning());
		for (int i = 0; i < Comm2.getNumberOfPartitions(); i++) {
			cout << "mem[" << i << "]: " << Comm2.getValidArray(i) << " ";
		}
		cout << "\nmemoryValid: " << memoryValid << "\n\n";

	// ----------------------------//
	// ----------------------------//
	// ----- InferenceRate -----//
	} else if ((strcmp(argv[5], "DN2PCIoTiR")) == 0) {
		InferenceRate InferenceRate(&g, r.getTargetNumberOfVertices(), singleMove, swaps, &t, g.getNumberOfVertices(), forcedInput, initPart, argv[6], argv[4], numberOfThreads);

		if ((strcmp(argv[6], "verb")) == 0) {
			t.printGraph();
		}
		cout << "\nsourceGraphFile: " << argv[1] << "\ntargetGraphFile: " << argv[2];
		cout << "\nSingleMoves: " << singleMove << ", swaps: " << swaps << ", locked input: " << forcedInput << ", initPart: " << initPart << ", DN2PCIoTiR, number of threads: " << numberOfThreads;
		cout << "\ninitPartFileName: " << argv[4];

		cout << "\n\nInitial partitioning:";

		InferenceRate.printPartitioning(InferenceRate.getInitialPartitioning());
		cout << "initialCost: " << -InferenceRate.computeCost(InferenceRate.getInitialPartitioning(), false) << "\n";
		int memoryValid = InferenceRate.validPartitioning(InferenceRate.getInitialPartitioning());
		for (int i = 0; i < InferenceRate.getNumberOfPartitions(); i++) {
			cout << "mem[" << i << "]: " << InferenceRate.getValidArray(i) << " ";
		}
		cout << "\nmemoryValid: " << memoryValid << "\n\n" << endl;
	
		InferenceRate.run(argv[3]);

		cout << "\nBest partitioning found:";
		InferenceRate.printPartitioning(InferenceRate.getBestPartitioning());
		cout << "bestCost: " << -InferenceRate.computeCost(InferenceRate.getBestPartitioning(), false) << "\n";
		memoryValid = InferenceRate.validPartitioning(InferenceRate.getBestPartitioning());
		for (int i = 0; i < InferenceRate.getNumberOfPartitions(); i++) {
			cout << "mem[" << i << "]: " << InferenceRate.getValidArray(i) << " ";
		}
		cout << "\nmemoryValid: " << memoryValid << "\n\n";
	}

	return 0;
}
