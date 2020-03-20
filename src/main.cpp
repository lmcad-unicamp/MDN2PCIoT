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
#include "ReadInputFiles.h"
#include "Comm.h"

using namespace std;

int main() {
	char fileName[29] = "source-graphs/LeNet-765v.grf";
	char targetName[27] = "target-graphs/machines.tgt";

	bool singleMove = true;
	bool swaps = true;
	int forcedInput = 0; // 256 to force input
	bool initPart = false;  //false: random initial partitioning
							// true: read initial partitioning from initPart.txt

	ReadInputFiles r(fileName, targetName);

	SourceGraph g(r.getNumberOfVertices(), r.getTargetNumberOfVertices());
	g = *r.getSourceGraph();
	g.printGraph();

	TargetGraph t(r.getTargetNumberOfVertices());
	t = *r.getTargetGraph();

	Comm Comm2(&g, r.getTargetNumberOfVertices(), singleMove, swaps, &t, g.getNumberOfVertices(), forcedInput, initPart);

	t.printGraph();
	cout << "\nSingleMoves: " << singleMove << ", swaps: " << swaps << ", locked input: " << forcedInput;

	cout << "\n\nInitial partitioning:";
	Comm2.printPartitioning(Comm2.getInitialPartitioning());
	cout << "initialCost: " << Comm2.computeCost(Comm2.getInitialPartitioning()) << "\n";
	int memoryValid = Comm2.validPartitioning(Comm2.getInitialPartitioning());
	for (int i = 0; i < Comm2.getNumberOfPartitions(); i++) {
		cout << "mem[" << i << "]: " << Comm2.getValidArray(i) << " ";
	}
	cout << "\nmemoryValid: " << memoryValid << "\n\n";
	
	Comm2.run();

	cout << "\nBest partitioning found:";
	Comm2.printPartitioning(Comm2.getBestPartitioning());
	cout << "bestCost: " << Comm2.computeCost(Comm2.getBestPartitioning()) << "\n";
	memoryValid = Comm2.validPartitioning(Comm2.getBestPartitioning());
	for (int i = 0; i < Comm2.getNumberOfPartitions(); i++) {
		cout << "mem[" << i << "]: " << Comm2.getValidArray(i) << " ";
	}
	cout << "\nmemoryValid: " << memoryValid << "\n\n";
	
	return 0;
}
