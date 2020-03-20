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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "ReadInputFiles.h"
#include "SourceGraph.h"
#include "TargetGraph.h"

using namespace std;

/*_________________________________________________________*/
// Graph definition (SourceGraph) 
// graph represented as an adjacency list
/*_________________________________________________________*/

// read grf file defined as the same input to SCOTCH (except for memory weights, which are not present in SCOTCH)
ReadInputFiles::ReadInputFiles(char *fileName, char *targetName) {
	FILE *fp=NULL;
	char flagAux, adjAux[6000], *fgetsAux, *fgetsEdge, *caux;
	int i, weightAux, enableVertexLabels, enableEdgeWeights, enableVertexWeights, enableMemoryWeights, aux, numVertices, numTarget;

	// read target file
	char isDeco[8];

	/* Opens SourceGraph file to read SourceGraph data*/
	fp = fopen(targetName, "r");
	if(fp == NULL){
		cout << "TargetGraph file could not be open! \n";
		exit(EXIT_FAILURE);
	}

	// read if the target is decomposition-defined
	caux = fgets(isDeco, 7, fp);
	if(strncmp(isDeco,"deco 0", 6)) {
		cout << "TargetGraph file is not decomposition-defined! \n";
		exit(EXIT_FAILURE);	
	}

	// read num_vertices
	i = fscanf(fp, "%d %d ", &numTarget, &weightAux);
	machines = new TargetGraph(numTarget);

	for(i = 0; i < numTarget; i++) {
		// read computational weights
		aux = fscanf(fp, "%d %d ", &weightAux, &weightAux);
	    //setComputationalWeight(i, weightAux);

		// read memory weights
		aux = fscanf(fp, "%d %d ", &weightAux, &weightAux);
		machines->setMemoryWeight(i, weightAux);
	}

	/* TODO: read connections between devices
	for(
		// read vertex adjacency list
		fgets(adj_aux, sizeof(adj_aux), fp);
		fgets_aux = strtok(adj_aux, " ");
		while(fgets_aux != NULL){
			GraphInsertArc(Source_graph, i, atoi(fgets_aux)); 
			fgets_aux = strtok(NULL, " ");
		}			
	} 
	arrumar */

	fclose(fp);
	fp = NULL;

	/* Opens SourceGraph file to read SourceGraph data*/
	fp = fopen(fileName, "r");
	if(fp == NULL){
		cout << "SourceGraph file could not be open! \n";
		exit(EXIT_FAILURE);
	}

	// read numVertices
	i = fscanf(fp, "%d\n%d ", &numVertices, &numVertices);
	inputGraph = new SourceGraph(numVertices, numTarget);

	// (not used, other graph parameters)
	i = fscanf(fp, "%d\n%d ", &i, &i);

	// read flags
	flagAux = fgetc(fp);
	enableMemoryWeights = atoi(&flagAux);
	flagAux = fgetc(fp);
	enableVertexLabels = atoi(&flagAux);
	flagAux = fgetc(fp);
	enableEdgeWeights = atoi(&flagAux);
	flagAux = fgetc(fp);
	enableVertexWeights = atoi(&flagAux);
	inputGraph->setFlags(enableVertexLabels, enableEdgeWeights, enableVertexWeights, enableMemoryWeights);

	for(i = 0; i < numVertices; i++) {
		// read memory weights if necessary
		if(enableMemoryWeights != 0) {
			aux = fscanf(fp, "%d ", &weightAux);
			inputGraph->setMemoryWeight(i, weightAux);
		} 

		// read vertex weights if necessary
		if(enableVertexWeights != 0) {
			aux = fscanf(fp, "%d ", &weightAux);
		}

		// read vertex degree (not used)
		aux = fscanf(fp, "%d ", &weightAux);
		if (weightAux == 0) continue;

		if(enableEdgeWeights == 0) {
			// read vertex adjacency list
			caux = fgets(adjAux, sizeof(adjAux), fp);
			fgetsAux = strtok(adjAux, " ");
			while(fgetsAux != NULL){
				//cout << "\ni: " << i << " edge: " << atoi(fgetsAux);
				inputGraph->insertArc(i, atoi(fgetsAux), 1);
				inputGraph->setAdjMatrix(i, atoi(fgetsAux), 1);
				fgetsAux = strtok(NULL, " ");
			}
		} else {
			// read edge weight before edge 
			// (edge weight[0] adjacency[0] edge weight[1] adjacency[1] ...)
			caux = fgets(adjAux, sizeof(adjAux), fp);
			fgetsEdge = strtok(adjAux, " ");
			while(fgetsEdge != NULL){
				fgetsAux = strtok(NULL, " ");
				//printf("\n%d %d %d", atoi(fgetsEdge), atoi(fgetsAux), i);
				inputGraph->insertArc(i, atoi(fgetsAux), atoi(fgetsEdge));
				inputGraph->setAdjMatrix(i, atoi(fgetsAux), atof(fgetsEdge));
				fgetsEdge = strtok(NULL, " ");
			}
		}
	} 

	fclose(fp);
	fp = NULL;
	//return Source_graph;
}

ReadInputFiles::~ReadInputFiles() {

}

SourceGraph* ReadInputFiles::getSourceGraph() {
	return inputGraph;
}

int ReadInputFiles::getNumberOfVertices() {
	return inputGraph->getNumberOfVertices();
}

int ReadInputFiles::getTargetNumberOfVertices() {
	return machines->getNumberOfVertices();
}

TargetGraph* ReadInputFiles::getTargetGraph() {
	return machines;
}

