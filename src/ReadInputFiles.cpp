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
ReadInputFiles::ReadInputFiles(char *fileName, char *targetName, int numberOfThreads) {
	FILE *fp=NULL;
	char flagAux, adjAux[600000], *fgetsAux, *fgetsEdge, *caux, comment[]="#";
	int i, j, weightAux, enableVertexLabels, enableEdgeWeights, enableVertexWeights, enableMemoryWeights, aux, numVertices, numTarget;

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
	machines = new TargetGraph(numTarget, numberOfThreads);

	for(i = 0; i < numTarget; i++) {
		// read computational weights
		aux = fscanf(fp, "%d %d ", &weightAux, &weightAux);
		machines->setComputationalWeight(i, weightAux);

		// read memory weights
		aux = fscanf(fp, "%d %d ", &weightAux, &weightAux);
		machines->setMemoryWeight(i, weightAux);
	}

	// TODO: jump to next line if it is a comment (#)
	/*aux = fscanf(fp, "%c", &flagAux);
	strcpy(adjAux, &flagAux);
	printf("\n %s", &flagAux);
	printf("\n %s", comment);
	while (strcmp(&flagAux,comment) == 0) {
		printf("\n #");
		fgetsAux = fgets(adjAux, 6000, fp);
		aux = fscanf(fp, "%c", &flagAux);
		strcpy(adjAux, &flagAux);
	}*/

	// read connections between devices
	for(i = 1; i < numTarget; i++) {
		for (j = 0; j < i; j++) {
			aux = fscanf(fp, "%d ", &weightAux);
			machines->setConnectionWeight(i, j, weightAux);
		}
	}

	fclose(fp);
	fp = NULL;

	/* __________________________________________________________________________________ */
	/* Opens SourceGraph file to read SourceGraph data*/
	fp = fopen(fileName, "r");
	if(fp == NULL){
		cout << "SourceGraph file could not be open! \n";
		exit(EXIT_FAILURE);
	}

	// read numVertices
	i = fscanf(fp, "%d\n%d ", &numVertices, &numVertices);

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

	// read number of layers
	int numberOfLayers;
	i = fscanf(fp, "%d", &numberOfLayers);

	// instantiate the source graph
	inputGraph = new SourceGraph(numVertices, numberOfLayers, numTarget);
	inputGraph->setFlags(enableVertexLabels, enableEdgeWeights, enableVertexWeights, enableMemoryWeights);
	
	// read positions at which each layer begins
	int layerInitialPos, sharedParam;
	for (i = 0; i < numberOfLayers; i++) {
		aux = fscanf(fp, "%d ", &layerInitialPos);
		inputGraph->setLayerInitialPos(i, layerInitialPos);
	}
	// read number of shared parameters per layer
	for (i = 0; i < numberOfLayers; i++) {
		aux = fscanf(fp, "%d ", &sharedParam);
		inputGraph->setSharedParam(i, sharedParam);
	}

	// read width and height of each layer
	int layerWidth, layerHeight;
	for (i = 0; i < numberOfLayers; i++) {
		aux = fscanf(fp, "%d ", &layerWidth);
		//cout << "\n lw: " << layerWidth;
		inputGraph->setLayerWidth(i, layerWidth);
	}
	for (i = 0; i < numberOfLayers; i++) {
		aux = fscanf(fp, "%d ", &layerHeight);
		inputGraph->setLayerHeight(i, layerHeight);
	}

	// read original depth of each layer
	int orDepth;
	for (i = 0; i < numberOfLayers; i++) {
		aux = fscanf(fp, "%d ", &orDepth);
		inputGraph->setOriginalDepth(i, orDepth);
	}	

	// read each vertex properties (memory, vertex (computational) weight, degree, edges, edge weights)
	for(i = 0; i < numVertices; i++) {
		// new way of accounting the memory of shared parameters (for multilevel KLP)
		inputGraph->setNodeSharedParam(i, inputGraph->getLayerOfVertex(i));

		// read memory weights if necessary 
		if(enableMemoryWeights != 0) {
			aux = fscanf(fp, "%d ", &weightAux);
			inputGraph->setMemoryWeight(i, weightAux);
		} 

		// read vertex weights if necessary
		if(enableVertexWeights != 0) {
			aux = fscanf(fp, "%d ", &weightAux);
			inputGraph->setVertexWeight(i, weightAux);
		}

		// read vertex degree
		aux = fscanf(fp, "%d ", &weightAux);
		inputGraph->setVertexDegree(i, weightAux);
		//cout << "\n i: " << i << ", deg: " << weightAux;
		if (weightAux == 0 || (i >= inputGraph->getLayerInitialPos(numberOfLayers - 1))) continue;

		if(enableEdgeWeights == 0) {
			// read vertex adjacency list
			caux = fgets(adjAux, sizeof(adjAux), fp);
			fgetsAux = strtok(adjAux, " ");
			while(fgetsAux != NULL){
				//cout << "\ni: " << i << " edge: " << atoi(fgetsAux);
				inputGraph->insertArc(i, atoi(fgetsAux), -2, 1);
				//inputGraph->setAdjMatrix(i, atoi(fgetsAux), 1);
				fgetsAux = strtok(NULL, " ");
			}
		} else {
			// read edge weight before edge 
			// (edge weight[0] adjacency[0] edge weight[1] adjacency[1] ...)
			caux = fgets(adjAux, sizeof(adjAux), fp);
			fgetsEdge = strtok(adjAux, " ");
			while(fgetsEdge != NULL){
				fgetsAux = strtok(NULL, " ");

				// the algorithms do not use this tag
				// define tag to calculate redundant neurons in each edge send
				int redN = -2222, iLayer = 0, adjLayer = 0;
				for (int j = 0; j < inputGraph->getNumberOfLayers(); j++) {
					if (iLayer == inputGraph->getNumberOfLayers() - 1)
						break;
					if (i < inputGraph->getLayerInitialPos(j + 1)) {
						break;
					} else {
						iLayer++;
					}
				}
				for (int j = 0; j < inputGraph->getNumberOfLayers(); j++) {
					if (adjLayer == inputGraph->getNumberOfLayers() - 1)
						break;
					//cout << "\n fgetsEdge: " << fgetsEdge << ", fgetsAux: " << fgetsAux << ", iLayer: " << iLayer << ", layerInitialPos: " << inputGraph->getLayerInitialPos(j + 1);
					if (atoi(fgetsAux) < inputGraph->getLayerInitialPos(j + 1)) {
						break;
					} else {
						adjLayer++;
					}
				}
				int k, j; // line, column of i
				k = (i - inputGraph->getLayerInitialPos(iLayer)) / inputGraph->getLayerWidth(iLayer);
				j = (i - inputGraph->getLayerInitialPos(iLayer)) % inputGraph->getLayerWidth(iLayer);
				int l, m; // line, column of atoi(fgetsAux) (adj)
				l = (atoi(fgetsAux) - inputGraph->getLayerInitialPos(adjLayer)) / inputGraph->getLayerWidth(iLayer);
				m = (atoi(fgetsAux) - inputGraph->getLayerInitialPos(adjLayer)) % inputGraph->getLayerWidth(iLayer);

				// not used
				if (adjLayer >= 5) {
					redN = -2222;
				} else if (k == l && j == m) {
					redN = 0;
				} else if (k + 1 == l && j + 1 == m) {
					redN = 33;
				} else if (k + 1 == l && j == m) {
					redN = 32;
				} else if (k == l && j + 1 == m) {
					redN = 1;
				} else if (k + 1 == l && j - 1 == m) {
					redN = 31;
				} else if (k == l && j - 1 == m) {
					redN = -1;
				} else if (k - 1 == l && j + 1 == m) {
					redN = -31;
				} else if (k - 1 == l && j == m) {
					redN = -32;
				} else if (k - 1 == l && j - 1 == m) {
					redN = -33;
				}

				//cout << "\n i: " << i << " iLayer: " << iLayer << " adj: " << atoi(fgetsAux) << " adjLayer: " << adjLayer << " redN: " << redN;

				//printf("\n%d %d %d", atoi(fgetsEdge), atoi(fgetsAux), i);
				inputGraph->insertArc(i, atoi(fgetsAux), atoi(fgetsEdge)/*, redN*/);
				//inputGraph->setAdjMatrix(i, atoi(fgetsAux), atof(fgetsEdge));
				fgetsEdge = strtok(NULL, " ");
			}
		}
	} 

	fclose(fp);
	fp = NULL;
	//return Source_graph;
}

const SourceGraph* ReadInputFiles::getSourceGraph() const {
	return inputGraph;
}

int ReadInputFiles::getNumberOfVertices() const {
	return inputGraph->getNumberOfVertices();
}

int ReadInputFiles::getTargetNumberOfVertices() const {
	return machines->getNumberOfVertices();
}

const TargetGraph* ReadInputFiles::getTargetGraph() const {
	return machines;
}

ReadInputFiles::~ReadInputFiles() {
	//cout << "\nReadInputFiles object destructor\n";
	delete inputGraph;
	delete machines;
}

