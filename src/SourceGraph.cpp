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

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include "SourceGraph.h"

using namespace std;

/* ADJACENCY LIST REPRESENTATION: the function GRAPHinit() builds a graph with vertices 0 1 .. V-1 and no arc. */
SourceGraph::SourceGraph(int n, int nLayers, int nPartitions, bool typ) : type(typ) { 
	vertex v;
	int odd=0;

	//cout << "\n type: "<< type;

	V = n;
	/*while ((V % nPartitions) != 0) {
		V++;
		odd = 1;
	}*/
	//cout << "\nV: " << V;
	A = 0;
	enableVertexLabels = 0;
	enableEdgeWeights = 0;
	enableVertexWeights = 0;
	enableMemory = 0;
	setNumberOfLayers(nLayers);
	layerInitialPos = (int *) malloc(nLayers * sizeof(int));
	if (layerInitialPos == NULL) {
		cout << "layerInitialPos could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	sharedParam = (int *) malloc(nLayers * sizeof(int));
	if (sharedParam == NULL) {
		cout << "sharedParam could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	layerWidth = (int *) malloc(nLayers * sizeof(int));
	if (layerWidth == NULL) {
		cout << "layerWidth could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	layerHeight = (int *) malloc(nLayers * sizeof(int));
	if (layerHeight == NULL) {
		cout << "layerHeight could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	originalDepth = (int *) malloc(nLayers * sizeof(int));
	if (originalDepth == NULL) {
		cout << "originalDepth could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	vertexWeight = (int *) malloc(V * sizeof(int));
	if (vertexWeight == NULL) {
		cout << "vertexWeight could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	memory = (int *) malloc(V * sizeof(int));
	if (memory == NULL) {
		cout << "memory could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	match = (int *) malloc(V * sizeof(int));
	if (match == NULL) {
		cout << "match could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	map = (int *) malloc(V * sizeof(int));
	if (map == NULL) {
		cout << "map could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	degree = (int *) malloc(V * sizeof(int));
	if (degree == NULL) {
		cout << "degree could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	sharedParamPerVertex = (bool **) malloc(V * sizeof(bool *));
	if (sharedParamPerVertex == NULL) {
		cout << "sharedParamPerVertex could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < V; i++) {
		sharedParamPerVertex[i] = (bool *) malloc(nLayers * sizeof(bool));
		if (sharedParamPerVertex[i] == NULL) {
			cout << "sharedParamPerVertex[i] could not be allocated! \n";
			exit(EXIT_FAILURE);
		}
	}
	/*adjMatrixSize = (int *) malloc(V * sizeof(int));
	if (adjMatrixSize == NULL) {
		cout << "adjMatrixSize could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	adjMatrix = (node_matrix **) malloc (V * sizeof (node_matrix *));
	if (adjMatrix == NULL) {
		cout << "adjMatrix could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < V; i++) {
		adjMatrix[i] = (node_matrix *) calloc (100000, sizeof (node_matrix));
		if (adjMatrix[i] == NULL) {
			cout << "adjMatrix[i] could not be allocated! \n";
			exit(EXIT_FAILURE);
		}
	}*/
	adj = (links *) malloc(V * sizeof (links));
	if (adj == NULL) {
		cout << "adj could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (v = 0; v < V; ++v) {
		vertexWeight[v] = 0;
		adj[v] = NULL;
		memory[v] = 1;
		match[v] = v;
		map[v] = -1;
		degree[v] = 0;
		//adjMatrixSize[v] = 0;
		for (int i = 0; i < nLayers; i++)
			sharedParamPerVertex[v][i] = false;
	}
	/*if (odd == 1) {
		for (int i = V - 1; i >= n ; i--) {	
			memory[i] = 0;
		}
	}*/
}

/* The function NEWnode() receives a vertex w and the address next of a node and returns the address of a new node such that a->w == w and a->next == next. */
//static link SourceGraph::NEWnode(vertex w, link next, int edgeW) { 
links SourceGraph::NEWnode(vertex w, links next, int edgeW, /*int redN,*/ int src, int srcMlvl) { 
	links a = (links) malloc(sizeof (struct node));

	a->w = w; 
	a->edgeWeight = edgeW;
	a->source = src;
	a->sourceMlvl = srcMlvl;
	//a->redundantNeurons = redN;
	a->matched = 0;
	a->next = next;     
	return a;                         
}

/* ADJACENCY LIST REPRESENTATION: the function GRAPHinsertArc() inserts an arc v-w in graph G. The function supposes that v and w are distinct, positive, and smaller than G->V. If the graph already has an arc v-w, the function does not do anything. */
void SourceGraph::insertArc(vertex v, vertex w, int edgeW, /*int redN,*/ int src) { 
	links a;

	if (w == -1) {
	   return;
	}
	for (a = adj[v]; a != NULL; a = a->next) 
	  if (a->w == w) return;
	adj[v] = NEWnode(w, adj[v], edgeW, /*redN,*/ src);
	A++;
	degree[v]++;
	//degree[w]++;
}

void SourceGraph::insertArcSum(vertex v, vertex w, int edgeW, int src) { 
	links a;

	if (w == -1) {
	   return;
	}
	for (a = adj[v]; a != NULL; a = a->next) 
		if (a->w == w) {
			a->edgeWeight += edgeW;
			return;
		}
	adj[v] = NEWnode(w, adj[v], edgeW, -2222, src);
	A++;
}

bool SourceGraph::insertArcSrcMlvl(vertex v, vertex w, int edgeW, /*int redN,*/ int src) { 
	links a;

	if (w == -1) {
		return false;
	}


	for (a = adj[v]; a != NULL; a = a->next) {
		if (a->sourceMlvl == src && a->w == w && a->edgeWeight == edgeW) {
			return false;
		}
	}	

	adj[v] = NEWnode(w, adj[v], edgeW, /*redN,*/ -1, src);
	A++;
	degree[v]++;
	//degree[w]++;
	return true;
}

bool SourceGraph::insertArcSrcMlvlPartitionGraph(vertex v, vertex w, int edgeW, /*int redN,*/ int src, int srcMlvl) { 
	links a;

	if (w == -1) {
		return false;
	}


	for (a = adj[v]; a != NULL; a = a->next) {
		if (a->source == src && a->w == w && a->edgeWeight == edgeW && a->sourceMlvl == srcMlvl) {
			return false;
		}
	}	

	adj[v] = NEWnode(w, adj[v], edgeW, /*redN,*/ src, srcMlvl);
	A++;
	degree[v]++;
	//degree[w]++;
	return true;
}

bool SourceGraph::insertArcSrc(vertex v, vertex w, int edgeW, /*int redN,*/ int src) { 
	links a;

	if (w == -1) {
		return false;
	}


	for (a = adj[v]; a != NULL; a = a->next) {
		if (a->source == src && a->w == w) {
			return false;
		}
	}	

	adj[v] = NEWnode(w, adj[v], edgeW, /*redN,*/ src);
	A++;
	return true;
}

bool SourceGraph::removeArcSrc(vertex v, vertex w, int src) { 
	links a = adj[v], prev=NULL;

	if (w == -1) {
		return false;
	}
	while (a != NULL) {
		//cout << "\n w: " << a->w << " src: " << a->source;
		if (a->w == w && a->source == src) {
			break;
		}
		prev = a;
		a = a->next;
	}
	if (prev == NULL) {
		adj[v] = adj[v]->next;
		free(a);
		A--;
		return true; 
	} else if (a != NULL) {
		prev->next = a->next;
		free(a);
		A--;
		return true;
	}
	return false;
}

/*void SourceGraph::setAdjMatrix(vertex v, vertex w, double edgeW) { 
	if (v == w || edgeW <= 0) {
		cout << "\n There may not be auto edges \n or there may not be zero or negative edge weights! \n";
		exit(EXIT_FAILURE);
	}
	adjMatrix[v][w].edgeWeight = edgeW;
	//adjMatrix[w][v] = edgeW;
}*/

/* Updates vertex weights */
void SourceGraph::setVertexWeight(vertex v, int vertexW) {
	vertexWeight[v] = vertexW;
}

/* Updates memory weights */
void SourceGraph::setMemoryWeight(vertex v, int memoryWeight) {
	memory[v] = memoryWeight;
}

/*int SourceGraph::getNumberOfVertices() const {
	return V;
}*/

/*link SourceGraph::getAdjOfVertex(int i) const {
	return adj[i];
}*/

int SourceGraph::getLayerOfVertex(int v) const {
	int layer = 0;
	
	for (int j = 0; j < getNumberOfLayers(); j++) {
		if (layer == getNumberOfLayers() - 1)
			break;
		if (v < getLayerInitialPos(j + 1)) {
			break;
		} else {
			layer++;
		}
	}
	return layer;
}

void SourceGraph::setFlags(int vertexLabel, int edgeWeight, int vertexWeight, int memory) {
	enableVertexLabels = vertexLabel;
	enableEdgeWeights = edgeWeight;
	enableVertexWeights = vertexWeight;
	enableMemory = memory;
}

void SourceGraph::printGraphHeader() const {
	cout << "\nSource: V=" << V << ", A=" << A << ", enableMemory=" << enableMemory << ", enableVertexLabels=" << enableVertexLabels << ", enableEdgeWeights=" << enableEdgeWeights << ", enableVertexWeights=" << enableVertexWeights << "\n";
}

void SourceGraph::printGraphSrc() const {
	int i=0;
	links aux;

	cout << "\nCost graph: V=" << V << ", A=" << A << "\n";

	for (i = 0; i < V; i++) {
		cout << "Vertex[" << i << "]: ";
		aux = adj[i];
		while(aux != NULL){
			cout << aux->w << " (" << aux->source << ") ";
			aux = aux->next;
		}
		if (enableEdgeWeights == 1) {
			cout << "Edge weights: ";
			aux = adj[i];
			while(aux != NULL) {
				cout << aux->edgeWeight << " ";
				aux = aux->next;
			}
		}
		cout << "\n";
	}
}


void SourceGraph::printGraph() const {
	int i=0;
	links aux;

	cout << "\nSource: V=" << V << ", A=" << A << ", enableMemory=" << enableMemory << ", enableVertexLabels=" << enableVertexLabels << ", enableEdgeWeights=" << enableEdgeWeights << ", enableVertexWeights=" << enableVertexWeights << "\n";

	cout << "Number of layers: " << numberOfLayers << "\nLayers start at: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << layerInitialPos[i] << " ";
	}
	cout << "\nNumber of shared parameters per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << sharedParam[i] << " ";
	}
	cout << "\nLayer width per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << layerWidth[i] << " ";
	}
	cout << "\nLayer height per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << layerHeight[i] << " ";
	}
	cout << "\nOriginal depth per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << originalDepth[i] << " ";
	}
	cout << "\n";

	for (i = 0; i < V; i++) {
		cout << "Vertex[" << i << "]: ";
		aux = adj[i];
		/*while(aux != NULL){
			cout << aux->w << " ";
			aux = aux->next;
		}
		if (enableEdgeWeights == 1) {
			cout << "Edge weights: ";
			aux = adj[i];
			while(aux != NULL) {
				cout << aux->edgeWeight << " ";
				aux = aux->next;
			}
		}*/
		cout << "Memory: " << memory[i] << " ";
		if (enableVertexWeights == 1) {
			printf("Computational weight: %d ", vertexWeight[i]);
		}
		cout << "Degree: " << degree[i] << " ";

		printf("Shared parameters: ");
		for (int l = 0; l < numberOfLayers; l++) {
			if (sharedParamPerVertex[i][l] == true) {
				printf("1 ");
			} else {
				printf("0 ");
			}
		}

		cout << "\n";
	}
}

void SourceGraph::printGraphMlvl() const {
	int i=0;
	links aux;

	cout << "\nSource: V=" << V << ", A=" << A << ", enableMemory=" << enableMemory << ", enableVertexLabels=" << enableVertexLabels << ", enableEdgeWeights=" << enableEdgeWeights << ", enableVertexWeights=" << enableVertexWeights << "\n";

	cout << "Number of layers: " << numberOfLayers << "\nLayers start at: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << layerInitialPos[i] << " ";
	}
	cout << "\nNumber of shared parameters per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << sharedParam[i] << " ";
	}
	cout << "\nLayer width per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << layerWidth[i] << " ";
	}
	cout << "\nLayer height per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << layerHeight[i] << " ";
	}
	cout << "\nOriginal depth per layer: ";
	for (i = 0; i < numberOfLayers; i++) {
		cout << originalDepth[i] << " ";
	}
	cout << "\n";

	for (i = 0; i < V; i++) {
		cout << "Vertex[" << i << "]: ";
		aux = adj[i];
		while(aux != NULL){
			cout << aux->w << " ";
			aux = aux->next;
		}
		if (enableEdgeWeights == 1) {
			cout << "Edge weights: ";
			aux = adj[i];
			while(aux != NULL) {
				cout << aux->edgeWeight << " ";
				aux = aux->next;
			}
		}
		cout << "Source: ";
		aux = adj[i];
		while(aux != NULL) {
			cout << aux->sourceMlvl << " ";
			aux = aux->next;
		}
		cout << "Memory: " << memory[i] << " ";
		if (enableVertexWeights == 1) {
			printf("Computational weight: %d ", vertexWeight[i]);
		}
		cout << "Degree: " << degree[i] << " ";

		printf("Shared parameters: ");
		for (int l = 0; l < numberOfLayers; l++) {
			if (sharedParamPerVertex[i][l] == true) {
				printf("1 ");
			} else {
				printf("0 ");
			}
		}

		cout << "\n";
	}
}

/*int SourceGraph::getEnableMemory() const {
	return enableMemory;
}*/

/*int SourceGraph::getMemory(vertex v) const {
	return memory[v];
}*/

void SourceGraph::setNumberOfLayers(int num) {
	numberOfLayers = num;
}

/*int SourceGraph::getNumberOfLayers() const {
	return numberOfLayers;
}*/

void SourceGraph::setLayerInitialPos(int i, int pos) {
	layerInitialPos[i] = pos;
}

/*int SourceGraph::getLayerInitialPos(int i) const {
	return layerInitialPos[i];
}*/

void SourceGraph::setSharedParam(int i, int numP) {
	sharedParam[i] = numP;
}

int SourceGraph::getSharedParams(int u) const {
	int shrdp = 0;

	for (int i = 0; i < numberOfLayers; i++) {
		shrdp += sharedParamPerVertex[u][i];
	}
	return shrdp;
}

void SourceGraph::setLayerWidth(int i, int pos) {
	layerWidth[i] = pos;
}

void SourceGraph::setLayerHeight(int i, int pos) {
	layerHeight[i] = pos;
}

void SourceGraph::setOriginalDepth(int i, int pos) {
	originalDepth[i] = pos;
}

/*int SourceGraph::getAdjMatrixEdgeWeight(int i, int j) {
	return adjMatrix[i][j].edgeWeight;
}*/

void SourceGraph::setMatch(int u, int v) {
	match[u] = v;
	match[v] = u;
}

void SourceGraph::setEdgeMatch(links node) {
	node->matched = 1;
}

void SourceGraph::setMap(int u, int v, int label) {
	map[u] = label;
	map[v] = label;
}

void SourceGraph::setVertexDegree(int v, int deg) {
	degree[v] = deg;
}

void SourceGraph::setNodeSharedParam(int v, int layer) {
	sharedParamPerVertex[v][layer] = true;
}

void SourceGraph::setNodeSharedParamArray(int v, bool *shrdppv) {
	for (int i = 0; i < numberOfLayers; i++) {
		if (sharedParamPerVertex[v][i] == false) 
			sharedParamPerVertex[v][i] = shrdppv[i];
	}
}

void SourceGraph::getSharedParamArray(int v, bool *shrdppv) const {
	for (int i = 0; i < numberOfLayers; i++) {
		shrdppv[i] = sharedParamPerVertex[v][i];
	}
}

void SourceGraph::getSharedParamPerLayer(int *shrdppl) const {
	for (int i = 0; i < numberOfLayers; i++) {
		shrdppl[i] = sharedParam[i];
	}
}

void SourceGraph::setSharedParamPerLayer(int *shrdppl) const {
	for (int i = 0; i < numberOfLayers; i++) {
		sharedParam[i] = shrdppl[i];
	}
}

// copy constructor
SourceGraph::SourceGraph(const SourceGraph &graphToCopy) 
	: V(graphToCopy.V), A(graphToCopy.A), 
	enableVertexLabels(graphToCopy.enableVertexLabels), 
	enableEdgeWeights(graphToCopy.enableEdgeWeights), 
	enableVertexWeights(graphToCopy.enableVertexWeights), 
	enableMemory(graphToCopy.enableMemory), 
	numberOfLayers(graphToCopy.numberOfLayers),
	type(graphToCopy.type) 
{
	//cout << "\nCopy constructor called\n";
	layerInitialPos = (int *) malloc(numberOfLayers * sizeof(int));
	if (layerInitialPos == NULL) {
		cout << "layerInitialPos could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	layerWidth = (int *) malloc(numberOfLayers * sizeof(int));
	if (layerWidth == NULL) {
		cout << "layerWidth could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	layerHeight = (int *) malloc(numberOfLayers * sizeof(int));
	if (layerHeight == NULL) {
		cout << "layerHeight could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	originalDepth = (int *) malloc(numberOfLayers * sizeof(int));
	if (originalDepth == NULL) {
		cout << "originalDepth could not be allocated! \n";
		exit(EXIT_FAILURE);
	}	
	sharedParam = (int *) malloc(numberOfLayers * sizeof(int));
	if (sharedParam == NULL) {
		cout << "sharedParam could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < numberOfLayers; i++) {
		layerInitialPos[i] = graphToCopy.layerInitialPos[i];
		sharedParam[i] = graphToCopy.sharedParam[i];
		layerWidth[i] = graphToCopy.layerWidth[i];
		layerHeight[i] = graphToCopy.layerHeight[i];
		originalDepth[i] = graphToCopy.originalDepth[i];
	}

	memory = (int *) malloc(V * sizeof(int));
	if (memory == NULL) {
		cout << "memory could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	degree = (int *) malloc(V * sizeof(int));
	if (degree == NULL) {
		cout << "degree could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	match = (int *) malloc(V * sizeof(int));
	if (match == NULL) {
		cout << "match could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	map = (int *) malloc(V * sizeof(int));
	if (map == NULL) {
		cout << "map could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	vertexWeight = (int *) malloc(V * sizeof(int));
	if (vertexWeight == NULL) {
		cout << "vertexWeight could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	sharedParamPerVertex = (bool **) malloc(V * sizeof(bool *));
	if (sharedParamPerVertex == NULL) {
		cout << "sharedParamPerVertex could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < V; i++) {
		sharedParamPerVertex[i] = (bool *) malloc(numberOfLayers * sizeof(bool));
		if (sharedParamPerVertex[i] == NULL) {
			cout << "sharedParamPerVertex[i] could not be allocated! \n";
			exit(EXIT_FAILURE);
		}
	}
	/*adjMatrix = (node_matrix **) malloc (V * sizeof (node_matrix *));
	if (adjMatrix == NULL) {
		cout << "adjMatrix could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < V; i++) {
		adjMatrix[i] = (node_matrix *) calloc (V, sizeof (node_matrix));
		if (adjMatrix[i] == NULL) {
			cout << "adjMatrix[i] could not be allocated! \n";
			exit(EXIT_FAILURE);
		}
	}*/
	adj = (links *) malloc(V * sizeof (links));
	if (adj == NULL) {
		cout << "adj could not be allocated! \n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < V; i++) {
		/*for (int j = 0; j < V; j++) {
			adjMatrix[i][j] = graphToCopy.adjMatrix[i][j];
		}*/
		adj[i] = NULL;
		memory[i] = graphToCopy.memory[i];
		match[i] = graphToCopy.match[i];
		map[i] = graphToCopy.map[i];
		degree[i] = graphToCopy.degree[i];
		vertexWeight[i] = graphToCopy.vertexWeight[i];
		for (int j = 0; j < numberOfLayers; j++) {
			sharedParamPerVertex[i][j] = graphToCopy.sharedParamPerVertex[i][j];
		}
		links no = graphToCopy.adj[i];
		while (no != NULL) {
			if (type == false) {
				insertArc(i, no->w, no->edgeWeight, /*no->redundantNeurons,*/ no->source);
			} else {
				insertArcSrc(i, no->w, no->edgeWeight, /*no->redundantNeurons,*/ no->source);
			}
			no = no->next;
		}
	}
}

// The parameter is a pointer to a pointer to the struct graph 
SourceGraph::~SourceGraph() {
	links no=NULL, aux=NULL;
	int i;

	if (layerInitialPos != NULL) {
		free(layerInitialPos);
		layerInitialPos = NULL;
	}
	if (layerWidth != NULL) {
		free(layerWidth);
		layerWidth = NULL;
	}
	if (layerHeight != NULL) {
		free(layerHeight);
		layerHeight = NULL;
	}
	if (originalDepth != NULL) {
		free(originalDepth);
		originalDepth = NULL;
	}
	if (sharedParam != NULL) {
		free(sharedParam);
		sharedParam = NULL;
	}
	if (memory != NULL) {
		free(memory);	
		memory = NULL;
	}
	if (match != NULL) {
		free(match);	
		match = NULL;
	}
	if (map != NULL) {
		free(map);	
		map = NULL;
	}
	if (degree != NULL) {
		free(degree);	
		degree = NULL;
	}
	if (vertexWeight != NULL) {
		free(vertexWeight);	
		vertexWeight = NULL;
	}
	if (sharedParamPerVertex != NULL) {
		for (i = 0; i < V; i++) {
			if (sharedParamPerVertex[i] != NULL) {
				free(sharedParamPerVertex[i]);
			}
		}
		free(sharedParamPerVertex);
		sharedParamPerVertex = NULL;
	}
	/*if (adjMatrix != NULL) {
		for (i = 0; i < V; i++) {
			if (adjMatrix[i] != NULL) {
				free(adjMatrix[i]);
			}
		}
		free(adjMatrix);
		adjMatrix = NULL;
	}*/
	if (adj != NULL) {
		for (i = 0; i < V; i++) {
			no = adj[i];
			while (no != NULL) {
				aux = no;
				no = no->next;
				free(aux);
				aux = NULL;
			}
		}
		no = NULL;
		free(adj);
		adj = NULL;
	}
}
