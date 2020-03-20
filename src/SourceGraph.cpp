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
SourceGraph::SourceGraph(int n, int nPartitions) { 
	vertex v;
	int odd=0;

	V = n;
	while ((V % nPartitions) != 0) {
		V++;
		odd = 1;
	}
	//cout << "\nV: " << V;
	A = 0;
	enableVertexLabels = 0;
	enableEdgeWeights = 0;
	enableVertexWeights = 0;
	/*vertexWeights = (int *) malloc(V * sizeof(int));
	if (vertexWeights == NULL) {
		cout << "vertexWeights could not be allocated! \n";
		exit(EXIT_FAILURE);
	}*/
	memory = (int *) malloc(V * sizeof(int));
	if (memory == NULL) {
		cout << "memory could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	adjMatrix = (double **) malloc (V * sizeof (double *));
	if (adjMatrix == NULL) {
		cout << "adjMatrix could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (int i = 0; i < V; i++) {
		adjMatrix[i] = (double *) calloc (V, sizeof (double));
		if (adjMatrix[i] == NULL) {
			cout << "adjMatrix[i] could not be allocated! \n";
			exit(EXIT_FAILURE);
		}
	}
	adj = (link *) malloc(V * sizeof (link));
	if (adj == NULL) {
		cout << "adj could not be allocated! \n";
		exit(EXIT_FAILURE);
	}
	for (v = 0; v < V; ++v) {
		//vertexWeights[v] = 1;
		adj[v] = NULL;
		memory[v] = 1;
	}
	if (odd == 1) {
		for (int i = V - 1; i >= n ; i--) {	
			memory[i] = 0;
		}
	}
}

/* The function NEWnode() receives a vertex w and the address next of a node and returns the address of a new node such that a->w == w and a->next == next. */
//static link SourceGraph::NEWnode(vertex w, link next, int edgeW) { 
link SourceGraph::NEWnode(vertex w, link next, int edgeW, int src) { 
	link a = (link) malloc(sizeof (struct node));

	a->w = w; 
	a->edgeWeight = edgeW;
	a->source = src;
	a->next = next;     
	return a;                         
}

/* ADJACENCY LIST REPRESENTATION: the function GRAPHinsertArc() inserts an arc v-w in graph G. The function supposes that v and w are distinct, positive, and smaller than G->V. If the graph already has an arc v-w, the function does not do anything. */
void SourceGraph::insertArc(vertex v, vertex w, int edgeW, int src) { 
	link a;

	if (w == -1) {
	   return;
	}
	for (a = adj[v]; a != NULL; a = a->next) 
	  if (a->w == w) return;
	adj[v] = NEWnode(w, adj[v], edgeW, src);
	A++;
}

bool SourceGraph::insertArcSrc(vertex v, vertex w, int edgeW, int src) { 
	link a;

	if (w == -1) {
		return false;
	}
	for (a = adj[v]; a != NULL; a = a->next) 
		if (a->source == src && a->w == w) {
			//cout << "\n v: " << v << " edge: " << a->w << " a->source: " << a->source; 
			return false; 
		}
	adj[v] = NEWnode(w, adj[v], edgeW, src);
	A++;
	return true;
}

bool SourceGraph::removeArcSrc(vertex v, vertex w, int src) { 
	link a = adj[v], prev=NULL;

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

void SourceGraph::setAdjMatrix(vertex v, vertex w, double edgeW) { 
	if (v == w || edgeW <= 0) {
		cout << "\n There may not be auto edges \n or there may not be zero or negative edge weights! \n";
		exit(EXIT_FAILURE);
	}
	adjMatrix[v][w] = edgeW;
	adjMatrix[w][v] = edgeW;
}

/* Updates vertex weights */
/*void GraphUpdateVertexWeight(Graph G, vertex v, int vertex_weight) {
	G->vertex_weights[v] = vertex_weight;
}*/

/* Updates memory weights */
void SourceGraph::setMemoryWeight(vertex v, int memoryWeight) {
	memory[v] = memoryWeight;
}

int SourceGraph::getNumberOfVertices() {
	return V;
}

link SourceGraph::getAdjOfVertex(int i) {
	return adj[i];
}

void SourceGraph::setFlags(int vertexLabel, int edgeWeight, int vertexWeight, int memory) {
	enableVertexLabels = vertexLabel;
	enableEdgeWeights = edgeWeight;
	enableVertexWeights = vertexWeight;
	enableMemory = memory;
}

void SourceGraph::printGraphSrc() {
	int i=0;
	link aux;

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
		//cout << "Memory: " << memory[i] << " ";
		cout << "\n";
	}
}


void SourceGraph::printGraph() {
	int i=0;
	link aux;

	cout << "\nSource: V=" << V << ", A=" << A << ", enableMemory=" << enableMemory << ", enableVertexLabels=" << enableVertexLabels << ", enableEdgeWeights=" << enableEdgeWeights << ", enableVertexWeights=" << enableVertexWeights << "\n";

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
		cout << "Memory: " << memory[i] << " ";
		cout << "\n";
	}
}

int SourceGraph::getEnableMemory(){
	return enableMemory;
}

int SourceGraph::getMemory(vertex v){
	return memory[v];
}

// The parameter is a pointer to a pointer to the struct graph 
SourceGraph::~SourceGraph() {
	link no=NULL, aux=NULL;
	int i;

	if (memory != NULL) {
		free(memory);	
		memory = NULL;
	}
	if (adjMatrix != NULL) {
		for (i = 0; i < V; i++) {
			if (adjMatrix[i] != NULL) {
				free(adjMatrix[i]);
			}
		}
		free(adjMatrix);
		adjMatrix = NULL;
	}
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

