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

// https://www.ime.usp.br/~pf/algoritmos_para_grafos/aulas/graphdatastructs.html
#ifndef GRAPH
#define GRAPH

/* Graph vertices are represented by objects of type vertex. */
#define vertex int

/* The adjacency list of a vertex v is composed by nodes of type node. Each node of the list corresponds to an arc and contains a neighbour w of v and the address of the next node in the list. A link is a pointer to a node. */
typedef struct node *link;

class SourceGraph {

public:
	/* ADJACENCY LIST REPRESENTATION: the function GRAPHinit() builds a graph with vertices 0 1 .. V-1 and no arc. */
	SourceGraph(int n, int nPartitions = 2);

	/* The function NEWnode() receives a vertex w and the address next of a node and returns the address of a new node such that a->w == w and a->next == next. */
	//static link NEWnode(vertex w, link next, int edgeW);
	link NEWnode(vertex w, link next, int edgeW, int src=1);
	/* ADJACENCY LIST REPRESENTATION: the function GRAPHinsertArc() inserts an arc v-w in graph G. The function supposes that v and w are distinct, positive, and smaller than G->V. If the graph already has an arc v-w, the function does not do anything. */
	void insertArc(vertex v, vertex w, int edgeW, int src=1);
	bool insertArcSrc(vertex v, vertex w, int edgeW, int src);
	/* Updates vertex weights */
	//void GraphUpdateVertexWeight(Graph G, vertex v, int vertex_weight);
	/* Updates memory weights */
	void setMemoryWeight(vertex v, int memoryWeight);
	/* Updates flags */
	int getNumberOfVertices();
	link getAdjOfVertex(int i);
	void setFlags(int vertexLabel, int edgeWeight, int vertexWeight, int memory);
	void printGraph();
	void printGraphSrc();
	int getEnableMemory();
	int getMemory(vertex v);
	bool removeArcSrc(vertex v, vertex w, int src);
	void setAdjMatrix(vertex v, vertex w, double edgeW);
	~SourceGraph();

private:
	int V;
	int A; 
	int enableVertexLabels;
	int enableEdgeWeights;
	int enableVertexWeights;
	int enableMemory;
	int *vertexWeights;
	int *memory;
	double **adjMatrix;
	link *adj;
};

#endif

#ifndef NODE
#define NODE
struct node { 
	vertex w; 
	int edgeWeight;
	int source;
	link next; 
};
#endif


