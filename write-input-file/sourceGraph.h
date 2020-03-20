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


/* Graph vertices are represented by objects of type vertex. */
#define vertex int

/* The adjacency list of a vertex v is composed by nodes of type node. Each node of the list corresponds to an arc and contains a neighbour w of v and the address of the next node in the list. A link is a pointer to a node. */
typedef struct node *link;

#ifndef GRAPH
#define GRAPH
struct graph {
   int V; 
   int A; 
   int enable_vertex_labels, enable_edge_weights, enable_vertex_weights, enable_memory;
   int *vertex_weights;
   int *memory;
   link *adj;
};
#endif

typedef struct graph *Graph;

#ifndef NODE
#define NODE
struct node { 
   vertex w; 
   int edge_weight;
   link next; 
};
#endif

static link NEWnode(vertex w, link next, int edge_weight);

Graph GraphInit(int V);

void GraphInsertArc(Graph G, vertex v, vertex w, int edge_weight);

/* Updates vertex weights */
void GraphUpdateVertexWeight(Graph G, vertex v, int vertex_weight);

/* Updates memory weights */
void GraphUpdateMemoryWeight(Graph G, vertex v, int memory_weight);

/* Updates flags */
void GraphUpdateFlags(Graph G, int enable_vertex_labels, int enable_edge_weights, int enable_vertex_weights, int enable_memory);

void PrintGraph(Graph G);

void FreeGraph(Graph *G);
