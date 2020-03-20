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
#include <malloc.h>
#include "sourceGraph.h"

/* The function NEWnode() receives a vertex w and the address next of a node and returns the address of a new node such that a->w == w and a->next == next. */
static link NEWnode(vertex w, link next, int edge_weight) { 
   link a = (link) malloc(sizeof (struct node));

   a->w = w; 
   a->edge_weight = edge_weight;
   a->next = next;     
   return a;                         
}

/* ADJACENCY LIST REPRESENTATION: the function GRAPHinit() builds a graph with vertices 0 1 .. V-1 and no arc. */
Graph GraphInit(int V) { 
   vertex v;
   Graph G=NULL;
   int odd=0;
 
   G = (Graph) malloc(sizeof *G);
   if (G == NULL) {
		printf("G could not be allocated! \n");
		exit(EXIT_FAILURE);
   }
   /*if ((V % 2) != 0) {
		V++;
		odd = 1;
   }*/
   G->V = V; 
   G->A = 0;
   G->enable_vertex_labels = 0;
   G->enable_edge_weights = 0;
   G->enable_vertex_weights = 0;
   G->vertex_weights = (int *) malloc(V * sizeof(int));
   if (G->vertex_weights == NULL) {
		printf("G->vertex_weights could not be allocated! \n");
		exit(EXIT_FAILURE);
   }
   G->memory = (int *) malloc(V * sizeof(int));
   if (G->memory == NULL) {
		printf("G->memory could not be allocated! \n");
		exit(EXIT_FAILURE);
   }
   G->adj = (link *) malloc(V * sizeof (link));
   if (G->adj == NULL) {
		printf("G->adj could not be allocated! \n");
		exit(EXIT_FAILURE);
   }
   for (v = 0; v < V; ++v) {
	  G->vertex_weights[v] = 1;
      G->adj[v] = NULL;
	  G->memory[v] = 1;
   }
   if (odd == 1) {
      G->memory[V-1] = 0;
   }
   return G;
}

/* ADJACENCY LIST REPRESENTATION: the function GRAPHinsertArc() inserts an arc v-w in graph G. The function supposes that v and w are distinct, positive, and smaller than G->V. If the graph already has an arc v-w, the function does not do anything. */
void GraphInsertArc(Graph G, vertex v, vertex w, int edge_weight) { 
   link a;

   if (w == -1) {
	   return;
   }
   for (a = G->adj[v]; a != NULL; a = a->next) 
      if (a->w == w) return;
   G->adj[v] = NEWnode(w, G->adj[v], edge_weight);
   G->A++;
}

/* Updates vertex weights */
void GraphUpdateVertexWeight(Graph G, vertex v, int vertex_weight) {
	G->vertex_weights[v] = vertex_weight;
}

/* Updates memory weights */
void GraphUpdateMemoryWeight(Graph G, vertex v, int memory_weight) {
	G->memory[v] = memory_weight;
}

void GraphUpdateFlags(Graph G, int vl, int ew, int vw, int m) {
	G->enable_vertex_labels = vl;
	G->enable_edge_weights = ew;
	G->enable_vertex_weights = vw;
	G->enable_memory = m;
}

void PrintGraph(Graph G) {
	int i=0;
	link aux;

	printf("\nSource: V=%d, A=%d, em=%d, evl=%d, eew=%d, evw=%d \n", G->V, G->A, G->enable_memory, G->enable_vertex_labels, G->enable_edge_weights, G->enable_vertex_weights);

	for (i = 0; i < G->V; i++) {
		printf("Vertex[%d]: ", i);
		aux = G->adj[i];
		while(aux != NULL){
			printf("%d ", aux->w);
			aux = aux->next;
		}
		printf("Memory: %d ", G->memory[i]);
		if (G->enable_vertex_weights == 1) {
			printf("Computational weight: %d ", G->vertex_weights[i]);
		}
		if (G->enable_edge_weights == 1) {
			printf("Edge weights: ");
			aux = G->adj[i];
			while(aux != NULL) {
				printf("%d ", aux->edge_weight);
				aux = aux->next;
			}
		}
		printf("\n");
	}
}

// The parameter is a pointer to a pointer to the struct graph
void FreeGraph(Graph *G) {
	link no=NULL, aux=NULL;
	int i;

	if (*G != NULL) { // *G represents the content to the double pointer, that is, the address pointed by G (which is a pointer to graph)
		if ((*G)->vertex_weights != NULL){
			free((*G)->vertex_weights);
			(*G)->vertex_weights = NULL;
		}
		if ((*G)->memory != NULL){
			free((*G)->memory);
			(*G)->memory = NULL;
		}
		if ((*G)->adj != NULL){		
			for (i = 0; i < (*G)->V; i++) {
				no = (*G)->adj[i];
				while (no != NULL) {
					aux = no;
					no = no->next;
					free(aux);
					aux = NULL;
				}
			}
			no = NULL;
			free((*G)->adj);
			(*G)->adj = NULL;
		}
		free(*G);
		*G = NULL; 
	}
}
