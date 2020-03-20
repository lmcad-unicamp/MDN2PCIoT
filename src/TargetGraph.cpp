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
#include "TargetGraph.h"

	/* ADJACENCY LIST REPRESENTATION: the function GRAPHinit() builds a graph with vertices 0 1 .. V-1 and no arc. */
TargetGraph::TargetGraph(int n) { 
   vertex v;

   V = n; 
   A = 0;
   computationalWeights = (int *) malloc(V * sizeof(int));
   if (computationalWeights == NULL) {
		std::cout << "computationalWeights could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   memory = (int *) malloc(V * sizeof(int));
   if (memory == NULL) {
		std::cout << "memory could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   assigned = (int *) malloc(V * sizeof(int));
   if (assigned == NULL) {
		std::cout << "assigned could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   adj = (targetLink *) malloc(V * sizeof (targetLink));
   if (adj == NULL) {
		std::cout << "adj could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   for (v = 0; v < V; ++v) {
	  computationalWeights[v] = 1;
      adj[v] = NULL;
	  memory[v] = 1;
	  assigned[v] = -1;
   }
}

/* Updates memory weights */
void TargetGraph::setMemoryWeight(vertex v, int memoryWeight) {
	memory[v] = memoryWeight;
}

int TargetGraph::getMemory(int machine) {
	return memory[machine];
}

void TargetGraph::printGraph() {
	int i=0;
	targetLink aux;

	std::cout << "\nTarget: V=" << V << "\n";

	for (i = 0; i < V; i++) {
		std::cout << "Vertex[" << i << "]: "; 
		aux = adj[i];
		while(aux != NULL){
			std::cout << aux->w << " ";
			aux = aux->next;
		}
		std::cout << "Computational power: " << computationalWeights[i] << ", Memory: " << memory[i];
		std::cout << "\n";
	}
}

int TargetGraph::getNumberOfVertices() {
	return V;
}

TargetGraph::~TargetGraph() {
	targetLink no=NULL, aux=NULL;
	int i;

	if (computationalWeights != NULL) {
		free(computationalWeights);
		computationalWeights = NULL;
	}
	if (memory != NULL) {
		free(memory);
		memory = NULL;
	}
	if (assigned != NULL) {
		free(assigned);
		assigned = NULL;
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
