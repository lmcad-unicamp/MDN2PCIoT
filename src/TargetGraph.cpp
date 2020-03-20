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
TargetGraph::TargetGraph(int n, int nThreads) { 
   vertex v;

   V = n; 
   A = 0;
   numberOfThreads = nThreads;
   computationalWeight = (int *) malloc(V * sizeof(int));
   if (computationalWeight == NULL) {
		std::cout << "computationalWeight could not be allocated! \n";
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
   assignedPerThread = (int **) malloc(numberOfThreads * sizeof(int **));
   if (assignedPerThread == NULL) {
		std::cout << "assignedPerThread could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   for (v = 0; v < numberOfThreads; v++) {
      assignedPerThread[v] = (int *) malloc(V * sizeof(int ));
      if (assignedPerThread[v] == NULL) {
         std::cout << "assignedPerThread[v] could not be allocated! \n";
         exit(EXIT_FAILURE);
   	  }
   }
   connectionMatrix = (int **) malloc(V * sizeof(int *));
   if (connectionMatrix == NULL) {
		std::cout << "connectionMatrix could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   for (v = 0; v < V; v++) {
	connectionMatrix[v] = (int *) malloc(V * sizeof(int ));
	if (connectionMatrix[v] == NULL) {
		std::cout << "connectionMatrix[v] could not be allocated! \n";
		exit(EXIT_FAILURE);
   	}
   }
   adj = (targetLink *) malloc(V * sizeof (targetLink));
   if (adj == NULL) {
		std::cout << "adj could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   for (v = 0; v < V; ++v) {
	computationalWeight[v] = 1;
	adj[v] = NULL;
	memory[v] = 1;
	assigned[v] = -1;
	for (int t = 0; t < nThreads; t++)
	    assignedPerThread[t][v] = -1;
	for (int i = 0; i < V; i++) {
		connectionMatrix[v][i] = 1;
		if (v == i) {
			connectionMatrix[v][i] = 0;
		}
	}
   }
}

void TargetGraph::insertArc(vertex v, vertex w) { 

}

/* Updates vertex weights */
/*void TargetGraph::targetsetComputationalWeight(vertex v, int computationalWeight) {
	computational_weights[v] = computationalWeight;
}*/

/* Updates memory weights */
void TargetGraph::setMemoryWeight(vertex v, int memoryWeight) {
	memory[v] = memoryWeight;
}

/* Updates computational weights */
void TargetGraph::setComputationalWeight(vertex v, int compWeight) {
	computationalWeight[v] = compWeight;
}

/* Updates connection weights */
void TargetGraph::setConnectionWeight(vertex src, vertex dst, int connWeight) {
	connectionMatrix[src][dst] = connWeight;
	connectionMatrix[dst][src] = connWeight;
}

int TargetGraph::getMemory(int machine) const {
	return memory[machine];
}

/* Updates assigned */
void TargetGraph::setAssigned(int partition, int device) {
	assigned[partition] = device;
	for (int t = 0; t < numberOfThreads; t++)
		assignedPerThread[t][partition] = device;
}

void TargetGraph::setAssignedPerThread(int partition, int device, int thread) {
	assignedPerThread[thread][partition] = device;
}

int TargetGraph::getComputationalWeight(int machine) const {
	return computationalWeight[machine];
}

int TargetGraph::getAssigned(int machine) const {
	return assigned[machine];
}

int TargetGraph::getConnectionWeight(int src, int dst) const {
	return connectionMatrix[src][dst];
}

void TargetGraph::printGraph() const {

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
		std::cout << "Computational power: " << computationalWeight[i] << ", Memory: " << memory[i];
		std::cout << "\n";
	}
	std::cout << "Connection matrix: \n";
	for (i = 1; i < V; i++) {
		for (int j = 0; j < i; j++) {
			std::cout << connectionMatrix[i][j] << " ";
		}
		std::cout << "\n";
	}	
}

int TargetGraph::getNumberOfVertices() const {
	return V;
}

// copy constructor
TargetGraph::TargetGraph(const TargetGraph &graphToCopy) 
	: V(graphToCopy.V), A(graphToCopy.A), numberOfThreads(graphToCopy.numberOfThreads) 
{
   computationalWeight = (int *) malloc(V * sizeof(int));
   if (computationalWeight == NULL) {
		std::cout << "computationalWeight could not be allocated! \n";
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
   assignedPerThread = (int **) malloc(numberOfThreads * sizeof(int **));
   if (assignedPerThread == NULL) {
		std::cout << "assignedPerThread could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   for (int v = 0; v < numberOfThreads; v++) {
      assignedPerThread[v] = (int *) malloc(V * sizeof(int ));
      if (assignedPerThread[v] == NULL) {
         std::cout << "assignedPerThread[v] could not be allocated! \n";
         exit(EXIT_FAILURE);
   	  }
   }
   connectionMatrix = (int **) malloc(V * sizeof(int *));
   if (connectionMatrix == NULL) {
		std::cout << "connectionMatrix could not be allocated! \n";
		exit(EXIT_FAILURE);
   }
   for (int v = 0; v < V; v++) {
	connectionMatrix[v] = (int *) malloc(V * sizeof(int ));
	if (connectionMatrix[v] == NULL) {
		std::cout << "connectionMatrix[v] could not be allocated! \n";
		exit(EXIT_FAILURE);
   	}
   }
   adj = (targetLink *) malloc(V * sizeof (targetLink));
   if (adj == NULL) {
		std::cout << "adj could not be allocated! \n";
		exit(EXIT_FAILURE);
   }

   for (int v = 0; v < V; ++v) {
		computationalWeight[v] = graphToCopy.computationalWeight[v];
		memory[v] = graphToCopy.memory[v];
		assigned[v] = graphToCopy.assigned[v];
		for (int t = 0; t < numberOfThreads; t++) 
			assignedPerThread[t][v] = graphToCopy.assignedPerThread[t][v];
		for (int i = 0; i < V; i++) {
			connectionMatrix[v][i] = graphToCopy.connectionMatrix[v][i];
		}
		adj[v] = NULL;
		targetLink no = graphToCopy.adj[v];
		while (no != NULL) {
			insertArc(v, no->w);
			no = no->next;
		}
   }
}

TargetGraph::~TargetGraph() {
	targetLink no=NULL, aux=NULL;
	int i;

	if (computationalWeight != NULL) {
		free(computationalWeight);
		computationalWeight = NULL;
	}
	if (memory != NULL) {
		free(memory);
		memory = NULL;
	}
	if (assigned != NULL) {
		free(assigned);
		assigned = NULL;
	}
	if (assignedPerThread != NULL) {
		for (int i = 0; i < numberOfThreads; i++) {
			if (assignedPerThread[i] != NULL) {
				free(assignedPerThread[i]);
			}
		}
		free(assignedPerThread);
		assignedPerThread = NULL;
	}
	if (connectionMatrix != NULL) {
		for (int i = 0; i < V; i++) {
			if (connectionMatrix[i] != NULL) {
				free(connectionMatrix[i]);
			}
		}
		free(connectionMatrix);
		connectionMatrix = NULL;
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
