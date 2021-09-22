/***************************************************************************
 *   Copyright (C) 2020 by Fabíola Martins Campos de Oliveira 		   *
 *   fabiola.bass@gmail.com			                           *
 *                       						   *
 *   This file is part of MDN²PCIoT.  					   *
 *                                      		   		   *
 *   MDN²PCIoT is free software: you can redistribute it and/or modify	   *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.				   *
 *									   *
 *   MDN²PCIoT is distributed in the hope that it will be useful,	   *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of	   *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	   *
 *   GNU General Public License for more details.			   *
 *									   *
 *   You should have received a copy of the GNU General Public License     *
 *   long with MDN²PCIoT.  If not, see <http://www.gnu.org/licenses/>.     *
 ***************************************************************************/

#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include "SourceGraph.h"
#include "CoarsenGraph.h"

using namespace std;

CoarsenGraph::CoarsenGraph(SourceGraph *sourceG, TargetGraph *tgtG, bool verb, const int *partitioning, int definedNumberOfCoarsenedGraphs) : target(tgtG), verb(verb) 
{
	numberOfCoarsenedGraphs = 1;
	int newNumberOfVertices = sourceG->getNumberOfVertices();

	// multilevel graph coarsening (for small number of vertices)
	if (sourceG->getNumberOfVertices() < 50) {
		while (newNumberOfVertices > tgtG->getNumberOfVertices()) {
			newNumberOfVertices /= 2;
			numberOfCoarsenedGraphs++;
		}
	} else if (sourceG->getNumberOfVertices() < 2500) {
		// multilevel graph coarsening (for large number of vertices)
		while (newNumberOfVertices > 50) {
			newNumberOfVertices /= 2;
			numberOfCoarsenedGraphs++;
		}
	} else {
		// number of graphs: c < numberOfCoarsenedGraphs - 1 
		numberOfCoarsenedGraphs = definedNumberOfCoarsenedGraphs;
	}
	//cout << "\nnumber; " << numberOfCoarsenedGraphs;

	// turn off multilevel
	//numberOfCoarsenedGraphs = 1;

	if (tgtG->getNumberOfVertices() >= 4 && tgtG->getNumberOfVertices() <= 11 && sourceG->getNumberOfVertices() < 700) {
		granularityBalanceFactor = 0.03125;
	} else if (sourceG->getNumberOfVertices() < 2500){
		granularityBalanceFactor = 0.25;		
	} else {
		granularityBalanceFactor = 0.25;
	}

	// 32 devices in the case of AlexNet
	if (tgtG->getNumberOfVertices() < 32)
		vertexGroupingLimit = granularityBalanceFactor * target->getMinMemory();
	else 
		// equivalent to METIS in the case of AlexNet (sum of the vertex sizes, which is equivalent to the sum of the memory required to store layer output data
		vertexGroupingLimit = 833265; 

	sourceCoarsenedGraph = (SourceGraph **) malloc(numberOfCoarsenedGraphs * sizeof(SourceGraph *));
	if (sourceCoarsenedGraph == NULL) {
		cout << "\n sourceCoarsenedGraph could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	sourceCoarsenedGraph[0] = sourceG;

	sortedDegree = (int *) malloc(sourceG->getNumberOfVertices() * sizeof(int));
	if (sortedDegree == NULL) {
		cout << "\n sortedDegree could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	sortedVerticesByDegree = (int *) malloc(sourceG->getNumberOfVertices() * sizeof(int));
	if (sortedVerticesByDegree == NULL) {
		cout << "\n sortedVerticesByDegree could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < sourceG->getNumberOfVertices(); i++) {
		sortedDegree[i] = sourceG->getVertexDegree(i);
		sortedVerticesByDegree[i] = i;
	}

	// Sort the graph vertices by the vertex degree
	sortDegree(0);

	p = (int *) malloc(sourceG->getNumberOfVertices() * sizeof(int));
	if (p == NULL) {
		cout << "\n p could not be allocated! \n";
		exit(EXIT_FAILURE);	
	}
	for (int i = 0; i < sourceG->getNumberOfVertices(); i++) {
		p[i] = partitioning[i];
	}

	// Coarsen the source graph
	coarsenSourceGraph();
}

// Sort the graph vertices by the vertex degree
void CoarsenGraph::sortDegree(int c) {
	mergeSortWithIndexes(sortedDegree, 1, sourceCoarsenedGraph[c]->getNumberOfVertices(), sortedVerticesByDegree);

	/*for (int i = 0; i < sourceCoarsenedGraph[0]->getNumberOfVertices(); i++) {
		cout << "\n sD: " << sortedDegree[i] << ", sVBD: " << sortedVerticesByDegree[i];
	}*/
}

// Seed the random number generator.
void CoarsenGraph::seed_rng(void)
{
    int fp = open("/dev/random", O_RDONLY);
    if (fp == -1) abort();
    unsigned seed;
    unsigned pos = 0;
    while (pos < sizeof(seed)) {
        int amt = read(fp, (char *) &seed + pos, sizeof(seed) - pos);
        if (amt <= 0) abort();
        pos += amt;
    }
    srand(seed);
    close(fp);
}

void CoarsenGraph::coarsenSourceGraph() {

	//match = new int[sourceG->getNumberOfVertices];

	for (int c = 0; c < numberOfCoarsenedGraphs - 1; c++) {

		cout << "\nGraph number:" << c;

		seed_rng();

		if (c > 0) {
			sortedDegree = (int *) malloc(sourceCoarsenedGraph[c]->getNumberOfVertices() * sizeof(int));
			if (sortedDegree == NULL) {
				cout << "\n sortedDegree could not be allocated! \n";
				exit(EXIT_FAILURE);	
			}
			sortedVerticesByDegree = (int *) malloc(sourceCoarsenedGraph[c]->getNumberOfVertices() * sizeof(int));
			if (sortedVerticesByDegree == NULL) {
				cout << "\n sortedVerticesByDegree could not be allocated! \n";
				exit(EXIT_FAILURE);	
			}
			for (int i = 0; i < sourceCoarsenedGraph[c]->getNumberOfVertices(); i++) {
				sortedDegree[i] = sourceCoarsenedGraph[c]->getVertexDegree(i);
				sortedVerticesByDegree[i] = i;
			}

			// Sort the graph vertices by the vertex degree
			sortDegree(c);
		}

		int map = 0, i, numberOfNonmatched = 0;

		// Use 2-hop matching first if you want initial partitionings more close to the greedy approach or best fit
		// try a 2-hop matching with the unmatched vertices in the case that the number of vertices that were not matched is greater than a threshold
		// Scan for unmatched vertices
		/*for (int k = 0; k < sourceCoarsenedGraph[c]->getNumberOfVertices(); k++) {
			// By increasing order of the degree of the vertices
			int i = k; //sortedVerticesByDegree[k];

			// check if i is not matched
			if (sourceCoarsenedGraph[c]->getMatch(i) == i) {

				// Scan for another unmatched vertex
				for (int m = 0; m < sourceCoarsenedGraph[c]->getNumberOfVertices(); m++) {
					// By increasing order of the degree of the vertices	
					int j = m; //sortedVerticesByDegree[m];

					if (i == j) continue;

					// check if j is not matched
					if (sourceCoarsenedGraph[c]->getMatch(j) == j /*&& p[i] == p[j]*//*) {

						// Scan for a common edge between vertices i and j
						for (links nodeAdji = sourceCoarsenedGraph[c]->getAdjOfVertex(i); nodeAdji != NULL; nodeAdji = nodeAdji->next) {
							for (links nodeAdjj = sourceCoarsenedGraph[c]->getAdjOfVertex(j); nodeAdjj != NULL; nodeAdjj = nodeAdjj->next) {
								if (nodeAdji->w == nodeAdjj->w) {
									// check if vertex i and j fits the maximum memory for a vertex
									// shared parameters
									int shared = 0;
									for (int layer = 0; layer < sourceCoarsenedGraph[c]->getNumberOfLayers(); layer++) {
										if (sourceCoarsenedGraph[c]->getNodeSharedParam(i, layer) || sourceCoarsenedGraph[c]->getNodeSharedParam(j, layer)) {
											shared += sourceCoarsenedGraph[c]->getSharedParam(layer); 
										}
									}

									// check if the matched vertex fits the memory of the device with the least amount of memory
									if (sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(j) + shared <= (vertexGroupingLimit)) {
										// match vertices i and j
										sourceCoarsenedGraph[c]->setMatch(i, j);
										sourceCoarsenedGraph[c]->setMap(i, j, map);
										if (verb == true) {
											cout << "\ni: " << i << ", j: " << j << ", map: " << map << ", mem: " << sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(j) + shared;
										}
										map++;
										break;
									}
								}
								// stop search if vertices i and j are matched
								if (sourceCoarsenedGraph[c]->getMatch(j) != j) {
									break;
								}
							}
							// stop search if vertices i and j are matched
							if (sourceCoarsenedGraph[c]->getMatch(i) != i) {
								break;
							}
						}
						// stop search if vertices i and j are matched
						if (sourceCoarsenedGraph[c]->getMatch(i) != i) {
							break;
						}
					}
				}
			}
		}*/	

		// find a maximal matching in the source graph
		for (int k = 0; k < sourceCoarsenedGraph[c]->getNumberOfVertices(); k++) {

			// random choice of matchings
			/*i = rand() % sourceCoarsenedGraph[c]->getNumberOfVertices();
			for (int j = 0; j < sourceCoarsenedGraph[c]->getNumberOfVertices(); j++) {	
				if (sourceCoarsenedGraph[c]->getMatch(i) != i) {
					if (i == sourceCoarsenedGraph[c]->getNumberOfVertices() - 1) {
						i = 0;
					} else {
						i++;
					}
				} else {
					break;
				}
			}*/
			// nonrandom
			//i = k;

			i = k; //sortedVerticesByDegree[k];	
			
			// find a maximal matching (with the largest edge weight for vertex i in the source graph)
			if (sourceCoarsenedGraph[c]->getMatch(i) == i) {
				int maxAdj = -1;
				int max = 0;
				for (links nodeAdj = sourceCoarsenedGraph[c]->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
					//cout << "\n i: " << i << ", p[i]: " << p[i] << ", w: " << nodeAdj->w << ", p[nodeAdj->w]: " << p[nodeAdj->w];
					if (sourceCoarsenedGraph[c]->getMatch(nodeAdj->w) == nodeAdj->w && p[i] == p[nodeAdj->w]) {
						// shared parameters
						int shared = 0;
						for (int layer = 0; layer < sourceCoarsenedGraph[c]->getNumberOfLayers(); layer++) {
							if (sourceCoarsenedGraph[c]->getNodeSharedParam(i, layer) || sourceCoarsenedGraph[c]->getNodeSharedParam(nodeAdj->w, layer)) {
								shared += sourceCoarsenedGraph[c]->getSharedParam(layer); 
							}
						}
//cout << "\n i: " << i << ", p[i]: " << p[i] << ", w: " << nodeAdj->w << ", p[nodeAdj->w]: " << p[nodeAdj->w] << ", memi: " << sourceCoarsenedGraph[c]->getMemory(i) << ", memadj: " << sourceCoarsenedGraph[c]->getMemory(nodeAdj->w) << ", shared: " << shared << ", cond: " << (vertexGroupingLimit);

						// check if the matched vertex fits the memory of the device with the least amount of memory
						if (sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(nodeAdj->w) + shared <= (vertexGroupingLimit)) {
							if (max < nodeAdj->edgeWeight) {
								maxAdj = nodeAdj->w;
								max = nodeAdj->edgeWeight;
							}
						}
					}
				}

				for (links nodeAdj = sourceCoarsenedGraph[c]->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
					if (sourceCoarsenedGraph[c]->getMatch(nodeAdj->w) == nodeAdj->w && p[i] == p[nodeAdj->w]) { // não comentar p aqui se quiser inicial para 11 subgrafos igual a 13672 e atual em 13544 no desktop
						// shared parameters
						int shared = 0;
						for (int layer = 0; layer < sourceCoarsenedGraph[c]->getNumberOfLayers(); layer++) {
							if (sourceCoarsenedGraph[c]->getNodeSharedParam(i, layer) || sourceCoarsenedGraph[c]->getNodeSharedParam(nodeAdj->w, layer)) {
								shared += sourceCoarsenedGraph[c]->getSharedParam(layer); 
							}
						}

						// check if the matched vertex fits the memory of the device with the least amount of memory
						if (sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(nodeAdj->w) + shared <= (vertexGroupingLimit)) {
							if (nodeAdj->w == maxAdj && nodeAdj->edgeWeight == max) {	
								sourceCoarsenedGraph[c]->setMatch(i, nodeAdj->w);
								sourceCoarsenedGraph[c]->setEdgeMatch(nodeAdj);
								sourceCoarsenedGraph[c]->setMap(i, nodeAdj->w, map);	
								if (verb == true) {
									cout << "\ni: " << i << ", adj: " << nodeAdj->w << ", edgeW: " << nodeAdj->edgeWeight << ", srcMlvl: " << nodeAdj->sourceMlvl << ", map: " << map << ", mem: " << sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(nodeAdj->w) + shared;
								}
								map++;
								break;
							}
						}
					}
				}
			}
			if (sourceCoarsenedGraph[c]->getMatch(i) == i) {
				numberOfNonmatched++;
			}
		}

		// Comment if you will use 2-hop-first matching
		// try a 2-hop matching with the unmatched vertices in the case that the number of vertices that were not matched is greater than a threshold
		if (numberOfNonmatched > 0.10 * sourceCoarsenedGraph[c]->getNumberOfVertices()) {
			// Scan for unmatched vertices
			for (int k = 0; k < sourceCoarsenedGraph[c]->getNumberOfVertices(); k++) {
				// By increasing order of the degree of the vertices
				int i = k; //sortedVerticesByDegree[k];

				// check if i is not matched
				if (sourceCoarsenedGraph[c]->getMatch(i) == i) {

					// Scan for another unmatched vertex
					for (int m = 0; m < sourceCoarsenedGraph[c]->getNumberOfVertices(); m++) {
						// By increasing order of the degree of the vertices	
						int j = m; //sortedVerticesByDegree[m];

						if (i == j) continue;

						// check if j is not matched
						if (sourceCoarsenedGraph[c]->getMatch(j) == j && p[i] == p[j]) {

							// Scan for a common edge between vertices i and j
							for (links nodeAdji = sourceCoarsenedGraph[c]->getAdjOfVertex(i); nodeAdji != NULL; nodeAdji = nodeAdji->next) {
								for (links nodeAdjj = sourceCoarsenedGraph[c]->getAdjOfVertex(j); nodeAdjj != NULL; nodeAdjj = nodeAdjj->next) {
									if (nodeAdji->w == nodeAdjj->w) {
										// check if vertex i and j fits the maximum memory for a vertex
										// shared parameters
										int shared = 0;
										for (int layer = 0; layer < sourceCoarsenedGraph[c]->getNumberOfLayers(); layer++) {
											if (sourceCoarsenedGraph[c]->getNodeSharedParam(i, layer) || sourceCoarsenedGraph[c]->getNodeSharedParam(j, layer)) {
												shared += sourceCoarsenedGraph[c]->getSharedParam(layer); 
											}
										}

										// check if the matched vertex fits the memory of the device with the least amount of memory
										if (sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(j) + shared <= (vertexGroupingLimit)) {
											// match vertices i and j
											sourceCoarsenedGraph[c]->setMatch(i, j);
											sourceCoarsenedGraph[c]->setMap(i, j, map);
											if (verb == true) {
												cout << "\ni: " << i << ", j: " << j << ", map: " << map << ", mem: " << sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(j) + shared;
											}
											map++;
											break;
										}
									}
									// stop search if vertices i and j are matched
									if (sourceCoarsenedGraph[c]->getMatch(j) != j) {
										break;
									}
								}
								// stop search if vertices i and j are matched
								if (sourceCoarsenedGraph[c]->getMatch(i) != i) {
									break;
								}
							}
							// stop search if vertices i and j are matched
							if (sourceCoarsenedGraph[c]->getMatch(i) != i) {
								break;
							}
						}
					}
				}
			}
		}	

		// include the unmatched vertices into map
		for (int i = 0; i < sourceCoarsenedGraph[c]->getNumberOfVertices(); i++) {
			if (sourceCoarsenedGraph[c]->getMatch(i) == i) {
				sourceCoarsenedGraph[c]->setMap(i, i, map);
				/*if (verb == true) {
					cout << "\ni: " << i << ", map: " << map << ", mem: " << sourceCoarsenedGraph[c]->getMemory(i);
				}*/
				map++;
			}
		}

		// build a coarsen graph
		sourceCoarsenedGraph[c + 1] = new SourceGraph(map, sourceCoarsenedGraph[c]->getNumberOfLayers(), target->getNumberOfVertices());
		sourceCoarsenedGraph[c + 1]->setFlags(sourceCoarsenedGraph[c]->getEnableVertexLabels(), sourceCoarsenedGraph[c]->getEnableEdgeWeights(), sourceCoarsenedGraph[c]->getEnableVertexWeights(), sourceCoarsenedGraph[c]->getEnableMemory());
		int shrdppl[sourceCoarsenedGraph[c]->getNumberOfLayers()];	
		sourceCoarsenedGraph[c]->getSharedParamPerLayer(shrdppl);
		sourceCoarsenedGraph[c + 1]->setSharedParamPerLayer(shrdppl);

		bool shrdppv[sourceCoarsenedGraph[c]->getNumberOfLayers()];
		for (int i = 0; i < sourceCoarsenedGraph[c]->getNumberOfVertices(); i++) {
			// unmatched vertices
			if (sourceCoarsenedGraph[c]->getMatch(i) == i) {
				/*cout << "\n i: " << i;
				cout << "\n2 memi: " << sourceCoarsenedGraph[c]->getMemory(i);*/
				sourceCoarsenedGraph[c + 1]->setMemoryWeight(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getMemory(i));
				sourceCoarsenedGraph[c + 1]->setVertexWeight(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getVertexWeight(i));

				sourceCoarsenedGraph[c]->getSharedParamArray(i, shrdppv);
				sourceCoarsenedGraph[c + 1]->setNodeSharedParamArray(sourceCoarsenedGraph[c]->getMap(i), shrdppv);

				p[sourceCoarsenedGraph[c]->getMap(i)] = p[i];
			}
			for (links nodeAdj = sourceCoarsenedGraph[c]->getAdjOfVertex(i); nodeAdj != NULL; nodeAdj = nodeAdj->next) {
				//cout << "\n i: " << i << ", adj: " << nodeAdj->w;
				//cout << ", matchi: " <<  sourceCoarsenedGraph[c]->getMatch(i) << ", matchadj: " << sourceCoarsenedGraph[c]->getMatch(nodeAdj->w);

				// for each matching edge
				// update hipervertex data
				if (sourceCoarsenedGraph[c]->getMatch(i) == nodeAdj->w && sourceCoarsenedGraph[c]->getMatch(nodeAdj->w) == i /*&& i < nodeAdj->w*/) {
					//cout << "\n1 memi: " << sourceCoarsenedGraph[c]->getMemory(i) << ", memadj: " << sourceCoarsenedGraph[c]->getMemory(nodeAdj->w);
					//cout << "\n1 mem: " << sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(nodeAdj->w);

					sourceCoarsenedGraph[c + 1]->setMemoryWeight(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(nodeAdj->w));
					sourceCoarsenedGraph[c + 1]->setVertexWeight(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getVertexWeight(i) + sourceCoarsenedGraph[c]->getVertexWeight(nodeAdj->w));

					sourceCoarsenedGraph[c]->getSharedParamArray(i, shrdppv);
					sourceCoarsenedGraph[c + 1]->setNodeSharedParamArray(sourceCoarsenedGraph[c]->getMap(i), shrdppv);
					sourceCoarsenedGraph[c]->getSharedParamArray(nodeAdj->w, shrdppv);
					sourceCoarsenedGraph[c + 1]->setNodeSharedParamArray(sourceCoarsenedGraph[c]->getMap(i), shrdppv);

					p[sourceCoarsenedGraph[c]->getMap(i)] = p[i];					

				// or insert the respective edges
				} else if (!(sourceCoarsenedGraph[c]->getMatch(i) == nodeAdj->w && sourceCoarsenedGraph[c]->getMatch(nodeAdj->w) == i && i > nodeAdj->w)) {

					int src;
					if (c == 0)
						src = i;
					else
						src = nodeAdj->sourceMlvl;

					//cout << "\n i: " << i << ", adj: " << nodeAdj->w;
					//cout << "\n3 weight: " << nodeAdj->edgeWeight << ", src: " << src;

					sourceCoarsenedGraph[c + 1]->insertArcSrcMlvl(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getMap(nodeAdj->w), nodeAdj->edgeWeight, src);
				}
			}

			// for 2-hop matchings
			for (int j = i + 1; j < sourceCoarsenedGraph[c]->getNumberOfVertices(); j++) {
				if (sourceCoarsenedGraph[c]->getMatch(i) == j) {
					sourceCoarsenedGraph[c + 1]->setMemoryWeight(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getMemory(i) + sourceCoarsenedGraph[c]->getMemory(j));
					sourceCoarsenedGraph[c + 1]->setVertexWeight(sourceCoarsenedGraph[c]->getMap(i), sourceCoarsenedGraph[c]->getVertexWeight(i) + sourceCoarsenedGraph[c]->getVertexWeight(j));

					sourceCoarsenedGraph[c]->getSharedParamArray(i, shrdppv);
					sourceCoarsenedGraph[c + 1]->setNodeSharedParamArray(sourceCoarsenedGraph[c]->getMap(i), shrdppv);
					sourceCoarsenedGraph[c]->getSharedParamArray(j, shrdppv);
					sourceCoarsenedGraph[c + 1]->setNodeSharedParamArray(sourceCoarsenedGraph[c]->getMap(i), shrdppv);
	
					p[sourceCoarsenedGraph[c]->getMap(i)] = p[i];
					break;
				}
			}
		}
		if (verb == true) {
			sourceCoarsenedGraph[c + 1]->printGraphMlvl();
		} else {
			sourceCoarsenedGraph[c + 1]->printGraphHeader();
		}

		if (sortedDegree != NULL)
			free(sortedDegree);
		if (sortedVerticesByDegree != NULL)
			free(sortedVerticesByDegree);

		/*cout << "\n";
		for (int i = 0; i < map; i++) {
			cout << " " << p[i];
		}*/
	}
}

void CoarsenGraph::mergeSortWithIndexes(int *A, int p, int r, int *indexes){
	int j;

	//static int **Ai = (int **) calloc(r, sizeof(int *));
	int **Ai = (int **) calloc(r, sizeof(int *));
	if (Ai == NULL) {
		printf("Ai could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
/**#	pragma omp parallel num_threads(getNumberOfThreads())
	{
#	pragma omp for*/
	for (j = 0; j < r; j++) {
		Ai[j] = (int *) calloc(2, sizeof(int));
		if (Ai[j] == NULL) {
			printf("Ai[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}

//#	pragma omp for
	for (j = 0; j < r; j++) {
		Ai[j][0] = j;
		Ai[j][1] = A[j];
	}
	//} // pragma

	/*for (j = p; j <= r; j++) {
		printf("\n%d %d", Ai[j - 1][0], Ai[j - 1][1]);
	}
	printf("\n");*/

	mergeSort(Ai, p, r);

/*#	pragma omp parallel num_threads(getNumberOfThreads())
	{
#	pragma omp for*/
	for (j = p; j <= r; j++) {
		A[j - 1] = Ai[j - 1][1];
		indexes[j - 1] = Ai[j - 1][0];
	}

	/*for (j = p; j <= r; j++) {
		printf("\n%d %d", Ai[j - 1][0], Ai[j - 1][1]);
	}
	printf("\n");*/

//#	pragma omp for
	for (j = 0; j < r; j++) {
		free(Ai[j]);
		Ai[j] = NULL;
	}	
	//} // pragma

	free(Ai);
	Ai=NULL;
}

void CoarsenGraph::mergeSort(int **A, int p, int r) {
	int j, q;

	if (p < r) {
		//q = floor((p+r)/2);
		q = (p + r)/2;
		mergeSort(A, p, q);
		mergeSort(A, q + 1, r);
		merge(A, p, q, r);
	}
}

void CoarsenGraph::merge(int **A, int p, int q, int r) {
	int n1, n2, **L=NULL, **R=NULL, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	L = (int **) calloc(n1 + 1, sizeof(int *));
	if (L == NULL) {
		printf("L could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < n1 + 1; j++) {
		L[j] = (int *) calloc(2, sizeof(int));
		if (L[j] == NULL) {
			printf("L[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}
	R = (int **) calloc(n2 + 1, sizeof(int *));
	if (R == NULL) {
		printf("R could not be allocated!\n");
		exit(EXIT_FAILURE);
	}
	for (j = 0; j < n2 + 1; j++) {
		R[j] = (int *) calloc(2, sizeof(int));
		if (R[j] == NULL) {
			printf("R[%d] could not be allocated!\n", j);
			exit(EXIT_FAILURE);
		}
	}

	for (i = 1; i < n1 + 1; i++) {
		L[i - 1][1] = A[p + i - 2][1];
		L[i - 1][0] = A[p + i - 2][0];
	}
	for (j = 1; j < n2 + 1; j++) {
		R[j - 1][1] = A[q + j - 1][1];
		R[j - 1][0] = A[q + j - 1][0];
	}

	L[n1][1] = 2147483647;
	L[n1][0] = -1;
	R[n2][1] = 2147483647;
	R[n2][0] = -1;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++) {
		if (L[i - 1][1] <= R[j - 1][1]) {
			A[k - 1][1] = L[i - 1][1];
			A[k - 1][0] = L[i - 1][0];
			i++;
		} else {
			A[k - 1][1] = R[j - 1][1];
			A[k - 1][0] = R[j - 1][0];
			j++;
		}
	}

	for (j = 0; j < n1 + 1; j++) {
		free(L[j]);
		L[j] = NULL;
	}	
	free(L);
	L=NULL;

	for (j = 0; j < n2 + 1; j++) {
		free(R[j]);
		R[j] = NULL;
	}	
	free(R);
	R=NULL;
}

CoarsenGraph::~CoarsenGraph() {
	/*delete target;
	if (sourceCoarsenedGraph != NULL) {
		for (int i = 0; i < numberOfCoarsenedGraphs; i++) {
			if (sourceCoarsenedGraph[i] != NULL) {
				delete sourceCoarsenedGraph[i];
			}
		}
		free(sourceCoarsenedGraph);
		sourceCoarsenedGraph = NULL;
	}*/
}
