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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <stdbool.h>
#include "sourceGraph.h"

// number of layers
#define nLayers 8
int initialFullyConnectedLayer = 5;
int wi[nLayers];
int hi[nLayers];
int di[nLayers];
int wiScale[nLayers];
int hiScale[nLayers];
int diScale[nLayers];
int eachEdge[nLayers];
int mem[nLayers];
int shared[nLayers];
int wf[nLayers];
int hf[nLayers];
int df[nLayers];
int initialPos[nLayers];
int num_vertices = 0;
int computationalWeight[nLayers];

void BuildGraph2DConvolutionalLayer(Graph Source_graph, int layer);

//void BuildGraph3DFilterConvolutionalLayer(Graph Source_graph, int scale);

void Build2DSubsamplingLayer(Graph Source_graph, int layer);

void BuildGraph2DConvolutionalLayerC3(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int layer);

void Build2DSubsamplingLayerS4(Graph Source_graph, int conv_starts, int width_input, int height_input, int layer);

void BuildGraph2DConvolutionalLayerC5(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int layer);

void BuildGraph2DFullConnectedLayerFC6(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts);

int WriteInputFile(Graph Source_graph, int num_vertices);

int main(int argc, char *argv[]) {
	Graph Source_graph=NULL;

	for (int i = 0; i < nLayers; i++) {
		wiScale[i] = 1;
		hiScale[i] = 1;
		diScale[i] = 1;
	}

	// 1st layer: input
	wi[0] = 32;
	hi[0] = 32;
	// TODO: accept DI greater than 1
	di[0] = 1;
	//wiScale[0] = 16;
	//hiScale[0] = 16;
	eachEdge[0] = wiScale[0] * hiScale[0] * diScale[0];
	mem[0] = wiScale[0] * hiScale[0] * diScale[0];
	shared[0] = 0;
	computationalWeight[0] = 0;

	// 2nd layer: conv
	wi[1] = 28;
	hi[1] = 28;
	di[1] = 6;
	//wiScale[1] = 14;
	//hiScale[1] = 14;
	diScale[1] = 6;
	eachEdge[1] = wiScale[1] * hiScale[1] * diScale[1];
	mem[1] = wiScale[1] * hiScale[1] * diScale[1];
	// shared parameters
	shared[1] = 156;
	// filter size
	wf[1] = 5;
	hf[1] = 5;
	df[1] = 6;
	computationalWeight[1] = (25 + 25 + 1) * 6;

	// 3th layer: subsampling
	wi[2] = 14;
	hi[2] = 14;
	di[2] = 6;
	//wiScale[2] = 7;
	//hiScale[2] = 7;
	diScale[2] = 6;
	eachEdge[2] = wiScale[2] * hiScale[2] * diScale[2];
	mem[2] = wiScale[2] * hiScale[2] * diScale[2];
	// shared parameters
	shared[2] = 12;
	// filter size
	wf[2] = 2;
	hf[2] = 2;
	df[2] = 6;
	computationalWeight[2] = (4 + 1 + 1) * 6;

	// 4th layer: conv
	wi[3] = 10;
	hi[3] = 10;
	di[3] = 16;
	//wiScale[3] = 5;
	//hiScale[3] = 5;
	diScale[3] = 16;
	eachEdge[3] = wiScale[3] * hiScale[3] * diScale[3];
	mem[3] = wiScale[3] * hiScale[3] * diScale[3];
	// shared parameters
	shared[3] = 1516;
	// filter size
	wf[3] = 5;
	hf[3] = 5;
	df[3] = 16;
	computationalWeight[3] = (25 + 25 + 1) * 15;

	// 5th layer: subsampling
	wi[4] = 5;
	hi[4] = 5;
	di[4] = 16;
	//wiScale[4] = 5;
	//hiScale[4] = 5;
	diScale[4] = 16;
	eachEdge[4] = wiScale[4] * hiScale[4] * diScale[4];
	mem[4] = wiScale[4] * hiScale[4] * diScale[4];
	// shared parameters
	shared[4] = 32;
	// filter size
	wf[4] = 2;
	hf[4] = 2;
	df[4] = 16;
	computationalWeight[4] = (4 + 1 + 1) * 16;

	// 6th layer: fully connected
	wi[5] = 1;
	hi[5] = 1;
	di[5] = 120;
	//diScale[5] = 24;
	eachEdge[5] = wiScale[5] * hiScale[5] * diScale[5];
	mem[5] = 402; // 1 + 401 * 1
	// shared parameters
	shared[5] = 0;
	wf[5] = 5;
	hf[5] = 5;
	df[5] = 16;
	computationalWeight[5] = 25 + 25 + 1;

	// 7th layer: fully connected
	wi[6] = 1;
	hi[6] = 1;
	di[6] = 84;
	//diScale[6] = 14;
	eachEdge[6] = 1;
	mem[6] = 122; // 1 + 121 * 1
	// shared parameters
	shared[6] = 0;
	computationalWeight[6] = 120 + 120;

	// 8th layer: fully connected
	wi[7] = 1;
	hi[7] = 1;
	di[7] = 10;
	//diScale[7] = 5;
	eachEdge[7] = 1;
	mem[7] = 86; // 1 + 85 * 1
	// shared parameters
	shared[7] = 0;
	computationalWeight[7] = 84 + 84;

	initialPos[0] = 0;
	for (int i = 1; i < nLayers; i++) {
		initialPos[i] += initialPos[i - 1] + (wi[i - 1] / wiScale[i - 1]) * (hi[i - 1] / hiScale[i - 1]) * (di[i - 1] / diScale[i - 1]); 
	}

	for (int i = 0; i < nLayers; i++) {
		num_vertices += (wi[i] / wiScale[i]) * (hi[i] / hiScale[i]) * (di[i] / diScale[i]);
	}

	Source_graph = GraphInit(num_vertices);
	Source_graph->enable_memory = 1;
	Source_graph->enable_vertex_labels = 0;
	Source_graph->enable_edge_weights = 1;
	Source_graph->enable_vertex_weights = 1;

	// conv C1
	//if (diScale[1] > 1) {
		BuildGraph2DConvolutionalLayer(Source_graph, 1);
	/*} else {
		BuildGraph3DFilterConvolutionalLayer(Source_graph, 1, 2, 2, 2, 2, 6);
	}*/

	// subsampling 2D S2
	Build2DSubsamplingLayer(Source_graph, 2); 

	// conv C3
	BuildGraph2DConvolutionalLayerC3(Source_graph, initialPos[3], wi[2], hi[2], initialPos[2], 3);

	// subsampling 2D S4
	Build2DSubsamplingLayerS4(Source_graph, initialPos[3], wi[3], hi[3], 4);

	// conv C5
	BuildGraph2DConvolutionalLayerC5(Source_graph, initialPos[5], wi[4], hi[4], initialPos[4], 5);	

	// FC6
	BuildGraph2DFullConnectedLayerFC6(Source_graph, di[5], di[6], initialPos[5], initialPos[6]);

	// Output
	BuildGraph2DFullConnectedLayerFC6(Source_graph, di[6], di[7], initialPos[6], initialPos[7]);

	// set each vertex memory
	for (int layer = 0; layer < nLayers; layer++) {
		int last;
		if (layer < nLayers - 1)
			last = initialPos[layer + 1];
		else
			last = num_vertices;
		for (int i = initialPos[layer]; i < last; i++) {
			GraphUpdateMemoryWeight(Source_graph, i, mem[layer]);
		}
	}

	WriteInputFile(Source_graph, num_vertices);
		
	PrintGraph(Source_graph);
	return 0;
}

/* __________________________________________________________*/
void BuildGraph2DFullConnectedLayerFC6(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts) {

	for (int i = 0; i < output; i++) {
		for (int j = 0; j < input; j++) {
			GraphInsertArc(Source_graph, fc_input_starts + j, fc_output_starts + i, 1);
			GraphInsertArc(Source_graph, fc_output_starts + i, fc_input_starts + j, 1);
		}
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayerC5(Graph Source_graph,
					  int conv_starts,
					  int width_input,
					  int height_input,
					  int input_i,
					  int layer) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = conv_starts;
	// conv result width
	conv_width = width_input;
	// conv result height
	conv_height = height_input;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + di[layer]; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < hf[layer]; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < wf[layer]; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d, %d %d %d", input_i + (inputi + filteri) * width_input + filterj, convi + convj, inputi, width_input, inputj);
					GraphInsertArc(Source_graph, input_i + (inputi + filteri) * width_input + filterj, convi + convj, mem[layer - 1]);
					GraphInsertArc(Source_graph, convi + convj, input_i + (inputi + filteri) * width_input + filterj, mem[layer - 1]);
				}	
			}
			inputj++;
		}
		inputj = 0;
		//inputi++;
	}
}

/* __________________________________________________________*/
void Build2DSubsamplingLayerS4(Graph Source_graph,
				   int conv_starts,
				   int width_input,
				   int height_input,
				   int layer) {
	// conv result index starts after input ends
	int conv_index = conv_starts;
	// conv result width
	int conv_width = wi[layer - 2] - wf[layer - 1] + 1;
	// conv result height
	int conv_height = hi[layer - 2] - hf[layer - 1] + 1;

	// subsampling result index start after 1st layer conv result ends
	int subsampIndex = initialPos[layer];
	// subsampling result width
	int subsamp_width = wi[layer];
	// subsampling result height
	int subsamp_height = hi[layer];

	// conv index
	int convi = conv_index;
	int convj = convi + conv_width;

	for (int subsampi = subsampIndex; subsampi < subsampIndex + subsamp_height * subsamp_width; subsampi += subsamp_width) {
		//printf("\n subsampi: %d, convi: %d", subsampi, convi);
		for (int subsampj = 0; subsampj < subsamp_width; subsampj++) {
			//printf(" subsampj: %d, sum: %d", subsampj, subsampi + subsampj);
			GraphInsertArc(Source_graph, convi, subsampi + subsampj, mem[layer - 1]);
			GraphInsertArc(Source_graph, convi + 1, subsampi + subsampj, mem[layer - 1]);
			GraphInsertArc(Source_graph, convj, subsampi + subsampj, mem[layer - 1]);
			GraphInsertArc(Source_graph, convj + 1, subsampi + subsampj, mem[layer - 1]);

			GraphInsertArc(Source_graph, subsampi + subsampj, convi, mem[layer - 1]);
			GraphInsertArc(Source_graph, subsampi + subsampj, convi + 1, mem[layer - 1]);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj, mem[layer - 1]);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj + 1, mem[layer - 1]);

			convi += 2;
			convj = convi + conv_width;
		}
		convi += conv_width;
		convj = convi + conv_width;
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayerC3(Graph Source_graph,
					  int conv_starts,
					  int width_input,
					  int height_input,
					  int input_i,
					  int layer) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = conv_starts;
	// conv result width
	conv_width = width_input - wf[layer] + 1;
	// conv result height
	conv_height = height_input - hf[layer] + 1;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + conv_width * conv_height; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < hf[layer]; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < wf[layer]; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d, %d %d %d", input_i + (inputi + filteri) * width_input + filterj + inputj, convi + convj, inputi, width_input, inputj);
					GraphInsertArc(Source_graph, input_i + (inputi + filteri) * width_input + filterj + inputj, convi + convj, mem[layer - 1]);
					GraphInsertArc(Source_graph, convi + convj, input_i + (inputi + filteri) * width_input + filterj + inputj, mem[layer - 1]);
				}	
			}
			inputj++;
		}
		inputj = 0;
		inputi++;
	}
}

/* __________________________________________________________*/
void Build2DSubsamplingLayer(Graph Source_graph,
				int layer) {
	// conv result index starts after input ends
	int conv_index = initialPos[layer - 1];
	// conv result width
	int conv_width = wi[layer - 2] - wf[layer - 1] + 1;
	// conv result height
	int conv_height = hi[layer - 2] - hf[layer - 1] + 1;

	// subsampling result index start after 1st layer conv result ends
	int subsampIndex = conv_index + conv_width * conv_height;
	// subsampling result width
	int subsamp_width = conv_width / 2;
	// subsampling result height
	int subsamp_height = conv_height / 2;

	// conv index
	int convi = conv_index;
	int convj = convi + conv_width;

	for (int subsampi = subsampIndex; subsampi < subsampIndex + subsamp_height * subsamp_width; subsampi += subsamp_width) {
		for (int subsampj = 0; subsampj < subsamp_width; subsampj++) {
			GraphInsertArc(Source_graph, convi, subsampi + subsampj, mem[layer -1]);
			GraphInsertArc(Source_graph, convi + 1, subsampi + subsampj, mem[layer -1]);
			GraphInsertArc(Source_graph, convj, subsampi + subsampj, mem[layer -1]);
			GraphInsertArc(Source_graph, convj + 1, subsampi + subsampj, mem[layer -1]);

			GraphInsertArc(Source_graph, subsampi + subsampj, convi, mem[layer -1]);
			GraphInsertArc(Source_graph, subsampi + subsampj, convi + 1, mem[layer -1]);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj, mem[layer -1]);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj + 1, mem[layer -1]);

			convi += 2;
			convj = convi + conv_width;
		}
		convi += conv_width;
		convj = convi + conv_width;
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayer(Graph Source_graph,
					int layer) {
	int filteri, filterj;
	int convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = initialPos[layer];
	// conv result width
	conv_width = wi[layer - 1] - wf[layer] + 1;
	// conv result height
	conv_height = hi[layer - 1] - hf[layer] + 1;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + conv_width * conv_height; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < hf[layer]; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < wf[layer]; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d", (inputi + filteri) * wi[layer - 1] + filterj + inputj, convi + convj);
					GraphInsertArc(Source_graph, (inputi + filteri) * wi[layer - 1] + filterj + inputj, convi + convj, eachEdge[layer - 1]);
					GraphInsertArc(Source_graph, convi + convj, (inputi + filteri) * wi[layer - 1] + filterj + inputj, eachEdge[layer - 1]);
				}	
			}
			inputj++;
		}
		inputj = 0;
		inputi++;
	}
}

/* __________________________________________________________*/
int WriteInputFile(Graph Source_graph,
				   int num_vertices) {
	FILE *wif;
	char filename[50];
	int vertice;
	link adj;

	// file name according to input parameters
	strcpy(filename, "LeNet-");
	filename[6] = (char) (num_vertices / 1000) + '0';
	filename[7] = (char) (num_vertices % 1000 / 100) + '0';
	filename[8] = (char) (num_vertices % 100 / 10) + '0';
	filename[9] = (char) (num_vertices % 10) + '0';
	strcat(filename, "vertices");
	strcat(filename, ".grf");
	filename[23] = '\0';

	wif = fopen(filename, "w");
	if (wif == NULL) {
		printf("\nCould not create output file!");
		return 1;
	}

	// write header information (number of vertices and flags)
	fprintf(wif, "0\n%d %d\n0  %d%d%d%d", Source_graph->V, Source_graph->A/2, Source_graph->enable_memory, Source_graph->enable_vertex_labels, Source_graph->enable_vertex_weights, Source_graph->enable_edge_weights);

	// write information about the size of layers and amount of shared memory for each layer (layer parameters)
	fprintf(wif, "\n%d\n", 8);
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", initialPos[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", shared[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", wi[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", hi[layer]);
	}
	fprintf(wif, "\n");
	/*for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", df[layer]);
	}*/
	fprintf(wif, "1 6 6 16 16 120 84 10");

	for (vertice = 0; vertice < Source_graph->V; vertice++) {
		// write amount of memory needed for the vertex, amount of computation and vertex degree (vertex degree not ready yet)
		// TODO: write correct vertex degree (writing 1 for now)
		int layer = 0;
		for (int j = 0; j < nLayers; j++) {
			if (layer == nLayers - 1)
				break;
			if (vertice < initialPos[j + 1]) {
				break;
			} else {
				layer++;
			}
		}
		fprintf(wif, "\n%d %d 1", Source_graph->memory[vertice], computationalWeight[layer]);
		// write list of adjacency of each vertex
		adj = Source_graph->adj[vertice]; 
		for (adj; adj != NULL; adj = adj->next) {
			fprintf(wif, " %d %d", adj->edge_weight, adj->w);
		}
	}

	fclose(wif);
	return 0;	
}
