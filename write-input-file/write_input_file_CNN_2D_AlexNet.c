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

// AlexNet //

// number of layers
#define nLayers 12
int initialFullyConnectedLayer = 9;
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
int stride[nLayers];
int padding[nLayers];

void BuildGraph2DConvolutionalLayer(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int layer);

void Build2DSubsamplingLayer(Graph Source_graph, int conv_starts, int layer);

void BuildGraph2DFullyConnectedLayer2Dinput(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int layer);

void BuildGraph2DFullyConnectedLayer(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts);

int WriteInputFile(Graph Source_graph, int num_vertices);

int main(int argc, char *argv[]) {
	Graph Source_graph=NULL;

	for (int i = 0; i < nLayers; i++) {
		wiScale[i] = 1;
		hiScale[i] = 1;
		diScale[i] = 1;
		stride[i] = 1;
		padding[i] = 0;
	}

	// layer 0: input
	wi[0] = 227;
	hi[0] = 227;
	// TODO: accept DI greater than 1 (depth greater than 1)
	di[0] = 3;
	//wiScale[0] = 16;
	//hiScale[0] = 16;
	diScale[0] = 3;
	eachEdge[0] = wiScale[0] * hiScale[0] * diScale[0];
	mem[0] = wiScale[0] * hiScale[0] * diScale[0];
	shared[0] = 0;
	computationalWeight[0] = 0;

	// layer 1: conv
	wi[1] = 55;
	hi[1] = 55;
	di[1] = 96;
	//wiScale[1] = 14;
	//hiScale[1] = 14;
	diScale[1] = 96;
	eachEdge[1] = wiScale[1] * hiScale[1] * diScale[1];
	mem[1] = wiScale[1] * hiScale[1] * diScale[1];
	// shared parameters
	shared[1] = 96 * (11 * 11 * 3 + 1);
	// filter size
	wf[1] = 11;
	hf[1] = 11;
	df[1] = 3;
	computationalWeight[1] = (121 * 3 * 2 + 1 + 1) * 96 + 13;
	stride[1] = 4;

	// layer 2: subsampling
	wi[2] = 27;
	hi[2] = 27;
	di[2] = 96;
	//wiScale[2] = 7;
	//hiScale[2] = 7;
	diScale[2] = 96;
	eachEdge[2] = wiScale[2] * hiScale[2] * diScale[2];
	mem[2] = wiScale[2] * hiScale[2] * diScale[2];
	// shared parameters
	shared[2] = 0;
	// filter size
	wf[2] = 3;
	hf[2] = 3;
	df[2] = 96;
	computationalWeight[2] = 9 * 96;
	stride[2] = 2;

	// layer 3: conv
	wi[3] = 27;
	hi[3] = 27;
	di[3] = 256;
	//wiScale[3] = 5;
	//hiScale[3] = 5;
	diScale[3] = 256;
	eachEdge[3] = wiScale[3] * hiScale[3] * diScale[3];
	mem[3] = wiScale[3] * hiScale[3] * diScale[3];
	// shared parameters
	shared[3] = 256 * (5 * 5 * 96 + 1);
	// filter size
	wf[3] = 5;
	hf[3] = 5;
	df[3] = 96;
	computationalWeight[3] = (25 + 26 + 1) * 96;
	padding[3] = 4;

	// layer 4: subsampling
	wi[4] = 13;
	hi[4] = 13;
	di[4] = 256;
	//wiScale[4] = 5;
	//hiScale[4] = 5;
	diScale[4] = 256;
	eachEdge[4] = wiScale[4] * hiScale[4] * diScale[4];
	mem[4] = wiScale[4] * hiScale[4] * diScale[4];
	// shared parameters
	shared[4] = 0;
	// filter size
	wf[4] = 3;
	hf[4] = 3;
	df[4] = 256;
	computationalWeight[4] = 9 * 256;
	stride[4] = 2;

	// layer 5: conv
	wi[5] = 13;
	hi[5] = 13;
	di[5] = 384;
	diScale[5] = 384;
	eachEdge[5] = wiScale[5] * hiScale[5] * diScale[5];
	mem[5] = wiScale[5] * hiScale[5] * diScale[5];
	// shared parameters
	shared[5] = (3 * 3 * 256 + 1) * 384;
	wf[5] = 3;
	hf[5] = 3;
	df[5] = 256;
	computationalWeight[5] = 256 * (9 + 10 + 1);
	padding[5] = 2;

	// layer 6: conv
	wi[6] = 13;
	hi[6] = 13;
	di[6] = 384;
	diScale[6] = 384;
	eachEdge[6] = wiScale[6] * hiScale[6] * diScale[6];
	mem[6] = wiScale[6] * hiScale[6] * diScale[6];
	// shared parameters
	shared[6] = (3 * 3 * 384 + 1) * 384;
	wf[6] = 3;
	hf[6] = 3;
	df[6] = 384;
	computationalWeight[6] = 384 * (9 + 10 + 1);
	padding[6] = 2;

	// layer 7: conv
	wi[7] = 13;
	hi[7] = 13;
	di[7] = 256;
	diScale[7] = 256;
	eachEdge[7] = wiScale[7] * hiScale[7] * diScale[7];
	mem[7] = wiScale[7] * hiScale[7] * diScale[7];
	// shared parameters
	shared[7] = (3 * 3 * 384 + 1) * 256;
	wf[7] = 3;
	hf[7] = 3;
	df[7] = 384;
	computationalWeight[7] = 384 * (9 + 10 + 1);
	padding[7] = 2;

	// layer 8: pooling
	wi[8] = 6;
	hi[8] = 6;
	di[8] = 256;
	diScale[8] = 256;
	eachEdge[8] = wiScale[8] * hiScale[8] * diScale[8];
	mem[8] = wiScale[8] * hiScale[8] * diScale[8];
	// shared parameters
	shared[8] = 0;
	wf[8] = 3;
	hf[8] = 3;
	df[8] = 256;
	computationalWeight[8] = 256 * 9;
	stride[8] = 2;

	// layer 9: fully connected
	wi[9] = 1;
	hi[9] = 1;
	di[9] = 4096;
	diScale[9] = 1;
	eachEdge[9] = wiScale[9] * hiScale[9] * diScale[9];
	mem[9] = 9216 + 1 + 1;
	// shared parameters
	shared[9] = 0;
	wf[9] = 6;
	hf[9] = 6;
	df[9] = 256;
	computationalWeight[9] = 9216 + 9217 + 1;

	// layer 10: fully connected
	wi[10] = 1;
	hi[10] = 1;
	di[10] = 4096;
	diScale[10] = 1;
	eachEdge[10] = wiScale[10] * hiScale[10] * diScale[10];
	mem[10] = 4096 + 1 + 1;
	// shared parameters
	shared[10] = 0;
	computationalWeight[10] = 4096 + 4097 + 1;

	// layer 11: softmax
	wi[11] = 1;
	hi[11] = 1;
	di[11] = 1000;
	diScale[11] = 1;
	eachEdge[11] = wiScale[11] * hiScale[11] * diScale[11];
	mem[11] = 4096 + 1 + 1;
	// shared parameters
	shared[11] = 0;
	computationalWeight[11] = 1 + 2 + 2 * (4096 - 1);

	initialPos[0] = 0;
	for (int i = 1; i < nLayers; i++) {
		initialPos[i] += initialPos[i - 1] + (wi[i - 1] / wiScale[i - 1]) * (hi[i - 1] / hiScale[i - 1]) * (di[i - 1] / diScale[i - 1]); 
		//printf("\ninitialPos[%d]: %d", i, initialPos[i]);
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
	BuildGraph2DConvolutionalLayer(Source_graph, initialPos[1], wi[0], hi[0], initialPos[0], 1);

	// subsampling S2
	Build2DSubsamplingLayer(Source_graph, initialPos[1], 2); 

	// conv C3
	BuildGraph2DConvolutionalLayer(Source_graph, initialPos[3], wi[2], hi[2], initialPos[2], 3);

	// subsampling S4
	Build2DSubsamplingLayer(Source_graph, initialPos[3], 4);

	// conv C5
	BuildGraph2DConvolutionalLayer(Source_graph, initialPos[5], wi[4], hi[4], initialPos[4], 5);	

	// conv C6
	BuildGraph2DConvolutionalLayer(Source_graph, initialPos[6], wi[5], hi[5], initialPos[5], 6);

	// conv C7
	BuildGraph2DConvolutionalLayer(Source_graph, initialPos[7], wi[6], hi[6], initialPos[6], 7);

	// pool P8
	Build2DSubsamplingLayer(Source_graph, initialPos[7], 8);

	// FC9
	BuildGraph2DFullyConnectedLayer2Dinput(Source_graph, initialPos[9], wi[8], hi[8], initialPos[8], 9);

	// FC10
	BuildGraph2DFullyConnectedLayer(Source_graph, di[9], di[10], initialPos[9], initialPos[10]);

	// Output (Softmax 11)
	BuildGraph2DFullyConnectedLayer(Source_graph, di[10], di[11], initialPos[10], initialPos[11]);

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
		
	//PrintGraph(Source_graph);
	return 0;
}

/* __________________________________________________________*/
void BuildGraph2DFullyConnectedLayer(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts) {

	for (int i = 0; i < output; i++) {
		for (int j = 0; j < input; j++) {
			GraphInsertArc(Source_graph, fc_input_starts + j, fc_output_starts + i, 1);
			GraphInsertArc(Source_graph, fc_output_starts + i, fc_input_starts + j, 1);
		}
	}
}

/* __________________________________________________________*/
void BuildGraph2DFullyConnectedLayer2Dinput(Graph Source_graph,
					  int fc_starts,
					  int width_input,
					  int height_input,
					  int input_i,
					  int layer) {
	int filteri, filterj, fcd, fc_index, fc_width, fc_height;

	// fc result index starts after input ends
	fc_index = fc_starts;
	// fc result width
	fc_width = width_input;
	// fc result height
	fc_height = height_input;

	for (fcd = fc_index; fcd < fc_index + di[layer]; fcd++) {
		//printf("\n fcd=%d", fcd); 
		for (filteri = 0; filteri < hf[layer]; filteri++) {
			//printf("\n filteri=%d", filteri);
			for (filterj = 0; filterj < wf[layer]; filterj++) { 
				//printf("\n filterj=%d", filterj);
				GraphInsertArc(Source_graph, input_i + (filteri) * width_input + filterj, fcd, eachEdge[layer - 1]);
				GraphInsertArc(Source_graph, fcd, input_i + (filteri) * width_input + filterj, eachEdge[layer - 1]);
			}	
		}
	}
}

/* __________________________________________________________*/
void Build2DSubsamplingLayer(Graph Source_graph,
				   int conv_starts,
				   int layer) {
	// conv result index starts after input ends
	int conv_index = conv_starts;
	// conv result width
	int conv_width = wi[layer - 1];
	// conv result height
	int conv_height = hi[layer - 1];

	// subsampling result index start after 1st layer conv result ends
	int subsampIndex = initialPos[layer];
	// subsampling result width
	int subsamp_width = wi[layer];
	// subsampling result height
	int subsamp_height = hi[layer];

	// conv index
	int convi = conv_index;
	//int convj = convi + conv_width;

	for (int subsampi = subsampIndex; subsampi < subsampIndex + subsamp_height * subsamp_width; subsampi += subsamp_width) {
		//printf("\n subsampi: %d, convi: %d", subsampi, convi);
		for (int subsampj = 0; subsampj < subsamp_width; subsampj++) {
			//printf(" subsampj: %d, sum: %d", subsampj, subsampi + subsampj);
			for (int filteri = 0; filteri < hf[layer]; filteri++) {
				for (int filterj = 0; filterj < wf[layer]; filterj++) {
					GraphInsertArc(Source_graph, convi + filterj + filteri * conv_width, subsampi + subsampj, eachEdge[layer - 1]);
	
					GraphInsertArc(Source_graph, subsampi + subsampj, convi + filterj + filteri * conv_width, eachEdge[layer - 1]);
				}
			}
			convi += stride[layer];
			//convj = convi + conv_width;
		}
		convi += wf[layer] - stride[layer] + conv_width;
		//convj = convi + conv_width;
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayer(Graph Source_graph,
					  int conv_starts,
					  int width_input,
					  int height_input,
					  int input_i,
					  int layer) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = conv_starts;
	// conv result width
	conv_width = (width_input - wf[layer] + padding[layer]) / stride[layer] + 1;
	// conv result height
	conv_height = (height_input - hf[layer] + padding[layer]) / stride[layer] + 1;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + conv_width * conv_height; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < hf[layer]; filteri++) {
				//printf("\n filteri=%d", filteri);
				if (padding[layer] != 0 && ((convi - conv_index)/conv_width < hf[layer] - padding[layer] / 2 - 1 || convi == conv_index + conv_width * conv_height - conv_width) && (filteri >= hf[layer] - padding[layer] / 2 + (convi - conv_index)/conv_width || filteri >= hf[layer] - padding[layer] / 2 + conv_height - 1 - (convi - conv_index)/conv_width))
					continue;
				for (filterj = 0; filterj < wf[layer]; filterj++) { 
					//printf("\n filterj=%d", filterj);
					if (padding[layer] != 0 && (convj < wf[layer] - padding[layer] / 2 - 1 || convj > conv_width - 1 - wf[layer] + padding[layer] / 2 + 1) && (filterj >= wf[layer] - padding[layer] / 2 + convj || filterj >= conv_width - 1 + wf[layer] - padding[layer] / 2 - convj))
						continue;
					//printf("\n%d %d, %d %d %d", input_i + (inputi + filteri) * width_input + filterj + inputj, convi + convj, inputi, width_input, inputj);
					GraphInsertArc(Source_graph, input_i + (inputi + filteri) * width_input + filterj + inputj, convi + convj, eachEdge[layer - 1]);
					GraphInsertArc(Source_graph, convi + convj, input_i + (inputi + filteri) * width_input + filterj + inputj, eachEdge[layer - 1]);
				}	
			}
			if (convj >= padding[layer] / 2) 
				inputj += stride[layer];
		}
		inputj = 0;
		if (convi >= conv_index + conv_width * padding[layer] / 2) 
			inputi += stride[layer];
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
	strcpy(filename, "AlexNet-");
	filename[8] = (char) (num_vertices / 10000) + '0';
	filename[9] = (char) (num_vertices % 10000 / 1000) + '0';
	filename[10] = (char) (num_vertices % 1000 / 100) + '0';
	filename[11] = (char) (num_vertices % 100 / 10) + '0';
	filename[12] = (char) (num_vertices % 10) + '0';
	strcat(filename, "vertices");
	strcat(filename, ".grf");
	filename[25] = '\0';

	wif = fopen(filename, "w");
	if (wif == NULL) {
		printf("\nCould not create output file!");
		return 1;
	}

	// write header information (number of vertices and flags)
	fprintf(wif, "0\n%d %d\n0  %d%d%d%d", Source_graph->V, Source_graph->A/2, Source_graph->enable_memory, Source_graph->enable_vertex_labels, Source_graph->enable_vertex_weights, Source_graph->enable_edge_weights);

	// write information about the size of layers and amount of shared memory for each layer (layer parameters)
	fprintf(wif, "\n%d\n", nLayers);
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
	/*fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", df[layer]);
	}*/
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", di[layer]);
	}

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
