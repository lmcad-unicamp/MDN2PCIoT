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
int computationalWeight[nLayers];

// input
int WI = 16;
int HI = 16;
// TODO: accept DI greater than 1
int DI = 1;

// 1st layer conv filter
int WF = 3;
int HF = 3;
int DF = 1;

// scale parameters for 1st conv
int scale = 4;
int originalDF = 6;
int trainParam = 156;
int compWeight = 51;

// 2nd layer subsampling
int scaleSubsamp = 4 * 6;
int trainParamSubsamp = 12;
int compWeightSubsamp = 6;

// subsequent computational (vertex) weights (FLOP/s)
int compWeightC2 = 48;
int compWeightS2 = 6;
int compWeightFC1 = 51;
int compWeightFC2 = 240;
int compWeightFC3 = 168;

void BuildGraph2DConvolutionalLayer(Graph Source_graph, int scale);

void BuildGraph3DFilterConvolutionalLayer(Graph Source_graph, int scale);

void Build2DSubsamplingLayer(Graph Source_graph, int scale);

void BuildGraph2DConvolutionalLayerC3(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int scale, int originalDF, int trainParam);

void Build2DSubsamplingLayerS4(Graph Source_graph, int conv_starts, int width_input, int height_input, int scale);

void BuildGraph2DConvolutionalLayerC5(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int scale, int originalDF, int trainParam);

void BuildGraph2DFullConnectedLayerFC6(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts, int trainParam, int compWeight);

int WriteInputFile(Graph Source_graph, int num_vertices);

int main(int argc, char *argv[]) {
    Graph Source_graph=NULL;
	int num_vertices;

	computationalWeight[0] = 1;
	computationalWeight[1] = (25 + 25 + 1) * 4;
	computationalWeight[2] = (4 + 1 + 1) * 4;
	computationalWeight[3] = (25 + 25 + 1) * 4;
	computationalWeight[4] = 4 + 1 + 1;
	computationalWeight[5] = 25 + 25 + 1;
	computationalWeight[6] = 120 + 120;
	computationalWeight[7] = 84 + 84;

	// input vertices + vertices with 1st layer conv result + vertices with 2nd layer subsampling result
	num_vertices = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 120 + /*FC6*/ 84 + /*Output*/ 10;

	Source_graph = GraphInit(num_vertices);
	Source_graph->enable_memory = 1;
	Source_graph->enable_vertex_labels = 0;
	Source_graph->enable_edge_weights = 1;
	Source_graph->enable_vertex_weights = 1;

	// conv C1
	if (DF == 1) {
		BuildGraph2DConvolutionalLayer(Source_graph, scale);
	} else {
		BuildGraph3DFilterConvolutionalLayer(Source_graph, scale);
	}

	// subsampling 2D S2
	Build2DSubsamplingLayer(Source_graph, scale); 

	// conv C3
	int conv_starts = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF + ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF;
	int width_input = ((WI - WF) + 1) / 2;
	int height_input = ((HI - HF) + 1) / 2;
	int inputi = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF;
	int originalDF = 16;
	int trainParam = 1516;
	BuildGraph2DConvolutionalLayerC3(Source_graph, conv_starts, width_input, height_input, inputi, scale, originalDF, 0);

	// subsampling 2D S4
	conv_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ (((WI - WF) + 1)/2 - WF + 1) * (((HI - HF) + 1)/2 - WF + 1);
	Build2DSubsamplingLayerS4(Source_graph, conv_starts, width_input, height_input, scale);

	width_input = 5;
	height_input = 5;

	// conv C5
	conv_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ (((WI - WF) + 1)/2 - WF + 1) * (((HI - HF) + 1)/2 - WF + 1) + width_input * height_input;
	width_input = 5;
	height_input = 5;
	inputi = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ (((WI - WF) + 1)/2 - WF + 1) * (((HI - HF) + 1)/2 - WF + 1);
	trainParam = 401;
	scale = 16;
	BuildGraph2DConvolutionalLayerC5(Source_graph, conv_starts, width_input, height_input, inputi, scale, originalDF, trainParam);	

	// FC6
	int output = 84;
	int input = 120;
	int fc_input_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF;
	int fc_output_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 120;
	BuildGraph2DFullConnectedLayerFC6(Source_graph, input, output, fc_input_starts, fc_output_starts, input, compWeightFC2);

	// Output
	input = 84;
	output = 10;
	fc_input_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 120;
	fc_output_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 120 + /*FC6*/ 84;
	BuildGraph2DFullConnectedLayerFC6(Source_graph, input, output, fc_input_starts, fc_output_starts, input - 1, compWeightFC3);

	WriteInputFile(Source_graph, num_vertices);
		
	//PrintGraph(Source_graph);
	return 0;
}

/* __________________________________________________________*/
void BuildGraph2DFullConnectedLayerFC6(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts, int trainParam, int compWeight) {

	for (int i = 0; i < output; i++) {
		for (int j = 0; j < input; j++) {
			GraphInsertArc(Source_graph, fc_input_starts + j, fc_output_starts + i, 1);
			GraphInsertArc(Source_graph, fc_output_starts + i, fc_input_starts + j, 1);
		}
	}

	for (int i = fc_output_starts; i < fc_output_starts + output; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, 1 + trainParam + 1);
		GraphUpdateVertexWeight(Source_graph, i, compWeight);
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayerC5(Graph Source_graph,
									  int conv_starts,
									  int width_input,
									  int height_input,
									  int input_i,
									  int scale,
									  int originalDF,
									  int trainParam) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = conv_starts;
	// conv result width
	conv_width = width_input;
	// conv result height
	conv_height = height_input;
	// conv filter
	int DF = 120;
	int HF = 5;
	int WF = 5;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + DF; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < HF; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < WF; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d, %d %d %d", input_i + (inputi + filteri) * width_input + filterj, convi + convj, inputi, width_input, inputj);
					GraphInsertArc(Source_graph, input_i + (inputi + filteri) * width_input + filterj, convi + convj, scale);
					GraphInsertArc(Source_graph, convi + convj, input_i + (inputi + filteri) * width_input + filterj, scale);
				}	
			}
			inputj++;
		}
		inputj = 0;
		//inputi++;
	}

	for (int i = conv_starts; i < conv_starts + 120; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, 1 + scale * width_input * height_input + 1);
		GraphUpdateVertexWeight(Source_graph, i, compWeightFC1);
	}
}


/* __________________________________________________________*/
void Build2DSubsamplingLayerS4(Graph Source_graph,
							   int conv_starts,
							   int width_input,
							   int height_input,
							   int scale) {
	// conv result index starts after input ends
	int conv_index = conv_starts;
	// conv result width
	int conv_width = width_input - WF + 1;
	// conv result height
	int conv_height = height_input - HF + 1;

	// subsampling result index start after 1st layer conv result ends
	int subsampIndex = conv_index;
	// subsampling result width
	int subsamp_width = conv_width;
	// subsampling result height
	int subsamp_height = conv_height;

	// conv index
	int convi = conv_index - conv_width * conv_height;
	//int convj = convi + conv_width;

	for (int subsampi = subsampIndex; subsampi < subsampIndex + subsamp_height * subsamp_width; subsampi += subsamp_width) {
		//printf("\n subsampi: %d, convi: %d", subsampi, convi);
		for (int subsampj = 0; subsampj < subsamp_width; subsampj++) {
			//printf(" subsampj: %d, sum: %d", subsampj, subsampi + subsampj);
			GraphInsertArc(Source_graph, convi, subsampi + subsampj, 64);
			/*GraphInsertArc(Source_graph, convi + 1, subsampi + subsampj, scaleSubsamp);
			GraphInsertArc(Source_graph, convj, subsampi + subsampj, scaleSubsamp);
			GraphInsertArc(Source_graph, convj + 1, subsampi + subsampj, scaleSubsamp);*/

			GraphInsertArc(Source_graph, subsampi + subsampj, convi, 64);
			/*GraphInsertArc(Source_graph, subsampi + subsampj, convi + 1, scaleSubsamp);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj, scaleSubsamp);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj + 1, scaleSubsamp);*/

			convi += 1;
			//convj = convi + conv_width;
		}
		//convi += conv_width;
		//convj = convi + conv_width;
	}

	for (int i = subsampIndex; i < subsampIndex + conv_width * conv_height; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, 16);
		GraphUpdateVertexWeight(Source_graph, i, 16 * compWeightS2);
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayerC3(Graph Source_graph,
									  int conv_starts,
									  int width_input,
									  int height_input,
									  int input_i,
									  int scale,
									  int originalDF,
									  int trainParam) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = conv_starts;
	// conv result width
	conv_width = width_input - WF + 1;
	// conv result height
	conv_height = height_input - HF + 1;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + conv_width * conv_height; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < HF; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < WF; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d, %d %d %d", input_i + (inputi + filteri) * width_input + filterj + inputj, convi + convj, inputi, width_input, inputj);
					GraphInsertArc(Source_graph, input_i + (inputi + filteri) * width_input + filterj + inputj, convi + convj, scale * 6);
					GraphInsertArc(Source_graph, convi + convj, input_i + (inputi + filteri) * width_input + filterj + inputj, scale * 6);
				}	
			}
			inputj++;
		}
		inputj = 0;
		inputi++;
	}

	for (int i = conv_starts; i < conv_starts + conv_width * conv_height; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale * originalDF);
		GraphUpdateVertexWeight(Source_graph, i, scale * originalDF * compWeightC2);
	}
}

/* __________________________________________________________*/
void Build2DSubsamplingLayer(Graph Source_graph,
							 int scale) {
	// conv result index starts after input ends
	int conv_index = WI * HI * DI;
	// conv result width
	int conv_width = WI - WF + 1;
	// conv result height
	int conv_height = HI - HF + 1;

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
			GraphInsertArc(Source_graph, convi, subsampi + subsampj, scaleSubsamp);
			GraphInsertArc(Source_graph, convi + 1, subsampi + subsampj, scaleSubsamp);
			GraphInsertArc(Source_graph, convj, subsampi + subsampj, scaleSubsamp);
			GraphInsertArc(Source_graph, convj + 1, subsampi + subsampj, scaleSubsamp);

			GraphInsertArc(Source_graph, subsampi + subsampj, convi, scaleSubsamp);
			GraphInsertArc(Source_graph, subsampi + subsampj, convi + 1, scaleSubsamp);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj, scaleSubsamp);
			GraphInsertArc(Source_graph, subsampi + subsampj, convj + 1, scaleSubsamp);

			convi += 2;
			convj = convi + conv_width;
		}
		convi += conv_width;
		convj = convi + conv_width;
	}

	int num_vertices = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF + ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF;
	for (int i = subsampIndex; i < num_vertices; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scaleSubsamp);
		GraphUpdateVertexWeight(Source_graph, i, scaleSubsamp * compWeightSubsamp);
	}
}

/* __________________________________________________________*/
void BuildGraph2DConvolutionalLayer(Graph Source_graph,
									int scale) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = WI * HI * DI;
	// conv result width
	conv_width = WI - WF + 1;
	// conv result height
	conv_height = HI - HF + 1;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + conv_width * conv_height; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < HF; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < WF; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d", (inputi + filteri) * WI + filterj + inputj, convi + convj);
					GraphInsertArc(Source_graph, (inputi + filteri) * WI + filterj + inputj, convi + convj, scale);
					GraphInsertArc(Source_graph, convi + convj, (inputi + filteri) * WI + filterj + inputj, scale);
				}	
			}
			inputj++;
		}
		inputj = 0;
		inputi++;
	}

	for (int i = 0; i < conv_index; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale);
		GraphUpdateVertexWeight(Source_graph, i, 0);
	}
	int num_vertices = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF;
	for (int i = conv_index; i < num_vertices; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale * originalDF);
		GraphUpdateVertexWeight(Source_graph, i, scale * originalDF * compWeight);
	}
}

/* __________________________________________________________*/
void BuildGraph3DFilterConvolutionalLayer(Graph Source_graph,
										  int scale) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height, convk;

	// conv result index starts after input ends
	conv_index = WI * HI * DI;
	// conv result width
	conv_width = WI - WF + 1;
	// conv result height
	conv_height = HI - HF + 1;
	// input index
	inputi = 0;
	inputj = 0;

	for (convi = conv_index; convi < conv_index + conv_width * conv_height; convi+=conv_width) {
		//printf("\n convi=%d", convi); 
		for (convj = 0; convj < conv_width; convj++) {
			//printf("\n convj=%d", convj);
			for (filteri = 0; filteri < HF; filteri++) {
				//printf("\n filteri=%d", filteri);
				for (filterj = 0; filterj < WF; filterj++) { 
					//printf("\n filterj=%d", filterj);
					//printf("\n%d %d", (inputi + filteri) * WI + filterj + inputj, convi + convj);
					for (convk = 0; convk < DF * conv_width * conv_height; convk += conv_width * conv_height) {
						GraphInsertArc(Source_graph, (inputi + filteri) * WI + filterj + inputj, convi + convj + convk, scale);
						GraphInsertArc(Source_graph, convi + convj + convk, (inputi + filteri) * WI + filterj + inputj, scale);
					}
				}	
			}
			inputj++;
		}
		inputj = 0;
		inputi++;
	}

	for (int i = 0; i < conv_index; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale);
		GraphUpdateVertexWeight(Source_graph, i, 0);
	}
	int num_vertices = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF;
	for (int i = conv_index; i < num_vertices; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale * originalDF + trainParam);
		GraphUpdateVertexWeight(Source_graph, i, scale * compWeight);
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
	fprintf(wif, "0\n%d %d\n0  %d%d%d%d", Source_graph->V, Source_graph->A, Source_graph->enable_memory, Source_graph->enable_vertex_labels, Source_graph->enable_vertex_weights, Source_graph->enable_edge_weights);

	// write information about the size of layers and amount of shared memory for each layer (layer parameters)
	fprintf(wif, "\n%d\n", 8);
	//for (int layer = 0; layer < num_layers; layer++) {
		fprintf(wif, "0 256 452 501 526 551 671 755");
	//}
	fprintf(wif, "\n");
	//for (int layer = 0; layer < num_layers; layer++) {
		fprintf(wif, "0 156 12 1516 32 0 0 0");
	//}

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
