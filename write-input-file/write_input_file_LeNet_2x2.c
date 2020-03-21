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

void BuildGraph2DConvolutionalLayer(Graph Source_graph, int scale);

void BuildGraph3DFilterConvolutionalLayer(Graph Source_graph, int scale);

void Build2DSubsamplingLayer(Graph Source_graph, int scale);

void BuildGraph2DConvolutionalLayerC3(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int scale, int originalDF, int trainParam);

void Build2DSubsamplingLayerS4(Graph Source_graph, int conv_starts, int width_input, int height_input, int scale);

void BuildGraph2DConvolutionalLayerC5(Graph Source_graph, int conv_starts, int width_input, int height_input, int input_i, int scale, int originalDF, int trainParam, int layer);

void BuildGraph2DFullConnectedLayerFC6(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts, int trainParam, int compWeight, int layer);

int WriteInputFile(Graph Source_graph, int num_vertices);

int main(int argc, char *argv[]) {
    Graph Source_graph=NULL;
	int num_vertices;

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
	wiScale[0] = 2;
	hiScale[0] = 2;
	eachEdge[0] = wiScale[0] * hiScale[0] * diScale[0];
	mem[0] = wiScale[0] * hiScale[0] * diScale[0];
	shared[0] = 0;
	computationalWeight[0] = 4;

	// 2nd layer: conv
	wi[1] = 28;
	hi[1] = 28;
	di[1] = 6;
	wiScale[1] = 2;
	hiScale[1] = 2;
	diScale[1] = 6;
	eachEdge[1] = wiScale[1] * hiScale[1] * diScale[1];
	mem[1] = wiScale[1] * hiScale[1] * diScale[1];
	// shared parameters
	shared[1] = 156;
	// filter size
	wf[1] = 5;
	hf[1] = 5;
	df[1] = 6;
	computationalWeight[1] = (25 + 25 + 1) * 24;

	// 3th layer: subsampling
	wi[2] = 14;
	hi[2] = 14;
	di[2] = 6;
	wiScale[2] = 2;
	hiScale[2] = 2;
	diScale[2] = 6;
	eachEdge[2] = wiScale[2] * hiScale[2] * diScale[2];
	mem[2] = wiScale[2] * hiScale[2] * diScale[2];
	// shared parameters
	shared[2] = 12;
	// filter size
	wf[2] = 2;
	hf[2] = 2;
	df[2] = 6;
	computationalWeight[2] = (4 + 1 + 1) * 24;

	// 4th layer: conv
	wi[3] = 10;
	hi[3] = 10;
	di[3] = 16;
	wiScale[3] = 2;
	hiScale[3] = 2;
	diScale[3] = 16;
	eachEdge[3] = wiScale[3] * hiScale[3] * diScale[3];
	mem[3] = wiScale[3] * hiScale[3] * diScale[3];
	// shared parameters
	shared[3] = 1516;
	// filter size
	wf[3] = 5;
	hf[3] = 5;
	df[3] = 16;
	computationalWeight[3] = (25 + 25 + 1) * 15 * 4;

	// 5th layer: subsampling
	wi[4] = 5;
	hi[4] = 5;
	di[4] = 16;
	wiScale[4] = 1;
	hiScale[4] = 1;
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
	diScale[5] = 4;
	eachEdge[5] = wiScale[5] * hiScale[5] * diScale[5];
	mem[5] = 1608; // 4 + 401 * 4
	// shared parameters
	shared[5] = 0;
	wf[5] = 5;
	hf[5] = 5;
	df[5] = 16;
	computationalWeight[5] = (25 + 25 + 1) * 4;

	// 7th layer: fully connected
	wi[6] = 1;
	hi[6] = 1;
	di[6] = 84;
	diScale[6] = 4;
	eachEdge[6] = 4;
	mem[6] = 488; // 4 + 121 * 4
	// shared parameters
	shared[6] = 0;
	computationalWeight[6] = (120 + 120) * 4;

	// 8th layer: fully connected
	wi[7] = 1;
	hi[7] = 1;
	di[7] = 10;
	diScale[7] = 5;
	eachEdge[7] = 4;
	mem[7] = 430; // 5 + 85 * 5
	// shared parameters
	shared[7] = 0;
	computationalWeight[7] = (84 + 84) * 5;

	initialPos[0] = 0;
	for (int i = 1; i < nLayers; i++) {
		initialPos[i] += initialPos[i - 1] + (wi[i - 1] / wiScale[i - 1]) * (hi[i - 1] / hiScale[i - 1]) * (di[i - 1] / diScale[i - 1]); 
	}

	// input vertices + vertices with 1st layer conv result + vertices with 2nd layer subsampling result
	num_vertices = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 30 + /*FC6*/ 21 + /*Output*/ 2;

	Source_graph = GraphInit(num_vertices);
	Source_graph->enable_memory = 1;
	Source_graph->enable_vertex_labels = 0;
	Source_graph->enable_edge_weights = 1;
	Source_graph->enable_vertex_weights = 1;

	for (int v = 0; v < num_vertices; v++) {
		int layer = 0;
		for (int j = 0; j < nLayers; j++) {
			if (layer == nLayers - 1)
				break;
			if (v < initialPos[j + 1]) {
				break;
			} else {
				layer++;
			}
		}		
		GraphUpdateVertexWeight(Source_graph, v, computationalWeight[layer]);
	}

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

	// conv C5
	width_input = 5;
	height_input = 5;
	inputi = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ (((WI - WF) + 1)/2 - WF + 1) * (((HI - HF) + 1)/2 - WF + 1);
	conv_starts = inputi + width_input * height_input;
	scale = 16;
	BuildGraph2DConvolutionalLayerC5(Source_graph, conv_starts, width_input, height_input, inputi, scale, originalDF, mem[5], 5);	

	// FC6
	int output = 21;
	int input = 30;
	int fc_input_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF;
	int fc_output_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 30;
	BuildGraph2DFullConnectedLayerFC6(Source_graph, input, output, fc_input_starts, fc_output_starts, mem[6], computationalWeight[6], 6);

	// Output
	input = 21;
	output = 2;
	fc_input_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 30;
	fc_output_starts = /*input*/ WI * HI * DI + /*C1*/ ((WI - WF) + 1) * ((HI - HF) + 1) * DF + /*S2*/ ((WI - WF) + 1)/2 * ((HI - HF) + 1)/2 * DF + /*C3*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*S4*/ ((((WI - WF) + 1)/2 - WF) + 1) * ((((HI - HF) + 1)/2 - HF) + 1) * DF + /*C5*/ 30 + /*FC6*/ 21;
	BuildGraph2DFullConnectedLayerFC6(Source_graph, input, output, fc_input_starts, fc_output_starts, mem[7], computationalWeight[7], 7);

	WriteInputFile(Source_graph, num_vertices);
		
	PrintGraph(Source_graph);
	return 0;
}

/* __________________________________________________________*/
void BuildGraph2DFullConnectedLayerFC6(Graph Source_graph, int input, int output, int fc_input_starts, int fc_output_starts, int trainParam, int compWeight, int layer) {

	for (int i = 0; i < output; i++) {
		for (int j = 0; j < input; j++) {
			GraphInsertArc(Source_graph, fc_input_starts + j, fc_output_starts + i, eachEdge[layer]);
			GraphInsertArc(Source_graph, fc_output_starts + i, fc_input_starts + j, eachEdge[layer]);
		}
	}

	for (int i = fc_output_starts; i < fc_output_starts + output; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, trainParam);
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
									  int trainParam,
					int layer) {
	int filteri, filterj, convi, convj, inputi, inputj, conv_index, conv_width, conv_height;

	// conv result index starts after input ends
	conv_index = conv_starts;
	// conv result width
	conv_width = width_input;
	// conv result height
	conv_height = height_input;
	// conv filter
	int DF = di[layer] / diScale[layer];
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

	for (int i = conv_starts; i < conv_starts + di[layer] / diScale[layer]; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, trainParam);
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
	}
	int num_vertices = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF;
	for (int i = conv_index; i < num_vertices; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale * originalDF);
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
		GraphUpdateVertexWeight(Source_graph, i, computationalWeight[0]);
	}
	int num_vertices = WI * HI * DI + ((WI - WF) + 1) * ((HI - HF) + 1) * DF;
	for (int i = conv_index; i < num_vertices; i++) {
		GraphUpdateMemoryWeight(Source_graph, i, scale * originalDF + trainParam);
		GraphUpdateVertexWeight(Source_graph, i, computationalWeight[1]);
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
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", initialPos[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", shared[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", wi[layer] / wiScale[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", hi[layer] / hiScale[layer]);
	}
	fprintf(wif, "\n");
	for (int layer = 0; layer < nLayers; layer++) {
		fprintf(wif, "%d ", di[layer]);
	}

	//int initialPos[8] = {0, 156, 12, 1516, 32, 0, 0, 0};

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
