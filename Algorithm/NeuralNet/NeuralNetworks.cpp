/*
* Source file of Portfolio Optimization Libraray, Version 1.0.0.0
* 
* Created at 2012-05 -17, 15:35 PM by Wilfred You
*
* NeuralNetworks.cpp - Implementation of single level neural networks algorithm
*/
#include "NeuralNetworks.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

void NetLayerAddNeurons(NetLayer* pLayer, double* receptors, double* synweights)
{
	int i;
	
	Neuron* pNeurons = (Neuron*)malloc(pLayer->listSize * sizeof(Neuron));

	for(i = 0; i < pLayer->listSize; i++)
	{
		pNeurons[i].id = i;
		pNeurons[i].number = pLayer->inchannel;
		pNeurons[i].soma = 0;
		pNeurons[i].bias = 0;
		pNeurons[i].receptors = receptors;
		pNeurons[i].synweights = &(synweights[i* pLayer->inchannel]);
	}
	
	pLayer->pNeuronList = pNeurons;
}

NeuralNetworks* ConstructNeuralNetworks(NeuroSetting* settings)
{
	int i;
	int numlayers;
	int totalsynweights, synweights, numneurons;
	double* pSynWeights;
	double* pCurLayerSynWeights;
	double* pReceptors;
	NeuralNetworks* pNeuralNetworks;
	NetLayer* pLayers;

	numlayers = settings->numlayers;
	pNeuralNetworks = (NeuralNetworks*)malloc(sizeof(NeuralNetworks));
	pLayers = (NetLayer*)malloc(numlayers * sizeof(NetLayer));

	totalsynweights = 0;
	synweights = settings->inchannels;
	for(i = 0; i < numlayers; i++)
	{
		numneurons = settings->numneurons[i];
		
		pLayers[i].id = i;
		pLayers[i].inchannel = synweights;
		pLayers[i].listSize = numneurons;
		
		totalsynweights += synweights * numneurons;
		synweights = numneurons;
	}

	pReceptors = (double*)malloc(settings->inchannels * sizeof(double));
	pSynWeights = (double*)malloc(totalsynweights * sizeof(double));
	
	pNeuralNetworks->iLayerCnts = numlayers;
	pNeuralNetworks->pLayers = pLayers;
	pNeuralNetworks->pReceptors = pReceptors;
	pNeuralNetworks->pSynWeights = pSynWeights;
	pNeuralNetworks->iReceptorCnts = settings->inchannels;
	pNeuralNetworks->iOutputsCnt = settings->ouchannels;
	pNeuralNetworks->iTotaSynWeights = totalsynweights;
	pNeuralNetworks->fnactivation = settings->fnactivation;
	pNeuralNetworks->fnPostHandleCallback = settings->fnPostHandleCallback;
	
	pCurLayerSynWeights = pSynWeights;
	for(i = 0; i < numlayers; i++)
	{		
		NetLayerAddNeurons(&pLayers[i], pReceptors, pCurLayerSynWeights);
		
		pLayers[i].pAxons = (double*)malloc(pLayers[i].listSize * sizeof(double));

		pCurLayerSynWeights += (pLayers[i].inchannel * pLayers[i].listSize);
		pReceptors = pLayers[i].pAxons;
	}

	return pNeuralNetworks;
}

void DestructNeuralNetworks(NeuralNetworks* pNeuralNetworks)
{
	if(NULL != pNeuralNetworks)
	{
		if(NULL != pNeuralNetworks->pLayers)
		{
			int i;

			for(i = 0; i < pNeuralNetworks->iLayerCnts; i++)
			{
				if(NULL != pNeuralNetworks->pLayers[i].pNeuronList)
				{
					free(pNeuralNetworks->pLayers[i].pNeuronList);
					pNeuralNetworks->pLayers[i].pNeuronList = NULL;
				}

				if(NULL != pNeuralNetworks->pLayers[i].pAxons)
				{
					free(pNeuralNetworks->pLayers[i].pAxons);
					pNeuralNetworks->pLayers[i].pAxons = NULL;
				}
			}

			free(pNeuralNetworks->pLayers);
			pNeuralNetworks->pLayers = NULL;
		}

		if(NULL != pNeuralNetworks->pReceptors)
		{
			free(pNeuralNetworks->pReceptors);
			pNeuralNetworks->pReceptors = NULL;
		}

		if(NULL != pNeuralNetworks->pSynWeights)
		{
			free(pNeuralNetworks->pSynWeights);	
			pNeuralNetworks->pSynWeights = NULL;
		}

		if(NULL != pNeuralNetworks->inTrainSignals)
		{
			free(pNeuralNetworks->inTrainSignals);	
			pNeuralNetworks->inTrainSignals = NULL;
		}

		if(NULL != pNeuralNetworks->ouLearnResults)
		{
			free(pNeuralNetworks->ouLearnResults);	
			pNeuralNetworks->ouLearnResults = NULL;
		}

		free(pNeuralNetworks);
		pNeuralNetworks = NULL;
	}
}


void LoadTrainingSignals2Net(NeuralNetworks* pNeuralNetworks, double* pSingals, int iLength)
{
	pNeuralNetworks->inTrainSignals = (double*)malloc(iLength * sizeof(double));
	pNeuralNetworks->ouLearnResults = (double*)malloc(iLength * pNeuralNetworks->iOutputsCnt * sizeof(double));

	memcpy(pNeuralNetworks->inTrainSignals, pSingals, iLength * sizeof(double));
	pNeuralNetworks->iTrainSingalCnts = iLength;
	pNeuralNetworks->iTrainingTimes = (iLength - pNeuralNetworks->iReceptorCnts) + 1;
}

void DownloadNeuralNetWeights(NeuralNetworks* pNeuralNetworks, double* pSynWeights, int iLength)
{
	int i;

	printf("\nBest Cost - Weights list:\n");
	for(i = 0; i < pNeuralNetworks->iTotaSynWeights; i++)
	{
		printf("\n weights[%d] = %lf", i, pNeuralNetworks->pSynWeights[i]);
	}

	memcpy(pSynWeights, pNeuralNetworks->pSynWeights, iLength*sizeof(double));
}

void UpdateNetworkNeurons(NeuralNetworks* pNeuralNetworks, double* pSynWeights, int iLength)
{
	memcpy(pNeuralNetworks->pSynWeights, pSynWeights, iLength*sizeof(double));
}

double PostTrainingNetworks(NeuralNetworks* pNeuralNetworks)
{
	return pNeuralNetworks->fnPostHandleCallback(pNeuralNetworks);
}

double ComputeNeuronBias(Neuron* pNeuron)
{
	int i;
	double bias = 0;

	for(i = 0; i < pNeuron->number; i++)
	{
		bias += (pNeuron->receptors[i] * pNeuron->synweights[i]);
//		printf("w[%d] = %lf, r[%d] = %lf\n",i, pNeuron->synweights[i], i, pNeuron->receptors[i]);
	}
	bias += pNeuron->soma;

	return bias;
}

void ComputeLayerAxons(NetLayer* pLayer, fnActivation fnactivation)
{
	int i;
	double bias;

	for(i = 0; i < pLayer->listSize; i++)
	{
//		printf("\nNeurons[%d]:\n", i);
		bias = ComputeNeuronBias(&(pLayer->pNeuronList[i]));
		pLayer->pAxons[i] = fnactivation(bias);
	}	
}

void TrainingNetworks(NeuralNetworks* pNeuralNetworks)
{
	int i, j;
	double* pMidResult;

	pMidResult = pNeuralNetworks->ouLearnResults;
	for(i = 0; i < pNeuralNetworks->iTrainingTimes; i++)
	{
		memcpy(pNeuralNetworks->pReceptors, &(pNeuralNetworks->inTrainSignals[i]), pNeuralNetworks->iReceptorCnts * sizeof(double));
		for(j = 0; j < pNeuralNetworks->iLayerCnts; j++)
		{	
			ComputeLayerAxons(&pNeuralNetworks->pLayers[j], pNeuralNetworks->fnactivation);
		}

		memcpy(pMidResult, (pNeuralNetworks->pLayers[pNeuralNetworks->iLayerCnts-1]).pAxons, pNeuralNetworks->iOutputsCnt * sizeof(double));
		pMidResult += pNeuralNetworks->iOutputsCnt;
	}
}

void StartNeuralNetWorks(NeuralNetworks* pNeuralNetworks)
{
	TrainingNetworks(pNeuralNetworks);
}

void DownloadNeuralNetOutputs(NeuralNetworks* pNeuralNetworks, double* ouResults, int ilength)
{
	memcpy(ouResults, pNeuralNetworks->ouLearnResults, ilength*sizeof(double));
}