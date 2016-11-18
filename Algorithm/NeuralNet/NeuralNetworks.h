#ifndef _H_NEURALNETWORKS_
#define _H_NEURALNETWORKS_

typedef double (*fnActivation)(double);
typedef double (*fnPostHandleCallback)(void*);

typedef struct _tagNeuron
{
	int		id;					/*Hold the index in layer*/
	int		number;				/*Hold the number of receptors*/
	double	bias;
	double	soma;
	double* receptors;
	double* synweights;
}Neuron;

typedef struct _tagNetLayer
{
	int		id;
	int		inchannel;
	int		listSize;
	double* pAxons;
	Neuron* pNeuronList;
}NetLayer;

typedef struct _tagNeuralNetworks
{
	NetLayer*	pLayers;
	int			iLayerCnts;

	double*		pReceptors;
	double*		pSynWeights;
	int			iReceptorCnts;
	int			iOutputsCnt;
	int			iTotaSynWeights;

	double*		inTrainSignals;
	double*		ouLearnResults;
	int			iTrainSingalCnts;
	int			iTrainingTimes;

	fnActivation fnactivation;
	fnPostHandleCallback fnPostHandleCallback;
}NeuralNetworks;

typedef struct _tagNeuroInfo
{
	double*	signals;
	int		length;

	int*	numneurons;
	int		numlayers;
	
	int		inchannels;
	int		ouchannels;
	
	fnActivation fnactivation;
	fnPostHandleCallback fnPostHandleCallback;
}NeuroSetting;

// Public Interfaces
NeuralNetworks* ConstructNeuralNetworks(NeuroSetting* settings);
void TrainingNetworks(NeuralNetworks* pNeuralNetworks);
void DestructNeuralNetworks(NeuralNetworks* pNeuralNetworks);
void LoadTrainingSignals2Net(NeuralNetworks* pNeuralNetworks, double* pSingals, int iLength);
void DownloadNeuralNetWeights(NeuralNetworks* pNeuralNetworks, double* pSynWeights, int iLength);
void UpdateNetworkNeurons(NeuralNetworks* pNeuralNetworks, double* pSynWeights, int iLength);
double PostTrainingNetworks(NeuralNetworks* pNeuralNetworks);
void StartNeuralNetWorks(NeuralNetworks* pNeuralNetworks);
void DownloadNeuralNetOutputs(NeuralNetworks* pNeuralNetworks, double* ouResults, int ilength);

#endif
