/*
* Source file of Portfolio Optimization Libraray, Version 1.0.0.0
* 
* Created at 2012-05 -17, 15:35 PM by Wilfred You
*
*/

// GA.cpp : Defines the entry point for the console application.
#include "..\..\Random\Random.h"
#include "GeneticAlgorithm.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <memory.h>
#include <time.h>
#include <Windows.h>
#include <assert.h>

/* Declaration of procedures used by this genetic algorithm */ 
static void Initialize(GAInfo*); 
static void Uninitialize(GAInfo*);
static void UpdateCost(Population * population, EvaluateFunc evfunc);
static void Reproduce(Population * population);
static void Crossover(Population * population);
static void Mutate(Population * population);
static void Swap(Chromosome** a, Chromosome** b);
static double Median(Chromosome** Chromosomes, int low, int high);
static void HalfInsertSort(Chromosome** Chromosomes, int length);
static void SortChromosomes(Chromosome** Chromosomes, int low, int high);
static void SortPopulation(Population * population);
static void Report(Population*); 
static void Resize(Population* population) ;
static void setGeneValue(Chromosome* chromosome, int index, double geValue);

#define TIME_ESTIMATE_START() {\
LARGE_INTEGER start;\
LARGE_INTEGER end ;\
LARGE_INTEGER frequency;\
if (!QueryPerformanceFrequency(&frequency)) \
{\
	return;\
}\
QueryPerformanceCounter(&start); 

#define TIME_ESTIMATE_END()	\
QueryPerformanceCounter(&end); \
printf("\nMain time cost:%f (s)", (double)(end.QuadPart - start.QuadPart) / (double)frequency.QuadPart); \
}

#define MAX_EVALUEATION_LENGTH	128

/*************************************************************/ 
/* setGeneValue: set the gene value */ 
/*************************************************************/ 
static void setGeneValue(Chromosome* chromosome, int index, double geValue)
{
	Gene* gene = &(chromosome->gene[index]);
	
    // Validates new value before assigning it.
    if (geValue < gene->geLower)
    {
        geValue = gene->geLower;
    }
    else if (geValue > gene->geUpper)
    {
        geValue = gene->geUpper;
    }

    // Checks for changes.
    if (geValue != gene->geValue)
    {
        gene->geValue = geValue;

        // Signals cost update is needed.
        chromosome->updated = false;
    }	
}

/***************************************************************/ 
/* Initialization function: Initializes the values of genes */ 
/* within the variables bounds. It also initializes (to zero) */ 
/* all cost values for each member of the population. It randomly */
/* generates values between these bounds for each gene of */ 
/* each genotype in the population. */ 
/***************************************************************/ 
static void Initialize(GAInfo* gaInfo) 
{ 
	int i; 
	int* mating_list;
	Chromosome* pChromList;
	Gene* pGeneList;

	pChromList = (Chromosome*)malloc((gaInfo->population.initsize)*sizeof(Chromosome));
	assert(NULL != pChromList);
	
	pGeneList = (Gene*)malloc((gaInfo->population.initsize)*(gaInfo->population.issuescope)*sizeof(Gene));
	assert(NULL != pGeneList);	

	mating_list = (int*)malloc((gaInfo->population.matingpoolsize)*sizeof(int));
	assert(NULL != mating_list);

	gaInfo->population.items = (Chromosome**)malloc((gaInfo->population.initsize)*sizeof(Chromosome*));
	assert(NULL != gaInfo->population.items);

	/* initialize variables within the bounds */ 
	for (i = 0; i < gaInfo->population.initsize; i++)
	{
		pChromList[i].cost = 0;
		pChromList[i].rCost = 0;
		pChromList[i].sProb = 0;
		pChromList[i].cProb = 0;
		pChromList[i].gene = &pGeneList[i * gaInfo->population.issuescope];
		gaInfo->gefunc(pChromList[i].gene, gaInfo->population.issuescope);
		pChromList[i].genecount = gaInfo->population.issuescope;
		pChromList[i].updated = false;
		pChromList[i].extId = i;
		gaInfo->population.items[i] = &pChromList[i];
	}
	//gaInfo->expifunc(pChromList);

	gaInfo->population.costListSize = MAX_EVALUEATION_LENGTH;
	gaInfo->population.minCostList = (double*)malloc((gaInfo->population.costListSize) * sizeof(double));
	assert(NULL != gaInfo->population.minCostList);

	gaInfo->population.genelist = pGeneList;
	gaInfo->population.chlist = pChromList;

	gaInfo->population.actualsize = gaInfo->population.initsize;
	gaInfo->population.mom_list = mating_list;
	gaInfo->population.dad_list = mating_list + gaInfo->population.matingpoolsize/2;
	gaInfo->population.generation = 0;
} 

/***************************************************************/ 
/* Uninitialize function: This function is to free all */
/* the resource allocated for GA processing */ 
/***************************************************************/ 
static void Uninitialize(GAInfo* gaInfo)
{
	if(NULL != gaInfo->population.items)
	{
		if(NULL != gaInfo->population.genelist)
		{
			free(gaInfo->population.genelist);
			gaInfo->population.genelist = NULL;
		}

		if(NULL != gaInfo->population.chlist)
		{
			free(gaInfo->population.chlist);
			gaInfo->population.chlist = NULL;
		}
		free(gaInfo->population.items);
		gaInfo->population.items = NULL;
	}

	if(NULL != gaInfo->population.mom_list)
	{
		free(gaInfo->population.mom_list);
		gaInfo->population.mom_list = NULL;
	}

	if(NULL != gaInfo->population.minCostList)
	{
		free(gaInfo->population.minCostList);
		gaInfo->population.minCostList = NULL;
	}
}

static void UpdateCost(Population * population, EvaluateFunc evfunc, void* context)
{
#if 0//ndef D_MULTI_THREAD
	int i;

	for (i = 0; i < population->actualsize; i++)
	{
		evfunc(population->items[i], context);
	}
#else
	evfunc(population->chlist, context);
#endif
}
/***************************************************************/ 
/* Reproduce: selects two parents that take part in */ 
/* the crossover. Implements a single points crossover */ 
/***************************************************************/ 
static void Reproduce(Population * population)
{
	// Declares local variables needed if can get this far.
	int i, mom_count, dad_count;
	int outstanding;
	double baseCost, normCost, totalCost, prob, thresholdProb;
	bool isDadNeeded = false;
	RANDINFO randinfo;
	
	// Safety check.
	if ((population->matingpoolsize%2 != 0) || (population->matingpoolsize <= 0) || (2*population->matingpoolsize > population->normalsize))
	{
		return;
	}

    // *************************************************************************
    // Select parents by cost or rank weighting depending on the value totalCost.
    // *************************************************************************
    // Lowest cost of the discarded chromosome.
	population->actualsize = population->normalsize;
    baseCost = population->items[population->matingpoolsize]->cost;

    // Normalized costs.
    totalCost = 0.0;
    for (i = 0; i < population->matingpoolsize; i++)
    {
        // Computes the normalized cost of the current chromosome.
        normCost = population->items[i]->cost - baseCost;

        // Saves the normalized cost for later use.
        //normCosts[i] = normCost;
		population->items[i]->rCost = normCost;

        // Sum of normalized costs.  This value will decide whether to use cost or
        // rank weighting later on.
        totalCost = totalCost + normCost;
    }

    // Computes the probability array.
    if (totalCost <= 0.000001)
    {
        // Computes probability by rank weighting scheme since totalCost = 0.
        // totalCost = 1 + 2 + 3 + ... + MatingPool
        totalCost = (1.0 + (double)population->matingpoolsize) * (double)population->matingpoolsize / 2.0;

        for (i = 0; i < population->matingpoolsize; i++)
        {
        	population->items[i]->sProb = (double)(population->matingpoolsize - i) / totalCost;
            //chromosome.sfitness[i] = (double)(population->matingpoolsize - i) / totalCost;
        }
    }
    else
    {
        // Computes probability by cost weighting scheme.
        for (i = 0; i < population->matingpoolsize; i++)
        {
            population->items[i]->sProb = abs(population->items[i]->rCost / totalCost);
        }
    }

    // Computes cumulative probabilities.
    population->items[0]->cProb = population->items[0]->sProb;
    for (i = 1; i < population->matingpoolsize; i++)
    {
        prob = population->items[i - 1]->cProb + population->items[i]->sProb;
        population->items[i]->cProb = prob;
    }

    // Selects matingPoolSize number of chromosomes for mating.
    outstanding = population->matingpoolsize;
    isDadNeeded = false; 
	mom_count = 0;
	dad_count = 0;
	InitRandomSpecial(&randinfo);
    while (outstanding > 0)
    {
        // Random threshold probability. 0 <= thresholdProb <= 1.
        thresholdProb = RandomUniform(&randinfo);

		// Finds the first cumulative probability that is greater than thresholdProb.
		for (i = 0; i < population->matingpoolsize; i++)
		{
            if (population->items[i]->cProb > thresholdProb)
            {
                if (isDadNeeded)
                {
					population->dad_list[dad_count++] = i;
                }
                else
                {
					population->mom_list[mom_count++] = i;
                }

                // Toggle flag so that mom and dad chromosomes are selected on an alternating basis.
                isDadNeeded = !isDadNeeded;

                // One less chromosome to be selected for pairing.
                outstanding--;
                break;
            }
        }
    }
}

/**************************************************************/ 
/* Crossover: performs crossover of each pair of selected parents. */ 
/**************************************************************/ 
static void Crossover(Population * population)
{
	int alpha, i, j, parents_count;
	int genesPerChromosome, discardCnt, matingPoolSize;
	Chromosome* mom; 
	Chromosome* dad;
	Chromosome* baby1;
	Chromosome*	baby2;
	double beta, newP1, newP2;
	Chromosome** babies;
	RANDINFO randinfo;
	
	// Deduces the mating pool size.
	matingPoolSize = population->matingpoolsize;
	
	// Safety check.			
	if ((matingPoolSize <= 0) || (2*matingPoolSize > population->normalsize) || (matingPoolSize % 2 != 0))
	{
		return;
	}
	
	genesPerChromosome = population->issuescope;
	
	// Calculates the number of chromosomes to be discarded.
	discardCnt = population->normalsize - matingPoolSize;
	
	// ***************************************************************************** //
	// * For speed and efficiency reasons, we should avoid continuously destroying
	// * and creating new chromosomes.	As such, existing chromosome objects should
	// * be re-used whenever possible.
	// ***************************************************************************** //
	babies = &population->items[matingPoolSize];
	
	// Creates random number generator to randomly select the parents for mating.
	InitRandomSpecial(&randinfo);
	parents_count = matingPoolSize/2;
	for (i = 0; i < parents_count; i++)
	{
		// Random crossover point.	0 <= alpha < genesPerChromosome.
		alpha = RandomNum(&randinfo, genesPerChromosome);
	
		// Blending factor. 0 <= beta < 1.
		beta = RandomUniform(&randinfo);
	
		// Retrieves the pair for mating.
		mom = population->items[population->mom_list[i]];
		dad = population->items[population->dad_list[i]];
	
		// Computes new parameter values.
		newP1 = mom->gene[alpha].geValue - beta * (mom->gene[alpha].geValue - dad->gene[alpha].geValue);
		newP2 = dad->gene[alpha].geValue + beta * (mom->gene[alpha].geValue - dad->gene[alpha].geValue);
	
		// Life continues...
		baby1 = babies[2 * i];
		baby2 = babies[2 * i + 1];

		setGeneValue(baby1, alpha, newP1);
		setGeneValue(baby2, alpha, newP2);
		//baby1.gene[alpha].geValue = newP1;
		//baby2.gene[alpha] = newP2;
		for (j = 0; j < alpha; j++)
		{
			setGeneValue(baby1, j, mom->gene[j].geValue);
			setGeneValue(baby2, j, dad->gene[j].geValue);
			//baby1[j] = mom[j];
			//baby2[j] = dad[j];
		}
	
		for (j = alpha+1; j < genesPerChromosome; j++)
		{
			setGeneValue(baby1, j, dad->gene[j].geValue);
			setGeneValue(baby2, j, mom->gene[j].geValue);
			//baby1[j] = dad[j];
			//baby2[j] = mom[j];
		}
	}			 
}

/**************************************************************/ 
/* Mutation: Random uniform mutation. A variable selected for */ 
/* mutation is replaced by a random value between lower and */ 
/* upper bounds of this variable */ 
/**************************************************************/ 
static void Mutate(Population * population)
{
	int i, chrmIdx, geneIdx;
	int genesPerChromosome;
	int mutateCnt;
	double upper, lower;
	RANDINFO randinfo;

  // Safety check.  Mutation is only possible if population is greater than 2.
    if (population->normalsize <= 2)
    {
        return;
    }

    // Short cut.
    genesPerChromosome = population->issuescope;

    // Total number of genes to be mutated
    mutateCnt = (int)(population->normalsize * genesPerChromosome * population->probmutation);

    // Random number generator for randomly selecting the genes to be mutated.
	InitRandomSpecial(&randinfo);

    // Mutates mutateCnt number of chromosomes.            
    for (i = 0; i < mutateCnt; i++)
    {
        // Selects the chromosome for mutation.  Note that
        //   1 <= chrmIdx < population.Size
        // in order to avoid the first(best) chromosome.
        chrmIdx = RandomNum(&randinfo, (population->normalsize - 2)) + 1;

        //  Selects the gene for mutation.  0 <= geneIdx < genesPerChromosome.
        geneIdx = RandomNum(&randinfo, genesPerChromosome);

        // Mutation takes place.  The new value is a random number within the
        // specified range of the gene.
		upper = population->items[chrmIdx]->gene[geneIdx].geUpper;
		lower = population->items[chrmIdx]->gene[geneIdx].geLower;	
		setGeneValue(population->items[chrmIdx], geneIdx, (RandomUniform(&randinfo)*(upper - lower)+lower));
    }
}

/*****************************************************************/
/* Swap function: swap two chromosome pointers*/
/*****************************************************************/
static void Swap(Chromosome** a, Chromosome** b)
{
    Chromosome* temp = *a;
    *a = *b;
    *b = temp;
}

/*****************************************************************/
/* Median function: using the median value for partion to support quicksort*/
/*****************************************************************/
static double Median(Chromosome** Chromosomes, int low, int high) 
{
    int mid = (low + high) / 2;

    if (Chromosomes[low]->cost > Chromosomes[mid]->cost)
        Swap(&Chromosomes[low], &Chromosomes[mid]);
    if (Chromosomes[low]->cost > Chromosomes[high]->cost)
        Swap(&Chromosomes[low], &Chromosomes[high]);
    if (Chromosomes[mid]->cost > Chromosomes[high]->cost)
        Swap(&Chromosomes[mid], &Chromosomes[high]);
        
    Swap(&Chromosomes[mid], &Chromosomes[low+1]); 

    return Chromosomes[low+1]->cost;
}

/*****************************************************************/
/* HalfInsertSort function: using HalfInsert to sort the chromosome list*/
/*****************************************************************/
static void HalfInsertSort(Chromosome** Chromosomes, int length) 
{
	int i, j;
	int low, high, mid;
	Chromosome* temp;
	
	for (i=1; i<length; i++)
	{
	    temp = Chromosomes[i];
	    low = 0;
	    high = i - 1;
	    while (low <= high) 
	    {
	          mid = (low + high) / 2;
	          if (Chromosomes[mid]->cost > temp->cost)
	                high = mid - 1;
	          else
	                low = mid + 1;
	    } 
		
	    for (j=i-1; j>high; j--)
	    {
			Chromosomes[j+1] = Chromosomes[j];
	    }
		
	    Chromosomes[high+1] = temp; 
	}
}

/*****************************************************************/
/* SortChromosomes function: using QuickSort to sort the chromosome list*/
/*****************************************************************/
static void SortChromosomes(Chromosome** Chromosomes, int low, int high)
{ 
	if (high-low >= 5)
	{
		double pivot = Median(Chromosomes, low, high); 
		int i = low + 1;
		int j = high;
		
		while (i < j)
		{
			while (Chromosomes[++i]->cost < pivot) {}
			while (Chromosomes[--j]->cost > pivot) {}
			if (i < j)
				Swap(&Chromosomes[i], &Chromosomes[j]);
		}
			
		Swap(&Chromosomes[j], &Chromosomes[low+1]); 

		SortChromosomes(Chromosomes, low, j-1);
		SortChromosomes(Chromosomes, j+1, high);
	}
	else
	{
		HalfInsertSort(&Chromosomes[low], high-low+1);
	}
}

/***************************************************************/ 
/* SortPopulation function: This function keeps sort the population*/ 
/* in ascending order by cost of the chromosome. Note that best */ 
/* individual will be re-older in the population*/ 
/***************************************************************/ 
static void SortPopulation(Population * population)
{
	SortChromosomes(population->items, 0, population->actualsize-1);
}

/***************************************************************/ 
/* Report function: Reports progress of the simulation. Data */ 
/* dumped into the output file are separated by commas */ 
/***************************************************************/ 
static void Report(Population* population) 
{ 
//	int i; 
	double minCost; /* min. population cost */ 
	double avg = 0.0; /* avg population fitness */ 
	double stddev = 0.0; /* std. deviation of population cost */ 
	double sum_square = 0.0; /* sum of square for std.  */ 
	double sum = 0.0; /* total population fitness */ 

#if 0
	for (i = 0; i < population->normalsize; i++) 
	{ 
		sum += population->items[i]->cost; 
	} 
	avg = sum/(double)population->normalsize; 

	for (i = 0; i < population->normalsize; i++) 
	{
		sum_square += ((population->items[i]->cost - avg)*(population->items[i]->cost - avg));
	}
	stddev = sqrt(sum_square /(population->normalsize - 1)); 
#endif
	minCost = population->items[0]->cost; 

	population->minCost = minCost;
	population->avg = avg;
	population->stddev = stddev;
} 

/**************************************************************/ 
/* Report function: Resize the population list according to */ 
static void Resize(Population* population, void* context) 
{

}

/**************************************************************/ 
/* Main function: Each generation involves selecting the best */ 
/* members, performing crossover & mutation and then */ 
/* evaluating the resulting population, until the terminating */ 
/* condition is satisfied */ 
/**************************************************************/ 
void GALaunch(GAInfo* gaInfo)
{
	Initialize(gaInfo);

	UpdateCost(&(gaInfo->population), gaInfo->evfunc, gaInfo->context);
	SortPopulation(&(gaInfo->population));
	while(gaInfo->population.generation < gaInfo->population.maxgeneration) 
	{ 
		gaInfo->population.generation++; 
		Reproduce(&(gaInfo->population)); 
		Crossover(&(gaInfo->population)); 
		Mutate(&(gaInfo->population)); 
		UpdateCost(&(gaInfo->population), gaInfo->evfunc, gaInfo->context);
		SortPopulation(&(gaInfo->population));
		Report(&(gaInfo->population)); 
		if(gaInfo->icfunc(&(gaInfo->population), gaInfo->context))
		{
			break;
		}
	} 
	
	Uninitialize(gaInfo);
}
