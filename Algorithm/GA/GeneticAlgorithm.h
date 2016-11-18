/*
* Source file of Portfolio Optimization Libraray, Version 1.0.0.0
* Created at 2012-05 -17, 15:35 PM by Wilfred You
*/

#ifndef _H_GENERATICALGORITHM_
#define _H_GENERATICALGORITHM_

/* Structures & Enum definition*/
typedef struct
{
	double geValue;			/* a string of variables */ 
	double geUpper;			/* GT's variables upper bound */ 
	double geLower;			/* GT's variables lower bound */ 
}Gene;

typedef struct genotype		/* genotype (GT), a member of the population */ 
{ 
	bool updated;			/* set the flag to update cost */
	int genecount;			/* gene numbers per chromosome*/
	Gene* gene;				/* Gene in the chromosome*/
	double cost;			/* GT's cost */ 
	double rCost;			/* relative cost */ 
	double sProb;			/* selection probability */
	double cProb;			/* cumulative probability */ 
	int extId;				/* extend parameters id */
}Chromosome; 

typedef struct
{
	int initsize;			/* initialize size of population*/
	int actualsize;			/* population for cost */
	int normalsize;			/* population size */
	int matingpoolsize;		/* mating pool size */
	int issuescope;			/* no. of problem variables*/
	int generation;			/* current generation no. */ 
	int maxgeneration;		/* max. number of generations */ 
	
	double probmutation;	/* probability of mutation */ 

	double minCost;			/* min. population cost */ 
	double avg;				/* avg population cost */ 
	double stddev;			/* std. deviation of population cost */

	double* minCostList;		/* Hold the min cost of latest generation */
	int costListSize;		/* Hold the size of min cost list */

	int* mom_list;			/* hold the mom index for mating */
	int* dad_list;			/* hold the dad index for mating */

	Gene* genelist;			/* hold the whole gene list for population */
	Chromosome* chlist;		/* backup the list of chromosome */
	Chromosome** items;		/* hold the pointer array of chromosome of population*/
}Population;

/* Callback function and Global Static Declarations */
typedef bool (*IsConvergedFunc)(Population*, void*);
typedef void (*EvaluateFunc)(Chromosome*, void*);
typedef void (*GeneCreateFunc)(Gene*, int);
//typedef void (*ExtParamInstallFunc)(Chromosome*);

/* Global Structure for GA */
typedef struct  
{
	IsConvergedFunc icfunc;			/* Specific the target function*/
	EvaluateFunc evfunc;			/* Specific population evaluation method*/
	GeneCreateFunc gefunc;			/* Specific gene generate method*/
//	ExtParamInstallFunc expifunc;	/* Specific the function to install the ext param*/
	void* context;					/* Hold the context for evaluation process */

	Population population;
}GAInfo;

/* Public interface for Genetic Algorithm*/
void GALaunch(GAInfo* gaInfo);
#endif
