/*
* Source file of Portfolio Optimization Libraray, Version 1.0.0.0
* Created at 2012-05 -17, 15:35 PM by Wilfred You
*/

#include "ProbabilityModel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int ProbModelCalculateBinCount(double* lpValues, int nCount, double dBinWidth, double *minValue)
{
	int i;
	double min = lpValues[0];
	double max = min;
	double value = 0.0;

	for (i = 1; i < nCount; i++)
	{
		value = lpValues[i];

		if (value < min)
		{
			min = value;
		}
		else if (value > max)
		{
			max = value;
		}
	}

	*minValue = min;

	// Calculates the number of bins needed.
	return (int)ceil((max - min) / dBinWidth) + 2;
}

void ProbModelCreateProDensity(PROBMODELINFO* probModelInfo)
{
	int i, totalOccurrence;
	int binIndex, binCount;
	double normFactor, mean, variance, diff, bin;

	binCount = probModelInfo->nBinCount;

	// Adjusts the minimum bin so that it is an exact multiple of binWidth.
	probModelInfo->dMinPrice = probModelInfo->dBinWidth * (int)floor(probModelInfo->dMinPrice / probModelInfo->dBinWidth);

	// Initializes bins and frequencies.
	bin = probModelInfo->dMinPrice;
	for (i = 0; i < probModelInfo->nBinCount; i++)
	{
		probModelInfo->dBins[i] = bin;
		probModelInfo->dFreqs[i] = 0.0;
		bin = bin + probModelInfo->dBinWidth;
	}

	// Constructs the histogram for the observations.			  
	totalOccurrence = 0;
	binIndex = 0;
	for (i = 0; i < probModelInfo->nLength; i++)
	{
		// Calculates bin index for the current observation.  Bin must be positive since it is the delta
		// from the min value.
		bin = abs(probModelInfo->inValues[i] - probModelInfo->dMinPrice);
		binIndex = (int)floor(bin / probModelInfo->dBinWidth);

		// Safety check.
		if (binIndex < 0)
		{
			char msg[128];
			sprintf_s(msg, 128, "Bin index cannot be negative.  Bin = {%f}, Bin index = {%d}", bin, binIndex);
		}

		// Increment its frequency by 1.
		probModelInfo->dFreqs[binIndex] = probModelInfo->dFreqs[binIndex] + 1;

		// Also increments total occurrence by 1.
		totalOccurrence = totalOccurrence + 1;
	}

	// Normalizes frequencies to make it a proper density.
	normFactor = totalOccurrence * probModelInfo->dBinWidth;
	for (i = 0; i < probModelInfo->nBinCount; i++)
	{
		probModelInfo->dFreqs[i] = probModelInfo->dFreqs[i] / normFactor;
	}

	// Calculates mean observation value. 
	mean = 0.0;
	for (i = 0; i < probModelInfo->nBinCount; i++)
	{
		mean = mean + probModelInfo->dBins[i] * probModelInfo->dFreqs[i] * probModelInfo->dBinWidth;
	}
	probModelInfo->dMean = mean;

	// Calculates variance of observation value. 
	variance = 0.0;
	diff = 0.0;
	for (i = 0; i < probModelInfo->nBinCount; i++)
	{
		diff = probModelInfo->dBins[i] - mean;
		variance = variance + diff * diff * probModelInfo->dFreqs[i] * probModelInfo->dBinWidth;
	}

	probModelInfo->dVariance = variance;
	probModelInfo->dStandardDev = sqrt(variance);	
}

void ProbModelCreateProDistribution(PROBMODELINFO* probModelInfo, PROPOINT* pPoint)
{
	int i;
	double oldProb = 0;
	double newProb;

	for (i = 0; i < probModelInfo->nBinCount; i++)
	{
		// Calculates the probability that the underlying random variable is <= density.GetBin[i];
		newProb = oldProb + probModelInfo->dFreqs[i] * probModelInfo->dBinWidth;

		// Safety check.  Probability can never exceed 1.0
		if (newProb > 1.0)
		{
			newProb = 1.0;
		}

		// Creates the corresponding point, re-use existing point if possible.  Points are automatically 
		// sorted in both x and y values in ascending order.
		pPoint[i].x = probModelInfo->dBins[i];
		pPoint[i].y = newProb;

		// To be used for next probability.
		oldProb = newProb;
	}
}
