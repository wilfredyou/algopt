/*
* Source file of Portfolio Optimization Libraray, Version 1.0.0.0
* Created at 2012-05 -17, 15:35 PM by Wilfred You
*/

#ifndef _H_PROBABILITY_
#define _H_PROBABILITY_

typedef struct _tagPOINT
{
	double x;
	double y;
}PROPOINT;

typedef struct _tagPROBMODELINFO
{
	double* inValues;
	int nLength;
	int nBinCount;
	double dMinPrice;
	double dBinWidth;

	double dMean;
	double dVariance;
	double dStandardDev;
	double* dBins;
	double* dFreqs;	
}PROBMODELINFO;

int ProbModelCalculateBinCount(double* lpValues, int nCount, double dBinWidth, double* minValue);

void ProbModelCreateProDensity(PROBMODELINFO* probModelInfo);

void ProbModelCreateProDistribution(PROBMODELINFO* probModelInfo, PROPOINT* pPoint);
#endif
