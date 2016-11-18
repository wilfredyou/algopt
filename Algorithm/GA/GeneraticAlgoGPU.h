#if !defined(_H_GENERATICALGOGPU_) 
#define _H_GENERATICALGOGPU_
#include "GeneraticAlgorithm.h"

extern "C" void UpdateCostGPU(Population * population, void* context);

extern "C" void UpdateOptmResultListGPU(Population * population, int bestIdx, void* context);

extern "C" void GenerateRandomCloud(void* context);

#endif
