#include "math.h"

#include "WolframRTL.h"

static WolframCompileLibrary_Functions funStructCompile;

static void * E0 = 0;

static void * E1 = 0;


static mbool initialize = 1;

#include "Lognormal_Mean.h"

DLLEXPORT int Initialize_Lognormal_Mean(WolframLibraryData libData)
{
if( initialize)
{
funStructCompile = libData->compileLibraryFunctions;
{
E0 = funStructCompile->getExpressionFunctionPointer(libData, "Hold[Function[List[mu, sig2], Power[E, Plus[mu, Times[Rational[1, 2], Power[sig2, 2]]]]]]");
}
if( E0 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
{
E1 = funStructCompile->getExpressionFunctionPointer(libData, "Hold[Function[List[mu, sig2], Times[Power[E, Plus[Times[2, mu], Power[sig2, 2]]], Plus[-1, Power[E, Power[sig2, 2]]]]]]");
}
if( E1 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
initialize = 0;
}
return 0;
}

DLLEXPORT void Uninitialize_Lognormal_Mean(WolframLibraryData libData)
{
if( !initialize)
{
initialize = 1;
}
}

DLLEXPORT int Lognormal_Mean(WolframLibraryData libData, mreal A1, mreal A2, MTensor *Res)
{
mreal R0_0;
mreal R0_1;
mreal R0_2;
mreal R0_3;
MTensor* T0_0;
MTensorInitializationData Tinit;
mreal *P0;
int err = 0;
Tinit = funStructCompile->GetInitializedMTensors(libData, 1);
T0_0 = MTensorInitializationData_getTensor(Tinit, 0);
R0_0 = A1;
R0_1 = A2;
{
int S0[2];
void * S1[2];
S0[0] = 3;
S1[0] = (void*) (&R0_0);
S0[1] = 3;
S1[1] = (void*) (&R0_1);
err = funStructCompile->evaluateFunctionExpression(libData, E0, 0, 0, 2, S0, S1, 3, 0, (void*) (&R0_2));
if( err)
{
goto error_label;
}
}
{
int S0[2];
void * S1[2];
S0[0] = 3;
S1[0] = (void*) (&R0_0);
S0[1] = 3;
S1[1] = (void*) (&R0_1);
err = funStructCompile->evaluateFunctionExpression(libData, E1, 0, 0, 2, S0, S1, 3, 0, (void*) (&R0_3));
if( err)
{
goto error_label;
}
}
{
mint S0[2] = {1, 2};
err = funStructCompile->MTensor_allocate(T0_0, 3, 2, S0);
if( err)
{
goto error_label;
}
P0 = MTensor_getRealDataMacro(*T0_0);
P0[0] = R0_2;
P0[1] = R0_3;
}
funStructCompile->MTensor_copy(Res, *T0_0);
error_label:
funStructCompile->ReleaseInitializedMTensors(Tinit);
funStructCompile->WolframLibraryData_cleanUp(libData, 1);
return err;
}

