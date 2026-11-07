#ifndef IRIT_OPT_CONFIG_H
#define IRIT_OPT_CONFIG_H

#include "inc_irit/cagd_lib.h"

#define IRIT_OPT_MAX_PATH 1024
#define IRIT_OPT_MAX_NAME 128

typedef struct OptConfigStruct {
    char MacroShapePath[IRIT_OPT_MAX_PATH];
    int SimulationType;
    int LatticeType;
    CagdRType SkinThickness;
    int LatticeCountU, LatticeCountV, LatticeCountW;
    int *LatticeCountLists[3];
    int LatticeCountListLengths[3];
    int BodyStackAxis;          /* 0=U, 1=V, 2=W, -1=infer from counts. */
    CagdBType UseMicroStructFFD;
    CagdRType InitialThickness;
    int MeshSample;
    int MaxIterations;
    int ConvergenceIterations;
    CagdRType ThicknessChangeTolerance;
    CagdRType TargetVolumeFraction;
    CagdRType TargetVolumeFractionTol;
    int OptimizationAlgorithm;
    CagdBType UseTargetVolume;
    CagdBType GeometryOnly;

    char BaseFolder[IRIT_OPT_MAX_PATH];
    char WorkDir[IRIT_OPT_MAX_PATH];
    char AnsysExe[IRIT_OPT_MAX_PATH];
    char AnsysProject[IRIT_OPT_MAX_PATH];
    char AnsysJournal[IRIT_OPT_MAX_PATH];
    char InpFile[IRIT_OPT_MAX_PATH];
    char ItdOutputFile[IRIT_OPT_MAX_PATH];
    char SedFile[IRIT_OPT_MAX_PATH];
    char DataNewFile[IRIT_OPT_MAX_PATH];
    char DataPrevFile[IRIT_OPT_MAX_PATH];
    char DesignParametersFile[IRIT_OPT_MAX_PATH];
    char SummaryFile[IRIT_OPT_MAX_PATH];
    char VolumeFile[IRIT_OPT_MAX_PATH];
    CagdRType Scale, Density, YoungsModulus, PoissonRatio, SedValue;
    CagdRType MinimumThickness, SafetyFactorSed;
    CagdBType RemoveThinTiles;
} OptConfigStruct;

int IritOptRun(const OptConfigStruct *Config);

#endif
