#ifndef IRIT_OPT_CONFIG_H
#define IRIT_OPT_CONFIG_H

#include "inc_irit/cagd_lib.h"

#define IRIT_OPT_MAX_PATH 1024
#define IRIT_OPT_MAX_NAME 128

typedef struct OptRunInputStruct {
    char MacroShapePath[IRIT_OPT_MAX_PATH];
    int SimulationType;
    int LatticeType;
    CagdRType SkinThickness;
    int LatticeCountU, LatticeCountV, LatticeCountW;
    int *LatticeCountLists[3];
    int LatticeCountListLengths[3];
    CagdRType InitialThickness;
    int MeshSample;
    int MaxIterations;
    int ConvergenceIterations;
    CagdRType ThicknessChangeTolerance;
    CagdRType TargetVolumeFraction;
    CagdRType TargetVolumeFractionTol;
    int OptimizationAlgorithm;
} OptRunInputStruct;

typedef struct OptConfigInputStruct {
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
} OptConfigInputStruct;

typedef struct OptConfigStruct {
    OptRunInputStruct Run;
    OptConfigInputStruct Environment;
} OptConfigStruct;

int IritOptRun(const OptConfigStruct *Config, CagdBType GeometryOnly);

#endif
