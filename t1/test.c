/***************************************************************************
*   test.c - Iterative lattice optimization and heterogeneous tile driver. *
****************************************************************************
* (C) Asaf Tzur / based on IRIT library infrastructure.                    *
****************************************************************************
* Written as a local IRIT-oriented optimization module.                    *
***************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/triv_lib.h"
#include "inc_irit/mvar_lib.h"
#include "inc_irit/user_lib.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <windows.h>
#include <nlopt.h>

#ifdef _WIN32
#   include <direct.h>
#   define OPT_PATH_SEP '\\'
#   define OPT_STRICMP   _stricmp
#   define OPT_STRNICMP  _strnicmp
#else
#   include <dirent.h>
#   include <sys/stat.h>
#   include <sys/types.h>
#   include <unistd.h>
#   include <errno.h>
#   include <strings.h>
#   define OPT_PATH_SEP '/'
#   define OPT_STRICMP   strcasecmp
#   define OPT_STRNICMP  strncasecmp
#endif

/***************************************************************************
*                              Global Limits.                              *
***************************************************************************/

#define OPT_MAX_PATH            1024
#define OPT_MAX_LINE            2048
#define OPT_MAX_NAME            128
#define OPT_MAX_ROWS            100000
#define OPT_MAX_COLS            100
#define OPT_MAX_NODE_PER_ELEM   32
#define OPT_INVALID_REAL        (-1.0e300)
#define OPT_MAX_TILE_THICKNESS  0.5
#define OPT_MIN_TILE_THICKNESS  0.001

/***************************************************************************
*                        Input File Structures.                            *
***************************************************************************/

typedef struct OptRunInputStruct {
    char MacroShapePath[OPT_MAX_PATH];

    int SimulationType;         /* 0: emulation, 1: ANSYS, 2: CalculiX. */

    int LatticeType;
    int LatticeCountU;
    int LatticeCountV;
    int LatticeCountW;
    int *BodyCounts[3];
    int BodyCountLengths[3];

    CagdRType InitialThickness;

    int MeshType;
    int MeshSample;

    int MaxIterations;
    int ConvergenceIterations;      /* Required consecutive small-change iterations. */
    CagdRType ThicknessChangeTolerance; /* Absolute max thickness-change tolerance. */
    CagdRType TargetVolumeFraction;
    CagdRType TargetVolumeFractionTol;
    //CagdRType TargetDeflection;
    //CagdRType TargetMass;
    //CagdRType DeflectionTolerance;
    //CagdRType MassTolerance;

    int OptimizationAlgorithm;
    int AlgoType;

    char RunFolderName[OPT_MAX_NAME];
} OptRunInputStruct;

typedef struct OptConfigInputStruct {
    char BaseFolder[OPT_MAX_PATH];
    char WorkDir[OPT_MAX_PATH];

    char AnsysExe[OPT_MAX_PATH];
    char AnsysProject[OPT_MAX_PATH];
    char AnsysJournal[OPT_MAX_PATH];
    char IritExe[OPT_MAX_PATH];

    char InpFile[OPT_MAX_PATH];
    char ItdOutputFile[OPT_MAX_PATH];
    char StressFile[OPT_MAX_PATH];
    char DeflectionFile[OPT_MAX_PATH];
    char SedFile[OPT_MAX_PATH];
    char DataNewFile[OPT_MAX_PATH];
    char DataPrevFile[OPT_MAX_PATH];
    char DesignParametersFile[OPT_MAX_PATH];
    char SummaryFile[OPT_MAX_PATH];
    char VolumeFile[OPT_MAX_PATH];

    CagdRType Scale;
    CagdRType Density;
    CagdRType YoungsModulus;
    CagdRType PoissonRatio;
    CagdRType SedValue;
    CagdRType YieldStress;

    CagdRType MinimumThickness;
    CagdRType SafetyFactorSed;
    CagdRType SafetyFactorStress;

    int NodesPerElement;

    int DefaultL;
    int DefaultLMin;
    int DefaultLDef;

    int ElementsPerTile;
    int MeshSampleTotal;

    /* 0 preserves the historical tile-0 ownership order; 1 uses the
     * component order verified from exported INP element centroids. */
    int Tile0OwnerMapping;

    CagdBType RemoveThinTiles;
} OptConfigInputStruct;

/***************************************************************************
*                   INP / Simulation Data Structures.                      *
*                                                                          *
* These structures hold the parsed mesh, simulation values, and the tile   *
* bookkeeping used by the optimization stages.                             *
***************************************************************************/

typedef struct OptNodeStruct {
    int Id;
    CagdRType X, Y, Z;
    CagdRType Stress;
} OptNodeStruct;

typedef struct OptElementStruct {
    int Id;
    int NodeIds[OPT_MAX_NODE_PER_ELEM];
    int NodeCount;
    CagdRType AvgStress;
    CagdRType MinStress;
    CagdRType MaxStress;
    CagdRType SED;
    CagdRType Volume;
    CagdRType StrainEnergy;
    CagdRType ComplianceContribution;
} OptElementStruct;

typedef struct OptSubTileStruct {
    OptElementStruct **Elements;
    int ElementCount;
    int ElementCapacity;

    CagdRType AvgStress;
    CagdRType MinStress;
    CagdRType MaxStress;
    CagdRType Thickness;
    CagdRType SED;
    CagdRType MaxSED;
    CagdRType Volume;
    CagdRType StrainEnergy;
    CagdRType ComplianceContribution;
    CagdRType Sensitivity;
    CagdRType Gain;
} OptSubTileStruct;

typedef struct OptTileStruct {
    char Name[OPT_MAX_NAME];
    OptSubTileStruct *SubTiles;
    int SubTileCount;

    CagdRType AvgStress;
    CagdRType MinStress;
    CagdRType MaxStress;
    CagdRType Volume;
    CagdRType StrainEnergy;
    CagdRType ComplianceContribution;
    CagdRType Gain;
} OptTileStruct;

typedef struct OptNodeArrayStruct {
    OptNodeStruct *Items;
    int Count;
    int Capacity;
} OptNodeArrayStruct;

typedef struct OptElementArrayStruct {
    OptElementStruct *Items;
    int Count;
    int Capacity;
} OptElementArrayStruct;

typedef struct OptScalarFieldStruct {
    CagdRType *Values;
    CagdBType *HasValue;
    int Capacity;
} OptScalarFieldStruct;

/***************************************************************************
*                  Legacy IRIT CSV / Microstructure State.                 *
*                                                                          *
* This section preserves the data layout consumed by the geometry          *
* generator.                                                               *
***************************************************************************/

typedef CagdRType (*LclThicknessFuncCBType)(CagdRType u,
                                            CagdRType v,
                                            CagdRType w);

typedef struct OptWorkflowStateStruct OptWorkflowStateStruct;

typedef struct OptBodyStruct {
    char Name[OPT_MAX_NAME];
    IritTrivTVStruct *Map;      /* Owned copy, one trivariate per body. */
    int Counts[3];
    int TileOffset, TileCount;
    CagdRType Min[3], Max[3];
} OptBodyStruct;

typedef struct OptMeshOwnerStruct {
    int Tile, SubTile;          /* One record per exported Bezier patch. */
} OptMeshOwnerStruct;

typedef struct UserMicroLocalDataStruct {
    IritTrivTVStruct *DefMap;
    CagdRType BndryThickness;
    LclThicknessFuncCBType ThicknessFuncCB;
    OptWorkflowStateStruct *Opt;
    OptBodyStruct *Body;
} UserMicroLocalDataStruct;

typedef struct OptLegacyStateStruct {
    float **Data;              /* Rows after the 4 header lines. */
    int Rows;
    int Cols;

    int LType;                 /* Legacy CSV row 1. */
    int Ua, Va, Wa;            /* Legacy CSV row 2. */
    int MeshSample;            /* Legacy CSV row 3. */
    int MeshType;              /* Legacy CSV row 4. */

} OptLegacyStateStruct;

typedef struct OptMMADataStruct {
    int Count;
    CagdRType *PrevX;
    CagdRType *ObjGrad;
    CagdRType *VolGrad;
    CagdRType ObjAtPrev;
    CagdRType VolAtPrev;
    CagdRType TargetVolume;
} OptMMADataStruct;


/***************************************************************************
*                         Workflow State Structure.                        *
***************************************************************************/

struct OptWorkflowStateStruct {
    OptRunInputStruct Run;
    OptConfigInputStruct Config;

    int L;
    int LMin;
    int LDef;
    int MeshSample;
    int ElementsPerTile;

    CagdRType MaxDef;
    CagdRType Mass;
    CagdRType SedAvg;
    CagdRType SedR;
    CagdRType DefMapVol;
    CagdRType LatticeVol;
    CagdRType LastThicknessChange;
    int ConsecutiveSmallChanges;

    char TimeIritStart[16];
    char TimeIritFinish[16];
    char TimeSimulationFinish[16];
    char FolderName[OPT_MAX_NAME];

    OptNodeArrayStruct Nodes;
    OptElementArrayStruct Elements;
    OptScalarFieldStruct StressData;
    OptScalarFieldStruct DeflectionData;
    OptScalarFieldStruct SedData;
    OptScalarFieldStruct VolumeData;

    OptBodyStruct *Bodies;
    int BodyCount, TotalTiles;
    OptMeshOwnerStruct *MeshOwners;
    int MeshOwnerCount, MeshOwnerCapacity;
    OptTileStruct *TileStructure;    /* Persistent tile/subtile topology.
                                       * Analysis fields and FE links are refreshed
                                       * each iteration; thickness is retained. */
    CagdRType *Results;              /* Flat array [iter][3]. */
    OptLegacyStateStruct Legacy;     /* Per-run IRIT CSV/callback state. */
};


#define RESULT_AT(Opt, iter, col) \
    ((Opt)->Results[(iter) * 3 + (col)])

/***************************************************************************
*                           Forward Declarations.                          *
***************************************************************************/

 /* ---------- config / parsing ---------- */
static void OptSetDefaultRunInputs(OptRunInputStruct *Run);
static void OptSetDefaultConfigInputs(OptConfigInputStruct *Cfg);
static int  LoadRunInputs(const char *FileName, OptRunInputStruct *Run);
static int  LoadConfigInputs(const char *FileName, OptConfigInputStruct *Cfg);
static int  ValidateInputs(const OptRunInputStruct *Run, const OptConfigInputStruct *Cfg);
static void ResolveLatticeLayout(OptWorkflowStateStruct *Opt);
static void InitializeWorkflowState(OptWorkflowStateStruct *Opt);
static void FreeWorkflowState(OptWorkflowStateStruct *Opt);

/* ---------- generic helpers ---------- */
static void OptFatal(const char *Msg);
static void OptFatalF(const char *Fmt, const char *Arg);
static char *OptTrim(char *Str);
static int OptIsBlank(const char *Str);
static void OptStripUtf8Bom(char *Line);
static void OptPathJoin(char *Dst,
                        size_t DstSize,
                        const char *A,
                        const char *B);
static const char *OptPathFileName(const char *Path);
static void OptMakeIndexedFileName(const char *FileName,
                                   int Index,
                                   char *Out,
                                   size_t OutSize);
static int OptFileExists(const char *FileName);
static void OptFormatNowHMS(char Buffer[16]);
static void OptMakeRunFolderName(char Buffer[OPT_MAX_NAME]);
static int OptWriteLine(const char *FileName, const char *Content, int Append);
static int OptParseInt(const char *Text, int *Value);
static int OptParseReal(const char *Text, CagdRType *Value);
static int OptParseBool(const char *Text);
static void *OptAllocZero(size_t Count, size_t Size);
static void *OptReallocCopy(void *OldData,
                            size_t OldBytes,
                            size_t NewBytes);
static int ParseKeyValueFile(const char *FileName,
                             int (*Handler)(const char *Key,
                                            const char *Value,
                                            void *Ctx),
                             void *Ctx);
static int RunInputHandler(const char *Key, const char *Value, void *Ctx);
static int ConfigInputHandler(const char *Key, const char *Value, void *Ctx);
static CagdRType ComputeObjVol(IritPrsrObjectStruct* PObj);
static int LatticeBy_VolumeFraction(OptWorkflowStateStruct *Opt);

/* ---------- dynamic containers ---------- */
static void NodeArrayReset(OptNodeArrayStruct *Arr);
static void ElementArrayReset(OptElementArrayStruct *Arr);
static void ScalarFieldReset(OptScalarFieldStruct *Field);
static int ScalarFieldEnsureCapacity(OptScalarFieldStruct *Field, int RequiredId);
static void ScalarFieldSet(OptScalarFieldStruct *Field, int Id, CagdRType Value);
static CagdBType ScalarFieldTryGet(const OptScalarFieldStruct *Field,
                                   int Id,
                                   CagdRType *Value);
static int NodeArrayPush(OptNodeArrayStruct *Arr, const OptNodeStruct *Node);
static int ElementArrayPush(OptElementArrayStruct *Arr,
                            const OptElementStruct *Elem);
static int SubTileAddElementPtr(OptSubTileStruct *SubTile,
                                OptElementStruct *Elem);
static void DestroyTileStructure(OptWorkflowStateStruct *Opt);
static int InitTileStructure(OptWorkflowStateStruct *Opt);
static void ClearTileAnalysisData(OptWorkflowStateStruct *Opt);

/* ---------- workflow core ---------- */
static int  InitLegacyDataFile(OptWorkflowStateStruct *Opt);
static int  Optimization_Main(int argc, char **argv);
static int  RunWorkflow(OptWorkflowStateStruct *Opt);
static int  RunIRITIntegrated(OptWorkflowStateStruct *Opt, int mode);
static int  RunSimulation(OptWorkflowStateStruct *Opt);
static int  RunSimulationEmulation(OptWorkflowStateStruct *Opt);
static int  RunAnsysBatch(OptWorkflowStateStruct *Opt);
static int  RunCalculixBatch(OptWorkflowStateStruct *Opt);
static int  ProcessIterationData(OptWorkflowStateStruct *Opt, int Iter);
static CagdBType EndConditionMet(const OptWorkflowStateStruct *Opt);

/* ---------- results / file operations ---------- */
static int OrganizeResults(OptWorkflowStateStruct *Opt,
                           const char *SourceFolder,
                           const char *TargetBaseFolder,
                           int Iter);
static int ReadLastRowOfDesignParameters(const char *FileName,
                                         CagdRType *Deflection,
                                         CagdRType *Mass);
static int CopyDirectoryRecursive(OptWorkflowStateStruct *Opt,
                                  const char *SourceDir,
                                  const char *TargetDir);
static int EnsureDirectoryRecursive(const char *DirName);
static int RunCommandLineBlocking(const char *CommandLine,
                                  const char *WorkingDir);
static int RunCommandLineBlockingToFile(const char *CommandLine,
                                        const char *WorkingDir,
                                        const char *LogFile);

/* ---------- simulation / data parsing ---------- */
static int ReadINP(OptWorkflowStateStruct *Opt, const char *INPPath);
static int ReadSimulationData(const char *FileName,
                              OptScalarFieldStruct *Field);
static CagdRType *ReadTileData(const OptWorkflowStateStruct *Opt, const char *FileName, int ExpectedCount);
static OptNodeStruct *FindNodeById(OptWorkflowStateStruct *Opt, int NodeId);
static void UpdateTileAnalysisData(OptWorkflowStateStruct *Opt);
static CagdRType DesignUpdate_OC(OptWorkflowStateStruct *Opt,
                                  int Count,
                                  const CagdRType *Prev,
                                  int Iter);
static CagdRType DesignUpdate_SED_Ratio(OptWorkflowStateStruct *Opt,
                             int Count,
                             const CagdRType *Prev,
                             int Iter);
static CagdRType DesignUpdate_ASAF(OptWorkflowStateStruct* Opt,
    int Count,
    const CagdRType* Prev,
    int Iter);
static CagdRType ComputeMaxThicknessChange(const OptWorkflowStateStruct *Opt,
                                           const CagdRType *Prev,
                                           int Count);
static CagdRType DesignUpdate_MMA(OptWorkflowStateStruct* Opt,
    int Count,
    const CagdRType* Prev,
    int Iter);
static double MMAObjective(unsigned n,
                           const double *x,
                           double *grad,
                           void *Data);
static double MMAVolumeConstraint(unsigned n,
                                  const double *x,
                                  double *grad,
                                  void *Data);

static void NeighbourMatching(OptWorkflowStateStruct *Opt);
static int CompareNeighbourEntries(const void *P1, const void *P2);
static void NeighbourMatchingB(OptWorkflowStateStruct *Opt,
                                const CagdRType Local[][3],
                                int LocalCount,
                                int LDefOffset);
static void OptMatchBodyBoundaries(OptWorkflowStateStruct *Opt,
                                   const CagdRType Local[][3],
                                   int LocalCount,
                                   int LDefOffset);
static void CreateLegacyDataFile(OptWorkflowStateStruct *Opt);
static CagdRType OptDesignMoveFraction(int Iter);
static CagdRType EnforceExactVolumeConstraint(OptWorkflowStateStruct *Opt,
                                                const CagdRType *Prev,
                                                int Count,
                                                int Iter);
static void InitOrResetResultsArray(OptWorkflowStateStruct *Opt);
static int OptLoadBodies(OptWorkflowStateStruct *Opt);
static int OptParseCounts(const char *Text, int **Values, int *Count);
static void OptWriteCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt);
static int OptReadCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt);
static int OptValidateMesh(OptWorkflowStateStruct *Opt);
static int OptBodyForTile(const OptWorkflowStateStruct *Opt, int Tile);
static void OptTagTile(OptWorkflowStateStruct *Opt, IritPrsrObjectStruct *Obj, int Tile, int *GeomSub);
static int OptCollectMesh(OptWorkflowStateStruct *Opt, IritPrsrObjectStruct *Obj,
                          IritPrsrObjectStruct *Out, int SaveOwners);
static IritPrsrObjectStruct *PostProcessTileBodies(IritPrsrObjectStruct *Tile,
                                      IritUserMicroPostProcessTileCBStruct *CB);

/* ---------- legacy IRIT geometry generation ---------- */
static float **ReadLegacyCsvToFloatArray(OptWorkflowStateStruct *Opt,
                                         const char *FileName,
                                         int *OutRows,
                                         int *OutCols);
static void FreeLegacyCsv(float **Data, int Rows);
static int LoadLegacyGeometryData(OptWorkflowStateStruct *Opt, const char *FileName);
static void UnloadLegacyGeometryData(OptWorkflowStateStruct *Opt);
static IritPrsrObjectStruct *PreProcessTile(
                                  IritPrsrObjectStruct* Tile,
                                  IritUserMicroPreProcessTileCBStruct* CBData);
static int PreProcessTile1FaceParam(
                    IritUserMicroTileBndryPrmStruct* BPrm,
                    CagdRType OuterRadius,
                    CagdBType Circ);
static int GenerateMicroStructures(OptWorkflowStateStruct *Opt, int mode);

#include "multibody_support.inc"

/***************************************************************************
*                             Helper Utilities.                            *
***************************************************************************/

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptFatal(const char *Msg)
{
    fprintf(stderr, "ERROR: %s\n", Msg);
    exit(EXIT_FAILURE);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptFatalF(const char *Fmt, const char *Arg)
{
    char Buffer[OPT_MAX_LINE];
    snprintf(Buffer, sizeof(Buffer), Fmt, Arg);
    OptFatal(Buffer);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static char *OptTrim(char *Str)
{
    char *End;

    if (Str == NULL)
        return NULL;

    while (*Str != '\0' && isspace((unsigned char)*Str))
        ++Str;

    End = Str + strlen(Str);
    while (End > Str && isspace((unsigned char)End[-1]))
        *--End = '\0';

    return Str;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int OptIsBlank(const char *Str)
{
    if (Str == NULL)
        return 1;

    while (*Str != '\0') {
        if (!isspace((unsigned char)*Str))
            return 0;
        ++Str;
    }
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptStripUtf8Bom(char *Line)
{
    unsigned char *Ptr = (unsigned char*)Line;
    if (Ptr != NULL && Ptr[0] == 0xEF && Ptr[1] == 0xBB && Ptr[2] == 0xBF) {
        memmove(Line, Line + 3, strlen(Line + 3) + 1);
    }
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int OptParseInt(const char *Text, int *Value)
{
    char *EndPtr;
    long V;

    if (OptIsBlank(Text))
        return 0;

    V = strtol(Text, &EndPtr, 10);
    while (EndPtr != NULL && *EndPtr != '\0' && isspace((unsigned char)*EndPtr))
        ++EndPtr;

    if (EndPtr == Text || (EndPtr != NULL && *EndPtr != '\0'))
        return 0;

    *Value = (int)V;
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int OptParseReal(const char *Text, CagdRType *Value)
{
    char *EndPtr;
    double V;

    if (OptIsBlank(Text))
        return 0;

    V = strtod(Text, &EndPtr);
    while (EndPtr != NULL && *EndPtr != '\0' && isspace((unsigned char)*EndPtr))
        ++EndPtr;

    if (EndPtr == Text || (EndPtr != NULL && *EndPtr != '\0'))
        return 0;

    *Value = (CagdRType)V;
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptPathJoin(char *Dst, size_t DstSize, const char *A, const char *B)
{
    size_t LenA;

    if (OptIsBlank(A)) {
        snprintf(Dst, DstSize, "%s", B != NULL ? B : "");
        return;
    }
    if (OptIsBlank(B)) {
        snprintf(Dst, DstSize, "%s", A);
        return;
    }

    LenA = strlen(A);
    if (A[LenA - 1] == '\\' || A[LenA - 1] == '/')
        snprintf(Dst, DstSize, "%s%s", A, B);
    else
        snprintf(Dst, DstSize, "%s%c%s", A, OPT_PATH_SEP, B);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static const char *OptPathFileName(const char *Path)
{
    const char *P1 = strrchr(Path, '\\');
    const char *P2 = strrchr(Path, '/');
    const char *P = P1 > P2 ? P1 : P2;
    return P != NULL ? P + 1 : Path;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptMakeIndexedFileName(const char *FileName,
                                   int Index,
                                   char *Out,
                                   size_t OutSize)
{
    const char *Slash1;
    const char *Slash2;
    const char *Slash;
    const char *Dot;
    size_t PrefixLen;

    if (FileName == NULL || Out == NULL || OutSize == 0)
        return;

    Slash1 = strrchr(FileName, '\\');
    Slash2 = strrchr(FileName, '/');
    Slash = Slash1 > Slash2 ? Slash1 : Slash2;
    Dot = strrchr(FileName, '.');

    if (Dot == NULL || (Slash != NULL && Dot < Slash)) {
        snprintf(Out, OutSize, "%s_%d", FileName, Index);
        return;
    }

    PrefixLen = (size_t)(Dot - FileName);
    if (PrefixLen >= OutSize)
        PrefixLen = OutSize - 1;

    memcpy(Out, FileName, PrefixLen);
    Out[PrefixLen] = '\0';
    snprintf(Out + PrefixLen, OutSize - PrefixLen, "_%d%s", Index, Dot);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int OptFileExists(const char *FileName)
{
    FILE *F = fopen(FileName, "rb");
    if (F != NULL) {
        fclose(F);
        return 1;
    }
    return 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptFormatNowHMS(char Buffer[16])
{
    time_t Now = time(NULL);
    struct tm TmBuf;
    struct tm* TmPtr = NULL;

#ifdef _WIN32
    if (localtime_s(&TmBuf, &Now) == 0)
        TmPtr = &TmBuf;
#else
    TmPtr = localtime_r(&Now, &TmBuf);
#endif

    if (TmPtr == NULL)
        snprintf(Buffer, 16, "00:00:00");
    else
        strftime(Buffer, 16, "%H:%M:%S", TmPtr);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptMakeRunFolderName(char Buffer[OPT_MAX_NAME])
{
    time_t Now = time(NULL);
    struct tm TmBuf;
    struct tm* TmPtr = NULL;

#ifdef _WIN32
    if (localtime_s(&TmBuf, &Now) == 0)
        TmPtr = &TmBuf;
#else
    TmPtr = localtime_r(&Now, &TmBuf);
#endif

    if (TmPtr == NULL)
        snprintf(Buffer, OPT_MAX_NAME, "run_unknown_time");
    else
        strftime(Buffer, OPT_MAX_NAME, "%Y%m%d_%H%M%S", TmPtr);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int OptWriteLine(const char *FileName, const char *Content, int Append)
{
    FILE *F = fopen(FileName, Append ? "a" : "w");
    if (F == NULL)
        return 0;

    fprintf(F, "%s\n", Content);
    fclose(F);
    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Allocates a zeroed memory block using the IRIT allocation mechanism.   *
*                                                                          *
* PARAMETERS:                                                              *
*   Count: Number of elements to allocate.                                 *
*   Size:  Size in bytes of one element.                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   void *: Allocated and zeroed memory block, or NULL on failure.         *
***************************************************************************/
static void *OptAllocZero(size_t Count, size_t Size)
{
    size_t Bytes;
    void *Data;

    if (Count == 0 || Size == 0)
        Bytes = 1;
    else {
        Bytes = Count * Size;
        if (Bytes / Count != Size)
            return NULL;
    }

    Data = IritMalloc(Bytes);
    if (Data != NULL)
        IRIT_ZAP_MEM(Data, Bytes);

    return Data;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reallocates a memory block by allocating a new block, copying the      *
* existing content, and freeing the old block.                             *
*                                                                          *
* PARAMETERS:                                                              *
*   OldData:  Existing block to replace.                                   *
*   OldBytes: Size in bytes of the old block.                              *
*   NewBytes: Size in bytes of the new block.                              *
*                                                                          *
* RETURN VALUE:                                                            *
*   void *: Reallocated memory block, or NULL on failure.                  *
***************************************************************************/
static void *OptReallocCopy(void *OldData, size_t OldBytes, size_t NewBytes)
{
    char *NewData;

    if (NewBytes == 0)
        NewBytes = 1;

    NewData = (char *)IritMalloc(NewBytes);
    if (NewData == NULL)
        return NULL;

    if (OldData != NULL && OldBytes > 0)
        memcpy(NewData, OldData, IRIT_MIN(OldBytes, NewBytes));
    if (NewBytes > OldBytes)
        memset(NewData + OldBytes, 0, NewBytes - OldBytes);
    if (OldData != NULL)
        IritFree(OldData);

    return NewData;
}

/***************************************************************************
*                            Input File Parsing.                           *
***************************************************************************/

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptSetDefaultRunInputs(OptRunInputStruct *Run)
{
    IRIT_ZAP_MEM(Run, sizeof(OptRunInputStruct));

    Run -> SimulationType = 1;
    Run -> LatticeType = -1;
    Run -> LatticeCountU = -1;
    Run -> LatticeCountV = -1;
    Run -> LatticeCountW = -1;
    Run -> InitialThickness = OPT_INVALID_REAL;
    Run -> MeshType = -1;
    Run -> MeshSample = -1;
    Run -> MaxIterations = -1;
    Run -> ConvergenceIterations = 3;
    Run -> ThicknessChangeTolerance = 5.0e-3 * OPT_MAX_TILE_THICKNESS;
    //Run -> TargetDeflection = OPT_INVALID_REAL;
    //Run -> TargetMass = OPT_INVALID_REAL;
    //Run -> DeflectionTolerance = OPT_INVALID_REAL;
    //Run -> MassTolerance = OPT_INVALID_REAL;
    Run -> TargetVolumeFraction = OPT_INVALID_REAL;
    Run -> TargetVolumeFractionTol = OPT_INVALID_REAL;
    Run -> OptimizationAlgorithm = 1;
    Run -> AlgoType = 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void OptSetDefaultConfigInputs(OptConfigInputStruct *Cfg)
{
    IRIT_ZAP_MEM(Cfg, sizeof(OptConfigInputStruct));

    Cfg -> Scale = 0.015;
    Cfg -> Density = 2700.0;
    Cfg -> YoungsModulus = 6.9e10;
    Cfg -> PoissonRatio = 0.33;
    Cfg -> SedValue = 0.15;
    Cfg -> YieldStress = 258.0;
    Cfg -> MinimumThickness = 1.0e-5;
    Cfg -> SafetyFactorSed = 1.1;
    Cfg -> SafetyFactorStress = 1.25;
    Cfg -> NodesPerElement = 8;
    Cfg -> DefaultL = 7;
    Cfg -> DefaultLMin = 6;
    Cfg -> DefaultLDef = 0;
    Cfg -> ElementsPerTile = -1;
    Cfg -> MeshSampleTotal = -1;
    Cfg -> Tile0OwnerMapping = 0;
    Cfg -> RemoveThinTiles = FALSE;
    OptPathJoin(Cfg->VolumeFile,
        sizeof(Cfg->VolumeFile),
        Cfg->WorkDir,
        "VOLUME.txt");
}

static int ParseKeyValueFile(const char *FileName,
    int (*Handler)(const char *Key, const char *Value, void *Ctx),
    void *Ctx)
{
    FILE *F;
    char Line[OPT_MAX_LINE];
    int LineNum = 0;

    F = fopen(FileName, "r");
    if (F == NULL) {
        fprintf(stderr, "Failed to open input file: %s\n", FileName);
        return 0;
    }

    while (fgets(Line, sizeof(Line), F) != NULL) {
        char *Eq;
        char *Key;
        char *Value;

        ++LineNum;
        OptStripUtf8Bom(Line);

        Key = OptTrim(Line);
        if (*Key == '\0' || *Key == '#')
            continue;

        Eq = strchr(Key, '=');
        if (Eq == NULL)
            continue;

        *Eq = '\0';
        Value = OptTrim(Eq + 1);
        Key = OptTrim(Key);

        if (!Handler(Key, Value, Ctx)) {
            fprintf(stderr, "Invalid key or value in %s at line %d: %s\n",
                FileName, LineNum, Key);
            fclose(F);
            return 0;
        }
    }

    fclose(F);
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int RunInputHandler(const char *Key, const char *Value, void *Ctx)
{
    OptRunInputStruct *Run = (OptRunInputStruct *)Ctx;

    if (strcmp(Key, "macro_shape_path") == 0) {
        snprintf(Run -> MacroShapePath, sizeof(Run -> MacroShapePath), "%s", Value);
    }
    else if (strcmp(Key, "simulation_type") == 0) {
        return OptParseInt(Value, &Run -> SimulationType);
    }
    else if (strcmp(Key, "lattice_type") == 0) {
        return OptParseInt(Value, &Run -> LatticeType);
    }
    else if (strcmp(Key, "lattice_count_u") == 0) {
        if (!OptParseCounts(Value, &Run->BodyCounts[0], &Run->BodyCountLengths[0]))
            return 0;
        Run->LatticeCountU = Run->BodyCounts[0][0];
        return 1;
    }
    else if (strcmp(Key, "lattice_count_v") == 0) {
        if (!OptParseCounts(Value, &Run->BodyCounts[1], &Run->BodyCountLengths[1]))
            return 0;
        Run->LatticeCountV = Run->BodyCounts[1][0];
        return 1;
    }
    else if ((strcmp(Key, "lattice_count_w") == 0 || strcmp(Key, "lattice count_w") == 0)) {
        if (!OptParseCounts(Value, &Run->BodyCounts[2], &Run->BodyCountLengths[2]))
            return 0;
        Run->LatticeCountW = Run->BodyCounts[2][0];
        return 1;
    }
    else if (strcmp(Key, "initial_thickness") == 0) {
        return OptParseReal(Value, &Run -> InitialThickness);
    }
    else if (strcmp(Key, "mesh_type") == 0) {
        return OptParseInt(Value, &Run -> MeshType);
    }
    else if (strcmp(Key, "mesh_sample") == 0) {
        return OptParseInt(Value, &Run -> MeshSample);
    }
    else if (strcmp(Key, "max_iterations") == 0) {
        return OptParseInt(Value, &Run -> MaxIterations);
    }
    else if (strcmp(Key, "convergence_iterations") == 0) {
        return OptParseInt(Value, &Run->ConvergenceIterations);
    }
    else if (strcmp(Key, "thickness_change_tolerance") == 0) {
        return OptParseReal(Value, &Run->ThicknessChangeTolerance);
    }
    /*
    else if (strcmp(Key, "target_deflection") == 0) {
        return OptParseReal(Value, &Run -> TargetDeflection);
    }
    else if (strcmp(Key, "target_mass") == 0) {
        if (!OptIsBlank(Value))
            return OptParseReal(Value, &Run -> TargetMass);
    }
    else if (strcmp(Key, "deflection_tolerance") == 0) {
        if (!OptIsBlank(Value))
            return OptParseReal(Value, &Run -> DeflectionTolerance);
    }
    else if (strcmp(Key, "mass_tolerance") == 0) {
        if (!OptIsBlank(Value))
            return OptParseReal(Value, &Run -> MassTolerance);
    }
    */
    else if (strcmp(Key, "target_volume_fraction") == 0) {
        if (!OptIsBlank(Value))
            return OptParseReal(Value, &Run->TargetVolumeFraction);
    }
    else if (strcmp(Key, "target_volume_fraction_tol") == 0) {
        if (!OptIsBlank(Value))
            return OptParseReal(Value, &Run->TargetVolumeFractionTol);
    }
    else if (strcmp(Key, "optimization_algorithm") == 0) {
        if (!OptIsBlank(Value))
            return OptParseInt(Value, &Run -> OptimizationAlgorithm);
    }
    else if (strcmp(Key, "run_folder_name") == 0) {
        snprintf(Run -> RunFolderName, sizeof(Run -> RunFolderName), "%s", Value);
    }
    else {
        return 0;
    }

    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int ConfigInputHandler(const char *Key, const char *Value, void *Ctx)
{
    OptConfigInputStruct *Cfg = (OptConfigInputStruct *)Ctx;

#define COPY_CFG_STR(FieldName, KeyName) \
    if (strcmp(Key, KeyName) == 0) { \
        snprintf(Cfg->FieldName, sizeof(Cfg->FieldName), "%s", Value); \
        return 1; \
    }

    COPY_CFG_STR(BaseFolder, "base_folder")
    COPY_CFG_STR(WorkDir, "work_dir")
    COPY_CFG_STR(AnsysExe, "ansys_exe")
    COPY_CFG_STR(AnsysProject, "ansys_project")
    COPY_CFG_STR(AnsysJournal, "ansys_journal")
    COPY_CFG_STR(IritExe, "irit_exe")
    COPY_CFG_STR(InpFile, "inp_file")
    COPY_CFG_STR(ItdOutputFile, "itd_output_file")
    COPY_CFG_STR(StressFile, "stress_file")
    COPY_CFG_STR(DeflectionFile, "deflection_file")
    COPY_CFG_STR(SedFile, "sed_file")
    COPY_CFG_STR(VolumeFile, "volume_file")
    COPY_CFG_STR(DataNewFile, "data_new_file")
    COPY_CFG_STR(DataPrevFile, "data_prev_file")
    COPY_CFG_STR(DesignParametersFile, "design_parameters_file")
    COPY_CFG_STR(SummaryFile, "summary_file")

#undef COPY_CFG_STR

    if (strcmp(Key, "scale") == 0)
        return OptParseReal(Value, &Cfg -> Scale);
    if (strcmp(Key, "density") == 0)
        return OptParseReal(Value, &Cfg -> Density);
    if (strcmp(Key, "youngs_modulus") == 0)
        return OptParseReal(Value, &Cfg -> YoungsModulus);
    if (strcmp(Key, "poisson_ratio") == 0)
        return OptParseReal(Value, &Cfg -> PoissonRatio);
    if (strcmp(Key, "sed_value") == 0)
        return OptParseReal(Value, &Cfg -> SedValue);
    if (strcmp(Key, "yield_stress") == 0)
        return OptParseReal(Value, &Cfg -> YieldStress);
    if (strcmp(Key, "minimum_thickness") == 0)
        return OptParseReal(Value, &Cfg -> MinimumThickness);
    if (strcmp(Key, "safety_factor_sed") == 0)
        return OptParseReal(Value, &Cfg -> SafetyFactorSed);
    if (strcmp(Key, "safety_factor_stress") == 0)
        return OptParseReal(Value, &Cfg -> SafetyFactorStress);
    if (strcmp(Key, "nodes_per_element") == 0)
        return OptParseInt(Value, &Cfg -> NodesPerElement);
    if (strcmp(Key, "default_L") == 0)
        return OptParseInt(Value, &Cfg -> DefaultL);
    if (strcmp(Key, "default_L_min") == 0)
        return OptParseInt(Value, &Cfg -> DefaultLMin);
    if (strcmp(Key, "default_L_def") == 0)
        return OptParseInt(Value, &Cfg -> DefaultLDef);
    if (strcmp(Key, "elements_per_tile") == 0) {
        if (!OptIsBlank(Value))
            return OptParseInt(Value, &Cfg -> ElementsPerTile);
        return 1;
    }
    if (strcmp(Key, "mesh_sample_total") == 0) {
        if (!OptIsBlank(Value))
            return OptParseInt(Value, &Cfg -> MeshSampleTotal);
        return 1;
    }
    if (strcmp(Key, "tile0_owner_mapping") == 0)
        return OptParseInt(Value, &Cfg -> Tile0OwnerMapping);

    if (strcmp(Key, "remove_thin_tiles") == 0) {
        Cfg -> RemoveThinTiles = OptParseBool(Value);
        return 1;
    }
    else if (OPT_STRICMP(Key, "volume_file") == 0) {
        snprintf(Cfg->VolumeFile, sizeof(Cfg->VolumeFile), "%s", Value);
    }

    return 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static CagdBType OptParseBool(const char *Value)
{
#ifdef _WIN32
    if (_stricmp(Value, "true") == 0 || _stricmp(Value, "yes") == 0 || strcmp(Value, "1") == 0)
        return TRUE;
    if (_stricmp(Value, "false") == 0 || _stricmp(Value, "no") == 0 || strcmp(Value, "0") == 0)
        return FALSE;
#else
    if (strcasecmp(Value, "true") == 0 || strcasecmp(Value, "yes") == 0 || strcmp(Value, "1") == 0)
        return TRUE;
    if (strcasecmp(Value, "false") == 0 || strcasecmp(Value, "no") == 0 || strcmp(Value, "0") == 0)
        return FALSE;
#endif

    fprintf(stderr, "Invalid boolean value: %s\n", Value);
    return FALSE;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int LoadRunInputs(const char *FileName, OptRunInputStruct *Run)
{
    OptSetDefaultRunInputs(Run);
    return ParseKeyValueFile(FileName, RunInputHandler, Run);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int LoadConfigInputs(const char *FileName, OptConfigInputStruct *Cfg)
{
    OptSetDefaultConfigInputs(Cfg);
    return ParseKeyValueFile(FileName, ConfigInputHandler, Cfg);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int ValidateInputs(const OptRunInputStruct *Run, const OptConfigInputStruct *Cfg)
{
    if (OptIsBlank(Run -> MacroShapePath)) {
        fprintf(stderr, "Missing macro_shape_path in opt_run.txt\n");
        return 0;
    }
    if (Run -> SimulationType < 0 || Run -> SimulationType > 2) {
        fprintf(stderr, "simulation_type must be 0 (emulation), 1 (ANSYS), or 2 (CalculiX).\n");
        return 0;
    }
    if (Run -> LatticeType < 0 || Run->LatticeType > 2 ||
        Run -> LatticeCountU <= 0 ||
        Run -> LatticeCountV <= 0 ||
        Run -> LatticeCountW <= 0) {
        fprintf(stderr, "Invalid lattice parameters in opt_run.txt\n");
        return 0;
    }
    if (Run->BodyCountLengths[0] != Run->BodyCountLengths[1] ||
        Run->BodyCountLengths[0] != Run->BodyCountLengths[2]) {
        fprintf(stderr, "U, V and W count lists must have equal lengths.\n");
        return 0;
    }
    if ((double)(Run->MeshSample - 1) * (Run->MeshSample - 1) * (Run->MeshSample - 1) > INT_MAX) {
        fprintf(stderr, "mesh_sample is too large.\n"); return 0;
    }
    if (Cfg->Scale <= 0.0 || !isfinite(Cfg->Scale)) return 0;
    if (Run -> InitialThickness <= 0.0 || Run -> MeshSample <= 1 || Run -> MaxIterations <= 0) {
        fprintf(stderr,
                "Invalid initial thickness / mesh sample / max iterations "
                "in opt_run.txt\n");
        return 0;
    }
    if (Run->ConvergenceIterations <= 0) {
        fprintf(stderr, "convergence_iterations must be >= 1 in opt_run.txt\n");
        return 0;
    }
    if (Run->ThicknessChangeTolerance <= 0.0 ||
        !isfinite(Run->ThicknessChangeTolerance)) {
        fprintf(stderr, "thickness_change_tolerance must be > 0 in opt_run.txt\n");
        return 0;
    }
    /*
    if (Run -> TargetDeflection == OPT_INVALID_REAL || Run -> TargetMass == OPT_INVALID_REAL ||
        Run -> DeflectionTolerance == OPT_INVALID_REAL ||
        Run -> MassTolerance == OPT_INVALID_REAL) {
        fprintf(stderr,
                "target_deflection, target_mass, deflection_tolerance "
                "and mass_tolerance must all be filled in opt_run.txt\n");
        return 0;
    }
    */
    if (OptIsBlank(Cfg -> BaseFolder) || OptIsBlank(Cfg -> WorkDir) ||
        OptIsBlank(Cfg -> InpFile) ||
        OptIsBlank(Cfg -> ItdOutputFile) || OptIsBlank(Cfg -> StressFile) ||
        OptIsBlank(Cfg -> DeflectionFile) || OptIsBlank(Cfg -> SedFile) ||
        OptIsBlank(Cfg -> DataNewFile) || OptIsBlank(Cfg -> DataPrevFile) ||
        OptIsBlank(Cfg -> DesignParametersFile) || OptIsBlank(Cfg -> SummaryFile)) {
        fprintf(stderr, "Missing required paths in opt_config.txt\n");
        return 0;
    }

    if (Run -> SimulationType == 1 &&
        (OptIsBlank(Cfg -> AnsysExe) || OptIsBlank(Cfg -> AnsysProject) ||
         OptIsBlank(Cfg -> AnsysJournal))) {
        fprintf(stderr, "ANSYS paths are required when simulation_type is 1.\n");
        return 0;
    }

    if (Cfg -> NodesPerElement <= 0 || Cfg -> NodesPerElement > OPT_MAX_NODE_PER_ELEM) {
        fprintf(stderr, "nodes_per_element in opt_config.txt is invalid for this program\n");
        return 0;
    }
    if (Cfg->Tile0OwnerMapping < 0 || Cfg->Tile0OwnerMapping > 1) {
        fprintf(stderr,
                "tile0_owner_mapping must be 0 (legacy) or 1 (centroid-verified).\n");
        return 0;
    }

    return 1;
}

/***************************************************************************
 *                      WORKFLOW STATE INITIALIZATION
 ***************************************************************************/

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void ResolveLatticeLayout(OptWorkflowStateStruct *Opt)
{
    switch (Opt->Run.LatticeType) {
    case 0: /* Face (3D Cross 1): 6 face arms. */
    case 1: /* Face (3D Cross 2): 6 face arms. */
        Opt->L = 7;
        Opt->LMin = 6;
        Opt->LDef = 0;
        break;
    case 2: /* Diagonal: 8 vertex arms. */
    case 3: /* Combined: 6 face + 8 vertex arms. */
        Opt->L = 24;
        Opt->LMin = 8;
        Opt->LDef = 8;
        break;
    default:
        Opt->L = Opt->Config.DefaultL;
        Opt->LMin = Opt->Config.DefaultLMin;
        Opt->LDef = Opt->Config.DefaultLDef;
        break;
    }

    Opt->MeshSample = (Opt->Run.MeshSample - 1) *
        (Opt->Run.MeshSample - 1) *
        (Opt->Run.MeshSample - 1);

    Opt->ElementsPerTile = Opt->L * Opt->MeshSample;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void InitOrResetResultsArray(OptWorkflowStateStruct *Opt)
{
    if (Opt->Results != NULL) {
        IritFree(Opt->Results);
        Opt->Results = NULL;
    }

    Opt->Results = (CagdRType *) OptAllocZero((size_t) Opt->Run.MaxIterations * 3,
                                              sizeof(CagdRType));
    if (Opt->Results == NULL)
        OptFatal("Failed to allocate results array");
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void InitializeWorkflowState(OptWorkflowStateStruct *Opt)
{
    IRIT_ZAP_MEM(&Opt->Nodes, sizeof(Opt->Nodes));
    IRIT_ZAP_MEM(&Opt->Elements, sizeof(Opt->Elements));
    IRIT_ZAP_MEM(&Opt->StressData, sizeof(Opt->StressData));
    IRIT_ZAP_MEM(&Opt->DeflectionData, sizeof(Opt->DeflectionData));
    IRIT_ZAP_MEM(&Opt->SedData, sizeof(Opt->SedData));

    Opt->LastThicknessChange = 1.0e30;
    Opt->ConsecutiveSmallChanges = 0;

    ResolveLatticeLayout(Opt);
    if (!OptLoadBodies(Opt))
        OptFatal("Cannot initialize body layout");

    /*
     * Build the persistent tile/subtile hierarchy ONCE.
     * Thickness starts at InitialThickness and is changed later only by
     * DesignUpdate_*() functions.  ANSYS-derived fields are refreshed after
     * every solve without rebuilding this structure.
     */
    if (!InitTileStructure(Opt))
        OptFatal("Cannot initialize persistent TileStructure");

    InitOrResetResultsArray(Opt);

    if (OptIsBlank(Opt->Run.RunFolderName))
        OptMakeRunFolderName(Opt->FolderName);
    else
        snprintf(Opt->FolderName, sizeof(Opt->FolderName), "%s", Opt->Run.RunFolderName);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void FreeWorkflowState(OptWorkflowStateStruct *Opt)
{
    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);
    ScalarFieldReset(&Opt->StressData);
    ScalarFieldReset(&Opt->DeflectionData);
    ScalarFieldReset(&Opt->SedData);
    DestroyTileStructure(Opt);
    {
        int B;
        for (B = 0; B < Opt->BodyCount; ++B)
            IritTrivTVFree(Opt->Bodies[B].Map);
        if (Opt->Bodies) IritFree(Opt->Bodies);
        if (Opt->MeshOwners) IritFree(Opt->MeshOwners);
        for (B = 0; B < 3; ++B)
            if (Opt->Run.BodyCounts[B]) IritFree(Opt->Run.BodyCounts[B]);
        ScalarFieldReset(&Opt->VolumeData);
    }

    if (Opt->Results != NULL) {
        IritFree(Opt->Results);
        Opt->Results = NULL;
    }
}

/***************************************************************************
 *                         DYNAMIC CONTAINER HELPERS
 ***************************************************************************/

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void NodeArrayReset(OptNodeArrayStruct *Arr)
{
    if (Arr -> Items != NULL)
        IritFree(Arr -> Items);
    Arr -> Items = NULL;
    Arr -> Count = 0;
    Arr -> Capacity = 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void ElementArrayReset(OptElementArrayStruct *Arr)
{
    if (Arr -> Items != NULL)
        IritFree(Arr -> Items);
    Arr -> Items = NULL;
    Arr -> Count = 0;
    Arr -> Capacity = 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int NodeArrayPush(OptNodeArrayStruct *Arr, const OptNodeStruct *Node)
{
    if (Arr -> Count >= Arr -> Capacity) {
        int NewCap = Arr -> Capacity == 0 ? 1024 : Arr -> Capacity * 2;
        OptNodeStruct *NewItems =
            (OptNodeStruct *) OptReallocCopy(Arr -> Items,
                                             (size_t) Arr -> Capacity *
                                                 sizeof(OptNodeStruct),
                                             (size_t) NewCap *
                                                 sizeof(OptNodeStruct));
        if (NewItems == NULL)
            return 0;
        Arr -> Items = NewItems;
        Arr -> Capacity = NewCap;
    }

    Arr -> Items[Arr -> Count++] = *Node;
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int ElementArrayPush(OptElementArrayStruct *Arr, const OptElementStruct *Elem)
{
    if (Arr -> Count >= Arr -> Capacity) {
        int NewCap = Arr -> Capacity == 0 ? 1024 : Arr -> Capacity * 2;
        OptElementStruct *NewItems =
            (OptElementStruct *) OptReallocCopy(Arr -> Items,
                                                (size_t) Arr -> Capacity *
                                                    sizeof(OptElementStruct),
                                                (size_t) NewCap *
                                                    sizeof(OptElementStruct));
        if (NewItems == NULL)
            return 0;
        Arr -> Items = NewItems;
        Arr -> Capacity = NewCap;
    }

    Arr -> Items[Arr -> Count++] = *Elem;
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void ScalarFieldReset(OptScalarFieldStruct *Field)
{
    if (Field -> Values != NULL)
        IritFree(Field -> Values);
    if (Field -> HasValue != NULL)
        IritFree(Field -> HasValue);

    Field -> Values = NULL;
    Field -> HasValue = NULL;
    Field -> Capacity = 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int ScalarFieldEnsureCapacity(OptScalarFieldStruct *Field, int RequiredId)
{
    int NewCap;
    CagdRType *NewValues;
    CagdBType* NewHasValue;
    int OldCap;

    if (RequiredId < Field -> Capacity)
        return 1;

    OldCap = Field -> Capacity;
    NewCap = Field -> Capacity == 0 ? 1024 : Field -> Capacity;
    while (NewCap <= RequiredId)
        NewCap *= 2;

    NewValues = (CagdRType *) OptReallocCopy(Field -> Values,
                                             (size_t) OldCap *
                                                 sizeof(CagdRType),
                                             (size_t) NewCap *
                                                 sizeof(CagdRType));
    if (NewValues == NULL)
        return 0;
    Field -> Values = NewValues;

    NewHasValue = (CagdBType *) OptReallocCopy(Field -> HasValue,
                                              (size_t) OldCap *
                                                  sizeof(CagdBType),
                                              (size_t) NewCap *
                                                  sizeof(CagdBType));
    if (NewHasValue == NULL)
        return 0;
    Field -> HasValue = NewHasValue;
    Field -> Capacity = NewCap;
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void ScalarFieldSet(OptScalarFieldStruct *Field, int Id, CagdRType Value)
{
    if (!ScalarFieldEnsureCapacity(Field, Id))
        OptFatal("Failed to grow scalar field");
    Field -> Values[Id] = Value;
    Field -> HasValue[Id] = TRUE;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static CagdBType ScalarFieldTryGet(const OptScalarFieldStruct *Field, int Id, CagdRType *Value)
{
    if (Id < 0 || Id >= Field -> Capacity || !Field -> HasValue[Id])
        return FALSE;

    *Value = Field -> Values[Id];
    return TRUE;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int SubTileAddElementPtr(OptSubTileStruct *SubTile, OptElementStruct *Elem)
{
    if (SubTile -> ElementCount >= SubTile -> ElementCapacity) {
        int NewCap = SubTile -> ElementCapacity == 0 ? 16 : SubTile -> ElementCapacity * 2;
        OptElementStruct **NewItems =
            (OptElementStruct **) OptReallocCopy(SubTile -> Elements,
                                                 (size_t) SubTile -> ElementCapacity *
                                                     sizeof(OptElementStruct *),
                                                 (size_t) NewCap *
                                                     sizeof(OptElementStruct *));
        if (NewItems == NULL)
            return 0;
        SubTile -> Elements = NewItems;
        SubTile -> ElementCapacity = NewCap;
    }

    SubTile -> Elements[SubTile -> ElementCount++] = Elem;
    return 1;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void DestroyTileStructure(OptWorkflowStateStruct *Opt)
{
    int TileCount, T, S;

    if (Opt->TileStructure == NULL)
        return;

    TileCount = Opt->TotalTiles;
    for (T = 0; T < TileCount; ++T) {
        if (Opt->TileStructure[T].SubTiles != NULL) {
            for (S = 0; S < Opt->TileStructure[T].SubTileCount; ++S) {
                if (Opt->TileStructure[T].SubTiles[S].Elements != NULL)
                    IritFree(Opt->TileStructure[T].SubTiles[S].Elements);
            }
            IritFree(Opt->TileStructure[T].SubTiles);
        }
    }

    IritFree(Opt->TileStructure);
    Opt->TileStructure = NULL;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int InitTileStructure(OptWorkflowStateStruct *Opt)
{
    int TileCount;
    int T;

    /*
     * This is a true one-time initialization.  Re-entering it would hide
     * a workflow error, so an existing structure is simply left untouched.
     */
    if (Opt->TileStructure != NULL)
        return 1;

    TileCount = Opt->TotalTiles;
    if (TileCount <= 0 || Opt->L <= 0)
        return 0;

    Opt->TileStructure =
        (OptTileStruct *) OptAllocZero((size_t) TileCount, sizeof(OptTileStruct));
    if (Opt->TileStructure == NULL)
        return 0;

    for (T = 0; T < TileCount; ++T) {
        OptTileStruct *Tile = &Opt->TileStructure[T];

        Tile->SubTileCount = Opt->L;
        Tile->SubTiles =
            (OptSubTileStruct *) OptAllocZero((size_t) Opt->L,
                                              sizeof(OptSubTileStruct));
        if (Tile->SubTiles == NULL) {
            DestroyTileStructure(Opt);
            return 0;
        }

        /* Tile identity/topology does not change during optimization. */
        snprintf(Tile->Name,
                 sizeof(Tile->Name),
                 "Body_%d_Tile_%d",
                 OptBodyForTile(Opt, T) + 1,
                 T);

        /*
         * Initial design state.  From this point onward Thickness is changed
         * only by DesignUpdate_*() (and explicit post-update corrections such
         * as neighbour matching / volume correction).
         */
        {
            int L;
            for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L)
                Tile->SubTiles[L].Thickness = Opt->Run.InitialThickness;
        }
    }

    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Clears only data that belongs to the previous ANSYS solution while     *
*   preserving the allocated tile/subtile hierarchy, element-pointer       *
*   buffers, capacities, names, and thickness design state.                 *
***************************************************************************/
static void ClearTileAnalysisData(OptWorkflowStateStruct *Opt)
{
    int T, S;

    if (Opt == NULL || Opt->TileStructure == NULL)
        return;

    for (T = 0; T < Opt->TotalTiles; ++T) {
        OptTileStruct *Tile = &Opt->TileStructure[T];

        Tile->AvgStress = 0.0;
        Tile->MinStress = 0.0;
        Tile->MaxStress = 0.0;
        Tile->Volume = 0.0;
        Tile->StrainEnergy = 0.0;
        Tile->ComplianceContribution = 0.0;
        Tile->Gain = 0.0;

        for (S = 0; S < Tile->SubTileCount; ++S) {
            OptSubTileStruct *SubTile = &Tile->SubTiles[S];

            /*
             * The old entries point into the old Opt->Elements array, which
             * ProcessIterationData() has already freed.  Null the entries but
             * KEEP the allocated pointer buffer and its capacity for reuse.
             */
            if (SubTile->Elements != NULL && SubTile->ElementCapacity > 0) {
                IRIT_ZAP_MEM(SubTile->Elements,
                             (size_t) SubTile->ElementCapacity *
                                 sizeof(OptElementStruct *));
            }
            SubTile->ElementCount = 0;

            /* Preserve SubTile->Thickness: it is design state. */
            SubTile->AvgStress = 0.0;
            SubTile->MinStress = 0.0;
            SubTile->MaxStress = 0.0;
            SubTile->SED = 0.0;
            SubTile->MaxSED = 0.0;
            SubTile->Volume = 0.0;
            SubTile->StrainEnergy = 0.0;
            SubTile->ComplianceContribution = 0.0;
            SubTile->Sensitivity = 0.0;
            SubTile->Gain = 0.0;
        }
    }
}

/***************************************************************************
 *                            WINDOWS FILE HELPERS
 ***************************************************************************/

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int EnsureDirectoryRecursive(const char *DirName)
{
#ifdef _WIN32
    char Tmp[OPT_MAX_PATH];
    char *P;

    snprintf(Tmp, sizeof(Tmp), "%s", DirName);
    for (P = Tmp; *P != '\0'; ++P) {
        if (*P == '\\' || *P == '/') {
            char Saved = *P;
            if (P != Tmp && P[-1] != ':') {
                *P = '\0';
                CreateDirectoryA(Tmp, NULL);
                *P = Saved;
            }
        }
    }
    CreateDirectoryA(Tmp, NULL);
    return 1;
#else
    (void)DirName;
    return 1;
#endif
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int CopyDirectoryRecursive(OptWorkflowStateStruct *Opt, const char *SourceDir, const char *TargetDir)
{
#ifdef _WIN32
    char Pattern[OPT_MAX_PATH];
    WIN32_FIND_DATAA Fd;
    HANDLE HFind;

    EnsureDirectoryRecursive(TargetDir);
    snprintf(Pattern, sizeof(Pattern), "%s\\*", SourceDir);

    HFind = FindFirstFileA(Pattern, &Fd);
    if (HFind == INVALID_HANDLE_VALUE)
        return 0;

    do {
        char SrcPath[OPT_MAX_PATH];
        char DstPath[OPT_MAX_PATH];
        const char *Name = Fd.cFileName;
        const char *SkipItd = OptPathFileName(Opt->Config.ItdOutputFile);
        const char *SkipExe = OptPathFileName(Opt->Config.IritExe);

        if (strcmp(Name, ".") == 0 || strcmp(Name, "..") == 0)
            continue;

        OptPathJoin(SrcPath, sizeof(SrcPath), SourceDir, Name);
        OptPathJoin(DstPath, sizeof(DstPath), TargetDir, Name);

        if (Fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            if (!CopyDirectoryRecursive(Opt, SrcPath, DstPath)) {
                FindClose(HFind);
                return 0;
            }
        }
        else {
            if (OPT_STRICMP(Name, SkipItd) == 0 || OPT_STRICMP(Name, SkipExe) == 0)
                continue;
            EnsureDirectoryRecursive(TargetDir);
            if (!CopyFileA(SrcPath, DstPath, FALSE)) {
                FindClose(HFind);
                return 0;
            }
        }
    } while (FindNextFileA(HFind, &Fd));

    FindClose(HFind);
    return 1;
#else
    (void)SourceDir;
    (void)TargetDir;
    return 0;
#endif
}

/***************************************************************************
 *                         PROCESS LAUNCHING (WINDOWS)
 ***************************************************************************/

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int RunCommandLineBlocking(const char *CommandLine, const char *WorkingDir)
{
#ifdef _WIN32
    STARTUPINFOA Si;
    PROCESS_INFORMATION Pi;
    DWORD ExitCode = (DWORD)-1;
    char CmdMutable[OPT_MAX_LINE * 2];

    IRIT_ZAP_MEM(&Si, sizeof(Si));
    IRIT_ZAP_MEM(&Pi, sizeof(Pi));
    Si.cb = sizeof(Si);
    Si.dwFlags = STARTF_USESHOWWINDOW;
    Si.wShowWindow = SW_HIDE;

    snprintf(CmdMutable, sizeof(CmdMutable), "%s", CommandLine);

    if (!CreateProcessA(NULL,
        CmdMutable,
        NULL,
        NULL,
        FALSE,
        CREATE_NO_WINDOW,
        NULL,
        WorkingDir,
        &Si,
        &Pi)) {
        fprintf(stderr, "CreateProcess failed for command:\n%s\n", CommandLine);
        return 0;
    }

    WaitForSingleObject(Pi.hProcess, INFINITE);
    GetExitCodeProcess(Pi.hProcess, &ExitCode);
    CloseHandle(Pi.hThread);
    CloseHandle(Pi.hProcess);

    return ExitCode == 0;
#else
    (void)WorkingDir;
    return system(CommandLine) == 0;
#endif
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Runs a child process synchronously and redirects both stdout and stderr *
*   to a fresh log file.  This is used for long CalculiX solves so progress *
*   and fatal solver messages are retained even though its console is       *
*   hidden.                                                                 *
***************************************************************************/
static int RunCommandLineBlockingToFile(const char *CommandLine,
                                        const char *WorkingDir,
                                        const char *LogFile)
{
#ifdef _WIN32
    SECURITY_ATTRIBUTES Sa;
    STARTUPINFOA Si;
    PROCESS_INFORMATION Pi;
    HANDLE LogHandle = INVALID_HANDLE_VALUE;
    DWORD ExitCode = (DWORD)-1;
    char CmdMutable[OPT_MAX_LINE * 2];

    IRIT_ZAP_MEM(&Sa, sizeof(Sa));
    IRIT_ZAP_MEM(&Si, sizeof(Si));
    IRIT_ZAP_MEM(&Pi, sizeof(Pi));
    Sa.nLength = sizeof(Sa);
    Sa.bInheritHandle = TRUE;

    LogHandle = CreateFileA(LogFile,
                            GENERIC_WRITE,
                            FILE_SHARE_READ | FILE_SHARE_WRITE,
                            &Sa,
                            CREATE_ALWAYS,
                            FILE_ATTRIBUTE_NORMAL,
                            NULL);
    if (LogHandle == INVALID_HANDLE_VALUE) {
        fprintf(stderr, "Cannot create process log: %s\n", LogFile);
        return FALSE;
    }

    Si.cb = sizeof(Si);
    Si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
    Si.wShowWindow = SW_HIDE;
    Si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    Si.hStdOutput = LogHandle;
    Si.hStdError = LogHandle;
    snprintf(CmdMutable, sizeof(CmdMutable), "%s", CommandLine);

    if (!CreateProcessA(NULL,
                        CmdMutable,
                        NULL,
                        NULL,
                        TRUE,
                        CREATE_NO_WINDOW,
                        NULL,
                        WorkingDir,
                        &Si,
                        &Pi)) {
        fprintf(stderr, "CreateProcess failed for command:\n%s\n", CommandLine);
        CloseHandle(LogHandle);
        return FALSE;
    }

    WaitForSingleObject(Pi.hProcess, INFINITE);
    GetExitCodeProcess(Pi.hProcess, &ExitCode);
    CloseHandle(Pi.hThread);
    CloseHandle(Pi.hProcess);
    CloseHandle(LogHandle);
    return ExitCode == 0;
#else
    char Redirected[OPT_MAX_LINE * 2];
    (void)WorkingDir;
    snprintf(Redirected, sizeof(Redirected), "%s > \"%s\" 2>&1",
             CommandLine, LogFile);
    return system(Redirected) == 0;
#endif
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Initializes the legacy thickness CSV file that is consumed by the      *
* geometry generation stage.                                               *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the file was created successfully, FALSE otherwise.       *
***************************************************************************/

static int InitLegacyDataFile(OptWorkflowStateStruct *Opt)
{
    FILE *F;
    int Times;
    int I;

    F = fopen(Opt->Config.DataNewFile, "w");
    if (F == NULL)
        return 0;

    OptWriteCsvHeader(F, Opt);

    Times = Opt->TotalTiles * Opt->LMin;
    for (I = 0; I < Times; ++I)
        fprintf(F, "%d,%.17g,0\n", I + 1, Opt->Run.InitialThickness);

    fclose(F);
    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Executes the full optimization workflow until convergence or until the *
* maximal number of iterations is reached.                                 *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the workflow completed, FALSE on failure.                 *
***************************************************************************/

static int RunWorkflow(OptWorkflowStateStruct *Opt)
{
    int Iter, Finished;
    char TargetBaseFolder[OPT_MAX_PATH];

    OptPathJoin(TargetBaseFolder,
                sizeof(TargetBaseFolder),
                Opt->Config.BaseFolder,
                Opt->FolderName);
    EnsureDirectoryRecursive(TargetBaseFolder);

    if (!InitLegacyDataFile(Opt))
        OptFatalF("Failed to initialize legacy data file: %s", Opt->Config.DataNewFile);

    if (!LatticeBy_VolumeFraction(Opt)) return 0;

    for (Iter = 1; Iter <= Opt->Run.MaxIterations; ++Iter) {
        OptFormatNowHMS(Opt->TimeIritStart);
        if (!RunIRITIntegrated(Opt, 0))
            OptFatal("IRIT generation step failed");

        OptFormatNowHMS(Opt->TimeIritFinish);
        if (!RunSimulation(Opt))
            OptFatal("Simulation execution step failed");

        OptFormatNowHMS(Opt->TimeSimulationFinish);
        Finished = ProcessIterationData(Opt, Iter);
        if (Opt->Run.SimulationType == 0) {
            if (!RunIRITIntegrated(Opt, 0))
                OptFatal("Post-update IRIT generation failed");
            fprintf(stderr, "Emulation finished after one design-update iteration.\n");
            break;
        }
        if (Finished) {
            fprintf(stderr, "Optimization finished at iteration %d\n", Iter);
            break;
        }
    }

    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Executes the integrated geometry generation stage and verifies that an *
* INP output file was created.                                             *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the generation stage succeeded, FALSE otherwise.          *
***************************************************************************/

static int RunIRITIntegrated(OptWorkflowStateStruct *Opt, int mode)
{
    if (!LoadLegacyGeometryData(Opt, Opt->Config.DataNewFile)) {
        fprintf(stderr, "Failed to load legacy geometry CSV: %s\n", Opt->Config.DataNewFile);
        return 0;
    }

    if (!GenerateMicroStructures(Opt, mode)) {
        UnloadLegacyGeometryData(Opt);
        return 0;
    }

    UnloadLegacyGeometryData(Opt);

    if(mode != 0)
	return 1;

    if (!OptFileExists(Opt->Config.InpFile)) {
        fprintf(stderr, "Expected INP output was not created: %s\n", Opt->Config.InpFile);
        return 0;
    }

    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Dispatches the selected simulation backend.                            *
***************************************************************************/
static int RunSimulation(OptWorkflowStateStruct *Opt)
{
    switch (Opt->Run.SimulationType) {
    case 0:
        return RunSimulationEmulation(Opt);
    case 1:
        return RunAnsysBatch(Opt);
    case 2:
        return RunCalculixBatch(Opt);
    default:
        fprintf(stderr, "Unsupported simulation_type: %d\n",
                Opt->Run.SimulationType);
        return FALSE;
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Placeholder for a simulation emulator.  The completed implementation   *
*   will synthesize the stress, deformation, SED, volume, and design-output *
*   files consumed by ProcessIterationData.                                *
***************************************************************************/
static int RunSimulationEmulation(OptWorkflowStateStruct *Opt)
{
    FILE *Sed = NULL, *Volume = NULL;
    CagdRType TargetSed = Opt->Config.SedValue * Opt->Config.SafetyFactorSed;
    int I;

    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);
    if (!ReadINP(Opt, Opt->Config.InpFile) || Opt->Elements.Count < 2)
        goto Failure;

    Sed = fopen(Opt->Config.SedFile, "w");
    Volume = fopen(Opt->Config.VolumeFile, "w");
    if (Sed == NULL || Volume == NULL)
        goto Failure;

    fprintf(Sed, "Element Number\tENERGYPOTENTIAL/VOLUME ()\n");
    fprintf(Volume, "Element Number\tVOLUME\n");
    for (I = 0; I < Opt->Elements.Count; ++I) {
        CagdRType T = (CagdRType) I / (Opt->Elements.Count - 1);
        CagdRType Value = TargetSed * (0.25 + 2.0 * T);
        fprintf(Sed, "%d\t%.17g\n", Opt->Elements.Items[I].Id, (double) Value);
        fprintf(Volume, "%d\t1\n", Opt->Elements.Items[I].Id);
    }

    fclose(Sed);
    fclose(Volume);
    fprintf(stderr, "Emulated %d elements: SED %.6g -> %.6g (target %.6g).\n",
            Opt->Elements.Count, (double) (0.25 * TargetSed),
            (double) (2.25 * TargetSed), (double) TargetSed);
    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);
    return TRUE;

Failure:
    if (Sed != NULL) fclose(Sed);
    if (Volume != NULL) fclose(Volume);
    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);
    return FALSE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Runs the ANSYS Workbench batch command for the current iteration.      *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the external command completed successfully, FALSE        *
*        otherwise.                                                        *
***************************************************************************/

static int RunAnsysBatch(OptWorkflowStateStruct *Opt)
{
    char CmdLine[OPT_MAX_LINE * 2];

    snprintf(CmdLine, sizeof(CmdLine),
        "\"%s\" -B -R \"%s\" -F \"%s\"",
        Opt->Config.AnsysExe,
        Opt->Config.AnsysJournal,
        Opt->Config.AnsysProject);

    return RunCommandLineBlocking(CmdLine, Opt->Config.WorkDir);
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Completes the mesh-only INP produced by IRIT, runs CalculiX, and        *
*   converts its ASCII .dat results to the small two-column files consumed *
*   by the existing optimization code.                                     *
*                                                                          *
*   CalculiX is run in SI units (m, N, Pa).  The generated text outputs use *
*   the same units as the former ANSYS exports: mm, MPa, mm^3 and           *
*   mJ/mm^3.  Numerically, 1 mJ/mm^3 equals 1 MPa.                          *
***************************************************************************/
static int RunCalculixBatch(OptWorkflowStateStruct *Opt)
{
    enum {
        CCX_SEC_NONE,
        CCX_SEC_DISP,
        CCX_SEC_STRESS,
        CCX_SEC_ENERGY,
        CCX_SEC_VOLUME
    } Section = CCX_SEC_NONE;
    const CagdRType TotalLoadY = -1000.0; /* Total vertical force, N. */
    const CagdRType PaToMPa = 1.0e-6;
    const CagdRType MToMm = 1.0e3;
    const CagdRType M3ToMm3 = 1.0e9;
    const char *EnvExe;
    char CcxExe[OPT_MAX_PATH];
    char TempInp[OPT_MAX_PATH];
    char JobName[OPT_MAX_NAME];
    char DatFile[OPT_MAX_PATH];
    char LogFile[OPT_MAX_PATH];
    char CmdLine[OPT_MAX_LINE * 2];
    char Line[OPT_MAX_LINE];
    FILE *Src = NULL, *Dst = NULL, *Dat = NULL;
    FILE *Sed = NULL, *VolumeOut = NULL, *StressOut = NULL;
    FILE *DeflectionOut = NULL, *DesignOut = NULL;
    CagdRType *ElementEnergy = NULL, *ElementVolume = NULL;
    CagdRType *ElementMises = NULL, *NodeDisp = NULL;
    CagdRType *NodeStressSum = NULL;
    int *NodeStressCount = NULL;
    int NodeMaxId = 0, ElementMaxId = 0;
    int FixedCount = 0, LoadedCount = 0;
    int DispCount = 0, EnergyCount = 0, VolumeCount = 0;
    int I, N, InNodeSection = FALSE, SkipGeneratedElset = FALSE;
    CagdRType MinZ = 1.0e300, MaxZ = -1.0e300, ZTolerance;
    CagdRType MaxDeflectionM = 0.0, TotalVolumeM3 = 0.0;
    int Ok = FALSE;

    /* Read the IRIT mesh first; this also gives us reliable node IDs and
     * coordinates for constructing the support and load node sets. */
    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);
    if (!ReadINP(Opt, Opt->Config.InpFile) ||
        Opt->Nodes.Count == 0 || Opt->Elements.Count == 0) {
        fprintf(stderr, "CalculiX: cannot read IRIT mesh %s\n",
                Opt->Config.InpFile);
        goto Cleanup;
    }

    for (I = 0; I < Opt->Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt->Nodes.Items[I];
        MinZ = IRIT_MIN(MinZ, Node->Z);
        MaxZ = IRIT_MAX(MaxZ, Node->Z);
        NodeMaxId = IRIT_MAX(NodeMaxId, Node->Id);
    }
    for (I = 0; I < Opt->Elements.Count; ++I)
        ElementMaxId = IRIT_MAX(ElementMaxId, Opt->Elements.Items[I].Id);

    ZTolerance = IRIT_MAX(1.0, fabs(MaxZ - MinZ)) * 1.0e-9;
    for (I = 0; I < Opt->Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt->Nodes.Items[I];
        if (fabs(Node->Z - MinZ) <= ZTolerance)
            ++FixedCount;
        if (fabs(Node->Z - MaxZ) <= ZTolerance &&
            Node->Y > -0.015 && Node->Y < 0.015)
            ++LoadedCount;
    }
    if (FixedCount == 0 || LoadedCount == 0) {
        fprintf(stderr,
                "CalculiX: failed to select BC nodes (fixed=%d, loaded=%d).\n",
                FixedCount, LoadedCount);
        goto Cleanup;
    }

    /* Rewrite node coordinates with bounded precision.  IRIT's full decimal
     * strings can exceed CalculiX's numeric input field width.  Connectivity
     * and IDs are copied unchanged. */
    snprintf(TempInp, sizeof(TempInp), "%s.ccx.tmp", Opt->Config.InpFile);
    Src = fopen(Opt->Config.InpFile, "r");
    Dst = fopen(TempInp, "w");
    if (Src == NULL || Dst == NULL)
        goto Cleanup;

    while (fgets(Line, sizeof(Line), Src) != NULL) {
        char Copy[OPT_MAX_LINE];
        char *P;

        snprintf(Copy, sizeof(Copy), "%s", Line);
        P = OptTrim(Copy);

        /* Makes the function safe if it is called twice on the same INP. */
        if (strcmp(P, "** BEGIN ADDED CALCULIX ANALYSIS DEFINITION") == 0)
            break;

        if (*P == '*') {
            /* IRIT emits one SET_n selection for every saved object.  These
             * sets are not used by this workflow and become very large for
             * dense lattices.  Remove only SET_n; retain EL-ALL and all BC
             * node sets used by the CalculiX analysis. */
            if (OPT_STRNICMP(P, "*ELSET, ELSET=SET_", 18) == 0) {
                SkipGeneratedElset = TRUE;
                InNodeSection = FALSE;
                continue;
            }
            SkipGeneratedElset = FALSE;
            InNodeSection = OPT_STRNICMP(P, "*NODE", 5) == 0;
            fputs(Line, Dst);
        }
        else if (SkipGeneratedElset) {
            continue;
        }
        else if (InNodeSection && *P != '\0' && strncmp(P, "**", 2) != 0) {
            int Id;
            double X, Y, Z;
            if (sscanf(P, "%d , %lf , %lf , %lf", &Id, &X, &Y, &Z) == 4)
                fprintf(Dst, "%d, %.11E, %.11E, %.11E\n", Id, X, Y, Z);
            else
                fputs(Line, Dst);
        }
        else {
            fputs(Line, Dst);
        }
    }
    fclose(Src);
    Src = NULL;

    fprintf(Dst, "** BEGIN ADDED CALCULIX ANALYSIS DEFINITION\n");
    fprintf(Dst, "** Support: minimum Z. Load: maximum Z and -0.015 < Y < 0.015.\n");
    fprintf(Dst, "*NSET, NSET=FIXED\n");
    N = 0;
    for (I = 0; I < Opt->Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt->Nodes.Items[I];
        if (fabs(Node->Z - MinZ) <= ZTolerance) {
            fprintf(Dst, "%d%s", Node->Id,
                    (++N % 16 == 0 || N == FixedCount) ? "\n" : ",");
        }
    }
    fprintf(Dst, "*NSET, NSET=LOADED\n");
    N = 0;
    for (I = 0; I < Opt->Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt->Nodes.Items[I];
        if (fabs(Node->Z - MaxZ) <= ZTolerance &&
            Node->Y > -0.015 && Node->Y < 0.015) {
            fprintf(Dst, "%d%s", Node->Id,
                    (++N % 16 == 0 || N == LoadedCount) ? "\n" : ",");
        }
    }

    /* Material, section, constraints, load, and all data needed by the
     * optimizer are kept in the same standalone INP file. */
    fprintf(Dst, "*MATERIAL, NAME=OPT_MATERIAL\n");
    fprintf(Dst, "*ELASTIC\n%.17g, %.17g\n",
            (double) Opt->Config.YoungsModulus,
            (double) Opt->Config.PoissonRatio);
    fprintf(Dst, "*DENSITY\n%.17g\n", (double) Opt->Config.Density);
    fprintf(Dst, "*SOLID SECTION, ELSET=EL-ALL, MATERIAL=OPT_MATERIAL\n");
    /* The iterative Cholesky solver is substantially more memory-efficient
     * than the default direct solver for multi-million-DOF lattice meshes. */
    fprintf(Dst, "*STEP\n*STATIC, SOLVER=ITERATIVECHOLESKY\n");
    fprintf(Dst, "*BOUNDARY\nFIXED, 1, 3, 0.0\n");
    fprintf(Dst, "*CLOAD\nLOADED, 2, %.17g\n",
            (double) (TotalLoadY / LoadedCount));
    /* Keep the .dat output intentionally small.  The optimizer needs nodal
     * displacement for DesignParameters and element energy/volume for SED;
     * it does not need a duplicate .frd file or integration-point stresses. */
    fprintf(Dst, "*NODE PRINT, NSET=N-ALL\nU\n");
    /* ELSE is whole-element strain energy; ENER would instead emit one
     * energy-density row per integration point and cannot be paired directly
     * with the whole-element EVOL value. */
    fprintf(Dst, "*EL PRINT, ELSET=EL-ALL\nELSE, EVOL\n");
    fprintf(Dst, "*END STEP\n");
    fprintf(Dst, "** END ADDED CALCULIX ANALYSIS DEFINITION\n");
    fclose(Dst);
    Dst = NULL;

#ifdef _WIN32
    if (!MoveFileExA(TempInp, Opt->Config.InpFile,
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_COPY_ALLOWED)) {
        fprintf(stderr, "CalculiX: cannot replace completed INP file.\n");
        goto Cleanup;
    }
#else
    if (rename(TempInp, Opt->Config.InpFile) != 0)
        goto Cleanup;
#endif

    /* CALCULIX_CCX permits another installation without changing the code. */
    EnvExe = getenv("CALCULIX_CCX");
    snprintf(CcxExe, sizeof(CcxExe), "%s",
             EnvExe != NULL && *EnvExe != '\0' ? EnvExe :
             "C:\\irit\\CalculiX-2.23.0-win-x64\\bin\\ccx_MT.exe");
    if (!OptFileExists(CcxExe)) {
        fprintf(stderr, "CalculiX executable not found: %s\n", CcxExe);
        goto Cleanup;
    }

#ifdef _WIN32
    {
        SYSTEM_INFO SysInfo;
        char ThreadText[32], Override[32];
        DWORD ThreadCount, OverrideLen;

        GetSystemInfo(&SysInfo);
        ThreadCount = SysInfo.dwNumberOfProcessors;
        /* Large CalculiX models are normally limited by memory bandwidth;
         * using every logical CPU can be slower and consumes more memory. */
        if (ThreadCount > 8)
            ThreadCount = 8;
        OverrideLen = GetEnvironmentVariableA("CALCULIX_THREADS", Override,
                                               sizeof(Override));
        if (OverrideLen > 0 && OverrideLen < sizeof(Override)) {
            int Requested = atoi(Override);
            if (Requested > 0)
                ThreadCount = (DWORD) Requested;
        }
        if (ThreadCount < 1)
            ThreadCount = 1;
        snprintf(ThreadText, sizeof(ThreadText), "%lu",
                 (unsigned long) ThreadCount);
        SetEnvironmentVariableA("OMP_NUM_THREADS", ThreadText);
        SetEnvironmentVariableA("CCX_NPROC_STIFFNESS", ThreadText);
        SetEnvironmentVariableA("CCX_NPROC_RESULTS", ThreadText);
        SetEnvironmentVariableA("CCX_NPROC_EQUATION_SOLVER", ThreadText);
        fprintf(stderr, "CalculiX thread count: %s\n", ThreadText);
    }
#endif

    snprintf(JobName, sizeof(JobName), "%s",
             OptPathFileName(Opt->Config.InpFile));
    {
        char *Dot = strrchr(JobName, '.');
        if (Dot != NULL)
            *Dot = '\0';
    }
    snprintf(CmdLine, sizeof(CmdLine), "\"%s\" -i \"%s\"", CcxExe, JobName);
    snprintf(LogFile, sizeof(LogFile), "%s%c%s_ccx.log",
             Opt->Config.WorkDir, OPT_PATH_SEP, JobName);
    fprintf(stderr, "CalculiX solver log: %s\n", LogFile);
    if (!RunCommandLineBlockingToFile(CmdLine, Opt->Config.WorkDir, LogFile)) {
        fprintf(stderr,
                "CalculiX solve failed for job %s; inspect %s\n",
                JobName, LogFile);
        goto Cleanup;
    }

    snprintf(DatFile, sizeof(DatFile), "%s%c%s.dat",
             Opt->Config.WorkDir, OPT_PATH_SEP, JobName);
    Dat = fopen(DatFile, "r");
    if (Dat == NULL) {
        fprintf(stderr, "CalculiX result file was not created: %s\n", DatFile);
        goto Cleanup;
    }

    ElementEnergy = (CagdRType *) OptAllocZero(ElementMaxId + 1,
                                                sizeof(CagdRType));
    ElementVolume = (CagdRType *) OptAllocZero(ElementMaxId + 1,
                                                sizeof(CagdRType));
    NodeDisp = (CagdRType *) OptAllocZero(NodeMaxId + 1, sizeof(CagdRType));
    if (ElementEnergy == NULL || ElementVolume == NULL || NodeDisp == NULL)
        goto Cleanup;

    /* The .dat headings make a much safer interface than fixed-column FRD.
     * ELSE and EVOL are whole-element values; no Gauss-point approximation
     * is needed for the element-average SED. */
    while (fgets(Line, sizeof(Line), Dat) != NULL) {
        char *P = OptTrim(Line);
        int Id, Ip;
        double A, B, C, D, E, F;

        if (strstr(P, "displacements (vx,vy,vz)") != NULL) {
            Section = CCX_SEC_DISP;
            continue;
        }
        if (strstr(P, "stresses (elem, integ.pnt.") != NULL) {
            Section = CCX_SEC_STRESS;
            continue;
        }
        if (strstr(P, "internal energy (element, energy)") != NULL) {
            Section = CCX_SEC_ENERGY;
            continue;
        }
        if (strstr(P, "volume (element, volume)") != NULL) {
            Section = CCX_SEC_VOLUME;
            continue;
        }
        if (strstr(P, "internal energy density") != NULL) {
            Section = CCX_SEC_NONE; /* ENER is diagnostic; ELSE/EVOL is exact. */
            continue;
        }

        if (Section == CCX_SEC_DISP &&
            sscanf(P, "%d %lf %lf %lf", &Id, &A, &B, &C) == 4 &&
            Id > 0 && Id <= NodeMaxId) {
            NodeDisp[Id] = (CagdRType) sqrt(A * A + B * B + C * C);
            MaxDeflectionM = IRIT_MAX(MaxDeflectionM, NodeDisp[Id]);
            ++DispCount;
        }
        else if (Section == CCX_SEC_ENERGY &&
                 sscanf(P, "%d %lf", &Id, &A) == 2 &&
                 Id > 0 && Id <= ElementMaxId) {
            ElementEnergy[Id] = (CagdRType) A;
            ++EnergyCount;
        }
        else if (Section == CCX_SEC_VOLUME &&
                 sscanf(P, "%d %lf", &Id, &A) == 2 &&
                 Id > 0 && Id <= ElementMaxId) {
            ElementVolume[Id] = (CagdRType) A;
            TotalVolumeM3 += (CagdRType) A;
            ++VolumeCount;
        }
    }
    fclose(Dat);
    Dat = NULL;

    if (DispCount < Opt->Nodes.Count ||
        EnergyCount != Opt->Elements.Count ||
        VolumeCount != Opt->Elements.Count) {
        fprintf(stderr,
                "CalculiX output count mismatch: disp=%d/%d, energy=%d/%d, volume=%d/%d.\n",
                DispCount, Opt->Nodes.Count, EnergyCount, Opt->Elements.Count,
                VolumeCount, Opt->Elements.Count);
        goto Cleanup;
    }

    Sed = fopen(Opt->Config.SedFile, "w");
    VolumeOut = fopen(Opt->Config.VolumeFile, "w");
    StressOut = fopen(Opt->Config.StressFile, "w");
    DeflectionOut = fopen(Opt->Config.DeflectionFile, "w");
    DesignOut = fopen(Opt->Config.DesignParametersFile, "w");
    if (Sed == NULL || VolumeOut == NULL || StressOut == NULL ||
        DeflectionOut == NULL || DesignOut == NULL)
        goto Cleanup;

    fprintf(Sed, "Element Number\tENERGYPOTENTIAL/VOLUME (mJ/mm^3)\n");
    fprintf(VolumeOut, "Element Number\tVOLUME (mm^3)\n");
    for (I = 0; I < Opt->Elements.Count; ++I) {
        int Id = Opt->Elements.Items[I].Id;
        if (ElementVolume[Id] <= 0.0)
            goto Cleanup;
        fprintf(Sed, "%d\t%.17g\n", Id,
                (double) (ElementEnergy[Id] / ElementVolume[Id] * PaToMPa));
        fprintf(VolumeOut, "%d\t%.17g\n", Id,
                (double) (ElementVolume[Id] * M3ToMm3));
    }

    /* Keep the configured legacy file present for ProcessIterationData.
     * Stress is deliberately not requested from CalculiX in optimization
     * mode because it is not used by any design-update algorithm. */
    fprintf(StressOut,
            "Node Number\tEquivalent Stress (not requested)\n");
    fprintf(DeflectionOut, "Node Number\tTotal Deformation (mm)\n");
    for (I = 0; I < Opt->Nodes.Count; ++I) {
        int Id = Opt->Nodes.Items[I].Id;
        fprintf(DeflectionOut, "%d\t%.17g\n", Id,
                (double) (NodeDisp[Id] * MToMm));
    }

    /* This deliberately matches ReadLastRowOfDesignParameters(): the first
     * numeric column is maximum deformation in mm and the second is mass kg. */
    fprintf(DesignOut, "Name,P1,P2\n");
    fprintf(DesignOut, "DP 0,%.17g,%.17g\n",
            (double) (MaxDeflectionM * MToMm),
            (double) (TotalVolumeM3 * Opt->Config.Density));

    fclose(Sed); Sed = NULL;
    fclose(VolumeOut); VolumeOut = NULL;
    fclose(StressOut); StressOut = NULL;
    fclose(DeflectionOut); DeflectionOut = NULL;
    fclose(DesignOut); DesignOut = NULL;

    fprintf(stderr,
            "CalculiX completed: %d elements, max deflection %.6g mm, mass %.6g kg.\n",
            Opt->Elements.Count, (double) (MaxDeflectionM * MToMm),
            (double) (TotalVolumeM3 * Opt->Config.Density));
    Ok = TRUE;

Cleanup:
    if (Src != NULL) fclose(Src);
    if (Dst != NULL) fclose(Dst);
    if (Dat != NULL) fclose(Dat);
    if (Sed != NULL) fclose(Sed);
    if (VolumeOut != NULL) fclose(VolumeOut);
    if (StressOut != NULL) fclose(StressOut);
    if (DeflectionOut != NULL) fclose(DeflectionOut);
    if (DesignOut != NULL) fclose(DesignOut);
    if (ElementEnergy != NULL) IritFree(ElementEnergy);
    if (ElementVolume != NULL) IritFree(ElementVolume);
    if (ElementMises != NULL) IritFree(ElementMises);
    if (NodeDisp != NULL) IritFree(NodeDisp);
    if (NodeStressSum != NULL) IritFree(NodeStressSum);
    if (NodeStressCount != NULL) IritFree(NodeStressCount);
    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);
    return Ok;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Parses the iteration outputs, updates the tile state, and writes the   *
* next legacy thickness data file.                                         *
*                                                                          *
* PARAMETERS:                                                              *
*   Iter: Iteration index to process.                                      *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the end condition was met, FALSE otherwise.               *
***************************************************************************/

static int ProcessIterationData(OptWorkflowStateStruct* Opt, int Iter)
{
    CagdRType* PrevThickness;
    int TileCount = Opt->TotalTiles;
    char TargetBaseFolder[OPT_MAX_PATH];

    NodeArrayReset(&Opt->Nodes);
    ElementArrayReset(&Opt->Elements);

    ScalarFieldReset(&Opt->StressData);
    ScalarFieldReset(&Opt->DeflectionData);
    ScalarFieldReset(&Opt->SedData);
    ScalarFieldReset(&Opt->VolumeData);

    OptPathJoin(TargetBaseFolder,
        sizeof(TargetBaseFolder),
        Opt->Config.BaseFolder,
        Opt->FolderName);

    if (!OrganizeResults(Opt, Opt->Config.WorkDir, TargetBaseFolder, Iter))
        return 0;

    if (!ReadINP(Opt, Opt->Config.InpFile))
        OptFatalF("Failed to parse INP file: %s", Opt->Config.InpFile);

    if (!ReadSimulationData(Opt->Config.StressFile, &Opt->StressData))
        OptFatalF("Failed to parse stress file: %s", Opt->Config.StressFile);

    if (!ReadSimulationData(Opt->Config.DeflectionFile, &Opt->DeflectionData))
        OptFatalF("Failed to parse deflection file: %s", Opt->Config.DeflectionFile);

    if (!ReadSimulationData(Opt->Config.SedFile, &Opt->SedData))
        OptFatalF("Failed to parse SED file: %s", Opt->Config.SedFile);

    if (Opt->Config.VolumeFile[0] == '\0') 
        OptFatal("Volume file path is empty. Check volume_file in opt_config.txt and ConfigInputHandler.");

    if (!ReadSimulationData(Opt->Config.VolumeFile, &Opt->VolumeData)) 
        OptFatalF("Failed to parse volume file: %s", Opt->Config.VolumeFile);

    PrevThickness = ReadTileData(Opt, Opt->Config.DataNewFile,
        TileCount * Opt->LMin);
    if (PrevThickness == NULL)
        OptFatalF("Failed to parse tile data file: %s",
            Opt->Config.DataNewFile);

    /*
     * Refresh ONLY the live analysis fields from the just-finished ANSYS
     * solution.  This does not change any lattice thickness.
     */
    UpdateTileAnalysisData(Opt);

    switch (Opt->Run.OptimizationAlgorithm) {
        case 0:
            Opt->LastThicknessChange =
                DesignUpdate_SED_Ratio(Opt, TileCount * Opt->LMin, PrevThickness, Iter);
            break;
        case 1:
            Opt->LastThicknessChange =
                DesignUpdate_MMA(Opt, TileCount * Opt->LMin, PrevThickness, Iter);
            break;
        case 2:
            Opt->LastThicknessChange =
                DesignUpdate_OC(Opt, TileCount * Opt->LMin, PrevThickness, Iter);
            break;
        case 3:
            Opt->LastThicknessChange =
                DesignUpdate_ASAF(Opt, TileCount * Opt->LMin, PrevThickness, Iter);
            break;
        default:
            fprintf(stderr, "Unknown optimization_algorithm = %d\n",
                    Opt->Run.OptimizationAlgorithm);
            IritFree(PrevThickness);
            return 0;
    }

    NeighbourMatching(Opt);

#ifdef _WIN32
    {
        char DataPrev1[OPT_MAX_PATH];
        char DataPrev2[OPT_MAX_PATH];

        OptMakeIndexedFileName(Opt->Config.DataPrevFile,
                               1,
                               DataPrev1,
                               sizeof(DataPrev1));
        OptMakeIndexedFileName(Opt->Config.DataPrevFile,
                               2,
                               DataPrev2,
                               sizeof(DataPrev2));

        if (OptFileExists(DataPrev1))
            CopyFileA(DataPrev1, DataPrev2, FALSE);

        CopyFileA(Opt->Config.DataNewFile, DataPrev1, FALSE);
        CopyFileA(Opt->Config.DataNewFile, Opt->Config.DataPrevFile, FALSE);
        DeleteFileA(Opt->Config.DataNewFile);
    }
#endif

    CreateLegacyDataFile(Opt);
    //(void)EnforceExactVolumeConstraint(Opt,PrevThickness,TileCount * Opt->LMin, Iter);

    /*
     * Measure convergence only after ALL thickness-changing operations
     * (design update, neighbour matching and exact-volume correction).
     */
    Opt->LastThicknessChange = ComputeMaxThicknessChange(Opt,
                                                         PrevThickness,
                                                         TileCount * Opt->LMin);

    fprintf(stderr,
            "Final max thickness change for iteration %d: %.8g (tol %.8g)\n",
            Iter,
            Opt->LastThicknessChange,
            Opt->Run.ThicknessChangeTolerance);

    if (Opt->LastThicknessChange <= Opt->Run.ThicknessChangeTolerance) {
        ++Opt->ConsecutiveSmallChanges;
        fprintf(stderr,
                "Small-change convergence count: %d/%d\n",
                Opt->ConsecutiveSmallChanges,
                Opt->Run.ConvergenceIterations);
    }
    else {
        Opt->ConsecutiveSmallChanges = 0;
    }

    IritFree(PrevThickness);

    if (Opt->ConsecutiveSmallChanges >= Opt->Run.ConvergenceIterations) {
        fprintf(stderr,
                "Optimization converged: max thickness change <= %.8g "
                "for %d consecutive iterations.\n",
                Opt->Run.ThicknessChangeTolerance,
                Opt->Run.ConvergenceIterations);
        return 1;
    }

    return 0;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static CagdBType EndConditionMet(const OptWorkflowStateStruct *Opt)
{
    /*
     * Reserved for future independent stopping criteria.
     * Thickness convergence is evaluated at the end of ProcessIterationData,
     * after all thickness-changing operations have completed.
     */
    (void)Opt;
    return FALSE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reads the design-parameter summary that was produced by the solver.    *
*                                                                          *
* PARAMETERS:                                                              *
*   FileName:   CSV file to parse.                                         *
*   Deflection: Location to store the last deflection value.               *
*   Mass:       Location to store the last mass value.                     *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the file was parsed successfully, FALSE otherwise.        *
***************************************************************************/

static int ReadLastRowOfDesignParameters(const char *FileName,
                                         CagdRType *Deflection,
                                         CagdRType *Mass)
{
    FILE *F;
    char Line[OPT_MAX_LINE];
    CagdRType LastDef = 0.0, LastMass = 0.0;

    F = fopen(FileName, "r");
    if (F == NULL)
        return 0;

    while (fgets(Line, sizeof(Line), F) != NULL) {
        char *P = OptTrim(Line);
        char *Tok;
        int Col = 0;
        char *Parts[4] = { NULL, NULL, NULL, NULL };

        if (*P == '\0' || *P == '#' || strncmp(P, "Name", 4) == 0)
            continue;

        Tok = strtok(P, ",");
        while (Tok != NULL && Col < 4) {
            Parts[Col++] = Tok;
            Tok = strtok(NULL, ",");
        }

        if (Col >= 3) {
            CagdRType TmpDef, TmpMass;
            if (OptParseReal(Parts[1], &TmpDef) && OptParseReal(Parts[2], &TmpMass)) {
                LastDef = TmpDef;
                LastMass = TmpMass;
            }
        }
    }

    fclose(F);
    *Deflection = LastDef;
    *Mass = LastMass;
    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Copies the iteration files into a run folder and appends a summary     *
* record for the iteration.                                                *
*                                                                          *
* PARAMETERS:                                                              *
*   SourceFolder:     Working directory to archive.                        *
*   TargetBaseFolder: Destination run folder.                              *
*   Iter:             Iteration index.                                     *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the results were organized successfully, FALSE otherwise. *
***************************************************************************/
static int OrganizeResults(OptWorkflowStateStruct *Opt,
                           const char *SourceFolder,
                           const char *TargetBaseFolder,
                           int Iter)
{
    char TargetFolder[OPT_MAX_PATH];
    char SummaryPath[OPT_MAX_PATH];
    char SummaryLine[OPT_MAX_LINE];

    snprintf(TargetFolder,
             sizeof(TargetFolder),
             "%s%cRun_%d",
             TargetBaseFolder,
             OPT_PATH_SEP,
             Iter);
    EnsureDirectoryRecursive(TargetFolder);

    if (!CopyDirectoryRecursive(Opt, SourceFolder, TargetFolder)) {
        fprintf(stderr, "Warning: failed to copy results folder to %s\n", TargetFolder);
    }

    if (!ReadLastRowOfDesignParameters(Opt->Config.DesignParametersFile,
                                       &Opt->MaxDef,
                                       &Opt->Mass)) {
        fprintf(stderr,
                "Failed to read design parameters CSV: %s\n",
                Opt->Config.DesignParametersFile);
        return 0;
    }

    RESULT_AT(Opt, Iter - 1, 0) = (CagdRType)Iter;
    RESULT_AT(Opt, Iter - 1, 1) = Opt->MaxDef;
    RESULT_AT(Opt, Iter - 1, 2) = Opt->Mass;

    OptPathJoin(SummaryPath, sizeof(SummaryPath), TargetBaseFolder, Opt->Config.SummaryFile);

    /* Write CSV header once, when the summary file is first created. */
    if (!OptFileExists(SummaryPath)) {
        OptWriteLine(SummaryPath,
            "Iteration,MaxDeflection,Mass,TimeIritStart,TimeIritFinish,TimeSimulationFinish",
            0);
    }

    snprintf(SummaryLine, sizeof(SummaryLine),
        "%d,%.4f,%.4f,%s,%s,%s",
        Iter,
        (double)Opt->MaxDef,
        (double)Opt->Mass,
        Opt->TimeIritStart,
        Opt->TimeIritFinish,
        Opt->TimeSimulationFinish);
    OptWriteLine(SummaryPath, SummaryLine, 1);

    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reads the generated INP mesh file into the native node and element     *
* containers.                                                              *
*                                                                          *
* PARAMETERS:                                                              *
*   INPPath: Input INP file to parse.                                      *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the file was parsed successfully, FALSE otherwise.        *
***************************************************************************/

static int ReadINP(OptWorkflowStateStruct *Opt, const char *INPPath)
{
    FILE *F;
    char Line[OPT_MAX_LINE];
    enum { SEC_NONE, SEC_NODE, SEC_ELEMENT } CurrentSection = SEC_NONE;

    F = fopen(INPPath, "r");
    if (F == NULL)
        return 0;

    while (fgets(Line, sizeof(Line), F) != NULL) {
        char *P = OptTrim(Line);

        if (*P == '\0' || strncmp(P, "**", 2) == 0)
            continue;

        if (*P == '*') {
            if (OPT_STRNICMP(P, "*NODE", 5) == 0)
                CurrentSection = SEC_NODE;
            else if (OPT_STRNICMP(P, "*ELEMENT", 8) == 0)
                CurrentSection = SEC_ELEMENT;
            else
                CurrentSection = SEC_NONE;
            continue;
        }

        if (CurrentSection == SEC_NODE) {
            char *Tok;
            int Col = 0;
            char *Parts[4] = { NULL, NULL, NULL, NULL };
            OptNodeStruct Node;

            IRIT_ZAP_MEM(&Node, sizeof(Node));

            Tok = strtok(P, ",");
            while (Tok != NULL && Col < 4) {
                Parts[Col++] = OptTrim(Tok);
                Tok = strtok(NULL, ",");
            }

            if (Col >= 4) {
                Node.Id = atoi(Parts[0]);
                Node.X = (CagdRType)atof(Parts[1]);
                Node.Y = (CagdRType)atof(Parts[2]);
                Node.Z = (CagdRType)atof(Parts[3]);
                if (!NodeArrayPush(&Opt->Nodes, &Node)) {
                    fclose(F);
                    return 0;
                }
            }
        }
        else if (CurrentSection == SEC_ELEMENT) {
            char *Tok;
            int Col = 0;
            OptElementStruct Elem;

            IRIT_ZAP_MEM(&Elem, sizeof(Elem));

            Tok = strtok(P, ",");
            while (Tok != NULL) {
                Tok = OptTrim(Tok);
                if (Col == 0)
                    Elem.Id = atoi(Tok);
                else if (Col - 1 < OPT_MAX_NODE_PER_ELEM)
                    Elem.NodeIds[Col - 1] = atoi(Tok);
                ++Col;
                Tok = strtok(NULL, ",");
            }

            Elem.NodeCount = Col - 1;
            if (Elem.Id > 0 && Elem.NodeCount > 0) {
                if (!ElementArrayPush(&Opt->Elements, &Elem)) {
                    fclose(F);
                    return 0;
                }
            }
        }
    }

    fclose(F);
    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reads a tab-separated scalar simulation file into an indexed field.    *
*                                                                          *
* PARAMETERS:                                                              *
*   FileName: Input file to parse.                                         *
*   Field:    Scalar field to populate.                                    *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the file was parsed successfully, FALSE otherwise.        *
***************************************************************************/

static int ReadSimulationData(const char *FileName, OptScalarFieldStruct *Field)
{
    FILE *F;
    char Line[OPT_MAX_LINE];
    int HeaderSkipped = 0;

    F = fopen(FileName, "r");
    if (F == NULL)
        return 0;

    while (fgets(Line, sizeof(Line), F) != NULL) {
        char *P = OptTrim(Line);
        char *IdTok;
        char *ValTok;
        int Id;
        CagdRType Val;

        if (!HeaderSkipped) {
            HeaderSkipped = 1;
            continue;
        }

        if (*P == '\0')
            continue;

        IdTok = strtok(P, "\t");
        ValTok = strtok(NULL, "\t");
        if (IdTok == NULL || ValTok == NULL)
            continue;

        Id = atoi(OptTrim(IdTok));
        Val = (CagdRType)atof(OptTrim(ValTok));
        ScalarFieldSet(Field, Id, Val);
    }

    fclose(F);
    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reads the legacy thickness CSV file into a dense array of thickness    *
* values.                                                                  *
*                                                                          *
* PARAMETERS:                                                              *
*   FileName:      Legacy thickness file to parse.                         *
*   ExpectedCount: Number of thickness entries expected in the file.       *
*                                                                          *
* RETURN VALUE:                                                            *
*   CagdRType *: Allocated thickness array, or NULL on failure.            *
***************************************************************************/

static CagdRType *ReadTileData(const OptWorkflowStateStruct *Opt,
                              const char *FileName, int ExpectedCount)
{
    FILE *F = fopen(FileName, "r");
    CagdRType *Prev = NULL;
    char Line[OPT_MAX_LINE], Extra;
    int N = 0, Id;
    double Thickness, Sed;
    if (!F) return NULL;
    if (!OptReadCsvHeader(F, Opt)) goto Done;
    Prev = (CagdRType *) OptAllocZero(ExpectedCount, sizeof(CagdRType));
    if (!Prev) goto Done;
    while (fgets(Line, sizeof(Line), F)) {
        if (OptIsBlank(Line)) continue;
        if (N >= ExpectedCount ||
            sscanf(Line, "%d,%lf,%lf %c", &Id, &Thickness, &Sed, &Extra) != 3 ||
            Id != N + 1 || !isfinite(Thickness) || Thickness <= 0.0 ||
            !isfinite(Sed)) goto Invalid;
        Prev[N++] = Thickness;
    }
    if (N != ExpectedCount) goto Invalid;
    goto Done;
Invalid:
    fprintf(stderr, "Invalid/missing/extra CSV design row at index %d in %s\n", N, FileName);
    IritFree(Prev); Prev = NULL;
Done:
    fclose(F);
    return Prev;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Finds a node by its identifier in the current node container.          *
*                                                                          *
* PARAMETERS:                                                              *
*   NodeId: Identifier of the node to locate.                              *
*                                                                          *
* RETURN VALUE:                                                            *
*   OptNodeStruct *: Matching node if found, NULL otherwise.               *
***************************************************************************/

static OptNodeStruct *FindNodeById(OptWorkflowStateStruct *Opt, int NodeId)
{
    int Guess = NodeId - 1;
    int I;

    if (Guess >= 0 && Guess < Opt->Nodes.Count && Opt->Nodes.Items[Guess].Id == NodeId)
        return &Opt->Nodes.Items[Guess];

    for (I = 0; I < Opt->Nodes.Count; ++I) {
        if (Opt->Nodes.Items[I].Id == NodeId)
            return &Opt->Nodes.Items[I];
    }

    return NULL;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Refreshes the persistent TileStructure with data from the current      *
*   ANSYS solution.  Topology and thickness are preserved.                  *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static void UpdateTileAnalysisData(OptWorkflowStateStruct* Opt)
{
    int I, J, K, L, X;
    int Count = 0;

    /*
     * TileStructure was created once in InitializeWorkflowState().
     * Clear only fields that belong to the previous ANSYS solution.
     * SubTile->Thickness is deliberately preserved.
     */
    if (Opt->TileStructure == NULL)
        OptFatal("TileStructure was not initialized before analysis update");
    ClearTileAnalysisData(Opt);

    if (Opt->Elements.Count != Opt->MeshOwnerCount * Opt->MeshSample) {
        fprintf(stderr,
                "Mesh ownership mismatch: elements=%d, owners=%d, configured elements/owner=%d.\n",
                Opt->Elements.Count, Opt->MeshOwnerCount, Opt->MeshSample);
        OptFatal("Mesh element count does not match exported ownership map");
    }
    for (I = 0; I < Opt->Elements.Count; ++I) {
        OptMeshOwnerStruct *Owner = &Opt->MeshOwners[I / Opt->MeshSample];
        if (Opt->Elements.Items[I].Id != I + 1)
            OptFatal("ANSYS INP element order/IDs differ from exported mesh");
        if (!SubTileAddElementPtr(&Opt->TileStructure[Owner->Tile].SubTiles[Owner->SubTile],
                                  &Opt->Elements.Items[I]))
            OptFatal("Cannot assign element to body/tile/subtile");
    }

    /* ------------------------------------------------------------
     * Assign nodal stress values.
     * ------------------------------------------------------------ */
    for (I = 0; I < Opt->Nodes.Count; ++I) {
        CagdRType Stress;

        if (ScalarFieldTryGet(&Opt->StressData,
            Opt->Nodes.Items[I].Id,
            &Stress)) {
            Opt->Nodes.Items[I].Stress = Stress;
        }
        else {
            Opt->Nodes.Items[I].Stress = -1.0;
        }
    }

    /* ------------------------------------------------------------
     * Assign element stress, SED, volume, strain energy, and
     * compliance contribution.
     * ------------------------------------------------------------ */
    for (I = 0; I < Opt->Elements.Count; ++I) {
        CagdRType Avg = 0.0;
        CagdRType Min = 50000.0;
        CagdRType Max = 0.0;
        CagdRType SedVal = 0.0;
        CagdRType VolVal = 0.0;
        int UsedNodeCount =
            IRIT_MIN(Opt->Config.NodesPerElement,
                Opt->Elements.Items[I].NodeCount);
        int N;

        for (N = 0; N < UsedNodeCount; ++N) {
            OptNodeStruct* Node =
                FindNodeById(Opt, Opt->Elements.Items[I].NodeIds[N]);
            CagdRType S = Node != NULL ? Node->Stress : -1.0;

            Avg += S;

            if (S > Max)
                Max = S;
            if (S < Min)
                Min = S;
        }

        Opt->Elements.Items[I].AvgStress =
            Avg / (UsedNodeCount > 0 ? UsedNodeCount : 1);
        Opt->Elements.Items[I].MinStress = Min;
        Opt->Elements.Items[I].MaxStress = Max;

        if (ScalarFieldTryGet(&Opt->SedData,
            Opt->Elements.Items[I].Id,
            &SedVal)) {
            Opt->Elements.Items[I].SED = SedVal;
        }
        else {
            Opt->Elements.Items[I].SED = 0.0;
        }

        if (ScalarFieldTryGet(&Opt->VolumeData,
            Opt->Elements.Items[I].Id,
            &VolVal)) {
            Opt->Elements.Items[I].Volume = VolVal;
        }
        else {
            Opt->Elements.Items[I].Volume = 0.0;
        }

        /*
         * Element strain energy:
         *     SE_e = SED_e * V_e
         *
         * Element compliance contribution:
         *     C_e ~= 2 * SE_e
         */
        Opt->Elements.Items[I].StrainEnergy =
            Opt->Elements.Items[I].SED *
            Opt->Elements.Items[I].Volume;

        Opt->Elements.Items[I].ComplianceContribution =
            2.0 * Opt->Elements.Items[I].StrainEnergy;
    }

    /* ------------------------------------------------------------
     * Map elements into tiles/subtiles and compute subtile values.
     * ------------------------------------------------------------ */
    for (I = 0; I < Opt->TotalTiles; ++I) {
        OptTileStruct* Tile = &Opt->TileStructure[I];

        Tile->Volume = 0.0;
        Tile->StrainEnergy = 0.0;
        Tile->ComplianceContribution = 0.0;
        Tile->Gain = 0.0;

        for (L = 0; L < Opt->L; ++L) {
            OptSubTileStruct* SubTile = &Tile->SubTiles[L];
            CagdRType SedSum = 0.0;
            CagdRType MaxSed = 0.0;
            CagdRType VolSum = 0.0;
            CagdRType SESum = 0.0;
            CagdRType CSum = 0.0;

            SubTile->Volume = 0.0;
            SubTile->StrainEnergy = 0.0;
            SubTile->ComplianceContribution = 0.0;
            SubTile->Gain = 0.0;

            for (X = 0; X < SubTile->ElementCount; ++X) {
                OptElementStruct* Elem;

                Elem = SubTile->Elements[X];

                SedSum += Elem->SED;

                if (Elem->SED > MaxSed)
                    MaxSed = Elem->SED;

                VolSum += Elem->Volume;
                SESum += Elem->StrainEnergy;
                CSum += Elem->ComplianceContribution;

                ++Count;
            }

            SubTile->Volume = VolSum;
            SubTile->StrainEnergy = SESum;
            SubTile->ComplianceContribution = CSum;

            /*
             * Prefer volume-weighted SED:
             *     SED_s = SE_s / V_s
             *
             * If volume is missing, fall back to the old arithmetic
             * average.
             */
            if (SubTile->Volume > 0.0) {
                SubTile->SED =
                    SubTile->StrainEnergy / SubTile->Volume;
            }
            else {
                SubTile->SED =
                    SedSum / (SubTile->ElementCount > 0 ?
                        SubTile->ElementCount : 1);
            }

            SubTile->MaxSED = MaxSed;

            /*
             * TO-like local update indicator:
             *     Gain = C_s / V_s
             *
             * Since C_s ~= 2 * SE_s, this is proportional to SED.
             */
            if (SubTile->Volume > 0.0) {
                SubTile->Gain =
                    SubTile->ComplianceContribution /
                    SubTile->Volume;
            }
            else {
                SubTile->Gain = 0.0;
            }

            Tile->Volume += SubTile->Volume;
            Tile->StrainEnergy += SubTile->StrainEnergy;
            Tile->ComplianceContribution +=
                SubTile->ComplianceContribution;
        }

        if (Tile->Volume > 0.0) {
            Tile->Gain =
                Tile->ComplianceContribution / Tile->Volume;
        }
        else {
            Tile->Gain = 0.0;
        }

        /* --------------------------------------------------------
         * Compute subtile stress values.
         * -------------------------------------------------------- */
        for (L = 0; L < Opt->L; ++L) {
            OptSubTileStruct* SubTile = &Tile->SubTiles[L];
            CagdRType Avg1 = 0.0;
            CagdRType Min1 = 50000.0;
            CagdRType Max1 = 0.0;

            for (X = 0; X < SubTile->ElementCount; ++X) {
                OptElementStruct* Elem = SubTile->Elements[X];

                Avg1 += Elem->AvgStress;

                if (Elem->MaxStress > Max1)
                    Max1 = Elem->MaxStress;
                if (Elem->MinStress < Min1)
                    Min1 = Elem->MinStress;
            }

            SubTile->AvgStress =
                Avg1 / (SubTile->ElementCount > 0 ?
                    SubTile->ElementCount : 1);
            SubTile->MinStress = Min1;
            SubTile->MaxStress = Max1;
        }

        /* --------------------------------------------------------
         * Compute tile stress values.
         * -------------------------------------------------------- */
        {
            CagdRType Avg2 = 0.0;
            CagdRType Min2 = 50000.0;
            CagdRType Max2 = 0.0;

            for (L = 0; L < Opt->L; ++L) {
                Avg2 += Tile->SubTiles[L].AvgStress;

                if (Tile->SubTiles[L].MaxStress > Max2)
                    Max2 = Tile->SubTiles[L].MaxStress;
                if (Tile->SubTiles[L].MinStress < Min2)
                    Min2 = Tile->SubTiles[L].MinStress;
            }

            Tile->AvgStress = Avg2 / Opt->L;
            Tile->MinStress = Min2;
            Tile->MaxStress = Max2;
        }
    }


}

/****************************************************************************
* DESCRIPTION:                                                              *
*   Updates the modifiable subtile thicknesses using an OC-like rule, based *
* on the VTS thickness optimization update.                                 *
*                                                                           *
*   The current SED value is used as a sensitivity proxy. The update uses   *
* damping, move limits, and lower/upper thickness bounds to avoid large     *
* jumps and oscillations. A cheap volume proxy is enforced by limiting the  *
* total sum of thickness values, instead of calling the expensive geometric *
* volume calculation every iteration.                                       *
*                                                                           *
*   The returned value is the maximum absolute thickness change in this     *
* iteration. This value is used as an optimization stopping criterion.      *
*                                                                           *
* PARAMETERS:                                                               *
*   Opt:   Optimization workflow state.                                     *
*   Count: Number of modifiable thickness entries.                          *
*   Prev:  Previous thickness values, read from the current data file.      *
*   Iter:  Current optimization iteration number.                           *
*                                                                           *
* RETURN VALUE:                                                             *
*   CagdRType: Maximum absolute thickness change in this iteration.         *
****************************************************************************/
static CagdRType DesignUpdate_OC(OptWorkflowStateStruct *Opt,
                              int Count,
                              const CagdRType *Prev,
                              int Iter)
{
    const CagdRType Eta = 0.5;
    const CagdRType TMin = Opt->Config.MinimumThickness;
    const CagdRType TMax = OPT_MAX_TILE_THICKNESS;
    const CagdRType Move = (TMax - TMin) / 5.0;
    const CagdRType Eps = 1.0e-12;
    const int MaxBisectIter = 60;

    int I, J, K, L, It, Counter;
    int SedCount = 0;
    CagdRType L1 = Eps, L2 = 1.0e12, Lambda = 1.0;
    CagdRType TargetSum, SumOld = 0.0, SumNew = 0.0;
    CagdRType MaxAbsChange = 0.0;
    CagdRType SedAvg = 0.0, SedSqAvg = 0.0;
    CagdRType SedTarget = Opt->Config.SedValue / Opt->Config.SafetyFactorSed;

    (void) Iter;

    if (SedTarget <= Eps)
        SedTarget = Eps;

    TargetSum = Count * Opt->Run.InitialThickness;
    if (TargetSum <= Eps) {
        for (I = 0; I < Count; ++I)
            TargetSum += Prev[I];
    }

    for (I = 0; I < Count; ++I)
        SumOld += Prev[I];

    for (It = 0; It < MaxBisectIter; ++It) {
        SumNew = 0.0;
        Counter = 0;
        Lambda = 0.5 * (L1 + L2);

        for (I = 0; I < Opt->TotalTiles; ++I) {
            OptTileStruct *Tile = &Opt->TileStructure[I];

            for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
                CagdRType TOld = Prev[Counter];
                CagdRType Sed = Tile->SubTiles[L].SED;
                CagdRType Sens = Sed > 0.0 ? Sed / SedTarget : Eps;
                CagdRType TCand = TOld * pow(Sens / Lambda, Eta);
                CagdRType TLow = TOld - Move;
                CagdRType THigh = TOld + Move;
                CagdRType TNew = TCand;

                if (TLow < TMin)
                    TLow = TMin;
                if (THigh > TMax)
                    THigh = TMax;

                if (TNew < TLow)
                    TNew = TLow;
                if (TNew > THigh)
                    TNew = THigh;

                SumNew += TNew;
                ++Counter;
            }
        }



        if (SumNew > TargetSum)
            L1 = Lambda;
        else
            L2 = Lambda;
    }

    Counter = 0;
    SumNew = 0.0;
    for (I = 0; I < Opt->TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt->TileStructure[I];

        for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
            CagdRType TOld = Prev[Counter];
            CagdRType Sed = Tile->SubTiles[L].SED;
            CagdRType Sens = Sed > 0.0 ? Sed / SedTarget : Eps;
            CagdRType TCand = TOld * pow(Sens / Lambda, Eta);
            CagdRType TLow = TOld - Move;
            CagdRType THigh = TOld + Move;
            CagdRType TNew = TCand;

            if (TLow < TMin)
                TLow = TMin;
            if (THigh > TMax)
                THigh = TMax;

            if (TNew < TLow)
                TNew = TLow;
            if (TNew > THigh)
                TNew = THigh;

            Tile->SubTiles[L].Thickness = TNew;
            SumNew += TNew;

            if (fabs(TNew - TOld) > MaxAbsChange)
                MaxAbsChange = fabs(TNew - TOld);

            if (TNew > TMin && Sed > 0.0) {
                SedAvg += Sed;
                SedSqAvg += Sed * Sed;
                ++SedCount;
            }

            ++Counter;
        }
    }



    if (SedCount > 0) {
        SedAvg /= SedCount;
        SedSqAvg /= SedCount;
        Opt->SedAvg = SedAvg;
        Opt->SedR = sqrt(fabs(SedSqAvg - SedAvg * SedAvg));
    }
    else {
        Opt->SedAvg = 0.0;
        Opt->SedR = 0.0;
    }

    fprintf(stderr,
            "DesignUpdate: Tsum old=%.8g, Tsum new=%.8g, target=%.8g, "
            "change=%.8g, SedAvg=%.8g, SedStd=%.8g\n",
            SumOld, SumNew, TargetSum, MaxAbsChange, Opt->SedAvg, Opt->SedR);

    return MaxAbsChange;
}

/***************************************************************************
* DESCRIPTION:
*   Computes the maximum absolute change between the final current design
*   and the previous iteration thickness vector.  This must be called only
*   after every operation that may modify thickness (including neighbour
*   matching and exact-volume correction).
***************************************************************************/
static CagdRType ComputeMaxThicknessChange(const OptWorkflowStateStruct *Opt,
                                           const CagdRType *Prev,
                                           int Count)
{
    int I, L;
    int Counter = 0;
    CagdRType MaxChange = 0.0;

    if (Opt == NULL || Prev == NULL || Count <= 0)
        return 0.0;

    for (I = 0; I < Opt->TotalTiles; ++I) {
        const OptTileStruct *Tile = &Opt->TileStructure[I];

        for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
            CagdRType Delta;

            if (Counter >= Count)
                OptFatal("Thickness count mismatch in ComputeMaxThicknessChange");

            Delta = fabs(Tile->SubTiles[L].Thickness - Prev[Counter]);
            if (Delta > MaxChange)
                MaxChange = Delta;

            ++Counter;
        }
    }

    if (Counter != Count)
        OptFatal("Thickness count mismatch in ComputeMaxThicknessChange");

    return MaxChange;
}
/***************************************************************************
* DESCRIPTION:                                                             *
*   Old direct SED-ratio thickness update. This version is preserved for   *
* comparison or rollback, but is no longer called by the main workflow.    *
*                                                                          *
* PARAMETERS:                                                              *
*   Count: Number of thickness entries in Prev.                            *
*   Prev:  Previous thickness values.                                      *
*   Iter:  Current iteration index.                                        *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static CagdRType DesignUpdate_SED_Ratio(OptWorkflowStateStruct* Opt,
    int Count,
    const CagdRType* Prev,
    int Iter)
{
    int IndCounter = 0;
    int SedrCounter = 0;
    int I, J, K, L;
    CagdRType Avg = 0.0;
    CagdRType Var = 0.0;
    CagdRType MaxThicknessChange = 0.0;
    CagdRType PowFactor = 0.35;
    CagdRType MaxChange = 0.1;
    CagdRType TargetSed;
    CagdRType* Sedr;



    if (Opt == NULL || Prev == NULL || Count <= 0)
        return 0.0;

    if (Iter == 1) {
        PowFactor = 0.4;
        MaxChange = 0.4;
    }

    /* Preserve the original V1 SED-ratio definition. */
    TargetSed = Opt->Config.SedValue * Opt->Config.SafetyFactorSed;
    if (TargetSed <= 0.0)
        OptFatal("Invalid target SED value in DesignUpdate_SED_Ratio");

    Sedr = (CagdRType*)OptAllocZero((size_t)Count, sizeof(CagdRType));
    if (Sedr == NULL)
        OptFatal("Failed to allocate SED scratch array");

    for (I = 0; I < Opt->TotalTiles; ++I) {
        OptTileStruct* Tile = &Opt->TileStructure[I];

        for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
            CagdRType TOld;
            CagdRType TNew;
            CagdRType SedFactor;
            CagdRType Change;
            CagdRType AbsThicknessChange;

            if (IndCounter >= Count) {
                fprintf(stderr,
                    "Warning: DesignUpdate_SED_Ratio reached Count limit early. "
                    "Count = %d\n",
                    Count);
                goto ComputeSedStats;
            }

            TOld = Prev[IndCounter];

            SedFactor = Tile->SubTiles[L].SED / TargetSed;

            /*
             * Safety against zero/negative/invalid SED values.
             * If SedFactor is very small, the tile will shrink,
             * but the minimum-thickness clamp below still protects it.
             */
            if (SedFactor < 0.0)
                SedFactor = 0.0;

            Change = pow(SedFactor, PowFactor);

            CagdRType RawChange = pow(SedFactor, PowFactor);

            fprintf(stderr,
                "Tile %d Sub %d: SED=%g Target=%g Ratio=%g RawChange=%g\n",
                I,
                L,
                Tile->SubTiles[L].SED,
                TargetSed,
                SedFactor,
                RawChange);

            if (Change > 1.0 + MaxChange)
                Change = 1.0 + MaxChange;

            if (Change < 1.0 - MaxChange)
                Change = 1.0 - MaxChange;

            TNew = TOld * Change;

            if (TNew < Opt->Config.MinimumThickness)
                TNew = Opt->Config.MinimumThickness;
            if (TNew > OPT_MAX_TILE_THICKNESS)
                TNew = OPT_MAX_TILE_THICKNESS;

            Tile->SubTiles[L].Thickness = TNew;

            AbsThicknessChange = fabs(TNew - TOld);
            if (AbsThicknessChange > MaxThicknessChange)
                MaxThicknessChange = AbsThicknessChange;

            if (TNew > Opt->Config.MinimumThickness) {
                Sedr[SedrCounter++] = Tile->SubTiles[L].SED;
            }

            ++IndCounter;
        }
    }



ComputeSedStats:

    if (IndCounter != Count) {
        fprintf(stderr,
            "Warning: DesignUpdate_SED_Ratio used %d thickness values, "
            "but Count was %d\n",
            IndCounter,
            Count);
    }

    if (SedrCounter > 0) {
        int I2;

        for (I2 = 0; I2 < SedrCounter; ++I2)
            Avg += Sedr[I2];

        Avg /= SedrCounter;

        for (I2 = 0; I2 < SedrCounter; ++I2)
            Var += IRIT_SQR(Sedr[I2] - Avg);

        Var /= SedrCounter;

        Opt->SedAvg = Avg;
        Opt->SedR = sqrt(Var);
    }
    else {
        Opt->SedAvg = 0.0;
        Opt->SedR = 0.0;
    }

    IritFree(Sedr);

    return MaxThicknessChange;
}

static CagdRType DesignUpdate_ASAF(OptWorkflowStateStruct* Opt,
    int Count,
    const CagdRType* Prev,
    int Iter)
{
    int Counter = 0;
    int SedrCounter = 0;
    int B, U, V, W, L;
    CagdRType Avg = 0.0;
    CagdRType Var = 0.0;
    CagdRType MaxThicknessChange = 0.0;
    CagdRType PowFactor = 0.35;
    CagdRType MaxChange = 0.15;
    CagdRType TargetSed;
    CagdRType* Sedr;

    if (Opt == NULL || Prev == NULL || Count <= 0)
        return 0.0;

    if (Iter == 1) {
        PowFactor = 0.4;
        MaxChange = 0.3;
    }

    TargetSed = Opt->Config.SedValue * Opt->Config.SafetyFactorSed;
    if (TargetSed <= 0.0)
        OptFatal("Invalid target SED value in DesignUpdate_ASAF");

    Sedr = (CagdRType*)OptAllocZero((size_t)Count, sizeof(CagdRType));
    if (Sedr == NULL)
        OptFatal("Failed to allocate SED scratch array");

    /*
     * Use the original physical traversal explicitly:
     *     body -> W -> V -> U -> SubTile.
     *
     * This keeps Prev[Counter], TileStructure and data_new.csv in the same
     * unambiguous order, including multi-body runs.
     */
    for (B = 0; B < Opt->BodyCount; ++B) {
        OptBodyStruct* Body = &Opt->Bodies[B];

        for (W = 0; W < Body->Counts[2]; ++W) {
            for (V = 0; V < Body->Counts[1]; ++V) {
                for (U = 0; U < Body->Counts[0]; ++U) {
                    int TileIndex = Body->TileOffset +
                        U + Body->Counts[0] * (V + Body->Counts[1] * W);
                    OptTileStruct* Tile = &Opt->TileStructure[TileIndex];

                    for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
                        CagdRType TOld;
                        CagdRType TNew;
                        CagdRType SedFactor;
                        CagdRType Change;
                        CagdRType AbsThicknessChange;

                        if (Counter >= Count)
                            OptFatal("DesignUpdate_ASAF thickness count mismatch");

                        TOld = Prev[Counter];
                        SedFactor = Tile->SubTiles[L].SED / TargetSed;
                        if (SedFactor < 0.0 || !isfinite(SedFactor))
                            SedFactor = 0.0;

                        Change = pow(SedFactor, PowFactor);

                        if (Change > 1.0 + MaxChange)
                            Change = 1.0 + MaxChange;
                        if (Change < 1.0 - MaxChange)
                            Change = 1.0 - MaxChange;

                        TNew = TOld * Change;

                        if (TNew < Opt->Config.MinimumThickness)
                            TNew = Opt->Config.MinimumThickness;
                        if (TNew > OPT_MAX_TILE_THICKNESS)
                            TNew = OPT_MAX_TILE_THICKNESS;

                        Tile->SubTiles[L].Thickness = TNew;

                        AbsThicknessChange = fabs(TNew - TOld);
                        if (AbsThicknessChange > MaxThicknessChange)
                            MaxThicknessChange = AbsThicknessChange;

                        if (TNew > Opt->Config.MinimumThickness)
                            Sedr[SedrCounter++] = Tile->SubTiles[L].SED;

                        fprintf(stderr,
                            "Body %d Tile[%d,%d,%d] GTile=%d Sub=%d "
                            "SED=%g Told=%g Tnew=%g Change=%g\n",
                            B + 1, U, V, W, TileIndex, L,
                            Tile->SubTiles[L].SED, TOld, TNew, Change);

                        ++Counter;
                    }
                }
            }
        }
    }

    if (Counter != Count)
        OptFatal("DesignUpdate_ASAF did not consume all thickness values");

    IritFree(Sedr);
    return MaxThicknessChange;
}

/****************************************************************************
* DESCRIPTION:                                                              *
*   Updates the modifiable subtile thicknesses using an OC-like rule, based *
* on the VTS thickness optimization update.                                 *
*                                                                           *
*   The current SED value is used as a sensitivity proxy. The update uses   *
* damping, move limits, and lower/upper thickness bounds to avoid large     *
* jumps and oscillations. A cheap volume proxy is enforced by limiting the  *
* total sum of thickness values, instead of calling the expensive geometric *
* volume calculation every iteration.                                       *
*                                                                           *
*   The returned value is the maximum absolute thickness change in this     *
* iteration. This value is used as an optimization stopping criterion.      *
*                                                                           *
* PARAMETERS:                                                               *
*   Opt:   Optimization workflow state.                                     *
*   Count: Number of modifiable thickness entries.                          *
*   Prev:  Previous thickness values, read from the current data file.      *
*   Iter:  Current optimization iteration number.                           *
*                                                                           *
* RETURN VALUE:                                                             *
*   CagdRType: Maximum absolute thickness change in this iteration.         *
****************************************************************************/
static double MMAObjective(unsigned n,
                           const double *x,
                           double *grad,
                           void *Data)
{
    OptMMADataStruct *MMAData = (OptMMADataStruct *) Data;
    double Value;
    unsigned I;

    Value = MMAData->ObjAtPrev;

    for (I = 0; I < n; ++I) {
        if (grad != NULL)
            grad[I] = MMAData->ObjGrad[I];

        Value += MMAData->ObjGrad[I] * (x[I] - MMAData->PrevX[I]);
    }

    return Value;
}

static double MMAVolumeConstraint(unsigned n,
                                  const double *x,
                                  double *grad,
                                  void *Data)
{
    OptMMADataStruct *MMAData = (OptMMADataStruct *) Data;
    double Value;
    unsigned I;

    Value = -MMAData->TargetVolume;

    for (I = 0; I < n; ++I) {
        if (grad != NULL)
            grad[I] = MMAData->VolGrad[I];

        Value += MMAData->VolGrad[I] * x[I];
    }

    return Value;
}

static CagdRType DesignUpdate_MMA(OptWorkflowStateStruct* Opt,
    int Count,
    const CagdRType* Prev,
    int Iter)
{
    const CagdRType Eps = 1.0e-12;
    const CagdRType TMin = Opt->Config.MinimumThickness;
    const CagdRType TMax = OPT_MAX_TILE_THICKNESS;
    const CagdRType MoveFrac = OptDesignMoveFraction(Iter);
    int I, J, K, L, Counter;
    CagdRType CurrentVolume = 0.0;
    CagdRType VolumeUnitScale = 1.0;
    CagdRType CurrentCompliance = 0.0;
    CagdRType MaxThicknessChange = 0.0;
    CagdRType SedAvg = 0.0;
    CagdRType SedSqAvg = 0.0;
    int SedCount = 0;
    double MinObj = 0.0;
    double *X, *LB, *UB;
    CagdRType *ObjGrad, *VolGrad;
    nlopt_opt MMAOpt;
    nlopt_result Result;
    OptMMADataStruct MMAData;

    if (Opt == NULL || Prev == NULL || Count <= 0)
        return 0.0;

    X = (double *) OptAllocZero((size_t) Count, sizeof(double));
    LB = (double *) OptAllocZero((size_t) Count, sizeof(double));
    UB = (double *) OptAllocZero((size_t) Count, sizeof(double));
    ObjGrad = (CagdRType *) OptAllocZero((size_t) Count, sizeof(CagdRType));
    VolGrad = (CagdRType *) OptAllocZero((size_t) Count, sizeof(CagdRType));

    if (X == NULL || LB == NULL || UB == NULL || ObjGrad == NULL || VolGrad == NULL)
        OptFatal("Failed to allocate MMA scratch arrays");

    Counter = 0;
    for (I = 0; I < Opt->TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt->TileStructure[I];

        for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
            OptSubTileStruct *SubTile = &Tile->SubTiles[L];
            CagdRType TOld = Iter <= 1 || Prev[Counter] <= TMin ?
                             Opt->Run.InitialThickness : Prev[Counter];
            CagdRType TLow = TOld * (1.0 - MoveFrac);
            CagdRType THigh = TOld * (1.0 + MoveFrac);
            CagdRType TSafe = TOld > Eps ? TOld : Eps;

            if (TLow < TMin)
                TLow = TMin;
            if (THigh > TMax)
                THigh = TMax;

            X[Counter] = TOld;
            LB[Counter] = TLow;
            UB[Counter] = THigh;

            CurrentVolume += SubTile->Volume;
            CurrentCompliance += SubTile->ComplianceContribution;

            ObjGrad[Counter] = -SubTile->ComplianceContribution / TSafe;
            VolGrad[Counter] = SubTile->Volume / TSafe;

            ++Counter;
        }
    }



    /*
     * Element volumes come from ANSYS and can be expressed in a different
     * unit system than the IRIT geometry. Preserve their relative
     * distribution for the sensitivity, but normalize their sum to the
     * exact current lattice volume used by the volume-fraction target.
     */
    if (CurrentVolume > Eps && Opt->LatticeVol > Eps) {
        VolumeUnitScale = Opt->LatticeVol / CurrentVolume;
        for (I = 0; I < Count; ++I)
            VolGrad[I] *= VolumeUnitScale;
        CurrentVolume = Opt->LatticeVol;
    }

    MMAData.Count = Count;
    MMAData.PrevX = (CagdRType *) Prev;
    MMAData.ObjGrad = ObjGrad;
    MMAData.VolGrad = VolGrad;
    MMAData.ObjAtPrev = CurrentCompliance;
    MMAData.TargetVolume = Opt->Run.TargetVolumeFraction * Opt->DefMapVol;

    if (MMAData.TargetVolume <= Eps)
        MMAData.TargetVolume = CurrentVolume;

    MMAOpt = nlopt_create(NLOPT_LD_MMA, (unsigned) Count);
    if (MMAOpt == NULL)
        OptFatal("Failed to create NLopt MMA optimizer");

    nlopt_set_lower_bounds(MMAOpt, LB);
    nlopt_set_upper_bounds(MMAOpt, UB);
    nlopt_set_min_objective(MMAOpt, MMAObjective, &MMAData);
    nlopt_add_inequality_constraint(MMAOpt,
                                    MMAVolumeConstraint,
                                    &MMAData,
                                    1.0e-8);
    nlopt_set_xtol_rel(MMAOpt, 1.0e-4);
    nlopt_set_maxeval(MMAOpt, 50);

    Result = nlopt_optimize(MMAOpt, X, &MinObj);

    if (Result < 0) {
        fprintf(stderr,
                "Warning: NLopt MMA failed with result %d. Keeping previous thicknesses.\n",
                (int) Result);
        for (I = 0; I < Count; ++I)
            X[I] = Prev[I];
    }

    Counter = 0;
    for (I = 0; I < Opt->TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt->TileStructure[I];

        for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
            CagdRType TNew = X[Counter];
            CagdRType Change = fabs(TNew - Prev[Counter]);

            if (TNew < TMin)
                TNew = TMin;
            if (TNew > TMax)
                TNew = TMax;

            Tile->SubTiles[L].Thickness = TNew;

            if (Change > MaxThicknessChange)
                MaxThicknessChange = Change;

            if (TNew > TMin && Tile->SubTiles[L].SED > 0.0) {
                SedAvg += Tile->SubTiles[L].SED;
                SedSqAvg += Tile->SubTiles[L].SED * Tile->SubTiles[L].SED;
                ++SedCount;
            }

            ++Counter;
        }
    }



    if (SedCount > 0) {
        SedAvg /= SedCount;
        SedSqAvg /= SedCount;
        Opt->SedAvg = SedAvg;
        Opt->SedR = sqrt(fabs(SedSqAvg - SedAvg * SedAvg));
    }
    else {
        Opt->SedAvg = 0.0;
        Opt->SedR = 0.0;
    }

    fprintf(stderr,
            "DesignUpdate_MMA: result=%d, compliance=%.8g, volume=%.8g, target=%.8g, unit-scale=%.8g, change=%.8g, SedAvg=%.8g, SedStd=%.8g\n",
            (int) Result,
            CurrentCompliance,
            CurrentVolume,
            MMAData.TargetVolume,
            VolumeUnitScale,
            MaxThicknessChange,
            Opt->SedAvg,
            Opt->SedR);

    nlopt_destroy(MMAOpt);
    IritFree(X);
    IritFree(LB);
    IritFree(UB);
    IritFree(ObjGrad);
    IritFree(VolGrad);

    return MaxThicknessChange;
}


/***************************************************************************
* DESCRIPTION:                                                             *
*   Adjusts initial lattice thickness until the lattice volume fraction    *
* is within the requested target tolerance. It generates the lattice,      *
* computes the lattice/deformation-map volume ratio, and updates the       *
* initial thickness accordingly.                                           *
*                                                                          *
* PARAMETERS:                                                              *
*   None.                                                                  *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the target volume fraction was reached, FALSE otherwise.  *
***************************************************************************/

static int LatticeBy_VolumeFraction(OptWorkflowStateStruct *Opt)
{
    int Iter;
    const int MaxIter = 50;

    for (Iter = 0; Iter < MaxIter; Iter++) {
        if (!RunIRITIntegrated(Opt, 1)) return 0;

        if (Opt->DefMapVol <= 0.0) {
            fprintf(stderr, "ERROR: Invalid DefMap volume: %.16g\n",
                Opt->DefMapVol);
            return 0;
        }
        {
            CagdRType VolFraction;
            CagdRType Target;
            CagdRType Tol;
            CagdRType LowerLimit;
            CagdRType UpperLimit;

            VolFraction = Opt->LatticeVol / Opt->DefMapVol;
            Target = Opt->Run.TargetVolumeFraction;
            Tol = Opt->Run.TargetVolumeFractionTol;

            LowerLimit = Target * (1.0 - Tol);
            UpperLimit = Target * (1.0 + Tol);

            if (VolFraction > UpperLimit) {
                fprintf(stdout,
                    "Finished iteration %d, volume ratio is: %.8g. "
                    "Target range is [%.8g, %.8g]. Redoing with smaller thickness.\n",
                    Iter + 1, VolFraction, LowerLimit, UpperLimit);

                Opt->Run.InitialThickness *= 0.9;
                fprintf(stdout, "New initial thickness: %.8g\n",
                    Opt->Run.InitialThickness);

                fflush(stdout);

                InitLegacyDataFile(Opt);
            }

            else if (VolFraction < LowerLimit) {
                fprintf(stdout,
                    "Finished iteration %d, volume ratio is: %.8g. "
                    "Target range is [%.8g, %.8g]. Redoing with larger thickness.\n",
                    Iter + 1, VolFraction, LowerLimit, UpperLimit);

                Opt->Run.InitialThickness *= 1.1;
                fprintf(stdout, "New initial thickness: %.8g\n",
                    Opt->Run.InitialThickness);

                fflush(stdout);

                InitLegacyDataFile(Opt);
            }

            else {
                fprintf(stdout,
                    "Finished iteration %d, volume ratio is: %.8g. "
                    "Target range is [%.8g, %.8g]. Finished.\n",
                    Iter + 1, VolFraction, LowerLimit, UpperLimit);

                fflush(stdout);

                return 1;
            }
        }
    }

    fprintf(stderr,
        "ERROR: Failed to reach target volume fraction after %d attempts.\n",
        MaxIter);

    return 0;
}
/***************************************************************************
* DESCRIPTION:                                                             *
*   Calculates the volume of an obejct, including a Trivar list.           *
* It then returns the total volume of the object and all its sub-objects.  *
*                                                                          *
* PARAMETERS:                                                              *
*   PObj - The object to compute the volume for.                           *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static CagdRType ComputeObjVol(IritPrsrObjectStruct* PObj)
{
    CagdRType Vol = 0.0;

    if (PObj == NULL) {
        return 0.0;
    }

    if (IRIT_PRSR_IS_OLST_OBJ(PObj)) {
        int i;
        IritPrsrObjectStruct* PObjTmp;

        for (i = 0; (PObjTmp = IritPrsrListObjectGet(PObj, i)) != NULL; i++) {
            Vol += ComputeObjVol(PObjTmp);
        }
    }
    else if (IRIT_PRSR_IS_TRIVAR_OBJ(PObj)) {
        IritTrivTVStruct* TV;

        for (TV = PObj->U.Trivars; TV != NULL; TV = TV->Pnext) {
            Vol += fabs(IritTrivTVVolume(TV, TRUE));
        }
    }

    return Vol;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Enforces equal thickness values on coincident tile interfaces.         *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static void NeighbourMatching(OptWorkflowStateStruct *Opt)
{
    static const CagdRType Local0Legacy[6][3] = {
        { 0.5, 0.5, 1.0 },
        { 0.5, 0.5, 0.0 },
        { 0.5, 1.0, 0.5 },
        { 0.5, 0.0, 0.5 },
        { 0.0, 0.5, 0.5 },
        { 1.0, 0.5, 0.5 }
    };

    /* Same order used by PreProcessTile for Cross 1:
     * W+, W-, V+, V-, U+, U-. */
    static const CagdRType Local0Verified[6][3] = {
        { 0.5, 0.5, 1.0 },
        { 0.5, 0.5, 0.0 },
        { 0.5, 1.0, 0.5 },
        { 0.5, 0.0, 0.5 },
        { 1.0, 0.5, 0.5 },
        { 0.0, 0.5, 0.5 }
    };

    static const CagdRType Local0_1[6][3] = {
    {0.5, 0.5, 0.0},  /* W- */
    {0.5, 0.5, 1.0},   /* W+ */
    {0.5, 0.0, 0.5},  /* V- */
    {0.5, 1.0, 0.5},  /* V+ */
    {0.0, 0.5, 0.5},  /* U- */
    {1.0, 0.5, 0.5}  /* U+ */
    };

    static const CagdRType Local2[8][3] = {
        { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 },
        { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }
    };

    switch (Opt->Run.LatticeType) {
    case 0:
        NeighbourMatchingB(Opt,
            Opt->Config.Tile0OwnerMapping == 1 ?
                Local0Verified : Local0Legacy,
            6, 0);
        break;
    case 1:
        NeighbourMatchingB(Opt, Local0_1, 6, 0);
        break;
    case 2:
        NeighbourMatchingB(Opt, Local2, 8, 8);
        break;
    default:
        NeighbourMatchingB(Opt, Local0Legacy, 6, 0);
        break;
    }
}

typedef struct OptNeighbourEntryStruct {
    int LinearId;
    CagdRType X, Y, Z;
} OptNeighbourEntryStruct;

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int CompareNeighbourEntries(const void *P1, const void *P2)
{
    const OptNeighbourEntryStruct* A = (const OptNeighbourEntryStruct*)P1;
    const OptNeighbourEntryStruct* B = (const OptNeighbourEntryStruct*)P2;

    if (A -> X < B -> X) return -1;
    if (A -> X > B -> X) return  1;
    if (A -> Y < B -> Y) return -1;
    if (A -> Y > B -> Y) return  1;
    if (A -> Z < B -> Z) return -1;
    if (A -> Z > B -> Z) return  1;
    return 0;
}


/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void NeighbourMatchingB(OptWorkflowStateStruct *Opt,
                         const CagdRType Local[][3], int LocalCount, int LDefOffset)
{
    int B, U, V, W, L, X = 0, I, J, G;
    int Total = Opt->TotalTiles * LocalCount;
    OptNeighbourEntryStruct *Entries = (OptNeighbourEntryStruct *)
        OptAllocZero(Total, sizeof(OptNeighbourEntryStruct));
    if (!Entries) OptFatal("Cannot allocate neighbour data");
    /* First match exact integer-grid interfaces separately within each body. */
    for (B = 0; B < Opt->BodyCount; ++B) {
        OptBodyStruct *Body = &Opt->Bodies[B];
        int Start = X;
        for (W = 0; W < Body->Counts[2]; ++W)
            for (V = 0; V < Body->Counts[1]; ++V)
                for (U = 0; U < Body->Counts[0]; ++U)
                    for (L = 0; L < LocalCount; ++L) {
                        int Tile = Body->TileOffset + U + Body->Counts[0] * (V + Body->Counts[1] * W);
                        Entries[X].LinearId = Tile * LocalCount + L;
                        Entries[X].X = U + Local[L][0];
                        Entries[X].Y = V + Local[L][1];
                        Entries[X].Z = W + Local[L][2]; ++X;
                    }
        qsort(Entries + Start, X - Start, sizeof(*Entries), CompareNeighbourEntries);
        for (I = Start; I < X; I = J) {
            CagdRType Avg = 0.0;
            for (J = I+1; J < X && CompareNeighbourEntries(&Entries[I], &Entries[J]) == 0; ++J);
            for (G = I; G < J; ++G) {
                int Id = Entries[G].LinearId;
                Avg += Opt->TileStructure[Id / LocalCount].SubTiles[Id % LocalCount + LDefOffset].Thickness;
            }
            Avg /= J-I;
            for (G = I; G < J; ++G) {
                int Id = Entries[G].LinearId;
                Opt->TileStructure[Id / LocalCount].SubTiles[Id % LocalCount + LDefOffset].Thickness = Avg;
            }
        }
    }
    IritFree(Entries);
    OptMatchBodyBoundaries(Opt, Local, LocalCount, LDefOffset);
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Writes the legacy thickness CSV file consumed by the geometry callback *
* code.                                                                    *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static void CreateLegacyDataFile(OptWorkflowStateStruct *Opt)
{
    FILE *F;
    int B, U, V, W, L;
    int CounterFile = 0;

    F = fopen(Opt->Config.DataNewFile, "w");
    if (F == NULL)
        OptFatalF("Failed to create data file: %s", Opt->Config.DataNewFile);

    OptWriteCsvHeader(F, Opt);

    /* Same ordering as DesignUpdate_ASAF: body -> W -> V -> U -> SubTile. */
    for (B = 0; B < Opt->BodyCount; ++B) {
        OptBodyStruct *Body = &Opt->Bodies[B];

        for (W = 0; W < Body->Counts[2]; ++W) {
            for (V = 0; V < Body->Counts[1]; ++V) {
                for (U = 0; U < Body->Counts[0]; ++U) {
                    int TileIndex = Body->TileOffset +
                        U + Body->Counts[0] * (V + Body->Counts[1] * W);

                    for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
                        fprintf(F,
                            "%d,%.17g,%.17g\n",
                            CounterFile + 1,
                            (double)Opt->TileStructure[TileIndex].SubTiles[L].Thickness,
                            (double)Opt->TileStructure[TileIndex].SubTiles[L].SED);
                        ++CounterFile;
                    }
                }
            }
        }
    }

    if (CounterFile != Opt->TotalTiles * Opt->LMin) {
        fclose(F);
        OptFatal("CreateLegacyDataFile thickness count mismatch");
    }

    fclose(F);
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Verifies the exact generated lattice volume after a design update.      *
* If the candidate violates the target volume fraction, thickness values    *
* are scaled down, neighbor matching is reapplied, and the geometry volume  *
* is recomputed before the next ANSYS solve is allowed to use the design.   *
*                                                                          *
* PARAMETERS:                                                              *
*   None                                                                   *
*                                                                          *
* RETURN VALUE:                                                            *
*   CagdRType: Maximal thickness correction applied by this volume gate.    *
***************************************************************************/

static CagdRType OptDesignMoveFraction(int Iter)
{
    return Iter <= 1 ? 0.05 : 0.15;
}

static CagdRType EnforceExactVolumeConstraint(OptWorkflowStateStruct *Opt,
                                                const CagdRType *Prev,
                                                int Count,
                                                int Iter)
{
    const int MaxIter = 12;
    const CagdRType Eps = 1.0e-12;
    const CagdRType TMin = Opt->Config.MinimumThickness;
    const CagdRType MoveFrac = OptDesignMoveFraction(Iter);
    CagdRType VolAtPrev;
    CagdRType TargetVolume;
    CagdRType MaxCorrection = 0.0;
    int It;

    if (Opt == NULL || Prev == NULL || Count <= 0 ||
        Opt->Run.TargetVolumeFraction <= 0.0)
        return 0.0;

    TargetVolume = Opt->Run.TargetVolumeFraction * Opt->DefMapVol;
    if (TargetVolume <= Eps)
        return 0.0;

    for (It = 0; It < MaxIter; ++It) {
        int I, J, K, L;
        int Counter = 0;
        int LimitedCount = 0;
        CagdRType Ratio;
        CagdRType Scale;

        CreateLegacyDataFile(Opt);

        if (!RunIRITIntegrated(Opt, 1))
            OptFatal("Failed to recompute geometry volume after design update");

        TargetVolume = Opt->Run.TargetVolumeFraction * Opt->DefMapVol;
        if (TargetVolume <= Eps)
            return MaxCorrection;

        Ratio = TargetVolume / Opt->LatticeVol;
        if (fabs(Ratio - 1.0) <= 1.0e-4) {
            fprintf(stderr,
                    "Exact volume check: lattice=%.8g, target=%.8g, correction=%.8g\n",
                    Opt->LatticeVol,
                    TargetVolume,
                    MaxCorrection);
            return MaxCorrection;
        }

        Scale = sqrt(Ratio);

        if (Scale <= 0.0)
            Scale = Ratio < 1.0 ? 0.95 : 1.05;

        for (I = 0; I < Opt->TotalTiles; ++I) {
            OptTileStruct *Tile = &Opt->TileStructure[I];

            for (L = Opt->LDef; L < Opt->LDef + Opt->LMin; ++L) {
                CagdRType OldT = Tile->SubTiles[L].Thickness;
                CagdRType PrevT = Iter <= 1 || Counter >= Count || Prev[Counter] <= TMin ?
                                  Opt->Run.InitialThickness : Prev[Counter];
                CagdRType TLow = PrevT * (1.0 - MoveFrac);
                CagdRType THigh = PrevT * (1.0 + MoveFrac);
                CagdRType NewT = OldT * Scale;
                CagdRType Delta;

                if (TLow < TMin)
                    TLow = TMin;
                if (THigh > OPT_MAX_TILE_THICKNESS)
                    THigh = OPT_MAX_TILE_THICKNESS;

                if (NewT < TLow) {
                    NewT = TLow;
                    ++LimitedCount;
                }
                else if (NewT > THigh) {
                    NewT = THigh;
                    ++LimitedCount;
                }

                Tile->SubTiles[L].Thickness = NewT;
                Delta = fabs(NewT - OldT);
                if (Delta > MaxCorrection)
                    MaxCorrection = Delta;

                ++Counter;
            }
        }



        NeighbourMatching(Opt);

        fprintf(stderr,
                "Exact volume check: lattice=%.8g target=%.8g; scaling thicknesses by %.8g, move-limited=%d\n",
                Opt->LatticeVol,
                TargetVolume,
                Scale,
                LimitedCount);

        if (LimitedCount >= Count) {
            fprintf(stderr,
                    "Exact volume check stopped at move limit before reaching target. Next iteration will continue correcting.\n");
            CreateLegacyDataFile(Opt);
            return MaxCorrection;
        }
    }

    CreateLegacyDataFile(Opt);
    if (!RunIRITIntegrated(Opt, 1))
        OptFatal("Failed to recompute geometry volume after final volume correction");

    if (fabs(Opt->LatticeVol / TargetVolume - 1.0) > 1.0e-4) {
        fprintf(stderr,
                "Warning: exact volume gate stopped away from target: lattice=%.8g, target=%.8g\n",
                Opt->LatticeVol,
                TargetVolume);
    }

    return MaxCorrection;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reads the legacy geometry CSV file and extracts its numeric content.   *
*                                                                          *
* PARAMETERS:                                                              *
*   FileName: Input legacy CSV file.                                       *
*   OutRows:  Location to store the number of data rows.                   *
*   OutCols:  Location to store the number of data columns.                *
*                                                                          *
* RETURN VALUE:                                                            *
*   float **: Matrix of parsed values, or NULL on failure.                 *
***************************************************************************/

static float **ReadLegacyCsvToFloatArray(OptWorkflowStateStruct *Opt,
                         const char *FileName, int *OutRows, int *OutCols)
{
    int I, N = Opt->TotalTiles * Opt->LMin;
    CagdRType *Values = ReadTileData(Opt, FileName, N);
    float **Data;
    if (!Values) return NULL;
    Data = (float **) OptAllocZero(N, sizeof(float *));
    if (!Data) { IritFree(Values); return NULL; }
    for (I = 0; I < N; ++I) {
        Data[I] = (float *) OptAllocZero(3, sizeof(float));
        if (!Data[I]) { FreeLegacyCsv(Data, I); IritFree(Values); return NULL; }
        Data[I][0] = (float) (I + 1);
        Data[I][1] = (float) Values[I];
    }
    IritFree(Values);
    Opt->Legacy.LType = Opt->Run.LatticeType;
    Opt->Legacy.MeshSample = Opt->Run.MeshSample;
    Opt->Legacy.MeshType = Opt->Run.MeshType;
    *OutRows = N; *OutCols = 3;
    return Data;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void FreeLegacyCsv(float **Data, int Rows)
{
    int I;
    if (Data == NULL)
        return;

    for (I = 0; I < Rows; ++I) {
        if (Data[I] != NULL)
            IritFree(Data[I]);
    }
    IritFree(Data);
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static int LoadLegacyGeometryData(OptWorkflowStateStruct *Opt, const char *FileName)
{
    Opt->Legacy.Data = ReadLegacyCsvToFloatArray(Opt, FileName, &Opt->Legacy.Rows, &Opt->Legacy.Cols);
    return Opt->Legacy.Data != NULL;
}

/***************************************************************************
* AUXILIARY:                                                               *
* Auxiliary function in the optimization workflow.                         *
***************************************************************************/

static void UnloadLegacyGeometryData(OptWorkflowStateStruct *Opt)
{
    FreeLegacyCsv(Opt->Legacy.Data, Opt->Legacy.Rows);
    Opt->Legacy.Data = NULL;
    Opt->Legacy.Rows = 0;
    Opt->Legacy.Cols = 0;
}

/***************************************************************************
 *                         IRIT CALLBACK CODE                               *
 ***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
*   Prepares the boundary parameters of one face of a synthesized tile.    *
*                                                                          *
* PARAMETERS:                                                              *
*   CBData: Callback data of the current tile.                             *
*   BPrm:   Boundary-parameter record to populate.                         *
*   Circ:   TRUE for circular cross sections, FALSE otherwise.             *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the parameter structure was filled successfully.          *
***************************************************************************/
static int PreProcessTile1FaceParam(
                    IritUserMicroTileBndryPrmStruct* BPrm,
                    CagdRType OuterRadius,
                    CagdBType Circ)
{
    IRIT_ZAP_MEM(BPrm, sizeof(IritUserMicroTileBndryPrmStruct));

    BPrm -> BndryShape = 0.0;
    BPrm -> CntrShape = 0.4;
    BPrm -> Circular = Circ;
    BPrm -> InnerRadius = 0.0;
    BPrm -> OuterRadius = OuterRadius;

    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Builds one tile geometry according to the current legacy thickness     *
* data and tile type.                                                      *
*                                                                          *
* PARAMETERS:                                                              *
*   Tile:   Tile geometry to preprocess. Expected to be NULL here.         *
*   CBData: Callback data of the current tile.                             *
*                                                                          *
* RETURN VALUE:                                                            *
*   IritPrsrObjectStruct *: Constructed tile object, or NULL on failure.   *
***************************************************************************/
static IritPrsrObjectStruct *PreProcessTile(
                                  IritPrsrObjectStruct* Tile,
                                  IritUserMicroPreProcessTileCBStruct* CBData)
{
    UserMicroLocalDataStruct* LclData = (UserMicroLocalDataStruct*)CBData -> CBFuncData;
    OptWorkflowStateStruct *Opt = LclData -> Opt;
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax,
        * LclMinDmn = CBData -> TileLclDmnMin,
        * LclMaxDmn = CBData -> TileLclDmnMax;
    IritGeomBBBboxStruct BBox;
    IritUserMicroTileBndryPrmStruct UMinPrms, UMaxPrms, VMinPrms, VMaxPrms, WMinPrms, WMaxPrms;
    int Index;
    int BaseIndex;
    int Arms = 6;
    CagdBType Circle = FALSE;
    //CagdBType RemoveThinTiles = LclData -> RemoveThinTiles;

    assert(Tile == NULL);
    IritTrivTVDomain(LclData -> DefMap, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);

    Index = (CBData -> TileIdxs[0]) +
        (CBData -> TileIdxs[1] * LclData->Body->Counts[0]) +
        (CBData -> TileIdxs[2] * LclData->Body->Counts[0] * LclData->Body->Counts[1]);

    fprintf(stderr,
        "Index:[%d], Tile[%d,%d,%d] from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.3f)\r",
        Index,
        CBData -> TileIdxs[0], CBData -> TileIdxs[1], CBData -> TileIdxs[2],
        LclMinDmn[0], LclMinDmn[1], LclMinDmn[2],
        LclMaxDmn[0], LclMaxDmn[1], LclMaxDmn[2]);

    switch (Opt->Legacy.LType) {
    case 0: /* Cross 1 */
    case 1: /* Cross 2 */
        Arms = 6;
        Circle = TRUE;
        break;
    case 2: /* Diagonal */
        Arms = 8;
        Circle = FALSE;
        break;
    case 3: /* Combined */
        Arms = 14;
        Circle = FALSE;
        break;
    default:
        Arms = 6;
        Circle = FALSE;
        break;
    }

    if (Index < 0 || Index >= LclData->Body->TileCount)
        OptFatal("Tile callback index outside body grid");
    BaseIndex = (LclData->Body->TileOffset + Index) * Arms;
    if (BaseIndex + Arms > Opt->Legacy.Rows)
        OptFatal("Tile callback exceeds thickness data");

    if (Opt->Config.RemoveThinTiles) {
        CagdRType Avg = 0.0;
        int I;
        for (I = 0; I < Arms; ++I)
            Avg += Opt->Legacy.Data[BaseIndex + I][1];
        Avg /= Arms;

        if (Avg < OPT_MIN_TILE_THICKNESS) {
            fprintf(stderr, "Warning: Removing thin tile at index %d with average thickness %.6f\n",
                Index, (double)Avg);
            return NULL;
        }
    }

    switch (Opt->Legacy.LType) {
    case 0:
    {
        if (!PreProcessTile1FaceParam(&WMaxPrms, Opt->Legacy.Data[BaseIndex + 0][1], Circle) ||
            !PreProcessTile1FaceParam(&WMinPrms, Opt->Legacy.Data[BaseIndex + 1][1], Circle) ||
            !PreProcessTile1FaceParam(&VMaxPrms, Opt->Legacy.Data[BaseIndex + 2][1], Circle) ||
            !PreProcessTile1FaceParam(&VMinPrms, Opt->Legacy.Data[BaseIndex + 3][1], Circle) ||
            !PreProcessTile1FaceParam(&UMaxPrms, Opt->Legacy.Data[BaseIndex + 4][1], Circle) ||
            !PreProcessTile1FaceParam(&UMinPrms, Opt->Legacy.Data[BaseIndex + 5][1], Circle)) {
            return NULL;
        }

        Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
            &VMinPrms, &VMaxPrms,
            &WMinPrms, &WMaxPrms,
            FALSE, NULL);

        IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
        if (BBox.Min[0] < -IRIT_UEPS || BBox.Max[0] > 1.0 + IRIT_UEPS ||
            BBox.Min[1] < -IRIT_UEPS || BBox.Max[1] > 1.0 + IRIT_UEPS ||
            BBox.Min[2] < -IRIT_UEPS || BBox.Max[2] > 1.0 + IRIT_UEPS) {
            fprintf(stderr, "Warning: Tile spans beyond the unit box.\n");
        }

        Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);
    }
    break;

    case 1:
    {
        static const int Map[6] = { 4, 5, 2, 3, 0, 1 };
        IrtRType OuterRadii[6];
        int i;
        CagdRType Avg = 0.0;

        for (i = 0; i < 6; ++i) {
            OuterRadii[i] = Opt->Legacy.Data[BaseIndex + Map[i]][1];
            Avg += OuterRadii[i];
        }
        Avg /= 6.0;

        Tile = IritUserMicro3DCrossTile2(FALSE,
            OuterRadii,
            NULL,
            2.5 * Avg,
            NULL);
        IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
        if (BBox.Min[0] < -IRIT_UEPS || BBox.Max[0] > 1.0 + IRIT_UEPS ||
            BBox.Min[1] < -IRIT_UEPS || BBox.Max[1] > 1.0 + IRIT_UEPS ||
            BBox.Min[2] < -IRIT_UEPS || BBox.Max[2] > 1.0 + IRIT_UEPS) {
            fprintf(stderr, "Warning: Tile spans beyond the unit box.\n");
        }
        Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);
    }
    break;

    case 2:
    {
        CagdRType CntrSize = 0.1;
        CagdRType CntrVertScl = 1.0;
        const CagdRType CrnrSizes[8] = {
            Opt->Legacy.Data[BaseIndex + 0][1],     Opt->Legacy.Data[BaseIndex + 1][1],
            Opt->Legacy.Data[BaseIndex + 2][1], Opt->Legacy.Data[BaseIndex + 3][1],
            Opt->Legacy.Data[BaseIndex + 4][1], Opt->Legacy.Data[BaseIndex + 5][1],
            Opt->Legacy.Data[BaseIndex + 6][1], Opt->Legacy.Data[BaseIndex + 7][1]
        };
        CagdRType CrnrVertScl[2] = { 0.4, 0.4 };
        const CagdRType *VertBars = NULL;
        int RefineVertRods2x2 = 1;
        CagdRType SmoothFactor = 0.2;

        (void)CntrVertScl;
        (void)CrnrVertScl;
        (void)VertBars;
        (void)RefineVertRods2x2;
        (void)SmoothFactor;

        Tile = IritUserMicroDiagTile1(CntrSize, 1.0,
            CrnrSizes, CrnrVertScl,
            NULL, 1, 0.2,
            NULL, NULL);
        IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
        if (BBox.Min[0] < -IRIT_UEPS || BBox.Max[0] > 1.0 + IRIT_UEPS ||
            BBox.Min[1] < -IRIT_UEPS || BBox.Max[1] > 1.0 + IRIT_UEPS ||
            BBox.Min[2] < -IRIT_UEPS || BBox.Max[2] > 1.0 + IRIT_UEPS) {
            fprintf(stderr, "Warning: Tile spans beyond the unit box.\n");
        }
        Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);
    }
    break;

    case 3:
    default:
        /* Combined tile construction is not implemented in this module. */
        break;
    }

    return Tile;
}

/***************************************************************************
* Adapt FFD callback data to the existing per-tile geometry builder.        *
***************************************************************************/
static IritPrsrObjectStruct *PreProcessTileFFD(
    IritPrsrObjectStruct *Tile,
    IritTrivFFDInfoStruct *FFDInfo)
{
    UserMicroLocalDataStruct *LclData =
        (UserMicroLocalDataStruct *) FFDInfo->_AuxInfo;
    IritUserMicroPreProcessTileCBStruct CBData;
    CagdRType Min[3], Max[3];
    int Counts[3], I;

    assert(Tile == NULL);
    IRIT_ZAP_MEM(&CBData, sizeof(CBData));
    CBData.CBFuncData = LclData;

    Counts[0] = LclData->Body->Counts[0];
    Counts[1] = LclData->Body->Counts[1];
    Counts[2] = LclData->Body->Counts[2];

    IritTrivTVDomain(LclData->DefMap,
                    &Min[0], &Max[0],
                    &Min[1], &Max[1],
                    &Min[2], &Max[2]);

    for (I = 0; I < 3; ++I) {
        /* Recover zero-based grid indices for the legacy thickness lookup. */
        CBData.TileIdxs[I] = (int) (
            (FFDInfo->_DmnMin[I] - Min[I]) *
            Counts[I] / (Max[I] - Min[I]) + 0.5);
        CBData.NumOfTiles[I] = Counts[I];
        CBData.TileLclDmnMin[I] =
            (FFDInfo->_DmnMin[I] - Min[I]) / (Max[I] - Min[I]);
        CBData.TileLclDmnMax[I] =
            (FFDInfo->_DmnMax[I] - Min[I]) / (Max[I] - Min[I]);
    }

    /* The existing builder scales the tile; FFD translates it afterwards. */
    IRIT_HMGN_MAT_COPY(CBData.Mat, FFDInfo->_SclMat);
    return PreProcessTile(NULL, &CBData);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Builds the heterogeneous micro-structure object for the current          *
* thickness state.  The function loads the deformation map, computes the     *
* deformation-map and lattice volumes, and optionally exports ITD and INP    *
* outputs.  Tile aggregation supports original composition, direct FFD,     *
* and composition with global placement and FFD approximation, selected     *
* using UseFFDTileAggregation and UseMicroStructFFD.                         *
*                                                                            *
* PARAMETERS:                                                                *
*   Opt:  Pointer to the active optimization workflow state.                 *
*   mode: Geometry-generation mode:                                          *
*         0 - generate geometry, compute volumes, and export ITD/INP.        *
*         1 - generate geometry and compute volumes only.                    *
*         2 - compute deformation-map volume only.                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: TRUE if the requested geometry-generation stage completed           *
*        successfully, FALSE otherwise.                                      *
*****************************************************************************/
static int GenerateMicroStructures(OptWorkflowStateStruct *Opt, int mode)
{
    /* UseMicroStructFFD takes priority; otherwise select original or direct FFD. */
    const int UseFFDTileAggregation = 0;
    const int UseMicroStructFFD = 1;
    IritPrsrObjectStruct *MS = NULL, *Scaled = NULL, *BodyMS = NULL;
    IritUserMicroParamStruct Param;
    IritUserMicroRegularParamStruct *Reg = &Param.U.RegularParam;
    IritMvarMVStruct *MV = NULL;
    CagdRType Steps[3][2];
    IritTrivFFDInfoStruct FFDInfo;
    IritPrsrINPSaveDfltFileParamsStruct Save;
    IrtHmgnMatType Mat;
    UserMicroLocalDataStruct LclData;
    CagdRType ScaleVol = Opt->Config.Scale * Opt->Config.Scale * Opt->Config.Scale;
    int B, D, Handler, Ok = FALSE;
    Opt->DefMapVol = 0.0;
    for (B = 0; B < Opt->BodyCount; ++B)
        Opt->DefMapVol += fabs(IritTrivTVVolume(Opt->Bodies[B].Map, TRUE)) * ScaleVol;
    if (mode == 2) return TRUE;
    MS = IritPrsrGenLISTObject(NULL);
    if (mode == 0) Opt->MeshOwnerCount = 0;
    for (B = 0; B < Opt->BodyCount; ++B) {
        OptBodyStruct *Body = &Opt->Bodies[B];
        IRIT_ZAP_MEM(&LclData, sizeof(LclData));
        LclData.Opt = Opt; LclData.Body = Body; LclData.DefMap = Body->Map;
        LclData.BndryThickness = 3.0;
        fprintf(stderr, "\nGenerating body %d/%d (%s), %d tiles\n", B+1, Opt->BodyCount, Body->Name, Body->TileCount);
        if (UseMicroStructFFD || !UseFFDTileAggregation) {
            IRIT_ZAP_MEM(&Param, sizeof(Param));
            MV = IritMvarCnvrtTVToMV(Body->Map);
            if (!MV) goto Cleanup;
            Param.TilingType = IRIT_USER_MICRO_TILE_REGULAR;
            Param.DeformMV = MV; Param.ApproxLowOrder = 2;
            Reg->TilingStepMode = TRUE;
            Reg->MaxPolyEdgeLen = 0.1;
            Reg->GlobalPlacement = UseMicroStructFFD ? TRUE : FALSE;
            Reg->FFDeformApprox = UseMicroStructFFD ? TRUE : FALSE;
            for (D = 0; D < 3; ++D) {
                Steps[D][0] = Body->Counts[D]; Steps[D][1] = 0;
                Reg->TilingSteps[D].TilesPerIntervals = Steps[D];
                Reg->TilingSteps[D].Len = 1;
            }
            Reg->PreProcessCBFunc = PreProcessTile;
            Reg->PostProcessCBFunc = PostProcessTileBodies;
            Reg->CBFuncData = &LclData;
            BodyMS = IritUserMicroStructComposition(&Param);
            IritMvarMVFree(MV); MV = NULL;
        }
        else {
            IRIT_ZAP_MEM(&FFDInfo, sizeof(FFDInfo));
            FFDInfo.FFDCB = PreProcessTileFFD; FFDInfo._AuxInfo = &LclData;
            BodyMS = IritTrivFFDTileObjectInTV(NULL, Body->Map,
                Body->Counts[0], Body->Counts[1], Body->Counts[2], FALSE,
                0.0, 0.0, NULL, &FFDInfo);
        }
        if (!BodyMS || !OptCollectMesh(Opt, BodyMS, MS, mode == 0)) goto Cleanup;
        IritPrsrFreeObject(BodyMS); BodyMS = NULL;
    }
    Opt->LatticeVol = ComputeObjVol(MS) * ScaleVol;
    if (mode == 0) {
        IritMiscMatGenMatUnifScale(Opt->Config.Scale, Mat);
        Scaled = IritGeomTransformObject(MS, Mat);
        if (!Scaled) goto Cleanup;
        IRIT_ZAP_MEM(&Save, sizeof(Save));
        for (D = 0; D < 3; ++D) Save.UVWSamples[D] = Opt->Run.MeshSample;
        Save.Epsilon = IRIT_EPS;
        Save.BzrOnly = FALSE; /* Already converted; ownership matches each patch. */
        strcpy(Save.RealFormat, "%.17g");
        Save.SaveDepth = 1;
        Save.GeometryType = IRIT_PRSR_FEM_ELEMENT_C3D8;
        Save.Messages = TRUE;
        if (!IritPrsrInpSaveFile(Scaled, Opt->Config.InpFile, &Save)) goto Cleanup;
        Handler = IritPrsrOpenDataFile(Opt->Config.ItdOutputFile, FALSE, 1);
        if (Handler < 0) goto Cleanup;
        IritPrsrPutObjectToHandler(Handler, MS);
        IritPrsrCloseStream(Handler, TRUE);
        if (!OptValidateMesh(Opt)) goto Cleanup;
    }
    Ok = TRUE;
Cleanup:
    if (MV) IritMvarMVFree(MV);
    if (BodyMS) IritPrsrFreeObject(BodyMS);
    if (Scaled) IritPrsrFreeObject(Scaled);
    if (MS) IritPrsrFreeObject(MS);
    return Ok;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Loads the input files and executes the optimization workflow.          *
*                                                                          *
* PARAMETERS:                                                              *
*   argc: Number of command-line arguments.                                *
*   argv: Command-line argument vector.                                    *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: Process exit status.                                              *
***************************************************************************/
static int Optimization_Main(int argc, char **argv)
{
    const char *RunFile = argc > 1 ? argv[1] : "opt_run.txt";
    const char *CfgFile = argc > 2 ? argv[2] : "opt_config.txt";
    OptWorkflowStateStruct OptData;
    OptWorkflowStateStruct *Opt = &OptData;

    IRIT_ZAP_MEM(Opt, sizeof(*Opt));

    if (!LoadRunInputs(RunFile, &Opt->Run))
        OptFatalF("Failed to parse run input file: %s", RunFile);

    if (!LoadConfigInputs(CfgFile, &Opt->Config))
        OptFatalF("Failed to parse config input file: %s", CfgFile);

    if (!ValidateInputs(&Opt->Run, &Opt->Config))
        OptFatal("Input validation failed");

    InitializeWorkflowState(Opt);

    if (argc > 3 && strcmp(argv[3], "--geometry-only") == 0) {
        int Ok = InitLegacyDataFile(Opt) && RunIRITIntegrated(Opt, 0);
        FreeWorkflowState(Opt);
        return Ok ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (!RunWorkflow(Opt))
        OptFatal("Workflow failed");

    FreeWorkflowState(Opt);
    return EXIT_SUCCESS;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Program entry point.                                                   *
*                                                                          *
* PARAMETERS:                                                              *
*   argc: Number of command-line arguments.                                *
*   argv: Command-line argument vector.                                    *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: Process exit status.                                              *
***************************************************************************/
int main(int argc, char **argv)
{
    return Optimization_Main(argc, argv);
}
