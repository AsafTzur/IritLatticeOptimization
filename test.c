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
#include "opt_config.h"

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
#define OPT_INVALID_REAL        (-1.0e300)
#define OPT_MAX_TILE_THICKNESS  0.5
#define OPT_MIN_TILE_THICKNESS  0.001
#define OPT_LEGACY_MESH_TYPE    1

/***************************************************************************
*                        Input File Structures.                            *
***************************************************************************/

/***************************************************************************
*                   INP / Simulation Data Structures.                      *
*                                                                          *
* These structures hold the parsed mesh, simulation values, and the tile   *
* bookkeeping used by the optimization stages.                             *
***************************************************************************/

typedef struct OptNodeStruct {
    int Id;
    CagdRType X, Y, Z;
} OptNodeStruct;

typedef struct OptElementStruct {
    int Id;
    CagdRType SED;
    CagdRType Volume;
    CagdRType StrainEnergy;
    CagdRType ComplianceContribution;
} OptElementStruct;

typedef struct OptSubTileStruct {
    OptElementStruct **Elements;
    int ElementCount;
    int ElementCapacity;

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

typedef CagdRType (*LclThicknessFuncType)(CagdRType u,
                                            CagdRType v,
                                            CagdRType w);
typedef int (*OptInputHandlerFuncType)(const char *Key,
                                       const char *Value,
                                       void *Ctx);

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

typedef struct OptNeighbourEntryStruct {
    int LinearId;
    CagdRType X, Y, Z;
} OptNeighbourEntryStruct;

typedef struct UserMicroLocalDataStruct {
    IritTrivTVStruct *DefMap;
    CagdRType BndryThickness;
    LclThicknessFuncType ThicknessFuncCB;
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

    CagdRType MaxDef;
    CagdRType Mass;
    CagdRType SedAvg;
    CagdRType SedR;
    CagdRType DefMapVol;
    CagdRType LatticeVol;
    CagdRType LastThicknessChange;
    int ConsecutiveSmallChanges;
    int CurrentIteration;
    CagdBType FinishReasonWritten;

    char TimeIritStart[16];
    char TimeIritFinish[16];
    char TimeSimulationFinish[16];
    char FolderName[OPT_MAX_NAME];

    OptNodeArrayStruct Nodes;
    OptElementArrayStruct Elements;
    OptScalarFieldStruct SedData;
    OptScalarFieldStruct VolumeData;

    OptBodyStruct *Bodies;
    int BodyCount, TotalTiles;
    OptMeshOwnerStruct *MeshOwners;
    int MeshOwnerCount, MeshOwnerCapacity;
    OptTileStruct *TileStructure;    /* Persistent tile/subtile topology.
                                       * Analysis fields and FE links are
                                       * refreshed each iteration; thickness
                                       * is retained. */
    CagdRType *Results;              /* Flat array [Iter][3]. */
    OptLegacyStateStruct Legacy;     /* Per-run IRIT CSV/callback state. */
};


#define RESULT_AT(Opt, Iter, Col) \
    ((Opt) -> Results[(Iter) * 3 + (Col)])

/***************************************************************************
*                           Forward Declarations                            *
***************************************************************************/

/* Workflow. */
static int OptimizationMain(int Argc, char **Argv);
int IritOptRun(const OptConfigStruct *Config, CagdBType GeometryOnly);
static int RunWorkflow(OptWorkflowStateStruct *Opt);
static int InitLegacyDataFile(OptWorkflowStateStruct *Opt);
static int RunIRITIntegrated(OptWorkflowStateStruct *Opt, int Mode);
static int RunSimulation(OptWorkflowStateStruct *Opt);
static int RunAnsysBatch(OptWorkflowStateStruct *Opt);
static int RunCalculixBatch(OptWorkflowStateStruct *Opt);
static int ProcessIterationData(OptWorkflowStateStruct *Opt, int Iter);
static int OrganizeResults(OptWorkflowStateStruct *Opt,
                           const char *SourceFolder,
                           const char *TargetBaseFolder,
                           int Iter);
static int OptWriteFinishReason(OptWorkflowStateStruct *Opt,
                                const char *Reason,
                                int Iter);
static int ReadLastRowOfDesignParameters(const char *FileName,
                                         CagdRType *Deflection,
                                         CagdRType *Mass);

/* Optimization. */
static void UpdateTileAnalysisData(OptWorkflowStateStruct *Opt);
static CagdRType DesignUpdateSEDRatio(OptWorkflowStateStruct *Opt,
                                      int Count,
                                      const CagdRType *Prev,
                                      int Iter);
static CagdRType DesignUpdateMMA(OptWorkflowStateStruct *Opt,
                                 int Count,
                                 const CagdRType *Prev,
                                 int Iter);
static double MMAObjective(unsigned n,
                           const double *x,
                           double *Grad,
                           void *Data);
static double MMAVolumeConstraint(unsigned n,
                                  const double *x,
                                  double *Grad,
                                  void *Data);
static CagdRType DesignUpdateOC(OptWorkflowStateStruct *Opt,
                                int Count,
                                const CagdRType *Prev);
static CagdRType DesignUpdateAsaf(OptWorkflowStateStruct *Opt,
                                  int Count,
                                  const CagdRType *Prev,
                                  int Iter);
static CagdRType ComputeMaxThicknessChange(const OptWorkflowStateStruct *Opt,
                                           const CagdRType *Prev,
                                           int Count);
static void NeighbourMatching(OptWorkflowStateStruct *Opt);
static int CompareNeighbourEntries(const void *P1, const void *P2);
static void NeighbourMatchingB(OptWorkflowStateStruct *Opt,
                               const CagdRType Local[][3],
                               int LocalCount,
                               int LDefOffset);
static CagdRType OptDesignMoveFraction(int Iter);
static int LatticeByVolumeFraction(OptWorkflowStateStruct *Opt);
static CagdRType ComputeObjVol(IritPrsrObjectStruct *PObj);
static void CreateLegacyDataFile(OptWorkflowStateStruct *Opt);

/* Geometry generation. */
static int LoadLegacyGeometryData(OptWorkflowStateStruct *Opt,
                                  const char *FileName);
static void UnloadLegacyGeometryData(OptWorkflowStateStruct *Opt);
static float **ReadLegacyCsvToFloatArray(OptWorkflowStateStruct *Opt,
                                         const char *FileName,
                                         int *OutRows,
                                         int *OutCols);
static void FreeLegacyCsv(float **Data, int Rows);
static int PreProcessTile1FaceParam(IritUserMicroTileBndryPrmStruct *BPrm,
                                    CagdRType OuterRadius,
                                    CagdBType Circ,
                                    const CagdRType SkinThickness[4]);
static IritPrsrObjectStruct *PreProcessTile(
                                  IritPrsrObjectStruct *Tile,
                                  IritUserMicroPreProcessTileCBStruct *CBData);
static IritPrsrObjectStruct *PreProcessTileFFD(
                                  IritPrsrObjectStruct *Tile,
                                  IritTrivFFDInfoStruct *FFDInfo);
static int GenerateMicroStructures(OptWorkflowStateStruct *Opt, int Mode);

/* Multi-body support, implemented in multibody_support.inc. */
static int OptParseCounts(const char *Text, int **Values, int *Count);
static int OptCountBodyTVs(const IritPrsrObjectStruct *Obj);
static int OptCopyBodyTVs(OptWorkflowStateStruct *Opt,
                          const IritPrsrObjectStruct *Obj,
                          int *BodyIndex,
                          int *TileOffset);
static int OptLoadBodies(OptWorkflowStateStruct *Opt);
static void OptWriteCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt);
static int OptReadCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt);
static int OptValidateMesh(OptWorkflowStateStruct *Opt);
static int OptBodyForTile(const OptWorkflowStateStruct *Opt, int Tile);
static void OptTagTile(IritPrsrObjectStruct *Obj,
                       int Tile,
                       int *GeomSub);
static int OptAppendMeshOwner(OptWorkflowStateStruct *Opt,
                              int Tile,
                              int SubTile);
static int OptCollectMesh(OptWorkflowStateStruct *Opt,
                          IritPrsrObjectStruct *Obj,
                          IritPrsrObjectStruct *Out,
                          int SaveOwners);
static IritPrsrObjectStruct *PostProcessTileBodies(
                          IritPrsrObjectStruct *Tile,
                          IritUserMicroPostProcessTileCBStruct *CB);

/* Simulation data. */
static int ReadINP(OptWorkflowStateStruct *Opt, const char *INPPath);
static int ReadSimulationData(const char *FileName,
                              OptScalarFieldStruct *Field);
static CagdRType *ReadTileData(const OptWorkflowStateStruct *Opt,
                               const char *FileName,
                               int ExpectedCount);

/* Workflow state and containers. */
static void InitializeWorkflowState(OptWorkflowStateStruct *Opt);
static void ResolveLatticeLayout(OptWorkflowStateStruct *Opt);
static void InitOrResetResultsArray(OptWorkflowStateStruct *Opt);
static int InitTileStructure(OptWorkflowStateStruct *Opt);
static void ClearTileAnalysisData(OptWorkflowStateStruct *Opt);
static void DestroyTileStructure(OptWorkflowStateStruct *Opt);
static void FreeWorkflowState(OptWorkflowStateStruct *Opt);
static void NodeArrayReset(OptNodeArrayStruct *Arr);
static void ElementArrayReset(OptElementArrayStruct *Arr);
static int NodeArrayPush(OptNodeArrayStruct *Arr,
                         const OptNodeStruct *Node);
static int ElementArrayPush(OptElementArrayStruct *Arr,
                            const OptElementStruct *Elem);
static void ScalarFieldReset(OptScalarFieldStruct *Field);
static int ScalarFieldEnsureCapacity(OptScalarFieldStruct *Field,
                                     int RequiredId);
static void ScalarFieldSet(OptScalarFieldStruct *Field,
                           int Id,
                           CagdRType Value);
static CagdBType ScalarFieldTryGet(const OptScalarFieldStruct *Field,
                                   int Id,
                                   CagdRType *Value);
static int SubTileAddElementPtr(OptSubTileStruct *SubTile,
                                OptElementStruct *Elem);

/* Unified input parsing. */
static void OptSetDefaultRunInputs(OptRunInputStruct *Run);
static void OptSetDefaultConfigInputs(OptConfigInputStruct *Cfg);
static void OptSetDefaultConfig(OptConfigStruct *Config);
static void OptFreeConfig(OptConfigStruct *Config);
static int OptResolveDefaultPaths(OptConfigStruct *Config);
static int OptCopyConfigToWorkflow(const OptConfigStruct *Config,
                                   OptWorkflowStateStruct *Opt);
static int LoadInputs(const char *FileName, OptConfigStruct *Config);
static int ParseKeyValueFile(const char *FileName,
                             OptInputHandlerFuncType Handler,
                             void *Ctx);
static int InputHandler(const char *Key, const char *Value, void *Ctx);
static CagdBType OptParseBool(const char *Value);
static int ValidateInputs(const OptRunInputStruct *Run,
                          const OptConfigInputStruct *Cfg);

/* File and process helpers. */
static int EnsureDirectoryRecursive(const char *DirName);
static int CopyDirectoryRecursive(OptWorkflowStateStruct *Opt,
                                  const char *SourceDir,
                                  const char *TargetDir);
static int RunCommandLineBlocking(const char *CommandLine,
                                  const char *WorkingDir);
static int RunCommandLineBlockingToFile(const char *CommandLine,
                                        const char *WorkingDir,
                                        const char *LogFile);
static int OptWriteLine(const char *FileName,
                        const char *Content,
                        int Append);

/* General helpers. */
static void OptFatal(const char *Msg);
static void OptFatalF(const char *Fmt, const char *Arg);
static char *OptTrim(char *Str);
static int OptIsBlank(const char *Str);
static void OptStripUtf8Bom(char *Line);
static int OptParseInt(const char *Text, int *Value);
static int OptParseReal(const char *Text, CagdRType *Value);
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
static void *OptAllocZero(size_t Count, size_t Size);
static void *OptReallocCopy(void *OldData,
                            size_t OldBytes,
                            size_t NewBytes);

/* Allows fatal errors raised deep in the workflow to record their reason. */
IRIT_STATIC_DATA OptWorkflowStateStruct *GlblActiveWorkflow = NULL;

/***************************************************************************
*                           Program and Workflow                            *
***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
*   Program entry point.                                                   *
*                                                                          *
* PARAMETERS:                                                              *
*   Argc: Number of command-line arguments.                                *
*   Argv: Command-line argument vector.                                    *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: Process exit status.                                              *
*                                                                          *
* KEYWORDS:                                                               *
*   main, lattice optimization                                             *
***************************************************************************/
int main(int Argc, char **Argv)
{
    return OptimizationMain(Argc, Argv);
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Loads the unified configuration file and executes the optimization     *
* workflow.                                                                *
*                                                                          *
* PARAMETERS:                                                              *
*   Argc: Number of command-line arguments.                                *
*   Argv: Command-line argument vector.                                    *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: Process exit status.                                              *
***************************************************************************/
static int OptimizationMain(int Argc, char **Argv)
{
    CagdBType GeometryOnly = FALSE;
    const char *ConfigFile = "opt_config.txt";
    OptConfigStruct Config;
    int i;

    OptSetDefaultConfig(&Config);

    for (i = 1; i < Argc; ++i) {
        if (strcmp(Argv[i], "--geometry-only") == 0)
            GeometryOnly = TRUE;
        else
            ConfigFile = Argv[i];
    }

    if (!LoadInputs(ConfigFile, &Config))
        OptFatalF("Failed to parse configuration file: %s", ConfigFile);

    if (!ValidateInputs(&Config.Run, &Config.Environment))
        OptFatal("Input validation failed");

    i = IritOptRun(&Config, GeometryOnly);
    OptFreeConfig(&Config);
    return i ? EXIT_SUCCESS : EXIT_FAILURE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Runs an optimization from one complete, caller-owned configuration.    *
*   This is the entry point intended for the future IRIT user-library      *
*   binding.                                                               *
***************************************************************************/
int IritOptRun(const OptConfigStruct *Config, CagdBType GeometryOnly)
{
    OptWorkflowStateStruct OptData;
    OptWorkflowStateStruct *Opt = &OptData;
    int Ok;

    if (Config == NULL ||
        !ValidateInputs(&Config -> Run, &Config -> Environment))
        return FALSE;
    IRIT_ZAP_MEM(Opt, sizeof(*Opt));
    if (!OptCopyConfigToWorkflow(Config, Opt))
        return FALSE;

    if (!GeometryOnly)
        GlblActiveWorkflow = Opt;
    InitializeWorkflowState(Opt);
    if (GeometryOnly)
        Ok = InitLegacyDataFile(Opt) && RunIRITIntegrated(Opt, 0);
    else {
        Ok = RunWorkflow(Opt);
        GlblActiveWorkflow = NULL;
    }
    FreeWorkflowState(Opt);
    return Ok;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Executes the full optimization workflow until convergence or until the *
* maximal number of iterations is reached.                                 *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the workflow completed, FALSE on failure.                 *
***************************************************************************/

static int RunWorkflow(OptWorkflowStateStruct *Opt)
{
    int Iter, Status;
    char TargetBaseFolder[OPT_MAX_PATH];

    OptPathJoin(TargetBaseFolder,
                sizeof(TargetBaseFolder),
                Opt -> Config.BaseFolder,
                Opt -> FolderName);
    EnsureDirectoryRecursive(TargetBaseFolder);

    if (!InitLegacyDataFile(Opt))
        OptFatalF("Failed to initialize legacy data file: %s",
                  Opt -> Config.DataNewFile);

    if (!LatticeByVolumeFraction(Opt)) {
        OptWriteFinishReason(Opt, "Initial volume-fraction setup failed", 0);
        return FALSE;
    }

    for (Iter = 1; Iter <= Opt -> Run.MaxIterations; ++Iter) {
        Opt -> CurrentIteration = Iter;
        OptFormatNowHMS(Opt -> TimeIritStart);
        if (!RunIRITIntegrated(Opt, 0))
            OptFatal("IRIT generation step failed");

        OptFormatNowHMS(Opt -> TimeIritFinish);
        if (!RunSimulation(Opt))
            OptFatal("Simulation execution step failed");

        OptFormatNowHMS(Opt -> TimeSimulationFinish);
        Status = ProcessIterationData(Opt, Iter);
        if (Status < 0) {
            OptWriteFinishReason(Opt, "Iteration data processing failed", Iter);
            return FALSE;
        }
        if (Status > 0) {
            fprintf(stderr, "Optimization finished at iteration %d\n", Iter);
            OptWriteFinishReason(Opt, "Converged", Iter);
            return TRUE;
        }
    }

    OptWriteFinishReason(Opt, "Maximum iteration count reached",
                         Opt -> Run.MaxIterations);
    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Initializes the legacy thickness CSV file that is consumed by the      *
* geometry generation stage.                                               *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the file was created successfully, FALSE otherwise.       *
***************************************************************************/

static int InitLegacyDataFile(OptWorkflowStateStruct *Opt)
{
    FILE *F;
    int Times;
    int I;

    F = fopen(Opt -> Config.DataNewFile, "w");
    if (F == NULL)
        return 0;

    OptWriteCsvHeader(F, Opt);

    Times = Opt -> TotalTiles * Opt -> LMin;
    for (I = 0; I < Times; ++I)
        fprintf(F, "%d,%.17g,0\n", I + 1, Opt -> Run.InitialThickness);

    fclose(F);
    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Appends the single final status record to the run summary.             *
***************************************************************************/
static int OptWriteFinishReason(OptWorkflowStateStruct *Opt,
                                const char *Reason,
                                int Iter)
{
    char RunFolder[OPT_MAX_PATH];
    char SummaryPath[OPT_MAX_PATH];
    char Line[OPT_MAX_LINE];

    if (Opt == NULL || Reason == NULL || Opt -> FinishReasonWritten ||
        OptIsBlank(Opt -> FolderName))
        return FALSE;

    OptPathJoin(RunFolder, sizeof(RunFolder),
                Opt -> Config.BaseFolder, Opt -> FolderName);
    if (!EnsureDirectoryRecursive(RunFolder))
        return FALSE;
    OptPathJoin(SummaryPath, sizeof(SummaryPath),
                RunFolder, Opt -> Config.SummaryFile);
    snprintf(Line, sizeof(Line),
             "FinishReason,%s,Iteration,%d", Reason, Iter);
    if (!OptWriteLine(SummaryPath, Line, 1))
        return FALSE;
    Opt -> FinishReasonWritten = TRUE;
    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Executes the integrated geometry generation stage and verifies that an *
* INP output file was created.                                             *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt:  Active optimization workflow state.                              *
*   Mode: Geometry-generation mode.                                        *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the generation stage succeeded, FALSE otherwise.          *
***************************************************************************/

static int RunIRITIntegrated(OptWorkflowStateStruct *Opt, int Mode)
{
    if (!LoadLegacyGeometryData(Opt, Opt -> Config.DataNewFile)) {
        fprintf(stderr,
                "Failed to load legacy geometry CSV: %s\n",
                Opt -> Config.DataNewFile);
        return FALSE;
    }

    if (!GenerateMicroStructures(Opt, Mode)) {
        UnloadLegacyGeometryData(Opt);
        return FALSE;
    }

    UnloadLegacyGeometryData(Opt);

    if (Mode != 0)
        return TRUE;

    if (!OptFileExists(Opt -> Config.InpFile)) {
        fprintf(stderr,
                "Expected INP output was not created: %s\n",
                Opt -> Config.InpFile);
        return FALSE;
    }

    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Dispatches the selected simulation backend.                            *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the selected solver completed successfully, FALSE         *
*        otherwise.                                                        *
***************************************************************************/
static int RunSimulation(OptWorkflowStateStruct *Opt)
{
    switch (Opt -> Run.SimulationType) {
        case 1:
            return RunAnsysBatch(Opt);
        case 2:
            return RunCalculixBatch(Opt);
        default:
            fprintf(stderr,
                    "Unsupported simulation_type: %d\n",
                    Opt -> Run.SimulationType);
            return FALSE;
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Runs the ANSYS Workbench batch command for the current iteration.      *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
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
        Opt -> Config.AnsysExe,
        Opt -> Config.AnsysJournal,
        Opt -> Config.AnsysProject);

    return RunCommandLineBlocking(CmdLine, Opt -> Config.WorkDir);
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Runs a child process synchronously and redirects stdout and stderr to   *
* a fresh log file.                                                        *
*                                                                          *
* PARAMETERS:                                                              *
*   CommandLine: Command line to execute.                                  *
*   WorkingDir: Working directory of the child process.                    *
*   LogFile:    File that receives stdout and stderr.                      *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the child process exits successfully, FALSE otherwise.    *
***************************************************************************/
static int RunCalculixBatch(OptWorkflowStateStruct *Opt)
{
    enum {
        CCX_SEC_NONE,
        CCX_SEC_DISP,
        CCX_SEC_ENERGY,
        CCX_SEC_VOLUME
    } Section = CCX_SEC_NONE;
    const CagdRType TotalLoadY = -1000.0; /* Total vertical force, N. */
    const CagdRType PaToMPa = 1.0e-6;
    const CagdRType MToMm = 1.0e3;
    const CagdRType M3ToMm3 = 1.0e9;
    const char *EnvExe = NULL;
    char CcxExe[OPT_MAX_PATH] = { 0 };
    char TempInp[OPT_MAX_PATH] = { 0 };
    char JobName[OPT_MAX_NAME] = { 0 };
    char DatFile[OPT_MAX_PATH] = { 0 };
    char LogFile[OPT_MAX_PATH] = { 0 };
    char CmdLine[OPT_MAX_LINE * 2] = { 0 };
    char Line[OPT_MAX_LINE] = { 0 };
    FILE *Src = NULL;
    FILE *Dst = NULL;
    FILE *Dat = NULL;
    FILE *Sed = NULL;
    FILE *VolumeOut = NULL;
    FILE *DesignOut = NULL;
    CagdRType *ElementEnergy = NULL;
    CagdRType *ElementVolume = NULL;
    CagdRType *NodeDisp = NULL;
    int NodeMaxId = 0;
    int ElementMaxId = 0;
    int FixedCount = 0;
    int LoadedCount = 0;
    int DispCount = 0;
    int EnergyCount = 0;
    int VolumeCount = 0;
    int I = 0, N = 0;
    int InNodeSection = FALSE;
    int SkipGeneratedElset = FALSE;
    CagdRType ZTolerance = 0.0;
    CagdRType MinZ = 1.0e300;
    CagdRType MaxZ = -1.0e300;
    CagdRType MaxDeflectionM = 0.0;
    CagdRType TotalVolumeM3 = 0.0;
    int Ok = FALSE;

    /* Read the IRIT mesh first; this also gives us reliable node IDs and
     * coordinates for constructing the support and load node sets. */
    NodeArrayReset(&Opt -> Nodes);
    ElementArrayReset(&Opt -> Elements);
    if (!ReadINP(Opt, Opt -> Config.InpFile) ||
        Opt -> Nodes.Count == 0 || Opt -> Elements.Count == 0) {
        fprintf(stderr, "CalculiX: cannot read IRIT mesh %s\n",
                Opt -> Config.InpFile);
        goto Cleanup;
    }

    for (I = 0; I < Opt -> Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt -> Nodes.Items[I];
        MinZ = IRIT_MIN(MinZ, Node -> Z);
        MaxZ = IRIT_MAX(MaxZ, Node -> Z);
        NodeMaxId = IRIT_MAX(NodeMaxId, Node -> Id);
    }
    for (I = 0; I < Opt -> Elements.Count; ++I)
        ElementMaxId = IRIT_MAX(ElementMaxId, Opt -> Elements.Items[I].Id);

    ZTolerance = IRIT_MAX(1.0, fabs(MaxZ - MinZ)) * 1.0e-9;
    for (I = 0; I < Opt -> Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt -> Nodes.Items[I];
        if (fabs(Node -> Z - MinZ) <= ZTolerance)
            ++FixedCount;
        if (fabs(Node -> Z - MaxZ) <= ZTolerance &&
            Node -> Y > -0.015 && Node -> Y < 0.015)
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
    snprintf(TempInp, sizeof(TempInp), "%s.ccx.tmp", Opt -> Config.InpFile);
    Src = fopen(Opt -> Config.InpFile, "r");
    Dst = fopen(TempInp, "w");
    if (Src == NULL || Dst == NULL)
        goto Cleanup;

    while (fgets(Line, sizeof(Line), Src) != NULL) {
        char Copy[OPT_MAX_LINE] = { 0 };
        char *P = NULL;

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
            int Id = 0;
            double X = 0.0, Y = 0.0, Z = 0.0;
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
    fprintf(Dst,
            "** Support: minimum Z. Load: maximum Z and "
            "-0.015 < Y < 0.015.\n");
    fprintf(Dst, "*NSET, NSET=FIXED\n");
    N = 0;
    for (I = 0; I < Opt -> Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt -> Nodes.Items[I];
        if (fabs(Node -> Z - MinZ) <= ZTolerance) {
            fprintf(Dst, "%d%s", Node -> Id,
                    (++N % 16 == 0 || N == FixedCount) ? "\n" : ",");
        }
    }
    fprintf(Dst, "*NSET, NSET=LOADED\n");
    N = 0;
    for (I = 0; I < Opt -> Nodes.Count; ++I) {
        OptNodeStruct *Node = &Opt -> Nodes.Items[I];
        if (fabs(Node -> Z - MaxZ) <= ZTolerance &&
            Node -> Y > -0.015 && Node -> Y < 0.015) {
            fprintf(Dst, "%d%s", Node -> Id,
                    (++N % 16 == 0 || N == LoadedCount) ? "\n" : ",");
        }
    }

    /* Material, section, constraints, load, and all data needed by the
     * optimizer are kept in the same standalone INP file. */
    fprintf(Dst, "*MATERIAL, NAME=OPT_MATERIAL\n");
    fprintf(Dst, "*ELASTIC\n%.17g, %.17g\n",
            (double) Opt -> Config.YoungsModulus,
            (double) Opt -> Config.PoissonRatio);
    fprintf(Dst, "*DENSITY\n%.17g\n", (double) Opt -> Config.Density);
    fprintf(Dst, "*SOLID SECTION, ELSET=EL-ALL, MATERIAL=OPT_MATERIAL\n");
    /* SPOOLES is currently selected.  To switch back to Iterative Cholesky,
     * comment out the active fprintf line and uncomment the alternative. */
    fprintf(Dst, "*STEP\n*STATIC, SOLVER=SPOOLES\n");
    /* fprintf(Dst, "*STEP\n*STATIC, SOLVER=ITERATIVECHOLESKY\n"); */
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
    if (!MoveFileExA(TempInp, Opt -> Config.InpFile,
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_COPY_ALLOWED)) {
        fprintf(stderr, "CalculiX: cannot replace completed INP file.\n");
        goto Cleanup;
    }
#else
    if (rename(TempInp, Opt -> Config.InpFile) != 0)
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
        SYSTEM_INFO SysInfo = { 0 };
        char ThreadText[32] = { 0 }, Override[32] = { 0 };
        DWORD ThreadCount = 0, OverrideLen = 0;

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
             OptPathFileName(Opt -> Config.InpFile));
    {
        char *Dot = strrchr(JobName, '.');
        if (Dot != NULL)
            *Dot = '\0';
    }
    snprintf(CmdLine, sizeof(CmdLine), "\"%s\" -i \"%s\"", CcxExe, JobName);
    snprintf(LogFile, sizeof(LogFile), "%s%c%s_ccx.log",
             Opt -> Config.WorkDir, OPT_PATH_SEP, JobName);
    fprintf(stderr, "CalculiX solver log: %s\n", LogFile);
    if (!RunCommandLineBlockingToFile(CmdLine,
                                      Opt -> Config.WorkDir,
                                      LogFile)) {
        fprintf(stderr,
                "CalculiX solve failed for job %s; inspect %s\n",
                JobName, LogFile);
        goto Cleanup;
    }

    snprintf(DatFile, sizeof(DatFile), "%s%c%s.dat",
             Opt -> Config.WorkDir, OPT_PATH_SEP, JobName);
    Dat = fopen(DatFile, "r");
    if (Dat == NULL) {
        fprintf(stderr, "CalculiX result file was not created: %s\n", DatFile);
        goto Cleanup;
    }

    ElementEnergy = (CagdRType *) OptAllocZero((size_t) ElementMaxId + 1u,
                                                sizeof(CagdRType));
    ElementVolume = (CagdRType *) OptAllocZero((size_t) ElementMaxId + 1u,
                                                sizeof(CagdRType));
    NodeDisp = (CagdRType *) OptAllocZero((size_t) NodeMaxId + 1u,
                                          sizeof(CagdRType));
    if (ElementEnergy == NULL || ElementVolume == NULL || NodeDisp == NULL)
        goto Cleanup;

    /* The .dat headings make a much safer interface than fixed-column FRD.
     * ELSE and EVOL are whole-element values; no Gauss-point approximation
     * is needed for the element-average SED. */
    while (fgets(Line, sizeof(Line), Dat) != NULL) {
        char *P = OptTrim(Line);
        int Id = 0;
        double A = 0.0, B = 0.0, C = 0.0;

        if (strstr(P, "displacements (vx,vy,vz)") != NULL) {
            Section = CCX_SEC_DISP;
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
            /* ENER is diagnostic; ELSE/EVOL is exact. */
            Section = CCX_SEC_NONE;
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

    if (DispCount < Opt -> Nodes.Count ||
        EnergyCount != Opt -> Elements.Count ||
        VolumeCount != Opt -> Elements.Count) {
        fprintf(stderr,
                "CalculiX output count mismatch: disp=%d/%d, energy=%d/%d, "
                "volume=%d/%d.\n",
                DispCount, Opt -> Nodes.Count,
                EnergyCount, Opt -> Elements.Count,
                VolumeCount, Opt -> Elements.Count);
        goto Cleanup;
    }

    Sed = fopen(Opt -> Config.SedFile, "w");
    VolumeOut = fopen(Opt -> Config.VolumeFile, "w");
    DesignOut = fopen(Opt -> Config.DesignParametersFile, "w");
    if (Sed == NULL || VolumeOut == NULL || DesignOut == NULL)
        goto Cleanup;

    fprintf(Sed, "Element Number\tENERGYPOTENTIAL/VOLUME (mJ/mm^3)\n");
    fprintf(VolumeOut, "Element Number\tVOLUME (mm^3)\n");
    for (I = 0; I < Opt -> Elements.Count; ++I) {
        int Id = Opt -> Elements.Items[I].Id;
        if (ElementVolume[Id] <= 0.0)
            goto Cleanup;
        fprintf(Sed, "%d\t%.17g\n", Id,
                (double) (ElementEnergy[Id] / ElementVolume[Id] * PaToMPa));
        fprintf(VolumeOut, "%d\t%.17g\n", Id,
                (double) (ElementVolume[Id] * M3ToMm3));
    }

    /* This deliberately matches ReadLastRowOfDesignParameters(): the first
     * numeric column is maximum deformation in mm and the second is mass kg. */
    fprintf(DesignOut, "Name,P1,P2\n");
    fprintf(DesignOut, "DP 0,%.17g,%.17g\n",
            (double) (MaxDeflectionM * MToMm),
            (double) (TotalVolumeM3 * Opt -> Config.Density));

    fclose(Sed);
    Sed = NULL;
    fclose(VolumeOut);
    VolumeOut = NULL;
    fclose(DesignOut);
    DesignOut = NULL;

    fprintf(stderr,
            "CalculiX completed: %d elements, max deflection %.6g mm, "
            "mass %.6g kg.\n",
            Opt -> Elements.Count, (double) (MaxDeflectionM * MToMm),
            (double) (TotalVolumeM3 * Opt -> Config.Density));
    Ok = TRUE;

Cleanup:
    if (Src != NULL)
        fclose(Src);
    if (Dst != NULL)
        fclose(Dst);
    if (Dat != NULL)
        fclose(Dat);
    if (Sed != NULL)
        fclose(Sed);
    if (VolumeOut != NULL)
        fclose(VolumeOut);
    if (DesignOut != NULL)
        fclose(DesignOut);
    if (ElementEnergy != NULL)
        IritFree(ElementEnergy);
    if (ElementVolume != NULL)
        IritFree(ElementVolume);
    if (NodeDisp != NULL)
        IritFree(NodeDisp);
    NodeArrayReset(&Opt -> Nodes);
    ElementArrayReset(&Opt -> Elements);
    return Ok;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Reads the current solver results, refreshes the persistent analysis     *
* data, applies one design update, matches neighbouring interfaces, and    *
* writes the thickness state for the next iteration.                       *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt:  Active optimization workflow state.                              *
*   Iter: Iteration index to process.                                      *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: 1 if converged, 0 to continue, and -1 on processing failure.      *
***************************************************************************/
static int ProcessIterationData(OptWorkflowStateStruct *Opt, int Iter)
{
    int TileCount = Opt -> TotalTiles;
    CagdRType *PrevThickness;
    char TargetBaseFolder[OPT_MAX_PATH];

    NodeArrayReset(&Opt -> Nodes);
    ElementArrayReset(&Opt -> Elements);
    ScalarFieldReset(&Opt -> SedData);
    ScalarFieldReset(&Opt -> VolumeData);

    OptPathJoin(TargetBaseFolder,
                sizeof(TargetBaseFolder),
                Opt -> Config.BaseFolder,
                Opt -> FolderName);

    if (!OrganizeResults(Opt,
                         Opt -> Config.WorkDir,
                         TargetBaseFolder,
                         Iter))
        return -1;

    if (!ReadINP(Opt, Opt -> Config.InpFile))
        OptFatalF("Failed to parse INP file: %s", Opt -> Config.InpFile);

    if (!ReadSimulationData(Opt -> Config.SedFile, &Opt -> SedData))
        OptFatalF("Failed to parse SED file: %s", Opt -> Config.SedFile);

    if (Opt -> Config.VolumeFile[0] == '\0')
        OptFatal("Volume file path is empty in the configuration file.");

    if (!ReadSimulationData(Opt -> Config.VolumeFile,
                            &Opt -> VolumeData))
        OptFatalF("Failed to parse volume file: %s",
                  Opt -> Config.VolumeFile);

    PrevThickness = ReadTileData(Opt,
                                 Opt -> Config.DataNewFile,
                                 TileCount * Opt -> LMin);
    if (PrevThickness == NULL)
        OptFatalF("Failed to parse tile data file: %s",
                  Opt -> Config.DataNewFile);

    /* Refresh analysis fields without changing the thickness design. */
    UpdateTileAnalysisData(Opt);

    switch (Opt -> Run.OptimizationAlgorithm) {
        case 0:
            Opt -> LastThicknessChange =
                DesignUpdateSEDRatio(Opt,
                                     TileCount * Opt -> LMin,
                                     PrevThickness,
                                     Iter);
            break;
        case 1:
            Opt -> LastThicknessChange =
                DesignUpdateMMA(Opt,
                                TileCount * Opt -> LMin,
                                PrevThickness,
                                Iter);
            break;
        case 2:
            Opt -> LastThicknessChange =
                DesignUpdateOC(Opt,
                               TileCount * Opt -> LMin,
                               PrevThickness);
            break;
        case 3:
            Opt -> LastThicknessChange =
                DesignUpdateAsaf(Opt,
                                 TileCount * Opt -> LMin,
                                 PrevThickness,
                                 Iter);
            break;
        default:
            fprintf(stderr,
                    "Unknown optimization_algorithm = %d\n",
                    Opt -> Run.OptimizationAlgorithm);
            IritFree(PrevThickness);
            return -1;
    }

    NeighbourMatching(Opt);

#ifdef _WIN32
    {
        char DataPrev1[OPT_MAX_PATH];
        char DataPrev2[OPT_MAX_PATH];

        OptMakeIndexedFileName(Opt -> Config.DataPrevFile,
                               1,
                               DataPrev1,
                               sizeof(DataPrev1));
        OptMakeIndexedFileName(Opt -> Config.DataPrevFile,
                               2,
                               DataPrev2,
                               sizeof(DataPrev2));

        if (OptFileExists(DataPrev1))
            CopyFileA(DataPrev1, DataPrev2, FALSE);

        CopyFileA(Opt -> Config.DataNewFile, DataPrev1, FALSE);
        CopyFileA(Opt -> Config.DataNewFile,
                  Opt -> Config.DataPrevFile,
                  FALSE);
        DeleteFileA(Opt -> Config.DataNewFile);
    }
#endif

    CreateLegacyDataFile(Opt);

    /* Evaluate convergence after every thickness-changing operation. */
    Opt -> LastThicknessChange =
        ComputeMaxThicknessChange(Opt,
                                  PrevThickness,
                                  TileCount * Opt -> LMin);

    fprintf(stderr,
            "Final max thickness change for iteration %d: %.8g "
            "(tol %.8g)\n",
            Iter,
            Opt -> LastThicknessChange,
            Opt -> Run.ThicknessChangeTolerance);

    if (Opt -> LastThicknessChange <= Opt -> Run.ThicknessChangeTolerance) {
        ++Opt -> ConsecutiveSmallChanges;
        fprintf(stderr,
                "Small-change convergence count: %d/%d\n",
                Opt -> ConsecutiveSmallChanges,
                Opt -> Run.ConvergenceIterations);
    }
    else {
        Opt -> ConsecutiveSmallChanges = 0;
    }

    IritFree(PrevThickness);

    if (Opt -> ConsecutiveSmallChanges >=
        Opt -> Run.ConvergenceIterations) {
        fprintf(stderr,
                "Optimization converged: max thickness change <= %.8g "
                "for %d consecutive iterations.\n",
                Opt -> Run.ThicknessChangeTolerance,
                Opt -> Run.ConvergenceIterations);
        return TRUE;
    }

    return FALSE;
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
        fprintf(stderr,
                "Warning: failed to copy results folder to %s\n",
                TargetFolder);
    }

    if (!ReadLastRowOfDesignParameters(Opt -> Config.DesignParametersFile,
                                       &Opt -> MaxDef,
                                       &Opt -> Mass)) {
        fprintf(stderr,
                "Failed to read design parameters CSV: %s\n",
                Opt -> Config.DesignParametersFile);
        return 0;
    }

    RESULT_AT(Opt, Iter - 1, 0) = (CagdRType)Iter;
    RESULT_AT(Opt, Iter - 1, 1) = Opt -> MaxDef;
    RESULT_AT(Opt, Iter - 1, 2) = Opt -> Mass;

    OptPathJoin(SummaryPath,
                sizeof(SummaryPath),
                TargetBaseFolder,
                Opt -> Config.SummaryFile);

    /* Write CSV header once, when the summary file is first created. */
    if (!OptFileExists(SummaryPath)) {
        OptWriteLine(SummaryPath,
            "Iteration,MaxDeflection,Mass,TimeIritStart,"
            "TimeIritFinish,TimeSimulationFinish",
            0);
    }

    snprintf(SummaryLine, sizeof(SummaryLine),
        "%d,%.4f,%.4f,%s,%s,%s",
        Iter,
        (double)Opt -> MaxDef,
        (double)Opt -> Mass,
        Opt -> TimeIritStart,
        Opt -> TimeIritFinish,
        Opt -> TimeSimulationFinish);
    OptWriteLine(SummaryPath, SummaryLine, 1);

    return 1;
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
    CagdRType LastDef = 0.0;
    CagdRType LastMass = 0.0;

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
            if (OptParseReal(Parts[1], &TmpDef) &&
                OptParseReal(Parts[2], &TmpMass)) {
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
*                      Optimization and Design Update                       *
***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
*   Refreshes the persistent TileStructure with data from the current      *
*   ANSYS solution.  Topology and thickness are preserved.                  *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static void UpdateTileAnalysisData(OptWorkflowStateStruct *Opt)
{
    int I, L, X;

    /*
     * TileStructure was created once in InitializeWorkflowState().
     * Clear only fields that belong to the previous ANSYS solution.
     * SubTile -> Thickness is deliberately preserved.
     */
    if (Opt -> TileStructure == NULL)
        OptFatal("TileStructure was not initialized before analysis update");
    ClearTileAnalysisData(Opt);

    if (Opt -> Elements.Count != Opt -> MeshOwnerCount * Opt -> MeshSample) {
        fprintf(stderr,
                "Mesh ownership mismatch: elements=%d, owners=%d, "
                "configured elements/owner=%d.\n",
                Opt -> Elements.Count,
                Opt -> MeshOwnerCount,
                Opt -> MeshSample);
        OptFatal("Mesh element count does not match exported ownership map");
    }
    for (I = 0; I < Opt -> Elements.Count; ++I) {
        OptMeshOwnerStruct *Owner = &Opt -> MeshOwners[I / Opt -> MeshSample];
        if (Opt -> Elements.Items[I].Id != I + 1)
            OptFatal("ANSYS INP element order/IDs differ from exported mesh");
        if (!SubTileAddElementPtr(
                    &Opt -> TileStructure[Owner -> Tile].
                        SubTiles[Owner -> SubTile],
                    &Opt -> Elements.Items[I]))
            OptFatal("Cannot assign element to body/tile/subtile");
    }

    /* ------------------------------------------------------------
     * Assign element SED, volume, strain energy, and compliance
     * contribution.
     * ------------------------------------------------------------ */
    for (I = 0; I < Opt -> Elements.Count; ++I) {
        CagdRType SedVal = 0.0;
        CagdRType VolVal = 0.0;

        if (ScalarFieldTryGet(&Opt -> SedData,
            Opt -> Elements.Items[I].Id,
            &SedVal)) {
            Opt -> Elements.Items[I].SED = SedVal;
        }
        else {
            Opt -> Elements.Items[I].SED = 0.0;
        }

        if (ScalarFieldTryGet(&Opt -> VolumeData,
            Opt -> Elements.Items[I].Id,
            &VolVal)) {
            Opt -> Elements.Items[I].Volume = VolVal;
        }
        else {
            Opt -> Elements.Items[I].Volume = 0.0;
        }

        /*
         * Element strain energy:
         *     SE_e = SED_e * V_e
         *
         * Element compliance contribution:
         *     C_e ~= 2 * SE_e
         */
        Opt -> Elements.Items[I].StrainEnergy =
            Opt -> Elements.Items[I].SED *
            Opt -> Elements.Items[I].Volume;

        Opt -> Elements.Items[I].ComplianceContribution =
            2.0 * Opt -> Elements.Items[I].StrainEnergy;
    }

    /* ------------------------------------------------------------
     * Map elements into tiles/subtiles and compute subtile values.
     * ------------------------------------------------------------ */
    for (I = 0; I < Opt -> TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt -> TileStructure[I];

        Tile -> Volume = 0.0;
        Tile -> StrainEnergy = 0.0;
        Tile -> ComplianceContribution = 0.0;
        Tile -> Gain = 0.0;

        for (L = 0; L < Opt -> L; ++L) {
            OptSubTileStruct *SubTile = &Tile -> SubTiles[L];
            CagdRType SedSum = 0.0;
            CagdRType MaxSed = 0.0;
            CagdRType VolSum = 0.0;
            CagdRType SESum = 0.0;
            CagdRType CSum = 0.0;

            SubTile -> Volume = 0.0;
            SubTile -> StrainEnergy = 0.0;
            SubTile -> ComplianceContribution = 0.0;
            SubTile -> Gain = 0.0;

            for (X = 0; X < SubTile -> ElementCount; ++X) {
                OptElementStruct *Elem;

                Elem = SubTile -> Elements[X];

                SedSum += Elem -> SED;

                if (Elem -> SED > MaxSed)
                    MaxSed = Elem -> SED;

                VolSum += Elem -> Volume;
                SESum += Elem -> StrainEnergy;
                CSum += Elem -> ComplianceContribution;

            }

            SubTile -> Volume = VolSum;
            SubTile -> StrainEnergy = SESum;
            SubTile -> ComplianceContribution = CSum;

            /*
             * Prefer volume-weighted SED:
             *     SED_s = SE_s / V_s
             *
             * If volume is missing, fall back to the old arithmetic
             * average.
             */
            if (SubTile -> Volume > 0.0) {
                SubTile -> SED =
                    SubTile -> StrainEnergy / SubTile -> Volume;
            }
            else {
                SubTile -> SED =
                    SedSum / (SubTile -> ElementCount > 0 ?
                        SubTile -> ElementCount : 1);
            }

            SubTile -> MaxSED = MaxSed;

            /*
             * TO-like local update indicator:
             *     Gain = C_s / V_s
             *
             * Since C_s ~= 2 * SE_s, this is proportional to SED.
             */
            if (SubTile -> Volume > 0.0) {
                SubTile -> Gain =
                    SubTile -> ComplianceContribution /
                    SubTile -> Volume;
            }
            else {
                SubTile -> Gain = 0.0;
            }

            Tile -> Volume += SubTile -> Volume;
            Tile -> StrainEnergy += SubTile -> StrainEnergy;
            Tile -> ComplianceContribution +=
                SubTile -> ComplianceContribution;
        }

        if (Tile -> Volume > 0.0) {
            Tile -> Gain =
                Tile -> ComplianceContribution / Tile -> Volume;
        }
        else {
            Tile -> Gain = 0.0;
        }
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Updates the modifiable subtile thicknesses with the direct SED-ratio  *
* rule used when optimization_algorithm is zero.                           *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt:   Active optimization workflow state.                             *
*   Count: Number of thickness entries in Prev.                            *
*   Prev:  Previous thickness values.                                      *
*   Iter:  Current iteration index.                                        *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static CagdRType DesignUpdateSEDRatio(OptWorkflowStateStruct *Opt,
    int Count,
    const CagdRType *Prev,
    int Iter)
{
    int IndCounter = 0;
    int SedrCounter = 0;
    int I, L;
    CagdRType Avg = 0.0;
    CagdRType Var = 0.0;
    CagdRType MaxThicknessChange = 0.0;
    CagdRType PowFactor = 0.35;
    CagdRType MaxChange = 0.1;
    CagdRType TargetSed;
    CagdRType *Sedr;


    if (Opt == NULL || Prev == NULL || Count <= 0)
        return 0.0;

    if (Iter == 1) {
        PowFactor = 0.4;
        MaxChange = 0.4;
    }

    /* Preserve the original V1 SED-ratio definition. */
    TargetSed = Opt -> Config.SedValue * Opt -> Config.SafetyFactorSed;
    if (TargetSed <= 0.0)
        OptFatal("Invalid target SED value in DesignUpdateSEDRatio");

    Sedr = (CagdRType *)OptAllocZero((size_t)Count, sizeof(CagdRType));
    if (Sedr == NULL)
        OptFatal("Failed to allocate SED scratch array");

    for (I = 0; I < Opt -> TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt -> TileStructure[I];

        for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
            CagdRType TOld;
            CagdRType TNew;
            CagdRType SedFactor;
            CagdRType Change;
            CagdRType AbsThicknessChange;

            if (IndCounter >= Count) {
                fprintf(stderr,
                    "Warning: DesignUpdateSEDRatio reached Count limit early. "
                    "Count = %d\n",
                    Count);
                goto ComputeSedStats;
            }

            TOld = Prev[IndCounter];

            SedFactor = Tile -> SubTiles[L].SED / TargetSed;

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
                Tile -> SubTiles[L].SED,
                TargetSed,
                SedFactor,
                RawChange);

            if (Change > 1.0 + MaxChange)
                Change = 1.0 + MaxChange;

            if (Change < 1.0 - MaxChange)
                Change = 1.0 - MaxChange;

            TNew = TOld * Change;

            if (TNew < Opt -> Config.MinimumThickness)
                TNew = Opt -> Config.MinimumThickness;
            if (TNew > OPT_MAX_TILE_THICKNESS)
                TNew = OPT_MAX_TILE_THICKNESS;

            Tile -> SubTiles[L].Thickness = TNew;

            AbsThicknessChange = fabs(TNew - TOld);
            if (AbsThicknessChange > MaxThicknessChange)
                MaxThicknessChange = AbsThicknessChange;

            if (TNew > Opt -> Config.MinimumThickness) {
                Sedr[SedrCounter++] = Tile -> SubTiles[L].SED;
            }

            ++IndCounter;
        }
    }


ComputeSedStats:

    if (IndCounter != Count) {
        fprintf(stderr,
            "Warning: DesignUpdateSEDRatio used %d thickness values, "
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

        Opt -> SedAvg = Avg;
        Opt -> SedR = sqrt(Var);
    }
    else {
        Opt -> SedAvg = 0.0;
        Opt -> SedR = 0.0;
    }

    IritFree(Sedr);

    return MaxThicknessChange;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Updates the modifiable subtile thicknesses with NLopt MMA using local  *
* compliance and volume sensitivities from the current finite-element      *
* solution. Move limits and thickness bounds constrain each design step.   *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt:   Active optimization workflow state.                             *
*   Count: Number of modifiable thickness entries.                         *
*   Prev:  Previous thickness values.                                      *
*   Iter:  Current optimization iteration number.                          *
*                                                                          *
* RETURN VALUE:                                                            *
*   CagdRType: Maximum absolute thickness change in this iteration.        *
***************************************************************************/
static CagdRType DesignUpdateMMA(OptWorkflowStateStruct *Opt,
                                 int Count,
                                 const CagdRType *Prev,
                                 int Iter)
{
    const CagdRType Eps = 1.0e-12;
    const CagdRType TMin = Opt -> Config.MinimumThickness;
    const CagdRType TMax = OPT_MAX_TILE_THICKNESS;
    const CagdRType MoveFrac = OptDesignMoveFraction(Iter);
    int I, L, Counter;
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

    if (X == NULL || LB == NULL || UB == NULL ||
        ObjGrad == NULL || VolGrad == NULL)
        OptFatal("Failed to allocate MMA scratch arrays");

    Counter = 0;
    for (I = 0; I < Opt -> TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt -> TileStructure[I];

        for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
            OptSubTileStruct *SubTile = &Tile -> SubTiles[L];
            CagdRType TOld = Iter <= 1 || Prev[Counter] <= TMin ?
                             Opt -> Run.InitialThickness : Prev[Counter];
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

            CurrentVolume += SubTile -> Volume;
            CurrentCompliance += SubTile -> ComplianceContribution;

            ObjGrad[Counter] = -SubTile -> ComplianceContribution / TSafe;
            VolGrad[Counter] = SubTile -> Volume / TSafe;

            ++Counter;
        }
    }


    /*
     * Element volumes come from ANSYS and can be expressed in a different
     * unit system than the IRIT geometry. Preserve their relative
     * distribution for the sensitivity, but normalize their sum to the
     * exact current lattice volume used by the volume-fraction target.
     */
    if (CurrentVolume > Eps && Opt -> LatticeVol > Eps) {
        VolumeUnitScale = Opt -> LatticeVol / CurrentVolume;
        for (I = 0; I < Count; ++I)
            VolGrad[I] *= VolumeUnitScale;
        CurrentVolume = Opt -> LatticeVol;
    }

    MMAData.Count = Count;
    MMAData.PrevX = (CagdRType *) Prev;
    MMAData.ObjGrad = ObjGrad;
    MMAData.VolGrad = VolGrad;
    MMAData.ObjAtPrev = CurrentCompliance;
    MMAData.TargetVolume = Opt -> Run.TargetVolumeFraction * Opt -> DefMapVol;

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
                "Warning: NLopt MMA failed with result %d. "
                "Keeping previous thicknesses.\n",
                (int) Result);
        for (I = 0; I < Count; ++I)
            X[I] = Prev[I];
    }

    Counter = 0;
    for (I = 0; I < Opt -> TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt -> TileStructure[I];

        for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
            CagdRType TNew = X[Counter];
            CagdRType Change = fabs(TNew - Prev[Counter]);

            if (TNew < TMin)
                TNew = TMin;
            if (TNew > TMax)
                TNew = TMax;

            Tile -> SubTiles[L].Thickness = TNew;

            if (Change > MaxThicknessChange)
                MaxThicknessChange = Change;

            if (TNew > TMin && Tile -> SubTiles[L].SED > 0.0) {
                SedAvg += Tile -> SubTiles[L].SED;
                SedSqAvg += Tile -> SubTiles[L].SED * Tile -> SubTiles[L].SED;
                ++SedCount;
            }

            ++Counter;
        }
    }


    if (SedCount > 0) {
        SedAvg /= SedCount;
        SedSqAvg /= SedCount;
        Opt -> SedAvg = SedAvg;
        Opt -> SedR = sqrt(fabs(SedSqAvg - SedAvg * SedAvg));
    }
    else {
        Opt -> SedAvg = 0.0;
        Opt -> SedR = 0.0;
    }

    fprintf(stderr,
            "DesignUpdateMMA: result=%d, compliance=%.8g, volume=%.8g, "
            "target=%.8g, unit-scale=%.8g, change=%.8g, SedAvg=%.8g, "
            "SedStd=%.8g\n",
            (int) Result,
            CurrentCompliance,
            CurrentVolume,
            MMAData.TargetVolume,
            VolumeUnitScale,
            MaxThicknessChange,
            Opt -> SedAvg,
            Opt -> SedR);

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
* Evaluates the linearized MMA compliance objective and its gradient.      *
***************************************************************************/
static double MMAObjective(unsigned n,
                           const double *x,
                           double *Grad,
                           void *Data)
{
    OptMMADataStruct *MMAData = (OptMMADataStruct *) Data;
    double Value;
    unsigned I;

    Value = MMAData -> ObjAtPrev;

    for (I = 0; I < n; ++I) {
        if (Grad != NULL)
            Grad[I] = MMAData -> ObjGrad[I];

        Value += MMAData -> ObjGrad[I] * (x[I] - MMAData -> PrevX[I]);
    }

    return Value;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Evaluates the linearized MMA volume constraint and its gradient.         *
***************************************************************************/
static double MMAVolumeConstraint(unsigned n,
                                  const double *x,
                                  double *Grad,
                                  void *Data)
{
    OptMMADataStruct *MMAData = (OptMMADataStruct *) Data;
    double Value;
    unsigned I;

    Value = -MMAData -> TargetVolume;

    for (I = 0; I < n; ++I) {
        if (Grad != NULL)
            Grad[I] = MMAData -> VolGrad[I];

        Value += MMAData -> VolGrad[I] * x[I];
    }

    return Value;
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
*                                                                           *
* RETURN VALUE:                                                             *
*   CagdRType: Maximum absolute thickness change in this iteration.         *
****************************************************************************/
static CagdRType DesignUpdateOC(OptWorkflowStateStruct *Opt,
                                int Count,
                                const CagdRType *Prev)
{
    const CagdRType Eta = 0.5;
    const CagdRType TMin = Opt -> Config.MinimumThickness;
    const CagdRType TMax = OPT_MAX_TILE_THICKNESS;
    const CagdRType Move = (TMax - TMin) / 5.0;
    const CagdRType Eps = 1.0e-12;
    const int MaxBisectIter = 60;

    int I, L, It, Counter;
    int SedCount = 0;
    CagdRType L1 = Eps;
    CagdRType L2 = 1.0e12;
    CagdRType Lambda = 1.0;
    CagdRType TargetSum;
    CagdRType SumOld = 0.0;
    CagdRType SumNew = 0.0;
    CagdRType MaxAbsChange = 0.0;
    CagdRType SedAvg = 0.0;
    CagdRType SedSqAvg = 0.0;
    CagdRType SedTarget = Opt -> Config.SedValue /
        Opt -> Config.SafetyFactorSed;

    if (SedTarget <= Eps)
        SedTarget = Eps;

    TargetSum = Count * Opt -> Run.InitialThickness;
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

        for (I = 0; I < Opt -> TotalTiles; ++I) {
            OptTileStruct *Tile = &Opt -> TileStructure[I];

            for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
                CagdRType TOld = Prev[Counter];
                CagdRType Sed = Tile -> SubTiles[L].SED;
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
    for (I = 0; I < Opt -> TotalTiles; ++I) {
        OptTileStruct *Tile = &Opt -> TileStructure[I];

        for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
            CagdRType TOld = Prev[Counter];
            CagdRType Sed = Tile -> SubTiles[L].SED;
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

            Tile -> SubTiles[L].Thickness = TNew;
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
        Opt -> SedAvg = SedAvg;
        Opt -> SedR = sqrt(fabs(SedSqAvg - SedAvg * SedAvg));
    }
    else {
        Opt -> SedAvg = 0.0;
        Opt -> SedR = 0.0;
    }

    fprintf(stderr,
            "DesignUpdate: Tsum old=%.8g, Tsum new=%.8g, target=%.8g, "
            "change=%.8g, SedAvg=%.8g, SedStd=%.8g\n",
            SumOld, SumNew, TargetSum, MaxAbsChange,
            Opt -> SedAvg, Opt -> SedR);

    return MaxAbsChange;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Updates thicknesses using the direct ASAF SED-ratio rule while         *
* preserving the physical body/W/V/U/subtile traversal order.              *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt:   Active optimization workflow state.                             *
*   Count: Number of modifiable thickness entries.                         *
*   Prev:  Previous thickness values.                                      *
*   Iter:  Current optimization iteration number.                          *
*                                                                          *
* RETURN VALUE:                                                            *
*   CagdRType: Maximum absolute thickness change in this iteration.        *
***************************************************************************/
static CagdRType DesignUpdateAsaf(OptWorkflowStateStruct *Opt,
                                  int Count,
                                  const CagdRType *Prev,
                                  int Iter)
{
    int Counter = 0;
    int B, U, V, W, L;
    CagdRType MaxThicknessChange = 0.0;
    CagdRType PowFactor = 0.35;
    CagdRType MaxChange = 0.3;
    CagdRType TargetSed;

    if (Opt == NULL || Prev == NULL || Count <= 0)
        return 0.0;

    if (Iter == 1) {
        PowFactor = 0.4;
        MaxChange = 0.3;
    }

    TargetSed = Opt -> Config.SedValue * Opt -> Config.SafetyFactorSed;
    if (TargetSed <= 0.0)
        OptFatal("Invalid target SED value in DesignUpdateAsaf");

    /*
     * Use the original physical traversal explicitly:
     *     body -> W -> V -> U -> SubTile.
     *
     * This keeps Prev[Counter], TileStructure and data_new.csv in the same
     * unambiguous order, including multi-body runs.
     */
    for (B = 0; B < Opt -> BodyCount; ++B) {
        OptBodyStruct *Body = &Opt -> Bodies[B];

        for (W = 0; W < Body -> Counts[2]; ++W) {
            for (V = 0; V < Body -> Counts[1]; ++V) {
                for (U = 0; U < Body -> Counts[0]; ++U) {
                    int TileIndex = Body -> TileOffset +
                        U + Body -> Counts[0] * (V + Body -> Counts[1] * W);
                    OptTileStruct *Tile = &Opt -> TileStructure[TileIndex];

                    for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
                        CagdRType TOld;
                        CagdRType TNew;
                        CagdRType SedFactor;
                        CagdRType Change;
                        CagdRType AbsThicknessChange;

                        if (Counter >= Count)
                            OptFatal(
                                "DesignUpdateAsaf thickness count mismatch");

                        TOld = Prev[Counter];
                        SedFactor = Tile -> SubTiles[L].SED / TargetSed;
                        if (SedFactor < 0.0 || !isfinite(SedFactor))
                            SedFactor = 0.0;

                        Change = pow(SedFactor, PowFactor);

                        if (Change > 1.0 + MaxChange)
                            Change = 1.0 + MaxChange;
                        if (Change < 1.0 - MaxChange)
                            Change = 1.0 - MaxChange;

                        TNew = TOld * Change;

                        if (TNew < Opt -> Config.MinimumThickness)
                            TNew = Opt -> Config.MinimumThickness;
                        if (TNew > OPT_MAX_TILE_THICKNESS)
                            TNew = OPT_MAX_TILE_THICKNESS;

                        Tile -> SubTiles[L].Thickness = TNew;

                        AbsThicknessChange = fabs(TNew - TOld);
                        if (AbsThicknessChange > MaxThicknessChange)
                            MaxThicknessChange = AbsThicknessChange;

                        fprintf(stderr,
                            "Body %d Tile[%d,%d,%d] GTile=%d Sub=%d "
                            "SED=%g Told=%g Tnew=%g Change=%g\n",
                            B + 1, U, V, W, TileIndex, L,
                            Tile -> SubTiles[L].SED, TOld, TNew, Change);

                        ++Counter;
                    }
                }
            }
        }
    }

    if (Counter != Count)
        OptFatal("DesignUpdateAsaf did not consume all thickness values");

    return MaxThicknessChange;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Computes the maximum absolute thickness change between the final       *
* current design and the previous iteration thickness vector. This is      *
* evaluated after all thickness-changing operations are complete.          *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt:   Active optimization workflow state.                             *
*   Prev:  Previous iteration thickness vector.                            *
*   Count: Number of entries expected in Prev.                             *
*                                                                          *
* RETURN VALUE:                                                            *
*   CagdRType: Maximum absolute thickness change.                          *
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

    for (I = 0; I < Opt -> TotalTiles; ++I) {
        const OptTileStruct *Tile = &Opt -> TileStructure[I];

        for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
            CagdRType Delta;

            if (Counter >= Count)
                OptFatal(
                    "Thickness count mismatch in ComputeMaxThicknessChange");

            Delta = fabs(Tile -> SubTiles[L].Thickness - Prev[Counter]);
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
*   Enforces equal thickness values on coincident tile interfaces.         *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static void NeighbourMatching(OptWorkflowStateStruct *Opt)
{
    static const CagdRType Local0[6][3] = {
        { 0.5, 0.5, 1.0 },  /* W+. */
        { 0.5, 0.5, 0.0 },  /* W-. */
        { 0.5, 1.0, 0.5 },  /* V+. */
        { 0.5, 0.0, 0.5 },  /* V-. */
        { 1.0, 0.5, 0.5 },  /* U+. */
        { 0.0, 0.5, 0.5 }   /* U-. */
    };

    static const CagdRType Local01[6][3] = {
        { 0.5, 0.5, 0.0 },  /* W-. */
        { 0.5, 0.5, 1.0 },  /* W+. */
        { 0.5, 0.0, 0.5 },  /* V-. */
        { 0.5, 1.0, 0.5 },  /* V+. */
        { 0.0, 0.5, 0.5 },  /* U-. */
        { 1.0, 0.5, 0.5 }   /* U+. */
    };

    static const CagdRType Local2[8][3] = {
        { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 },
        { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }
    };

    switch (Opt -> Run.LatticeType) {
        case 0:
            NeighbourMatchingB(Opt, Local0, 6, 0);
            break;
        case 1:
            NeighbourMatchingB(Opt, Local01, 6, 0);
            break;
        case 2:
            NeighbourMatchingB(Opt, Local2, 8, 8);
            break;
        default:
            NeighbourMatchingB(Opt, Local0, 6, 0);
            break;
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Orders neighbouring sample points lexicographically by X, Y, then Z.     *
***************************************************************************/

static int CompareNeighbourEntries(const void *P1, const void *P2)
{
    const OptNeighbourEntryStruct *A = (const OptNeighbourEntryStruct *)P1;
    const OptNeighbourEntryStruct *B = (const OptNeighbourEntryStruct *)P2;

    if (A -> X < B -> X)
        return -1;
    if (A -> X > B -> X)
        return 1;
    if (A -> Y < B -> Y)
        return -1;
    if (A -> Y > B -> Y)
        return 1;
    if (A -> Z < B -> Z)
        return -1;
    if (A -> Z > B -> Z)
        return 1;
    return 0;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Matches coincident arm-end thicknesses independently within each body.   *
***************************************************************************/
static void NeighbourMatchingB(OptWorkflowStateStruct *Opt,
                               const CagdRType Local[][3],
                               int LocalCount,
                               int LDefOffset)
{
    int B, U, V, W, L, I, J, G;
    int X = 0;
    int Total = Opt -> TotalTiles * LocalCount;
    OptNeighbourEntryStruct *Entries;

    Entries = (OptNeighbourEntryStruct *)
        OptAllocZero(Total, sizeof(OptNeighbourEntryStruct));
    if (Entries == NULL)
        OptFatal("Cannot allocate neighbour data");

    /* Match exact integer-grid interfaces separately within each body. */
    for (B = 0; B < Opt -> BodyCount; ++B) {
        int Start = X;
        OptBodyStruct *Body = &Opt -> Bodies[B];

        for (W = 0; W < Body -> Counts[2]; ++W) {
            for (V = 0; V < Body -> Counts[1]; ++V) {
                for (U = 0; U < Body -> Counts[0]; ++U) {
                    for (L = 0; L < LocalCount; ++L) {
                        int Tile = Body -> TileOffset + U +
                            Body -> Counts[0] *
                            (V + Body -> Counts[1] * W);

                        Entries[X].LinearId = Tile * LocalCount + L;
                        Entries[X].X = U + Local[L][0];
                        Entries[X].Y = V + Local[L][1];
                        Entries[X].Z = W + Local[L][2];
                        ++X;
                    }
                }
            }
        }

        qsort(Entries + Start,
              X - Start,
              sizeof(*Entries),
              CompareNeighbourEntries);

        for (I = Start; I < X; I = J) {
            CagdRType Avg = 0.0;

            for (J = I + 1;
                 J < X &&
                 CompareNeighbourEntries(&Entries[I], &Entries[J]) == 0;
                 ++J)
                ;

            for (G = I; G < J; ++G) {
                int Id = Entries[G].LinearId;
                OptSubTileStruct *SubTile =
                    &Opt -> TileStructure[Id / LocalCount].SubTiles[
                        Id % LocalCount + LDefOffset];

                Avg += SubTile -> Thickness;
            }

            Avg /= J - I;
            for (G = I; G < J; ++G) {
                int Id = Entries[G].LinearId;
                OptSubTileStruct *SubTile =
                    &Opt -> TileStructure[Id / LocalCount].SubTiles[
                        Id % LocalCount + LDefOffset];

                SubTile -> Thickness = Avg;
            }
        }
    }

    IritFree(Entries);
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Returns the relative move limit used by the MMA thickness update.        *
***************************************************************************/

static CagdRType OptDesignMoveFraction(int Iter)
{
    return Iter <= 1 ? 0.05 : 0.15;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Adjusts initial lattice thickness until the lattice volume fraction    *
* is within the requested target tolerance. It generates the lattice,      *
* computes the lattice/deformation-map volume ratio, and updates the       *
* initial thickness accordingly.                                           *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the target volume fraction was reached, FALSE otherwise.  *
***************************************************************************/

static int LatticeByVolumeFraction(OptWorkflowStateStruct *Opt)
{
    int Iter;
    const int MaxIter = 50;

    for (Iter = 0; Iter < MaxIter; Iter++) {
        if (!RunIRITIntegrated(Opt, 1))
            return 0;

        if (Opt -> DefMapVol <= 0.0) {
            fprintf(stderr, "ERROR: Invalid DefMap volume: %.16g\n",
                Opt -> DefMapVol);
            return 0;
        }
        {
            CagdRType VolFraction;
            CagdRType Target;
            CagdRType Tol;
            CagdRType LowerLimit;
            CagdRType UpperLimit;

            VolFraction = Opt -> LatticeVol / Opt -> DefMapVol;
            Target = Opt -> Run.TargetVolumeFraction;
            Tol = Opt -> Run.TargetVolumeFractionTol;

            LowerLimit = Target * (1.0 - Tol);
            UpperLimit = Target * (1.0 + Tol);

            if (VolFraction > UpperLimit) {
                fprintf(stdout,
                    "Finished iteration %d, volume ratio is: %.8g. "
                    "Target range is [%.8g, %.8g]. Redoing with smaller "
                    "thickness.\n",
                    Iter + 1, VolFraction, LowerLimit, UpperLimit);

                Opt -> Run.InitialThickness *= 0.9;
                fprintf(stdout, "New initial thickness: %.8g\n",
                    Opt -> Run.InitialThickness);

                fflush(stdout);

                InitLegacyDataFile(Opt);
            }

            else if (VolFraction < LowerLimit) {
                fprintf(stdout,
                    "Finished iteration %d, volume ratio is: %.8g. "
                    "Target range is [%.8g, %.8g]. Redoing with larger "
                    "thickness.\n",
                    Iter + 1, VolFraction, LowerLimit, UpperLimit);

                Opt -> Run.InitialThickness *= 1.1;
                fprintf(stdout, "New initial thickness: %.8g\n",
                    Opt -> Run.InitialThickness);

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

static CagdRType ComputeObjVol(IritPrsrObjectStruct *PObj)
{
    CagdRType Vol = 0.0;

    if (PObj == NULL) {
        return 0.0;
    }

    if (IRIT_PRSR_IS_OLST_OBJ(PObj)) {
        int i;
        IritPrsrObjectStruct *PObjTmp;

        for (i = 0; (PObjTmp = IritPrsrListObjectGet(PObj, i)) != NULL; i++) {
            Vol += ComputeObjVol(PObjTmp);
        }
    }
    else if (IRIT_PRSR_IS_TRIVAR_OBJ(PObj)) {
        IritTrivTVStruct *TV;

        for (TV = PObj -> U.Trivars; TV != NULL; TV = TV -> Pnext) {
            Vol += fabs(IritTrivTVVolume(TV, TRUE));
        }
    }

    return Vol;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Writes the legacy thickness CSV file consumed by the geometry callback *
* code.                                                                    *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   void                                                                   *
***************************************************************************/

static void CreateLegacyDataFile(OptWorkflowStateStruct *Opt)
{
    FILE *F;
    int B, U, V, W, L;
    int CounterFile = 0;

    F = fopen(Opt -> Config.DataNewFile, "w");
    if (F == NULL)
        OptFatalF("Failed to create data file: %s", Opt -> Config.DataNewFile);

    OptWriteCsvHeader(F, Opt);

    /* Same ordering as DesignUpdateAsaf: body -> W -> V -> U -> SubTile. */
    for (B = 0; B < Opt -> BodyCount; ++B) {
        OptBodyStruct *Body = &Opt -> Bodies[B];

        for (W = 0; W < Body -> Counts[2]; ++W) {
            for (V = 0; V < Body -> Counts[1]; ++V) {
                for (U = 0; U < Body -> Counts[0]; ++U) {
                    int TileIndex = Body -> TileOffset +
                        U + Body -> Counts[0] * (V + Body -> Counts[1] * W);

                    for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L) {
                        fprintf(F,
                            "%d,%.17g,%.17g\n",
                            CounterFile + 1,
                            (double) Opt -> TileStructure[TileIndex].
                                SubTiles[L].Thickness,
                            (double) Opt -> TileStructure[TileIndex].
                                SubTiles[L].SED);
                        ++CounterFile;
                    }
                }
            }
        }
    }

    if (CounterFile != Opt -> TotalTiles * Opt -> LMin) {
        fclose(F);
        OptFatal("CreateLegacyDataFile thickness count mismatch");
    }

    fclose(F);
}

/***************************************************************************
*                            Geometry Generation                            *
***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
* Loads the legacy thickness CSV used by geometry-generation callbacks.    *
***************************************************************************/

static int LoadLegacyGeometryData(OptWorkflowStateStruct *Opt,
                                  const char *FileName)
{
    Opt -> Legacy.Data = ReadLegacyCsvToFloatArray(Opt,
                                                    FileName,
                                                    &Opt -> Legacy.Rows,
                                                    &Opt -> Legacy.Cols);
    return Opt -> Legacy.Data != NULL;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Releases the legacy CSV table and clears its dimensions.                 *
***************************************************************************/

static void UnloadLegacyGeometryData(OptWorkflowStateStruct *Opt)
{
    FreeLegacyCsv(Opt -> Legacy.Data, Opt -> Legacy.Rows);
    Opt -> Legacy.Data = NULL;
    Opt -> Legacy.Rows = 0;
    Opt -> Legacy.Cols = 0;
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
    int I;
    int N = Opt -> TotalTiles * Opt -> LMin;
    CagdRType *Values = ReadTileData(Opt, FileName, N);
    float **Data;
    if (!Values)
        return NULL;
    Data = (float **) OptAllocZero(N, sizeof(float *));
    if (Data == NULL) {
        IritFree(Values);
        return NULL;
    }
    for (I = 0; I < N; ++I) {
        Data[I] = (float *) OptAllocZero(3, sizeof(float));
        if (Data[I] == NULL) {
            FreeLegacyCsv(Data, I);
            IritFree(Values);
            return NULL;
        }
        Data[I][0] = (float) (I + 1);
        Data[I][1] = (float) Values[I];
    }
    IritFree(Values);
    Opt -> Legacy.LType = Opt -> Run.LatticeType;
    Opt -> Legacy.MeshSample = Opt -> Run.MeshSample;
    *OutRows = N;
    *OutCols = 3;
    return Data;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Frees a row-oriented legacy CSV table allocated by the parser.           *
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
* DESCRIPTION:                                                             *
*   Prepares the boundary parameters of one face of a synthesized tile.    *
*                                                                          *
* PARAMETERS:                                                              *
*   BPrm:       Boundary-parameter record to populate.                     *
*   OuterRadius: Requested outer radius.                                   *
*   Circ:       TRUE for circular cross sections, FALSE otherwise.         *
*   SkinThickness: Skin thickness at the four face corners; zeros disable. *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the parameter structure was filled successfully.          *
***************************************************************************/
static int PreProcessTile1FaceParam(
                    IritUserMicroTileBndryPrmStruct *BPrm,
                    CagdRType OuterRadius,
                    CagdBType Circ,
                    const CagdRType SkinThickness[4])
{
    IRIT_ZAP_MEM(BPrm, sizeof(IritUserMicroTileBndryPrmStruct));

    BPrm -> BndryShape = 0.1;
    BPrm -> CntrShape = 0.4;
    BPrm -> Circular = Circ;
    BPrm -> InnerRadius = 0.0;
    BPrm -> OuterRadius = OuterRadius;
    memcpy(BPrm->Bndry, SkinThickness, sizeof(BPrm->Bndry));
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
                                  IritPrsrObjectStruct *Tile,
                                  IritUserMicroPreProcessTileCBStruct *CBData)
{
    UserMicroLocalDataStruct *LclData =
        (UserMicroLocalDataStruct *) CBData -> CBFuncData;
    OptWorkflowStateStruct *Opt = LclData -> Opt;
    CagdRType *LclMinDmn = CBData -> TileLclDmnMin;
    CagdRType *LclMaxDmn = CBData -> TileLclDmnMax;
    IritGeomBBBboxStruct BBox;
    IritUserMicroTileBndryPrmStruct UMinPrms, UMaxPrms;
    IritUserMicroTileBndryPrmStruct VMinPrms, VMaxPrms;
    IritUserMicroTileBndryPrmStruct WMinPrms, WMaxPrms;
    int Index;
    int BaseIndex;
    int Arms = 6;
    CagdBType Circle = FALSE;

    assert(Tile == NULL);

    Index = CBData -> TileIdxs[0] +
        CBData -> TileIdxs[1] * LclData -> Body -> Counts[0] +
        CBData -> TileIdxs[2] * LclData -> Body -> Counts[0] *
            LclData -> Body -> Counts[1];

    fprintf(stderr,
            "Index:[%d], Tile[%d,%d,%d] from "
            "(%.3f, %.3f %.3f) to (%.3f, %.3f, %.3f)\r",
            Index,
            CBData -> TileIdxs[0],
            CBData -> TileIdxs[1],
            CBData -> TileIdxs[2],
            LclMinDmn[0], LclMinDmn[1], LclMinDmn[2],
            LclMaxDmn[0], LclMaxDmn[1], LclMaxDmn[2]);

    switch (Opt -> Legacy.LType) {
        case 0: /* Cross 1. */
        case 1: /* Cross 2. */
            Arms = 6;
            Circle = TRUE;
            break;
        case 2: /* Diagonal. */
            Arms = 8;
            Circle = FALSE;
            break;
        case 3: /* Combined. */
            Arms = 14;
            Circle = FALSE;
            break;
        default:
            Arms = 6;
            Circle = FALSE;
            break;
    }

    if (Index < 0 || Index >= LclData -> Body -> TileCount)
        OptFatal("Tile callback index outside body grid");

    BaseIndex = (LclData -> Body -> TileOffset + Index) * Arms;
    if (BaseIndex + Arms > Opt -> Legacy.Rows)
        OptFatal("Tile callback exceeds thickness data");

    if (Opt -> Config.RemoveThinTiles) {
        CagdRType Avg = 0.0;
        int i;

        for (i = 0; i < Arms; ++i)
            Avg += Opt -> Legacy.Data[BaseIndex + i][1];
        Avg /= Arms;

        if (Avg < OPT_MIN_TILE_THICKNESS) {
            fprintf(stderr,
                    "Warning: Removing thin tile at index %d with average "
                    "thickness %.6f\n",
                    Index,
                    (double) Avg);
            return NULL;
        }
    }

    switch (Opt -> Legacy.LType) {
        case 0: {
            static const CagdRType NoSkin[4] =
                { 0.0, 0.0, 0.0, 0.0 };
            CagdRType UniformSkin[4];
            int SkinIndex;

            for (SkinIndex = 0; SkinIndex < 4; ++SkinIndex)
                UniformSkin[SkinIndex] = Opt -> Run.SkinThickness;

            /* Currently applied to both U faces of every type-0 tile.
             * Restricting it to exterior macro-shape faces is a separate
             * topology-aware step. */
            if (!PreProcessTile1FaceParam(
                                      &WMaxPrms,
                                      Opt -> Legacy.Data[BaseIndex + 0][1],
                                      Circle, NoSkin) ||
                !PreProcessTile1FaceParam(
                                      &WMinPrms,
                                      Opt -> Legacy.Data[BaseIndex + 1][1],
                                      Circle, NoSkin) ||
                !PreProcessTile1FaceParam(
                                      &VMaxPrms,
                                      Opt -> Legacy.Data[BaseIndex + 2][1],
                                      Circle, NoSkin) ||
                !PreProcessTile1FaceParam(
                                      &VMinPrms,
                                      Opt -> Legacy.Data[BaseIndex + 3][1],
                                      Circle, NoSkin) ||
                !PreProcessTile1FaceParam(
                                      &UMaxPrms,
                                      Opt -> Legacy.Data[BaseIndex + 4][1],
                                      Circle, UniformSkin) ||
                !PreProcessTile1FaceParam(
                                      &UMinPrms,
                                      Opt -> Legacy.Data[BaseIndex + 5][1],
                                      Circle, UniformSkin)) {
                return NULL;
            }

            Tile = IritUserMicro3DCrossTile(&UMinPrms,
                                            &UMaxPrms,
                                            &VMinPrms,
                                            &VMaxPrms,
                                            &WMinPrms,
                                            &WMaxPrms,
                                            FALSE,
                                            NULL);

            IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
            if (BBox.Min[0] < -IRIT_UEPS ||
                BBox.Max[0] > 1.0 + IRIT_UEPS ||
                BBox.Min[1] < -IRIT_UEPS ||
                BBox.Max[1] > 1.0 + IRIT_UEPS ||
                BBox.Min[2] < -IRIT_UEPS ||
                BBox.Max[2] > 1.0 + IRIT_UEPS) {
                fprintf(stderr,
                        "Warning: Tile spans beyond the unit box.\n");
            }

            Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);
            break;
        }
        case 1: {
            static const int Map[6] = { 4, 5, 2, 3, 0, 1 };
            IrtRType OuterRadii[6];
            CagdRType Avg = 0.0;
            int i;

            for (i = 0; i < 6; ++i) {
                OuterRadii[i] =
                    Opt -> Legacy.Data[BaseIndex + Map[i]][1];
                Avg += OuterRadii[i];
            }
            Avg /= 6.0;

            Tile = IritUserMicro3DCrossTile2(FALSE,
                                             OuterRadii,
                                             NULL,
                                             2.5 * Avg,
                                             NULL);
            IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
            if (BBox.Min[0] < -IRIT_UEPS ||
                BBox.Max[0] > 1.0 + IRIT_UEPS ||
                BBox.Min[1] < -IRIT_UEPS ||
                BBox.Max[1] > 1.0 + IRIT_UEPS ||
                BBox.Min[2] < -IRIT_UEPS ||
                BBox.Max[2] > 1.0 + IRIT_UEPS) {
                fprintf(stderr,
                        "Warning: Tile spans beyond the unit box.\n");
            }

            Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);
            break;
        }
        case 2: {
            const CagdRType CntrSize = 0.1;
            const CagdRType CrnrSizes[8] = {
                Opt -> Legacy.Data[BaseIndex + 0][1],
                Opt -> Legacy.Data[BaseIndex + 1][1],
                Opt -> Legacy.Data[BaseIndex + 2][1],
                Opt -> Legacy.Data[BaseIndex + 3][1],
                Opt -> Legacy.Data[BaseIndex + 4][1],
                Opt -> Legacy.Data[BaseIndex + 5][1],
                Opt -> Legacy.Data[BaseIndex + 6][1],
                Opt -> Legacy.Data[BaseIndex + 7][1]
            };
            CagdRType CrnrVertScl[2] = { 0.4, 0.4 };

            Tile = IritUserMicroDiagTile1(CntrSize,
                                          1.0,
                                          CrnrSizes,
                                          CrnrVertScl,
                                          NULL,
                                          1,
                                          0.2,
                                          NULL,
                                          NULL);
            IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
            if (BBox.Min[0] < -IRIT_UEPS ||
                BBox.Max[0] > 1.0 + IRIT_UEPS ||
                BBox.Min[1] < -IRIT_UEPS ||
                BBox.Max[1] > 1.0 + IRIT_UEPS ||
                BBox.Min[2] < -IRIT_UEPS ||
                BBox.Max[2] > 1.0 + IRIT_UEPS) {
                fprintf(stderr,
                        "Warning: Tile spans beyond the unit box.\n");
            }

            Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);
            break;
        }
        case 3:
        default:
            /* Combined tile construction is not implemented here. */
            break;
    }

    return Tile;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Adapts FFD callback data to the existing per-tile geometry builder.      *
***************************************************************************/
static IritPrsrObjectStruct *PreProcessTileFFD(
    IritPrsrObjectStruct *Tile,
    IritTrivFFDInfoStruct *FFDInfo)
{
    UserMicroLocalDataStruct *LclData =
        (UserMicroLocalDataStruct *) FFDInfo -> _AuxInfo;
    IritUserMicroPreProcessTileCBStruct CBData;
    CagdRType Min[3], Max[3];
    int Counts[3], I;

    assert(Tile == NULL);
    IRIT_ZAP_MEM(&CBData, sizeof(CBData));
    CBData.CBFuncData = LclData;

    Counts[0] = LclData -> Body -> Counts[0];
    Counts[1] = LclData -> Body -> Counts[1];
    Counts[2] = LclData -> Body -> Counts[2];

    IritTrivTVDomain(LclData -> DefMap,
                    &Min[0], &Max[0],
                    &Min[1], &Max[1],
                    &Min[2], &Max[2]);

    for (I = 0; I < 3; ++I) {
        /* Recover zero-based grid indices for the legacy thickness lookup. */
        CBData.TileIdxs[I] = (int) (
            (FFDInfo -> _DmnMin[I] - Min[I]) *
            Counts[I] / (Max[I] - Min[I]) + 0.5);
        CBData.NumOfTiles[I] = Counts[I];
        CBData.TileLclDmnMin[I] =
            (FFDInfo -> _DmnMin[I] - Min[I]) / (Max[I] - Min[I]);
        CBData.TileLclDmnMax[I] =
            (FFDInfo -> _DmnMax[I] - Min[I]) / (Max[I] - Min[I]);
    }

    /* The existing builder scales the tile; FFD translates it afterwards. */
    IRIT_HMGN_MAT_COPY(CBData.Mat, FFDInfo -> _SclMat);
    return PreProcessTile(NULL, &CBData);
}

/****************************************************************************
* DESCRIPTION:                                                              *
*   Builds the heterogeneous microstructure for the current thickness       *
* state. It computes deformation-map and lattice volumes and, in export     *
* mode, writes the ITD and INP geometry files.                              *
*                                                                           *
* PARAMETERS:                                                               *
*   Opt:  Active optimization workflow state.                               *
*   Mode: Geometry-generation mode:                                         *
*         0 - Generate geometry and export ITD/INP.                         *
*         1 - Generate geometry and compute volumes only.                   *
*         2 - Compute deformation-map volume only.                          *
*                                                                           *
* RETURN VALUE:                                                             *
*   int: TRUE if the requested stage completed successfully, FALSE          *
*        otherwise.                                                         *
****************************************************************************/
static int GenerateMicroStructures(OptWorkflowStateStruct *Opt, int Mode)
{
    const int UseFFDTileAggregation = 0;
    /* Exact local deformation is required for full boundary skins on highly
     * curved maps; the approximate/global FFD path can invert skin TVs. */
    const int UseMicroStructFFD = 0;
    int B, D, Handler;
    int Ok = FALSE;
    CagdRType Steps[3][2];
    CagdRType ScaleVol = Opt -> Config.Scale * Opt -> Config.Scale *
                         Opt -> Config.Scale;
    IrtHmgnMatType Mat;
    IritMvarMVStruct *MV = NULL;
    IritPrsrObjectStruct *MS = NULL;
    IritPrsrObjectStruct *Scaled = NULL;
    IritPrsrObjectStruct *BodyMS = NULL;
    IritUserMicroParamStruct Param;
    IritUserMicroRegularParamStruct *Reg = &Param.U.RegularParam;
    IritTrivFFDInfoStruct FFDInfo;
    IritPrsrINPSaveDfltFileParamsStruct Save;
    UserMicroLocalDataStruct LclData;

    Opt -> DefMapVol = 0.0;
    for (B = 0; B < Opt -> BodyCount; ++B) {
        Opt -> DefMapVol +=
            fabs(IritTrivTVVolume(Opt -> Bodies[B].Map, TRUE)) * ScaleVol;
    }

    if (Mode == 2)
        return TRUE;

    MS = IritPrsrGenLISTObject(NULL);
    if (MS == NULL)
        return FALSE;

    if (Mode == 0)
        Opt -> MeshOwnerCount = 0;

    for (B = 0; B < Opt -> BodyCount; ++B) {
        OptBodyStruct *Body = &Opt -> Bodies[B];

        IRIT_ZAP_MEM(&LclData, sizeof(LclData));
        LclData.Opt = Opt;
        LclData.Body = Body;
        LclData.DefMap = Body -> Map;
        LclData.BndryThickness = Opt -> Run.SkinThickness;

        fprintf(stderr,
                "\nGenerating body %d/%d (%s), %d tiles\n",
                B + 1,
                Opt -> BodyCount,
                Body -> Name,
                Body -> TileCount);

        if (UseMicroStructFFD || !UseFFDTileAggregation) {
            IRIT_ZAP_MEM(&Param, sizeof(Param));
            MV = IritMvarCnvrtTVToMV(Body -> Map);
            if (MV == NULL)
                goto Cleanup;

            Param.TilingType = IRIT_USER_MICRO_TILE_REGULAR;
            Param.DeformMV = MV;
            Param.ApproxLowOrder = 2;
            Reg -> TilingStepMode = TRUE;
            Reg -> MaxPolyEdgeLen = 0.1;
            Reg -> GlobalPlacement = UseMicroStructFFD ? TRUE : FALSE;
            Reg -> FFDeformApprox = UseMicroStructFFD ? TRUE : FALSE;

            for (D = 0; D < 3; ++D) {
                Steps[D][0] = Body -> Counts[D];
                Steps[D][1] = 0;
                Reg -> TilingSteps[D].TilesPerIntervals = Steps[D];
                Reg -> TilingSteps[D].Len = 1;
            }

            Reg -> PreProcessCBFunc = PreProcessTile;
            Reg -> PostProcessCBFunc = PostProcessTileBodies;
            Reg -> CBFuncData = &LclData;
            BodyMS = IritUserMicroStructComposition(&Param);

            IritMvarMVFree(MV);
            MV = NULL;
        }
        else {
            IRIT_ZAP_MEM(&FFDInfo, sizeof(FFDInfo));
            FFDInfo.FFDCB = PreProcessTileFFD;
            FFDInfo._AuxInfo = &LclData;
            BodyMS = IritTrivFFDTileObjectInTV(NULL,
                                               Body -> Map,
                                               Body -> Counts[0],
                                               Body -> Counts[1],
                                               Body -> Counts[2],
                                               FALSE,
                                               0.0,
                                               0.0,
                                               NULL,
                                               &FFDInfo);
        }

        if (BodyMS == NULL ||
            !OptCollectMesh(Opt, BodyMS, MS, Mode == 0))
            goto Cleanup;

        IritPrsrFreeObject(BodyMS);
        BodyMS = NULL;
    }

    Opt -> LatticeVol = ComputeObjVol(MS) * ScaleVol;
    if (Mode == 0) {
        IritMiscMatGenMatUnifScale(Opt -> Config.Scale, Mat);
        Scaled = IritGeomTransformObject(MS, Mat);
        if (Scaled == NULL)
            goto Cleanup;

        IRIT_ZAP_MEM(&Save, sizeof(Save));
        for (D = 0; D < 3; ++D)
            Save.UVWSamples[D] = Opt -> Run.MeshSample;

        Save.Epsilon = IRIT_EPS;
        Save.BzrOnly = FALSE;
        strcpy(Save.RealFormat, "%.17g");
        Save.SaveDepth = 1;
        Save.GeometryType = IRIT_PRSR_FEM_ELEMENT_C3D8;
        Save.Messages = TRUE;

        if (!IritPrsrInpSaveFile(Scaled, Opt -> Config.InpFile, &Save))
            goto Cleanup;

        Handler = IritPrsrOpenDataFile(Opt -> Config.ItdOutputFile,
                                       FALSE,
                                       1);
        if (Handler < 0)
            goto Cleanup;

        IritPrsrPutObjectToHandler(Handler, MS);
        IritPrsrCloseStream(Handler, TRUE);
        if (!OptValidateMesh(Opt))
            goto Cleanup;
    }

    Ok = TRUE;

Cleanup:
    if (MV != NULL)
        IritMvarMVFree(MV);
    if (BodyMS != NULL)
        IritPrsrFreeObject(BodyMS);
    if (Scaled != NULL)
        IritPrsrFreeObject(Scaled);
    if (MS != NULL)
        IritPrsrFreeObject(MS);

    return Ok;
}

/***************************************************************************
*                        Multi-Body Geometry Support                        *
***************************************************************************/

/* Multi-body support used by the optimization driver.  This file is kept
 * separate because it is also convenient to regenerate it independently of
 * the main experimental driver. */

#define OPT_ATTR_TILE     "_opt_tile"
#define OPT_ATTR_SUBTILE  "_opt_subtile"

IRIT_STATIC_DATA const int
    OptCross1OwnerOrder[7] = { 6, 5, 4, 3, 2, 1, 0 };
IRIT_STATIC_DATA const int
    OptCross2OwnerOrder[7] = { 6, 4, 5, 2, 3, 0, 1 };

/***************************************************************************
* DESCRIPTION:                                                             *
* Parses a positive-integer lattice-count list.                           *
***************************************************************************/
static int OptParseCounts(const char *Text, int **Values, int *Count)
{
    const char *P = Text;
    int N = 0;
    int Capacity = 0;
    int *Result = NULL;

    if (Text == NULL || Values == NULL || Count == NULL)
        return FALSE;

    while (*P != '\0') {
        char *End;
        long Value;
        int *NewResult;

        while (isspace((unsigned char) *P) || *P == ',' || *P == ';' ||
               *P == 'x' || *P == 'X')
            ++P;
        if (*P == '\0')
            break;

        errno = 0;
        Value = strtol(P, &End, 10);
        if (End == P || errno == ERANGE || Value <= 0 || Value > INT_MAX)
            goto Error;
        P = End;
        if (*P != '\0' && !isspace((unsigned char) *P) && *P != ',' &&
            *P != ';' && *P != 'x' && *P != 'X')
            goto Error;

        if (N == Capacity) {
            int NewCapacity = Capacity == 0 ? 4 : Capacity * 2;
            NewResult = (int *) OptReallocCopy(Result,
                                                (size_t) Capacity * sizeof(int),
                                                (size_t) NewCapacity *
                                                    sizeof(int));
            if (NewResult == NULL)
                goto Error;
            Result = NewResult;
            Capacity = NewCapacity;
        }
        Result[N++] = (int) Value;
    }

    if (N == 0)
        goto Error;
    if (*Values != NULL)
        IritFree(*Values);
    *Values = Result;
    *Count = N;
    return TRUE;

Error:
    if (Result != NULL)
        IritFree(Result);
    return FALSE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Counts trivariates that represent macro-shape bodies.                   *
***************************************************************************/
static int OptCountBodyTVs(const IritPrsrObjectStruct *Obj)
{
    int I, Count = 0;
    const IritTrivTVStruct *TV;

    for (; Obj != NULL; Obj = Obj -> Pnext) {
        if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
            for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
                Count += OptCountBodyTVs(IritPrsrListObjectGet(Obj, I));
        }
        else if (IRIT_PRSR_IS_TRIVAR_OBJ(Obj)) {
            for (TV = Obj -> U.Trivars; TV != NULL; TV = TV -> Pnext)
                ++Count;
        }
    }
    return Count;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Copies macro-shape trivariates into the body table.                     *
***************************************************************************/
static int OptCopyBodyTVs(OptWorkflowStateStruct *Opt,
                          const IritPrsrObjectStruct *Obj,
                          int *BodyIndex,
                          int *TileOffset)
{
    int I;
    const IritTrivTVStruct *TV;

    for (; Obj != NULL; Obj = Obj -> Pnext) {
        if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
            for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
                if (!OptCopyBodyTVs(Opt, IritPrsrListObjectGet(Obj, I),
                                    BodyIndex, TileOffset))
                    return FALSE;
        }
        else if (IRIT_PRSR_IS_TRIVAR_OBJ(Obj)) {
            for (TV = Obj -> U.Trivars; TV != NULL; TV = TV -> Pnext) {
                OptBodyStruct *Body = &Opt -> Bodies[*BodyIndex];
                long long TileCount;

                Body -> Map = IritTrivTVCopy(TV);
                if (Body -> Map == NULL)
                    return FALSE;
                Body -> Map -> Pnext = NULL;
                snprintf(Body -> Name, sizeof(Body -> Name), "%s%s%d",
                         Obj -> ObjName != NULL ? Obj -> ObjName : "Body_",
                         Obj -> ObjName != NULL ? "_" : "", *BodyIndex + 1);
                for (I = 0; I < 3; ++I)
                    Body -> Counts[I] = Opt -> Run.BodyCounts[I][*BodyIndex];
                TileCount = (long long) Body -> Counts[0] * Body -> Counts[1] *
                            Body -> Counts[2];
                if (TileCount <= 0 || TileCount > INT_MAX - *TileOffset)
                    return FALSE;
                Body -> TileOffset = *TileOffset;
                Body -> TileCount = (int) TileCount;
                IritTrivTVDomain(Body -> Map,
                                 &Body -> Min[0], &Body -> Max[0],
                                 &Body -> Min[1], &Body -> Max[1],
                                 &Body -> Min[2], &Body -> Max[2]);
                *TileOffset += Body -> TileCount;
                ++*BodyIndex;
            }
        }
    }
    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Loads and initializes the multi-body macro-shape layout.                *
***************************************************************************/
static int OptLoadBodies(OptWorkflowStateStruct *Opt)
{
    const char *FileName = Opt -> Run.MacroShapePath;
    IritPrsrObjectStruct *Objects;
    int D;
    int B = 0;
    int TileOffset = 0;

    IritPrsrSetFlattenObjects(FALSE);
    Objects = IritPrsrGetDataFiles(&FileName, 1, TRUE, FALSE);
    if (Objects == NULL)
        return FALSE;
    Opt -> BodyCount = OptCountBodyTVs(Objects);
    if (Opt -> BodyCount <= 0) {
        IritPrsrFreeObjectList(Objects);
        return FALSE;
    }
    for (D = 0; D < 3; ++D) {
        if (Opt -> Run.BodyCountLengths[D] != Opt -> BodyCount) {
            fprintf(stderr,
                    "Lattice count list %d has %d values; macro shape has "
                    "%d bodies.\n",
                    D,
                    Opt -> Run.BodyCountLengths[D],
                    Opt -> BodyCount);
            IritPrsrFreeObjectList(Objects);
            return FALSE;
        }
    }
    Opt -> Bodies = (OptBodyStruct *) OptAllocZero((size_t) Opt -> BodyCount,
                                                 sizeof(OptBodyStruct));
    if (Opt -> Bodies == NULL ||
        !OptCopyBodyTVs(Opt, Objects, &B, &TileOffset)) {
        IritPrsrFreeObjectList(Objects);
        return FALSE;
    }
    IritPrsrFreeObjectList(Objects);
    Opt -> TotalTiles = TileOffset;
    return B == Opt -> BodyCount;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Writes the legacy lattice CSV header.                                   *
***************************************************************************/
static void OptWriteCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt)
{
    const OptBodyStruct *Body = &Opt -> Bodies[0];
    fprintf(F, "lattice_type,%d\n", Opt -> Run.LatticeType);
    fprintf(F, "lattice_counts,%d,%d,%d\n",
            Body -> Counts[0], Body -> Counts[1], Body -> Counts[2]);
    fprintf(F, "mesh_sample,%d\n", Opt -> Run.MeshSample);
    fprintf(F, "mesh_type,%d\n", OPT_LEGACY_MESH_TYPE);
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Reads and validates the legacy lattice CSV header.                      *
***************************************************************************/
static int OptReadCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt)
{
    char Line[OPT_MAX_LINE], Extra;
    int LatticeType, U, V, W, MeshSample, MeshType;
    if (fgets(Line, sizeof(Line), F) == NULL ||
        sscanf(Line, "lattice_type,%d %c", &LatticeType, &Extra) != 1 ||
        fgets(Line, sizeof(Line), F) == NULL ||
        sscanf(Line, "lattice_counts,%d,%d,%d %c", &U, &V, &W, &Extra) != 3 ||
        fgets(Line, sizeof(Line), F) == NULL ||
        sscanf(Line, "mesh_sample,%d %c", &MeshSample, &Extra) != 1 ||
        fgets(Line, sizeof(Line), F) == NULL ||
        sscanf(Line, "mesh_type,%d %c", &MeshType, &Extra) != 1)
        return FALSE;
    return LatticeType == Opt -> Run.LatticeType &&
           U == Opt -> Bodies[0].Counts[0] && V == Opt -> Bodies[0].Counts[1] &&
           W == Opt -> Bodies[0].Counts[2] &&
           MeshSample == Opt -> Run.MeshSample &&
           MeshType == OPT_LEGACY_MESH_TYPE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Returns the body index that owns a global tile index.                   *
***************************************************************************/
static int OptBodyForTile(const OptWorkflowStateStruct *Opt, int Tile)
{
    int B;
    for (B = 0; B < Opt -> BodyCount; ++B)
        if (Tile >= Opt -> Bodies[B].TileOffset &&
            Tile < Opt -> Bodies[B].TileOffset + Opt -> Bodies[B].TileCount)
            return B;
    return -1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Tags generated tile geometry with tile ownership attributes.            *
***************************************************************************/
static void OptTagTile(IritPrsrObjectStruct *Obj,
                       int Tile, int *GeomSub)
{
    int I;
    if (Obj == NULL)
        return;
    if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
        for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
            OptTagTile(IritPrsrListObjectGet(Obj, I), Tile, GeomSub);
    }
    else {
        IritMiscAttrSetObjectIntAttrib(Obj, OPT_ATTR_TILE, Tile);
        IritMiscAttrSetObjectIntAttrib(Obj, OPT_ATTR_SUBTILE, (*GeomSub)++);
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Appends one exported mesh-patch ownership record.                       *
***************************************************************************/
static int OptAppendMeshOwner(OptWorkflowStateStruct *Opt,
                              int Tile,
                              int SubTile)
{
    if (Opt -> MeshOwnerCount == Opt -> MeshOwnerCapacity) {
        int NewCapacity = Opt -> MeshOwnerCapacity == 0 ? 256 :
                          Opt -> MeshOwnerCapacity * 2;
        OptMeshOwnerStruct *NewOwners = (OptMeshOwnerStruct *)
            OptReallocCopy(Opt -> MeshOwners,
                           (size_t) Opt -> MeshOwnerCapacity *
                               sizeof(*NewOwners),
                           (size_t) NewCapacity * sizeof(*NewOwners));
        if (NewOwners == NULL)
            return FALSE;
        Opt -> MeshOwners = NewOwners;
        Opt -> MeshOwnerCapacity = NewCapacity;
    }
    Opt -> MeshOwners[Opt -> MeshOwnerCount].Tile = Tile;
    Opt -> MeshOwners[Opt -> MeshOwnerCount].SubTile = SubTile;
    ++Opt -> MeshOwnerCount;
    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Collects generated geometry and records patch ownership.                *
***************************************************************************/
static int OptCollectMesh(OptWorkflowStateStruct *Opt,
                          IritPrsrObjectStruct *Obj,
                          IritPrsrObjectStruct *Out,
                          int SaveOwners)
{
    int I, Tile, SubTile;
    IritTrivTVStruct *TV;
    if (Obj == NULL)
        return TRUE;
    if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
        for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
            if (!OptCollectMesh(Opt, IritPrsrListObjectGet(Obj, I), Out,
                                SaveOwners))
                return FALSE;
        return TRUE;
    }
    Tile = IritMiscAttrGetObjectIntAttrib(Obj, OPT_ATTR_TILE);
    SubTile = IritMiscAttrGetObjectIntAttrib(Obj, OPT_ATTR_SUBTILE);

    if (SaveOwners && IRIT_PRSR_IS_TRIVAR_OBJ(Obj) &&
        (Opt -> Run.LatticeType == 0 || Opt -> Run.LatticeType == 1)) {
        /* Cross1 export order is core,U-,U+,V-,V+,W-,W+, while its
         * design order is W+,W-,V+,V-,U+,U-,core. Cross2 exports
         * core,U-,U+,V-,V+,W-,W+. Map both to design indices. */
        const int *Order = Opt -> Run.LatticeType == 0 ?
            OptCross1OwnerOrder : OptCross2OwnerOrder;

        if (Tile < 0 || Tile >= Opt -> TotalTiles)
            return FALSE;
        for (I = 0, TV = Obj -> U.Trivars; TV != NULL; TV = TV -> Pnext, ++I) {
            if (I >= 7 || !OptAppendMeshOwner(Opt, Tile, Order[I]))
                return FALSE;
            IritPrsrListObjectAppend(
                Out, IritPrsrGenTRIVARObject(IritTrivTVCopy(TV)));
        }
        return I == 7;
    }

    if (SaveOwners && IRIT_PRSR_IS_TRIVAR_OBJ(Obj) &&
        Opt -> Run.LatticeType == 2) {
        /* DiagTile1 returns three eight-TV objects in this order:
         * center joints, diagonal rods, corner joints. Each internal list
         * was built with IRIT_LIST_PUSH, so vertex order is 7 down to 0.
         * Design variables occupy rod indices 8..15 (LDef == 8). */
        int GroupBase;

        if (Tile < 0 || Tile >= Opt -> TotalTiles ||
            SubTile < 0 || SubTile > 2)
            return FALSE;
        GroupBase = SubTile * 8;
        for (I = 0, TV = Obj -> U.Trivars; TV != NULL; TV = TV -> Pnext, ++I) {
            if (I >= 8 ||
                !OptAppendMeshOwner(Opt, Tile, GroupBase + (7 - I)))
                return FALSE;
            IritPrsrListObjectAppend(
                Out, IritPrsrGenTRIVARObject(IritTrivTVCopy(TV)));
        }
        return I == 8;
    }

    if (SaveOwners && (Tile < 0 || Tile >= Opt -> TotalTiles ||
                       SubTile < 0 || SubTile >= Opt -> L ||
                       !OptAppendMeshOwner(Opt, Tile, SubTile)))
        return FALSE;
    IritPrsrListObjectAppend(Out, IritPrsrCopyObject(NULL, Obj, TRUE));
    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Tags a generated tile after microstructure placement.                   *
***************************************************************************/
static IritPrsrObjectStruct *PostProcessTileBodies(
                            IritPrsrObjectStruct *Tile,
                            IritUserMicroPostProcessTileCBStruct *CB)
{
    UserMicroLocalDataStruct *Local =
        (UserMicroLocalDataStruct *) CB -> CBFuncData;
    OptBodyStruct *Body = Local -> Body;
    int SubTile = 0;
    int TileIndex = Body -> TileOffset + CB -> TileIdxs[0] +
        Body -> Counts[0] *
            (CB -> TileIdxs[1] + Body -> Counts[1] * CB -> TileIdxs[2]);
    OptTagTile(Tile, TileIndex, &SubTile);
    return Tile;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Validates all exported mesh ownership records.                          *
***************************************************************************/
static int OptValidateMesh(OptWorkflowStateStruct *Opt)
{
    int I;
    if (Opt -> MeshOwnerCount <= 0) {
        fprintf(stderr, "Generated mesh contains no owned geometry patches.\n");
        return FALSE;
    }
    for (I = 0; I < Opt -> MeshOwnerCount; ++I)
        if (Opt -> MeshOwners[I].Tile < 0 ||
            Opt -> MeshOwners[I].Tile >= Opt -> TotalTiles ||
            Opt -> MeshOwners[I].SubTile < 0 ||
            Opt -> MeshOwners[I].SubTile >= Opt -> L)
            return FALSE;
    return TRUE;
}



/***************************************************************************
*                     Simulation Data and Mesh Parsing                      *
***************************************************************************/

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
                if (!NodeArrayPush(&Opt -> Nodes, &Node)) {
                    fclose(F);
                    return 0;
                }
            }
        }
        else if (CurrentSection == SEC_ELEMENT) {
            OptElementStruct Elem;

            IRIT_ZAP_MEM(&Elem, sizeof(Elem));
            Elem.Id = atoi(P);
            if (Elem.Id > 0) {
                if (!ElementArrayPush(&Opt -> Elements, &Elem)) {
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
    int Id;
    int N = 0;
    double Thickness, Sed;
    if (!F)
        return NULL;
    if (!OptReadCsvHeader(F, Opt))
        goto Done;
    Prev = (CagdRType *) OptAllocZero(ExpectedCount, sizeof(CagdRType));
    if (!Prev)
        goto Done;
    while (fgets(Line, sizeof(Line), F)) {
        if (OptIsBlank(Line))
            continue;
        if (N >= ExpectedCount ||
            sscanf(Line, "%d,%lf,%lf %c", &Id, &Thickness, &Sed, &Extra) != 3 ||
            Id != N + 1 || !isfinite(Thickness) || Thickness <= 0.0 ||
            !isfinite(Sed)) goto Invalid;
        Prev[N++] = Thickness;
    }
    if (N != ExpectedCount)
        goto Invalid;
    goto Done;
Invalid:
    fprintf(stderr,
            "Invalid/missing/extra CSV design row at index %d in %s\n",
            N,
            FileName);
    IritFree(Prev);
    Prev = NULL;
Done:
    fclose(F);
    return Prev;
}

/***************************************************************************
*                       Workflow State and Containers                       *
***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
* Initializes per-run containers, body layout, tile topology, and results. *
***************************************************************************/

static void InitializeWorkflowState(OptWorkflowStateStruct *Opt)
{
    OptMakeRunFolderName(Opt -> FolderName);
    IRIT_ZAP_MEM(&Opt -> Nodes, sizeof(Opt -> Nodes));
    IRIT_ZAP_MEM(&Opt -> Elements, sizeof(Opt -> Elements));
    IRIT_ZAP_MEM(&Opt -> SedData, sizeof(Opt -> SedData));

    Opt -> LastThicknessChange = 1.0e30;
    Opt -> ConsecutiveSmallChanges = 0;

    ResolveLatticeLayout(Opt);
    if (!OptLoadBodies(Opt))
        OptFatal("Cannot initialize body layout");

    /*
     * Build the persistent tile/subtile hierarchy ONCE.
     * Thickness starts at InitialThickness and is changed later only by
     * design-update functions. ANSYS-derived fields are refreshed after
     * every solve without rebuilding this structure.
     */
    if (!InitTileStructure(Opt))
        OptFatal("Cannot initialize persistent TileStructure");

    InitOrResetResultsArray(Opt);

}

/***************************************************************************
* DESCRIPTION:                                                             *
* Maps the selected lattice type to its tile and design-variable layout.   *
***************************************************************************/

static void ResolveLatticeLayout(OptWorkflowStateStruct *Opt)
{
    switch (Opt -> Run.LatticeType) {
        case 0: /* Face (3D Cross 1): 6 face arms. */
        case 1: /* Face (3D Cross 2): 6 face arms. */
            Opt -> L = 7;
            Opt -> LMin = 6;
            Opt -> LDef = 0;
            break;
        case 2: /* Diagonal: 8 vertex arms. */
        case 3: /* Combined: 6 face + 8 vertex arms. */
            Opt -> L = 24;
            Opt -> LMin = 8;
            Opt -> LDef = 8;
            break;
        default:
            OptFatal("Unsupported lattice type");
            return;
    }

    Opt -> MeshSample = (Opt -> Run.MeshSample - 1) *
        (Opt -> Run.MeshSample - 1) *
        (Opt -> Run.MeshSample - 1);
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Allocates a cleared three-column result table for all iterations.        *
***************************************************************************/

static void InitOrResetResultsArray(OptWorkflowStateStruct *Opt)
{
    if (Opt -> Results != NULL) {
        IritFree(Opt -> Results);
        Opt -> Results = NULL;
    }

    Opt -> Results =
        (CagdRType *) OptAllocZero((size_t) Opt -> Run.MaxIterations * 3,
                                   sizeof(CagdRType));
    if (Opt -> Results == NULL)
        OptFatal("Failed to allocate results array");
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Creates the persistent tile/subtile hierarchy and initial thicknesses.   *
***************************************************************************/

static int InitTileStructure(OptWorkflowStateStruct *Opt)
{
    int TileCount;
    int T;

    /*
     * This is a true one-time initialization.  Re-entering it would hide
     * a workflow error, so an existing structure is simply left untouched.
     */
    if (Opt -> TileStructure != NULL)
        return 1;

    TileCount = Opt -> TotalTiles;
    if (TileCount <= 0 || Opt -> L <= 0)
        return 0;

    Opt -> TileStructure =
        (OptTileStruct *) OptAllocZero((size_t) TileCount,
                                       sizeof(OptTileStruct));
    if (Opt -> TileStructure == NULL)
        return 0;

    for (T = 0; T < TileCount; ++T) {
        OptTileStruct *Tile = &Opt -> TileStructure[T];

        Tile -> SubTileCount = Opt -> L;
        Tile -> SubTiles =
            (OptSubTileStruct *) OptAllocZero((size_t) Opt -> L,
                                              sizeof(OptSubTileStruct));
        if (Tile -> SubTiles == NULL) {
            DestroyTileStructure(Opt);
            return 0;
        }

        /* Tile identity/topology does not change during optimization. */
        snprintf(Tile -> Name,
                 sizeof(Tile -> Name),
                 "Body_%d_Tile_%d",
                 OptBodyForTile(Opt, T) + 1,
                 T);

        /*
         * Initial design state.  From this point onward Thickness is changed
         * only by design-update functions and interface matching such
         * as neighbour matching / volume correction).
         */
        {
            int L;
            for (L = Opt -> LDef; L < Opt -> LDef + Opt -> LMin; ++L)
                Tile -> SubTiles[L].Thickness = Opt -> Run.InitialThickness;
        }
    }

    return 1;
}

/***************************************************************************
* DESCRIPTION:                                                             *
*   Clears analysis values from the previous solver result while          *
* preserving the persistent tile/subtile topology, allocated element      *
* pointer buffers, and thickness design state.                             *
*                                                                          *
* PARAMETERS:                                                              *
*   Opt: Active optimization workflow state.                               *
*                                                                          *
* RETURN VALUE:                                                            *
*   void: No return value.                                                 *
***************************************************************************/
static void ClearTileAnalysisData(OptWorkflowStateStruct *Opt)
{
    int T, S;

    if (Opt == NULL || Opt -> TileStructure == NULL)
        return;

    for (T = 0; T < Opt -> TotalTiles; ++T) {
        OptTileStruct *Tile = &Opt -> TileStructure[T];

        Tile -> Volume = 0.0;
        Tile -> StrainEnergy = 0.0;
        Tile -> ComplianceContribution = 0.0;
        Tile -> Gain = 0.0;

        for (S = 0; S < Tile -> SubTileCount; ++S) {
            OptSubTileStruct *SubTile = &Tile -> SubTiles[S];

            /*
             * The old entries point into the old Opt -> Elements array, which
             * ProcessIterationData() has already freed.  Null the entries but
             * KEEP the allocated pointer buffer and its capacity for reuse.
             */
            if (SubTile -> Elements != NULL && SubTile -> ElementCapacity > 0) {
                IRIT_ZAP_MEM(SubTile -> Elements,
                             (size_t) SubTile -> ElementCapacity *
                                 sizeof(OptElementStruct *));
            }
            SubTile -> ElementCount = 0;

            /* Preserve SubTile -> Thickness: it is design state. */
            SubTile -> SED = 0.0;
            SubTile -> MaxSED = 0.0;
            SubTile -> Volume = 0.0;
            SubTile -> StrainEnergy = 0.0;
            SubTile -> ComplianceContribution = 0.0;
            SubTile -> Sensitivity = 0.0;
            SubTile -> Gain = 0.0;
        }
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Releases every tile, subtile, and finite-element ownership list.         *
***************************************************************************/

static void DestroyTileStructure(OptWorkflowStateStruct *Opt)
{
    int TileCount, T, S;

    if (Opt -> TileStructure == NULL)
        return;

    TileCount = Opt -> TotalTiles;
    for (T = 0; T < TileCount; ++T) {
        if (Opt -> TileStructure[T].SubTiles != NULL) {
            for (S = 0; S < Opt -> TileStructure[T].SubTileCount; ++S) {
                if (Opt -> TileStructure[T].SubTiles[S].Elements != NULL)
                    IritFree(Opt -> TileStructure[T].SubTiles[S].Elements);
            }
            IritFree(Opt -> TileStructure[T].SubTiles);
        }
    }

    IritFree(Opt -> TileStructure);
    Opt -> TileStructure = NULL;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Releases all allocations owned by one optimization workflow.             *
***************************************************************************/

static void FreeWorkflowState(OptWorkflowStateStruct *Opt)
{
    NodeArrayReset(&Opt -> Nodes);
    ElementArrayReset(&Opt -> Elements);
    ScalarFieldReset(&Opt -> SedData);
    DestroyTileStructure(Opt);
    {
        int B;
        for (B = 0; B < Opt -> BodyCount; ++B)
            IritTrivTVFree(Opt -> Bodies[B].Map);
        if (Opt -> Bodies)
            IritFree(Opt -> Bodies);
        if (Opt -> MeshOwners)
            IritFree(Opt -> MeshOwners);
        for (B = 0; B < 3; ++B)
            if (Opt -> Run.BodyCounts[B])
                IritFree(Opt -> Run.BodyCounts[B]);
        ScalarFieldReset(&Opt -> VolumeData);
    }

    if (Opt -> Results != NULL) {
        IritFree(Opt -> Results);
        Opt -> Results = NULL;
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Releases the node array and restores an empty reusable state.            *
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
* DESCRIPTION:                                                             *
* Releases the element array and restores an empty reusable state.         *
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
* DESCRIPTION:                                                             *
* Appends a node, growing the contiguous node array when required.         *
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
* DESCRIPTION:                                                             *
* Appends an element, growing the contiguous element array as needed.      *
***************************************************************************/

static int ElementArrayPush(OptElementArrayStruct *Arr,
                            const OptElementStruct *Elem)
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
* DESCRIPTION:                                                             *
* Releases scalar values and their parallel validity bitmap.               *
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
* DESCRIPTION:                                                             *
* Grows a sparse ID-indexed scalar field while preserving stored values.   *
***************************************************************************/

static int ScalarFieldEnsureCapacity(OptScalarFieldStruct *Field,
                                     int RequiredId)
{
    int NewCap;
    CagdRType *NewValues;
    CagdBType *NewHasValue;
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
* DESCRIPTION:                                                             *
* Stores a scalar by solver ID and marks that entry as available.          *
***************************************************************************/

static void ScalarFieldSet(OptScalarFieldStruct *Field, int Id, CagdRType Value)
{
    if (!ScalarFieldEnsureCapacity(Field, Id))
        OptFatal("Failed to grow scalar field");
    Field -> Values[Id] = Value;
    Field -> HasValue[Id] = TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Retrieves a scalar only when the requested solver ID has a value.        *
***************************************************************************/

static CagdBType ScalarFieldTryGet(const OptScalarFieldStruct *Field,
                                   int Id,
                                   CagdRType *Value)
{
    if (Id < 0 || Id >= Field -> Capacity || !Field -> HasValue[Id])
        return FALSE;

    *Value = Field -> Values[Id];
    return TRUE;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Associates one finite element with a subtile, growing its pointer list.  *
***************************************************************************/

static int SubTileAddElementPtr(OptSubTileStruct *SubTile,
                                OptElementStruct *Elem)
{
    if (SubTile -> ElementCount >= SubTile -> ElementCapacity) {
        int NewCap = SubTile -> ElementCapacity == 0 ?
            16 : SubTile -> ElementCapacity * 2;
        OptElementStruct **NewItems =
            (OptElementStruct **) OptReallocCopy(
                                        SubTile -> Elements,
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
*                           Unified Input Parsing                           *
***************************************************************************/

/****************************************************************************
* DESCRIPTION:                                                              *
* Initializes defaults for run-specific inputs before the configuration     *
* file is parsed.                                                           *
****************************************************************************/
static void OptSetDefaultRunInputs(OptRunInputStruct *Run)
{
    IRIT_ZAP_MEM(Run, sizeof(OptRunInputStruct));

    Run -> SimulationType = 1;
    Run -> LatticeType = -1;
    Run -> SkinThickness = 0.0;
    Run -> LatticeCountU = -1;
    Run -> LatticeCountV = -1;
    Run -> LatticeCountW = -1;
    Run -> InitialThickness = OPT_INVALID_REAL;
    Run -> MeshSample = -1;
    Run -> MaxIterations = -1;
    Run -> ConvergenceIterations = 3;
    Run -> ThicknessChangeTolerance =
        5.0e-3 * OPT_MAX_TILE_THICKNESS;
    Run -> TargetVolumeFraction = OPT_INVALID_REAL;
    Run -> TargetVolumeFractionTol = OPT_INVALID_REAL;
    Run -> OptimizationAlgorithm = 1;
}

/****************************************************************************
* DESCRIPTION:                                                              *
* Initializes defaults for machine, material, and optimization settings.    *
****************************************************************************/
static void OptSetDefaultConfigInputs(OptConfigInputStruct *Cfg)
{
    IRIT_ZAP_MEM(Cfg, sizeof(OptConfigInputStruct));

    Cfg -> Scale = 0.015;
    Cfg -> Density = 2700.0;
    Cfg -> YoungsModulus = 6.9e10;
    Cfg -> PoissonRatio = 0.33;
    Cfg -> SedValue = 0.15;
    Cfg -> MinimumThickness = 1.0e-5;
    Cfg -> SafetyFactorSed = 1.1;
    Cfg -> RemoveThinTiles = FALSE;
}

static void OptSetDefaultConfig(OptConfigStruct *Config)
{
    IRIT_ZAP_MEM(Config, sizeof(*Config));
    OptSetDefaultRunInputs(&Config -> Run);
    OptSetDefaultConfigInputs(&Config -> Environment);
    snprintf(Config -> Environment.AnsysExe,
             sizeof(Config -> Environment.AnsysExe),
             "%s",
             "C:\\Program Files\\ANSYS Inc\\v222\\Framework\\bin\\Win64\\RunWB2.exe");
}

static void OptFreeConfig(OptConfigStruct *Config)
{
    int D;

    for (D = 0; D < 3; ++D) {
        if (Config -> Run.BodyCounts[D] != NULL)
            IritFree(Config -> Run.BodyCounts[D]);
        Config -> Run.BodyCounts[D] = NULL;
        Config -> Run.BodyCountLengths[D] = 0;
    }
}

/* Expands the single configured root into the conventional run layout.
 * Explicit nonempty values remain valid overrides. */
static int OptResolveDefaultPaths(OptConfigStruct *Config)
{
    OptRunInputStruct *Run = &Config -> Run;
    OptConfigInputStruct *Cfg = &Config -> Environment;
    char AnsysDir[OPT_MAX_PATH];

    if (OptIsBlank(Cfg -> BaseFolder)) {
        fprintf(stderr, "Missing root_folder in configuration file.\n");
        return FALSE;
    }
    if (OptIsBlank(Cfg -> WorkDir))
        OptPathJoin(Cfg -> WorkDir, sizeof(Cfg -> WorkDir),
                    Cfg -> BaseFolder, "TEMP_SET");
    OptPathJoin(AnsysDir, sizeof(AnsysDir), Cfg -> BaseFolder, "ANSYS_SET");

#define OPT_DEFAULT_PATH(Field, Parent, Name) \
    do { if (OptIsBlank(Field)) \
        OptPathJoin(Field, sizeof(Field), Parent, Name); } while (0)
    if (OptIsBlank(Run -> MacroShapePath))
        snprintf(Run -> MacroShapePath,
                 sizeof(Run -> MacroShapePath),
                 "%s",
                 "DefMap.itd");
    {
        char MacroShapeFile[OPT_MAX_PATH];

        snprintf(MacroShapeFile, sizeof(MacroShapeFile), "%s",
                 Run -> MacroShapePath);
        OptPathJoin(Run -> MacroShapePath,
                    sizeof(Run -> MacroShapePath),
                    Cfg -> WorkDir,
                    MacroShapeFile);
    }
    OPT_DEFAULT_PATH(Cfg -> AnsysProject, AnsysDir, "WB2.wbpj");
    OPT_DEFAULT_PATH(Cfg -> AnsysJournal, AnsysDir, "RunWB.wbjn");
    /* The mesh exporter owns this file; it is not a user setting. */
    OptPathJoin(Cfg -> InpFile, sizeof(Cfg -> InpFile),
                Cfg -> WorkDir, "optimized_geo.inp");
    OPT_DEFAULT_PATH(Cfg -> ItdOutputFile, Cfg -> WorkDir, "Lever_new.itd");
    OPT_DEFAULT_PATH(Cfg -> SedFile, Cfg -> WorkDir, "SED.txt");
    OPT_DEFAULT_PATH(Cfg -> VolumeFile, Cfg -> WorkDir, "VOLUME.txt");
    OPT_DEFAULT_PATH(Cfg -> DataNewFile, Cfg -> WorkDir, "data_new.csv");
    OPT_DEFAULT_PATH(Cfg -> DataPrevFile, Cfg -> WorkDir, "data_prev.csv");
    OPT_DEFAULT_PATH(Cfg -> DesignParametersFile, Cfg -> WorkDir,
                     "DesignParameters.csv");
    if (OptIsBlank(Cfg -> SummaryFile))
        snprintf(Cfg -> SummaryFile, sizeof(Cfg -> SummaryFile),
                 "%s", "Summary.csv");
#undef OPT_DEFAULT_PATH
    return TRUE;
}

static int OptCopyConfigToWorkflow(const OptConfigStruct *Config,
                                   OptWorkflowStateStruct *Opt)
{
    int D;

    Opt -> Run = Config -> Run;
    Opt -> Config = Config -> Environment;
    for (D = 0; D < 3; ++D) {
        size_t Bytes;

        Opt -> Run.BodyCounts[D] = NULL;
        Bytes = (size_t) Config -> Run.BodyCountLengths[D] * sizeof(int);
        if (Bytes == 0)
            continue;
        Opt -> Run.BodyCounts[D] = (int *) OptAllocZero(1, Bytes);
        if (Opt -> Run.BodyCounts[D] == NULL) {
            FreeWorkflowState(Opt);
            return FALSE;
        }
        memcpy(Opt -> Run.BodyCounts[D], Config -> Run.BodyCounts[D], Bytes);
    }
    return TRUE;
}

/****************************************************************************
* DESCRIPTION:                                                              *
* Loads all run, optimization, and machine settings from one configuration   *
* file.                                                                     *
****************************************************************************/
static int LoadInputs(const char *FileName, OptConfigStruct *Config)
{
    if (!ParseKeyValueFile(FileName, InputHandler, Config))
        return FALSE;
    return OptResolveDefaultPaths(Config);
}

/****************************************************************************
* DESCRIPTION:                                                              *
* Reads a key-value configuration file and dispatches each parsed entry to   *
* the supplied handler.                                                     *
****************************************************************************/
static int ParseKeyValueFile(const char *FileName,
                             OptInputHandlerFuncType Handler,
                             void *Ctx)
{
    int LineNum = 0;
    FILE *F;
    char Line[OPT_MAX_LINE];

    F = fopen(FileName, "r");
    if (F == NULL) {
        fprintf(stderr, "Failed to open input file: %s\n", FileName);
        return FALSE;
    }

    while (fgets(Line, sizeof(Line), F) != NULL) {
        char *Eq, *Key, *Value;

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
            fprintf(stderr,
                    "Invalid key or value in %s at line %d: %s\n",
                    FileName,
                    LineNum,
                    Key);
            fclose(F);
            return FALSE;
        }
    }

    fclose(F);
    return TRUE;
}

/****************************************************************************
* DESCRIPTION:                                                              *
* Maps one key from the unified configuration file into the run or fixed     *
* configuration structure.                                                  *
****************************************************************************/
static int InputHandler(const char *Key, const char *Value, void *Ctx)
{
    OptConfigStruct *Config = (OptConfigStruct *) Ctx;
    OptRunInputStruct *Run = &Config -> Run;
    OptConfigInputStruct *Cfg = &Config -> Environment;

    if (strcmp(Key, "macro_shape_file") == 0) {
        snprintf(Run -> MacroShapePath,
                 sizeof(Run -> MacroShapePath),
                 "%s",
                 Value);
        return TRUE;
    }
    if (strcmp(Key, "simulation_type") == 0)
        return OptParseInt(Value, &Run -> SimulationType);
    if (strcmp(Key, "lattice_type") == 0)
        return OptParseInt(Value, &Run -> LatticeType);
    if (strcmp(Key, "skin_thickness") == 0)
        return OptParseReal(Value, &Run -> SkinThickness);
    if (strcmp(Key, "lattice_count_u") == 0) {
        if (!OptParseCounts(Value,
                            &Run -> BodyCounts[0],
                            &Run -> BodyCountLengths[0]))
            return FALSE;
        Run -> LatticeCountU = Run -> BodyCounts[0][0];
        return TRUE;
    }
    if (strcmp(Key, "lattice_count_v") == 0) {
        if (!OptParseCounts(Value,
                            &Run -> BodyCounts[1],
                            &Run -> BodyCountLengths[1]))
            return FALSE;
        Run -> LatticeCountV = Run -> BodyCounts[1][0];
        return TRUE;
    }
    if (strcmp(Key, "lattice_count_w") == 0) {
        if (!OptParseCounts(Value,
                            &Run -> BodyCounts[2],
                            &Run -> BodyCountLengths[2]))
            return FALSE;
        Run -> LatticeCountW = Run -> BodyCounts[2][0];
        return TRUE;
    }
    if (strcmp(Key, "initial_thickness") == 0)
        return OptParseReal(Value, &Run -> InitialThickness);
    if (strcmp(Key, "mesh_sample") == 0)
        return OptParseInt(Value, &Run -> MeshSample);
    if (strcmp(Key, "max_iterations") == 0)
        return OptParseInt(Value, &Run -> MaxIterations);
    if (strcmp(Key, "convergence_iterations") == 0)
        return OptParseInt(Value, &Run -> ConvergenceIterations);
    if (strcmp(Key, "thickness_change_tolerance") == 0)
        return OptParseReal(Value, &Run -> ThicknessChangeTolerance);
    if (strcmp(Key, "target_volume_fraction") == 0) {
        if (OptIsBlank(Value))
            return TRUE;
        return OptParseReal(Value, &Run -> TargetVolumeFraction);
    }
    if (strcmp(Key, "target_volume_fraction_tol") == 0) {
        if (OptIsBlank(Value))
            return TRUE;
        return OptParseReal(Value, &Run -> TargetVolumeFractionTol);
    }
    if (strcmp(Key, "optimization_algorithm") == 0) {
        if (OptIsBlank(Value))
            return TRUE;
        return OptParseInt(Value, &Run -> OptimizationAlgorithm);
    }
    if (strcmp(Key, "root_folder") == 0 ||
        strcmp(Key, "base_folder") == 0) {
        snprintf(Cfg -> BaseFolder, sizeof(Cfg -> BaseFolder), "%s", Value);
        return TRUE;
    }
    if (strcmp(Key, "work_dir") == 0) {
        snprintf(Cfg -> WorkDir, sizeof(Cfg -> WorkDir), "%s", Value);
        return TRUE;
    }
    if (strcmp(Key, "ansys_exe") == 0) {
        snprintf(Cfg -> AnsysExe, sizeof(Cfg -> AnsysExe), "%s", Value);
        return TRUE;
    }
    if (strcmp(Key, "ansys_project") == 0) {
        snprintf(Cfg -> AnsysProject,
                 sizeof(Cfg -> AnsysProject),
                 "%s",
                 Value);
        return TRUE;
    }
    if (strcmp(Key, "ansys_journal") == 0) {
        snprintf(Cfg -> AnsysJournal,
                 sizeof(Cfg -> AnsysJournal),
                 "%s",
                 Value);
        return TRUE;
    }
    if (strcmp(Key, "itd_output_file") == 0) {
        snprintf(Cfg -> ItdOutputFile,
                 sizeof(Cfg -> ItdOutputFile),
                 "%s",
                 Value);
        return TRUE;
    }
    if (strcmp(Key, "sed_file") == 0) {
        snprintf(Cfg -> SedFile, sizeof(Cfg -> SedFile), "%s", Value);
        return TRUE;
    }
    if (strcmp(Key, "volume_file") == 0) {
        snprintf(Cfg -> VolumeFile, sizeof(Cfg -> VolumeFile), "%s", Value);
        return TRUE;
    }
    if (strcmp(Key, "data_new_file") == 0) {
        snprintf(Cfg -> DataNewFile, sizeof(Cfg -> DataNewFile), "%s", Value);
        return TRUE;
    }
    if (strcmp(Key, "data_prev_file") == 0) {
        snprintf(Cfg -> DataPrevFile,
                 sizeof(Cfg -> DataPrevFile),
                 "%s",
                 Value);
        return TRUE;
    }
    if (strcmp(Key, "design_parameters_file") == 0) {
        snprintf(Cfg -> DesignParametersFile,
                 sizeof(Cfg -> DesignParametersFile),
                 "%s",
                 Value);
        return TRUE;
    }
    if (strcmp(Key, "summary_file") == 0) {
        snprintf(Cfg -> SummaryFile,
                 sizeof(Cfg -> SummaryFile),
                 "%s",
                 Value);
        return TRUE;
    }

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
    if (strcmp(Key, "minimum_thickness") == 0)
        return OptParseReal(Value, &Cfg -> MinimumThickness);
    if (strcmp(Key, "safety_factor_sed") == 0)
        return OptParseReal(Value, &Cfg -> SafetyFactorSed);
    if (strcmp(Key, "remove_thin_tiles") == 0) {
        Cfg -> RemoveThinTiles = OptParseBool(Value);
        return TRUE;
    }

    return FALSE;
}

/****************************************************************************
* DESCRIPTION:                                                              *
* Parses a boolean configuration value.                                     *
****************************************************************************/
static CagdBType OptParseBool(const char *Value)
{
    if (Value == NULL)
        return FALSE;
    if (OPT_STRICMP(Value, "true") == 0 || strcmp(Value, "1") == 0 ||
        OPT_STRICMP(Value, "yes") == 0 || OPT_STRICMP(Value, "on") == 0)
        return TRUE;

    return FALSE;
}

/****************************************************************************
* DESCRIPTION:                                                              *
* Verifies the unified configuration after it has been parsed.              *
****************************************************************************/
static int ValidateInputs(const OptRunInputStruct *Run,
                          const OptConfigInputStruct *Cfg)
{
    if (OptIsBlank(Run -> MacroShapePath)) {
        fprintf(stderr, "Missing macro_shape_file in configuration file.\n");
        return FALSE;
    }
    {
        size_t Len = strlen(Run -> MacroShapePath);

        if (Len < 4 || OPT_STRICMP(Run -> MacroShapePath + Len - 4,
                                   ".itd") != 0) {
            fprintf(stderr, "macro_shape_file must name an .itd file.\n");
            return FALSE;
        }
    }
    if (Run -> SimulationType != 1 && Run -> SimulationType != 2) {
        fprintf(stderr,
                "simulation_type must be 1 (ANSYS) or 2 (CalculiX).\n");
        return FALSE;
    }
    if (Run -> LatticeType < 0 || Run -> LatticeType > 2 ||
        Run -> LatticeCountU <= 0 || Run -> LatticeCountV <= 0 ||
        Run -> LatticeCountW <= 0) {
        fprintf(stderr, "Invalid lattice parameters in configuration file.\n");
        return FALSE;
    }
    if (!isfinite(Run -> SkinThickness) || Run -> SkinThickness < 0.0 ||
        Run -> SkinThickness > 0.5) {
        fprintf(stderr, "skin_thickness must be between 0 and 0.5.\n");
        return FALSE;
    }
    if (Run -> SkinThickness > 0.0 && Run -> LatticeType != 0) {
        fprintf(stderr,
                "skin_thickness is supported only for lattice_type 0.\n");
        return FALSE;
    }
    if (Run -> BodyCountLengths[0] != Run -> BodyCountLengths[1] ||
        Run -> BodyCountLengths[0] != Run -> BodyCountLengths[2]) {
        fprintf(stderr, "U, V and W count lists must have equal lengths.\n");
        return FALSE;
    }
    if ((double) (Run -> MeshSample - 1) * (Run -> MeshSample - 1) *
        (Run -> MeshSample - 1) > INT_MAX) {
        fprintf(stderr, "mesh_sample is too large.\n");
        return FALSE;
    }
    if (Cfg -> Scale <= 0.0 || !isfinite(Cfg -> Scale))
        return FALSE;
    if (Run -> InitialThickness <= 0.0 || Run -> MeshSample <= 1 ||
        Run -> MaxIterations <= 0) {
        fprintf(stderr,
                "Invalid initial thickness, mesh sample, or max iterations "
                "in configuration file.\n");
        return FALSE;
    }
    if (Run -> ConvergenceIterations <= 0) {
        fprintf(stderr, "convergence_iterations must be at least 1.\n");
        return FALSE;
    }
    if (Run -> ThicknessChangeTolerance <= 0.0 ||
        !isfinite(Run -> ThicknessChangeTolerance)) {
        fprintf(stderr, "thickness_change_tolerance must be positive.\n");
        return FALSE;
    }
    if (OptIsBlank(Cfg -> BaseFolder) || OptIsBlank(Cfg -> WorkDir) ||
        OptIsBlank(Cfg -> InpFile) || OptIsBlank(Cfg -> ItdOutputFile) ||
        OptIsBlank(Cfg -> SedFile) || OptIsBlank(Cfg -> VolumeFile) ||
        OptIsBlank(Cfg -> DataNewFile) || OptIsBlank(Cfg -> DataPrevFile) ||
        OptIsBlank(Cfg -> DesignParametersFile) ||
        OptIsBlank(Cfg -> SummaryFile)) {
        fprintf(stderr, "Missing required paths in configuration file.\n");
        return FALSE;
    }
    if (Run -> SimulationType == 1 &&
        (OptIsBlank(Cfg -> AnsysExe) || OptIsBlank(Cfg -> AnsysProject) ||
         OptIsBlank(Cfg -> AnsysJournal))) {
        fprintf(stderr,
                "ANSYS paths are required when simulation_type is 1.\n");
        return FALSE;
    }
    return TRUE;
}

/***************************************************************************
*                         File and Process Helpers                          *
***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
* Creates every component of a directory path if it does not exist.       *
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
* DESCRIPTION:                                                             *
* Recursively copies iteration output while excluding active run files.    *
***************************************************************************/

static int CopyDirectoryRecursive(OptWorkflowStateStruct *Opt,
                                  const char *SourceDir,
                                  const char *TargetDir)
{
#ifdef _WIN32
    char Pattern[OPT_MAX_PATH];
    char ExePath[OPT_MAX_PATH];
    const char *SkipExe = NULL;
    DWORD ExeLen;
    WIN32_FIND_DATAA Fd;
    HANDLE HFind;

    EnsureDirectoryRecursive(TargetDir);
    snprintf(Pattern, sizeof(Pattern), "%s\\*", SourceDir);

    HFind = FindFirstFileA(Pattern, &Fd);
    if (HFind == INVALID_HANDLE_VALUE)
        return 0;

    ExeLen = GetModuleFileNameA(NULL,
                                ExePath,
                                (DWORD) sizeof(ExePath));
    if (ExeLen > 0 && ExeLen < (DWORD) sizeof(ExePath))
        SkipExe = OptPathFileName(ExePath);

    do {
        char SrcPath[OPT_MAX_PATH];
        char DstPath[OPT_MAX_PATH];
        const char *Name = Fd.cFileName;
        const char *SkipItd = OptPathFileName(Opt -> Config.ItdOutputFile);

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
            if (OPT_STRICMP(Name, SkipItd) == 0 ||
                (SkipExe != NULL && OPT_STRICMP(Name, SkipExe) == 0))
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
* DESCRIPTION:                                                             *
* Executes an external solver command and waits for its exit status.       *
***************************************************************************/

static int RunCommandLineBlocking(const char *CommandLine,
                                  const char *WorkingDir)
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
*   Runs a child process synchronously and redirects stdout and stderr to  *
* a fresh log file.                                                        *
*                                                                          *
* PARAMETERS:                                                              *
*   CommandLine: Command line to execute.                                  *
*   WorkingDir: Working directory of the child process.                    *
*   LogFile:    File that receives stdout and stderr.                      *
*                                                                          *
* RETURN VALUE:                                                            *
*   int: TRUE if the child process exits successfully, FALSE otherwise.    *
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
* Writes one text line, either replacing or appending to the target file.  *
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
*                              General Helpers                              *
***************************************************************************/

/***************************************************************************
* DESCRIPTION:                                                             *
* Records a fatal workflow reason, reports it, and terminates execution.   *
***************************************************************************/

static void OptFatal(const char *Msg)
{
    if (GlblActiveWorkflow != NULL)
        OptWriteFinishReason(GlblActiveWorkflow,
                             Msg,
                             GlblActiveWorkflow -> CurrentIteration);
    fprintf(stderr, "ERROR: %s\n", Msg);
    exit(EXIT_FAILURE);
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Formats a one-string-argument error before invoking the fatal handler.   *
***************************************************************************/

static void OptFatalF(const char *Fmt, const char *Arg)
{
    char Buffer[OPT_MAX_LINE];
    snprintf(Buffer, sizeof(Buffer), Fmt, Arg);
    OptFatal(Buffer);
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Trims leading and trailing ASCII whitespace from a mutable string.       *
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
* DESCRIPTION:                                                             *
* Tests whether a string is null, empty, or contains only whitespace.      *
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
* DESCRIPTION:                                                             *
* Removes an optional UTF-8 byte-order mark from the first input line.     *
***************************************************************************/

static void OptStripUtf8Bom(char *Line)
{
    unsigned char *Ptr = (unsigned char *)Line;
    if (Ptr != NULL && Ptr[0] == 0xEF && Ptr[1] == 0xBB && Ptr[2] == 0xBF) {
        memmove(Line, Line + 3, strlen(Line + 3) + 1);
    }
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Parses a complete base-10 integer and rejects overflow or trailing text. *
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
* DESCRIPTION:                                                             *
* Parses a finite real value and rejects overflow or trailing text.        *
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
* DESCRIPTION:                                                             *
* Joins two path components using the platform directory separator.       *
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
* DESCRIPTION:                                                             *
* Returns the final filename component of a Windows or POSIX path.         *
***************************************************************************/

static const char *OptPathFileName(const char *Path)
{
    const char *P1 = strrchr(Path, '\\');
    const char *P2 = strrchr(Path, '/');
    const char *P = P1 > P2 ? P1 : P2;
    return P != NULL ? P + 1 : Path;
}

/***************************************************************************
* DESCRIPTION:                                                             *
* Inserts an iteration index before a filename extension.                  *
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
* DESCRIPTION:                                                             *
* Checks whether a regular path can be opened for binary reading.          *
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
* DESCRIPTION:                                                             *
* Formats the current local time for an iteration summary record.          *
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
* DESCRIPTION:                                                             *
* Generates the timestamp-based directory name for one optimization run.  *
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

    if (Bytes > UINT_MAX)
        return NULL;

    Data = IritMalloc((unsigned int) Bytes);
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

    if (NewBytes > UINT_MAX)
        return NULL;

    NewData = (char *)IritMalloc((unsigned int) NewBytes);
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
