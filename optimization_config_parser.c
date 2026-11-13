#include "optimization_config_parser.h"

#include "inc_irit/irit_sm.h"
#include "inc_irit/misc_lib.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#   define CONFIG_PATH_SEP '\\'
#   define CONFIG_STRICMP _stricmp
#else
#   include <strings.h>
#   define CONFIG_PATH_SEP '/'
#   define CONFIG_STRICMP strcasecmp
#endif

#define CONFIG_LINE_MAX 2048
#define CONFIG_INVALID_REAL (-1.0e300)

typedef int (*ConfigEntryHandlerType)(const char *Key,
                                      const char *Value,
                                      void *Context);

static char *ConfigTrim(char *Text)
{
    char *End;

    while (isspace((unsigned char) *Text))
        ++Text;
    End = Text + strlen(Text);
    while (End > Text && isspace((unsigned char) End[-1]))
        --End;
    *End = '\0';
    return Text;
}

static int ConfigIsBlank(const char *Text)
{
    if (Text == NULL)
        return TRUE;
    while (*Text != '\0') {
        if (!isspace((unsigned char) *Text))
            return FALSE;
        ++Text;
    }
    return TRUE;
}

static void ConfigStripUtf8Bom(char *Line)
{
    const unsigned char *Bytes = (const unsigned char *) Line;

    if (Bytes[0] == 0xef && Bytes[1] == 0xbb && Bytes[2] == 0xbf)
        memmove(Line, Line + 3, strlen(Line + 3) + 1);
}

static int ConfigParseInt(const char *Text, int *Value)
{
    char *End;
    long Parsed;

    if (Text == NULL || Value == NULL)
        return FALSE;
    errno = 0;
    Parsed = strtol(Text, &End, 10);
    while (isspace((unsigned char) *End))
        ++End;
    if (errno == ERANGE || End == Text || *End != '\0' ||
        Parsed < INT_MIN || Parsed > INT_MAX)
        return FALSE;
    *Value = (int) Parsed;
    return TRUE;
}

static int ConfigParseReal(const char *Text, CagdRType *Value)
{
    char *End;
    double Parsed;

    if (Text == NULL || Value == NULL)
        return FALSE;
    errno = 0;
    Parsed = strtod(Text, &End);
    while (isspace((unsigned char) *End))
        ++End;
    if (errno == ERANGE || End == Text || *End != '\0' || !isfinite(Parsed))
        return FALSE;
    *Value = (CagdRType) Parsed;
    return TRUE;
}

static CagdBType ConfigParseBool(const char *Value)
{
    return Value != NULL &&
        (CONFIG_STRICMP(Value, "true") == 0 || strcmp(Value, "1") == 0 ||
         CONFIG_STRICMP(Value, "yes") == 0 ||
         CONFIG_STRICMP(Value, "on") == 0);
}

static void ConfigPathJoin(char *Destination,
                           size_t DestinationSize,
                           const char *Parent,
                           const char *Child)
{
    size_t ParentLength;

    if (DestinationSize == 0)
        return;
    if (ConfigIsBlank(Parent)) {
        snprintf(Destination, DestinationSize, "%s", Child != NULL ? Child : "");
        return;
    }
    if (ConfigIsBlank(Child)) {
        snprintf(Destination, DestinationSize, "%s", Parent);
        return;
    }
    ParentLength = strlen(Parent);
    if (Parent[ParentLength - 1] == '/' || Parent[ParentLength - 1] == '\\')
        snprintf(Destination, DestinationSize, "%s%s", Parent, Child);
    else
        snprintf(Destination, DestinationSize, "%s%c%s",
                 Parent, CONFIG_PATH_SEP, Child);
}

static void *ConfigReallocCopy(void *OldData,
                               size_t OldBytes,
                               size_t NewBytes)
{
    void *NewData;

    if (NewBytes > UINT_MAX)
        return NULL;
    NewData = IritMalloc((unsigned int) NewBytes);
    if (NewData == NULL)
        return NULL;
    if (OldData != NULL && OldBytes > 0)
        memcpy(NewData, OldData, OldBytes);
    if (OldData != NULL)
        IritFree(OldData);
    return NewData;
}

static int ConfigParseCountList(const char *Text, int **Values, int *Count)
{
    const char *Cursor = Text;
    int Length = 0, Capacity = 0;
    int *Result = NULL;

    if (Text == NULL || Values == NULL || Count == NULL)
        return FALSE;
    while (*Cursor != '\0') {
        char *End;
        long Parsed;

        while (isspace((unsigned char) *Cursor) || *Cursor == ',' ||
               *Cursor == ';' || *Cursor == 'x' || *Cursor == 'X')
            ++Cursor;
        if (*Cursor == '\0')
            break;
        errno = 0;
        Parsed = strtol(Cursor, &End, 10);
        if (End == Cursor || errno == ERANGE || Parsed <= 0 || Parsed > INT_MAX) {
            if (Result != NULL)
                IritFree(Result);
            return FALSE;
        }
        Cursor = End;
        if (*Cursor != '\0' && !isspace((unsigned char) *Cursor) &&
            *Cursor != ',' && *Cursor != ';' && *Cursor != 'x' &&
            *Cursor != 'X') {
            if (Result != NULL)
                IritFree(Result);
            return FALSE;
        }
        if (Length == Capacity) {
            int NewCapacity = Capacity == 0 ? 4 : Capacity * 2;
            int *NewResult = (int *) ConfigReallocCopy(
                Result, (size_t) Capacity * sizeof(int),
                (size_t) NewCapacity * sizeof(int));
            if (NewResult == NULL) {
                if (Result != NULL)
                    IritFree(Result);
                return FALSE;
            }
            Result = NewResult;
            Capacity = NewCapacity;
        }
        Result[Length++] = (int) Parsed;
    }
    if (Length == 0) {
        if (Result != NULL)
            IritFree(Result);
        return FALSE;
    }
    if (*Values != NULL)
        IritFree(*Values);
    *Values = Result;
    *Count = Length;
    return TRUE;
}

void OptimizationConfigInit(OptConfigStruct *Config)
{
    IRIT_ZAP_MEM(Config, sizeof(*Config));
    Config->SimulationType = 1;
    Config->LatticeType = -1;
    Config->LatticeCountU = -1;
    Config->LatticeCountV = -1;
    Config->LatticeCountW = -1;
    Config->BodyStackAxis = -1;
    Config->UseMicroStructFFD = TRUE;
    Config->InitialThickness = CONFIG_INVALID_REAL;
    Config->MeshSample = -1;
    Config->MaxIterations = -1;
    Config->ConvergenceIterations = 3;
    Config->ThicknessChangeTolerance = 2.5e-3;
    Config->TargetVolumeFraction = CONFIG_INVALID_REAL;
    Config->TargetVolumeFractionTol = CONFIG_INVALID_REAL;
    Config->OptimizationAlgorithm = 1;
    Config->UseTargetVolume = TRUE;
    Config->Scale = 0.015;
    Config->Density = 2700.0;
    Config->YoungsModulus = 6.9e10;
    Config->PoissonRatio = 0.33;
    Config->SedValue = 0.15;
    Config->MinimumThickness = 1.0e-5;
    Config->SafetyFactorSed = 1.1;
    snprintf(Config->AnsysExe, sizeof(Config->AnsysExe), "%s",
             "C:\\Program Files\\ANSYS Inc\\v222\\Framework\\bin\\Win64\\RunWB2.exe");
}

void OptimizationConfigFree(OptConfigStruct *Config)
{
    int Direction;

    if (Config == NULL)
        return;
    for (Direction = 0; Direction < 3; ++Direction) {
        if (Config->LatticeCountLists[Direction] != NULL)
            IritFree(Config->LatticeCountLists[Direction]);
        Config->LatticeCountLists[Direction] = NULL;
        Config->LatticeCountListLengths[Direction] = 0;
    }
}

static int ConfigResolvePaths(OptConfigStruct *Config)
{
    char AnsysDirectory[IRIT_OPT_MAX_PATH];
    char MacroShapeFile[IRIT_OPT_MAX_PATH];

    if (ConfigIsBlank(Config->BaseFolder)) {
        fprintf(stderr, "Missing root_folder in configuration file.\n");
        return FALSE;
    }
    if (ConfigIsBlank(Config->WorkDir))
        ConfigPathJoin(Config->WorkDir, sizeof(Config->WorkDir),
                       Config->BaseFolder, "TEMP_SET");
    ConfigPathJoin(AnsysDirectory, sizeof(AnsysDirectory),
                   Config->BaseFolder, "ANSYS_SET");
    if (ConfigIsBlank(Config->MacroShapePath))
        snprintf(Config->MacroShapePath, sizeof(Config->MacroShapePath),
                 "%s", "DefMap.itd");
    snprintf(MacroShapeFile, sizeof(MacroShapeFile), "%s", Config->MacroShapePath);
    ConfigPathJoin(Config->MacroShapePath, sizeof(Config->MacroShapePath),
                   Config->WorkDir, MacroShapeFile);

#define CONFIG_DEFAULT_PATH(Field, Parent, Name) \
    do { if (ConfigIsBlank(Field)) \
        ConfigPathJoin(Field, sizeof(Field), Parent, Name); } while (0)
    CONFIG_DEFAULT_PATH(Config->AnsysProject, AnsysDirectory, "WB2.wbpj");
    CONFIG_DEFAULT_PATH(Config->AnsysJournal, AnsysDirectory, "RunWB.wbjn");
    ConfigPathJoin(Config->InpFile, sizeof(Config->InpFile),
                   Config->WorkDir, "optimized_geo.inp");
    CONFIG_DEFAULT_PATH(Config->ItdOutputFile, Config->WorkDir, "optimized_geo.itd");
    CONFIG_DEFAULT_PATH(Config->SedFile, Config->WorkDir, "SED.txt");
    CONFIG_DEFAULT_PATH(Config->VolumeFile, Config->WorkDir, "VOLUME.txt");
    CONFIG_DEFAULT_PATH(Config->DataNewFile, Config->WorkDir, "data_new.csv");
    CONFIG_DEFAULT_PATH(Config->DataPrevFile, Config->WorkDir, "data_prev.csv");
    CONFIG_DEFAULT_PATH(Config->DesignParametersFile, Config->WorkDir,
                        "DesignParameters.csv");
    if (ConfigIsBlank(Config->SummaryFile))
        snprintf(Config->SummaryFile, sizeof(Config->SummaryFile), "%s", "Summary.csv");
#undef CONFIG_DEFAULT_PATH
    return TRUE;
}

static int ConfigHandleEntry(const char *Key, const char *Value, void *Context)
{
    OptConfigStruct *Config = (OptConfigStruct *) Context;

#define CONFIG_STRING(KeyName, Field) \
    if (strcmp(Key, KeyName) == 0) { \
        snprintf(Config->Field, sizeof(Config->Field), "%s", Value); return TRUE; }
#define CONFIG_INT(KeyName, Field) \
    if (strcmp(Key, KeyName) == 0) return ConfigParseInt(Value, &Config->Field)
#define CONFIG_REAL(KeyName, Field) \
    if (strcmp(Key, KeyName) == 0) return ConfigParseReal(Value, &Config->Field)

    CONFIG_STRING("macro_shape_file", MacroShapePath);
    CONFIG_INT("simulation_type", SimulationType);
    CONFIG_INT("lattice_type", LatticeType);
    CONFIG_REAL("skin_thickness", SkinThickness);
    if (strcmp(Key, "lattice_count_u") == 0) {
        if (!ConfigParseCountList(Value, &Config->LatticeCountLists[0],
                                  &Config->LatticeCountListLengths[0])) return FALSE;
        Config->LatticeCountU = Config->LatticeCountLists[0][0]; return TRUE;
    }
    if (strcmp(Key, "lattice_count_v") == 0) {
        if (!ConfigParseCountList(Value, &Config->LatticeCountLists[1],
                                  &Config->LatticeCountListLengths[1])) return FALSE;
        Config->LatticeCountV = Config->LatticeCountLists[1][0]; return TRUE;
    }
    if (strcmp(Key, "lattice_count_w") == 0) {
        if (!ConfigParseCountList(Value, &Config->LatticeCountLists[2],
                                  &Config->LatticeCountListLengths[2])) return FALSE;
        Config->LatticeCountW = Config->LatticeCountLists[2][0]; return TRUE;
    }
    if (strcmp(Key, "body_stack_axis") == 0) {
        if (CONFIG_STRICMP(Value, "u") == 0) Config->BodyStackAxis = 0;
        else if (CONFIG_STRICMP(Value, "v") == 0) Config->BodyStackAxis = 1;
        else if (CONFIG_STRICMP(Value, "w") == 0) Config->BodyStackAxis = 2;
        else return FALSE;
        return TRUE;
    }
    if (strcmp(Key, "use_microstructure_ffd") == 0) {
        Config->UseMicroStructFFD = ConfigParseBool(Value); return TRUE;
    }
    CONFIG_REAL("initial_thickness", InitialThickness);
    CONFIG_INT("mesh_sample", MeshSample);
    CONFIG_INT("max_iterations", MaxIterations);
    CONFIG_INT("convergence_iterations", ConvergenceIterations);
    CONFIG_REAL("thickness_change_tolerance", ThicknessChangeTolerance);
    if (strcmp(Key, "use_target_volume") == 0) {
        Config->UseTargetVolume = ConfigParseBool(Value); return TRUE;
    }
    if (strcmp(Key, "target_volume_fraction") == 0)
        return ConfigIsBlank(Value) || ConfigParseReal(Value, &Config->TargetVolumeFraction);
    if (strcmp(Key, "target_volume_fraction_tol") == 0)
        return ConfigIsBlank(Value) || ConfigParseReal(Value, &Config->TargetVolumeFractionTol);
    if (strcmp(Key, "optimization_algorithm") == 0)
        return ConfigIsBlank(Value) || ConfigParseInt(Value, &Config->OptimizationAlgorithm);
    if (strcmp(Key, "root_folder") == 0 || strcmp(Key, "base_folder") == 0) {
        snprintf(Config->BaseFolder, sizeof(Config->BaseFolder), "%s", Value); return TRUE;
    }
    CONFIG_STRING("work_dir", WorkDir);
    CONFIG_STRING("ansys_exe", AnsysExe);
    CONFIG_STRING("ansys_project", AnsysProject);
    CONFIG_STRING("ansys_journal", AnsysJournal);
    CONFIG_STRING("itd_output_file", ItdOutputFile);
    CONFIG_STRING("sed_file", SedFile);
    CONFIG_STRING("volume_file", VolumeFile);
    CONFIG_STRING("data_new_file", DataNewFile);
    CONFIG_STRING("data_prev_file", DataPrevFile);
    CONFIG_STRING("design_parameters_file", DesignParametersFile);
    CONFIG_STRING("summary_file", SummaryFile);
    CONFIG_REAL("scale", Scale);
    CONFIG_REAL("density", Density);
    CONFIG_REAL("youngs_modulus", YoungsModulus);
    CONFIG_REAL("poisson_ratio", PoissonRatio);
    CONFIG_REAL("sed_value", SedValue);
    CONFIG_REAL("minimum_thickness", MinimumThickness);
    CONFIG_REAL("safety_factor_sed", SafetyFactorSed);
    if (strcmp(Key, "remove_thin_tiles") == 0) {
        Config->RemoveThinTiles = ConfigParseBool(Value); return TRUE;
    }
#undef CONFIG_STRING
#undef CONFIG_INT
#undef CONFIG_REAL
    return FALSE;
}

static int ConfigParseFile(const char *FileName,
                           ConfigEntryHandlerType Handler,
                           void *Context)
{
    FILE *File = fopen(FileName, "r");
    char Line[CONFIG_LINE_MAX];
    int LineNumber = 0;

    if (File == NULL)
        return FALSE;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        char *Equals, *Key, *Value;
        ++LineNumber;
        ConfigStripUtf8Bom(Line);
        Key = ConfigTrim(Line);
        if (*Key == '\0' || *Key == '#')
            continue;
        Equals = strchr(Key, '=');
        if (Equals == NULL)
            continue;
        *Equals = '\0';
        Value = ConfigTrim(Equals + 1);
        Key = ConfigTrim(Key);
        if (!Handler(Key, Value, Context)) {
            fprintf(stderr, "Invalid key or value in %s at line %d: %s\n",
                    FileName, LineNumber, Key);
            fclose(File);
            return FALSE;
        }
    }
    fclose(File);
    return TRUE;
}

int OptimizationConfigValidate(const OptConfigStruct *Config)
{
    size_t Length;

    if (Config == NULL || Config->BodyStackAxis < -1 || Config->BodyStackAxis > 2)
        return FALSE;
    if (ConfigIsBlank(Config->MacroShapePath))
        return FALSE;
    Length = strlen(Config->MacroShapePath);
    if (Length < 4 || CONFIG_STRICMP(Config->MacroShapePath + Length - 4, ".itd") != 0)
        return FALSE;
    if (Config->SimulationType != 1 && Config->SimulationType != 2)
        return FALSE;
    if (Config->LatticeType < 0 || Config->LatticeType > 2 ||
        Config->LatticeCountU <= 0 || Config->LatticeCountV <= 0 ||
        Config->LatticeCountW <= 0)
        return FALSE;
    if (!isfinite(Config->SkinThickness) || Config->SkinThickness < 0.0 ||
        Config->SkinThickness > 0.5 ||
        (Config->SkinThickness > 0.0 && Config->LatticeType != 0))
        return FALSE;
    if (Config->LatticeCountListLengths[0] != Config->LatticeCountListLengths[1] ||
        Config->LatticeCountListLengths[0] != Config->LatticeCountListLengths[2])
        return FALSE;
    if (Config->MeshSample <= 1 ||
        (double) (Config->MeshSample - 1) * (Config->MeshSample - 1) *
            (Config->MeshSample - 1) > INT_MAX)
        return FALSE;
    if (Config->Scale <= 0.0 || !isfinite(Config->Scale) ||
        Config->RemoveThinTiles || Config->InitialThickness <= 0.0 ||
        Config->MaxIterations <= 0 || Config->ConvergenceIterations <= 0 ||
        Config->ThicknessChangeTolerance <= 0.0 ||
        !isfinite(Config->ThicknessChangeTolerance))
        return FALSE;
    if (ConfigIsBlank(Config->BaseFolder) || ConfigIsBlank(Config->WorkDir) ||
        ConfigIsBlank(Config->InpFile) || ConfigIsBlank(Config->ItdOutputFile) ||
        ConfigIsBlank(Config->SedFile) || ConfigIsBlank(Config->VolumeFile) ||
        ConfigIsBlank(Config->DataNewFile) || ConfigIsBlank(Config->DataPrevFile) ||
        ConfigIsBlank(Config->DesignParametersFile) ||
        ConfigIsBlank(Config->SummaryFile))
        return FALSE;
    if (Config->SimulationType == 1 &&
        (ConfigIsBlank(Config->AnsysExe) || ConfigIsBlank(Config->AnsysProject) ||
         ConfigIsBlank(Config->AnsysJournal)))
        return FALSE;
    return TRUE;
}

int OptimizationConfigLoadFile(const char *FileName, OptConfigStruct *Config)
{
    if (FileName == NULL || Config == NULL)
        return FALSE;
    OptimizationConfigInit(Config);
    if (!ConfigParseFile(FileName, ConfigHandleEntry, Config) ||
        !ConfigResolvePaths(Config) || !OptimizationConfigValidate(Config)) {
        OptimizationConfigFree(Config);
        return FALSE;
    }
    return TRUE;
}
