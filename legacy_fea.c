/* Legacy FEA compatibility implementation.
 *
 * This file is intentionally included by test.c after its private workflow
 * types and helper declarations.  Do not compile it as a separate translation
 * unit.  Delete this file and the include when CalculiX and temporary INP
 * augmentation are retired. */

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

    if (!CreateProcessA(NULL, CmdMutable, NULL, NULL, FALSE,
                        CREATE_NO_WINDOW, NULL, WorkingDir, &Si, &Pi)) {
        fprintf(stderr, "CreateProcess failed for command:\n%s\n",
                CommandLine);
        return FALSE;
    }
    WaitForSingleObject(Pi.hProcess, INFINITE);
    GetExitCodeProcess(Pi.hProcess, &ExitCode);
    CloseHandle(Pi.hThread);
    CloseHandle(Pi.hProcess);
    return ExitCode == 0;
#else
    (void) WorkingDir;
    return system(CommandLine) == 0;
#endif
}

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
    LogHandle = CreateFileA(LogFile, GENERIC_WRITE,
                            FILE_SHARE_READ | FILE_SHARE_WRITE, &Sa,
                            CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
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
    if (!CreateProcessA(NULL, CmdMutable, NULL, NULL, TRUE,
                        CREATE_NO_WINDOW, NULL, WorkingDir, &Si, &Pi)) {
        fprintf(stderr, "CreateProcess failed for command:\n%s\n",
                CommandLine);
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
    (void) WorkingDir;
    snprintf(Redirected, sizeof(Redirected), "%s > \"%s\" 2>&1",
             CommandLine, LogFile);
    return system(Redirected) == 0;
#endif
}

static void DefineINPFaceSelections(IritPrsrObjectStruct *Tile,
                                    int Support, int SupportSide,
                                    int Pressure)
{
    int I;
    IritTrivTVStruct *TV;

    if (Tile == NULL)
        return;
    if (IRIT_PRSR_IS_OLST_OBJ(Tile)) {
        IritPrsrObjectStruct *Child;

        for (I = 0; (Child = IritPrsrListObjectGet(Tile, I)) != NULL; ++I)
            DefineINPFaceSelections(Child, Support, SupportSide, Pressure);
    }
    else if (IRIT_PRSR_IS_TRIVAR_OBJ(Tile)) {
        for (TV = Tile -> U.Trivars; TV != NULL; TV = TV -> Pnext) {
            IritCagdSrfStruct *Faces[6];
            char Bndry[256] = { 0 };

            IritTrivBndrySrfsFromTVToData(TV, FALSE, FALSE, Faces);
            for (I = 0; I < 6; ++I) {
                IritCagdBBoxStruct BBox;
                const char *Name = NULL;

                IritCagdSrfBBox(Faces[I], &BBox);
                if (Support) {
                    CagdBType IsSupportFace = FALSE;

                    if (SupportSide == OPT_SUPPORT_SIDE_VMIN) {
                        IsSupportFace =
                            fabs(BBox.Min[1]) <= IRIT_UEPS &&
                            fabs(BBox.Max[1]) <= IRIT_UEPS;
                    }
                    else if (SupportSide == OPT_SUPPORT_SIDE_VMAX) {
                        IsSupportFace =
                            fabs(BBox.Min[1] - 1.0) <= IRIT_UEPS &&
                            fabs(BBox.Max[1] - 1.0) <= IRIT_UEPS;
                    }
                    if (IsSupportFace)
                        Name = "SUPPORT";
                }

                if (Name == NULL &&
                    Pressure && fabs(BBox.Min[2]) <= IRIT_UEPS &&
                    fabs(BBox.Max[2]) <= IRIT_UEPS)
                    Name = "PRESSURE";
                if (Name != NULL)
                    sprintf(Bndry + strlen(Bndry), "%s%d::%s",
                            Bndry[0] ? "," : "", I, Name);
                IritCagdSrfFree(Faces[I]);
            }
            IritMiscAttrIDSetStrAttrib(&TV -> Attr,
                                      IRIT_ATTR_CREATE_ID(Bndry), Bndry);
        }
    }
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
    switch (Opt -> Config.SimulationType) {
        case 1:
            /* Winglet analysis definition.  To return to the cantilever test,
             * replace the call below with AddAnsysAnalysisDefinition(Opt). */
            if (!AddAnsysAnalysisDefinitionWinglet(Opt)) {
                fprintf(stderr,
                    "Failed to add ANSYS winglet analysis definition to INP.\n");
                return FALSE;
            }
            return RunAnsysBatch(Opt);
        case 2:
            return RunCalculixBatch(Opt);
        default:
            fprintf(stderr,
                    "Unsupported simulation_type: %d\n",
                    Opt -> Config.SimulationType);
            return FALSE;
    }
}

/***************************************************************************
* DESCRIPTION:
*   Adds an ANSYS External-Model-compatible ABAQUS analysis definition
*   to the IRIT-generated INP file.
*
*   The original IRIT mesh is preserved.  The function adds:
*       - FIXED node set
*       - LOADED node set
*       - Material definition
*       - Material assignment to EL-ALL
*       - Boundary conditions
*       - Nodal force
*       - One analysis/load step
*
*   The generated INP remains an ABAQUS-format INP that can be imported
*   through ANSYS Workbench External Model.
*
* RETURN VALUE:
*   TRUE on success, FALSE on failure.
***************************************************************************/
static int AddAnsysAnalysisDefinition(OptWorkflowStateStruct* Opt)
{
    const CagdRType TotalLoadY = -1000.0;  /* Total load [N]. */

    char TempInp[OPT_MAX_PATH];
    char Line[OPT_MAX_LINE];

    FILE* Src = NULL;
    FILE* Dst = NULL;

    CagdRType MinZ = 1.0e300;
    CagdRType MaxZ = -1.0e300;
    CagdRType ZTolerance;

    int FixedCount = 0;
    int LoadedCount = 0;

    int I;
    int N;

    if (Opt == NULL)
        return FALSE;

    if (!Opt->Lattice.FEMTopologyBound ||
        Opt->FEM.NodeCount <= 0) {
        fprintf(stderr,
            "ANSYS INP: FEM topology is unavailable.\n");
        return FALSE;
    }

    /***********************************************************************
     * Determine model Z limits.
     ***********************************************************************/
    for (I = 0; I < Opt->FEM.NodeCount; ++I) {
        OptNodeStruct* Node = &Opt->FEM.Nodes[I];

        MinZ = IRIT_MIN(MinZ, Node->Z);
        MaxZ = IRIT_MAX(MaxZ, Node->Z);
    }

    ZTolerance =
        IRIT_MAX(1.0, fabs(MaxZ - MinZ)) * 1.0e-9;

    /***********************************************************************
     * Count support and loaded nodes.
     *
     * Same physical definition currently used by the CalculiX solver:
     *
     *   FIXED  = minimum-Z face.
     *
     *   LOADED = maximum-Z face,
     *            restricted to -0.015 < Y < 0.015.
     ***********************************************************************/
    for (I = 0; I < Opt->FEM.NodeCount; ++I) {
        OptNodeStruct* Node = &Opt->FEM.Nodes[I];

        if (fabs(Node->Z - MinZ) <= ZTolerance)
            ++FixedCount;

        if (fabs(Node->Z - MaxZ) <= ZTolerance &&
            Node->Y > -0.015 &&
            Node->Y < 0.015)
            ++LoadedCount;
    }

    if (FixedCount == 0 || LoadedCount == 0) {
        fprintf(stderr,
            "ANSYS INP: failed to select BC nodes "
            "(fixed=%d, loaded=%d).\n",
            FixedCount,
            LoadedCount);

        return FALSE;
    }

    /***********************************************************************
     * Rewrite through a temporary file.
     *
     * This also makes the function safe if it is accidentally called twice:
     * an existing ANSYS block is removed before a new one is written.
     ***********************************************************************/
    snprintf(TempInp,
        sizeof(TempInp),
        "%s.ansys.tmp",
        Opt->Config.InpFile);

    Src = fopen(Opt->Config.InpFile, "r");
    if (Src == NULL) {
        fprintf(stderr,
            "ANSYS INP: cannot open input file: %s\n",
            Opt->Config.InpFile);
        return FALSE;
    }

    Dst = fopen(TempInp, "w");
    if (Dst == NULL) {
        fprintf(stderr,
            "ANSYS INP: cannot create temporary file: %s\n",
            TempInp);
        fclose(Src);
        remove(TempInp);
        return FALSE;
    }

    /***********************************************************************
     * Copy original IRIT INP.
     ***********************************************************************/
    while (fgets(Line, sizeof(Line), Src) != NULL) {
        char Copy[OPT_MAX_LINE];
        char* P;

        snprintf(Copy, sizeof(Copy), "%s", Line);
        P = OptTrim(Copy);

        /*
         * If the file already contains our block, discard the old block.
         * The block is always appended at the end of the file.
         */
        if (strcmp(P,
            "** BEGIN ADDED ANSYS ANALYSIS DEFINITION") == 0)
            break;

        fputs(Line, Dst);
    }

    fclose(Src);
    Src = NULL;

    /***********************************************************************
     * Begin ANSYS-compatible ABAQUS analysis data.
     ***********************************************************************/
    fprintf(Dst,
        "\n"
        "** BEGIN ADDED ANSYS ANALYSIS DEFINITION\n");

    fprintf(Dst,
        "** Generated automatically for ANSYS External Model\n");

    fprintf(Dst,
        "** FIXED: minimum Z\n"
        "** LOADED: maximum Z and -0.015 < Y < 0.015\n");

    /***********************************************************************
     * Fixed node set.
     ***********************************************************************/
    fprintf(Dst, "*NSET, NSET=FIXED\n");

    N = 0;

    for (I = 0; I < Opt->FEM.NodeCount; ++I) {
        OptNodeStruct* Node = &Opt->FEM.Nodes[I];

        if (fabs(Node->Z - MinZ) <= ZTolerance) {

            ++N;

            fprintf(Dst,
                "%d%s",
                Node->Id,
                (N % 16 == 0 || N == FixedCount)
                ? "\n"
                : ",");
        }
    }

    /***********************************************************************
     * Loaded node set.
     ***********************************************************************/
    fprintf(Dst, "*NSET, NSET=LOADED\n");

    N = 0;

    for (I = 0; I < Opt->FEM.NodeCount; ++I) {
        OptNodeStruct* Node = &Opt->FEM.Nodes[I];

        if (fabs(Node->Z - MaxZ) <= ZTolerance &&
            Node->Y > -0.015 &&
            Node->Y < 0.015) {

            ++N;

            fprintf(Dst,
                "%d%s",
                Node->Id,
                (N % 16 == 0 || N == LoadedCount)
                ? "\n"
                : ",");
        }
    }

    /***********************************************************************
     * Material.
     ***********************************************************************/
    fprintf(Dst,
        "*MATERIAL, NAME=OPT_MATERIAL\n");

    fprintf(Dst,
        "*ELASTIC\n"
        "%.17g, %.17g\n",
        (double)Opt->Config.YoungsModulus,
        (double)Opt->Config.PoissonRatio);

    fprintf(Dst,
        "*DENSITY\n"
        "%.17g\n",
        (double)Opt->Config.Density);

    /***********************************************************************
     * Assign material to all IRIT elements.
     *
     * IRIT already creates:
     *
     *   *ELEMENT, TYPE=C3D8, ELSET=EL-ALL
     ***********************************************************************/
    fprintf(Dst,
        "*SOLID SECTION, ELSET=EL-ALL, MATERIAL=OPT_MATERIAL\n");

    /*
     * Standard ABAQUS solid-section data line.
     */
    fprintf(Dst, ",\n");

    /***********************************************************************
     * Analysis/load step.
     *
     * DO NOT add "*STATIC, SOLVER=SPOOLES".
     *
     * ANSYS External Model recognizes STEP and imports the BC/load data;
     * the actual analysis environment is your Workbench Static Structural
     * system.
     ***********************************************************************/
    fprintf(Dst, "*STEP\n");

    /***********************************************************************
     * Fully constrain the minimum-Z nodes:
     *
     * DOF 1 = X
     * DOF 2 = Y
     * DOF 3 = Z
     ***********************************************************************/
    fprintf(Dst,
        "*BOUNDARY\n"
        "FIXED, 1, 3, 0.0\n");

    /***********************************************************************
     * Apply total Y load, distributed equally over all loaded nodes.
     ***********************************************************************/
    fprintf(Dst,
        "*CLOAD\n"
        "LOADED, 2, %.17g\n",
        (double)(TotalLoadY /
            (CagdRType)LoadedCount));

    fprintf(Dst, "*END STEP\n");

    fprintf(Dst,
        "** END ADDED ANSYS ANALYSIS DEFINITION\n");

    fclose(Dst);
    Dst = NULL;

    /***********************************************************************
     * Replace original INP with completed INP.
     ***********************************************************************/
#ifdef _WIN32

    if (!MoveFileExA(TempInp,
        Opt->Config.InpFile,
        MOVEFILE_REPLACE_EXISTING |
        MOVEFILE_COPY_ALLOWED)) {

        fprintf(stderr,
            "ANSYS INP: failed to replace original INP.\n");

        remove(TempInp);
        return FALSE;
    }

#else

    if (rename(TempInp,
        Opt->Config.InpFile) != 0) {

        fprintf(stderr,
            "ANSYS INP: failed to replace original INP.\n");

        remove(TempInp);
        return FALSE;
    }

#endif

    fprintf(stderr,
        "ANSYS analysis definition added to INP "
        "(fixed=%d, loaded=%d, total Fy=%.6g N).\n",
        FixedCount,
        LoadedCount,
        (double)TotalLoadY);

    return TRUE;
}

/***************************************************************************
* DESCRIPTION:
*   Adds the hard-coded ANSYS External-Model analysis definition used by
*   the winglet test case to the IRIT-generated ABAQUS-format INP file.
*
*   The IRIT geometry callback is responsible for exporting the named
*   selections used here:
*
*       SUPPORT  - root boundary at the selected V-min face.
*       PRESSURE - loaded cover at the selected W-min face.
*
*   This function deliberately does not reconstruct those selections from
*   deformed XYZ coordinates.  It only references the UVW-based selections
*   already embedded in the INP by DefineINPFaceSelections().
*
*   Hard-coded winglet test properties:
*       PA 2200 elastic modulus     = 1.8 GPa
*       Poisson ratio               = 0.39
*       Density                     = 930 kg/m^3
*       Reference strength          = 45 MPa (comment only; linear solve)
*       Applied normal pressure     = 10 kPa
*
*   The pressure value is a test value and can be changed at the top of
*   this function.  Positive ABAQUS pressure acts normal to the selected
*   surface according to its surface orientation.
*
* RETURN VALUE:
*   TRUE on success, FALSE on failure.
***************************************************************************/
static int AbortTemporaryInp(FILE *Src, FILE *Dst, void *Allocation,
                             const char *TempInp)
{
    if (Src != NULL)
        fclose(Src);
    if (Dst != NULL)
        fclose(Dst);
    if (Allocation != NULL)
        IritFree(Allocation);
    remove(TempInp);
    return FALSE;
}

static int AddAnsysAnalysisDefinitionWinglet(OptWorkflowStateStruct* Opt)
{
    const CagdRType YoungsModulusPa = 1.8e9;
    const CagdRType PoissonRatio = 0.39;
    const CagdRType DensityKgM3 = 930.0;
    const CagdRType ReferenceStrengthPa = 45.0e6;
    const CagdRType PressurePa = 1.0e4;  /* 10 kPa test pressure. */
    const char *SupportSideLabel =
        OPT_SUPPORT_SIDE == OPT_SUPPORT_SIDE_VMAX ? "V-max" : "V-min";
    static const int C3D8FaceNodes[6][4] = {
        { 0, 1, 2, 3 },
        { 4, 7, 6, 5 },
        { 0, 4, 5, 1 },
        { 1, 5, 6, 2 },
        { 2, 6, 7, 3 },
        { 3, 7, 4, 0 }
    };

    char TempInp[OPT_MAX_PATH];
    char Line[OPT_MAX_LINE];
    FILE* Src = NULL;
    FILE* Dst = NULL;
    int* SupportNodeIds = NULL;
    int SupportNodeCount = 0;
    int SupportFaceCount = 0;
    int InSupportSurface = FALSE;
    int HasSupportSurface = FALSE;
    int HasPressureSurface = FALSE;
    int I, J, K;

    if (Opt == NULL || OptIsBlank(Opt->Config.InpFile))
        return FALSE;

    /* Convert the exported SUPPORT element faces through the persistent
     * C3D8 connectivity, independently of skin thickness and XYZ shape. */
    if (!Opt->Lattice.FEMTopologyBound ||
        Opt->FEM.NodeCount <= 0 ||
        Opt->FEM.ElementCount <= 0) {
        fprintf(stderr,
                "ANSYS winglet INP: FEM topology is unavailable.\n");
        return FALSE;
    }

    SupportNodeIds = (int*) OptAllocZero((size_t) Opt->FEM.NodeCount,
                                         sizeof(int));
    if (SupportNodeIds == NULL) {
        fprintf(stderr,
                "ANSYS winglet INP: failed to allocate SUPPORT node list.\n");
        return FALSE;
    }

    snprintf(TempInp,
             sizeof(TempInp),
             "%s.ansys_winglet.tmp",
             Opt->Config.InpFile);

    Src = fopen(Opt->Config.InpFile, "r");
    if (Src == NULL) {
        fprintf(stderr,
                "ANSYS winglet INP: cannot open input file: %s\n",
                Opt->Config.InpFile);
        return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
    }

    Dst = fopen(TempInp, "w");
    if (Dst == NULL) {
        fprintf(stderr,
                "ANSYS winglet INP: cannot create temporary file: %s\n",
                TempInp);
        return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
    }

    /* Copy the IRIT file while converting its UVW-generated SUPPORT surface
     * into the node set required by *BOUNDARY.  PRESSURE remains a surface. */
    while (fgets(Line, sizeof(Line), Src) != NULL) {
        char Copy[OPT_MAX_LINE];
        char* P;

        snprintf(Copy, sizeof(Copy), "%s", Line);
        P = OptTrim(Copy);

        if (strcmp(P, "** BEGIN ADDED ANSYS WINGLET ANALYSIS DEFINITION") == 0 ||
            strcmp(P, "** BEGIN ADDED ANSYS ANALYSIS DEFINITION") == 0)
            break;

        if (strncmp(P, "**", 2) == 0) {
            fputs(Line, Dst);
            continue;
        }

        if (*P == '*') {
            InSupportSurface = FALSE;
            if (OPT_STRNICMP(P, "*SURFACE", 8) == 0) {
                if (strstr(P, "SUPPORT") != NULL) {
                    HasSupportSurface = TRUE;
                    InSupportSurface = TRUE;
                }
                if (strstr(P, "PRESSURE") != NULL)
                    HasPressureSurface = TRUE;
            }
        }
        else if (InSupportSurface && *P != '\0') {
            int ElementId = 0;
            int FaceNumber = 0;
            char FacePrefix = '\0';
            OptElementStruct* Elem = NULL;

            if (sscanf(P, "%d , %c%d",
                       &ElementId, &FacePrefix, &FaceNumber) != 3 ||
                toupper((unsigned char) FacePrefix) != 'S' ||
                FaceNumber < 1 || FaceNumber > 6) {
                fprintf(stderr,
                        "ANSYS winglet INP: invalid SUPPORT surface row: %s\n",
                        P);
                return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
            }

            for (I = 0; I < Opt->FEM.ElementCount; ++I) {
                if (Opt->FEM.Elements[I].Id == ElementId) {
                    Elem = &Opt->FEM.Elements[I];
                    break;
                }
            }
            if (Elem == NULL || Elem->NodeCount != OPT_C3D8_NODE_COUNT) {
                fprintf(stderr,
                        "ANSYS winglet INP: SUPPORT references unknown/non-C3D8 "
                        "element %d.\n",
                        ElementId);
                return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
            }

            ++SupportFaceCount;
            for (K = 0; K < 4; ++K) {
                int NodeId = Elem->NodeIds[
                    C3D8FaceNodes[FaceNumber - 1][K]];
                int AlreadyPresent = FALSE;

                for (J = 0; J < SupportNodeCount; ++J) {
                    if (SupportNodeIds[J] == NodeId) {
                        AlreadyPresent = TRUE;
                        break;
                    }
                }
                if (!AlreadyPresent) {
                    if (SupportNodeCount >= Opt->FEM.NodeCount) {
                        fprintf(stderr,
                                "ANSYS winglet INP: SUPPORT node list overflow.\n");
                        return AbortTemporaryInp(Src, Dst, SupportNodeIds,
                                                 TempInp);
                    }
                    SupportNodeIds[SupportNodeCount++] = NodeId;
                }
            }
        }

        fputs(Line, Dst);
    }

    fclose(Src);
    Src = NULL;

    if (!HasSupportSurface || SupportFaceCount <= 0 || SupportNodeCount <= 0) {
        fprintf(stderr,
                "ANSYS winglet INP: named element surface SUPPORT was not "
                "found or contained no C3D8 faces.\n"
                "Check DefineINPFaceSelections() and the exported INP.\n");
        return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
    }

    if (!HasPressureSurface) {
        fprintf(stderr,
                "ANSYS winglet INP: named surface PRESSURE was not found.\n"
                "Check DefineINPFaceSelections() and the exported INP.\n");
        return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
    }

    fprintf(Dst, "\n** BEGIN ADDED ANSYS WINGLET ANALYSIS DEFINITION\n");
    fprintf(Dst,
            "** Winglet test: SUPPORT=%s of body %d, PRESSURE=W-min\n",
            SupportSideLabel,
            OPT_SUPPORT_BODY_INDEX);
    fprintf(Dst,
            "** SUPPORT_NODES derived from %d exported SUPPORT C3D8 faces "
            "(%d unique nodes)\n",
            SupportFaceCount,
            SupportNodeCount);
    fprintf(Dst, "** Material: EOS PA 2200 / Nylon 12 test definition\n");
    fprintf(Dst,
            "** Reference strength for post-processing only: %.17g Pa\n",
            (double) ReferenceStrengthPa);
    fprintf(Dst,
            "** Applied normal pressure: %.17g Pa\n",
            (double) PressurePa);

    fprintf(Dst, "*NSET, NSET=SUPPORT_NODES\n");
    for (I = 0; I < SupportNodeCount; ++I) {
        fprintf(Dst,
                "%d%s",
                SupportNodeIds[I],
                ((I + 1) % 16 == 0 || I + 1 == SupportNodeCount)
                    ? "\n"
                    : ",");
    }

    /* Linear-elastic PA 2200 material used by the static winglet test. */
    fprintf(Dst, "*MATERIAL, NAME=NYLON2200\n");
    fprintf(Dst,
            "*ELASTIC\n"
            "%.17g, %.17g\n",
            (double) YoungsModulusPa,
            (double) PoissonRatio);
    fprintf(Dst,
            "*DENSITY\n"
            "%.17g\n",
            (double) DensityKgM3);

    /* IRIT exports every C3D8 element into EL-ALL. */
    fprintf(Dst,
            "*SOLID SECTION, ELSET=EL-ALL, MATERIAL=NYLON2200\n");
    fprintf(Dst, ",\n");

    /* ANSYS External Model reads the first ABAQUS load step. */
    fprintf(Dst, "*STEP\n");

    /* Fully clamp all nodes belonging to the selected SUPPORT faces. */
    fprintf(Dst,
            "*BOUNDARY\n"
            "SUPPORT_NODES, 1, 3, 0.0\n");

    /* Apply a true surface-normal pressure to the exported PRESSURE surface.
     * ANSYS External Model supports ABAQUS *DSLOAD entries with load type P. */
    fprintf(Dst,
            "*DSLOAD\n"
            "PRESSURE, P, %.17g\n",
            (double) PressurePa);

    fprintf(Dst, "*END STEP\n");
    fprintf(Dst, "** END ADDED ANSYS WINGLET ANALYSIS DEFINITION\n");

    fclose(Dst);
    Dst = NULL;

#ifdef _WIN32
    if (!MoveFileExA(TempInp,
                     Opt->Config.InpFile,
                     MOVEFILE_REPLACE_EXISTING |
                     MOVEFILE_COPY_ALLOWED)) {
        fprintf(stderr,
                "ANSYS winglet INP: failed to replace original INP.\n");
        return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
    }
#else
    if (rename(TempInp, Opt->Config.InpFile) != 0) {
        fprintf(stderr,
                "ANSYS winglet INP: failed to replace original INP.\n");
        return AbortTemporaryInp(Src, Dst, SupportNodeIds, TempInp);
    }
#endif

    fprintf(stderr,
            "ANSYS winglet analysis definition added: "
            "support=body%d/%s, support_faces=%d, support_nodes=%d, "
            "PA2200 E=%.6g GPa, pressure=%.6g kPa.\n",
            OPT_SUPPORT_BODY_INDEX,
            SupportSideLabel,
            SupportFaceCount,
            SupportNodeCount,
            (double) (YoungsModulusPa * 1.0e-9),
            (double) (PressurePa * 1.0e-3));

    IritFree(SupportNodeIds);
    return TRUE;
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
static int FinishCalculixBatch(FILE *Src, FILE *Dst, FILE *Dat, FILE *Sed,
                               FILE *VolumeOut, FILE *DesignOut,
                               CagdRType *ElementEnergy,
                               CagdRType *ElementVolume,
                               CagdRType *NodeDisp,
                               int Result)
{
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
    return Result;
}

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

    /* Reuse the FEM topology parsed once from the first exported INP. */
    if (!Opt -> Lattice.FEMTopologyBound ||
        Opt -> FEM.NodeCount == 0 || Opt -> FEM.ElementCount == 0) {
        fprintf(stderr, "CalculiX: persistent FEM topology is unavailable.\n");
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
    }

    for (I = 0; I < Opt -> FEM.NodeCount; ++I) {
        OptNodeStruct *Node = &Opt -> FEM.Nodes[I];
        MinZ = IRIT_MIN(MinZ, Node -> Z);
        MaxZ = IRIT_MAX(MaxZ, Node -> Z);
        NodeMaxId = IRIT_MAX(NodeMaxId, Node -> Id);
    }
    for (I = 0; I < Opt -> FEM.ElementCount; ++I)
        ElementMaxId = IRIT_MAX(ElementMaxId, Opt -> FEM.Elements[I].Id);

    ZTolerance = IRIT_MAX(1.0, fabs(MaxZ - MinZ)) * 1.0e-9;
    for (I = 0; I < Opt -> FEM.NodeCount; ++I) {
        OptNodeStruct *Node = &Opt -> FEM.Nodes[I];
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
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
    }

    /* Rewrite node coordinates with bounded precision.  IRIT's full decimal
     * strings can exceed CalculiX's numeric input field width.  Connectivity
     * and IDs are copied unchanged. */
    snprintf(TempInp, sizeof(TempInp), "%s.ccx.tmp", Opt -> Config.InpFile);
    Src = fopen(Opt -> Config.InpFile, "r");
    Dst = fopen(TempInp, "w");
    if (Src == NULL || Dst == NULL)
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);

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
    for (I = 0; I < Opt -> FEM.NodeCount; ++I) {
        OptNodeStruct *Node = &Opt -> FEM.Nodes[I];
        if (fabs(Node -> Z - MinZ) <= ZTolerance) {
            fprintf(Dst, "%d%s", Node -> Id,
                    (++N % 16 == 0 || N == FixedCount) ? "\n" : ",");
        }
    }
    fprintf(Dst, "*NSET, NSET=LOADED\n");
    N = 0;
    for (I = 0; I < Opt -> FEM.NodeCount; ++I) {
        OptNodeStruct *Node = &Opt -> FEM.Nodes[I];
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
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
    }
#else
    if (rename(TempInp, Opt -> Config.InpFile) != 0)
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
#endif

    /* CALCULIX_CCX permits another installation without changing the code. */
    EnvExe = getenv("CALCULIX_CCX");
    snprintf(CcxExe, sizeof(CcxExe), "%s",
             EnvExe != NULL && *EnvExe != '\0' ? EnvExe :
             "C:\\irit\\CalculiX-2.23.0-win-x64\\bin\\ccx_MT.exe");
    if (!OptFileExists(CcxExe)) {
        fprintf(stderr, "CalculiX executable not found: %s\n", CcxExe);
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
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
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
    }

    snprintf(DatFile, sizeof(DatFile), "%s%c%s.dat",
             Opt -> Config.WorkDir, OPT_PATH_SEP, JobName);
    Dat = fopen(DatFile, "r");
    if (Dat == NULL) {
        fprintf(stderr, "CalculiX result file was not created: %s\n", DatFile);
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
    }

    ElementEnergy = (CagdRType *) OptAllocZero((size_t) ElementMaxId + 1u,
                                                sizeof(CagdRType));
    ElementVolume = (CagdRType *) OptAllocZero((size_t) ElementMaxId + 1u,
                                                sizeof(CagdRType));
    NodeDisp = (CagdRType *) OptAllocZero((size_t) NodeMaxId + 1u,
                                          sizeof(CagdRType));
    if (ElementEnergy == NULL || ElementVolume == NULL || NodeDisp == NULL)
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);

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

    if (DispCount < Opt -> FEM.NodeCount ||
        EnergyCount != Opt -> FEM.ElementCount ||
        VolumeCount != Opt -> FEM.ElementCount) {
        fprintf(stderr,
                "CalculiX output count mismatch: disp=%d/%d, energy=%d/%d, "
                "volume=%d/%d.\n",
                DispCount, Opt -> FEM.NodeCount,
                EnergyCount, Opt -> FEM.ElementCount,
                VolumeCount, Opt -> FEM.ElementCount);
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);
    }

    Sed = fopen(Opt -> Config.SedFile, "w");
    VolumeOut = fopen(Opt -> Config.VolumeFile, "w");
    DesignOut = fopen(Opt -> Config.DesignParametersFile, "w");
    if (Sed == NULL || VolumeOut == NULL || DesignOut == NULL)
        return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                                   ElementEnergy, ElementVolume, NodeDisp,
                                   FALSE);

    fprintf(Sed, "Element Number\tENERGYPOTENTIAL/VOLUME (mJ/mm^3)\n");
    fprintf(VolumeOut, "Element Number\tVOLUME (mm^3)\n");
    for (I = 0; I < Opt -> FEM.ElementCount; ++I) {
        int Id = Opt -> FEM.Elements[I].Id;
        if (ElementVolume[Id] <= 0.0)
            return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut,
                                       DesignOut, ElementEnergy, ElementVolume,
                                       NodeDisp, FALSE);
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
            Opt -> FEM.ElementCount, (double) (MaxDeflectionM * MToMm),
            (double) (TotalVolumeM3 * Opt -> Config.Density));
    Ok = TRUE;
    return FinishCalculixBatch(Src, Dst, Dat, Sed, VolumeOut, DesignOut,
                               ElementEnergy, ElementVolume, NodeDisp, Ok);
}
