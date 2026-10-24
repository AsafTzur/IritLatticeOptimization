/* Multi-body support used by the optimization driver.  This file is kept
 * separate because it is also convenient to regenerate it independently of
 * the main experimental driver. */

#define OPT_ATTR_TILE     "_opt_tile"
#define OPT_ATTR_SUBTILE  "_opt_subtile"

static int OptParseCounts(const char *Text, int **Values, int *Count)
{
    const char *P = Text;
    int *Result = NULL, N = 0, Capacity = 0;

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
                                                (size_t) NewCapacity * sizeof(int));
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

static int OptCountBodyTVs(const IritPrsrObjectStruct *Obj)
{
    int I, Count = 0;
    const IritTrivTVStruct *TV;

    for (; Obj != NULL; Obj = Obj->Pnext) {
        if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
            for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
                Count += OptCountBodyTVs(IritPrsrListObjectGet(Obj, I));
        }
        else if (IRIT_PRSR_IS_TRIVAR_OBJ(Obj)) {
            for (TV = Obj->U.Trivars; TV != NULL; TV = TV->Pnext)
                ++Count;
        }
    }
    return Count;
}

static int OptCopyBodyTVs(OptWorkflowStateStruct *Opt,
                          const IritPrsrObjectStruct *Obj,
                          int *BodyIndex,
                          int *TileOffset)
{
    int I;
    const IritTrivTVStruct *TV;

    for (; Obj != NULL; Obj = Obj->Pnext) {
        if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
            for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
                if (!OptCopyBodyTVs(Opt, IritPrsrListObjectGet(Obj, I),
                                    BodyIndex, TileOffset))
                    return FALSE;
        }
        else if (IRIT_PRSR_IS_TRIVAR_OBJ(Obj)) {
            for (TV = Obj->U.Trivars; TV != NULL; TV = TV->Pnext) {
                OptBodyStruct *Body = &Opt->Bodies[*BodyIndex];
                long long TileCount;

                Body->Map = IritTrivTVCopy(TV);
                if (Body->Map == NULL)
                    return FALSE;
                Body->Map->Pnext = NULL;
                snprintf(Body->Name, sizeof(Body->Name), "%s%s%d",
                         Obj->ObjName != NULL ? Obj->ObjName : "Body_",
                         Obj->ObjName != NULL ? "_" : "", *BodyIndex + 1);
                for (I = 0; I < 3; ++I)
                    Body->Counts[I] = Opt->Run.BodyCounts[I][*BodyIndex];
                TileCount = (long long) Body->Counts[0] * Body->Counts[1] *
                            Body->Counts[2];
                if (TileCount <= 0 || TileCount > INT_MAX - *TileOffset)
                    return FALSE;
                Body->TileOffset = *TileOffset;
                Body->TileCount = (int) TileCount;
                IritTrivTVDomain(Body->Map,
                                 &Body->Min[0], &Body->Max[0],
                                 &Body->Min[1], &Body->Max[1],
                                 &Body->Min[2], &Body->Max[2]);
                *TileOffset += Body->TileCount;
                ++*BodyIndex;
            }
        }
    }
    return TRUE;
}

static int OptLoadBodies(OptWorkflowStateStruct *Opt)
{
    const char *FileName = Opt->Run.MacroShapePath;
    IritPrsrObjectStruct *Objects;
    int B = 0, TileOffset = 0, D;

    IritPrsrSetFlattenObjects(FALSE);
    Objects = IritPrsrGetDataFiles(&FileName, 1, TRUE, FALSE);
    if (Objects == NULL)
        return FALSE;
    Opt->BodyCount = OptCountBodyTVs(Objects);
    if (Opt->BodyCount <= 0) {
        IritPrsrFreeObjectList(Objects);
        return FALSE;
    }
    for (D = 0; D < 3; ++D) {
        if (Opt->Run.BodyCountLengths[D] != Opt->BodyCount) {
            fprintf(stderr, "Lattice count list %d has %d values; macro shape has %d bodies.\n",
                    D, Opt->Run.BodyCountLengths[D], Opt->BodyCount);
            IritPrsrFreeObjectList(Objects);
            return FALSE;
        }
    }
    Opt->Bodies = (OptBodyStruct *) OptAllocZero((size_t) Opt->BodyCount,
                                                 sizeof(OptBodyStruct));
    if (Opt->Bodies == NULL ||
        !OptCopyBodyTVs(Opt, Objects, &B, &TileOffset)) {
        IritPrsrFreeObjectList(Objects);
        return FALSE;
    }
    IritPrsrFreeObjectList(Objects);
    Opt->TotalTiles = TileOffset;
    return B == Opt->BodyCount;
}

static void OptWriteCsvHeader(FILE *F, const OptWorkflowStateStruct *Opt)
{
    const OptBodyStruct *Body = &Opt->Bodies[0];
    fprintf(F, "lattice_type,%d\n", Opt->Run.LatticeType);
    fprintf(F, "lattice_counts,%d,%d,%d\n",
            Body->Counts[0], Body->Counts[1], Body->Counts[2]);
    fprintf(F, "mesh_sample,%d\n", Opt->Run.MeshSample);
    fprintf(F, "mesh_type,%d\n", Opt->Run.MeshType);
}

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
    return LatticeType == Opt->Run.LatticeType &&
           U == Opt->Bodies[0].Counts[0] && V == Opt->Bodies[0].Counts[1] &&
           W == Opt->Bodies[0].Counts[2] &&
           MeshSample == Opt->Run.MeshSample && MeshType == Opt->Run.MeshType;
}

static int OptBodyForTile(const OptWorkflowStateStruct *Opt, int Tile)
{
    int B;
    for (B = 0; B < Opt->BodyCount; ++B)
        if (Tile >= Opt->Bodies[B].TileOffset &&
            Tile < Opt->Bodies[B].TileOffset + Opt->Bodies[B].TileCount)
            return B;
    return -1;
}

static void OptTagTile(OptWorkflowStateStruct *Opt, IritPrsrObjectStruct *Obj,
                       int Tile, int *GeomSub)
{
    int I;
    (void) Opt;
    if (Obj == NULL)
        return;
    if (IRIT_PRSR_IS_OLST_OBJ(Obj)) {
        for (I = 0; I < IritPrsrListObjectLength(Obj); ++I)
            OptTagTile(Opt, IritPrsrListObjectGet(Obj, I), Tile, GeomSub);
    }
    else {
        IritMiscAttrSetObjectIntAttrib(Obj, OPT_ATTR_TILE, Tile);
        IritMiscAttrSetObjectIntAttrib(Obj, OPT_ATTR_SUBTILE, (*GeomSub)++);
    }
}

static int OptAppendMeshOwner(OptWorkflowStateStruct *Opt, int Tile, int SubTile)
{
    if (Opt->MeshOwnerCount == Opt->MeshOwnerCapacity) {
        int NewCapacity = Opt->MeshOwnerCapacity == 0 ? 256 :
                          Opt->MeshOwnerCapacity * 2;
        OptMeshOwnerStruct *NewOwners = (OptMeshOwnerStruct *)
            OptReallocCopy(Opt->MeshOwners,
                           (size_t) Opt->MeshOwnerCapacity * sizeof(*NewOwners),
                           (size_t) NewCapacity * sizeof(*NewOwners));
        if (NewOwners == NULL)
            return FALSE;
        Opt->MeshOwners = NewOwners;
        Opt->MeshOwnerCapacity = NewCapacity;
    }
    Opt->MeshOwners[Opt->MeshOwnerCount].Tile = Tile;
    Opt->MeshOwners[Opt->MeshOwnerCount].SubTile = SubTile;
    ++Opt->MeshOwnerCount;
    return TRUE;
}

static int OptCollectMesh(OptWorkflowStateStruct *Opt, IritPrsrObjectStruct *Obj,
                          IritPrsrObjectStruct *Out, int SaveOwners)
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
        (Opt->Run.LatticeType == 0 || Opt->Run.LatticeType == 1)) {
        /* Cross1 emits W+,W-,V+,V-,U+,U-,core. Cross2 emits
         * core,U-,U+,V-,V+,W-,W+. Map both to the design indices. */
        static const int Cross1LegacyOrder[7] = { 0, 1, 2, 3, 4, 5, 6 };
        /* Actual Cross 1 export order, verified from INP centroids, is:
         * core,U-,U+,V-,V+,W-,W+.  PreProcessTile design order is:
         * W+,W-,V+,V-,U+,U-,core. */
        static const int Cross1VerifiedOrder[7] = { 6, 5, 4, 3, 2, 1, 0 };
        static const int Cross2Order[7] = { 6, 4, 5, 2, 3, 0, 1 };
        const int *Order = Opt->Run.LatticeType == 0 ?
            (Opt->Config.Tile0OwnerMapping == 1 ?
                Cross1VerifiedOrder : Cross1LegacyOrder) : Cross2Order;

        if (Tile < 0 || Tile >= Opt->TotalTiles)
            return FALSE;
        for (I = 0, TV = Obj->U.Trivars; TV != NULL; TV = TV->Pnext, ++I) {
            if (I >= 7 || !OptAppendMeshOwner(Opt, Tile, Order[I]))
                return FALSE;
            IritPrsrListObjectAppend(
                Out, IritPrsrGenTRIVARObject(IritTrivTVCopy(TV)));
        }
        return I == 7;
    }

    if (SaveOwners && IRIT_PRSR_IS_TRIVAR_OBJ(Obj) &&
        Opt->Run.LatticeType == 2) {
        /* DiagTile1 returns three eight-TV objects in this order:
         * center joints, diagonal rods, corner joints. Each internal list
         * was built with IRIT_LIST_PUSH, so vertex order is 7 down to 0.
         * Design variables occupy rod indices 8..15 (LDef == 8). */
        int GroupBase;

        if (Tile < 0 || Tile >= Opt->TotalTiles ||
            SubTile < 0 || SubTile > 2)
            return FALSE;
        GroupBase = SubTile * 8;
        for (I = 0, TV = Obj->U.Trivars; TV != NULL; TV = TV->Pnext, ++I) {
            if (I >= 8 ||
                !OptAppendMeshOwner(Opt, Tile, GroupBase + (7 - I)))
                return FALSE;
            IritPrsrListObjectAppend(
                Out, IritPrsrGenTRIVARObject(IritTrivTVCopy(TV)));
        }
        return I == 8;
    }

    if (SaveOwners && (Tile < 0 || Tile >= Opt->TotalTiles ||
                       SubTile < 0 || SubTile >= Opt->L ||
                       !OptAppendMeshOwner(Opt, Tile, SubTile)))
        return FALSE;
    IritPrsrListObjectAppend(Out, IritPrsrCopyObject(NULL, Obj, TRUE));
    return TRUE;
}

static IritPrsrObjectStruct *PostProcessTileBodies(
                            IritPrsrObjectStruct *Tile,
                            IritUserMicroPostProcessTileCBStruct *CB)
{
    UserMicroLocalDataStruct *Local =
        (UserMicroLocalDataStruct *) CB->CBFuncData;
    OptBodyStruct *Body = Local->Body;
    int SubTile = 0;
    int TileIndex = Body->TileOffset + CB->TileIdxs[0] +
        Body->Counts[0] * (CB->TileIdxs[1] + Body->Counts[1] * CB->TileIdxs[2]);
    OptTagTile(Local->Opt, Tile, TileIndex, &SubTile);
    return Tile;
}

static int OptValidateMesh(OptWorkflowStateStruct *Opt)
{
    int I;
    if (Opt->MeshOwnerCount <= 0) {
        fprintf(stderr, "Generated mesh contains no owned geometry patches.\n");
        return FALSE;
    }
    for (I = 0; I < Opt->MeshOwnerCount; ++I)
        if (Opt->MeshOwners[I].Tile < 0 ||
            Opt->MeshOwners[I].Tile >= Opt->TotalTiles ||
            Opt->MeshOwners[I].SubTile < 0 ||
            Opt->MeshOwners[I].SubTile >= Opt->L)
            return FALSE;
    return TRUE;
}

static void OptMatchBodyBoundaries(OptWorkflowStateStruct *Opt,
                                   const CagdRType Local[][3],
                                   int LocalCount,
                                   int LDefOffset)
{
    /* Within-body interfaces were matched immediately before this call.
     * Cross-body topology is not implicit in separate trivariates, so leave
     * unrelated bodies independent instead of joining them by index. */
    (void) Opt;
    (void) Local;
    (void) LocalCount;
    (void) LDefOffset;
}
