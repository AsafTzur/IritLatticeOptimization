/*****************************************************************************
*   Constructs locally varying trivariate tiles	of minimal accumulated depth *
* in microstructure, using a call back function, toward	better CT approval.  *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, Dec 2017    *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
    IritTrivTVStruct *DefMap;
} UserMicroLocalDataStruct;

static int PreProcessTile1FaceParam(
				  IritUserMicroPreProcessTileCBStruct *CBData,
				  int X,
				  int Y,
				  int Z,
				  IritTrivTVBndryType Bndry,
				  IritUserMicroTileBndryPrmStruct *BPrm);
static IritPrsrObjectStruct *PreProcessTile(
				  IritPrsrObjectStruct *Tile,
				  IritUserMicroPreProcessTileCBStruct *CBData);

/*****************************************************************************
* DESCRIPTION:								     *
*   Prepare one	face parameters	for the	3D grid	tile synthesized on the	fly. *
*									     *
* PARAMETERS:								     *
*   CBData:	The call back data.					     *
*   X,Y, Z:	Indices	of the tile in deformation function.		     *
*   BPrm:	Parameters of this face	to update.			     *
*									     *
* RETURN VALUE:								     *
*   int:								     *
*****************************************************************************/
static int PreProcessTile1FaceParam(
				  IritUserMicroPreProcessTileCBStruct *CBData,
				  int X,
				  int Y,
				  int Z,
				  IritTrivTVBndryType Bndry,
				  IritUserMicroTileBndryPrmStruct *BPrm)
{
    int i, XYZ2,
        X2 = X,
        Y2 = Y,
        Z2 = Z,
        XYZ = ((X * 17) % CBData -> NumOfTiles[0] + Y * 157);
    IRIT_STATIC_DATA int
	// QueenPlacementIndcs[8] = { 5, 2, 4, 6, 0, 3, 1, 7 };
        QueenPlacementIndcs[8] = { 1, 4, 6, 0, 3, 7, 5, 2 };

    IRIT_ZAP_MEM(BPrm, sizeof(IritUserMicroTileBndryPrmStruct));

    BPrm -> OuterRadius = 0.04;
    BPrm -> BndryShape = 0.5;

    switch (Bndry) {
        case IRIT_TRIV_U_MIN_BNDRY:
	    X2 -= 1;
	    break;
        case IRIT_TRIV_U_MAX_BNDRY:
	    X2 += 1;
	    break;
        case IRIT_TRIV_V_MIN_BNDRY:
	    Y2 -= 1;
	    break;
        case IRIT_TRIV_V_MAX_BNDRY:
	    Y2 += 1;
	    break;
        case IRIT_TRIV_W_MIN_BNDRY:
	    Z2 -= 1;
	    break;
        case IRIT_TRIV_W_MAX_BNDRY:
	    Z2 += 1;
	    break;
        default:
	    assert(0);
	    return 0;
    }
    XYZ2 = ((X2 * 17) % CBData -> NumOfTiles[0] + Y2 * 157);

    /* Compute center location in [0,1]. */
    BPrm -> FaceCenter[0] = Bndry == IRIT_TRIV_U_MIN_BNDRY ||
		            Bndry == IRIT_TRIV_U_MAX_BNDRY ?
			    0.5 : (QueenPlacementIndcs[XYZ % 8] +
				   QueenPlacementIndcs[XYZ2 % 8]) / 14.0;
    BPrm -> FaceCenter[1] = Bndry == IRIT_TRIV_V_MIN_BNDRY ||
		            Bndry == IRIT_TRIV_V_MAX_BNDRY ?
			    0.5 : (QueenPlacementIndcs[(XYZ + 2) % 8] +
				   QueenPlacementIndcs[(XYZ2 + 2) % 8]) / 14.0;
    BPrm -> FaceCenter[2] = 0.5;


    /* Map center location to [-0.4 to 0.4]. */
    for (i = 0; i < 3; i++)
      BPrm -> FaceCenter[i] = BPrm -> FaceCenter[i] * 0.8 - 0.4;

    if (Bndry == IRIT_TRIV_W_MIN_BNDRY || Bndry == IRIT_TRIV_W_MAX_BNDRY) {
	/* Force no translation vertically. */
	BPrm -> FaceCenter[2] = -IRIT_INFNTY;
    }

    fprintf(stderr, "Face %d :: %d  %d  %d ::  %9.5f  %9.5f  %9.5f\n",
	    Bndry, X, Y, Z,
	    BPrm -> FaceCenter[0], BPrm -> FaceCenter[1], BPrm -> FaceCenter[2]);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Call back function from the	microstructure code to synthesize a tile.    *
*									     *
* PARAMETERS:								     *
*   Tile:   Tile to preprocess.	 Here we expect	a NULL as we build tile	     *
*	    from scratch.						     *
*   CBData: The	call back data.						     *
*									     *
* RETURN VALUE:								     *
*   IritPrsrObjectStruct *: processed tile.				     *
*****************************************************************************/
static IritPrsrObjectStruct *PreProcessTile(
				  IritPrsrObjectStruct *Tile,
				  IritUserMicroPreProcessTileCBStruct *CBData)
{
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;
    CagdRType
        *LclMinDmn = CBData -> TileLclDmnMin,
        *LclMaxDmn = CBData -> TileLclDmnMax;
    IritGeomBBBboxStruct BBox;
    IritUserMicroTileBndryPrmStruct UMinPrms, UMaxPrms,
			            VMinPrms, VMaxPrms, 
                                    WMinPrms, WMaxPrms;

    assert(Tile == NULL);                         /* We build tiles here... */

    fprintf(stderr, "Tile[%d,%d,%d] from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.3f)\n",
	    CBData -> TileIdxs[0],
	    CBData -> TileIdxs[1],
	    CBData -> TileIdxs[2],
	    LclMinDmn[0],
	    LclMinDmn[1],
	    LclMinDmn[2],
	    LclMaxDmn[0],
	    LclMaxDmn[1],
	    LclMaxDmn[2]);

    if (!PreProcessTile1FaceParam(CBData,
				  CBData -> TileIdxs[0],
				  CBData -> TileIdxs[1],
				  CBData -> TileIdxs[2],
				  IRIT_TRIV_U_MIN_BNDRY, &UMinPrms) ||
	!PreProcessTile1FaceParam(CBData,
				  CBData -> TileIdxs[0],
				  CBData -> TileIdxs[1],
				  CBData -> TileIdxs[2],
				  IRIT_TRIV_U_MAX_BNDRY, &UMaxPrms) ||
	!PreProcessTile1FaceParam(CBData,
				  CBData -> TileIdxs[0],
				  CBData -> TileIdxs[1],
				  CBData -> TileIdxs[2],
				  IRIT_TRIV_V_MIN_BNDRY, &VMinPrms) ||
	!PreProcessTile1FaceParam(CBData,
				  CBData -> TileIdxs[0],
				  CBData -> TileIdxs[1],
				  CBData -> TileIdxs[2],
				  IRIT_TRIV_V_MAX_BNDRY, &VMaxPrms) ||
	!PreProcessTile1FaceParam(CBData,
				  CBData -> TileIdxs[0],
				  CBData -> TileIdxs[1],
				  CBData -> TileIdxs[2],
				  IRIT_TRIV_W_MIN_BNDRY, &WMinPrms) ||
	!PreProcessTile1FaceParam(CBData,
				  CBData -> TileIdxs[0],
				  CBData -> TileIdxs[1],
				  CBData -> TileIdxs[2],
				  IRIT_TRIV_W_MAX_BNDRY, &WMaxPrms)) {
        return NULL;
    }

    Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms, &VMinPrms, &VMaxPrms,
				    &WMinPrms, &WMaxPrms, FALSE, NULL);
    assert(Tile != NULL);
    IritGeomBBComputeBboxObject(Tile, &BBox, FALSE);
    if (BBox.Min[0] < 0.0 || BBox.Max[0] > 1.0 ||
	BBox.Min[1] < 0.0 || BBox.Max[1] > 1.0 ||
	BBox.Min[2] < 0.0 || BBox.Max[2] > 1.0) {
        fprintf(stderr, "Warning: Tile spans beyond the unit box.\n");
    }

#define DEBUG_VERIFY_JACOBIAN
#ifdef DEBUG_VERIFY_JACOBIAN
    {
	IritTrivTVStruct *TV;

	for (TV = Tile -> U.Trivars; TV != NULL; TV = TV -> Pnext) {
	    IritMvarMVStruct
		*J = IritMvarCalculateTVJacobian(TV);
	    IritCagdBBoxStruct BBox;
	    
	    IritMvarMVBBox(J, &BBox);
	    IritMvarMVFree(J);
	    if (BBox.Min[0] * BBox.Max[0] < 0.0) {
		fprintf(stderr, "Negative Jacobian tile found\n");
		IritTrivDbg(TV);
	    }
	}
    }
#endif /* DEBUG_VERIFY_JACOBIAN */

    Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);

#define DEBUG_USER_MS_MAKE_TV_OBJS
#ifdef DEBUG_USER_MS_MAKE_TV_OBJS
    {
        IritPrsrObjectStruct *PTmp;
	IritTrivTVStruct *TV, *BTV,
	    *BzrTVs = NULL;	

	for (TV = Tile -> U.Trivars; TV != NULL; TV = TV -> Pnext) {
	    if (IRIT_TRIV_IS_BEZIER_TV(TV)) {
		BTV = IritTrivTVCopy(TV);
		IRIT_LIST_PUSH(BTV, BzrTVs);
	    }
	    else
		BzrTVs = IritCagdListAppend(IritTrivCnvrtBsp2BzrTV(TV),
					    BzrTVs);
	}
	PTmp = IritPrsrLnkListToListObject(BzrTVs, IRIT_PRSR_OBJ_TRIVAR);
	IritPrsrFreeObject(Tile);
	Tile = PTmp;
    }
#endif /* DEBUG_USER_MS_MAKE_TV_OBJS */

    return Tile;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Create a micro structure with a varying-in-size tiling example.	     *
*									     *
* PARAMETERS:								     *
*   None								     *
*									     *
* RETURN VALUE:								     *
*   void								     *
*****************************************************************************/
static void GenerateMicroStructure(void)
{
    const char
        *InputDefMap = "DefMap.itd";
    int i, Handler;    
    IritPrsrObjectStruct *MS, *DefMapPObj;
    IritMvarMVStruct *DeformMV;
    IritTrivTVStruct *TVMap;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    DefMapPObj = IritPrsrGetDataFiles(&InputDefMap, 1, FALSE, FALSE);
    if (DefMapPObj == NULL) {
        fprintf(stderr, "Failed to load the deformation function.\n");
        return;
    }
    assert(IRIT_PRSR_IS_TRIVAR_OBJ(DefMapPObj));
    TVMap = DefMapPObj -> U.Trivars;
    DeformMV = IritMvarCnvrtTVToMV(TVMap);

    /* Create the structure to be passed to the call back function. */
    LclData.DefMap = TVMap;

    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.TilingType = IRIT_USER_MICRO_TILE_REGULAR;
    MSParam.DeformMV = DeformMV;
    MSParam.ApproxLowOrder = 4;

    /* Sets boundary end conditions on the geometry - cap the tiles in all  */
    /* boundaries and as a side effect color trivar tiles on boundaries,    */
    /* so one can set boundary conditions (i.e. toward analysis).           */
    MSParam.ShellCapBits = 0;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;       /* Tile is synthesized on the fly. */
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i].TilesPerIntervals = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i].Len = 1;
    }
    MSRegularParam -> TilingSteps[0].TilesPerIntervals[0] = 8;
    MSRegularParam -> TilingSteps[1].TilesPerIntervals[0] = 8;
    MSRegularParam -> TilingSteps[2].TilesPerIntervals[0] = 8;


    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessTile;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */

    /* Generate the MS with shifted centers. */
    fprintf(stderr, "Generating MS with shifted centers...\n");
    
    MS = IritUserMicroStructComposition(&MSParam);/* Cnstrct microstructure.*/

    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    Handler = IritPrsrOpenDataFile("MSShiftedCenters.itd", FALSE, 1); 
    if (MS != NULL) {
        IritPrsrPutObjectToHandler(Handler, MS);
	IritPrsrFreeObject(MS);
    }
    IritPrsrCloseStream(Handler, TRUE);

    /** End **/
    IritMvarMVFree(DeformMV);
    IritUserMicroTileFree(MSRegularParam -> Tile);

    /* Free the structure for the call back function. */
    IritTrivTVFree(TVMap);
}

int main(int argc, char **argv)
{ 
    GenerateMicroStructure();
    
    return 0;
}
