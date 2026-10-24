/*****************************************************************************
*   Constructs locally varying tiles in	microstructure constructions using   *
* a call back function.	 Tile is loaded	from a file.			     *
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

typedef struct UserMicroLocalDataStruct {/* User specific data in CB funcs. */
    IritTrivTVStruct *DefMap;
    IritTrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
} UserMicroLocalDataStruct;

static IritPrsrObjectStruct *PreProcessTile(
				  IritPrsrObjectStruct *Tile,
				  IritUserMicroPreProcessTileCBStruct *CBData)
{
    char const *TileFile = "/temp/tile.itd";
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax;
    CagdPType DfDw;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    IritTrivTVDomain(LclData -> DefMap,
		     &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);
    IRIT_TRIV_TV_EVAL_E3(LclData -> DwDefMap,
		         CBData -> TileLclDmnMin[0],
		         CBData -> TileLclDmnMin[1],
		         CBData -> TileLclDmnMin[2], DfDw);

    fprintf(stderr, "Tile[%d,%d,%d] from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f), DfDw = (%.4f %.4f %.4f)\n",
	    CBData -> TileIdxs[0],
	    CBData -> TileIdxs[1],
	    CBData -> TileIdxs[2],
	    CBData -> TileLclDmnMin[0],
	    CBData -> TileLclDmnMin[1],
	    CBData -> TileLclDmnMin[2],
	    CBData -> TileLclDmnMax[0],
	    CBData -> TileLclDmnMax[1],
	    CBData -> TileLclDmnMax[2],
	    DfDw[0], DfDw[1], DfDw[2]);

    IritPrsrSetFlattenObjects(FALSE);
    Tile = IritPrsrGetDataFiles(&TileFile, 1, 1, 1);
    assert(Tile != NULL);
    Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);

#ifdef DEBUG_PRE_PROCESS_RETURN_NULL
    {
        static ii = 0;
	if (++ii % 3 == 0)
	    return NULL;
    }
#endif /* DEBUG_PRE_PROCESS_RETURN_NULL */

    return Tile;
}

static IritPrsrObjectStruct *PostProcessTile(IritPrsrObjectStruct *Tile,
				       IritUserMicroPostProcessTileCBStruct
				                                    *CBData)
{
    // IritPrsrStdoutObject(Tile, FALSE);

    // IritPrsrFreeObject(Tile);
    // return NULL;

#ifdef DEBUG_POST_PROCESS_RETURN_NULL
    {
        static ii = 0;
	if (++ii % 3 == 0)
	    return NULL;
    }
#endif /* DEBUG_POST_PROCESS_RETURN_NULL */

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
    char *ErrStr;
    int i, ErrLine, Handler;    
    IritPrsrObjectStruct *MS;
    IritMvarMVStruct *DeformMV;
    IritTrivTVStruct *TVMap;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    TVMap = IritTrivTVReadFromFile("/temp/DefMap.itd", &ErrStr, &ErrLine);
    DeformMV = IritMvarCnvrtTVToMV(TVMap);

    /* Create the structure to be passed to the call back function. */
    LclData.DefMap = TVMap;
    LclData.DuDefMap = IritTrivTVDerive(TVMap, IRIT_TRIV_CONST_U_DIR);
    LclData.DvDefMap = IritTrivTVDerive(TVMap, IRIT_TRIV_CONST_V_DIR);
    LclData.DwDefMap = IritTrivTVDerive(TVMap, IRIT_TRIV_CONST_W_DIR);

    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.TilingType = IRIT_USER_MICRO_TILE_REGULAR;

    /* Sets boundary end conditions on the geometry - cap the tiles in all  */
    /* boundaries and as a side effect color trivar tiles on boundaries,    */
    /* so one can set boundary conditions (i.e. toward analysis).           */
    /*
    MSParam.ShellCapBits = IRIT_USER_MICRO_BIT_CAP_UMIN |
                           IRIT_USER_MICRO_BIT_CAP_UMAX |
                           IRIT_USER_MICRO_BIT_CAP_VMIN |
                           IRIT_USER_MICRO_BIT_CAP_VMAX |
                           IRIT_USER_MICRO_BIT_CAP_WMIN |
			   IRIT_USER_MICRO_BIT_CAP_WMAX;
    */
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;
    MSParam.ApproxLowOrder = 4;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;       /* Tile is synthesized on the fly. */
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i].TilesPerIntervals = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i].Len = 1;
    }
    MSRegularParam -> TilingSteps[0].TilesPerIntervals[0] = 7;
    MSRegularParam -> TilingSteps[1].TilesPerIntervals[0] = 5;
    MSRegularParam -> TilingSteps[2].TilesPerIntervals[0] = 3;

    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessTile;
    MSRegularParam -> PostProcessCBFunc = PostProcessTile;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */
    
    MS = IritUserMicroStructComposition(&MSParam);/* Cnstrct microstructure.*/

    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    IritMvarMVFree(DeformMV);

    IritUserMicroTileFree(MSRegularParam -> Tile);

    Handler = IritPrsrOpenDataFile("/temp/MS.itd", FALSE, 1); 
    if (MS != NULL) {
        IritPrsrPutObjectToHandler(Handler, MS);
	IritPrsrFreeObject(MS);
    }
    IritPrsrCloseStream(Handler, TRUE);

    /* Free the structure for the call back function. */
    IritTrivTVFree(TVMap);
    IritTrivTVFree(LclData.DuDefMap);
    IritTrivTVFree(LclData.DvDefMap);
    IritTrivTVFree(LclData.DwDefMap);
}

int main(int argc, char **argv)
{ 
    GenerateMicroStructure();
    
    return 0;
}
