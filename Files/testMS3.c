/*****************************************************************************
*   Constructs locally varying tiles in	microstructure constructions using   *
* a call back function.	Tile is	loaded from a file.			     *
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
    IritTrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
} UserMicroLocalDataStruct;


static IritPrsrObjectStruct *PreProcessTile(
				  IritPrsrObjectStruct *Tile,
				  IritUserMicroPreProcessTileCBStruct *CBData)
{
    static IritPrsrObjectStruct
	*ReadTile = NULL;
    const char *tileFile = "tile.itd";
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

    fprintf(stderr, "Tile[%d,%d,%d] locally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n\t   globally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n\t   DfDw = (%.4f %.4f %.4f)\n",
	    CBData -> TileIdxs[0],
	    CBData -> TileIdxs[1],
	    CBData -> TileIdxs[2],
	    CBData -> TileLclDmnMin[0],
	    CBData -> TileLclDmnMin[1],
	    CBData -> TileLclDmnMin[2],
	    CBData -> TileLclDmnMax[0],
	    CBData -> TileLclDmnMax[1],
	    CBData -> TileLclDmnMax[2],
	    CBData -> DefMapDmnMin[0],
	    CBData -> DefMapDmnMin[1],
	    CBData -> DefMapDmnMin[2],
	    CBData -> DefMapDmnMax[0],
	    CBData -> DefMapDmnMax[1],
	    CBData -> DefMapDmnMax[2],
	    DfDw[0], DfDw[1], DfDw[2]);

    IritPrsrFreeObject(Tile);                  /* Free the old (dummy) tile. */
    
    IritPrsrSetFlattenObjects(FALSE);
    if (ReadTile == NULL) {
	if ((ReadTile = IritPrsrGetDataFiles(&tileFile, 1, 1, 1)) == NULL) {
	    fprintf(stderr, "Failed to load the tile\n");
	    exit(1);
	}
    }
    Tile = IritGeomTransformObject(ReadTile, CBData -> Mat);

    return Tile;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Converts a hierarchy of geometries into a linear list of geometries.     *
*									     *
*									     *
* PARAMETERS:								     *
*   MS:	    Input hierarchy.						     *
*   LinMS:  The	resulting equivalent linear geometry.			     *
*									     *
* RETURN VALUE:								     *
*   void								     *
*****************************************************************************/
static void FlattenHierarchy2LinearList(const IritPrsrObjectStruct *MS,
					IritPrsrObjectStruct *LinMS)
{
    int i;

    switch (MS -> ObjType) {
	case IRIT_PRSR_OBJ_LIST_OBJ:
	    for (i = 0; i < IritPrsrListObjectLength(MS); i++)
		FlattenHierarchy2LinearList(IritPrsrListObjectGet(MS, i),
					    LinMS);
	    break;
	default:
	    IritPrsrListObjectAppend(LinMS, IritPrsrCopyObject(NULL, MS,
							       TRUE));
	    break;
    }
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
    IritPrsrObjectStruct *MS, *LinMS;
    IritMvarMVStruct *DeformMV;
    IritTrivTVStruct *TVMap;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    if ((TVMap = IritTrivTVReadFromFile("DefMap.itd",
				     &ErrStr, &ErrLine)) == NULL) {
	 fprintf(stderr, "Failed to load deformation map\n");
	 exit(1);
    }
    DeformMV = IritMvarCnvrtTVToMV(TVMap);

    /* Create the structure to be passed to the call back function. */
    LclData.DefMap = TVMap;
    LclData.DuDefMap = IritTrivTVDerive(TVMap, IRIT_TRIV_CONST_U_DIR);
    LclData.DvDefMap = IritTrivTVDerive(TVMap, IRIT_TRIV_CONST_V_DIR);
    LclData.DwDefMap = IritTrivTVDerive(TVMap, IRIT_TRIV_CONST_W_DIR);

    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.TilingType = IRIT_USER_MICRO_TILE_REGULAR;

    /* Sets boundary end conditions on the geometry - cap the tiles in five */
    /* boundaries and do a solid shell base at the (WMin) bottom.	    */
    /* No capping/shelling for now.
    MSParam.ShellCapBits = IRIT_USER_MICRO_BIT_CAP_UMIN |
                           IRIT_USER_MICRO_BIT_CAP_UMAX |
                           IRIT_USER_MICRO_BIT_CAP_VMIN |
                           IRIT_USER_MICRO_BIT_CAP_VMAX |
                           IRIT_USER_MICRO_BIT_SHELL_WMIN |
                           IRIT_USER_MICRO_BIT_CAP_WMAX;
    */
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;/* Tile is synthesized by call back func. */
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i].TilesPerIntervals = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i].Len = 1;
	MSRegularParam -> TilingSteps[i].TilesPerIntervals[0] = 2;
    }

    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessTile;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */
    
    MS = IritUserMicroStructComposition(&MSParam);/* Cnstrct microstructure.*/

    LinMS = IritPrsrGenLISTObject(NULL);/* Create an empty list to populate.*/
    FlattenHierarchy2LinearList(MS, LinMS);  /* Convert into a linear list. */
    IritPrsrFreeObject(MS);
    MS = LinMS;
    fprintf(stderr, "%d surfaces created\n", IritPrsrListObjectLength(MS));

    IritMvarMVFree(DeformMV);
    IritUserMicroTileFree(MSRegularParam -> Tile);
    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    Handler = IritPrsrOpenDataFile("MS.itd", FALSE, 1); 
    IritPrsrPutObjectToHandler(Handler, MS);
    IritPrsrFreeObject(MS);
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
