/*****************************************************************************
*   Constructs locally varying tiles in	microstructure constructions using   *
* a call back function.							     *
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

typedef enum {
    MS_TILE_CROSS,
    MS_TILE_DIAGONAL,
    MS_TILE_CROSS_DIAGONAL
} MSTileType;

typedef struct IritMicroStrctCBStruct {
    IrtRType *Graded;
    IrtRType *Params;
    MSTileType TileType;  
} IritMicroStrctCBStruct;

static IritPrsrObjectStruct *MSGetTileAux(MSTileType TileType,
					  IrtRType *Params,
					  IrtRType *Graded,
					  char ** const ErrStr);
static IritPrsrObjectStruct *MSGradingTileCB(
				 IritPrsrObjectStruct *Tile,
				 IritUserMicroPreProcessTileCBStruct *CBData);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates and return the desired tile.                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   Tile:     Tile type to create.                                           *
*   Params:   For Cross tile: (Outer Radius, Inner Radius)                   *
*             For Diagonal tiles: (Center Size, Corner Size, Smooth Factor). *
*             For Cross-Diagonal tile: (Cross Radius, Diagonal Radius)       *
*   Graded:   We be used to create a graded lattice in the ruled direction   *
*             (variable arm thicknesses).				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IritPrsrObjectStruct *:    The created tile.                             *
*****************************************************************************/
static IritPrsrObjectStruct *MSGetTileAux(MSTileType TileType,
					  IrtRType *Params,
					  IrtRType *Graded,
					  char ** const ErrStr)
{
    int i;
    IrtRType
        GradedMid = IRIT_BLEND(Graded[0], Graded[1], 0.5);
    IritPrsrObjectStruct
        *TileObj = NULL;

    *ErrStr = NULL;

    switch (TileType) {
        default:
        case MS_TILE_CROSS:
	{
	    IrtRType OuterRadii[6], InnerRadii[6];

	    for (i = 0; i < 6; i++) {
	        IrtRType 
		    Gr = i == 0 ? Graded[0] 
			        : (i == 5 ? Graded[1] : GradedMid);
	        OuterRadii[i] = Params[0] * Gr;
		InnerRadii[i] = Params[1] * Gr;
	    }

	    TileObj = IritUserMicro3DCrossTile2(FALSE, OuterRadii,
						Params[1] == 0.0 ?
						    NULL :
						    InnerRadii,
						2 * Params[0], ErrStr);
	    break;
	}
	case MS_TILE_DIAGONAL:
	{
	    IrtRType CrnrSizes[8],
	        CrnrVertScl[2] = { 1.0, 1.0 };

	    for (i = 0; i < 8; i++)
	        CrnrSizes[i] = Params[1] * (i < 4 ? Graded[0] : Graded[1]);

	    TileObj = IritUserMicroDiagTile1(Params[0], 1.0, CrnrSizes,
					     CrnrVertScl, NULL, FALSE,
					     Params[2], NULL, ErrStr);
	    break;
	}
        case MS_TILE_CROSS_DIAGONAL:
	{
	    IrtRType CrossRadii[6], DiagRadii[8];

	    for (i = 0; i < 6; i++)
	        CrossRadii[i] = Params[0];
	    for (i = 0; i < 8; i++)
	      DiagRadii[i] = Params[1] * (i < 4 ? Graded[0] : Graded[1]);

	    TileObj = IritUserMicro3DCrossDiagTile(TRUE, CrossRadii,
						   DiagRadii, 2 * Params[0],
						   ErrStr);
	    break;
	}
    }

    return TileObj;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Given t, between zero and one, returns the required grading at that t.   *
*									     *
* PARAMETERS:								     *
*   CBdata:  An	optional call back data	that can be passed by the main	     *
*	     function to these call back functions.			     *
*   t:	     Parameter in the third dimension, in [0, 1].		     *
*									     *
* RETURN VALUE:								     *
*   IritPrsrObjectStruct:  Created graded tile.				     *
*****************************************************************************/
static IritPrsrObjectStruct *PreProcessTile(
				  IritPrsrObjectStruct *Tile,
				  IritUserMicroPreProcessTileCBStruct *CBData)
{
    char *ErrStr;
    IritMicroStrctCBStruct
	*LclData = (IritMicroStrctCBStruct *) CBData -> CBFuncData;
    CagdRType Graded[2],
	WMin = CBData -> DefMapDmnMin[2],
	Dw = CBData -> DefMapDmnMax[2] - WMin,
        *LclMinDmn = CBData -> TileLclDmnMin,
        *LclMaxDmn = CBData -> TileLclDmnMax;

    assert(Tile == NULL);                         /* We build tiles here... */
    
    Graded[0] = IRIT_BLEND(LclData -> Graded[1], LclData -> Graded[0],
			   WMin + Dw * LclMinDmn[2]);
    Graded[1] = IRIT_BLEND(LclData -> Graded[1], LclData -> Graded[0],
			   WMin + Dw * LclMaxDmn[2]);

    Tile = MSGetTileAux(LclData -> TileType, LclData -> Params,
			Graded, &ErrStr);

    if (Tile == NULL)
        Tile = IritPrsrGenNUMValObject(0.0);         /* Return something... */
    else
        Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);

    return Tile;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Create a micro structure with a varying-in-size tiling example.	     *
*									     *
* PARAMETERS:								     *
*   TileType:  Type of tile to use.					     *
*									     *
* RETURN VALUE:								     *
*   void								     *
*****************************************************************************/
static void GenerateMicroStructure(MSTileType TileType)
{
    int i, Handler;
    IrtRType TileParams[3],
	Graded[2] = { 0.5, 1.2 };
    IritPrsrObjectStruct *MS,
	*DefMap = NULL;
    IritMvarMVStruct *DeformMV;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRegularParamStruct *MSRegularParam;
    IritMicroStrctCBStruct LclData;

    if ((Handler = IritPrsrOpenDataFile("TVDef.itd", TRUE, TRUE)) >= 0) {
	DefMap = IritPrsrGetObjects(Handler);
	IritPrsrCloseStream(Handler, TRUE);
    }

    if (DefMap == NULL) {
	 fprintf(stderr, "Failed to load deformation map\n");
	 exit(1);
    }
    assert(IRIT_PRSR_IS_TRIVAR_OBJ(DefMap));
    DeformMV = IritMvarCnvrtTVToMV(DefMap -> U.Trivars);

    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.TilingType = IRIT_USER_MICRO_TILE_REGULAR;

    MSParam.DeformMV = DeformMV;
    MSParam.ApproxLowOrder = 4;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;/* Tile is synthesized by call back func. */
    MSRegularParam -> TilingStepMode = TRUE;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i].TilesPerIntervals = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i].Len = 1;
	MSRegularParam -> TilingSteps[i].TilesPerIntervals[0] = 1;
    }

    switch (TileType) {
	default:
	    assert(0);
	case MS_TILE_CROSS:
	    TileParams[0] = 0.2;
	    TileParams[1] = 0.0;
	    break;
	case MS_TILE_DIAGONAL:
	    TileParams[0] = 0.2;
	    TileParams[1] = 0.1;
	    TileParams[2] = 0.4;
	    break;
	case MS_TILE_CROSS_DIAGONAL:
	    TileParams[0] = 0.05;
	    TileParams[1] = 3.5;
	    break;
    }
    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessTile;
    LclData.Params = TileParams;
    LclData.Graded = Graded;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */

    MS = IritUserMicroStructComposition(&MSParam);/* Cnstrct microstructure.*/

    IritMvarMVFree(DeformMV);
    IritUserMicroTileFree(MSRegularParam -> Tile);
    for (i = 0; i < 3; ++i)
        IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    IritPrsrPutObjectToFile3("MSTest.itd", MS, 0); 

    IritPrsrFreeObject(MS);
    IritPrsrFreeObject(DefMap);
}

int main(int argc, char **argv)
{
    GenerateMicroStructure(MS_TILE_CROSS);
    
    return 0;
}
