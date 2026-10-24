/*****************************************************************************
*   Constructs locally varying trivariate tiles	in FFD constructions	     *
* using a call back function.						     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, Feb 2025    *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

typedef CagdRType (*ConeThicknessCBFuncType)(CagdRType u,
					     CagdRType v,
					     CagdRType w);

typedef struct ConeInfoStruct {
    IritTrivTVStruct *TV;
    char *MSFileName;
    CagdRType X, Y, Z;  /* Number of tiles in each axis. */
    CagdRType SmoothFactor;
    IritTrivTVDirType VarDir;
} ConeInfoStruct;

typedef struct AuxInfoCBStruct {
    int VarDir;
    CagdRType  SmoothFactor;
} AuxInfoCBStruct;

static CagdRType ThicknessVaryingUVW1(CagdRType u, CagdRType v, CagdRType w);
static CagdRType ThicknessVaryingUVW2(CagdRType u, CagdRType v, CagdRType w);
static CagdRType ThicknessVaryingUVW3(CagdRType u, CagdRType v, CagdRType w);
static IritPrsrObjectStruct *PreProcessTileCB(IritPrsrObjectStruct *Tile,
					      IritTrivFFDInfoStruct *CBInfo);
static void GenerateMicroStructures(const ConeInfoStruct *ConeInfo);

/*****************************************************************************
* DESCRIPTION:								     *
*   A u varying tiler.							     *
*									     *
* PARAMETERS:								     *
*   u, v, w:	UVW location to	estimate thickness at, in the def. map.	     *
*									     *
* RETURN VALUE:								     *
*   CagdRType:	the desired thickness.					     *
*****************************************************************************/
static CagdRType ThicknessVaryingUVW1(CagdRType u, CagdRType v, CagdRType w)
{
    return 0.05 + u * 0.25;                  /* Maps [0, 1] to [0.05, 0.3] */
}

/*****************************************************************************
* DESCRIPTION:								     *
*   A v varying tiler.							     *
*									     *
* PARAMETERS:								     *
*   u, v, w:	UVW location to	estimate thickness at, in the def. map.	     *
*									     *
* RETURN VALUE:								     *
*   CagdRType:	the desired thickness.					     *
*****************************************************************************/
static CagdRType ThicknessVaryingUVW2(CagdRType u, CagdRType v, CagdRType w)
{
    return 0.01 + v * 0.29;		     /* Maps [0, 1] to [0.01, 0.3] */
}

/*****************************************************************************
* DESCRIPTION:								     *
*   A w varying tiler.							     *
*									     *
* PARAMETERS:								     *
*   u, v, w:	UVW location to	estimate thickness at, in the def. map.	     *
*									     *
* RETURN VALUE:								     *
*   CagdRType:	the desired thickness.					     *
*****************************************************************************/
static CagdRType ThicknessVaryingUVW3(CagdRType u, CagdRType v, CagdRType w)
{
    return 0.3 - w * 0.29;                   /* Maps [0, 1] to [0.3, 0.01] */
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Call back function to synthesize the tiles.				     *
*									     *
* PARAMETERS:								     *
*   Tile:   Tile to pre-process.  Here we expect a NULL as we build tile     *
*	    from scratch.						     *
*   FFDInfo: The call back info.					     *
*									     *
* RETURN VALUE:								     *
*   IritPrsrObjectStruct *:	Created tile.				     *
*****************************************************************************/
static IritPrsrObjectStruct *PreProcessTileCB(IritPrsrObjectStruct *Tile,
					      IritTrivFFDInfoStruct *FFDInfo)
{
    int DirIdx;
    AuxInfoCBStruct
        *AuxInfo = (AuxInfoCBStruct *) FFDInfo -> _AuxInfo;
    CagdRType CenterSize, CrnrSizes[8], CrnrVertScl[2], *Min, *Max;
    ConeThicknessCBFuncType VarCallFunc;

    assert(Tile == NULL);

    /* Using diagonal tile. */
    switch (AuxInfo -> VarDir) {
        default:
	    assert(0);
        case IRIT_TRIV_CONST_U_DIR:
	    VarCallFunc = ThicknessVaryingUVW1;
	    DirIdx = 0;
	    break;
        case IRIT_TRIV_CONST_V_DIR:
	    VarCallFunc = ThicknessVaryingUVW2;
	    DirIdx = 1;
	    break;
        case IRIT_TRIV_CONST_W_DIR:
	    VarCallFunc = ThicknessVaryingUVW3;
	    DirIdx = 2;
	    break;
    }

    Min = FFDInfo -> _DmnMin;
    Max = FFDInfo -> _DmnMax;

    /* for face based tiles - like a cross.
    FaceSizes[0] = VarCallFunc(FFDInfo -> _DmnMin[DirIdx]);
    FaceSizes[1] = VarCallFunc(FFDInfo -> _DmnMin[DirIdx]);
    FaceSizes[2] = VarCallFunc(FFDInfo -> _DmnMin[DirIdx]);
    FaceSizes[3] = VarCallFunc(FFDInfo -> _DmnMin[DirIdx]);
    FaceSizes[4] = VarCallFunc(FFDInfo -> _DmnMax[DirIdx]);
    CrnrSizes[5] = VarCallFunc(FFDInfo -> _DmnMax[DirIdx]);
    */

    CrnrSizes[0] = VarCallFunc(Min[0], Min[1], Min[2]);
    CrnrSizes[1] = VarCallFunc(Max[0], Min[1], Min[2]);
    CrnrSizes[2] = VarCallFunc(Min[0], Max[1], Min[2]);
    CrnrSizes[3] = VarCallFunc(Max[0], Max[1], Min[2]);
    CrnrSizes[4] = VarCallFunc(Min[0], Min[1], Max[2]);
    CrnrSizes[5] = VarCallFunc(Max[0], Min[1], Max[2]);
    CrnrSizes[6] = VarCallFunc(Min[0], Max[1], Max[2]);
    CrnrSizes[7] = VarCallFunc(Max[0], Max[1], Max[2]);


    CenterSize = VarCallFunc((FFDInfo -> _DmnMin[0] + 
			      FFDInfo -> _DmnMax[0]) * 0.5,
			     (FFDInfo -> _DmnMin[1] + 
			      FFDInfo -> _DmnMax[1]) * 0.5,
			     (FFDInfo -> _DmnMin[2] + 
			      FFDInfo -> _DmnMax[2]) * 0.5);

    CrnrVertScl[0] = CrnrVertScl[1] = 0.5;

    Tile = IritUserMicroDiagTile1(CenterSize, 1.0, CrnrSizes, CrnrVertScl, 
				  NULL, FALSE, AuxInfo -> SmoothFactor,
				  NULL, NULL);
    Tile = IritGeomTransformObjectInPlace(Tile, FFDInfo -> _SclMat);

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
static void GenerateMicroStructures(const ConeInfoStruct *ConeInfo)
{
    char MSFileName[IRIT_LINE_LEN];
    IritPrsrObjectStruct *MS, *Tile;
    AuxInfoCBStruct AuxInfo;
    const IritTrivTVStruct *DefMap;
    IritTrivFFDInfoStruct FFDInfo;

    sprintf(MSFileName, "%s_%d%d%d.itd", ConeInfo -> MSFileName,
	    (int) ConeInfo -> X, (int) ConeInfo -> Y, (int) ConeInfo -> Z);

    IRIT_ZAP_MEM(&FFDInfo, sizeof(IritTrivFFDInfoStruct));
    FFDInfo.FFDCB = PreProcessTileCB;

    AuxInfo.VarDir = ConeInfo -> VarDir;
    AuxInfo.SmoothFactor = ConeInfo -> SmoothFactor;
    FFDInfo._AuxInfo = &AuxInfo;

    DefMap = ConeInfo -> TV;
    if (DefMap == NULL) {
        fprintf(stderr, "Failed to get the deformation function.\n");
        return;
    }

    if (0) {
        CagdRType
	    CrnrSizes[8] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 },
	    CrnrVertScl[2] = { 0.5, 0.5 };

	Tile = IritUserMicroDiagTile1(0.2, 1.0, CrnrSizes, CrnrVertScl,
				      NULL, FALSE, 0.0, NULL, NULL);
    }

    MS = IritTrivFFDTileObjectInTV(NULL, DefMap,
				   ConeInfo -> X, ConeInfo -> Y, ConeInfo -> Z,
				   0, 0.0, 0.0, NULL,
				   &FFDInfo);

    IritPrsrPutObjectToFile3(MSFileName, MS, 0);

    IritPrsrFreeObject(MS);

}

int main(int argc, char **argv)
{
    int i;
    CagdVType
        Center = { 0.0, 0.0, 0.0 };
    IritTrivTVStruct
        *ConeTV = IritTrivPrimCone2(Center, 1.0, 0.5, 1.0, FALSE);
    const ConeInfoStruct ConesInfo[] =
    {
        { ConeTV,  "testFFDCB1", 2, 2, 4, 0.0,  IRIT_TRIV_CONST_W_DIR },
        { ConeTV,  "testFFDCB2", 2, 2, 6, 0.1,  IRIT_TRIV_CONST_U_DIR },
        { ConeTV,  "testFFDCB3", 3, 3, 8, 0.15, IRIT_TRIV_CONST_V_DIR },

        { NULL, NULL, 0, 0, 0, IRIT_TRIV_CONST_W_DIR }
    };

    for (i = 0; ConesInfo[i].TV != NULL; i++) {
        GenerateMicroStructures(&ConesInfo[i]);
    }

    IritTrivTVFree(ConeTV);
    
    return 0;
}
