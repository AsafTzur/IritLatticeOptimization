/*****************************************************************************
*   Constructs locally varying surface tiles in	microstructure constructions *
* using	a call back function.						     *
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

/*****************************************************************************
* DESCRIPTION:								     *
*   Maps a parameter between zero and one to a radius of the cone to use.    *
*									     *
* PARAMETERS:								     *
*   Param:	     Parameter between zero and	one.			     *
*   MinRad, MaxRad:  Radii amplitudes to allow the sine	wave in	between.     *
*									     *
* RETURN VALUE:								     *
*   CagdRType *:    Radius of cone.					     *
*****************************************************************************/
static CagdRType MapParamToAmplitude(CagdRType Param,
				     CagdRType MinRad,
				     CagdRType MaxRad)
{
    /* Map [0, 1] in Param to a sine period but start and end at MinRad: */
    return (sin(Param * 2 * M_PI - M_PI / 2) + 1.0) * 0.5
						* (MaxRad - MinRad) + MinRad;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Constructs a cone parallel to the Z	axis from Z = 0	(radius	Rad1) to     *
* Z = 1	(and radius Rad2).  Cone will have no bases.			     *
*									     *
* PARAMETERS:								     *
*   Rad1, Rad2:	 Radii of two cone bases.				     *
*									     *
* RETURN VALUE:								     *
*   IritCagdSrfStruct *:  The constructed cone.				     *
*****************************************************************************/
static IritCagdSrfStruct *MakeZConeTile(CagdRType Rad1, CagdRType Rad2)
{
    static const IritCagdPtStruct
	Center = { NULL, NULL, { 0.0, 0.0, 0.0 } };
    IritCagdCrvStruct *TCrv,
        *Circ1 = IritCagdBspCrvCreatePCircle(&Center, Rad1),
        *Circ2 = IritCagdBspCrvCreatePCircle(&Center, Rad2);
    IritCagdSrfStruct *Srf;
    IrtHmgnMatType Mat;

    IritMiscMatGenMatTrans(0.0, 0.0, 1.0, Mat);
    TCrv = IritCagdCrvMatTransform(Circ2, Mat);
    IritCagdCrvFree(Circ2);
    Circ2 = TCrv;

    Srf = IritCagdRuledSrf(Circ1, Circ2, 2, 2);

    IritCagdCrvFree(Circ1);
    IritCagdCrvFree(Circ2);

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Constructs a tile as three open cones in [0, 1]^3 that locally fit the   *
* micro-structure in [UMin, UMax, VMin,	VMax, WMIn, WMax], and assuming	the  *
* entire trivariate mapping function domains' is [0, 1]^3.		     *
*   The	tiles are sized	so the entire structure	will create a sine wave	     *
* approximation.							     *
*									     *
* PARAMETERS:								     *
*   UMin, UMax,	VMin, VMax, WMin, WMax:	 Position of this tile in the	     *
*	deformation trivariate that is assumed to have domain [0, 1]^3.	     *
*   MinRad, MaxRad:  Radii amplitudes to allow the sine	wave in	between.     *
*									     *
* RETURN VALUE:								     *
*   IritPrsrObjectStruct *: Tile with three ruled surfaces spanning [0, 1]^3 *
*			 as this local tile for this UVW location.	     *
*****************************************************************************/
static IritPrsrObjectStruct *MakeSineCrossTile(CagdRType UMin,
					       CagdRType UMax,
					       CagdRType VMin,
					       CagdRType VMax,
					       CagdRType WMin,
					       CagdRType WMax,
					       CagdRType MinRad,
					       CagdRType MaxRad)
{
    IrtHmgnMatType Mat1, Mat2;
    IritCagdSrfStruct *TSrf,
        *UCone = MakeZConeTile(MapParamToAmplitude(UMin, MinRad, MaxRad),
			       MapParamToAmplitude(UMax, MinRad, MaxRad)),
        *VCone = MakeZConeTile(MapParamToAmplitude(VMin, MinRad, MaxRad),
			       MapParamToAmplitude(VMax, MinRad, MaxRad)),
        *WCone = MakeZConeTile(MapParamToAmplitude(WMin, MinRad, MaxRad),
			       MapParamToAmplitude(WMax, MinRad, MaxRad));
    IritPrsrObjectStruct *PObj;

    /* Transform the cones to their correct position and orientation. */
    IritMiscMatGenMatRotY1(M_PI / 2.0, Mat1);
    IritMiscMatGenMatTrans(0.0, 0.5, 0.5, Mat2);
    IritMiscMatMultTwo4by4(Mat1, Mat1, Mat2);
    TSrf = IritCagdSrfMatTransform(UCone, Mat1);
    IritCagdSrfFree(UCone);
    UCone = TSrf;

    IritMiscMatGenMatRotX1(-M_PI / 2.0, Mat1);
    IritMiscMatGenMatTrans(0.5, 0.0, 0.5, Mat2);
    IritMiscMatMultTwo4by4(Mat1, Mat1, Mat2);
    TSrf = IritCagdSrfMatTransform(VCone, Mat1);
    IritCagdSrfFree(VCone);
    VCone = TSrf;

    IritMiscMatGenMatTrans(0.5, 0.5, 0.0, Mat1);
    TSrf = IritCagdSrfMatTransform(WCone, Mat1);
    IritCagdSrfFree(WCone);
    WCone = TSrf;

    /* Build a tile from the three surfaces. */
    PObj = IritPrsrGenLISTObject(IritPrsrGenSRFObject(UCone));
    IritPrsrListObjectAppend(PObj, IritPrsrGenSRFObject(VCone));
    IritPrsrListObjectAppend(PObj, IritPrsrGenSRFObject(WCone));

    return PObj;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Called back	data just before the tile is mapped - we will substitute in  *
* our local tile and free the old one.					     *
*									     *
* PARAMETERS:								     *
*   Tile:   Canonical tile as provided by the user.			     *
*   CBData: Call back data as we provide in IritUserMicroStructComposition   *
* below.								     *
*									     *
* RETURN VALUE:								     *
*   IritPrsrObjectStruct *:	New tile to use	here.			     *
*****************************************************************************/
static IritPrsrObjectStruct *PreProcessSineTile(
				   IritPrsrObjectStruct *Tile,
				   IritUserMicroPreProcessTileCBStruct *CBData)
{
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax, Du, Dv, Dw;
    CagdPType DfDw;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    IritTrivTVDomain(LclData -> DefMap,
		     &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);
    IRIT_TRIV_TV_EVAL_E3(LclData -> DwDefMap,
		         CBData -> TileLclDmnMin[0],
		         CBData -> TileLclDmnMin[1],
		         CBData -> TileLclDmnMin[2], DfDw);
    Du = UMax - UMin;
    Dv = UMax - UMin;
    Dw = UMax - UMin;

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

    IritPrsrFreeObject(Tile);                       /* Free the old (dummy) tile. */
    /* and create new tile as part of 3 orthogonal sine tubes, in [0, 1]^3: */
    Tile = MakeSineCrossTile((CBData -> TileLclDmnMin[0] - UMin) / Du,
			     (CBData -> TileLclDmnMax[0] - UMin) / Du,
			     (CBData -> TileLclDmnMin[1] - VMin) / Dv,
			     (CBData -> TileLclDmnMax[1] - VMin) / Dv,
			     (CBData -> TileLclDmnMin[2] - WMin) / Dw,
			     (CBData -> TileLclDmnMax[2] - WMin) / Dw,
			     0.1, 0.3);
    /* Map new tile to its proper place in deformation function UVW space. */
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
static void GenerateSineMicroStructure(void)
{
    char *ErrStr;
    int ErrLine, Handler, i;    
    IritPrsrObjectStruct *MS;
    IritMvarMVStruct *DeformMV;
    IritTrivTVStruct *TVMap;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    TVMap = IritTrivTVReadFromFile("c:/temp/DefMap.itd", &ErrStr, &ErrLine);
    if (TVMap == NULL) {
        fprintf(stderr, "Failed to load the deformation function.\n");
        return;
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
    MSParam.ShellCapBits = IRIT_USER_MICRO_BIT_CAP_UMIN |
                           IRIT_USER_MICRO_BIT_CAP_UMAX |
                           IRIT_USER_MICRO_BIT_CAP_VMIN |
                           IRIT_USER_MICRO_BIT_CAP_VMAX |
                           IRIT_USER_MICRO_BIT_SHELL_WMIN |
                           IRIT_USER_MICRO_BIT_CAP_WMAX;
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;      /* Dummy tile - will not be used... */
    MSRegularParam -> TilingStepMode = TRUE;

    MSRegularParam -> TilingSteps[0].TilesPerIntervals = 
			      (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
    MSRegularParam -> TilingSteps[0].Len = 1;
    MSRegularParam -> TilingSteps[0].TilesPerIntervals[0] = 7;

    MSRegularParam -> TilingSteps[1].TilesPerIntervals = 
			      (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
    MSRegularParam -> TilingSteps[1].Len = 1;
    MSRegularParam -> TilingSteps[1].TilesPerIntervals[0] = 5;

    MSRegularParam -> TilingSteps[2].TilesPerIntervals = 
			      (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
    MSRegularParam -> TilingSteps[2].Len = 1;
    MSRegularParam -> TilingSteps[2].TilesPerIntervals[0] = 2;

    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessSineTile;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */
    
    MS = IritUserMicroStructComposition(&MSParam);/* Cnstrct microstructure.*/

    IritMvarMVFree(DeformMV);
    IritUserMicroTileFree(MSRegularParam -> Tile);
    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    Handler = IritPrsrOpenDataFile("c:/temp/MS2.itd", FALSE, 1); 
    IritPrsrPutObjectToHandler(Handler, MS);
    IritPrsrFreeObject(MS);
    IritPrsrCloseStream(Handler, TRUE);

    for (i = 0; i < 3; i++) 
	IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    /* Free the structure for the call back function. */
    IritTrivTVFree(TVMap);
    IritTrivTVFree(LclData.DuDefMap);
    IritTrivTVFree(LclData.DvDefMap);
    IritTrivTVFree(LclData.DwDefMap);
}

int main(int argc, char **argv)
{ 
    GenerateSineMicroStructure();
    
    return 0;
}
