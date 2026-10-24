/*****************************************************************************
* Microstructures test example - direct	MS building from a given tile.	     *
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

/* Create a tile from a list of surfaces.		*/
IritUserMicroTileStruct *UserMicroCreateTileFrmSrfs(IritCagdSrfStruct *SrfList)
{
    IritPrsrObjectStruct *IPObj,
	*IPListObj = IritPrsrGenLISTObject(NULL);
    IritUserMicroTileStruct *Tile;
    IritCagdSrfStruct *Srf, *SrfNext;

    for (Srf = SrfList; Srf != NULL; Srf = SrfNext) {
	IPObj = IritPrsrGenSRFObject(Srf);
	IritPrsrListObjectAppend(IPListObj, IPObj);
	SrfNext = Srf -> Pnext;
	Srf -> Pnext = NULL;
    }

    Tile = IritUserMicroTileNew(IPListObj);
    return Tile;
}

/* Create a tile from a list of trimmed surfaces.	*/
IritUserMicroTileStruct *UserMicroCreateTileFrmTrimSrfs(
					        IritTrimSrfStruct *TrmSrfList)
{
    IritPrsrObjectStruct *IPObj,
	*IPListObj = IritPrsrGenLISTObject(NULL);
    IritUserMicroTileStruct *Tile;
    IritTrimSrfStruct *TrmSrf, *TrmSrfNext;

    for (TrmSrf = TrmSrfList; TrmSrf != NULL; TrmSrf = TrmSrfNext) {
	IPObj = IritPrsrGenTRIMSRFObject(TrmSrf);
	IritPrsrListObjectAppend(IPListObj, IPObj);
	TrmSrfNext = TrmSrf -> Pnext;
	TrmSrf -> Pnext = NULL;
    }

    Tile = IritUserMicroTileNew(IPListObj);

    return Tile;
}

/* Generate a list of two surfaces. */
IritCagdSrfStruct *GenSrfList(void) 
{
    IritCagdSrfStruct *Srf, *Srf2;
    CagdVType Center;

    Center[0] = Center[1] = 0.5;
    Center[2] = 0;

    Srf = IritCagdPrimCylinderSrf(Center, 0.5, 1, FALSE,
				  IRIT_CAGD_PRIM_CAPS_NONE);
    Srf2 = IritCagdPrimCylinderSrf(Center, 0.25, 1, FALSE,
				   IRIT_CAGD_PRIM_CAPS_NONE);
    Srf -> Pnext = Srf2;

    return Srf;
}

/* Generate a list of two trimmed-surfaces.	*/
IritTrimSrfStruct *GenTrimSrfList(void)
{
    IritTrimSrfStruct *TrmSrf, *TrmSrf2;
    IritTrimCrvStruct *TrCrv;
    IritTrimCrvSegStruct *TrCrvSeg;
    IritCagdSrfStruct *Srf;
    CagdVType Center;
    IritCagdPtStruct CircCenter;
    IritCagdCrvStruct *TCrv;

    Center[0] = Center[1] = 0.5;
    Center[2] = 0;

    IRIT_PT_COPY(CircCenter.Pt, Center);

    Srf = IritCagdPrimCylinderSrf(Center, 0.5, 1, FALSE,
				  IRIT_CAGD_PRIM_CAPS_NONE);
    TCrv = IritCagdBspCrvCreatePCircle(&CircCenter, 0.45);    
    TrCrvSeg = IritTrimCrvSegNew(TCrv, NULL);
    TrCrv = IritTrimCrvNew(TrCrvSeg);
    TrmSrf = IritTrimSrfNew(Srf, TrCrv, FALSE);


    Srf = IritCagdPrimCylinderSrf(Center, 0.25, 1, FALSE,
				  IRIT_CAGD_PRIM_CAPS_NONE);
    TCrv = IritCagdBspCrvCreatePCircle(&CircCenter, 0.25);    
    TrCrvSeg = IritTrimCrvSegNew(TCrv, NULL);
    TrCrv = IritTrimCrvNew(TrCrvSeg);
    TrmSrf2 = IritTrimSrfNew(Srf, TrCrv, FALSE);

    TrmSrf -> Pnext = TrmSrf2;

    return TrmSrf;
}

/* A tiling example.			*/
void TestFunction(CagdBType TrimFlag)
{
    char *ErrStr;
    int i, ErrLine, Handler;    
    IritCagdSrfStruct *SrfList;
    IritTrimSrfStruct *TrimSrfList;
    IritPrsrObjectStruct *MS;
    IritUserMicroTileStruct *Tile;
    IritMvarMVStruct *DeformMV;
    IritTrivTVStruct *TVMap1;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRegularParamStruct *MSRegularParam;

#ifdef DEBUG_TEST_DYN_MEM
    IritMiscDynMemoryDbgCheckMark(&Start, &KeepStackStep, &TrackAllocID);
#endif /* DEBUG_TEST_DYN_MEM */

    IritPrsrSetFlattenObjects(FALSE);

    TVMap1 = IritTrivTVReadFromFile("/temp/DefMap.itd", &ErrStr, &ErrLine);
    if (TVMap1 == NULL) {
        fprintf(stderr, "Failed to load the deformation function\n");
	exit(1);
    }

    if (TrimFlag == 0) {
	/*	Surface list example.		*/
	SrfList = GenSrfList();
	Tile = UserMicroCreateTileFrmSrfs(SrfList);
    }
    else if (TrimFlag == 1){
	/*	Trimmed surface list example.    */
	TrimSrfList = GenTrimSrfList();
	Tile = UserMicroCreateTileFrmTrimSrfs(TrimSrfList);
    }
    else {
        /*      Load a tile.                     */
        const char
	    *TileFileName = "/temp/ModelTile.itd";
        IritPrsrObjectStruct
	    *PObjTile = IritPrsrGetDataFiles(&TileFileName, 1, TRUE, TRUE);

	if (PObjTile == NULL) {
	    fprintf(stderr, "Failed to load the tile\n");
	    exit(1);
	}
	Tile = IritUserMicroTileNew(PObjTile);
    }

    DeformMV = IritMvarCnvrtTVToMV(TVMap1);
    IritTrivTVFree(TVMap1);

    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.DeformMV = DeformMV;
    MSParam.TilingType = IRIT_USER_MICRO_TILE_REGULAR;
    MSParam.ShellThickness = 0;
    MSParam.ShellCapBits = 0;
    MSParam.ApproxLowOrder = 0;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = Tile;
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> C0DiscontWScale = 0;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i].TilesPerIntervals = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i].Len = 1;
	MSRegularParam -> TilingSteps[i].TilesPerIntervals[0] = i + 2;
    }

    MS = IritUserMicroStructComposition(&MSParam);

    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    Handler = IritPrsrOpenDataFile("/temp/MS.itd", FALSE, 1); 
    IritPrsrPutObjectToHandler(Handler, MS);
    IritPrsrFreeObject(MS);
    IritUserMicroTileFree(Tile);
    IritPrsrCloseStream(Handler, TRUE);

#ifdef DEBUG_TEST_DYN_MEM
    Start = 0;
    IritMiscDynMemoryDbgCheckMark(&Start, &KeepStackStep, &TrackAllocID);
#endif /* DEBUG_TEST_DYN_MEM */
}

int main(int argc, char **argv)
{
    TestFunction(2);
    
    return 0;
}
