/*****************************************************************************
*   Constructs locally varying tiles in	microstructure constructions using   *
* a call back function.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, Dec 2017    *
*****************************************************************************/

/* DefMap.itd created using:

Disk = boolone( pcircle( vector( 0, 0, 0 ), 1 ) );

TV = tFromSrfs( list( Disk * tz( -3 ),
		      Disk * tz( -2.3 ),
		      Disk * sc( 3 ) * tz( -0.8 ),
	              Disk * sc( 3 ) * tz( 0.8 ),
		      Disk * tz( 2.3 ),
		      Disk * tz( 3 ) ),
	        3, kv_open ) * sz( 1.5 );

view( TV, 1 );

save( "FatTV", TV );

*/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

static IrtRType
    GlblBase = 0.015,
    GlblAmpl = 0.1,
    GlblFreq = 1.0;

typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
    IritTrivTVStruct *DefMap;
    IritTrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
} UserMicroLocalDataStruct;


static IritPrsrObjectStruct *PreProcessTile(
				  IritPrsrObjectStruct *Tile,
				  IritUserMicroPreProcessTileCBStruct *CBData)
{
    IritUserMicroTileBndryPrmStruct
	UMinPrms, UMaxPrms, VMinPrms, VMaxPrms, WMinPrms, WMaxPrms;
    CagdRType t, Dw, WMid;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    Dw = CBData -> DefMapDmnMax[2] - CBData -> DefMapDmnMin[2];

    IRIT_ZAP_MEM(&UMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&UMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));

    UMinPrms.Circular = UMaxPrms.Circular = TRUE;
    VMinPrms.Circular = VMaxPrms.Circular = TRUE;
    WMinPrms.Circular = WMaxPrms.Circular = TRUE;

    /* Here is the function that controls the radii of a local tile. */
    WMid = (CBData -> TileLclDmnMin[2] + CBData -> TileLclDmnMax[2]) * 0.5;

    t = CBData -> DefMapDmnMin[2] + Dw * WMid;
    UMinPrms.OuterRadius =
        UMaxPrms.OuterRadius =
            VMinPrms.OuterRadius =
                VMaxPrms.OuterRadius =
		    GlblAmpl * IRIT_SQR(cos(GlblFreq * t * M_PI)) + GlblBase;

    t = CBData -> DefMapDmnMin[2] + Dw * CBData -> TileLclDmnMin[2];
    WMinPrms.OuterRadius =
                    GlblAmpl * IRIT_SQR(cos(GlblFreq * t * M_PI)) + GlblBase;

    t = CBData -> DefMapDmnMin[2] + Dw * CBData -> TileLclDmnMax[2];
    WMaxPrms.OuterRadius =
		    GlblAmpl * IRIT_SQR(cos(GlblFreq * t * M_PI)) + GlblBase;

    Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				    &VMinPrms, &VMaxPrms,
				    &WMinPrms, &WMaxPrms, FALSE, NULL);
    assert(Tile != NULL);
    Tile = IritGeomTransformObjectInPlace(Tile, CBData -> Mat);

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
*   OutputName:	 File name to save result in.				     *
*									     *
* RETURN VALUE:								     *
*   void								     *
*****************************************************************************/
static void GenerateMicroStructure(const char *OutputName)
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

    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.TilingType = IRIT_USER_MICRO_TILE_REGULAR;

    /* Sets boundary end conditions on the geometry - cap the tiles in five */
    /* boundaries and do a solid shell base at the (WMin) bottom.	    */
    /* No capping/shelling for now.
    MSParam.ShellCapBits = USER_MICRO_BIT_CAP_UMIN |
                           USER_MICRO_BIT_CAP_UMAX |
                           USER_MICRO_BIT_CAP_VMIN |
                           USER_MICRO_BIT_CAP_VMAX |
                           USER_MICRO_BIT_SHELL_WMIN |
                           USER_MICRO_BIT_CAP_WMAX;
    */
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;
    MSParam.ApproxLowOrder = 4;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;/* Tile is synthesized by call back func. */
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i].TilesPerIntervals = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i].Len = 1;
	MSRegularParam -> TilingSteps[i].TilesPerIntervals[0] = 4;
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
    fprintf(stderr, "%d tiles created\n", IritPrsrListObjectLength(MS));

    IritMvarMVFree(DeformMV);
    IritUserMicroTileFree(MSRegularParam -> Tile);
    for (i = 0; i < 3; ++i)
        IritFree(MSRegularParam -> TilingSteps[i].TilesPerIntervals);

    Handler = IritPrsrOpenDataFile(OutputName, FALSE, 1); 
    IritPrsrPutObjectToHandler(Handler, MS);
    IritPrsrFreeObject(MS);
    IritPrsrCloseStream(Handler, TRUE);

    /* Free the structure for the call back function. */
    IritTrivTVFree(TVMap);
}

int main(int argc, char **argv)
{
    if (argc != 5 ||
        sscanf(argv[1], "%lf", &GlblAmpl) != 1 ||
	sscanf(argv[2], "%lf", &GlblBase) != 1 ||
	sscanf(argv[3], "%lf", &GlblFreq) != 1) {
        fprintf(stderr, "Usage:  testMS6 Amplitude Base Freq Output.itd\nExample: testMS6 0.1 0.015 1 MS.itd");
        exit(1);
    }

    fprintf(stderr, "Running with: Amplitude %f, Base %f, Freq %f, into \"%s\"\n",
	    GlblAmpl, GlblBase, GlblFreq, argv[4]);


    GenerateMicroStructure(argv[4]);
    
    return 0;
}
