/*****************************************************************************
*   Constructs bifurcation tiles.					     *
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
    int i, ErrLine;    
    IritMvarMVStruct *DeformMV;
    IritTrivTVStruct *TVMap;
    IritUserMicroParamStruct MSParam;
    IritUserMicroRandomBiFurcationParamStruct *MSRandomBFParams;

    if ((TVMap = IritTrivTVReadFromFile("/temp/DefMap.itd", 
				    &ErrStr, &ErrLine)) == NULL) {
	fprintf(stderr, "Failed to read def. map TV.\n");
	exit(1);
    }
    DeformMV = IritMvarCnvrtTVToMV(TVMap);

    /* Create the structure to be passed to the call back function. */
    IRIT_ZAP_MEM(&MSParam, sizeof(IritUserMicroParamStruct));
    MSParam.TilingType = IRIT_USER_MICRO_TILE_BIFURCATION_RANDOM;
    MSParam.ShellThickness = -1.0;
    MSParam.ShellCapBits = IRIT_USER_MICRO_BIT_CAP_ALL;
    MSParam.DeformMV = DeformMV;

    MSRandomBFParams = &MSParam.U.RandomBiFurcationParam;

    for (i = 0; i < 3; ++i) {
        MSRandomBFParams -> Orders[i] = 3;
        MSRandomBFParams -> NumCPInTile[i] = 8;
    }
    MSRandomBFParams -> SubdThreshold = 0.2;
    MSRandomBFParams -> RandomFactor  = 0.0;

    IritTrivTVSetDomain(TVMap, 0, 1, 0,1, 0, 1);

#define SAMPLE_POINTS_IN_MS
#ifdef SAMPLE_POINTS_IN_MS
{
    CagdRType x, y, z,
	Step = 0.05;
    IritUserMicro2TilingStructPtr Tiling;
    CagdPType UVWPt;

    Tiling = IritUserMicroRandomBifurcationTiling(
					    MSParam.DeformMV, 
					    MSRandomBFParams -> Orders,
					    MSRandomBFParams -> NumCPInTile,
					    MSRandomBFParams -> SubdThreshold,
					    MSRandomBFParams -> RandomFactor,
					    MSParam.ShellCapBits,
					    MSParam.ShellThickness);


    /* Sample points in the tiling. */
    for (z = Step * 0.5; z < 1.0; z += Step) {
        UVWPt[2] = z;

	for (y = Step * 0.5; y < 1.0; y += Step) {
	    UVWPt[1] = y;

	    for (x = Step * 0.5; x < 1.0; x += Step) {
	        CagdRType ResValue;
		IritPrsrObjectStruct *PPt;

		UVWPt[0] = x;

		if (!IritUserMicroFunctionalEvaluateUV(Tiling, UVWPt,
						       &ResValue)) {
		    fprintf(stderr, "Failed to evaluate TV.\n");
		    exit(1);
		}

		PPt = IritPrsrGenPTObject(&UVWPt[0], &UVWPt[1], &UVWPt[2]);
		if (ResValue > 0) {
		    CagdRType
			J = IritTrivTVEvalJacobian(
						TVMap,
						UVWPt[0], UVWPt[1], UVWPt[2]);

		    IritMiscAttrIDSetObjectRealAttrib(
					  PPt,
					  IRIT_ATTR_CREATE_ID(Jacobian), J);
		    IritMiscAttrIDSetObjectRGBColor(PPt, 0, 255, 0);
		}
		else
		    IritMiscAttrIDSetObjectRGBColor(PPt, 255, 0, 0);

		if (ResValue > 0) /* Del this 'if' to also get outside pts. */
		    IritPrsrStdoutObject(PPt, FALSE);
		IritPrsrFreeObject(PPt);
	    }
	}
    }

    IritUserMicroFunctionalFreeTiling(Tiling);
}
#else
{
    IritPrsrObjectStruct
        *MS = IritUserMicroStructComposition(&MSParam);/* Cnstr. microstructure. */
    int Handler = IritPrsrOpenDataFile("/temp/MS.itd", FALSE, 1); 

    if (MS != NULL) {
        IritPrsrPutObjectToHandler(Handler, MS);
	IritPrsrFreeObject(MS);
    }
    IritPrsrCloseStream(Handler, TRUE);
}
#endif /* SAMPLE_POINTS_IN_MS */

    /* Free the structure for the call back function. */
    IritTrivTVFree(TVMap);
}

int main(int argc, char **argv)
{ 
    GenerateMicroStructure();
    
    return 0;
}
