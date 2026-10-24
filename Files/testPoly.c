/*****************************************************************************
* Filter the data.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"

int main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(TRUE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    /* Process the geometry. */
	    if (IRIT_PRSR_IS_POLY_OBJ(PObj)) {
	        IritPrsrPolygonStruct
		    *Pl = IritGeomMergeClosedLoopHoles(PObj -> U.Pl,
						       PObj -> U.Pl -> Pnext);

		IritPrsrStderrObject(IritPrsrGenPOLYObject(Pl));
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: PolyArea geom.itd\n");
	exit(2);
    }

    return 0;
}

int main2(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(TRUE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct *PTmp,
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    /* Filter the geometry. */
	    for (PTmp = PObj; PTmp != NULL; PTmp = PTmp -> Pnext) {
	        if (IRIT_PRSR_IS_TRIMSRF_OBJ(PTmp)) {
		    //IritPrsrObjectStruct
		    //   *E3Crv = IritMiscAttrGetObjectObjAttrib(PTmp, "E3Crv");
		    //
		    //if (E3Crv != NULL)
			IritPrsrStdoutObject(PTmp, FALSE);
		}
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: PolyArea geom.itd\n");
	exit(2);
    }

    return 0;
}

int main1(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
	IritPrsrSetFlattenObjects(TRUE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    if (IRIT_PRSR_IS_CRV_OBJ(PObj)) {
		IritCagdCrvStruct
		    *MCrvs = IritTrimMergeTrimmingCurves2Loops2(PObj -> U.Crvs,
							    IRIT_EPS);
		IritCagdDbg(MCrvs);
		IritCagdCrvFreeList(MCrvs);
	    }
	}
    }
    return 0;
}
