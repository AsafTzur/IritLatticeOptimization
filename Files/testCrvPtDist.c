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

int main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(TRUE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct *PTmp,
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    if (IRIT_PRSR_IS_CRV_OBJ(PObj)) {
		int n = 0;
		CagdRType x, y,
		    Step = 0.1;
		IritCagdBBoxStruct BBox;
		IritCagdCrvStruct
		    *Crv = PObj -> U.Crvs;
		void
		    *CrvPtDistCache = IritSymbDistCrvPointPrep(Crv);

		IritCagdCrvBBox(Crv, &BBox);

		Handler = IritPrsrOpenDataFile("/temp/inclusion.itd",
					       FALSE, TRUE);

		for (x = BBox.Min[0] - 0.5; x < BBox.Max[0] + 0.5; x += Step) {
		    IritCagdPtStruct Pt;

		    Pt.Pt[0] = x;
		    Pt.Pt[2] = 0.0;
		    for (y = BBox.Min[1] - 0.5;
			 y < BBox.Max[1] + 0.5;
			 y += Step) {
			n++;
			Pt.Pt[1] = y;
			if (IritSymbCrvPointInclusion(Crv, Pt.Pt,
						      IRIT_EPS) == 1) {
			    /* Dump the point in green. */
			    PTmp = IritPrsrGenPTObject(&Pt.Pt[0],
						       &Pt.Pt[1], &Pt.Pt[2]);
			    IritMiscAttrSetObjectRGBColor(PTmp, 0, 255, 0);
			    IritPrsrPutObjectToHandler(Handler, PTmp);
			    IritPrsrFreeObject(PTmp);
			}
			else {
			    CagdRType t;
			    IritCagdPtStruct PtE2;
			    IritCagdCrvStruct *Ln;

			    /* Dump the point in red. */
			    PTmp = IritPrsrGenPTObject(&Pt.Pt[0],
						       &Pt.Pt[1], &Pt.Pt[2]);
			    IritMiscAttrSetObjectRGBColor(PTmp, 255, 0, 0);
			    IritPrsrPutObjectToHandler(Handler, PTmp);
			    IritPrsrFreeObject(PTmp);

			    /* Dump line to closest point on crv (distance). */
			    t = IritSymbDistCrvPoint(Crv, CrvPtDistCache,
						     Pt.Pt, TRUE, IRIT_EPS,
						     IRIT_CAGD_PT_NONE);
			    PtE2.Pt[2] = 0.0;
			    IRIT_CAGD_CRV_EVAL_E2(Crv, t, PtE2.Pt);
			    Ln = IritCagdMergePtPt(&Pt, &PtE2);
			    PTmp = IritPrsrGenCRVObject(Ln);
			    IritMiscAttrSetObjectRGBColor(PTmp, 0, 0, 255);
			    IritPrsrPutObjectToHandler(Handler, PTmp);
			    IritPrsrFreeObject(PTmp);
			}
		    }
		}

		IritPrsrCloseStream(Handler, TRUE);
		IritSymbDistCrvPointFree(CrvPtDistCache);
		fprintf(stderr, "%d points were tested...\n", n);
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: test geom.itd\n");
	exit(2);
    }

    return 0;
}
