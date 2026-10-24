/*****************************************************************************
* Compute the minimal distance to a C^1	closed surface and inside/out.	     *
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

#define SRF_PT_DIST_SUBDIV_TOL 1e-3
#define SRF_PT_DIST_NUMER_TOL  1e-8

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

	    if (IRIT_PRSR_IS_SRF_OBJ(PObj)) {
		int n = 0;
		CagdRType x, y, z,
		    Step = 0.3;
		IritCagdBBoxStruct BBox;
		IritCagdSrfStruct *NrmlSrf,
		    *Srf = PObj -> U.Srfs;
		void
		    *SrfPtDistCache = IritMvarDistSrfPointPrep(Srf);

		IritCagdSrfBBox(Srf, &BBox);

		Handler = IritPrsrOpenDataFile("/temp/inclusion.itd",
					       FALSE, TRUE);

		NrmlSrf = IritSymbSrfNormalSrf(Srf);

		for (x = BBox.Min[0] - 0.5; x < BBox.Max[0] + 0.5; x += Step) {
		    IritCagdPtStruct Pt, Pt2;

		    Pt.Pt[0] = x;
		    for (y = BBox.Min[1] - 0.5;
			 y < BBox.Max[1] + 0.5;
			 y += Step) {
		        Pt.Pt[1] = y;
			for (z = BBox.Min[1] - 0.5;
			     z < BBox.Max[1] + 0.5;
			     z += Step) {
			    CagdUVType UV;
			    CagdVType Nrml, Vec;
			    CagdPType Pos;

			    n++;
			    Pt.Pt[2] = z;

			    if (IritMvarDistSrfPoint(Srf, SrfPtDistCache,
						     Pt.Pt, NULL, TRUE,
						     SRF_PT_DIST_SUBDIV_TOL,
						     SRF_PT_DIST_NUMER_TOL,
						     NULL, UV,
						     IRIT_CAGD_PT_NONE) ==
								        NULL) {
				fprintf(stderr, "Failed to compute surface distance for \"[%g %g %g]\"\n",
				        Pt.Pt[0], Pt.Pt[1], Pt.Pt[2]);
				exit(1);
			    }

			    IRIT_CAGD_SRF_EVAL_E3(NrmlSrf, UV[0], UV[1], Nrml);
			    IRIT_CAGD_SRF_EVAL_E3(Srf, UV[0], UV[1], Pos);
			    IRIT_VEC_SUB(Vec, Pt.Pt, Pos);

			    if (IRIT_DOT_PROD(Vec, Nrml) > 0.0) {
			        /* Dump the point in green. */
			        PTmp = IritPrsrGenPTObject(&Pt.Pt[0],
						     &Pt.Pt[1], &Pt.Pt[2]);
				IritMiscAttrSetObjectRGBColor(PTmp, 0, 255, 0);
				IritPrsrPutObjectToHandler(Handler, PTmp);
				IritPrsrFreeObject(PTmp);
			    }
			    else {
				IritCagdCrvStruct *Ln;
				
				/* Dump the point in red. */
				PTmp = IritPrsrGenPTObject(&Pt.Pt[0],
						     &Pt.Pt[1], &Pt.Pt[2]);
				IritMiscAttrSetObjectRGBColor(PTmp, 255, 0, 0);
				IritPrsrPutObjectToHandler(Handler, PTmp);
				IritPrsrFreeObject(PTmp);

				/* Dump a line to closest pt on surface. */
				IRIT_PT_COPY(Pt2.Pt, Pos);
				Ln = IritCagdMergePtPt(&Pt, &Pt2);
				PTmp = IritPrsrGenCRVObject(Ln);
				IritMiscAttrSetObjectRGBColor(PTmp, 0, 0, 255);
				IritPrsrPutObjectToHandler(Handler, PTmp);
				IritPrsrFreeObject(PTmp);
			    }
			}
		    }
		}

		IritPrsrCloseStream(Handler, TRUE);
		IritMvarDistSrfPointFree(SrfPtDistCache);
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
