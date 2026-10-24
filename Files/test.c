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

void main(int argc, char **argv)
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
	        if (IRIT_PRSR_IS_POLY_OBJ(PTmp)) {
		    IritPrsrPolygonStruct *Pl;

		    for (Pl = PTmp -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
		        int n = IritPrsrVrtxListLen(Pl -> PVertex);
			IritPrsrVertexStruct *V;
			IrtPlnType Plane;

			if (n > 4) {
			    fprintf(stderr,
				    "More than 4 vertices (n = %d)\n", n);
#			    ifdef DEBUG
				IritPrsrVrtxDbg(Pl -> PVertex);
#			    endif /* DEBUG */
			}

			IRIT_PLANE_COPY(Plane, Pl -> Plane);
			if (IritPrsrUpdatePolyPlane(Pl)) {
			    if (IRIT_DOT_PROD(Plane, Pl -> Plane) < 0.0) {
				fprintf(stderr, "Poly normal is flipped\n");
#				ifdef DEBUG
				    IritPrsrVrtxDbg(Pl -> PVertex);
#				endif /* DEBUG */
			    }
			}
			else {
			    fprintf(stderr,
				    "Failed to compute plane for poly\n");
#			    ifdef DEBUG
				IritPrsrVrtxDbg(Pl -> PVertex);
#			    endif /* DEBUG */
			}

			/* Disable all normals. */
			V = Pl -> PVertex;
			do {
			    IRIT_PRSR_RST_NORMAL_VRTX(V);
			    V = V -> Pnext;
			}
			while (V != NULL && V != Pl -> PVertex);
		    }
		}
	    }

	    IritPrsrPutObjectToFile3("/temp/TestPolys.itd", PObj, 0);
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

    exit(0);
}
