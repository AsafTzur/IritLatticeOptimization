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
	    if (IRIT_PRSR_IS_TRIVAR_OBJ(PObj)) {
		int IDir;
		CagdRType t;
		IritTrivTVDirType Dir;
		IritMvarMVStruct *MV;

		Dir = IRIT_TRIV_NO_DIR;
		t = IRIT_INFNTY;
		if (IritTrivTVKnotHasC0Discont(PObj -> U.Trivars, &Dir, &t)) {
		    fprintf(stderr, "TV C0 Discont:  %d  %f: %smesh discont.\n",
			    Dir - IRIT_TRIV_CONST_U_DIR, t,
			    IritTrivTVIsMeshC0DiscontAt(PObj -> U.Trivars,
							Dir, t)
			    ? "" : "no ");

		}
		else {
		    Dir = IRIT_TRIV_NO_DIR;
		    t = IRIT_INFNTY;
		    if (IritTrivTVKnotHasC1Discont(PObj -> U.Trivars,
						   &Dir, &t)) {
			fprintf(stderr, "TV C1 Discont:  %d  %f: %smesh discont.\n",
				Dir - IRIT_TRIV_CONST_U_DIR, t,
				IritTrivTVIsMeshC1DiscontAt(PObj -> U.Trivars,
							    Dir, t)
				? "" : "no ");
		    }
		}

		MV = IritMvarCnvrtTVToMV(PObj -> U.Trivars);
		IDir = -1;
		t = IRIT_INFNTY;
		if (IritMvarMVKnotHasC0Discont(MV, &IDir, &t)) {
		    fprintf(stderr, "MV C0 Discont:  %d  %f: %smesh discont.\n",
			    IDir, t,
			    IritMvarMVIsMeshC0DiscontAt(MV, IDir, t)
			    ? "" : "no ");

		}
		else {
		    IDir = -1;
		    t = IRIT_INFNTY;
		    if (IritMvarMVKnotHasC1Discont(MV, &IDir, &t)) {
			fprintf(stderr, "MV C1 Discont:  %d  %f: %smesh discont.\n",
				IDir, t,
				IritMvarMVIsMeshC1DiscontAt(MV, IDir, t)
				? "" : "no ");
		    }
		}

	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: test tv.itd\n");
	exit(2);
    }

    exit(0);
}

int main4(int argc, char **argv)
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
	    if (IRIT_PRSR_IS_CRV_OBJ(PObj)) {
	        int n;
		CagdRType Pt[3];
		IritCagdCrvStruct
		    *Crv = PObj -> U.Crvs;
		void
		    *PrepPtCrv = IritSymbDistCrvPointPrep(Crv);

		IritMiscCPUTime(TRUE);
 
		Pt[2] = 0.0;
		for (n = 0; n < 100000; n++) {
		    Pt[0] = IritMiscRandom(-2, 2);
		    Pt[1] = IritMiscRandom(-2, 2);
		    IritSymbDistCrvPoint(Crv, PrepPtCrv, Pt,
					 TRUE, IRIT_EPS, IRIT_CAGD_PT_NONE);
		}

		IritSymbDistCrvPointFree(PrepPtCrv);
		fprintf(stderr, "Time = %lf\n", IritMiscCPUTime(FALSE));
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: test crv.itd\n");
	exit(2);
    }

    return 0;
}

int main3(int argc, char **argv)
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
	    if (IRIT_PRSR_IS_SRF_OBJ(PObj)) {
	        int n;
		CagdRType UMin, UMax, VMin, VMax, u, v, k1, k2;
		CagdVType D1, D2;
		IritCagdSrfStruct
		    *Srf = PObj -> U.Srfs;
		
		IritCagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
		IritSymbEvalSrfCurvPrep(Srf, TRUE);
		for (n = 0; n < 1000000; n++) {
		    u = UMin + IritMiscRandom(0, 1) * (UMax - UMin);
		    v = VMin + IritMiscRandom(0, 1) * (VMax - VMin);

		    IritSymbEvalSrfCurvature(Srf, u, v, TRUE,
					     &k1, &k2, D1, D2);
		}
		IritSymbEvalSrfCurvPrep(Srf, FALSE);
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
		    //    *E3Crv = IritMiscAttrGetObjectObjAttrib(PTmp, "E3Crv");
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
