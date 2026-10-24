/*****************************************************************************
* Point	surface	distance computation.					     *
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

#define PERTURB_EPS 5e-2
#define NUMERIC_TOL 1e-10

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

	    if (IRIT_PRSR_IS_SRF_OBJ(PObj)) {
		int i;
		CagdRType UMin, UMax, VMin, VMax, Tm;
	        IritCagdSrfStruct
		    *Srf = PObj -> U.Srfs;
		void
		    *SrfPtHandle = IritMvarDistSrfPointPrep(Srf);

		IritCagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

		fprintf(stderr, "Starting numeric solution, 200000 points...");
		Tm = IritMiscCPUTime(TRUE);

		/* Guess UV values, evaluate to E3 and call SrfPtDist to   */
		/* find back the UV and verify.				   */
		for (i = 0; i < 200000; i++) {
		    CagdUVType UV;
		    CagdPType PtE3;
		    IritMvarPtStruct *AllInvUV, *InvUV;
		    IritMvarPtStruct
			*PerturbUV = IritMvarPtNew(2);

		    UV[0] = IritMiscRandom(UMin, UMax);
		    UV[1] = IritMiscRandom(VMin, VMax);

		    IRIT_CAGD_SRF_EVAL_E3(Srf, UV[0], UV[1], PtE3);

		    /* Perturb the location by Eps. */
		    PerturbUV -> Pt[0] =
				 UV[0] + IritMiscRandom(-PERTURB_EPS,
							PERTURB_EPS);
		    PerturbUV -> Pt[1] =
				 UV[1] + IritMiscRandom(-PERTURB_EPS,
							PERTURB_EPS);
		    PerturbUV -> Pt[0] =
				   IRIT_BOUND(PerturbUV -> Pt[0], UMin, UMax);
		    PerturbUV -> Pt[1] = 
				   IRIT_BOUND(PerturbUV -> Pt[1], VMin, VMax);

		    AllInvUV = IritMvarLclDistSrfPoint(Srf, SrfPtHandle, PtE3,
						       &PerturbUV, 0.01,
						       NUMERIC_TOL);

		    if (AllInvUV == NULL)
			fprintf(stderr, "Failed to compute inverse for UV = (%f, %f)\n",
			UV[0], UV[1]);
		    else {
			CagdRType d,
			    MinDist = IRIT_INFNTY;

			for (InvUV = AllInvUV;
			     InvUV != NULL;
			     InvUV = InvUV -> Pnext) {
			    d = IRIT_PT2D_DIST_SQR(UV, AllInvUV -> Pt);
			    if (MinDist > d)
				MinDist = d;
			}

			/* Compare distance squares! */
			if (MinDist > IRIT_SQR(NUMERIC_TOL)) {
			    fprintf(stderr, "\nImprecise solution: UV = (%f, %f), error = %18.16g\n",
				    UV[0], UV[1], sqrt(MinDist));
			}
		    }

		    IritMvarPtFreeList(AllInvUV);
		}

		Tm = IritMiscCPUTime(FALSE);
		fprintf(stderr, "done (time = %f)\n", Tm);
		Tm = IritMiscCPUTime(TRUE);
		fprintf(stderr, "Starting full blown solution, 20000 points...");

		/* Guess UV values, evaluate to E3 and call SrfPtDist to   */
		/* find back the UV and verify.				   */
		for (i = 0; i < 20000; i++) {
		    CagdUVType UV;
		    CagdPType PtE3;
		    IritMvarPtStruct *AllInvUV, *InvUV;

		    UV[0] = IritMiscRandom(UMin, UMax);
		    UV[1] = IritMiscRandom(VMin, VMax);

		    IRIT_CAGD_SRF_EVAL_E3(Srf, UV[0], UV[1], PtE3);

		    AllInvUV = IritMvarLclDistSrfPoint(Srf, SrfPtHandle, PtE3,
						       NULL, 0.01, NUMERIC_TOL);

		    if (AllInvUV == NULL)
		        fprintf(stderr, "Failed to compute inverse for UV = (%f, %f)\n",
				UV[0], UV[1]);
		    else {
			for (InvUV = AllInvUV;
			     InvUV != NULL;
			     InvUV = InvUV -> Pnext) {
			    if (IRIT_PT_APX_EQ_E2_EPS(UV, InvUV -> Pt,
						      NUMERIC_TOL))
				break;
			}
			if (InvUV == NULL) {
		            fprintf(stderr, "Failed to find precise solution for UV = (%f, %f)\n",
				    UV[0], UV[1]);
			}
		    }
		}

		Tm = IritMiscCPUTime(FALSE);
		fprintf(stderr, "done (time = %f)\n", Tm);

		IritMvarDistSrfPointFree(SrfPtHandle);		
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
