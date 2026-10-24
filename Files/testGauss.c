/*****************************************************************************
* Some naive exhaustive	attempts to search for interesting surfaces...	     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/misc_lib.h"

#define MAX_COEF_VAL	3
#define MAX_DEGREE	2
#define NUM_EXPRESSIONS 5

const int ExprIdx[NUM_EXPRESSIONS] = { 1, 2, 3, 4, 6 };

int AdvanceVec(CagdRType *Vec)
{
    int i;

    for (i = 0; i < NUM_EXPRESSIONS; i++) {
        if (Vec[ExprIdx[i]] == MAX_COEF_VAL)
	    Vec[ExprIdx[i]] = -MAX_COEF_VAL;
	else {
	    Vec[ExprIdx[i]]++;
	    break;
	}
    }

    return i >= NUM_EXPRESSIONS;
}

int main(int argc, char **argv)
{
    int i,
	Count = 0;
    IritCagdSrfStruct
        *PwrSrf = IritCagdPwrSrfNew(MAX_DEGREE + 1, MAX_DEGREE + 1,
				    IRIT_CAGD_PT_E3_TYPE);
    CagdRType
        **Pts = PwrSrf -> Points;

    /* Reset the coefficients' values. */
    IRIT_ZAP_MEM(Pts[3], sizeof(CagdRType) * IRIT_SQR(MAX_DEGREE + 1));
    for (i = 0; i < NUM_EXPRESSIONS; i++)
        Pts[3][ExprIdx[i]] = -MAX_COEF_VAL;

    do { /* Z axis */
        IRIT_GEN_COPY(Pts[2], Pts[3],
		      sizeof(CagdRType) * IRIT_SQR(MAX_DEGREE + 1));

        do { /* Y axis */
	    IRIT_GEN_COPY(Pts[1], Pts[2],
			  sizeof(CagdRType) * IRIT_SQR(MAX_DEGREE + 1));

	    do { /* X axis */
	        /* Examine this surface. */
	        IritCagdBBoxStruct BBox;
	        IritCagdSrfStruct
		    *BzrSrf = IritCagdCnvrtPwr2BzrSrf(PwrSrf),
		    *Numer = IritSymbSrfGaussCurvature(BzrSrf, TRUE);
		  // *Numer = IritSymbSrfMeanNumer(BzrSrf);

		if (Count == -50547)
		    printf("err");

		IritCagdSrfBBox(Numer, &BBox);
		IritCagdSrfFree(Numer);

		if (0 & IRIT_APX_EQ(BBox.Min[0], 0.0) &&
		    IRIT_APX_EQ(BBox.Max[0], 0.0) &&
		    IritCagdSrfIsCoplanarCtlMesh(BzrSrf) > IRIT_EPS) {
		    printf("Found a developable surface (%d):\n", Count);
		    IritPrsrStdoutObject(IritPrsrGenSRFObject(BzrSrf), FALSE);
		    IritPrsrStdoutObject(IritPrsrGenSRFObject(PwrSrf), FALSE);
		}

		IritCagdSrfFree(BzrSrf);
		if (++Count % 1000000 == 0)
		    fprintf(stderr, "Computed %d samples\n", Count);
	    }
	    while (!AdvanceVec(Pts[1]));
	}
	while (!AdvanceVec(Pts[2]));
    }
    while (!AdvanceVec(Pts[3]));

    printf("Done\n");
    getchar();

    return 0;
}
