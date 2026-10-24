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
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

void main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(FALSE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct *Obj1, *Obj2, *Obj3, *Obj4,
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    /* PObj should hold a list object : */
	    Obj1 = IritPrsrListObjectGet(PObj, 0);
	    Obj2 = IritPrsrListObjectGet(PObj, 1);
	    Obj3 = IritPrsrListObjectGet(PObj, 2);
	    Obj4 = IritPrsrListObjectGet(PObj, 3);

	    /* Computing matching UV curves on the two given surfaces:      */
	    /* MvarMakeCrvsOnSrfsMatchSpeed is obsolete -		    */
	    /*             use IritTrimMatch2ndCrvLenSpeedAs1stCrv instead. */
	    MvarMakeCrvsOnSrfsMatchSpeed(Obj1 -> U.Srfs, Obj2 -> U.Srfs,
					 &Obj3 -> U.Crvs, &Obj4 -> U.Crvs);

	    /* Verify the result. */
	    assert(Obj3 -> U.Crvs -> Length == Obj4 -> U.Crvs -> Length);
	    if (Obj3 -> U.Crvs -> Order == 2 && Obj4 -> U.Crvs -> Order == 2) {
		int i;

	        for (i = 0; i < Obj3 -> U.Crvs -> Length; i++) {
		    CagdRType *R;
		    CagdPType Pt1, Pt2;

		    IRIT_CAGD_SRF_EVAL_E3(Obj1 -> U.Srfs,
				          Obj3 -> U.Crvs -> Points[1][i],
				          Obj3 -> U.Crvs -> Points[2][i], Pt1);
		    IRIT_CAGD_SRF_EVAL_E3(Obj2 -> U.Srfs,
				          Obj4 -> U.Crvs -> Points[1][i],
				          Obj4 -> U.Crvs -> Points[2][i], Pt2);
		    assert(IRIT_PT_APX_EQ_EPS(Pt1, Pt2, IRIT_EPS * 0.1));
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

    exit(0);
}
