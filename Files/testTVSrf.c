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

int main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(FALSE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct *PTV, *PSrf, *PTSrfs,
		*PObj = IritPrsrGetObjects(Handler);
	    IritTrimSrfStruct *AllTSrfs;

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    /* PObj should hold a list object with a TV and a Srf: */
	    PTV = IritPrsrListObjectGet(PObj, 0);
	    PSrf = IritPrsrListObjectGet(PObj, 1);

	    if (!IRIT_PRSR_IS_TRIVAR_OBJ(PTV) || !IRIT_PRSR_IS_SRF_OBJ(PSrf)) {
	        fprintf(stderr, "Expects a trivar and a surface\n");
		exit(1);
	    }

	    IritMdlBooleanSetTolerances(0.05, 1e-10, 0.025);

	    /* The key function - divide Srfs at all knot planes of TV: */
	    AllTSrfs = IritUserDivideSrfsAtAllTVInteriorKnot(PSrf -> U.Srfs,
							     PTV -> U.Trivars);

	    PTSrfs = IritPrsrGenTRIMSRFObject(AllTSrfs);
	    IritPrsrStdoutObject(PTSrfs, FALSE);
	    IritPrsrFreeObject(PTSrfs);
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
