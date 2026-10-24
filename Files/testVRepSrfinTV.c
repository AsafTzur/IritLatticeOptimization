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
	    IritPrsrObjectStruct
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    if (IRIT_PRSR_IS_VMODEL_OBJ(PObj)) {
		const IritVMdlVModelStruct
		    *VMdl = PObj -> U.VMdls;
		IritMdlModelStruct *Mdl,
		    *Mdls = IritVMdlCnvrtVMdls2Mdls(NULL, VMdl);
		IritMdlTrimSrfStruct *TSrf;

		for (Mdl = Mdls; Mdl != NULL; Mdl = Mdl -> Pnext) {
		    for (TSrf = Mdl -> TrimSrfList;
			 TSrf != NULL;
			 TSrf = TSrf -> Pnext) {
		        IritTrivTVBndryType Bndry;
			const IritTrivTVStruct *TV;

			TV = IritVMdlGetTVBoundarySrf(VMdl, TSrf -> Srf,
						  &Bndry, NULL);
			if (TV == NULL) {
			    fprintf(stderr, "Failed to find TV\n");
			}
			else {
			    fprintf(stderr,
				    "Found TV - Srf is boundary %d\n",
				    Bndry);
			}
		    }
		}
	    }
	}
    }

    return 0;
}
