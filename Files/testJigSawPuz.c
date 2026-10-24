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

IritPrsrObjectStruct *UserPuz3DTest(IritPrsrObjectStruct *PObj, 
			      const int *StepsMin,
			      const int *StepsMax,
			      CagdVType PolyApproxInfo);

void main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(TRUE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    //IritMiscDynMemoryDbgInitTest2(0x1F, 373);

	    {
	        CagdVType
		    PolyApproxInfo = { 0.02, 0.001, 0.0 };
		IritPrsrObjectStruct
		    *PTiles = UserPuz3DTest(PObj, NULL, NULL,
					    PolyApproxInfo);

		IritPrsrPutObjectToFile3("/temp/Tiles.itd", PTiles, 0);
	    }
	}
    }

    exit(0);
}
