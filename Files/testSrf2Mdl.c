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
#include "inc_irit/iritprsr.h"
#include <time.h>

void main(int argc, char **argv)
{   
    int Handler, Is2Manifold, StitchedCrvs,
	DRes = 0;
    IrtRType TimeTaken;
    const char
	*FileName = "/temp/MS1.itd";
    IritCagdSrfStruct *Srfs;
    IritPrsrObjectStruct *IPObj;
    IritMdlModelStruct* Mdl;
    clock_t t;

    Handler = IritPrsrOpenDataFile(FileName, TRUE, TRUE);
    IPObj = IritPrsrGetObjects(Handler);
    Srfs = IPObj -> U.Srfs;
    t = clock();
    Mdl = IritMdlPrimListOfSrfs2Model2(Srfs, &Is2Manifold, &StitchedCrvs);
    t = clock() - t;
    TimeTaken = ((IrtRType) t) / CLOCKS_PER_SEC; // in seconds
    printf( "Stitching computed in %lf secs.\n", TimeTaken);

    if (Mdl == NULL) {
	fprintf(stderr, "Failed to build a model from the surfaces\n");
    }

//    DRes = IritMdlDebugVerify(Mdl, TRUE, TRUE);
    fprintf(stderr, "Model is 2-Manifold and is %s\nStitched Crvs Num %d\nVerification Res %d\n",
	   Is2Manifold == 1 ? "open" : "closed", StitchedCrvs, DRes);
}
