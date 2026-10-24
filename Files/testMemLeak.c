/*****************************************************************************
* Filter the data.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#define DEBUG_IRIT_MALLOC

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"

int main(int argc, char **argv)
{
    /* START */
    int i;
    char *Name = "c:/temp/trim_srf.itd";
    IritPrsrObjectStruct 
	*trim_srf = IritPrsrGetDataFiles((const char * const *) &Name,
					 1, 1, 1);


    {
	IrtRType ri, rj, rk;

	ri = 1;
	rj = 30;
	rk = 0;
	IritMiscDynMemoryDbgCheckMark(&ri, &rj, &rk);
    }

    for (i = 0; i < 100; i++) {
	IritCagdSrfStruct
	    *untrim_srfs_phys = IritTrimSrfCnvrt2BzrRglrSrf2(
						       trim_srf -> U.TrimSrfs,
						       TRUE,TRUE,1.0e-15);
	IritCagdSrfFreeList(untrim_srfs_phys);
    }


    {
	IrtRType ri, rj, rk;

	ri = 0;
	rj = 30;
	rk = 0;
	IritMiscDynMemoryDbgCheckMark(&ri, &rj, &rk);
    }
    IritPrsrFreeObject(trim_srf);

    return 0;    /* END */
}
