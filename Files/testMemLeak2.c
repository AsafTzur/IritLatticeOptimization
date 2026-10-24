/*****************************************************************************
* Filter the data.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#ifndef DEBUG_IRIT_MALLOC
#define DEBUG_IRIT_MALLOC
#endif /* DEBUG_IRIT_MALLOC */

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/geom_lib.h"

int main(int argc, char **argv)
{
    char *p, *q;
    double *r;
    double Start, KeepStep, TrackID;

    Start = 1;
    KeepStep = 0;
    TrackID = 0;
    IritMiscDynMemoryDbgCheckMark(&Start, &KeepStep, &TrackID);

    p = IritMalloc(sizeof(int));
    q = IritMalloc(sizeof(p));
    r = IritMalloc(sizeof(double));

    /* p[-1] = 9; */
    p[5] = 17;
    q[-1] = 0;
    r[0] = 0;
    r[1] = 0;

    IritFree(p);
    IritFree(r);

    Start = 0;
    IritMiscDynMemoryDbgCheckMark(&Start, &KeepStep, &TrackID);

    return 0;
}
