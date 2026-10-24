/*****************************************************************************
* Trivariate Boolean sum example.					     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, Dec 2017    *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

/* A boolean sum between 6 planes forming a unit [0, 1]^3 cube. */
IritTrivTVStruct *BoolSum1(void)
{
    IritCagdSrfStruct
	*Srf1 = IritCagdPrimPlaneSrf(0, 0, 1, 1, 1), /* Z = 1. */
	*Srf2 = IritCagdPrimPlaneSrf(0, 0, 1, 1, 0), /* Z = 0. */
	*Srf3 = IritCagdPrimPlaneXZSrf(0, 0, 1, 1, 1), /* Y = 1. */
	*Srf4 = IritCagdPrimPlaneXZSrf(0, 0, 1, 1, 0), /* Y = 0. */
	*Srf5 = IritCagdPrimPlaneYZSrf(0, 0, 1, 1, 1), /* X = 1. */
	*Srf6 = IritCagdPrimPlaneYZSrf(0, 0, 1, 1, 0); /* X = 0. */
    IritTrivTVStruct
	*TV = IritMvarTrivarBoolSum3(Srf1, Srf2, Srf3, Srf4, Srf5, Srf6);

    IritCagdSrfFree(Srf1);
    IritCagdSrfFree(Srf2);
    IritCagdSrfFree(Srf3);
    IritCagdSrfFree(Srf4);
    IritCagdSrfFree(Srf5);
    IritCagdSrfFree(Srf6);

    return TV;
}

/* Extract four boundary faces from an open cylinder (no caps) and use them */
/* to construct a trivariate cylinder with caps computed automatically.     */
IritTrivTVStruct *BoolSum2(void)
{
    CagdVType 
	Center = { -1.0, -1.0, 0.0 };
    IritCagdSrfStruct
	*CylSrf = IritCagdPrimCylinderSrf(Center, 1.0, 1.0, FALSE, FALSE),
	*Srf1 = IritCagdSrfRegionFromSrf(CylSrf, 0.0, 1.0,
					 IRIT_CAGD_CONST_U_DIR),
	*Srf2 = IritCagdSrfRegionFromSrf(CylSrf, 1.0, 2.0,
					 IRIT_CAGD_CONST_U_DIR),
	*Srf3 = IritCagdSrfRegionFromSrf(CylSrf, 2.0, 3.0,
					 IRIT_CAGD_CONST_U_DIR),
	*Srf4 = IritCagdSrfRegionFromSrf(CylSrf, 3.0, 4.0,
					 IRIT_CAGD_CONST_U_DIR);
    IritTrivTVStruct
	*TV = IritMvarTrivarBoolSum3(Srf1, Srf2, Srf3, Srf4, NULL, NULL);

    IritCagdSrfFree(Srf1);
    IritCagdSrfFree(Srf2);
    IritCagdSrfFree(Srf3);
    IritCagdSrfFree(Srf4);
    IritCagdSrfFree(CylSrf);

    return TV;
}

int main(int argc, char **argv)
{
    int FHandler;
    IritPrsrObjectStruct
	*TV1Obj = IritPrsrGenTRIVARObject(BoolSum1()),
	*TV2Obj = IritPrsrGenTRIVARObject(BoolSum2());

    FHandler = IritPrsrOpenDataFile("/temp/TVs.itd", FALSE, 1); 
    IritPrsrPutObjectToHandler(FHandler, TV1Obj);
    IritPrsrPutObjectToHandler(FHandler, TV2Obj);
    IritPrsrFreeObject(TV1Obj);
    IritPrsrFreeObject(TV1Obj);
    IritPrsrCloseStream(FHandler, TRUE);

    return 0;
}
