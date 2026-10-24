/*****************************************************************************
*   Constructs locally varying tiles in	microstructure constructions using   *
* a call back function.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, Dec 2017    *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/mdl_lib.h"

int main(int argc, char **argv)
{
    CagdVType
        SphereCenter = { 0.0, 0.0, 0.0 },
        ConeCenter = { 0.2, 0.0, -2.0 };
    IritCagdSrfStruct
        *Srf1 = IritCagdPrimSphereSrf(SphereCenter, 1.0, FALSE),
        *Srf2 = IritCagdPrimConeSrf(ConeCenter, 1.0, 4, FALSE, TRUE);
    IritMdlModelStruct
        *Mdl1 = IritMdlCnvrtSrf2Mdl(Srf1, TRUE),
        *Mdl2 = IritMdlCnvrtSrf2Mdl(Srf2, TRUE);
    IritPrsrObjectStruct *Inter, *Union, *Subtr;

    IritCagdSrfFree(Srf1);
    IritCagdSrfFree(Srf2);

    Inter = IritMdlBooleanIntersection(Mdl1, Mdl2, NULL, NULL);
    Subtr = IritMdlBooleanSubtraction(Mdl1, Mdl2, NULL, NULL);
    Union = IritMdlBooleanUnion(Mdl1, Mdl2, NULL, NULL);

    IritMdlModelFree(Mdl1);
    IritMdlModelFree(Mdl2);

    IritPrsrPutObjectToFile3("InterMdlBool.itd", Inter, 0);
    IritPrsrPutObjectToFile3("SubtrMdlBool.itd", Subtr, 0);
    IritPrsrPutObjectToFile3("UnionMdlBool.itd", Union, 0);

    IritPrsrFreeObject(Inter);
    IritPrsrFreeObject(Subtr);
    IritPrsrFreeObject(Union);
}
