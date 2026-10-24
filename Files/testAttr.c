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

static void AddAttributes(IritPrsrObjectStruct *PObj,
			  IrtHmgnMatType Mat,
			  void *Data)
{
    switch (PObj -> ObjType) {
	case IRIT_PRSR_OBJ_CURVE:
	    IritMiscAttrSetIntAttrib(&PObj -> U.Crvs -> Attr,
				     "IntAttr", 10000);
	    break;
	case IRIT_PRSR_OBJ_SURFACE:
	    IritMiscAttrSetIntAttrib(&PObj -> U.Srfs -> Attr,
				     "IntAttr", 10001);
	    break;
	case IRIT_PRSR_OBJ_TRISRF:
	    IritMiscAttrSetIntAttrib(&PObj -> U.TriSrfs -> Attr,
				     "IntAttr", 10002);
	    break;
	case IRIT_PRSR_OBJ_TRIMSRF:
	    IritMiscAttrSetIntAttrib(&PObj -> U.TrimSrfs -> Attr,
				     "IntAttr", 10010);
	    IritMiscAttrSetIntAttrib(&PObj -> U.TrimSrfs -> Srf -> Attr,
				     "IntAttr",
				     10011);
	    IritMiscAttrSetIntAttrib(
			      &PObj -> U.TrimSrfs -> TrimCrvList -> Attr,
			     "IntAttr", 10012);
	    IritMiscAttrSetIntAttrib(&PObj -> U.TrimSrfs -> TrimCrvList ->
							TrimCrvSegList -> Attr,
			     "IntAttr", 10013); 
	    IritMiscAttrSetIntAttrib(&PObj -> U.TrimSrfs -> TrimCrvList ->
					       TrimCrvSegList -> UVCrv -> Attr,
			     "IntAttr", 10014);
	    break;
	case IRIT_PRSR_OBJ_TRIVAR:
	    IritMiscAttrSetIntAttrib(&PObj -> U.Srfs -> Attr,
				     "IntAttr", 10020);
	    break;
	case IRIT_PRSR_OBJ_MODEL:
	    IritMiscAttrSetIntAttrib(&PObj -> U.Mdls -> Attr,
				     "IntAttr", 10030);
	    IritMiscAttrSetIntAttrib(&PObj -> U.Mdls -> TrimSegList -> Attr,
			     "IntAttr", 10031);
	    IritMiscAttrSetIntAttrib(
			 &PObj -> U.Mdls -> TrimSegList -> UVCrvFirst -> Attr,
			 "IntAttr", 10032);
	    IritMiscAttrSetIntAttrib(&PObj -> U.Mdls -> TrimSrfList -> Attr,
			     "IntAttr", 10033);
	    IritMiscAttrSetIntAttrib(
			     &PObj -> U.Mdls -> TrimSrfList -> Srf -> Attr,
			     "IntAttr", 10034);
	    IritMiscAttrSetIntAttrib(
			    &PObj -> U.Mdls -> TrimSrfList -> LoopList -> Attr,
			    "IntAttr", 10035);
	    break;
	case IRIT_PRSR_OBJ_MULTIVAR:
	    IritMiscAttrSetIntAttrib(&PObj -> U.MultiVars -> Attr, 
				     "IntAttr", 10040);
	    break;
        default:
	    break;
    }
}

int main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IritPrsrSetFlattenObjects(FALSE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    const char
		*TmpName = "/temp/AttrObjs.itd";
	    IritPrsrObjectStruct *PObj2, *PObj;
	    IritPrsrTraverseObjHierarchyStruct TraversState;

	    PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    IritPrsrTraverseObjHierarchyInitState(&TraversState);
	    TraversState.ApplyFunc = AddAttributes;
	    IritPrsrTraverseObjListHierarchy(PObj, &TraversState);

	    /* Save data with attributes. */
	    IritPrsrPutObjectToFile3(TmpName, PObj, 0);

	    /* And reload. */
	    PObj2 = IritPrsrGetDataFiles(&TmpName, 1, TRUE, TRUE);
	    IritPrsrStdoutObject(PObj2, FALSE);

	    /* Now read the data. */
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: IOAttr geom.itd\n");
	exit(2);
    }

    return 0;
}
