/*****************************************************************************
* Rocket fuel synthesis	using V-rep micro elements.			     *
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
#include "inc_irit/user_lib.h"

typedef struct RocketFuelCallBackStruct {
    const IritCagdCrvStruct *ThrustCurve;/* Scalar pos. func. in t = [0, 1]. */
} RocketFuelCallBackStruct;

static void UserRocketFuelElementDump(int LayerNum,
				      const int NumElements[],
				      int x,
				      int y,
				      CagdRType LayerArea,
				      CagdRType ElementArea,
				      CagdRType ElementDensity);
static CagdRType RocketFuelThrustProfileCB(void *CBData, CagdRType t);
static CagdRType EvalLclFuelThrust(void *CBData,
				   IritCagdSrfStruct *LclElement,
				   CagdRType LclDepth,
				   IritCagdSrfStruct *GlblFront,
				   const int Indices[],
				   const int NumElements[],
				   CagdRType GlblFrontArea);
static IritPrsrObjectStruct *RocketFuelDesign(
			    const IritPrsrObjectStruct *FuelGeometry,
			    const int MeshSize[3],
			    const IritCagdCrvStruct *ThrustCurve,
			    IritUserRocketFuelElementCBFuncType ElementCBFunc);

/*****************************************************************************
* DESCRIPTION:								     M
*   Main.  Executed as "test { -s Nx,Ny,Nz } FuelGeom.itd FuelProfile.itd"   M
*  where:								     M
*  1. (Nx, Ny, Nz) defines the discretization sizes of the fuel	geometry     M
*     (in parameter space, Nz is depth). Default is (10, 10, 10).	     M
*  2. FuelGeom.itd specifies a geometry	of a surface to	be rotated along Z   M
*     to form a	trivariate, or a trivariate directly.			     M
*  3  The geometry of the burning profile.  Only last coordinate (Y for	an   M
*     XY curve)	is evaluated with respect to the parameter time	between	zero M
*     and one.								     M
*									     *
* PARAMETERS:								     M
*   argc: N.S.F.I.							     M
*   argv: N.S.F.I.							     M
*									     *
* RETURN VALUE:								     M
*   void								     M
*									     *
* SEE ALSO:								     M
*									     M
*									     *
* KEYWORDS:								     M
*   main								     M
*****************************************************************************/
int main(int argc, char **argv)
{
    if (argc >= 3) {
	int MeshSize[3] = { 10, 10, 10 };
        IritPrsrObjectStruct *MS, *FuelObj, *BurnProfileObj;

	if (strcmp(argv[1], "-s") == 0) {
	    if (sscanf(argv[2], "%d,%d,%d",
		       &MeshSize[0], &MeshSize[1], &MeshSize[2]) != 3) {
	        fprintf(stderr, "Failed to fetch grid size");
		exit(1);	        
	    }

	    argv += 2;
	    argc -= 2;
	}
	fprintf(stderr, "Grid size is set to (%d, %d, %d)",
		MeshSize[0], MeshSize[1], MeshSize[2]);

	FuelObj = IritPrsrGetDataFiles((const char * const *) &argv[1], 1,
				 TRUE, TRUE);
	if (FuelObj == NULL ||
	    (!IRIT_PRSR_IS_SRF_OBJ(FuelObj) && !IRIT_PRSR_IS_TRIVAR_OBJ(FuelObj))) {
	    fprintf(stderr, "Failed to open fuel (surface or trivariate) file \"%s\"\n", argv[1]);
	    exit(1);
	}

	BurnProfileObj = IritPrsrGetDataFiles((const char * const *) &argv[2],
					1, TRUE, TRUE);
	if (BurnProfileObj == NULL || !IRIT_PRSR_IS_CRV_OBJ(BurnProfileObj)) {
	    fprintf(stderr, "Failed to open thrust profile file \"%s\"\n",
		    argv[2]);
	    exit(1);
	}

	MS = RocketFuelDesign(FuelObj, MeshSize, BurnProfileObj -> U.Crvs,
			      UserRocketFuelElementDump);

	IritPrsrStdoutObject(MS, FALSE);
	IritPrsrFreeObject(MS);
    }
    else {
	fprintf(stderr, "Usage: test { -s Nx,Ny,Nz } FeulGeom.itd FuelProfileCrv.itd.\n");
	exit(2);
    }

    return 0;
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Dumps to stdout the	information on all individual elements,	layer by     *
* layer.								     *
*									     *
* PARAMETERS:								     *
*   LayerNum:	  Index	of later.  Starting from zero.			     *
*   NumElements:  Number of elements in	a layer, in x &	y (array of size 2). *
*   x, y:	  Indices of this element in the current layer.		     *
*   LayerArea:	  The area of the whole	layer.				     *
*   ElementArea:  The area of this element.				     *
*   ElementDensity: The	density	of this	element.			     *
*									     *
* RETURN VALUE:								     *
*   void								     *
*****************************************************************************/
static void UserRocketFuelElementDump(int LayerNum,
				      const int	NumElements[],
				      int x,
				      int y,
				      CagdRType	LayerArea,
				      CagdRType	ElementArea,
				      CagdRType	ElementDensity)
{
    static int
        Layer = -1;
    static CagdRType
        AccumArea = 0.0;

    if (Layer != LayerNum) {
        if (LayerNum != 0) {
	    fprintf(stderr, "Layer %d accumulated area is %f\n",
		    LayerNum - 1, AccumArea);
	}
	AccumArea = 0.0;
	Layer = LayerNum;
    }
    AccumArea += ElementArea;

    fprintf(stderr,
	   "Layer %d (%d x %d): Area = %f, Element (%d, %d): Area = %f, Density = %f\n",
	    LayerNum, NumElements[0], NumElements[1], LayerArea,
	    x, y, ElementArea,ElementDensity);
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Given time,	between	zero and one, returns the required thrust at that    *
* time.	 Thrust	is extracted from the highest dimensions in the	curve.	     *
*									     *
* PARAMETERS:								     *
*   CBdata:  An	optional call back data	that can be passed by the main	     *
*	     function to these call back functions.			     *
*   t:	     Time.  Normalized to zero -> one.				     *
*									     *
* RETURN VALUE:								     *
*   CagdRType:	  Required Thrust.					     *
*****************************************************************************/
static CagdRType RocketFuelThrustProfileCB(void *CBData, CagdRType t)
{
    CagdPType Thrust;
    RocketFuelCallBackStruct
	*RFCB = (RocketFuelCallBackStruct *) CBData;
    const IritCagdCrvStruct 
	*Crv = RFCB -> ThrustCurve;
    int Dim = IRIT_CAGD_NUM_OF_PT_COORD(Crv -> PType);

    if (Dim > 3) {
	fprintf(stderr, "RocketFuel: Expected a thrust curve in 3D, not more...\n");
	return 1.0; /* Does not support this... */
    }

    IRIT_CAGD_CRV_EVAL_E3(Crv, t, Thrust);

    if (Thrust[Dim - 1] < IRIT_EPS) {
	fprintf(stderr, "RocketFuel: Warning:  Thrust at t = %f is %f\n",
		t, Thrust[Dim - 1]);
    }

    fprintf(stderr, "RocketFuel: Thrust at t = %f is %f\n",
	    t, Thrust[Dim - 1]);


    return Thrust[Dim - 1];
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Given the local micro-element front	surface	geometry, and the entire     *
* front's surface geometry, returns the	thrust that this small element will  *
* generate.  A default function	in case	no such	function is provided.	     *
*									     *
* PARAMETERS:								     *
*   CBdata:	   An optional call back data that can be passed by the	main *
*		   function to these call back functions.		     *
*   LclElement:	   The front surface of	the current micro element.	     *
*   LclDepth:	   The local layer depth of this LclElement.		     *
*   GlblFront:	   The entire front surface.				     *
*   Indices:	   2D indices of this element in layer.			     *
*		   Between 0 and < NumElements[], as (IdxX, IdxY).	     *
*   NumElements:   How many elements to	place in each layer as (Nx, Ny).     *
*   GlblFrontArea: The surface area of GlblFront.			     *
*									     *
* RETURN VALUE:								     *
*   CagdRType:	   Contribution	thrust of the local element.		     *
*****************************************************************************/
static CagdRType EvalLclFuelThrust(void *CBData,
				   IritCagdSrfStruct *LclElement,
				   CagdRType LclDepth,
				   IritCagdSrfStruct *GlblFront,
				   const int Indices[],
				   const int NumElements[],
				   CagdRType GlblFrontArea)
{
    return LclDepth * GlblFrontArea / (NumElements[0] * NumElements[1]);
}

/*****************************************************************************
* DESCRIPTION:								     *
*   Tool to design rocket function out of heterogeneous	material following   *
* burring profile as specified by ThrustCurve.				     *
*									     *
* PARAMETERS:								     *
*   FuelGeometry: Either A 2D surface cross section to revolve around Z	to   *
*		  create the basic trivariate fuel (as a volume	of	     *
*		  revolution CrossSection (V must be the expected burning    *
*		  direction then, from VMin, to	VMax), or		     *
*		  a 3D full trivariate,	in which case the burning is	     *
*		  expected from	WMin to	WMax.				     *
*   MeshSize:	  A vector holding (u, v, w) trivariate	division into	     *
		  micro-elements, w bing the depth burning (WMin to WMax).   *
*   ThrustCurve:  A thrust function (using the largest axis in the curve)    *
*		  over domain [0, 1].  Can be nil() for	uniform	burn.	     *
*   ElementCBFunc:  A call back	function to be invoked on every	element.     *
*									     *
* RETURN VALUE:								     *
*   IritPrsrObjectStruct *:  The computed fuel's heterogeneous geometry	     *
*****************************************************************************/
static IritPrsrObjectStruct *RocketFuelDesign(
			    const IritPrsrObjectStruct *FuelGeometry,
			    const int MeshSize[3],
			    const IritCagdCrvStruct *ThrustCurve,
			    IritUserRocketFuelElementCBFuncType ElementCBFunc)
{
    int BndrySrfs, SliceThrough;
    CagdRType
	ApplyRGB = 1.0;
    IritTrivTVStruct *TV;
    IritUserRocketFuelThrustProfileFuncType ThrustProfileFunc;
    IritPrsrObjectStruct *FuelGeom, *FuelObj;
    RocketFuelCallBackStruct RFCB;

    IRIT_ZAP_MEM(&RFCB, sizeof(RocketFuelCallBackStruct));

    if (IRIT_PRSR_IS_SRF_OBJ(FuelGeometry)) {
	if (IRIT_CAGD_IS_RATIONAL_SRF(FuelGeometry -> U.Srfs))
	    TV = IritTrivTVOfRev(FuelGeometry -> U.Srfs);
	else
	    TV = IritTrivTVOfRevPolynomialApprox(FuelGeometry -> U.Srfs);
    }
    else {
	assert(IRIT_PRSR_IS_TRIVAR_OBJ(FuelGeometry));
	TV = IritTrivTVCopy(FuelGeometry -> U.Trivars);
    }

    if (ThrustCurve != NULL) {
        RFCB.ThrustCurve = ThrustCurve;
	ThrustProfileFunc = RocketFuelThrustProfileCB;
    }
    else
	ThrustProfileFunc = NULL;

    /* BndrtSrf - TRUE for only surfaces,				     */
    /*            FALSE for trivariate data (slower).			     */
    /* SliceThru - TRUE to slice out a pie of solid fuel to see a cross      */
    /*              section.  FALSE to provide the full 360 degrees of fuel. */
    BndrySrfs = TRUE;
    SliceThrough = TRUE;
    FuelGeom = IritUserRocketFuelDesign3D(TV, MeshSize[2], MeshSize,
				      SliceThrough, BndrySrfs, &ApplyRGB,
				      ThrustProfileFunc, EvalLclFuelThrust,
				      NULL, NULL, UserRocketFuelElementDump,
				      &RFCB);

    /* Convert to object per geometric entity. */
    if (BndrySrfs) {
	FuelObj = IritPrsrLnkListToListObject(FuelGeom -> U.Srfs,
					      IRIT_PRSR_OBJ_SURFACE);
	FuelGeom -> U.Srfs = NULL;
    }
    else {
	FuelObj = IritPrsrLnkListToListObject(FuelGeom -> U.Trivars,
					      IRIT_PRSR_OBJ_TRIVAR);
	FuelGeom -> U.Trivars = NULL;
    }
    IritPrsrFreeObject(FuelGeom);

    return FuelObj;
}
