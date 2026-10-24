/******************************************************************************
* testTopo.c - Finding zeroes of a Bezier curve	using factoring.	     *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
*******************************************************************************
* Written by Jinesh Machchhar, July 2015.				     *
******************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/symb_lib.h"
#include "inc_irit/user_lib.h"
#include "inc_irit/mvar_lib.h"
#include "inc_irit/usertopo.h"



/* Create a vector of points. */
static IritUserTopoUnstrctGeomPtStruct *CreatePtVec(CagdPType Origin,
						    CagdRType DistXYZ[3],
						    int NumPtXYZ[3], 
						    int FirstID)
{
    int x, y, z, 
	i = 0;
    CagdPType Pt;
    IritUserTopoUnstrctGeomPtStruct
	*PtVec = (IritUserTopoUnstrctGeomPtStruct *)
	    IritMalloc(sizeof(IritUserTopoUnstrctGeomPtStruct) * (NumPtXYZ[0] *
						  NumPtXYZ[1] * NumPtXYZ[2]));

    IRIT_PT_COPY(Pt, Origin);

    for (z = 0, Pt[2] = Origin[2]; z < NumPtXYZ[2]; z++, Pt[2] += DistXYZ[2]) {
	for (y = 0, Pt[1] = Origin[1]; y < NumPtXYZ[1]; y++, Pt[1] += DistXYZ[1]) {
	    for (x = 0, Pt[0] = Origin[0];
		 x < NumPtXYZ[0];
		 x++, Pt[0] += DistXYZ[0]) {
		IRIT_PT_COPY(PtVec[i].Pt, Pt);
		PtVec[i].ID = FirstID;
		PtVec[i].Attr = NULL;
		PtVec[i].Attr = NULL;
		FirstID++;
		i++;
	    }
	}
    }

    return PtVec;
}

/* Create a linear curve and add to the grid. */
static void AddCrvToGrid(IritUserTopoUnstrctGeomStruct *Ud,
		         int FirstPtId,
		         int NumPtXYZ[3],
		         int LineType,
		         int LineValue1,
		         int LineValue2)
{
    int ID, Ind, NumPt, *PtIdVec,
	x = 0,
	y = 0,
	z = 0,
	CrvInd = 0;
    IritCagdCrvStruct *Crv;
    IritPrsrObjectStruct *IPObj;

    switch (LineType) {
	case 0: y = LineValue1;
	    z = LineValue2;	   
	    NumPt = NumPtXYZ[0];
	    Crv = IritCagdBzrCrvNew(NumPt, IRIT_CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (x = 0; x < NumPtXYZ[0]; x++) {
		Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) +
		                  y * NumPtXYZ[0] + x;
		PtIdVec[CrvInd++] = Ind;
	    }
	    break;
	case 1: x = LineValue1;
	    z = LineValue2;
	    NumPt = NumPtXYZ[1];
	    Crv = IritCagdBzrCrvNew(NumPt, IRIT_CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    Crv = IritCagdBzrCrvNew(NumPt, IRIT_CAGD_PT_E3_TYPE);
	    for (y = 0; y < NumPtXYZ[1]; y++) {
		Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) +
				  y * NumPtXYZ[0] + x;
		PtIdVec[CrvInd++] = Ind;
	    }
	    break;
	case 2: x = LineValue1;
	    y = LineValue2;
	    NumPt = NumPtXYZ[2];
	    Crv = IritCagdBzrCrvNew(NumPt, IRIT_CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    Crv = IritCagdBzrCrvNew(NumPt, IRIT_CAGD_PT_E3_TYPE);
	    for (z = 0; z < NumPtXYZ[2]; z++) {
		Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) +
				  y * NumPtXYZ[0] + x;
		PtIdVec[CrvInd++] = Ind;
	    }
	    break;
	default:
	    assert(0);
    }

    IPObj = IritPrsrGenCRVObject(Crv);    
    IritUserTopoAddCell(Ud, PtIdVec, NumPt, IPObj, &ID);
    IritFree(PtIdVec);
}


/* create a planar surface and add to the grid. */
static void AddSrfToGrid(IritUserTopoUnstrctGeomStruct *Ud,
		         int FirstPtId,
		         int NumPtXYZ[3],
		         int PlaneType,
		         int PlaneValue)
{
    int ID, Ind, NumPt, *PtIdVec,
	x = 0, 
	y = 0, 
	z = 0,
	SrfInd = 0;
    IritCagdSrfStruct *Srf;
    IritPrsrObjectStruct *IPObj;

    assert(PlaneValue >= 0 && PlaneValue < NumPtXYZ[PlaneType]);

    switch (PlaneType) {
	case 0:	x = PlaneValue;
	    NumPt = NumPtXYZ[1] * NumPtXYZ[2];
	    Srf = IritCagdBzrSrfNew(NumPtXYZ[1], NumPtXYZ[2],
				    IRIT_CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (z = 0; z < NumPtXYZ[2]; z++) {
		for (y = 0; y < NumPtXYZ[1]; y++) {
		    Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) +
			              y * NumPtXYZ[0] + x;
		    PtIdVec[SrfInd++] = Ind;
		}
	    }
	    break;
	case 1:	y = PlaneValue;
	    NumPt = NumPtXYZ[0] * NumPtXYZ[2];
	    Srf = IritCagdBzrSrfNew(NumPtXYZ[0], NumPtXYZ[2],
				    IRIT_CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (z = 0; z < NumPtXYZ[2]; z++) {
		for (x = 0; x < NumPtXYZ[0]; x++) {
		    Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) +
			              y * NumPtXYZ[0] + x;
		    PtIdVec[SrfInd++] = Ind;
		}
	    }
	    break;
	case 2: z = PlaneValue;
	    NumPt = NumPtXYZ[0] * NumPtXYZ[1];
	    Srf = IritCagdBzrSrfNew(NumPtXYZ[0], NumPtXYZ[1],
				    IRIT_CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (y = 0; y < NumPtXYZ[1]; y++) {
		for (x = 0; x < NumPtXYZ[0]; x++) {
		    Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) +
				      y * NumPtXYZ[0] + x;
		    PtIdVec[SrfInd++] = Ind;
		}
	    }
	    break;
	default:
	    assert(0);
    }   
    
    IPObj = IritPrsrGenSRFObject(Srf);    
    IritUserTopoAddCell(Ud, PtIdVec, NumPt, IPObj, &ID);
    IritFree(PtIdVec);
}


/* Print all adjacency relations amongst entities of the grid. */
static void PrintAdjacencies(IritUserTopoUnstrctGeomStruct *Ud)
{
    int NumEntAdj, i, EntId, *AdjList;
    IritPrsrObjectStruct *Ent;

    for (Ent = Ud -> CellList; Ent != NULL; Ent = Ent -> Pnext) {
	EntId = IritUserTopoObjectToId(Ud, Ent);
	NumEntAdj = IritUserTopoCellsAdjacentToCell(Ud, EntId, &AdjList);
	printf("%d: ", EntId);
	for (i = 0; i < NumEntAdj; i++) {
	    printf("%d, ", AdjList[i]);
	}
/*	IritFree(AdjList);			dbgheap	*/
	printf("\n");
    }
    
}

/* Render all the entities of the grid onto the given handle. */
static void RenderGridEntities(const IritUserTopoUnstrctGeomStruct *Ud,
			       int ClearScreen,
			       int Handle)
{
    IritPrsrObjectStruct *Obj,
	*ClrObj = IritPrsrGenStrObject("command_", "clear", NULL);

    if (ClearScreen) {
	IritPrsrSocWriteOneObject(Handle, ClrObj);
    }

    assert(Ud);

    for (Obj = Ud -> CellList; Obj != NULL; Obj = Obj -> Pnext) {
	IritPrsrSocWriteOneObject(Handle, Obj);
    }

    IritPrsrFreeObject(ClrObj);
}

/* An example of a call-back function which keeps first 500 points and discards the rest. */
int *UserTopoFilterPtCBFunction(const IritUserTopoUnstrctGeomStruct *Ud)
{
    int i,
	*FilterPtIndx = (int *) IritMalloc(sizeof(int) * Ud -> NumPts);

    for (i = 0; i < Ud -> NumPts; i++) {
	if (i >= 0 && i < 500) {
	    FilterPtIndx[i] = 1;
	}
	else {
	    FilterPtIndx[i] = 0;
	}
    }
    return FilterPtIndx;
}


/* Test case for unstructured grid. */
static void TestCase1(int IOHandle) 
{
    CagdPType Origin;
    CagdRType DistXYZ[3];
    int NumPtXYZ[3], NumPt, *IdVec, *CloneMap,
	FirstPtId = 0;
    IritUserTopoUnstrctGeomPtStruct *PtVec, *PtVec2;
    IritUserTopoUnstrctGeomStruct *Ud1, *Ud2, *Ud12, *Ud12B2, *Ud12B1,
	*Ud12B12, *Ud12F;
    IritPrsrObjectStruct *Ent;

    Origin[0] = Origin[1] = Origin[2] = 0;
    DistXYZ[0] = DistXYZ[1] = DistXYZ[2] = 0.1;
    NumPtXYZ[0] = NumPtXYZ[1] = NumPtXYZ[2] = 10;
    NumPt = NumPtXYZ[0] * NumPtXYZ[1] * NumPtXYZ[2];

    PtVec = CreatePtVec(Origin, DistXYZ, NumPtXYZ, FirstPtId);
    Ud1 = IritUserTopoUnstrctGeomNew();
    IritUserTopoSetPoints(Ud1, PtVec, NumPt); 

    AddCrvToGrid(Ud1, FirstPtId, NumPtXYZ, 0, 1, 1);		   /* Ent 0. */
    Ent = IritUserTopoIdToObject(Ud1, 0);
    IritUserTopoSetCellIntAttr(Ud1, Ent, "color", 1);

    AddCrvToGrid(Ud1, FirstPtId, NumPtXYZ, 1, NumPtXYZ[0] - 1,
		 NumPtXYZ[2] - 1);/* Ent 1. */
    Ent = IritUserTopoIdToObject(Ud1, 1);
    IritUserTopoSetCellIntAttr(Ud1, Ent, "color", 1);

    AddCrvToGrid(Ud1, FirstPtId, NumPtXYZ, 2, NumPtXYZ[0] - 1,
		 NumPtXYZ[1] - 1);/* Ent 2. */
    Ent = IritUserTopoIdToObject(Ud1, 2);
    IritUserTopoSetCellIntAttr(Ud1, Ent, "color", 3);

    AddSrfToGrid(Ud1, FirstPtId, NumPtXYZ, 2, 0);		   /* Ent 3. */
    Ent = IritUserTopoIdToObject(Ud1, 3);
    IritUserTopoSetCellIntAttr(Ud1, Ent, "color", 10);

    IritUserTopoUnstrctGeomUpdate(&Ud1, IRIT_EPS);
    printf("Ud1 adj\n");
    PrintAdjacencies(Ud1);  
    RenderGridEntities(Ud1, TRUE, IOHandle);

    Origin[2] = 0.9;
    NumPtXYZ[2] = 2;
    NumPt = NumPtXYZ[0] * NumPtXYZ[1] * NumPtXYZ[2];
    PtVec2 = CreatePtVec(Origin, DistXYZ, NumPtXYZ, FirstPtId);
    Ud2 = IritUserTopoUnstrctGeomNew();
    IritUserTopoSetPoints(Ud2, PtVec2, NumPt); 

    AddCrvToGrid(Ud2, FirstPtId, NumPtXYZ, 0,
		 NumPtXYZ[1] - 1, 0);/* Ent 0, will become Ent 4. */
    Ent = IritUserTopoIdToObject(Ud2, 0);
    IritUserTopoSetCellIntAttr(Ud2, Ent, "color", 4);
    
    IritUserTopoUnstrctGeomUpdate(&Ud2, IRIT_EPS);
    printf("Ud2 adj\n");
    PrintAdjacencies(Ud2);  
    RenderGridEntities(Ud2, TRUE, IOHandle);
    
    Ud12 = IritUserTopoAppendUnstrctGeoms(Ud1, Ud2, IRIT_EPS,
					  &CloneMap);/* 100 points common in Ud1 and Ud2. */
    printf("Ud12 adj\n");
    PrintAdjacencies(Ud12);
    RenderGridEntities(Ud12, TRUE, IOHandle);
    
    Ud12B1 = IritUserTopoCrvBndryFilter(Ud12);
    printf("Ud12B1 adj\n");
    PrintAdjacencies(Ud12B1);
    RenderGridEntities(Ud12B1, TRUE, IOHandle);

    Ud12B2 = IritUserTopoSrfBndryFilter(Ud12);
    printf("Ud12B2 adj\n");
    PrintAdjacencies(Ud12B2);
    RenderGridEntities(Ud12B2, TRUE, IOHandle);

    Ud12B12 = IritUserTopoAppendUnstrctGeoms(Ud12B2, Ud12B1, IRIT_EPS, &CloneMap);
    printf("Ud12B12 adj\n");
    PrintAdjacencies(Ud12B12);
    RenderGridEntities(Ud12B12, TRUE, IOHandle);

    IritUserTopoGetCellAttrThreshold(Ud12, "color", 3, 4, &IdVec);
    IritFree(IdVec);

    IritUserTopoSetFilterGridCallBackFunc(UserTopoFilterPtCBFunction);
    Ud12F = IritUserTopoApplyFilterToGrid(Ud12, TRUE);
    printf("Ud12F adj\n");
    PrintAdjacencies(Ud12F);
    RenderGridEntities(Ud12F, TRUE, IOHandle);

    IritUserTopoUnstrctGeomFree(Ud1);
    IritUserTopoUnstrctGeomFree(Ud2);
    IritUserTopoUnstrctGeomFree(Ud12);
    IritUserTopoUnstrctGeomFree(Ud12B1);
    IritUserTopoUnstrctGeomFree(Ud12B2);
    IritUserTopoUnstrctGeomFree(Ud12B12);
    IritUserTopoUnstrctGeomFree(Ud12F);
}

int main(int argc, char **argv)
{
    int PrgmIO = 0;
    char *Program = getenv("IRIT_DISPLAY");
    
    IritPrsrSocSrvrInit();            /* Initialize the listen socket for clients. */

#ifdef __WINNT__
    if (Program == NULL)
	Program = "wntgdrvsD -s-";
#endif /* __WINNT__ */
#ifdef __UNIX__
    if (Program == NULL)
	Program = "x11drvs -s-";
#endif /* __UNIX__ */
     

    if ((PrgmIO = IritPrsrSocExecAndConnect(Program, getenv("IRIT_BIN_IPC") !=
	NULL)) >= 0)
    {
	int i;

	TestCase1(PrgmIO);
	scanf("%d", &i);	
    }    

    return 0;
}
