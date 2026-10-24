/*****************************************************************************
*   Exemplifies the synthesis of 3D cross parameteric tiles (for the         *
* microstructure constructions using a call back function).		     *
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

int main(int argc, char **argv)
{
    int i;
    IritUserMicroTileBndryPrmStruct
	UMinPrms, UMaxPrms, VMinPrms, VMaxPrms, WMinPrms, WMaxPrms;
    IritPrsrObjectStruct *Tile;

    IRIT_ZAP_MEM(&UMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&UMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));

    UMinPrms.Circular = TRUE;
    VMinPrms.Circular = TRUE;
    WMaxPrms.Circular = TRUE;
    for (i = 0; i < 4; i++)
	WMinPrms.Bndry[i] = 0.2;
    for (i = 0; i < 4; i++)
	WMaxPrms.Bndry[i] = 0.1;
    UMinPrms.OuterRadius = 0.12;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.12;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.12;
    Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms, FALSE, NULL);
    IritPrsrPutObjectToFile3("c:/temp/TileSolid1.itd", Tile, 0);
    IritPrsrFreeObject(Tile);

    for (i = 0; i < 4; i++)
	VMinPrms.Bndry[i] = 0.14;
    for (i = 0; i < 4; i++)
	VMaxPrms.Bndry[i] = 0.05;
    UMinPrms.OuterRadius = 0.13;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.11;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.1;
    Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms, FALSE, NULL);
    IritPrsrPutObjectToFile3("c:/temp/TileSolid1diff.itd", Tile, 0);
    IritPrsrFreeObject(Tile);

    for (i = 0; i < 4; i++)
	UMinPrms.Bndry[i] = 0.24;
    for (i = 0; i < 4; i++)
	UMaxPrms.Bndry[i] = 0.15;
    UMinPrms.OuterRadius = 0.12;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.12;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.12;
    for (i = 0; i < 4; i++)
	WMinPrms.Bndry[i] = 0.12;
    for (i = 0; i < 4; i++)
	WMaxPrms.Bndry[i] = 0.1;
    UMinPrms.InnerRadius = 0.08;
    UMaxPrms.InnerRadius = 0.08;
    VMinPrms.InnerRadius = 0.08;
    VMaxPrms.InnerRadius = 0.08;
    WMinPrms.InnerRadius = 0.08;
    WMaxPrms.InnerRadius = 0.08;
    Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms, FALSE, NULL);
    IritPrsrPutObjectToFile3("c:/temp/TileHollowed1.itd", Tile, 0);
    IritPrsrFreeObject(Tile);

    for (i = 0; i < 4; i++)
	WMinPrms.Bndry[i] = 0.12;
    for (i = 0; i < 4; i++)
	WMaxPrms.Bndry[i] = 0.1;
    UMinPrms.OuterRadius = 0.13;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.11;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.1;
    Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms, FALSE, NULL);
    IritPrsrPutObjectToFile3("c:/temp/TileHollowed1diff.itd", Tile, 0);
    IritPrsrFreeObject(Tile);

    IRIT_ZAP_MEM(&UMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&UMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMaxPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    UMinPrms.OuterRadius = 0.13;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.11;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.1;

#ifndef DEBUG
    fprintf(stderr, "starting 50000 solids...");
    for (i = 0; i < 50000; i++) {
	Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				    &VMinPrms, &VMaxPrms,
				    &WMinPrms, &WMaxPrms, FALSE, NULL);
        IritPrsrFreeObject(Tile);
    }
    fprintf(stderr, "done\n.");

    UMinPrms.InnerRadius = 0.1;
    UMaxPrms.InnerRadius = 0.11;
    VMinPrms.InnerRadius = 0.1;
    VMaxPrms.InnerRadius = 0.085;
    WMinPrms.InnerRadius = 0.09;
    WMaxPrms.InnerRadius = 0.08;
    for (i = 0; i < 4; i++)
	UMinPrms.Bndry[i] = 0.24;

    fprintf(stderr, "starting 10000 hollowed (and one bndry)...");
    for (i = 0; i < 10000; i++) {
	Tile = IritUserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				    &VMinPrms, &VMaxPrms,
				    &WMinPrms, &WMaxPrms, FALSE, NULL);
        IritPrsrFreeObject(Tile);
    }
    fprintf(stderr, "done\n.");
#endif /* DEBUG */

    exit(0);
}
