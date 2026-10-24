/*****************************************************************************
*   Exemplifes the synthesis of	bifurcatiob parameteric	tiles (for the	     *
* microstructure constructions using a call back function).		     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 2018   *
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
    IritUserMicroTileBndryPrmStruct WMinPrms, WMax1Prms, WMax2Prms;
    IritPrsrObjectStruct *Tile;

    IRIT_ZAP_MEM(&WMinPrms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMax1Prms, sizeof(IritUserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMax2Prms, sizeof(IritUserMicroTileBndryPrmStruct));

    WMinPrms.Circular = TRUE;
    WMax1Prms.Circular = TRUE;
    WMax2Prms.Circular = TRUE;

    WMinPrms.OuterRadius = 0.3;
    WMax1Prms.OuterRadius = 0.12;
    WMax2Prms.OuterRadius = 0.08;

    Tile = IritUserMicroBifurcate1to2Tile(&WMinPrms, &WMax1Prms, &WMax2Prms,
				      0.5, 0.05, FALSE, NULL);
    IritPrsrPutObjectToFile3("c:/temp/BifTileSrfSolid1.itd", Tile, 0);
    IritPrsrFreeObject(Tile);

    Tile = IritUserMicroBifurcate1to2Tile(&WMinPrms, &WMax1Prms, &WMax2Prms,
				      0.5, 0.05, TRUE, NULL);
    IritPrsrPutObjectToFile3("c:/temp/BifTileTVSolid1.itd", Tile, 0);
    IritPrsrFreeObject(Tile);

    exit(0);
}
