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

/*****************************************************************************
* DESCRIPTION:								     M
*   Slice volumetric geometry.						     M
*									     *
* PARAMETERS:								     M
*   VModel:    The vmodel to slice with	a plane.			     M
*									     *
* RETURN VALUE:								     M
*   void								     M
*									     *
* SEE ALSO:								     M
*									     M
*									     *
* KEYWORDS:								     M
*   IrtVMODELPlanarSlicer						     M
*****************************************************************************/
void IrtVMODELPlanarSlicer(const IritPrsrObjectStruct *VModel)
{
    int Idx, x, y, SliceSize[2];
    IritCagdBBoxStruct BBox;
    IritVMdlParamsStruct Params;
    struct IritVMdlSlicerInfoStruct *Info;
    IritMiscImgPixelStruct *ImageLine, **Image;
    IritMiscWriteGenInfoStructPtr GI;
    IritVMdlVolumeElementStruct *VolumeElement;

    IritVMdlVModelBBox(NULL, VModel -> U.VMdls, &BBox, TRUE);

    IritVMdlDefaultParams(&Params);
    Params.SlicerParams.Min[0] = BBox.Min[0];
    Params.SlicerParams.Min[1] = BBox.Min[1];
    Params.SlicerParams.Max[0] = BBox.Max[0];
    Params.SlicerParams.Max[1] = BBox.Max[0];
    Params.SlicerParams.ZRes = 0.01;
    Params.SlicerParams.XYRes[0] = 0.01;
    Params.SlicerParams.XYRes[1] = 0.01;

    SliceSize[0] = SliceSize[1] = 200;

    Info = IritVMdlSlicerInitVModel(&Params, VModel -> U.VMdls);
    IritVMdlSlicerGetSliceSize(Info, SliceSize);    /* Fetch the slice size. */

    /* Allocate Image info. */
    Image = (IritMiscImgPixelStruct **)
	          IritMalloc(sizeof(IritMiscImgPixelStruct *) * SliceSize[1]);
    for (y = 0; y < SliceSize[1]; y++) {
        Image[y] = ImageLine = (IritMiscImgPixelStruct *)
	            IritMalloc(sizeof(IritMiscImgPixelStruct) * SliceSize[0]);
	for (x = 0; x < SliceSize[0]; x++)	         /* Init background. */
	    ImageLine[x].r = ImageLine[x].g = ImageLine[x].b = 255;
    }

    /* Traverse the VCells in the VModel. */
    for (Idx = 0,VolumeElement = VModel -> U.VMdls -> VolumeElements;
	 VolumeElement != NULL;
	 Idx++, VolumeElement = VolumeElement -> Pnext) {
        IritVMdlSlicerInitVModelVElement(&Params, Info, VModel -> U.VMdls,
				     VolumeElement);

	/* Compute the slice. */
	IritVMdlSlicerSliceAtZLevel(NULL, Info, 0.0, 0);

	for (y = 0; y < SliceSize[1]; y++) {
	    ImageLine = Image[y];

	    for (x = 0; x < SliceSize[0]; x++) {
		CagdPType UVWParams, Pos;

		if (IritVMdlSlicerGetCurrSliceXY(Info, x, y, UVWParams,
						 Pos, NULL)) {
		    ImageLine[x].r = Idx & 0x01 ? 50 : 200;
		    ImageLine[x].g = Idx & 0x02 ? 50 : 200;
		    ImageLine[x].b = Idx & 0x04 ? 50 : 200;
		}
	    }
	}
    }

    /* Save the image. */
    if ((GI = IritMiscImgWriteOpenFile(NULL, "c:/temp/slice.ppm",
				       IRIT_MISC_IMAGE_UNKNOWN_TYPE,
				       FALSE, SliceSize[0],
				       SliceSize[1])) == NULL)
        return;

    for (y = 0; y < SliceSize[0]; y++) {
        IritMiscImgWritePutLine(GI, NULL, Image[y]);
	IritFree(Image[y]);
    }
    IritFree(Image);

    IritMiscImgWriteCloseFile(GI);

    IritVMdlSlicerFree(Info);
}

int main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
	IritPrsrSetFlattenObjects(TRUE);
	if ((Handler = IritPrsrOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IritPrsrObjectStruct *PTmp,
		*PObj = IritPrsrGetObjects(Handler);

	    /* Done with file - close it. */
	    IritPrsrCloseStream(Handler, TRUE);

	    /* Filter the geometry. */
	    for (PTmp = PObj; PTmp != NULL; PTmp = PTmp -> Pnext) {
		if (IRIT_PRSR_IS_VMODEL_OBJ(PTmp)) {
		    IrtVMODELPlanarSlicer(PTmp);
		}
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: PolyArea geom.itd\n");
	exit(2);
    }

    return 0;
}
