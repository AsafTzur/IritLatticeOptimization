/*****************************************************************************
* Filter the data.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include <stdlib.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/vmdl_lib.h"

/***************** Declarations of user-defined functions. *****************/

static void PrintChar(IritVMdlSlicerOutputImageStruct *Image,
		      const char *Output,
		      const int *Resolution);
static IritCagdPointType Matrix3x3PropertyFunction(IritVMdlBlendPrepStruct
					                      *BlendCoordPrep,
					           const int ObjIndex,
					           void *Output);
static IritCagdPointType CharPropertyFunction(IritVMdlBlendPrepStruct
					                      *BlendCoordPrep,
					      const int ObjIndex,
					      void *Output);
static void *Matrix3x3BlendingRule(
			      IritVMdlBlendPrepStruct *BlendCoordPrep,
			      const IrtRType Pt[3],
			      void *OutputPtr,
			      IritVMdlBlendCoordinatesFuncType Coordinates,
			      IritVMdlBlendPropertyFuncType PropertyFunction);
static void *CharBlendingRule(IritVMdlBlendPrepStruct *BlendCoordPrep,
			      const IrtRType Pt[3],
			      void *OutputPtr,
			      IritVMdlBlendCoordinatesFuncType Coordinates,
			      IritVMdlBlendPropertyFuncType PropertyFunction);
static void PrintMatrix(IritVMdlSlicerOutputImageStruct *Image,
                        const char *Output,
                        const int *Resolution);

/***************** Definitions of user-defined functions. ******************/
/*       A contact Matrix per ObjIndex (no use of BlendCoordPrep).         */
static IritCagdPointType Matrix3x3PropertyFunction(IritVMdlBlendPrepStruct
					                    *BlendCoordPrep,
					           const int ObjIndex,
					           void *Output)
{
    if (Output != NULL) {
        int i, j;
        IrtRType
            (*Ptr)[3][3] = (IrtRType(*)[3][3]) Output,
            (*Matrix)[3] = *Ptr,
            Counter[3] = { 0.0,0.0,0.0 };

        for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j) {
                Matrix[i][j] = (IrtRType) rand();
                Counter[i] += Matrix[i][j];
            }

        for (i = 0; i < 3; ++i) {
            if (ObjIndex == 0) {
                if (i == 0)
                    for (j = 0; j < 3; ++j)
                        Matrix[i][j] *= 255.0 / Counter[i];
                else if (i == 1)
                    for (j = 0; j < 3; ++j)
                        Matrix[i][j] *= 0.0;
                else
                    for (j = 0; j < 3; ++j)
                        Matrix[i][j] *= 100.0 / Counter[i];
            }
            else {
                if (i == 0)
                    for (j = 0; j < 3; ++j)
                        Matrix[i][j] *= 10.0 / Counter[i];
                else if (i == 1)
                    for (j = 0; j < 3; ++j)
                        Matrix[i][j] *= 255.0 / Counter[i];
                else
                    for (j = 0; j < 3; ++j)
                        Matrix[i][j] *= 125.0 / Counter[i];
            }
        }
    }

    return IRIT_CAGD_MAKE_PT_TYPE(FALSE, 0);
}

/***************** Definitions of user-defined functions. ******************/
/*        A contact char per ObjIndex (no use of BlendCoordPrep).          */
static IritCagdPointType CharPropertyFunction(IritVMdlBlendPrepStruct
					                  *BlendCoordPrep,
					      const int ObjIndex,
					      void *Output)
{
    if (Output != NULL) {
        char
            *Ptr = (char *) Output;

        if (ObjIndex == 0)
            *Ptr = 'a';
        else
            *Ptr = 'g';
    }

    return IRIT_CAGD_MAKE_PT_TYPE(FALSE, 0);
}

/***************** Definitions of user-defined functions. ******************/
/*       A contact vector per ObjIndex (no use of BlendCoordPrep).         */
#ifdef IRIT_NOT_USED
static IritCagdPointType Vector3PropertyFunction(IritVMdlBlendPrepStruct
					                  *BlendCoordPrep,
					         const int ObjIndex,
					         void *Output)
{
    if (Output != NULL) {
        int i;
        IrtRType
        (*Ptr)[3] = (IrtRType(*)[3]) Output,
            *Vector = *Ptr;

        for (i = 0; i < 3; ++i) {
            if (ObjIndex == 0) {
                if (i == 0)
                    Vector[i] = 255.0;
                else if (i == 1)
                    Vector[i] = 0.0;
                else
                    Vector[i] = 100.0;
            }
            else {
                if (i == 0)
                    Vector[i] = 10.0;
                else if (i == 1)
                    Vector[i] = 255.0;
                else
                    Vector[i] = 125.0;
            }
        }
    }

    return IRIT_CAGD_MAKE_PT_TYPE(FALSE, 3);
}
#endif /* IRIT_NOT_USED */

/************************ Blending for matrices. ***************************/
static void *Matrix3x3BlendingRule(
			     IritVMdlBlendPrepStruct *BlendCoordPrep,
			     const IrtRType Pt[3],
			     void *OutputPtr,
			     IritVMdlBlendCoordinatesFuncType Coordinates,
			     IritVMdlBlendPropertyFuncType PropertyFunction)
{
    const size_t
        SizeTopology = BlendCoordPrep -> NumVPrims;
    int inTriv, i, j,
        ObjIndex = -1;  /* If pt belongs to only one VMdl, this is its idx. */
    IrtRType
        StandardValue = 255.0 / 3;  /* What to assign in case nothing else. */
    IritVMdlBlendTopologyStruct *TopolStructIter,
        *TopolStructStart = BlendCoordPrep -> TopolStruct;
    IritVMdlBlendPointStruct
        *Point = BlendCoordPrep -> Point;
    IrtRType
        (*PropertyAux)[3][3] = NULL,
        (*Ptr)[3][3] = (IrtRType(*)[3][3]) OutputPtr;
    IrtRType
        (*Matrix)[3] = *Ptr;

    if (!Coordinates)
        Coordinates = IritVMdlBlendComputeDistCoordinates;
    if (!PropertyFunction)
        PropertyFunction = IritVMdlBlendDfltRGBPropertyFunction;
    if (!Pt) {
        for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                Matrix[i][j] = StandardValue;

        return OutputPtr;
    }

    inTriv = IritVMdlBlendInIntersection(BlendCoordPrep, &ObjIndex);

    /* P belongs to at least a TV. */
    if (inTriv == 0) {
        IRIT_ZAP_MEM(Point -> BlendWeights, SizeTopology * sizeof(CagdRType));
        Point -> BlendWeights[ObjIndex] = 1.0;

        PropertyFunction(BlendCoordPrep, ObjIndex, OutputPtr);

        return OutputPtr;
    }
    else if (inTriv < 0) {
        /* P is outside the V-model. */

        for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                Matrix[i][j] = StandardValue;

        return OutputPtr;
    }
    else {
        /* P belongs to a transition region. */
        if (Point -> BlendWeights[0] < 0)
            Coordinates(BlendCoordPrep, Pt);
        PropertyAux = (IrtRType(*)[3][3])
            IritMalloc(sizeof(*PropertyAux));
    }

    /* Blending the values of the relevant properties. */
    for (i = 0; i < 3; ++i)
        IRIT_ZAP_MEM(Matrix[i], 3 * sizeof(IrtRType));

    for (TopolStructIter = TopolStructStart;
        TopolStructIter != NULL;
        TopolStructIter = TopolStructIter -> Pnext) {
        size_t
            IndexObj = TopolStructIter -> Index;
        CagdRType
            BlendWeight = Point -> BlendWeights[IndexObj];
        IrtRType
	    (*MatrixAux)[3] = *PropertyAux;

        if (!Point -> InclusionVPrims[IndexObj])
            continue;

        PropertyFunction(BlendCoordPrep, (int) IndexObj,
			 (void *) PropertyAux);

        for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                Matrix[i][j] += BlendWeight * MatrixAux[i][j];
    }

    if (PropertyAux)
        IritFree(PropertyAux);

    return OutputPtr;
}

/************************** Blending for chars. *****************************/
static void *CharBlendingRule(IritVMdlBlendPrepStruct *BlendCoordPrep,
			      const IrtRType Pt[3],
			      void *OutputPtr,
			      IritVMdlBlendCoordinatesFuncType Coordinates,
			      IritVMdlBlendPropertyFuncType PropertyFunction)
{
    const size_t
        SizeTopology = BlendCoordPrep -> NumVPrims;
    int inTriv,
        ObjIndex = -1;  /* If pt belongs to only one VMdl, this is its idx. */
    char
        StandardValue = 'x';        /* What to assign in case nothing else. */
    IritVMdlBlendPointStruct
        *Point = BlendCoordPrep -> Point;
    char
        *PropertyAux = NULL,
        *Ptr = (char *) OutputPtr;

    if (!Coordinates)
        Coordinates = IritVMdlBlendComputeDistCoordinates;
    if (!PropertyFunction)
        PropertyFunction = IritVMdlBlendDfltRGBPropertyFunction;
    if (!Pt) {
        *Ptr = StandardValue;

        return OutputPtr;
    }

    inTriv = IritVMdlBlendInIntersection(BlendCoordPrep, &ObjIndex);

    /* P belongs to at least a TV. */
    if (inTriv == 0) {
        IRIT_ZAP_MEM(Point -> BlendWeights, SizeTopology * sizeof(CagdRType));
        Point -> BlendWeights[ObjIndex] = 1.0;

        PropertyFunction(BlendCoordPrep, ObjIndex, OutputPtr);

        return OutputPtr;
    }
    else if (inTriv < 0) {
        /* P is outside the V-model. */
        return OutputPtr;
    }
    else {
        /* P belongs to a transition region. */
        if (Point -> BlendWeights[0] < 0)
            Coordinates(BlendCoordPrep, Pt);
        PropertyAux = (char *)
            IritMalloc(sizeof(char));
    }

    /* Blending the values of the relevant properties. */
    if (0 < Point -> BlendWeights[1] &&
        Point -> BlendWeights[1] <= 0.2)
        (*PropertyAux) = 'b';
    else if (0.2 < Point -> BlendWeights[1] &&
	     Point -> BlendWeights[1]  <= 0.4)
        (*PropertyAux) = 'c';
    else if (0.4 < Point -> BlendWeights[1] &&
	     Point -> BlendWeights[1] <= 0.6)
        (*PropertyAux) = 'd';
    else if (0.6 < Point -> BlendWeights[1] &&
	     Point -> BlendWeights[1] <= 0.8)
        (*PropertyAux) = 'e';
    else if (0.8 < Point -> BlendWeights[1] &&
	     Point -> BlendWeights[1] < 1)
        (*PropertyAux) = 'f';

    *Ptr = (*PropertyAux);

    if (PropertyAux)
        IritFree(PropertyAux);

    return OutputPtr;
}

/************************* Blending for vectors. ****************************/
#ifdef IRIT_NOT_USED
static void *Vector3BlendingRule(
			      IritVMdlBlendPrepStruct *BlendCoordPrep,
			      const IrtRType Pt[3],
			      void *OutputPtr,
			      IritVMdlBlendCoordinatesFuncType Coordinates,
			      IritVMdlBlendPropertyFuncType PropertyFunction)
{
    const size_t
        SizeTopology = BlendCoordPrep -> NumVPrims;
    int inTriv, i,
        ObjIndex = -1;  /* If pt belongs to only one VMdl, this is its idx. */
    IrtRType
        StandardValue = 255.0;      /* What to assign in case nothing else. */
    IritVMdlBlendTopologyStruct *TopolStructIter,
        *TopolStructStart = BlendCoordPrep -> TopolStruct;
    IritVMdlBlendPointStruct
        *Point = BlendCoordPrep -> Point;
    IrtRType
    (*PropertyAux)[3] = NULL,
        (*Ptr)[3] = (IrtRType(*)[3]) OutputPtr;
    IrtRType
        *Vector = *Ptr;

    if (!Coordinates)
        Coordinates = IritVMdlBlendComputeDistCoordinates;
    if (!PropertyFunction)
        PropertyFunction = IritVMdlBlendDfltRGBPropertyFunction;
    if (!Pt) {
        for (i = 0; i < 3; ++i)
            Vector[i] = StandardValue;

        return OutputPtr;
    }

    inTriv = IritVMdlBlendInIntersection(BlendCoordPrep, &ObjIndex);

    /* P belongs to at least a TV. */
    if (inTriv == 0) {
        IRIT_ZAP_MEM(Point -> BlendWeights, SizeTopology * sizeof(CagdRType));
        Point -> BlendWeights[ObjIndex] = 1.0;

        PropertyFunction(BlendCoordPrep, ObjIndex, OutputPtr);

        return OutputPtr;
    }
    else if (inTriv < 0) {
        /* P is outside the V-model. */
        for (i = 0; i < 3; ++i)
            Vector[i] = StandardValue;

        return OutputPtr;
    }
    else {
        /* P belongs to a transition region. */
        if (Point -> BlendWeights[0] < 0)
            Coordinates(BlendCoordPrep, Pt);
        PropertyAux = (IrtRType(*)[3]) IritMalloc(sizeof(*PropertyAux));
    }

    /* Blending the values of the relevant properties. */
    IRIT_ZAP_MEM(Vector, 3 * sizeof(IrtRType));

    for (TopolStructIter = TopolStructStart;
        TopolStructIter != NULL;
        TopolStructIter = TopolStructIter -> Pnext) {
        size_t
            IndexObj = TopolStructIter -> Index;
        CagdRType
            BlendWeight = Point -> BlendWeights[IndexObj];
        IrtRType
            *VectorAux = *PropertyAux;

        if (!Point -> InclusionVPrims[IndexObj])
            continue;

        PropertyFunction(BlendCoordPrep, (int) IndexObj,
			 (void *) PropertyAux);

        for (i = 0; i < 3; ++i)
            Vector[i] += BlendWeight * VectorAux[i];
    }

    if (PropertyAux)
        IritFree(PropertyAux);

    return OutputPtr;
}
#endif /* IRIT_NOT_USED */

static void PrintMatrix(IritVMdlSlicerOutputImageStruct *Image,
                        const char *Output,
                        const int *Resolution)
{
    if (Output) {
        unsigned int
            Size = (unsigned int) (Resolution[0] + 1);
        int i, j;
        IritMiscImgPixelStruct
            *ImageLine = (IritMiscImgPixelStruct *)
		            IritMalloc(sizeof(IritMiscImgPixelStruct) * Size);
        IritMiscWriteGenInfoStructPtr GI;

        if ((GI = IritMiscImgWriteOpenFile(NULL, Output,
					   IRIT_MISC_IMAGE_UNKNOWN_TYPE,
					   FALSE, Resolution[0] + 1,
					   Resolution[1] + 1)) == NULL)
            return;

        for (i = 0; i <= Resolution[1]; ++i) {

            IRIT_ZAP_MEM(ImageLine, Size * sizeof(IritMiscImgPixelStruct));

            for (j = 0; j < (int) Size; ++j) {
                const int
                    k = j * (Resolution[1] + 1) + i;
                void
		    *Point = &((IritMiscImgPixelStruct *) Image -> Points)[k];
                IrtRType
                    (*Ptr)[3][3] = (IrtRType(*)[3][3]) Point,
                    (*Matrix)[3] = *Ptr,
                    Value1 = Matrix[0][0] + Matrix[0][1] + Matrix[0][2],
                    Value2 = Matrix[1][0] + Matrix[1][1] + Matrix[1][2],
                    Value3 = Matrix[2][0] + Matrix[2][1] + Matrix[2][2];

                ImageLine[j].r = Value1 > 255.0 ?
                    (IrtBType) 255.0 :
                    (Value1 < 0.0 ?
                        (IrtBType) 0.0 :
                        (IrtBType) (Value1));
                ImageLine[j].g = Value2 > 255.0 ?
                    (IrtBType) 255.0 :
                    (Value2 < 0.0 ?
                        (IrtBType) 0.0 :
                        (IrtBType) (Value2));
                ImageLine[j].b = Value3 > 255.0 ?
                    (IrtBType) 255.0 :
                    (Value3 < 0.0 ?
                        (IrtBType) 0.0 :
                        (IrtBType) (Value3));
            }

            IritMiscImgWritePutLine(GI, NULL, ImageLine);
        }
        IritMiscImgWriteCloseFile(GI);
        IritFree(ImageLine);
    }
    else
        fprintf(stderr, "Output is NULL. Nothing saved.\n");
}

static void PrintChar(IritVMdlSlicerOutputImageStruct *Image,
		      const char *Output,
		      const int *Resolution)
{
    if (Output) {
        unsigned int
            Size = (unsigned int) (Resolution[0] + 1);
        int i, j;
        IritMiscImgPixelStruct
            *ImageLine = (IritMiscImgPixelStruct *)
			    IritMalloc(sizeof(IritMiscImgPixelStruct) * Size);
        IritMiscWriteGenInfoStructPtr GI;

        if ((GI = IritMiscImgWriteOpenFile(NULL, Output,
					   IRIT_MISC_IMAGE_UNKNOWN_TYPE,
					   FALSE, Resolution[0] + 1,
					   Resolution[1] + 1)) == NULL)
            return;

        for (i = 0; i <= Resolution[1]; ++i) {
            IRIT_ZAP_MEM(ImageLine, Size * sizeof(IritMiscImgPixelStruct));

            for (j = 0; j < (int) Size; ++j) {
                const int
                    k = j * (Resolution[1] + 1) + i;
                void
                    *Point = &((IritMiscImgPixelStruct *) Image -> Points)[k];
                char
                    *Ptr = (char *) Point,
                    Char = (*Ptr);

                if (Char == 'a') {
                    ImageLine[j].r = (IrtBType) 255.0;
                    ImageLine[j].g = (IrtBType) 0.0;
                    ImageLine[j].b = (IrtBType) 100.0;
                }
                else if (Char == 'b') {
                    ImageLine[j].r = (IrtBType) (255.0 * 0.8 + 10.0 * 0.2);
                    ImageLine[j].g = (IrtBType) (0.2 * 255.0);
                    ImageLine[j].b = (IrtBType) (100.0 * 0.8 + 125.0 * 0.2);
                }
                else if (Char == 'c') {
                    ImageLine[j].r = (IrtBType) (255.0 * 0.6 + 10.0 * 0.4);
                    ImageLine[j].g = (IrtBType) (0.4 * 255.0);
                    ImageLine[j].b = (IrtBType) (100.0 * 0.6 + 125.0 * 0.4);
                }
                else if (Char == 'd') {
                    ImageLine[j].r = (IrtBType) (255.0 * 0.4 + 10.0 * 0.6);
                    ImageLine[j].g = (IrtBType) (0.6 * 255.0);
                    ImageLine[j].b = (IrtBType) (100.0 * 0.4 + 125.0 * 0.6);
                }
                else if (Char == 'e') {
                    ImageLine[j].r = (IrtBType) (255.0 * 0.2 + 10.0 * 0.8);
                    ImageLine[j].g = (IrtBType) (0.8 * 255.0);
                    ImageLine[j].b = (IrtBType) (100.0 * 0.2 + 125.0 * 0.8);
                }
                else if (Char == 'f') {
                    ImageLine[j].r = (IrtBType) (255.0 * 0.1 + 10.0 * 0.9);
                    ImageLine[j].g = (IrtBType) (0.9 * 255.0);
                    ImageLine[j].b = (IrtBType) (100.0 * 0.1 + 125.0 * 0.9);
                }
                else if (Char == 'g')  {
                    ImageLine[j].r = (IrtBType) 10.0;
                    ImageLine[j].g = (IrtBType) 255.0;
                    ImageLine[j].b = (IrtBType) 125.0;
                }
                else if (Char == 'x') {
                    ImageLine[j].r = (IrtBType) 255.0;
                    ImageLine[j].g = (IrtBType) 255.0;
                    ImageLine[j].b = (IrtBType) 255.0;
                }
            }

            IritMiscImgWritePutLine(GI, NULL, ImageLine);
        }
        IritMiscImgWriteCloseFile(GI);
        IritFree(ImageLine);
    }
    else
        fprintf(stderr, "Output is NULL. Nothing saved.\n");
}

#ifdef IRIT_NOT_USED
static IritMiscImgPixelStruct **ConvertVector2Img(
				       IritVMdlSlicerOutputImageStruct *Image,
				       const int *Resolution)
{
    int i, j;
    IritMiscImgPixelStruct
        **Output = (IritMiscImgPixelStruct **)
            IritMalloc(((size_t) Resolution[1] + 1) *
		                            sizeof(IritMiscImgPixelStruct *));

    IRIT_ZAP_MEM(Output,
        ((size_t) Resolution[1] + 1) * sizeof(IritMiscImgPixelStruct *));

    for (i = 0; i < Resolution[1] + 1; ++i) {
        Output[i] = (IritMiscImgPixelStruct *)
	      IritMalloc(((size_t) Resolution[0] + 1) *
					     sizeof(IritMiscImgPixelStruct));
        IRIT_ZAP_MEM(Output[i], ((size_t) Resolution[0] + 1) *
		     sizeof(IritMiscImgPixelStruct));
    }

    for (i = 0; i < Resolution[0] + 1; ++i) {
        for (j = 0; j < Resolution[1] + 1; ++j) {
            const int
                k = i * (Resolution[1] + 1) + j;
            void
                *Ptr = (char *) (Image -> Points) +
                k * sizeof(IrtRType[3]);
            IrtRType
                (*Point)[3] = (IrtRType(*)[3]) Ptr,
                Entry1 = (*Point)[0],
                Entry2 = (*Point)[1],
                Entry3 = (*Point)[2];

            Output[j][i].r = (IrtBType) Entry1;
            Output[j][i].g = (IrtBType) Entry2;
            Output[j][i].b = (IrtBType) Entry3;
        }
    }
    return Output;
}
#endif /* IRIT_NOT_USED */

int main(int argc, char **argv)
{
    int Handler,
        ResolutionMatrix[2] = { 100, 57 },
        ResolutionChar[2] = { 100, 57 };
    const IrtRType
        z = 0.1;
    IritVMdlBlendCoordinatesFuncType
        Coordinates = IritVMdlBlendComputeDistCoordinates;
    IritVMdlBlendPropertyFuncType
        PropertyFunctionMatrix = Matrix3x3PropertyFunction,
        PropertyFunctionChar = CharPropertyFunction;
    IritVMdlBlendRuleFuncType
        BlendingRuleMatrix = Matrix3x3BlendingRule,
        BlendingRuleChar = CharBlendingRule;

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
                    IritVMdlVModelStruct
                        *VMdl = PTmp -> U.VMdls;
                    IritVMdlBlendPrepStruct
		        *BlendCoordPrep = IritVMdlBlendCoordPrep(NULL, VMdl);

                    /* Blending matrices over a whole slice. */
                    IritVMdlSlicerOutputImageStruct
                        *OutputImageChar =
		            IritVMdlBlendPropertiesSlice(
						     BlendCoordPrep, z,
						     ResolutionChar,
						     Coordinates,
						     PropertyFunctionChar,
						     BlendingRuleChar,
						     sizeof(char)),
                        *OutputImageMatrix =
		            IritVMdlBlendPropertiesSlice(
						     BlendCoordPrep, z,
						     ResolutionMatrix,
						     Coordinates,
						     PropertyFunctionMatrix,
						     BlendingRuleMatrix,
						     sizeof(IrtRType[3][3]));

                    /* Convert matrices to a printable format and print it. */
                    PrintChar(OutputImageChar, "char.ppm", ResolutionChar);
                    PrintMatrix(OutputImageMatrix, "matrices.ppm",
				ResolutionMatrix);

                    IritVMdlSlicerOutputImageFree(OutputImageChar);
                    IritVMdlSlicerOutputImageFree(OutputImageMatrix);
                    IritVMdlBlendCoordPrepFree(BlendCoordPrep);
                }
            }

            IritPrsrPutObjectToFile3("/temp/TestOut.itd", PObj, 0);
        }
        else {
            fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
            exit(1);
        }
    }
    else {
        fprintf(stderr, "Usage: test VModel.itd\n");
        exit(2);
    }

    return 0;
}
