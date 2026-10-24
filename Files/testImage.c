/*****************************************************************************
* Filter the data.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute	of Technology		     *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/misc_lib.h"

void main(int argc, char **argv)
{
    int x = -1,
        y = -1;

    /* Do we have XY pixels to check> */
    if (argc > 2 && strcmp(argv[1], "-p") == 0) {
        assert (argc > 4);
	argc--;
	argv++;

	if (sscanf(argv[1], "%d", &x) != 1) {
	    fprintf(stderr,
		    "Expected an X image coordinate after \"-p\" option\n");
	    exit(2);
	}
	argc--;
	argv++;

	if (sscanf(argv[1], "%d", &y) != 1)
	    y = -1;
	else {
	    argc--;
	    argv++;
	}
    }

    if (argv[1] == NULL || (x < 0 && y < 0)) {
	fprintf(stderr, "Usage: testImage -p x { y } image_file\n");
	exit(2);
    }

    fprintf(stderr, "processing image \"%s\", examined pixel is (%d, %d)\n",
	    argv[1], x, y);

    if (argc == 2) {
        int MaxX, MaxY,
	    Alpha = TRUE;
        IritMiscImgRGBAPxlStruct *RGBAImage;
        IritMiscImgPixelStruct
	    *RGBImage = IritMiscImgReadImage(argv[1], &MaxX, &MaxY, &Alpha);

	if (RGBImage == NULL) {
	    fprintf(stderr, "Failed to open image file \"%s\"\n", argv[1]);
	    exit(1);
	}
	if (!Alpha) {
	    fprintf(stderr, "No alpha in the given image \"%s\"\n", argv[1]);
	    exit(4);
	}

	RGBAImage = (IritMiscImgRGBAPxlStruct *) RGBImage;
	RGBImage = NULL;

	/* Verify the desired coordinates of the pixels? */
	if (x >= 0 && x > MaxX) {
	    fprintf(stderr, "Input X value > image max X value (%d)\n", MaxX);
	    IritFree(RGBAImage);
	    exit(5);
	}
	if (y >= 0 && y > MaxY) {
	    fprintf(stderr, "Input Y value > image max Y value (%d)\n", MaxY);
	    IritFree(RGBAImage);
	    exit(6);
	}

	if (x >= 0 || y >= 0) {
	    /* Dump pixel values to stderr. */
	    if (x < 0) {
	        /* Dump an entire row. */
	        fprintf(stderr, "Pixel [   X,    Y]:    R   G   B    A\n");
	        for (x = 0; x <= MaxX; x++) {
		    IritMiscImgRGBAPxlStruct
		        *p = &RGBAImage[x + (MaxX - y) * (MaxX + 1)];

		    fprintf(stderr, "Pixel [%4d, %4d]:  %3d %3d %3d  %3d\n",
			    x, y, p -> r, p -> g, p -> b, p -> a);
		}
	    }
	    else if (y < 0) {
	        /* Dump the entire column. */
	        fprintf(stderr, "Pixel [   X,    Y]:   R   G   B   A\n");
	        for (y = 0; y <= MaxY; y++) {
		    IritMiscImgRGBAPxlStruct
		        *p = &RGBAImage[x + (MaxX - y) * (MaxX + 1)];

		    fprintf(stderr, "Pixel [%4d, %4d]:  %3d %3d %3d %3d\n",
			    x, y, p -> r, p -> g, p -> b, p -> a);
		}
	    }
	    else { /* Dump a single pixel. */
	        IritMiscImgRGBAPxlStruct
		    *p = &RGBAImage[x + (MaxX - y) * (MaxX + 1)];

		fprintf(stderr, "Pixel [   X,    Y]:   R   G   B   A\n");
		fprintf(stderr, "Pixel [%4d, %4d]:  %3d %3d %3d %3d\n",
			x, y, p -> r, p -> g, p -> b, p -> a);
	    }
	}

	IritFree(RGBAImage);
    }
    else {
	fprintf(stderr, "Usage: testImage { -p x { y } } image_file\n");
	exit(2);
    }

    exit(0);
}
