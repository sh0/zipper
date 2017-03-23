/*

Raw cyberware files.

---------------------------------------------------------------

Copyright (c) 1994 The Board of Trustees of The Leland Stanford
Junior University.  All rights reserved.

Permission to use, copy, modify and distribute this software and its
documentation for any purpose is hereby granted without fee, provided
that the above copyright notice and this permission notice appear in
all copies of this software and that you do not sell the software.

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

*/

#include <limits.h>

typedef struct RawData {
    int nlg;      /* number of range columns (samples in x) */
    int nlt;      /* number of rows (samples in y) */
    long lgincr;      /* distance between columns, in micrometers */
    long ltincr;      /* distance between rows, in micrometers */
    float* y;     /* vertical range positions */
    float* z;     /* depth range positions */
} RawData;

/* range data from a PLY file */
typedef struct RangeData {
    int nlg;      /* number of range columns (samples in x) */
    int nlt;      /* number of rows (samples in y) */
    int interlaced;       /* is it interlaced data? */
    int num_points;       /* number of range points */
    Vector* points;   /* the range points */
    float* confidence;    /* confidence */
    float* intensity;     /* intensity */
    unsigned char* red;   /* red part of color */
    unsigned char* grn;   /* green part of color */
    unsigned char* blu;   /* blue part of color */
    int* pnt_indices; /* row x col indices to the positions */
    int has_color;        /* color information? */
    int has_intensity;    /* intensity information? */
    int has_confidence;   /* confidence information? */
    int mult_confidence;  /* multiply confidence? */
    char obj_info[50][PATH_MAX];  /* Beware the arbitrary "50"!! */
    int num_obj_info;
} RangeData;

