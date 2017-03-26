/*
 * Copyright (c) 1995-2017, Stanford University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of Stanford University nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY STANFORD UNIVERSITY ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL STANFORD UNIVERSITY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ZIPPER_RAW_H
#define ZIPPER_RAW_H

// Internal
#include "zipper.h"
#include "matrix.h"

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

// Declarations
RawData* read_raw_geom(char* name);
int get_raw_coord(Scan* scan, int lt, int lg, Vector vec);

#endif
