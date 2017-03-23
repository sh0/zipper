/*

Read data from raw cyberware files.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "malloc.h"
#include "cyfile.h"
#include "zipper.h"
#include "matrix.h"
#include "raw.h"


extern int global_chew_count;


/******************************************************************************
Read geometry from a raw file.

Entry:
  name - name of file to read from

Exit:
  returns pointer to data, or NULL if it couldn't read from file
******************************************************************************/

RawData* read_raw_geom(name)
char* name;
{
    int i;
    FILE* fp;
    char filename[80];
    RawData* rawdata;
    int nlg;
    int nlt;

    strcpy(filename, name);
    if (strlen(filename) < 4 ||
        strcmp(filename + strlen(filename) - 4, ".raw") != 0)
        strcat(filename, ".raw");

    fp = fopen(filename, "r");

    if (fp == NULL) {
        fprintf(stderr, "Can't read file %s\n", filename);
        return (NULL);
    }

    if (fp == NULL) return NULL;

    rawdata = (RawData*) myalloc(sizeof(RawData));

    fread(&nlg, sizeof(int), 1, fp);
    fread(&nlt, sizeof(int), 1, fp);
    fread(&rawdata->lgincr, sizeof(long), 1, fp);
    fread(&rawdata->ltincr, sizeof(long), 1, fp);
    rawdata->nlg = nlg;
    rawdata->nlt = nlt;

    /*
    printf ("nlg nlt: %d %d\n", nlg, nlt);
    printf ("lgincr ltincr: %d %d\n", rawdata->lgincr, rawdata->ltincr);
    */

    rawdata->y = (float*) myalloc(sizeof(float) * nlg * nlt);
    rawdata->z = (float*) myalloc(sizeof(float) * nlg * nlt);

    for (i = 0; i < nlg; i++) {
        fread(&rawdata->z[i * nlt], sizeof(float) * nlt, 1, fp);
        fread(&rawdata->y[i * nlt], sizeof(float) * nlt, 1, fp);
    }

    fclose(fp);

#if 0
    /* maybe chew back the edges of a mesh */
    if (global_chew_count)
        chew_raw_edges(rawdata);
#endif

    return (rawdata);
}

#define REALVOID (-HUGE)
#define IS_REALVOID(x) ((x) < REALVOID/1000.0)

/******************************************************************************
Return a range point from the raw geometry, given (i,j) indices.

Entry:
  scan  - scan that holds the range data
  lt,lg - indices to fetch the data from

Exit:
  vec - the returned range position
  returns 0 if data point exists, 1 if it is void
******************************************************************************/

get_raw_coord(scan, lt, lg, vec)
Scan* scan;
int lt, lg;
Vector vec;
{
    RawData* rawdata = scan->raw_geom;
    float to_meters = 1.0e-6;

    vec[X] = to_meters * rawdata->lgincr * (lg - rawdata->nlg / 2);
    if (lt % 2 == 1)
        vec[X] += to_meters * rawdata->lgincr * 0.5;

    vec[Z] = rawdata->y[lg * rawdata->nlt + lt];
    vec[Y] = rawdata->z[lg * rawdata->nlt + lt];

    vec[Y] -= to_meters * rawdata->ltincr * (rawdata->nlt / 2);

    if (IS_REALVOID(vec[Z]))
        return (1);

    return (0);
}


/******************************************************************************
Chew back the edges of a raw range image.
******************************************************************************/

chew_raw_edges(rawdata)
RawData* rawdata;
{
    int i, k;
    int lt, lg;
    int nlt, nlg;
    unsigned char* chew;
    float z, nz;
    float tolerance;
    int index;
    int count;

    nlt = rawdata->nlt;
    nlg = rawdata->nlg;

    /* 1 cm tolerance for jumps */
    tolerance = 0.01;

    chew = (unsigned char*) myalloc(sizeof(unsigned char) * nlt * nlg);

    /* chew back the edge of the range image a number of times */
    for (k = 0; k < global_chew_count; k++) {

        count = 0;

        /* zero out the chew flags */
        for (i = 0; i < nlt * nlg; i++)
            chew[i] = 0;

        /* find where to chew */

        for (lt = 1; lt < nlt - 1; lt++)
            for (lg = 1; lg < nlg - 1; lg++) {

                index = lg * nlt + lt;
                z = rawdata->y[index];

                if (IS_REALVOID(z))
                    continue;

                nz = rawdata->y[index + 1];
                if (z - nz > tolerance) chew[index] = 1;
                nz = rawdata->y[index - 1];
                if (z - nz > tolerance) chew[index] = 1;
                nz = rawdata->y[index + nlt];
                if (z - nz > tolerance) chew[index] = 1;
                nz = rawdata->y[index - nlt];
                if (z - nz > tolerance) chew[index] = 1;
            }

        /* make those chewed places into void places */
        for (lt = 1; lt < nlt - 1; lt++)
            for (lg = 1; lg < nlg - 1; lg++) {
                index = lg * nlt + lt;
                if (chew[index]) {
                    rawdata->y[index] = REALVOID;
                    count++;
                }
            }

        printf("%d chewed\n", count);

    }

    free(chew);
}

