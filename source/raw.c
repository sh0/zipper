/*
 * Read data from raw cyberware files.
 *
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

// External
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Intermal
#include "raw.h"

/******************************************************************************
Read geometry from a raw file.

Entry:
  name - name of file to read from

Exit:
  returns pointer to data, or NULL if it couldn't read from file
******************************************************************************/
RawData* read_raw_geom(char* name)
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

    rawdata = (RawData*) malloc(sizeof(RawData));

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

    rawdata->y = (float*) malloc(sizeof(float) * nlg * nlt);
    rawdata->z = (float*) malloc(sizeof(float) * nlg * nlt);

    for (i = 0; i < nlg; i++) {
        fread(&rawdata->z[i * nlt], sizeof(float) * nlt, 1, fp);
        fread(&rawdata->y[i * nlt], sizeof(float) * nlt, 1, fp);
    }

    fclose(fp);

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
int get_raw_coord(Scan* scan, int lt, int lg, Vector vec)
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
