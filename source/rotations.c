/*

Routines for rotation of objects.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include "strings.h"
#include "malloc.h"
#include "cyfile.h"
#include "zipper.h"


/******************************************************************************
Build the rotation matrix for an object.

Entry:
  sc - object to build matrix for
******************************************************************************/

build_rotmat(sc)
Scan* sc;
{
    int i, j, k;
    float theta;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            sc->rotmat[i][j] = (i == j);

    theta = - sc->rotate * M_PI / 180.0;

    sc->rotmat[X][X] = cos(theta);
    sc->rotmat[X][Z] = sin(theta);
    sc->rotmat[Z][X] = -sin(theta);
    sc->rotmat[Z][Z] = cos(theta);
}

