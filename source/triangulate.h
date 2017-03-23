/*

Header defining a polygon type.

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

#define BORDER       1
#define INTERIOR     2
#define CREATED      3
#define OLD_CREATED  4

#define NPMAX 40

typedef struct Npoly {
    int nverts;           /* number of vertices */
    unsigned char border[NPMAX];  /* if vertex is part of border polygon */
    int index[NPMAX];     /* index of vertex */
    float x[NPMAX];
    float y[NPMAX];
    int which_tri;        /* index of which triangle it belongs to */
} Npoly;


