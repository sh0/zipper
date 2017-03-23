/*

Matrix and Vector header.

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

#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef float Matrix[4][4];
typedef float Vector[3];
typedef float Vector4[4];
typedef float Plane[4];
typedef float Quaternion[4];

#define X 0
#define Y 1
#define Z 2
#define W 3

extern float vnorm();
extern float vdot();
extern float vlen();
extern void vset(float*, float, float, float);
extern void vscale(float*, float);

#endif /* __MATRIX_H__ */

