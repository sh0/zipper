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

#ifndef ZIPPER_MATRIX_H
#define ZIPPER_MATRIX_H

#include <math.h>

typedef float Matrix[4][4];
typedef float Vector[3];
typedef float Vector4[4];
typedef float Plane[4];
typedef float Quaternion[4];

#define X 0
#define Y 1
#define Z 2
#define W 3

// Old vector operations
void vcopy(const float*, float*);
void vset(float*, float, float, float);
float vlength(const float*);
void vscale(float*, float);
void vadd(const float*, const float*, float*);
void vsub(const float*, const float*, float*);
float vdot(const float*, const float*);
void vcross(const float*, const float*, float*);

// Vector operations
float vnorm(Vector v);
float vlen(Vector a);
void vapply(Matrix m, Vector a, Vector b);
void vsub2(Vector a, Vector b, Vector c);

// Matrix operations
void mat_apply(Matrix m, Vector v);
void mat_ident(Matrix m);
void mat_translate(Matrix m, float x, float y, float z);
void mat_mult(Matrix prod, Matrix a, Matrix b);
void mat_copy(Matrix dest, Matrix source);

#endif
