/*

Matrix module.

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

#include "matrix.h"

// Old vector operations
void vcopy(const float* v1, float* v2)
{
    int i;
    for (i = 0 ; i < 3 ; i++)
        v2[i] = v1[i];
}

void vset(float* v, float x, float y, float z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

float vlength(const float* v)
{
    return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void vscale(float* v, float div)
{
    v[0] *= div;
    v[1] *= div;
    v[2] *= div;
}

void vadd(const float* src1, const float* src2, float* dst)
{
    dst[0] = src1[0] + src2[0];
    dst[1] = src1[1] + src2[1];
    dst[2] = src1[2] + src2[2];
}

void vsub(const float* src1, const float* src2, float* dst)
{
    dst[0] = src1[0] - src2[0];
    dst[1] = src1[1] - src2[1];
    dst[2] = src1[2] - src2[2];
}

float vdot(const float* v1, const float* v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void vcross(const float* v1, const float* v2, float* cross)
{
    float temp[3];
    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    vcopy(temp, cross);
}

// Vector operations
float vnorm(Vector v)
{
    float len;
    len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len == 0)
        return 0.0f;

    v[0] /= len;
    v[1] /= len;
    v[2] /= len;
    return len;
}

float vlen(Vector a)
{
    return (sqrt(a [0] * a [0] + a [1] * a [1] + a [2] * a [2]));
}

void vapply(Matrix m, Vector a, Vector b)
{
    int j;
    Vector t;

    for (j = 0; j <= 2; j++)
        t [j] = a [0] * m [0] [j] + a [1] * m [1] [j] +
                a [2] * m [2] [j] + m [3] [j];
    for (j = 0; j <= 2; j++)
        b [j] = t [j];
}

void vsub2(Vector a, Vector b, Vector c)
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

// Matrix operations
void mat_apply(Matrix m, Vector v)
{
    int j;
    Vector t;
    for (j = 0; j <= 2; j++)
        t [j] = v[0] * m[0][j] + v[1] * m[1][j] + v[2] * m[2][j] + m[3][j];
    v[0] = t[0];
    v[1] = t[1];
    v[2] = t[2];
}

void mat_ident(Matrix m)
{
    int i;
    for (i = 0; i <= 3; i++) {
        m[i][0] = 0.0;
        m[i][1] = 0.0;
        m[i][2] = 0.0;
        m[i][3] = 0.0;
        m[i][i] = 1.0;
    }
}

void mat_translate(Matrix m, float x, float y, float z)
{
    mat_ident(m);
    m[3][0] = x;
    m[3][1] = y;
    m[3][2] = z;
}

void mat_mult(Matrix prod, Matrix a, Matrix b)
{
    int i, j;
    Matrix m;

    for (i = 0; i <= 3; i++) {
        for (j = 0; j <= 3; j++) {
            m [i] [j] = a [0] [j] * b [i] [0] + a [1] [j] * b [i] [1] +
                        a [2] [j] * b [i] [2] + a [3] [j] * b [i] [3];
        }
    }
    for (i = 0; i <= 3; i++) {
        for (j = 0; j <= 3; j++)
            prod [i] [j] = m [i] [j];
    }
}

void mat_copy(Matrix dest, Matrix source)
{
    int i, j;
    for (i = 0; i <= 3; i++) {
        for (j = 0; j <= 3; j++)
            dest[i][j] = source[i][j];
    }
}
