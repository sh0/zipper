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

#ifndef ZIPPER_MATRIX_H
#define ZIPPER_MATRIX_H

// External
#include <math.h>

// Types
typedef float Matrix[4][4];
typedef float Vector[3];
typedef float Vector4[4];
typedef float Plane[4];
typedef float Quaternion[4];

// Constants
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
