/*

Vector routines.

Greg Turk - written a long time ago

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
#include <matrix.h>
#include <math.h>

extern Matrix mat_form;                /* from the matrix routines */


/******************************************************************************
Adds two vectors.
c = a + b
******************************************************************************/

my_vadd (a,b,c)
Vector a,b,c;
{
  c [0] = a [0] + b [0];
  c [1] = a [1] + b [1];
  c [2] = a [2] + b [2];
}


/******************************************************************************
Subtract two vectors.
c = a - b
******************************************************************************/

my_vsub (a,b,c)
Vector a,b,c;
{
  c [0] = a [0] - b [0];
  c [1] = a [1] - b [1];
  c [2] = a [2] - b [2];
}


/******************************************************************************
Cross product of two vectors.
c = a cross b
******************************************************************************/

my_vcross (c,a,b)
  Vector a,b,c;
{
  float x,y;

  x    = a[Y] * b[Z] - a[Z] * b[Y];
  y    = a[Z] * b[X] - a[X] * b[Z];
  c[Z] = a[X] * b[Y] - a[Y] * b[X];

  c[X] = x;
  c[Y] = y;
}


/******************************************************************************
Normalize a vector.
******************************************************************************/

float vnorm(v)
  Vector v;
{
  float len;

  len = sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (len == 0) {
    /*
    printf ("vnorm: zero length vector\n");
    */
    return (0.0);
  }

  v[0] /= len;
  v[1] /= len;
  v[2] /= len;

  return (len);
}


/******************************************************************************
Returns dot product of two vectors.
******************************************************************************/

float my_vdot (a,b)
Vector a,b;
{
  return (a [0] * b [0] + a [1] * b [1] + a [2] * b [2]);
}


/******************************************************************************
Scale a vector.
******************************************************************************/

my_vscale (a,scale)
Vector a;
float scale;
{
  a [0] = a [0] * scale;
  a [1] = a [1] * scale;
  a [2] = a [2] * scale;
}


/******************************************************************************
Return the length of a vector.
******************************************************************************/

float vlen (a)
Vector a;
{
  return (sqrt (a [0] * a [0] + a [1] * a [1] + a [2] * a [2]));
}


/******************************************************************************
Create a vector from three real numbers.
******************************************************************************/

vcompose (x,y,z,a)
float x,y,z;
Vector a;
{
  a [0] = x;
  a [1] = y;
  a [2] = z;
}


/******************************************************************************
Copy a vector.

Entry:
  source - vector to copy from

Exit:
  dest - vector to copy to
******************************************************************************/

my_vcopy (source,dest)
Vector source,dest;
{
  dest [0] = source [0];
  dest [1] = source [1];
  dest [2] = source [2];
}


/******************************************************************************
Print out a vector.
******************************************************************************/

print_vector (a)
Vector a;
{
  printf ("%f  %f  %f\n", a [0], a [1], a [2]);
}


/******************************************************************************
Multiply a vector by a transformation matrix.  Assumes the last row of the
matrix is [0 0 0 1].

b = m * a
******************************************************************************/

vapply (m,a,b)
Matrix m;
Vector a,b;
{
  int j;
  Vector t;

  for (j = 0; j <= 2; j++)
    t [j] = a [0] * m [0] [j] + a [1] * m [1] [j] +
            a [2] * m [2] [j] + m [3] [j];

  for (j = 0; j <= 2; j++)
    b [j] = t [j];
}


/******************************************************************************
Transform a point by the current transformation matrix.

Entry:
  x,y,z - point to transform

Exit:
  x,y,z - transformed point
******************************************************************************/

transform_point (x,y,z)
float *x,*y,*z;
{
  float xx,yy,zz;

  xx = *x * mat_form [0] [0] + *y * mat_form [1] [0] +
       *z * mat_form [2] [0] + mat_form [3] [0];

  yy = *x * mat_form [0] [1] + *y * mat_form [1] [1] +
       *z * mat_form [2] [1] + mat_form [3] [1];

  zz = *x * mat_form [0] [2] + *y * mat_form [1] [2] +
       *z * mat_form [2] [2] + mat_form [3] [2];

  *x = xx;
  *y = yy;
  *z = zz;
}


/******************************************************************************
Transform a plane equation by the current transformation matrix.

Entry:
  a,b,c,d - planar coefficients

Exit:
  a,b,c,d - transformed coefficients
******************************************************************************/

transform_plane (a,b,c,d)
float *a,*b,*c,*d;
{
  float aa,bb,cc,dd;

  aa = *a * mat_form [0] [0] + *b * mat_form [1] [0] +
       *c * mat_form [2] [0] + mat_form [3] [0];

  bb = *a * mat_form [0] [1] + *b * mat_form [1] [1] +
       *c * mat_form [2] [1] + mat_form [3] [1];

  cc = *a * mat_form [0] [2] + *b * mat_form [1] [2] +
       *c * mat_form [2] [2] + mat_form [3] [2];

  dd = *a * mat_form [0] [3] + *b * mat_form [1] [3] +
       *c * mat_form [2] [3] + mat_form [3] [3];

  *a = aa;
  *b = bb;
  *c = cc;
  *d = dd;
}


/******************************************************************************
Rotate a vector onto a plane.
Axes 1, 2, and 3 correspond to x, y, and z.

Entry:
  v      - vector to rotate
  around - axis to rotate around
  onto   - with 'around' defines positive half plane to rotate 'v' onto

Exit:
  v - rotated vector
  m - matrix that performs rotation
******************************************************************************/

rotate_vector_to_plane (around,onto,v,m)
int around,onto;
Vector v;
Matrix m;
{
  float len;
  float a,b;
  int p1,p2;
  Matrix mat;

  p1 = onto - 1;

  /*  figure out the third axis  */

  if (onto + around == 3)
    p2 = 2;
  else if (onto + around == 5)
    p2 = 0;
  else
    p2 = 1;

  a = v [p1];
  b = v [p2];

  if ((a != 0) || (b != 0)) {

    len = sqrt (a*a + b*b);

    a = a / len;
    b = b / len;

    mat_ident (mat);

    mat [p1] [p1] =  a;   mat [p2] [p1] = b;
    mat [p1] [p2] = -b;   mat [p2] [p2] = a;

    vapply (mat, v, v);
    mat_mult (m, mat, m);
  }
}


/******************************************************************************
Alters the current transformation matrix to represent a specified viewpoint.
Should be made to check if the 'up' vector is paralell to the line joining
the 'from' and 'at' vectors.

Entry:
  ax,ay,az - what is being looked at; placed along positive z-axis
  fx,fy,fz - the eye position; will be placed at (0,0,0)
  ux,uy,uz - which direction is up; placed on positive y part of yz-plane

Exit:
  modifies the current transformation matrix 'mat_form'
******************************************************************************/

#if 0
view (ax,ay,az,fx,fy,fz,ux,uy,uz)
float ax,ay,az,fx,fy,fz,ux,uy,uz;
{
  Vector from,at,up;
  Vector v;
  Matrix mat;
  int i,j;

  vcompose (ax, ay, az, at);
  vcompose (fx, fy, fz, from);
  vcompose (ux, uy, uz, up);

  push();
  mat_ident (mat_form);

  translate (-from [0], -from [1], -from [2]);
  vapply (mat_form, at, v);

  rotate_vector_to_plane (3, 2, v, mat_form);
  rotate_vector_to_plane (1, 3, v, mat_form);

  vadd (from, up, v);
  vapply (mat_form, v, v);
  rotate_vector_to_plane (3, 2, v, mat_form);

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      mat [i] [j] = mat_form [i] [j];

  pop();
  mat_mult (mat_form, mat_form, mat);
}
#endif
