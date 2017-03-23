/*

Matrix inversion module.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define  MAT_SIZE  4               /* size of matrix */
#define  DOUBLE_SIZE  8
#define  BOUND  0.0000001

typedef float Matrix [MAT_SIZE] [MAT_SIZE];
typedef float Double_mat [DOUBLE_SIZE] [MAT_SIZE];

/*
This matrix inversion routine is based on transforming the matrix we wish
to invert into the identity matrix.  To keep track of the operations that
we
use to transform the matrix, we perform those same operations on an identity
matrix.  We will make use of a double matrix:

[ m11 m21 m31 m41  1   0   0   0 ]
[ m12 m22 m32 m42  0   1   0   0 ]
[ m13 m23 m33 m43  0   0   1   0 ]
[ m14 m24 m34 m44  0   0   0   1 ]

The operations we may perform are these:

   1)  switch two rows
   2)  subtract a multiple of one row from another row
   3)  multiply each element of a row by a number

When we are done, the double matrix should look like this:

[ 1   0   0   0  i11 i21 i31 i41 ]
[ 0   1   0   0  i12 i22 i32 i42 ]
[ 0   0   1   0  i13 i23 i33 i43 ]
[ 0   0   0   1  i14 i24 i34 i44 ]

The right half of the double matrix is  the inverse.
*/


/******************************************************************************
Invert a 4 by 4 matrix.
Returns the determinant of the matrix.
******************************************************************************/

float mat_invert (mat)
Matrix mat;
{
  Double_mat dmat;
  int a,b;
  int best;
  float det = 1.0;

  invert_to_double (mat, dmat);

  /*  convert to upper triangular, with one's on diagonal  */

  for (a = 0; a < MAT_SIZE; a++) {

    /* find the row with the largest absolute value in column "a" */

    best = invert_best_lead (a, dmat);

    /* if no row has good leading value, the matrix  */
    /* is not invertible and the determinant is zero */

    if (best == -1)
      return (0.0);

    if (best != a)  {
      invert_switch_rows (a, best, dmat);
      det = -1.0 * det;
    }

    det = det * dmat [a] [a];
    invert_multiply_row (a, 1 / dmat [a] [a], dmat);

    for (b = a + 1; b < MAT_SIZE; b++)
      invert_subtract_rows (a, b, dmat [a] [b], dmat);
  }

  /*  now get rid of the upper triangle  */

  for (a = MAT_SIZE - 1; a >= 1; a--)
    for (b = a - 1; b >= 0; b--)
      invert_subtract_rows (a, b, dmat [a] [b], dmat);

  invert_from_double (mat, dmat);

  return (det);
}


/******************************************************************************
Place a matrix into the left side of a double matrix and place the identity
matrix on the right side.
******************************************************************************/

invert_to_double (mat,dmat)
Matrix mat;
Double_mat dmat;
{
  int i,j;

  for (i = 0; i < MAT_SIZE; i ++)
    for (j = 0; j < MAT_SIZE; j++) {
      dmat [i] [j] = mat [i] [j];
      dmat [i + MAT_SIZE] [j] = 0.0;
    }

  for (j = 0; j < MAT_SIZE; j++)
    dmat [j + MAT_SIZE] [j] = 1.0;
}


/******************************************************************************
Returns the right half of a double matrix.
******************************************************************************/

invert_from_double (mat,dmat)
Matrix mat;
Double_mat dmat;
{
  int i,j;

  for (i = 0; i < MAT_SIZE; i++)
    for (j = 0; j < MAT_SIZE; j++)
      mat [i] [j] = dmat [i + MAT_SIZE] [j];
}


/******************************************************************************
Switch two rows of a double matrix.
******************************************************************************/

invert_switch_rows (row1,row2,dmat)
int row1,row2;
Double_mat dmat;
{
  int i;
  float t;

  for (i = 0; i < DOUBLE_SIZE; i++) {
    t = dmat [i] [row1];
    dmat [i] [row1] = dmat [i] [row2];
    dmat [i] [row2] = t;
  }
}


/******************************************************************************
Subtract a multiple of row 1 from row 2.
******************************************************************************/

invert_subtract_rows (row1,row2,m,dmat)
int row1,row2;
float m;
Double_mat dmat;
{
  int i;

  for (i = 0; i < DOUBLE_SIZE; i++)
    dmat [i] [row2] = dmat [i] [row2] - m * dmat [i] [row1];
}


/******************************************************************************
Multiply each element of a row by a number.
******************************************************************************/

invert_multiply_row (row,m,dmat)
int row;
float m;
Double_mat dmat;
{
  int i;

  for (i = 0; i < DOUBLE_SIZE; i++)
    dmat [i] [row] = m * dmat [i] [row];
}


/******************************************************************************
Returns the number of the row whose leading value is largest in absolute
value.

Entry:
  row  - look from this row on down
  dmat - matrix to look in

Exit:
  returns best row number, or -1 if there is no good row
******************************************************************************/

invert_best_lead(row,dmat)
int row;
Double_mat dmat;
{
  float max;
  int j,pos;

  max = -1.0;

  for (j = row; j < MAT_SIZE; j++)
    if (fabs (dmat [row] [j]) > max) {
      max = fabs (dmat [row] [j]);
      pos = j;
    }

  if (max < BOUND)
    return (-1);
  else
    return (pos);
}


/******************************************************************************
Print out a double matrix. For debugging purposes.
******************************************************************************/

static print_double (dmat)
Double_mat dmat;
{
  int i,j;

  printf ("\n");

  for (j = 0; j < MAT_SIZE; j++) {
    for (i = 0; i < DOUBLE_SIZE; i++)
      printf ("%f  ", dmat [i] [j]);
    printf ("\n");
  }

  printf ("\n");
}
