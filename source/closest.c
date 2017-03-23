/*

Find translation and rotation to best match up pairs of points.

This method of registering point sets is from the paper:

  "A Method for Registration of 3-D Shapes"
  Paul J. Besl and Neil D. McKay
  IEEE Transactions on Pattern Analysis and Machine Intelligence
  Vol. 14, No. 2, February 1992

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
#include <event.h>
#include <matrix.h>
#include <matrix1.h>
#include <matrix2.h>
#include <myalloc.h>

typedef float Matrix3[3][3];

typedef struct Pair {
  Vector p1,p2;
  float weight;
} Pair;

static Pair **pairs = NULL;
static int npairs = 0;
static int max_pairs;


/******************************************************************************
Initialize the list of pairs.
******************************************************************************/

init_pairs()
{
  int i;

  if (pairs != NULL) {
    for (i = 0; i < npairs; i++)
      free (pairs[i]);
    free (pairs);
  }

  npairs = 0;

  max_pairs = 5000;
  pairs = (Pair **) myalloc (sizeof (Pair *) * max_pairs);
}


/******************************************************************************
Add a pair to the list of pairs.

Entry:
  p1,p2  - pair of corresponding points
  weight - how much this pair should be weighted with respect to other pairs
******************************************************************************/

add_pair(p1,p2,weight)
  Vector p1,p2;
  float weight;
{
  Pair *new;

  if (npairs >= max_pairs) {
    max_pairs *= 2;
    pairs = (Pair **) realloc (pairs, sizeof (Pair *) * max_pairs);
  }

  new = (Pair *) myalloc (sizeof (Pair));
  pairs[npairs] = new;
  npairs++;

  new->p1[X] = p1[X];
  new->p1[Y] = p1[Y];
  new->p1[Z] = p1[Z];

  new->p2[X] = p2[X];
  new->p2[Y] = p2[Y];
  new->p2[Z] = p2[Z];

  new->weight = weight;
}


/******************************************************************************
Find the best rotation and translation to match two collections of points.
The set of all first elements from the pairs is one collection of points,
and the second elements forms the other collection.

Exit:
  rotmat   - matrix with the best rotation and translation for the match
  quat_rot - quaternion form of the rotation
  trans    - translational component
******************************************************************************/

match_pairs(rotmat,quat_rot,trans)
  Matrix rotmat;
  Quaternion quat_rot;
  Vector trans;
{
  int i,j,k;
  Vector c1,c2;		/* center of mass for two point collections */
  Vector v1,v2;
  float recip;
  float tr;
  MAT *m, *q;
  VEC *v;
  Matrix3 cov;
  Matrix3 aij;
  Quaternion quat;
  float weight;
  float weight_sum;

  /* find the center of mass for the two collections */

  weight_sum = 0;
  c1[X] = c1[Y] = c1[Z] = 0;
  c2[X] = c2[Y] = c2[Z] = 0;
  for (i = 0; i < npairs; i++) {

    weight = pairs[i]->weight;
    weight_sum += weight;

    c1[X] += pairs[i]->p1[X] * weight;
    c1[Y] += pairs[i]->p1[Y] * weight;
    c1[Z] += pairs[i]->p1[Z] * weight;

    c2[X] += pairs[i]->p2[X] * weight;
    c2[Y] += pairs[i]->p2[Y] * weight;
    c2[Z] += pairs[i]->p2[Z] * weight;
  }

  /*
  recip = 1.0 / npairs;
  */
  recip = 1.0 / weight_sum;

  c1[X] *= recip;
  c1[Y] *= recip;
  c1[Z] *= recip;

  c2[X] *= recip;
  c2[Y] *= recip;
  c2[Z] *= recip;

  /* create the cross-covariance matrix */

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      cov[i][j] = 0;

  for (k = 0; k < npairs; k++) {

    weight = pairs[k]->weight;

    v1[X] = (pairs[k]->p1[X] - c1[X]) * weight;
    v1[Y] = (pairs[k]->p1[Y] - c1[Y]) * weight;
    v1[Z] = (pairs[k]->p1[Z] - c1[Z]) * weight;

    v2[X] = pairs[k]->p2[X] - c2[X];
    v2[Y] = pairs[k]->p2[Y] - c2[Y];
    v2[Z] = pairs[k]->p2[Z] - c2[Z];

    cov[X][X] += v1[X] * v2[X];
    cov[X][Y] += v1[X] * v2[Y];
    cov[X][Z] += v1[X] * v2[Z];

    cov[Y][X] += v1[Y] * v2[X];
    cov[Y][Y] += v1[Y] * v2[Y];
    cov[Y][Z] += v1[Y] * v2[Z];

    cov[Z][X] += v1[Z] * v2[X];
    cov[Z][Y] += v1[Z] * v2[Y];
    cov[Z][Z] += v1[Z] * v2[Z];
  }

#if 0
  /* I don't think this is necessary */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      cov[i][j] *= recip;
#endif

  /* aij = cov - transpose(cov) */

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      aij[i][j] = cov[i][j] - cov[j][i];

  /* find the trace of the covariance matrix */

  tr = cov[X][X] + cov[Y][Y] + cov[Z][Z];

  /* create the necessary 4x4 matrix */

  m = m_get(4,4);
  q = m_get(4,4);
  v = v_get(4);

  m->me[0][0] = tr;

  m->me[1][0] = m->me[0][1] = aij[1][2];
  m->me[2][0] = m->me[0][2] = aij[2][0];
  m->me[3][0] = m->me[0][3] = aij[0][1];

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      m->me[i+1][j+1] = cov[i][j] + cov[j][i] - (i == j) * tr;

  /* find the eigenvector corresponding to the largest eigenvalue */
  /* of this matrix */

  v = symmeig( m, q, v );

  if( v->ve[0] > v->ve[1] ){
      if( v->ve[0] > v->ve[2] ){
	  if( v->ve[0] > v->ve[3] )
	    i = 0;
	  else
	    i = 3;
      }
      else{
	  if( v->ve[2] > v->ve[3] )
	    i = 2;
	  else
	    i =3;
      }
  }
  else{
      if( v->ve[1] > v->ve[2] ){
	  if( v->ve[1] > v->ve[3] )
	    i = 1;
	  else
	    i = 3;
      }
      else{
	  if( v->ve[2] > v->ve[3] )
	    i = 2;
	  else
	    i =3;
      }
  }

  quat[0] = q->me[0][i];
  quat[1] = q->me[1][i];
  quat[2] = q->me[2][i];
  quat[3] = q->me[3][i];

  convert_quat_to_mat (quat, rotmat);

  /* determine best translation */
  vapply (rotmat, c2, c2);
  rotmat[3][0] = c1[0] - c2[0];
  rotmat[3][1] = c1[1] - c2[1];
  rotmat[3][2] = c1[2] - c2[2];

  /* return rotation and translation info */
  quat_rot[0] = quat[0];
  quat_rot[1] = quat[1];
  quat_rot[2] = quat[2];
  quat_rot[3] = quat[3];
  trans[X] = rotmat[3][0];
  trans[Y] = rotmat[3][1];
  trans[Z] = rotmat[3][2];
}


/******************************************************************************
Convert a quaternion into a rotation matrix.

Entry:
  q - quaternion = (x,y,z,w)

Exit:
  mat - rotation matrix
******************************************************************************/

convert_quat_to_mat(q,mat)
  Quaternion q;
  Matrix mat;
{
 float q00,q01,q02,q03;
 float q11,q12,q13;
 float q22,q23;
 float q33;

 q00 = q[0] * q[0];
 q01 = q[0] * q[1];
 q02 = q[0] * q[2];
 q03 = q[0] * q[3];

 q11 = q[1] * q[1];
 q12 = q[1] * q[2];
 q13 = q[1] * q[3];

 q22 = q[2] * q[2];
 q23 = q[2] * q[3];

 q33 = q[3] * q[3];

 mat[X][X] = q00 + q11 - q22 - q33;
 mat[X][Y] = 2 * (q12 - q03);
 mat[X][Z] = 2 * (q13 + q02);
 mat[X][W] = 0;

 mat[Y][X] =  2 * (q12 + q03);
 mat[Y][Y] =  q00 + q22 - q11 - q33;
 mat[Y][Z] = 2 * (q23 - q01);
 mat[Y][W] = 0;

 mat[Z][X] = 2 * (q13 - q02);
 mat[Z][Y] = 2 * (q23 + q01);
 mat[Z][Z] = q00 + q33 - q11 - q22;
 mat[Z][W] = 0;

 mat[W][X] = 0;
 mat[W][Y] = 0;
 mat[W][Z] = 0;
 mat[W][W] = 1;
}

