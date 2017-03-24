/*

Draw various objects.

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

#include "zipper.h"

int backface_tri(Triangle* tri, Matrix mat, Matrix timat);
void mesh_to_world(Scan* sc, Vector invec, Vector outvec);
void world_to_mesh(Scan* sc, Vector invec, Vector outvec);
void world_to_mesh_normal(Scan* sc, Vector in_norm, Vector out_norm);
void mesh_to_world_normal(Scan* sc, Vector in_norm, Vector out_norm);

/******************************************************************************
Say if a triangle is front or back facing.

Entry:
  tri   - triangle to determine facing direction
  mat   - current transformation matrix
  timat - inverse transpose of current transformation

Exit:
  returns 1 if triangle is backfacing, 0 if it is frontfacing
******************************************************************************/
int backface_tri(Triangle* tri, Matrix mat, Matrix timat)
{
    Vector norm;
    Vector pnt;
    Vertex** v = tri->verts;
    float dot;

    /* transform the normal of the triangle */
    norm[X] = -tri->aa;
    norm[Y] = -tri->bb;
    norm[Z] = -tri->cc;
    mat_apply(timat, norm);


    if (0) { // is_orthographic
        if (norm[Z] < 0)
            return (1);
        else
            return (0);
    };


    /* transform the centriod of the triangle */
    pnt[X] = (v[0]->coord[X] + v[1]->coord[X] + v[2]->coord[X]) * 0.333333;
    pnt[Y] = (v[0]->coord[Y] + v[1]->coord[Y] + v[2]->coord[Y]) * 0.333333;
    pnt[Z] = (v[0]->coord[Z] + v[1]->coord[Z] + v[2]->coord[Z]) * 0.333333;
    mat_apply(mat, pnt);

    dot = vdot(norm, pnt);

    if (dot > 0)
        return (1);
    else
        return (0);
}

/******************************************************************************
Transform point from meshes space to world space.

Entry:
  sc    - scan of mesh
  invec - point to transform

Exit:
  outvec - transformed point
******************************************************************************/
void mesh_to_world(Scan* sc, Vector invec, Vector outvec)
{
    int i;
    Vector v;

    for (i = 0; i < 3; i++)
        v[i] = invec[X] * sc->rotmat[X][i] + invec[Y] * sc->rotmat[Y][i] +
               invec[Z] * sc->rotmat[Z][i];

    outvec[X] = v[X] + sc->xtrans;
    outvec[Y] = v[Y] + sc->ytrans;
    outvec[Z] = v[Z] + sc->ztrans;
}

/******************************************************************************
Transform point from world space to a meshes space.

Entry:
  sc    - scan of mesh
  invec - point to transform

Exit:
  outvec - transformed point
******************************************************************************/
void world_to_mesh(Scan* sc, Vector invec, Vector outvec)
{
    int i;
    Vector v;

    v[X] = invec[X] - sc->xtrans;
    v[Y] = invec[Y] - sc->ytrans;
    v[Z] = invec[Z] - sc->ztrans;

    for (i = 0; i < 3; i++)
        outvec[i] = v[X] * sc->rotmat[i][X] + v[Y] * sc->rotmat[i][Y] +
                    v[Z] * sc->rotmat[i][Z];
}

/******************************************************************************
Transform normal from world space to a meshes space.

Entry:
  sc      - scan of mesh
  in_norm - surface normal to transform

Exit:
  out_norm - transformed point
******************************************************************************/
void world_to_mesh_normal(Scan* sc, Vector in_norm, Vector out_norm)
{
    int i;
    Vector v;

    v[X] = in_norm[X];
    v[Y] = in_norm[Y];
    v[Z] = in_norm[Z];

    for (i = 0; i < 3; i++)
        out_norm[i] = v[X] * sc->rotmat[i][X] + v[Y] * sc->rotmat[i][Y] +
                      v[Z] * sc->rotmat[i][Z];
}

/******************************************************************************
Transform normal from meshes space to world space.

Entry:
  sc      - scan of mesh
  in_norm - point to transform

Exit:
  out_norm - transformed point
******************************************************************************/
void mesh_to_world_normal(Scan* sc, Vector in_norm, Vector out_norm)
{
    int i;
    Vector v;

    for (i = 0; i < 3; i++)
        v[i] = in_norm[X] * sc->rotmat[X][i] + in_norm[Y] * sc->rotmat[Y][i] +
               in_norm[Z] * sc->rotmat[Z][i];

    out_norm[X] = v[X];
    out_norm[Y] = v[Y];
    out_norm[Z] = v[Z];
}
