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
#include <fcntl.h>
#include <strings.h>

#include "zipper.h"

void left_proc();
void middle_proc();
void right_proc();
void quit_proc();
void refresh_drawing();

int edge_flag = 0;      /* draw polygons on the edge of the mesh? */
int bound_flag = 0;     /* draw edges on boundary of mesh? */
int pick_flag = 0;      /* pick a graphical primitive? */
int back_cull = 1;      /* cull backfacing polygons? */
int back_red = 0;               /* background color */
int back_grn = 0;
int back_blu = 0;
int copied_to_back = 0;         /* has the image been copied to backbuffer? */

/* This should be a variable; a function of hardware performance */
#define MAX_TRIANGLES_FOR_INTERACTION 5000

#define MOVE_WORLD -1
#define MOVE_LIGHT -2
int move_num = MOVE_WORLD;      /* move an object? */

/* "extra" lines to be displayed each frame */
static Vector* line1 = NULL;
static Vector* line2 = NULL;
static unsigned int* line_colors = NULL;
static int line_num;
static int line_max;


/******************************************************************************
Say if a triangle is front or back facing.

Entry:
  tri   - triangle to determine facing direction
  mat   - current transformation matrix
  timat - inverse transpose of current transformation

Exit:
  returns 1 if triangle is backfacing, 0 if it is frontfacing
******************************************************************************/

int backface_tri(tri, mat, timat)
Triangle* tri;
Matrix mat, timat;
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

mesh_to_world(sc, invec, outvec)
Scan* sc;
Vector invec, outvec;
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

world_to_mesh(sc, invec, outvec)
Scan* sc;
Vector invec, outvec;
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

world_to_mesh_normal(sc, in_norm, out_norm)
Scan* sc;
Vector in_norm, out_norm;
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

mesh_to_world_normal(sc, in_norm, out_norm)
Scan* sc;
Vector in_norm, out_norm;
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


/******************************************************************************
Initialize the data structures for "extra" lines to be displayed each frame.
******************************************************************************/

init_extra_lines()
{
    /* maybe free up old space */
    if (line1 != NULL) {
        free(line1);
        free(line2);
        free(line_colors);
        line1 = NULL;
        line2 = NULL;
        line_colors = NULL;
    }

    line_num = 0;
}


/******************************************************************************
Add a line to the "extra" lines to be displayed each frame.

Entry:
  p1,p2 - endpoints of line segment to be displayed
  color - line color
******************************************************************************/

add_extra_line(p1, p2, color)
Vector p1, p2;
unsigned int color;
{
    /* create new space for lines if necessary */
    if (line1 == NULL) {
        line_max = 200;
        line1 = (Vector*) malloc(sizeof(Vector) * line_max);
        line2 = (Vector*) malloc(sizeof(Vector) * line_max);
        line_colors = (unsigned int*) malloc(sizeof(unsigned int) * line_max);
    }

    /* reallocate more space if necessary */

    if (line_num == line_max) {
        line_max *= 2;
        line1 = (Vector*) realloc(line1, sizeof(Vector) * line_max);
        line2 = (Vector*) realloc(line2, sizeof(Vector) * line_max);
        line_colors = (unsigned int*)
                      realloc(line_colors, sizeof(unsigned int) * line_max);
    }

    /* add the new line to the list */
    vcopy(p1, line1[line_num]);
    vcopy(p2, line2[line_num]);
    line_colors[line_num] = color;
    line_num++;
}
