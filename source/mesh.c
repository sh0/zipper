/*

Mesh building routines.

Greg Turk, September 1993

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
#include <assert.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>

#include "zipper.h"
#include "raw.h"

static int CONF_EDGE_ZERO;
static float CONF_EDGE_COUNT_FACTOR;
static float CONF_ANGLE;
static float CONF_EXPONENT;


set_conf_edge_count_factor(factor)
float factor;
{
    CONF_EDGE_COUNT_FACTOR = factor;
}


float
get_conf_angle()
{
    return CONF_ANGLE;
}


set_conf_angle(factor)
float factor;
{
    CONF_ANGLE = factor;
}


float
get_conf_exponent()
{
    return CONF_EXPONENT;
}


set_conf_exponent(factor)
float factor;
{
    CONF_EXPONENT = factor;
}


float
get_conf_edge_count_factor()
{
    return CONF_EDGE_COUNT_FACTOR;
}


set_conf_edge_zero(set)
int set;
{
    CONF_EDGE_ZERO = set != 0;
}


int
get_conf_edge_zero()
{
    return CONF_EDGE_ZERO;
}


/******************************************************************************
Place into a scan a newly-created triangle mesh of a given vertex spacing.

Entry:
  sc    - scan data to make into triangle mesh
  level - level of mesh (how detailed), in [0,1,2,3]
******************************************************************************/

create_scan_mesh(sc, level)
Scan* sc;
int level;
{
    Mesh* mesh;
    Mesh* make_mesh_raw();
    Mesh* make_mesh_ply();

    /* return if we've already created this mesh */
    if (sc->meshes[level] != NULL)
        return;

    /* make sure all higher-level (less-detailed) meshes are created */
    if (level < 3)
        create_scan_mesh(sc, level + 1);

    /* create the mesh */

    if (sc->file_type == CYFILE)
        assert(0);
        //mesh = make_mesh(sc, level, TABLE_DIST * level_to_inc(level));
    else if (sc->file_type == RAWFILE)
        mesh = make_mesh_raw(sc, level, TABLE_DIST * level_to_inc(level));
    else if (sc->file_type == PLYRANGEFILE)
        mesh = make_mesh_ply(sc, level, TABLE_DIST * level_to_inc(level));
    else {
        printf("create_scan_mesh: wrong file type\n");
        return;
    }

    sc->meshes[level] = mesh;

    printf("%s (level %d): %d verts, %d tris.\n",
           sc->name, level + 1, mesh->nverts, mesh->ntris);
}


/******************************************************************************
Create a triangle mesh from scan data.

Entry:
  sc         - scan data to make into triangle mesh
  level      - level of mesh (how detailed), in [0,1,2,3]
  table_dist - distance for spatial subdivision hash table

Exit:
  returns pointer to newly-created mesh
******************************************************************************/
#if 0
Mesh* make_mesh(sc, level, table_dist)
Scan* sc;
int level;
float table_dist;
{
    int i, j;
    int ii, jj;
    int a, b;
    int inc;
    Vector vec;
    GSPEC* gs = sc->gs;
    int* vert_index;
    int result;
    int index;
    int in1, in2, in3, in4;
    Vertex* vt1, *vt2, *vt3, *vt4;
    int count;
    float max_length;
    int max_lt, max_lg;
    int rotational_flag;
    extern float edge_length_max();
    Mesh* mesh;
    Triangle* make_triangle();
    Triangle* t1, *t2;

    /* decide if this is rotational or translational scan */

    if (gs->flags & FLAG_CARTESIAN)
        rotational_flag = 0;
    else
        rotational_flag = 1;

    /* pick how far apart the mesh samples are, based on the level of */
    /* detail requested */

    inc = level_to_inc(level);

    /* get maximum okay length of a triangle edge */

    max_length = edge_length_max(level);

    /* create mesh */

    mesh = (Mesh*) malloc(sizeof(Mesh));

    /* allocate space for new triangles and vertices */

    mesh->nverts = 0;
    max_lt = (gs->nlt - 1) / inc + 1;
    max_lg = (gs->nlg - 1) / inc + 1;
    mesh->max_verts = max_lt * max_lg;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    mesh->ntris = 0;
    mesh->max_tris = mesh->max_verts * 2;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 200;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 200;
    mesh->parent_scan = sc;

    /* create table saying whether a vertex is okay and where it is */
    /* in the vertex list */

    vert_index = (int*) malloc(sizeof(int) * max_lt * max_lg);
    for (i = 0; i < max_lt * max_lg; i++)
        vert_index[i] = -1;

    /* create the vertices */

    for (i = 0, a = 0; i < gs->nlt; i += inc, a++) {
        for (j = 0, b = 0; j < gs->nlg; j += inc, b++) {
            result = get_gs_coord(sc, i, j, vec);
            if (result == 0) {        /* coordinate okay */
                index = make_vertex(mesh, vec);
                vert_index[a + b * max_lt] = index;
                /* tuck longitude into confidence for vertex error finding later */
                mesh->verts[index]->confidence = j;
            } else {          /* coordinate void */
                vert_index[a + b * max_lt] = -1;
            }
        }
    }

    /* create the triangles */

    for (i = 0; i < max_lt - 1; i ++)
        for (j = 0; j < max_lg - 1 + rotational_flag; j ++) {

            ii = (i + 1) % max_lt;
            jj = (j + 1) % max_lg;

            /* count the number of good vertices */
            in1 = vert_index[ i +  j * max_lt];
            in2 = vert_index[ii +  j * max_lt];
            in3 = vert_index[ii + jj * max_lt];
            in4 = vert_index[ i + jj * max_lt];
            count = (in1 >= 0) + (in2 >= 0) + (in3 >= 0) + (in4 >= 0);
            if (in1 >= 0)
                vt1 = mesh->verts[in1];
            if (in2 >= 0)
                vt2 = mesh->verts[in2];
            if (in3 >= 0)
                vt3 = mesh->verts[in3];
            if (in4 >= 0)
                vt4 = mesh->verts[in4];

            t1 = t2 = NULL;

            if (count == 4) {     /* all 4 vertices okay, so make 2 tris */
                float len1, len2;
                Vector v1, v2;

                /* compute lengths of cross-edges */
                vsub(vt1->coord, vt3->coord, v1);
                vsub(vt2->coord, vt4->coord, v2);
                len1 = vlen(v1);
                len2 = vlen(v2);

                /* make triangles that minimize the cross-length */
                if (len1 < len2) {
                    t1 = make_triangle(mesh, vt1, vt2, vt3, max_length);
                    t2 = make_triangle(mesh, vt1, vt3, vt4, max_length);
                } else {
                    t1 = make_triangle(mesh, vt2, vt3, vt4, max_length);
                    t2 = make_triangle(mesh, vt2, vt4, vt1, max_length);
                }
            } else if (count == 3) {  /* only 3 vertices okay, so make 1 tri */
                if (in1 == -1) {
                    t1 = make_triangle(mesh, vt2, vt3, vt4, max_length);
                } else if (in2 == -1) {
                    t1 = make_triangle(mesh, vt1, vt3, vt4, max_length);
                } else if (in3 == -1) {
                    t1 = make_triangle(mesh, vt1, vt2, vt4, max_length);
                } else { /* in4 == -1 */
                    t1 = make_triangle(mesh, vt1, vt2, vt3, max_length);
                }
            }

            /* see if triangles are backwards */
            if ((count == 3 || count == 4) && rotational_flag) {
                float z1, z2;

                if (in1 && in3) {
                    z1 = GETR(gs, i * inc, j * inc);
                    z2 = GETR(gs, ii * inc, jj * inc);
                } else {
                    z1 = GETR(gs, ii * inc, j * inc);
                    z2 = GETR(gs, i * inc, jj * inc);
                }

                if (z1 < 0 && z2 < 0) {
                    if (t1) {
                        /*
                        t1->normal[X] *= -1;
                        t1->normal[Y] *= -1;
                        t1->normal[Z] *= -1;
                        */
                        t1->aa *= -1;
                        t1->bb *= -1;
                        t1->cc *= -1;
                        t1->dd *= -1;
                    }
                    if (t2) {
                        /*
                        t2->normal[X] *= -1;
                        t2->normal[Y] *= -1;
                        t2->normal[Z] *= -1;
                        */
                        t2->aa *= -1;
                        t2->bb *= -1;
                        t2->cc *= -1;
                        t2->dd *= -1;
                    }
                }

            }

        }

    /* compute vertex normals */
    find_vertex_normals(mesh);

    /* compute certainty about vertex positions */
    vertex_errors(mesh, sc, rotational_flag);

    /* free the space used by the index table */
    free(vert_index);

    /* initialize hash table for vertices in mesh */
    init_table(mesh, table_dist);

    /* find the edges of the mesh */
    find_mesh_edges(mesh);

    /* return pointer to new mesh */
    return (mesh);
}
#endif

/******************************************************************************
Create a triangle mesh from scan data.

Entry:
  sc         - scan data to make into triangle mesh
  level      - level of mesh (how detailed), in [0,1,2,3]
  table_dist - distance for spatial subdivision hash table

Exit:
  returns pointer to newly-created mesh
******************************************************************************/

Mesh* make_mesh_raw(sc, level, table_dist)
Scan* sc;
int level;
float table_dist;
{
    int i, j;
    int ii, jj;
    int a, b;
    int inc;
    Vector vec;
    RawData* rawdata = sc->raw_geom;
    int* vert_index;
    int result;
    int index;
    int in1, in2, in3, in4;
    Vertex* vt1, *vt2, *vt3, *vt4;
    int count;
    float max_length;
    int max_lt, max_lg;
    extern float edge_length_max();
    Mesh* mesh;

    /* pick how far apart the mesh samples are, based on the level of */
    /* detail requested */

    inc = level_to_inc(level);

    /* get maximum okay length of a triangle edge */

    max_length = edge_length_max(level);

    /* create mesh */

    mesh = (Mesh*) malloc(sizeof(Mesh));

    /* allocate space for new triangles and vertices */

    mesh->nverts = 0;
    max_lt = (rawdata->nlt - 1) / inc + 1;
    max_lg = (rawdata->nlg - 1) / inc + 1;
    mesh->max_verts = max_lt * max_lg;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    mesh->ntris = 0;
    mesh->max_tris = mesh->max_verts * 2;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 200;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 200;
    mesh->parent_scan = sc;

    /* create table saying whether a vertex is okay and where it is */
    /* in the vertex list */

    vert_index = (int*) malloc(sizeof(int) * max_lt * max_lg);
    for (i = 0; i < max_lt * max_lg; i++)
        vert_index[i] = -1;

    /* create the vertices */

    for (i = 0, a = 0; i < rawdata->nlt; i += inc, a++) {
        for (j = 0, b = 0; j < rawdata->nlg; j += inc, b++) {
            result = get_raw_coord(sc, i, j, vec);
            if (result == 0) {        /* coordinate okay */
                index = make_vertex(mesh, vec);
                vert_index[a + b * max_lt] = index;
                /* tuck longitude into confidence for vertex error finding later */
                mesh->verts[index]->confidence = j;
            } else {          /* coordinate void */
                vert_index[a + b * max_lt] = -1;
            }
        }
    }

    /* create the triangles */

    for (i = 0; i < max_lt - 1; i ++)
        for (j = 0; j < max_lg - 1; j ++) {

            ii = (i + 1) % max_lt;
            jj = (j + 1) % max_lg;

            /* count the number of good vertices */
            in1 = vert_index[ i +  j * max_lt];
            in2 = vert_index[ii +  j * max_lt];
            in3 = vert_index[ii + jj * max_lt];
            in4 = vert_index[ i + jj * max_lt];
            count = (in1 >= 0) + (in2 >= 0) + (in3 >= 0) + (in4 >= 0);
            if (in1 >= 0)
                vt1 = mesh->verts[in1];
            if (in2 >= 0)
                vt2 = mesh->verts[in2];
            if (in3 >= 0)
                vt3 = mesh->verts[in3];
            if (in4 >= 0)
                vt4 = mesh->verts[in4];

            if (count == 4) {     /* all 4 vertices okay, so make 2 tris */
                float len1, len2;
                Vector v1, v2;

                /* compute lengths of cross-edges */
                vsub(vt1->coord, vt3->coord, v1);
                vsub(vt2->coord, vt4->coord, v2);
                len1 = vlen(v1);
                len2 = vlen(v2);

                /* make triangles that minimize the cross-length */
                if (len1 < len2) {
                    make_triangle(mesh, vt1, vt2, vt3, max_length);
                    make_triangle(mesh, vt1, vt3, vt4, max_length);
                } else {
                    make_triangle(mesh, vt2, vt3, vt4, max_length);
                    make_triangle(mesh, vt2, vt4, vt1, max_length);
                }
            } else if (count == 3) {  /* only 3 vertices okay, so make 1 tri */
                if (in1 == -1) {
                    make_triangle(mesh, vt2, vt3, vt4, max_length);
                } else if (in2 == -1) {
                    make_triangle(mesh, vt1, vt3, vt4, max_length);
                } else if (in3 == -1) {
                    make_triangle(mesh, vt1, vt2, vt4, max_length);
                } else { /* in4 == -1 */
                    make_triangle(mesh, vt1, vt2, vt3, max_length);
                }
            }

        }

    /* compute vertex normals */
    find_vertex_normals(mesh);

    /* compute certainty about vertex positions */
    vertex_errors(mesh, sc, 0);

    /* free the space used by the index table */
    free(vert_index);

    /* initialize hash table for vertices in mesh */
    init_table(mesh, table_dist);

    /* find the edges of the mesh */
    find_mesh_edges(mesh);

    /* lower the confidence at edges */
    lower_edge_confidence(mesh, level);

    /* return pointer to new mesh */
    return (mesh);
}


/******************************************************************************
Create a triangle mesh from scan data.

Entry:
  sc         - scan data to make into triangle mesh
  level      - level of mesh (how detailed), in [0,1,2,3]
  table_dist - distance for spatial subdivision hash table

Exit:
  returns pointer to newly-created mesh
******************************************************************************/

Mesh* make_mesh_ply(sc, level, table_dist)
Scan* sc;
int level;
float table_dist;
{
    int i, j;
    int ii, jj;
    int a, b;
    int inc;
    Vector vec;
    RangeData* plydata = sc->ply_geom;
    int result;
    int index;
    int in1, in2, in3, in4;
    Vertex* vt1, *vt2, *vt3, *vt4;
    int count;
    float max_length;
    int max_lt, max_lg;
    extern float edge_length_max();
    Mesh* mesh;
    int nlt, nlg;
    int* vert_index;

    /* pick how far apart the mesh samples are, based on the level of */
    /* detail requested */

    inc = level_to_inc(level);

    /* get maximum okay length of a triangle edge */

    max_length = edge_length_max(level);

    /* create mesh */

    mesh = (Mesh*) malloc(sizeof(Mesh));

    /* allocate space for new triangles and vertices */

    mesh->nverts = 0;
    max_lt = (plydata->nlt - 1) / inc + 1;
    max_lg = (plydata->nlg - 1) / inc + 1;
    nlt = plydata->nlt;
    nlg = plydata->nlg;
    mesh->max_verts = max_lt * max_lg;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    mesh->ntris = 0;
    mesh->max_tris = mesh->max_verts * 2;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 200;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 200;
    mesh->parent_scan = sc;

    /* create a list saying whether a vertex is going to be used */

    vert_index = (int*) malloc(sizeof(int) * plydata->num_points);
    for (i = 0; i < plydata->num_points; i++)
        vert_index[i] = -1;

    /* see which vertices will be used in a triangle */

    for (i = 0; i <= nlt - inc; i += inc)
        for (j = 0; j <= nlg - inc; j += inc) {
            in1 = plydata->pnt_indices[i + j * nlt];
            if (in1 >= 0)
                vert_index[in1] = 1;
        }

    /* create the vertices */

    for (i = 0; i < plydata->num_points; i++) {
        if (vert_index[i] == -1)
            continue;
        vec[X] = plydata->points[i][X];
        vec[Y] = plydata->points[i][Y];
        vec[Z] = plydata->points[i][Z];
        index = make_vertex(mesh, vec);
        vert_index[i] = index;
        mesh->verts[index]->confidence = plydata->confidence[i];
        mesh->verts[index]->intensity = plydata->intensity[i];
        mesh->verts[index]->red = plydata->red[i];
        mesh->verts[index]->grn = plydata->grn[i];
        mesh->verts[index]->blu = plydata->blu[i];
    }

    /* create the triangles */

    for (i = 0; i < nlt - inc; i += inc)
        for (j = 0; j < nlg - inc; j += inc) {

            ii = (i + inc) % nlt;
            jj = (j + inc) % nlg;

            /* count the number of good vertices */
            in1 = plydata->pnt_indices[ i +  j * nlt];
            in2 = plydata->pnt_indices[ i + jj * nlt];
            in3 = plydata->pnt_indices[ii + jj * nlt];
            in4 = plydata->pnt_indices[ii +  j * nlt];
            count = (in1 >= 0) + (in2 >= 0) + (in3 >= 0) + (in4 >= 0);
            if (in1 >= 0)
                vt1 = mesh->verts[vert_index[in1]];
            if (in2 >= 0)
                vt2 = mesh->verts[vert_index[in2]];
            if (in3 >= 0)
                vt3 = mesh->verts[vert_index[in3]];
            if (in4 >= 0)
                vt4 = mesh->verts[vert_index[in4]];

            if (count == 4) {     /* all 4 vertices okay, so make 2 tris */
                float len1, len2;
                Vector v1, v2;

                /* compute lengths of cross-edges */
                vsub(vt1->coord, vt3->coord, v1);
                vsub(vt2->coord, vt4->coord, v2);
                len1 = vlen(v1);
                len2 = vlen(v2);

                /* make triangles that minimize the cross-length */
                if (len1 < len2) {
                    make_triangle(mesh, vt1, vt2, vt3, max_length);
                    make_triangle(mesh, vt1, vt3, vt4, max_length);
                } else {
                    make_triangle(mesh, vt2, vt3, vt4, max_length);
                    make_triangle(mesh, vt2, vt4, vt1, max_length);
                }
            } else if (count == 3) {  /* only 3 vertices okay, so make 1 tri */
                if (in1 == -1) {
                    make_triangle(mesh, vt2, vt3, vt4, max_length);
                } else if (in2 == -1) {
                    make_triangle(mesh, vt1, vt3, vt4, max_length);
                } else if (in3 == -1) {
                    make_triangle(mesh, vt1, vt2, vt4, max_length);
                } else { /* in4 == -1 */
                    make_triangle(mesh, vt1, vt2, vt3, max_length);
                }
            }

        }

    /* free up the vertex index list */
    free(vert_index);

    /* compute vertex normals */
    find_vertex_normals(mesh);

    /* compute certainty about vertex positions */
    if (!plydata->has_confidence)
        vertex_errors(mesh, sc, 0, 0);
    else if (plydata->mult_confidence)
        vertex_errors(mesh, sc, 0, 1);


    /* initialize hash table for vertices in mesh */
    /*
      printf("\nSkipping hash table initialization!!!\n");
    */
    init_table(mesh, table_dist);

    /* find the edges of the mesh */
    find_mesh_edges(mesh);

    /* lower the confidence at edges */
#if 0
    if (!plydata->has_confidence)
#endif
        lower_edge_confidence(mesh, level);

    /* return pointer to new mesh */
    return (mesh);
}


/******************************************************************************
Compute how much error should be associated with the vertices of a mesh.

Entry:
  mesh     - mesh holding the vertices
  scan     - scan containing mesh
  rot_flag - mesh from a rotational scan? 1 = rotational scan, 0 = linear scan
  mult     - multiply confidence by current value?
******************************************************************************/

vertex_errors(mesh, scan, rot_flag, mult)
Mesh* mesh;
Scan* scan;
int rot_flag, mult;
{
    int i, j;
    float val;
    Vertex* vert;
    int lg;
    float c, s;

    if (rot_flag)
        /* rotational scans */
        for (i = 0; i < mesh->nverts; i++) {

            vert = mesh->verts[i];

            /* we stashed the longitude in confidence for now */
            lg = (int) vert->confidence;

            /* (s, 0, c) is laser direction */
            c = scan->cos_theta[lg] * 1.0e6;
            s = scan->sin_theta[lg] * 1.0e6;

            /* dot product of laser direction and surface normal */
            val = s * vert->normal[X] + c * vert->normal[Z];

            /* clamp the value to [0,1] */
            if (val < 0)
                val = 0;
            else if (val > 1)
                val = 1;

            if (mult)
                vert->confidence *= val;
            else
                vert->confidence = val;
        } else
        /* for linear scans */
        for (i = 0; i < mesh->nverts; i++) {

            vert = mesh->verts[i];

            /* dot product of laser direction and surface normal */
            /* (laser direction is (0, 0, 1)) */
            val = vert->normal[Z];

            val = vert->normal[X] * sin(CONF_ANGLE * M_PI / 180) +
                  vert->normal[Z] * cos(CONF_ANGLE * M_PI / 180);

            if (CONF_EXPONENT != 1.0)
                val = pow(fabs(val), CONF_EXPONENT);

            /* clamp the value to [0,1] */
            if (val < 0)
                val = 0;
            else if (val > 1)
                val = 1;

            if (mult)
                vert->confidence *= val;
            else
                vert->confidence = val;
        }
}


/******************************************************************************
Lower the confidence value on edges.

Entry:
  mesh  - mesh on which to lower the edge confidence
  level - level of mesh detail
******************************************************************************/

lower_edge_confidence(mesh, level)
Mesh* mesh;
int level;
{
    int i, j, k;
    int pass;
    int val;
    float recip;
    Vertex* v;
    int chew_count;

    switch (level) {
        case 0:
            chew_count = (int)(8 * CONF_EDGE_COUNT_FACTOR + 0.5);
            break;
        case 1:
            chew_count = (int)(4 * CONF_EDGE_COUNT_FACTOR + 0.5);
            break;
        case 2:
            chew_count = (int)(2 * CONF_EDGE_COUNT_FACTOR + 0.5);
            break;
        case 3:
            chew_count = (int)(1 * CONF_EDGE_COUNT_FACTOR + 0.5);
            break;
        default:
            fprintf(stderr, "lower_edge_confidence: bad switch %d\n", level);
            exit(-1);
    }

    /* make several passes through the vertices */
    for (pass = 1; pass < chew_count; pass++) {

        /* propagate higher on-edge values away from edges */
        for (i = 0; i < mesh->nverts; i++) {

            v = mesh->verts[i];
            if (v->on_edge != 0)
                continue;

            for (j = 0; j < v->nverts; j++)
                if (v->verts[j]->on_edge == pass) {
                    v->on_edge = pass + 1;
                    break;
                }
        }

    }

    /* lower the confidences on the edge */

    if (CONF_EDGE_ZERO) {
        recip = 1.0 / (chew_count);

        for (i = 0; i < mesh->nverts; i++) {
            v = mesh->verts[i];
            val = v->on_edge;
            if (val) {
                v->confidence *= (val - 1) * recip;
                if (val > 1)
                    v->on_edge = 0;
            }
        }
    } else {
        recip = 1.0 / (chew_count + 1);

        for (i = 0; i < mesh->nverts; i++) {
            v = mesh->verts[i];
            val = v->on_edge;
            if (val) {
                v->confidence *= val * recip;
                if (val > 1)
                    v->on_edge = 0;
            }
        }
    }
}


/******************************************************************************
Free up all the memory used in a mesh.

Entry:
  mesh - mesh to clear out
******************************************************************************/

clear_mesh(mesh)
Mesh* mesh;
{
    int i;
    Vertex* v;
    Triangle* t;

    /* free the hash table */
    if (mesh->table->verts != NULL)
        free(mesh->table->verts);

    /* free the vertices */
    if (mesh->verts != NULL) {
        for (i = 0; i < mesh->nverts; i++) {
            v = mesh->verts[i];
            free(v->verts);
            free(v->tris);
            if (v->edges)
                free(v->edges);
            free(v);
        }
        free(mesh->verts);
    }

    /* free the triangles */
    if (mesh->tris != NULL) {
        for (i = 0; i < mesh->ntris; i++) {
            t = mesh->tris[i];
            if (t->more)
                free(t->more);
            free(t);
        }
        free(mesh->tris);
    }

    /* free the edges */
    if (mesh->edges != NULL) {
        for (i = 0; i < mesh->nedges; i++) {
            free(mesh->edges[i]);
        }
        free(mesh->edges);
    }

    mesh->ntris = 0;
    mesh->nverts = 0;
    mesh->nedges = 0;
}


/******************************************************************************
Create a new vertex and add it to the list of vertices in a mesh.

Entry:
  mesh - mesh to add to
  vec  - coordinate of new vertex

Exit:
  returns index within list to new vertex
******************************************************************************/

int make_vertex(mesh, vec)
Mesh* mesh;
Vector vec;
{
    Vertex* vert;

    /* maybe make room for more vertices */

    if (mesh->nverts == mesh->max_verts) {
        mesh->max_verts = (int)(mesh->max_verts * 1.5);
        mesh->verts = (Vertex**)
                      realloc(mesh->verts, sizeof(Vertex*) * mesh->max_verts);
    }

    /* create new vertex and add it to the list */

    vert = (Vertex*) malloc(sizeof(Vertex));
    vert->coord[X] = vec[X];
    vert->coord[Y] = vec[Y];
    vert->coord[Z] = vec[Z];

    vert->normal[X] = 0;
    vert->normal[Y] = 0;
    vert->normal[Z] = 0;

    vert->ntris = 0;
    vert->max_tris = 8;
    vert->index = mesh->nverts;
    vert->moving = 0;
    vert->cinfo = NULL;
    vert->old_mesh = mesh;
    vert->confidence = 0;
    vert->tris = (Triangle**) malloc(sizeof(Triangle*) * vert->max_tris);

    vert->nverts = 0;
    vert->max_verts = 8;
    vert->verts = (Vertex**) malloc(sizeof(Vertex*) * vert->max_verts);

    vert->nedges = 0;
    vert->max_edges = 0;
    vert->edges = NULL;

    mesh->verts[mesh->nverts] = vert;
    mesh->nverts++;

    /* return index to the new vertex */

    return (mesh->nverts - 1);
}


/******************************************************************************
Create a new triangle and add it to the list in a mesh.

Entry:
  mesh        - mesh to add to
  vt1,vt2,vt3 - vertices to make triangle from
  max_len     - maximum allowed length of a triangle edge

Exit:
  returns pointer to newly-created triangle, or NULL if triangle was too big
******************************************************************************/

Triangle* make_triangle(mesh, vt1, vt2, vt3, max_len)
Mesh* mesh;
Vertex* vt1, *vt2, *vt3;
float max_len;
{
    int i, j;
    Triangle* tri;
    Vector v1, v2, v3;
    Vector cross;

    /* check edge lengths of triangle */
    vsub(vt3, vt1, v1);
    vsub(vt2, vt1, v2);
    vsub(vt3, vt2, v3);
    /* return if any edge is too long */
    if (vlen(v1) > max_len || vlen(v2) > max_len || vlen(v3) > max_len)
        return (NULL);

    /* maybe make room for more triangles */
    if (mesh->ntris >= mesh->max_tris) {
        mesh->max_tris = (int)(mesh->max_tris * 1.5);
        mesh->tris = (Triangle**)
                     realloc(mesh->tris, sizeof(Triangle*) * mesh->max_tris);
    }

    /* compute surface normal to the triangle */
    vcross(v1, v2, cross);
    vnorm(cross);

    /* create new triangle and add it to the list */

    tri = (Triangle*) malloc(sizeof(Triangle));
    tri->verts[0] = vt1;
    tri->verts[1] = vt2;
    tri->verts[2] = vt3;
    tri->mark = 0;
    tri->eat_mark = 0;
    tri->index = mesh->ntris;
    tri->more = NULL;
    tri->clips = NULL;
    tri->dont_touch = 0;
    if (set_triangle_geometry(tri) == 1) {
        free(tri);
        return NULL;
    }
    mesh->tris[mesh->ntris] = tri;
    mesh->ntris++;

    /* add this new triangle to each of its vertices lists */

    add_tri_to_vert(vt1, tri);
    add_tri_to_vert(vt2, tri);
    add_tri_to_vert(vt3, tri);

    /* add vertices to other vertices' lists */
    for (i = 0; i < 3; i++) {

        int found1 = 0;
        int found2 = 0;
        Vertex* vert = tri->verts[i];

        /* see if either of the other two vertices are already in vert's list */
        for (j = 0; j < vert->nverts; j++) {
            if (vert->verts[j] == tri->verts[(i + 1) % 3])
                found1 = 1;
            if (vert->verts[j] == tri->verts[(i + 2) % 3])
                found2 = 1;
        }

        /* add the ones not already on the list to vert's list */

        if (!found1) {
            /* maybe allocate more room for vertex list */
            if (vert->nverts >= vert->max_verts) {
                vert->max_verts += 4;
                vert->verts =
                    (Vertex**) realloc(vert->verts, sizeof(Vertex*) * vert->max_verts);
            }
            vert->verts[vert->nverts++] = tri->verts[(i + 1) % 3];
        }

        if (!found2) {
            /* maybe allocate more room for vertex list */
            if (vert->nverts >= vert->max_verts) {
                vert->max_verts += 4;
                vert->verts =
                    (Vertex**) realloc(vert->verts, sizeof(Vertex*) * vert->max_verts);
            }
            vert->verts[vert->nverts++] = tri->verts[(i + 2) % 3];
        }
    }

    /* return pointer to new triangle */
    return (tri);
}


/******************************************************************************
Remove a triangle from a mesh.

Warning: Since this routine changes the order of triangles in the mesh list,
care must be taken when looping through the list of triangles to delete
some of them.  The loop should proceed from high to low index.

Entry:
  tri    - triangle to remove
  mesh   - mesh triangle belongs to
  dverts - flag saying whether to delete un-used vertices (1 = delete them,
       0 = leave them alone)
******************************************************************************/

delete_triangle(tri, mesh, dverts)
Triangle* tri;
Mesh* mesh;
int dverts;
{
    int i, j;
    int index;
    Vertex* v[3];

    /* remove mention of this triangle from its vertices */
    /* (deleting the vertices if they belong to no other triangles) */

    for (i = 0; i < 3; i++)
        remove_tri_from_vert(tri->verts[i], tri, i, mesh, dverts);

    /* check the index of the triangle */
    index = tri->index;
    if (mesh->tris[index] != tri) {
        fprintf(stderr, "delete_triangle: triangle index wrong\n");
        exit(-1);
    }

    /* overwrite tri with the last triangle in the list */
    mesh->tris[index] = mesh->tris[--mesh->ntris];
    mesh->tris[index]->index = index;

    /* free up triangle's memory */
    if (tri->more)
        free(tri->more);
    free(tri);
}


/******************************************************************************
Remove a triangle from the list of triangles of a given vertex.

Entry:
  vert   - the vertex whose list is to be shrunk
  tri    - the triangle to remove to the list
  num    - index of the vertex in the triangle
  mesh   - mesh the triangle is from
  dverts - flag saying whether to delete un-used vertices (1 = delete them,
       0 = leave them alone)
******************************************************************************/

remove_tri_from_vert(vert, tri, num, mesh, dverts)
Vertex* vert;
Triangle* tri;
int num;
Mesh* mesh;
int dverts;
{
    int i;
    int index;
    int found = 0;
    Vertex* vert2;

    /* find reference to this triangle */
    for (i = 0; i < vert->ntris; i++) {
        if (vert->tris[i] == tri) {
            found = 1;
            index = i;
            break;
        }
    }

    /* check for error */
    if (!found) {
        fprintf(stderr, "remove_tri_from_vert: can't find reference to tri\n");
        exit(-1);
    }

    /* remove the triangle from the list by overwriting it with the last */
    /* triangle in the list */
    vert->tris[index] = vert->tris[--vert->ntris];

    /* see if the vertex's vertex list needs changing */
    /* (beware, this is icky) */

    /* check one of the other vertices */
    vert2 = tri->verts[(num + 1) % 3];
    found = 0;
    /* look for vert2 in the triangles of vert */
    for (i = 0; i < vert->ntris; i++)
        if (vert->tris[i]->verts[0] == vert2 ||
            vert->tris[i]->verts[1] == vert2 ||
            vert->tris[i]->verts[2] == vert2) {
            found = 1;
            break;
        }

    /* if none of the triangles use this vertex, see where it is in list */
    /* and delete it */
    if (!found) {
        /* find vert2 in vert's list */
        found = 0;
        for (i = 0; i < vert->nverts; i++) {
            if (vert->verts[i] == vert2) {
                found = 1;
                index = i;
                break;
            }
        }
        /* error check */
        if (!found) {
            fprintf(stderr, "remove_tri_from_vert: inconsistant vertex\n");
            return;
            /*
            exit (-1);
            */
        }
        /* actually remove the reference */
        vert->verts[index] = vert->verts[--vert->nverts];
    }

    /* now check the other */
    vert2 = tri->verts[(num + 2) % 3];
    found = 0;
    /* look for vert2 in the triangles of vert */
    for (i = 0; i < vert->ntris; i++)
        if (vert->tris[i]->verts[0] == vert2 ||
            vert->tris[i]->verts[1] == vert2 ||
            vert->tris[i]->verts[2] == vert2) {
            found = 1;
            break;
        }

    /* if none of the triangles use this vertex, see where it is in list */
    /* and delete it */
    if (!found) {
        /* find vert2 in vert's list */
        found = 0;
        for (i = 0; i < vert->nverts; i++) {
            if (vert->verts[i] == vert2) {
                found = 1;
                index = i;
                break;
            }
        }
        /* error check */
        if (!found) {
            fprintf(stderr, "remove_tri_from_vert: inconsistant vertex\n");
            return;
            /*
            exit (-1);
            */
        }
        /* actually remove the reference */
        vert->verts[index] = vert->verts[--vert->nverts];
    }

    /* mark this vertex as being on the edge of the mesh */
    vert->on_edge = 1;

    /* if there are no more triangles listed for this vertex and if */
    /* the flag says it is okay to, go ahead and delete this vertex */

    if (dverts && vert->ntris == 0)
        delete_vertex(vert, mesh);
}


/******************************************************************************
Remove a vertex from a mesh.

Entry:
  vert - vertex to remove
  mesh - mesh to remove vertex from
******************************************************************************/

delete_vertex(vert, mesh)
Vertex* vert;
Mesh* mesh;
{
    int index;

    /* free up triangle and vertex list */
    free(vert->tris);
    free(vert->verts);

    /* remove this vertex from the hash table */
    remove_from_hash(vert, mesh);

    /* overwrite vertex with the last vertex in the list */
    index = vert->index;
    mesh->verts[index] = mesh->verts[--mesh->nverts];
    mesh->verts[index]->index = index;
}


/******************************************************************************
Remove all the vertices of a mesh that are used by zero triangles.

Entry:
  mesh - mesh to remove un-used triangles from
******************************************************************************/

remove_unused_verts(mesh)
Mesh* mesh;
{
    int i;
    Vertex* v;
    int count = 0;

    /* examine each vertex of the mesh */

    for (i = mesh->nverts - 1; i >= 0; i--) {

        v = mesh->verts[i];

        /* delete a vertex if it has no vertices and no triangles */
        if (v->ntris == 0 || v->nverts == 0) {
            if (v->nverts == 0 && v->nverts == 0) {
                delete_vertex(v, mesh);
                count++;
            } else { /* error check */
                fprintf(stderr, "remove_unused_verts: inconsistancy at vertex %d\n",
                        v->index);
                exit(-1);
            }
        }
    }

    /*
      printf ("deleted %d un-used vertices\n", count);
    */
}


/******************************************************************************
Check to see that a proposed triangle is okay to add to a mesh.

Entry:
  v1,v2,v3 - vertices of proposed triangle, in the order they will appear
         around the potential triangle

Exit:
  returns 1 if proposed triangle is okay, 0 if it is not
******************************************************************************/

int check_proposed_tri(v1, v2, v3)
Vertex* v1, *v2, *v3;
{
    /* check each edge of proposed triangle */

    if (check_proposed_edge(v1, v2) == 0)
        return (0);

    if (check_proposed_edge(v2, v3) == 0)
        return (0);

    if (check_proposed_edge(v3, v1) == 0)
        return (0);

    /* if we have made it here, it is okay to add the triangle */

    /* Not true!  Can have degenerate triangle (set_triangle_geometry)! */
    /* - B. Curless, 9/14/95 */

    return (1);
}


/******************************************************************************
Check to see that a proposed new edge (for a new triangle) is okay to add.
This is for avoiding edges shared by more than two triangles, and also
to make sure that new triangles attach on the opposite side of an edge from
a triangle that already uses an edge.

Entry:
  v1,v2 - vertices of proposed edge, in the order they will appear around
      the potential triangle

Exit:
  returns 1 if proposed edge is okay, 0 if it is not
******************************************************************************/

int check_proposed_edge(v1, v2)
Vertex* v1, *v2;
{
    int i, j, k;
    Triangle* tri;
    Triangle* uses_v1;
    int use_count;
    int found;
    int index, found_index;

    /* find any triangles around v1 that already use the proposed edge */
    use_count = 0;
    for (i = 0; i < v1->ntris; i++) {

        tri = v1->tris[i];

        /* find which vertex of this triangle is v1 */
        found = 0;
        for (j = 0; j < 3; j++) {
            if (tri->verts[j] == v1) {
                index = j;
                found = 1;
                break;
            }
        }

        /* error check */
        if (!found) {
            fprintf(stderr,
                    "check_proposed_edge: triangle doesn't point to correct vertex, #1\n");
            exit(-1);
        }

        /* count this triangle as using the potential edge, and remember it */
        if (tri->verts[(index + 1) % 3] == v2 || tri->verts[(index + 2) % 3] == v2) {
            uses_v1 = tri;
            found_index = index;
            use_count++;
        }
    }

    /* it's okay to add the proposed edge if no triangles use it yet */
    if (use_count == 0)
        return (1);

    /* if we've got exactly one triangle using this edge, we need to */
    /* see which orientation it is */
    if (use_count == 1) {
        if (uses_v1->verts[(found_index + 1) % 3] == v2) {
#ifdef DEBUG_CLIP
            printf("same order\n");
#endif
            return (0);  /* same order as proposed edge */
        } else
            return (1);  /* opposite order, which is okay */
    }

    /* we can't add the proposed edge if two triangles already use it */
    if (use_count == 2) {
#ifdef DEBUG_CLIP
        printf("two already\n");
#endif
        return (0);
    }

    /* if we get here, we've got more than two triangles using the */
    /* edge, and this is WRONG! */

#ifdef DEBUG_CLIP
    fprintf(stderr,
            "check_proposed_edge: too many edges already: %d\n", use_count);
#endif

    return (0);
}


/******************************************************************************
Add a triangle to the list of triangles of a given vertex.

Entry:
  vert - the vertex whose list is to be expanded
  tri  - the triangle to add to the list
******************************************************************************/

add_tri_to_vert(vert, tri)
Vertex* vert;
Triangle* tri;
{
    /* maybe allocate more room for triangle list */
    if (vert->ntris >= vert->max_tris) {
        vert->max_tris += 4;
        vert->tris =
            (Triangle**) realloc(vert->tris, sizeof(Triangle*) * vert->max_tris);
    }

    /* add the triangle to the list */
    vert->tris[vert->ntris++] = tri;
}


/******************************************************************************
Compute the geometric information of a triangle from its vertices.
******************************************************************************/

int
set_triangle_geometry(tri)
Triangle* tri;
{
    int result;

    result = plane_thru_vectors(tri->verts[0], tri->verts[1], tri->verts[2],
                                &tri->aa, &tri->bb, &tri->cc, &tri->dd);

    if (result == 1) {
        return result;
    }

    /*
      tri->normal[X] = -tri->aa;
      tri->normal[Y] = -tri->bb;
      tri->normal[Z] = -tri->cc;
    */

    result = compute_edge_planes(tri);

    return result;
}


/******************************************************************************
Determine equation of the plane containing three vectors.

Entry:
  v0,v1,v2 - vectors

Exit:
  aa,bb,cc,dd - describes plane
  returns 0 if completed normally, 1 if degenerate plane
******************************************************************************/

int plane_thru_vectors(v0, v1, v2, aa, bb, cc, dd)
Vector v0, v1, v2;
float* aa, *bb, *cc, *dd;
{
    double a, b, c, d;
    double len;
    double recip;

    a =  v0[Y] * (v1[Z] - v2[Z]) + v1[Y] * (v2[Z] - v0[Z]) +
         v2[Y] * (v0[Z] - v1[Z]);
    b = -v0[X] * (v1[Z] - v2[Z]) - v1[X] * (v2[Z] - v0[Z]) -
        v2[X] * (v0[Z] - v1[Z]);
    c =  v0[X] * (v1[Y] - v2[Y]) + v1[X] * (v2[Y] - v0[Y]) +
         v2[X] * (v0[Y] - v1[Y]);
    d =  -a * v0[X] - b * v0[Y] - c * v0[Z];

    len = sqrt(a * a + b * b + c * c);
    if (len == 0) {
        fprintf(stderr, "plane_thru_vectors: degenerate point configuration\n");
        /*
        printf ("v0: %f %f %f\n", v0[X], v0[Y], v0[Z]);
        printf ("v1: %f %f %f\n", v1[X], v1[Y], v1[Z]);
        printf ("v2: %f %f %f\n", v2[X], v2[Y], v2[Z]);
        */
        *aa = *bb = *cc = *dd = 0;
        return (1);
    }

    recip = 1.0 / len;
    *aa = a * recip;
    *bb = b * recip;
    *cc = c * recip;
    *dd = d * recip;

    return (0);
}


/******************************************************************************
Compute planes that pass through edges of a triangle and are perpendicular
to the plane containing the triangle.

Entry:
  tri - triangle to compute edge planes for
******************************************************************************/

int
compute_edge_planes(tri)
Triangle* tri;
{
    int i;
    Vector v;
    float* v0, *v1, *v2;
    float a, b, c, d;

    /* find plane through two vertices and perpendicular to plane of triangle */

    for (i = 0; i < 3; i++) {

        v0 = tri->verts[i]->coord;
        v1 = tri->verts[(i + 1) % 3]->coord;
        v2 = tri->verts[(i + 2) % 3]->coord;

        /* make a vertex out of v0 perpendicular to plane of polygon */
        v[X] = v0[X] + tri->aa;
        v[Y] = v0[Y] + tri->bb;
        v[Z] = v0[Z] + tri->cc;

        /* find plane through these three vertices */
        if (plane_thru_vectors(v, v0, v1, &a, &b, &c, &d) == 1) {
            return 1;
        }

        /* points in polygon are positive when plugged into plane equation */
        if (v2[X] * a + v2[Y] * b + v2[Z] * c + d > 0) {
            tri->a[i] = a;
            tri->b[i] = b;
            tri->c[i] = c;
            tri->d[i] = d;
        } else {
            tri->a[i] = -a;
            tri->b[i] = -b;
            tri->c[i] = -c;
            tri->d[i] = -d;
        }
    }
    return 0;
}


/******************************************************************************
Calculate vertex normals by averaging the normals of the vertice's triangles.

Entry:
  mesh - mesh to find vertex normals for
******************************************************************************/

find_vertex_normals(mesh)
Mesh* mesh;
{
    int i;

    /* go through all vertices of the mesh */

    for (i = 0; i < mesh->nverts; i++)
        find_vertex_normal(mesh->verts[i]);
}


/******************************************************************************
Find normal at a vertex.

Entry:
  vert - vertex at which to find surface normal
******************************************************************************/

find_vertex_normal(vert)
Vertex* vert;
{
    int j;
    Vector sum;

    /* special case for vertices with NO triangles */
    if (vert->ntris == 0) {
        vert->normal[X] = 0;
        vert->normal[Y] = 0;
        vert->normal[Z] = 0;
        return;
    }

    sum[X] = sum[Y] = sum[Z] = 0;

    /* add surface normal of all triangles of the vertex */
    for (j = 0; j < vert->ntris; j++) {
        /*
        vadd (vert->tris[j]->normal, sum, sum);
        */
        sum[X] -= vert->tris[j]->aa;
        sum[Y] -= vert->tris[j]->bb;
        sum[Z] -= vert->tris[j]->cc;
    }

    /* normalize this sum of normals and save it away */
    vnorm(sum);
    vcopy(sum, vert->normal);
}


/******************************************************************************
Find which vertices of a mesh are on the edge of the mesh.

Entry:
  mesh - mesh to find edges of
******************************************************************************/

find_mesh_edges(mesh)
Mesh* mesh;
{
    int i;

    /* examine each vertex of a mesh to see if it's on the edge */

    for (i = 0; i < mesh->nverts; i++)
        vertex_edge_test(mesh->verts[i]);
}


/******************************************************************************
Mark a vertex to say whether it is on the edge of a mesh or not.

Entry:
  vert: vertex to mark as on edge or not

Exit:
  returns 0 if everything was okay, 1 if the mesh is build funny
******************************************************************************/

int vertex_edge_test(vert)
Vertex* vert;
{
    int i;
    Triangle* tri;
    unsigned char on_edge;
    static int been_here = 0;
    int bad_mesh = 0;

    /* initialize list that counts how many times an edge has been marked */
    /* (each edge should be marked twice, once for each triangle, */
    /*  otherwise the vertex is on an edge) */

    for (i = 0; i < vert->nverts; i++)
        vert->verts[i]->count = 0;

    /* go through each triangle of vertex, counting how many times an */
    /* adjacent vertex has been used */

    for (i = 0; i < vert->ntris; i++) {
        tri = vert->tris[i];

        /* count the participation of each vertex of the triangle */
        tri->verts[0]->count++;
        tri->verts[1]->count++;
        tri->verts[2]->count++;
    }

    /* examine the counts of the neighboring vertices to see if they */
    /* prove the vertex to be an edge.  also check for mesh consistancy */

    on_edge = 0;  /* assume we're not on an edge */

    for (i = 0; i < vert->nverts; i++) {
        if (vert->verts[i]->count == 1) {
            on_edge = 1;
        } else if (vert->verts[i]->count == 2) {
            /* this is okay */
        } else {
            /* if we're here, our mesh is built wrong */
            bad_mesh = 1;
            if (!been_here) {
                been_here = 1;
                fprintf(stderr, "vertex_edge_test: %d count on an edge\n",
                        vert->verts[i]->count);
                /*
                fprintf (stderr, "(You'll only get this message once, so beware!)\n");
                       */
            }
        }
    }

    vert->on_edge = on_edge;

    /* say if the mesh is improperly built */
    return (bad_mesh);
}


/******************************************************************************
Do all sorts of checks to make sure mesh is consistant.  Remove vertices and
triangles that cause problems.

Entry:
  scan - scan containing the mesh to clean up
******************************************************************************/

clean_up_mesh(scan)
Scan* scan;
{
    int i, j, k;
    Mesh* mesh;
    int count1, count2;
    Vertex* vert;
    Vertex* v;
    Vertex* v1, *v2;
    Triangle* tri, *tri2;
    Triangle* shared_triangle();
    int abnormal, num_adj;
    int nverts;
    int num;
    static Vertex* near_verts[100];
    static Triangle* near_tris[3];

    mesh = scan->meshes[mesh_level];

    /* compute which vertices are on the boundary */
    find_mesh_edges(mesh);




    /* delete vertices where 4 edges meet */
    /* and eliminate triple edges (or greater) */

    count1 = 0;
    count2 = 0;
    for (i = 0; i < mesh->nverts; i++) {

        vert = mesh->verts[i];

        /* ignore vertices not on the mesh edge */
        if (!vert->on_edge)
            continue;

        /* initialize list that counts how many times an edge has been marked */
        /* (each edge should be marked twice, once for each triangle, */
        /*  otherwise the vertex is on an edge) */

        for (j = 0; j < vert->nverts; j++)
            vert->verts[j]->count = 0;

        /* go through each triangle of vertex, counting how many times an */
        /* adjacent vertex has been used */

        for (j = 0; j < vert->ntris; j++) {
            tri = vert->tris[j];
            tri->verts[0]->count++;
            tri->verts[1]->count++;
            tri->verts[2]->count++;
        }

        /* examine the counts of the neighboring vertices to see if they */
        /* tell if the vertex has more than 2 edges meeting or if some edges */
        /* are triple edges or worse */

        abnormal = 0;
        num_adj = 0;
        for (j = 0; j < vert->nverts; j++) {
            v = vert->verts[j];
            /* see if this vertex has and edge shared by too many triangles */
            if (v->count != 1 && v->count != 2) {
                abnormal = 1;
                break;
            }
            /* count the number of edges meeting at the vertex */
            if (v->count == 1) {
                num_adj++;
            }
        }

        /* delete the triangles from any abnormal vertex */
        if (abnormal || num_adj > 2) {

            if (abnormal)
                count1++;
            else
                count2++;

            /* save a list of adjacent vertices */
            nverts = vert->nverts;
            for (j = 0; j < nverts; j++)
                near_verts[j] = vert->verts[j];

            /* delete the vertice's triangles */
            for (j = vert->ntris - 1; j >= 0; j--)
                delete_triangle(vert->tris[j], mesh, 0);

            /* make sure there are no adjacent vertices left */
            if (vert->nverts != 0) {
                fprintf(stderr, "clean_up_mesh: still vertices left in vertex!\n");
            }

            /* re-check to see if the other vertices are on an edge */
            for (j = 0; j < nverts; j++)
                vertex_edge_test(near_verts[j]);
        }
    }

#ifdef VERBOSE

    if (count1)
        printf("Got rid of %d vertices that have more than 2 edges meeting\n",
               count1);

    if (count2)
        printf("Deleted %d vertices where there were triple edges.\n", count2);

#endif


    /* get rid of pairs of triangles that share more than one edge */
    count1 = 0;
    count2 = 0;
    for (i = mesh->ntris - 1; i >= 0; i--) {

        tri = mesh->tris[i];
        abnormal = 0;

        /* collect adjacent triangles in near_tris[] */
        for (j = 0; j < 3; j++) {

            num = edges_shared_count(tri->verts[j], tri->verts[(j + 1) % 3]);

            if (num > 2) {
                abnormal = 1;
                break;
            }

            if (num == 2) {
                if (shared_triangle(0) == tri)
                    near_tris[j] = shared_triangle(1);
                else if (shared_triangle(1) == tri)
                    near_tris[j] = shared_triangle(0);
                else {
                    abnormal = 1;
                }
            } else {
                near_tris[j] = NULL;
            }
        }

        /* check to see if any of the adjacent triangles are duplicated */
        if (!abnormal) {
            if (near_tris[0] && near_tris[1] && near_tris[0] == near_tris[1])
                abnormal = 1;
            if (near_tris[0] && near_tris[2] && near_tris[0] == near_tris[2])
                abnormal = 1;
            if (near_tris[1] && near_tris[2] && near_tris[1] == near_tris[2])
                abnormal = 1;
        }

        /* delete the triangle if it is abnormal */
        if (abnormal) {
            count1++;
            delete_triangle(tri, mesh, 0);
            continue;
        }

        /* check if any adjacent triangles have their edges going the other way */
        for (j = 0; j < 3; j++) {
            if (near_tris[j] == NULL)
                continue;
            v1 = tri->verts[j];
            v2 = tri->verts[(j + 1) % 3];
            tri2 = near_tris[j];
            if (tri2->verts[0] == v1 && tri2->verts[1] == v2)
                abnormal = 1;
            if (tri2->verts[1] == v1 && tri2->verts[2] == v2)
                abnormal = 1;
            if (tri2->verts[2] == v1 && tri2->verts[0] == v2)
                abnormal = 1;
        }

        /* delete the triangle if it is abnormal */
        if (abnormal) {
            count2++;
            delete_triangle(tri, mesh, 0);
        }

    }

#ifdef VERBOSE

    if (count1) {
        printf("Deleted %d tris that shared more than one edge with a neighbor\n",
               count1);
    }

    if (count2) {
        printf("Deleted %d tris with vertex order opposite that of a neighbor\n",
               count2);
    }

#endif


    /* check for degenerate triangles (those where two or more */
    /* vertices are the same */

    count1 = 0;
    for (i = mesh->ntris - 1; i >= 0; i--) {
        tri = mesh->tris[i];
        if (tri->verts[0] == tri->verts[1] ||
            tri->verts[0] == tri->verts[2] ||
            tri->verts[1] == tri->verts[2]) {
            count1++;
            delete_triangle(tri, mesh, 0);
        }
    }

#ifdef VERBOSE

    if (count1)
        printf("Found and deleted %d degenerate triangles.\n", count1);

#endif

    /* get rid of un-used vertices */
    remove_unused_verts(mesh);

    /* re-compute which vertices are on the boundary */
    find_mesh_edges(mesh);

    /* re-compute all vertex normals */
    find_vertex_normals(mesh);
}

