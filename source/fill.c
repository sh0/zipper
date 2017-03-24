/*

Nice boundary loop filling.

Greg Turk, December 1993

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

#include "zipper.h"
#include "matrix.h"
#include "meshops.h"

/* vertex that is part of a hole that is being filled */
typedef struct FillVertex {
    Vertex* vert;         /* vertex helping cover the hole */
    unsigned char on_edge;    /* on the edge of the hole? */
    Vertex* e1, *e2;      /* adjacent vertices on edge, if any */
} FillVertex;

/* triangle used to fill a hole */
typedef struct FillTri {
    Triangle* tri;        /* triangle used to fill hole */
    int index;            /* index into list containing triangle */
} FillTri;

static FillVertex** fverts = NULL;
static int num_fverts = 0;
static int max_fverts = 50;

static FillTri** ftris = NULL;
static int num_ftris = 0;
static int max_ftris = 50;

// Parameters
static float FILL_EDGE_LENGTH_FACTOR;
static float FILL_EDGE_LENGTH;

void update_fill_resolution()
{
    FILL_EDGE_LENGTH = ZIPPER_RESOLUTION * FILL_EDGE_LENGTH_FACTOR;
}

void set_fill_edge_length_factor(float factor)
{
    FILL_EDGE_LENGTH_FACTOR = factor;
    FILL_EDGE_LENGTH = ZIPPER_RESOLUTION * FILL_EDGE_LENGTH_FACTOR;
}

float get_fill_edge_length_factor()
{
    return FILL_EDGE_LENGTH_FACTOR;
}

void better_fill_loop(int loop, Scan* scan);
int fix_fill_size(int loop, Scan* scan, float max_len, float* size);
int maybe_split_edge(Triangle* tri, Scan* scan, int idx1, int idx2);
void init_fill_lists();
void new_ftri(Triangle* tri);
void delete_ftri(FillTri* ftri);
void new_fvert(Vertex* vert, int on_edge, Vertex* v1, Vertex* v2);
void fill_tri_split2(FillTri* ftri, Mesh* mesh, int index);
void fill_tri_split3(FillTri* ftri, Mesh* mesh, int index);
void fill_tri_split4(FillTri* ftri, Mesh* mesh, int index);
void smooth_hole_vertices(Scan* scan);
void swap_hole_edges(Scan* sc);

/******************************************************************************
Fill in a loop of a mesh.

Entry:
  loop - index of loop to fill
  scan - scan containing loop
******************************************************************************/
void better_fill_loop(int loop, Scan* scan)
{
    int i;
    Mesh* mesh;
    Edge* e, *fedge;
    int result;
    int been_around;
    Vertex** vlist;
    int vcount = 0;
    int index;
    int self_intersect;
    int p1, p2, p3;
    Vector norm;
    Vector vec;
    Triangle* new_tri;
    Vertex* vert;
    int num;
    float init_size, size;
    int inc;
    float max_len;
    int iter_count;

    mesh = scan->meshes[mesh_level];

    /* go around the entire loop to see if all edges are oriented properly */
    /* and also come up with rough "normal" */

    vset(norm, 0.0, 0.0, 0.0);
    fedge = mesh->looplist.loops[loop];
    been_around = 0;
    for (e = fedge; e != fedge || !been_around; e = e->next) {
        been_around = 1;
        result = check_proposed_edge(e->v2, e->v1);
        if (!result) {
            fprintf(stderr, "can't fill this loop\n");
            return;
        }
        vadd(norm, e->v1->normal, norm);
        vcount++;
    }

    /* find coefficients for plane equation */
    vnorm(norm);

    /* initialize the polygon splitter */
    result = init_splitter(norm[X], norm[Y], norm[Z], 0.0);
    if (result) {
        fprintf(stderr, "can't fill this loop\n");
        return;
    }

    /* make a list of vertices around the loop */
    vlist = (Vertex**) malloc(sizeof(Vertex*) * vcount);
    index = 0;
    been_around = 0;
    for (e = fedge; e != fedge || !been_around; e = e->prev) {
        been_around = 1;
        vlist[index++] = e->v1;
    }

    /* transform the loop's vertices to the xy-plane */
    /* and send them to the splitter */

    for (i = 0; i < vcount; i++) {
        vcopy(vlist[i]->coord, vec);
        add_boundary_point(vec[X], vec[Y], vec[Z], i);
    }

    /* call the splitter */
    self_intersect = greedy_connect();

    /* check to see if we successfully filled the loop */
    if (self_intersect) {
        fprintf(stderr, "can't fill this loop\n");
        return;
    }

    /* mark all vertices in mesh as not part of the hole */
    for (i = 0; i < mesh->nverts; i++)
        mesh->verts[i]->move_to = NULL;

    /* initialize list of fill triangles and vertices */
    init_fill_lists();

    /* create the list of vertices around the hole */
    for (i = 0; i < vcount; i++) {
        new_fvert(vlist[i], 1, vlist[(i + vcount - 1) % vcount], vlist[(i + 1) % vcount]);
    }

    /* create the new triangles */
    for (i = 0; i < get_ntris(); i++) {
        get_triangle(i, &p1, &p2, &p3);
        if (check_proposed_tri(vlist[p1], vlist[p2], vlist[p3])) {
            new_tri = make_triangle(mesh, vlist[p1], vlist[p2], vlist[p3], 1e20);

            /* check_proposed_tri() is inadequate!  - B. Curless  9/14/95 */
            if (new_tri != NULL)
                new_ftri(new_tri);
        }
    }

    /* re-calculate vertex info */
    for (i = 0; i < vcount; i++) {
        vertex_edge_test(vlist[i]);
        find_vertex_normal(vlist[i]);
    }

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;

    free(vlist);

    /* determine maximum allowable edge length */
    inc = level_to_inc(mesh_level);
    max_len = inc * FILL_EDGE_LENGTH;

    /* make the triangles fit a maximum size requirement */

#if 0
    for (i = 0; i < 3; i++) {
        printf("iteration %d\n", i);
        fix_fill_size(loop, scan, max_len, &size);
        swap_hole_edges(scan);
        smooth_hole_vertices(scan);
    }
#endif

    swap_hole_edges(scan);
    fix_fill_size(loop, scan, max_len, &init_size);

    /*
      printf ("initial size = %f\n", init_size);
      printf ("ratio = %f\n", init_size / max_len);
    */

    iter_count = 1;
    while ((num = fix_fill_size(loop, scan, max_len, &size)) > 0 &&
           iter_count < init_size / max_len) {
        /*
        printf ("size = %f\n", size);
        printf ("%d new vertices\n", num);
        */
        swap_hole_edges(scan);
        smooth_hole_vertices(scan);
        smooth_hole_vertices(scan);
        swap_hole_edges(scan);
        smooth_hole_vertices(scan);
        smooth_hole_vertices(scan);
        swap_hole_edges(scan);
        iter_count++;
    }

    /* fix the geometry of all the triangles that fill the hole */
    for (i = 0; i < mesh->ntris; i++)
        set_triangle_geometry(mesh->tris[i]);

    /* compute normals and edge conditions for all new vertices */

    for (i = 0; i < num_fverts; i++) {
        if (fverts[i]->on_edge)
            continue;
        vert = fverts[i]->vert;

        /* normals and edge conditions */
        find_vertex_normal(vert);
        vertex_edge_test(vert);

        /* add new vertices to the hash table */
        add_to_hash(vert, mesh);
    }

    /* remove the "more" triangle info */

    for (i = 0; i < num_ftris; i++)
        if (ftris[i]->tri->more) {
            free(ftris[i]->tri->more);
            ftris[i]->tri->more = NULL;
        }

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;

    printf("done filling hole\n");
}

/******************************************************************************
Make sure that the triangles that fill a hole are no larger than a particular
maximum size.

Entry:
  loop    - index of loop to fill
  scan    - scan containing loop
  max_len - maximum allowable edge length (split any that are larger)

Exit:
  size - size of largest edge
  returns number of new vertices created
******************************************************************************/
int fix_fill_size(int loop, Scan* scan, float max_len, float* size)
{
    int i, j;
    Triangle* tri;
    FillVertex* verts[3];
    Vertex* v1, *v2;
    Vector diff;
    float len;
    int n1, n2, n3;
    int code;
    Mesh* mesh = scan->meshes[mesh_level];
    FillTri* ftri;
    int num_created = 0;
    float max_size = -1e20;

    /* check the size of all the triangles filling the hole */
    for (i = 0; i < num_ftris; i++) {
        ftri = ftris[i];
        tri = ftris[i]->tri;

        /* get pointers to the fill vertices of this triangle */
        verts[0] = (FillVertex*) tri->verts[0]->move_to;
        verts[1] = (FillVertex*) tri->verts[1]->move_to;
        verts[2] = (FillVertex*) tri->verts[2]->move_to;

        /* examine each edge of the triangle to see if it should be split */
        for (j = 0; j < 3; j++) {

            /* find length of edge */
            v1 = verts[j]->vert;
            v2 = verts[(j + 1) % 3]->vert;
            vsub(v1->coord, v2->coord, diff);
            len = vlen(diff);

            /* see if this is new maximum length edge */
            if (len > max_size)
                max_size = len;

            /* if edge is too long, mark it for splitting */
#if 0
            printf("max len: %f %f\n", max_len, len);
#endif
            if (len > max_len) {
                num_created += maybe_split_edge(tri, scan, j, (j + 1) % 3);
            }
        }

    }

    /* check each triangle and see if it needs to be split */
    for (i = num_ftris - 1; i >= 0; i--) {

        ftri = ftris[i];
        tri = ftris[i]->tri;

        n1 = (tri->more->mids[0] != NULL);
        n2 = (tri->more->mids[1] != NULL);
        n3 = (tri->more->mids[2] != NULL);
        code = (n3 << 2) | (n2 << 1) | n1;

        switch (code) {
            case 0:
                break;
            case 1:
                fill_tri_split2(ftri, mesh, 0);
                break;
            case 2:
                fill_tri_split2(ftri, mesh, 1);
                break;
            case 3:
                fill_tri_split3(ftri, mesh, 0);
                break;
            case 4:
                fill_tri_split2(ftri, mesh, 2);
                break;
            case 5:
                fill_tri_split3(ftri, mesh, 2);
                break;
            case 6:
                fill_tri_split3(ftri, mesh, 1);
                break;
            case 7:
                fill_tri_split4(ftri, mesh, 0);
                break;
        }
    }

    /* return the number of new vertices */
    *size = max_size;
    return (num_created);
}


/******************************************************************************
Create a new vertex at the midpoint of an edge, if this has not already
been done.

Entry:
  tri       - triangle whose edge will be split
  scan      - mesh that tri belongs to
  idx1,idx2 - indices of vertices that define the edge

Exit:
  returns 1 if it split the edge, 0 if not
******************************************************************************/
int maybe_split_edge(Triangle* tri, Scan* scan, int idx1, int idx2)
{
    int shared;
    Vertex* v1, *v2;
    Triangle* shared_triangle();
    Triangle* t1, *t2;
    FillVertex* fv1, *fv2;
    Triangle* other_tri;
    int other_index;
    int index;
    Vertex* new_vert;
    Vector pos;
    Mesh* mesh;

    mesh = scan->meshes[mesh_level];

    v1 = tri->verts[idx1];
    v2 = tri->verts[idx2];

    /* see if this edge is on the hole's boundary, and if so */
    /* then don't split the edge */

    fv1 = (FillVertex*) v1->move_to;
    fv2 = (FillVertex*) v2->move_to;
    if (fv1->on_edge && fv2->on_edge) {
        if (fv1->e1 == v2 || fv1->e2 == v2 ||
            fv2->e1 == v1 || fv2->e2 == v1)
            return (0);
    }

    /* see if we've already split this edge */

    if (tri->more->mids[idx1] != NULL) {
        return (0);
    }

    /* look for which triangles share this edge */

    shared = edges_shared_count(v1, v2);

    /* consistancy check */
    if (shared != 1 && shared != 2) {
        fprintf(stderr, "maybe_split_tri: wrong number of shared tris = %d\n",
                shared);
        exit(-1);
    }

    /* if the edge is shared, find which triangle shares it */

    other_tri = NULL;

    if (shared == 2) {

        t1 = shared_triangle(0);
        t2 = shared_triangle(1);

        if (t1 == tri)
            other_tri = t2;
        else if (t2 == tri)
            other_tri = t1;
        else {
            fprintf(stderr, "maybe_split_tri: can't find triangle\n");
            exit(-1);
        }

        /* find out indices for the edge vertices */
        if (other_tri->verts[0] == v2 && other_tri->verts[1] == v1)
            other_index = 0;
        else if (other_tri->verts[1] == v2 && other_tri->verts[2] == v1)
            other_index = 1;
        else if (other_tri->verts[2] == v2 && other_tri->verts[0] == v1)
            other_index = 2;
        else {
            fprintf(stderr, "maybe_split_tri: can't find other edge\n");
            exit(-1);
        }
    }

    /* compute the edge's midpoint */
    pos[X] = 0.5 * (v1->coord[X] + v2->coord[X]);
    pos[Y] = 0.5 * (v1->coord[Y] + v2->coord[Y]);
    pos[Z] = 0.5 * (v1->coord[Z] + v2->coord[Z]);

    /* make the vertex */
    index = make_vertex(mesh, pos);
    new_vert = mesh->verts[index];

    /* add the new vertex to the appropriate triangle(s) */
    tri->more->mids[idx1] = new_vert;
    if (shared == 2) {
        if (other_tri->more->mids[other_index] != NULL)
            fprintf(stderr, "maybe_split_tri: already a vertex here\n");
        other_tri->more->mids[other_index] = new_vert;
    }

    /* add the vertex to the list of hole-filling vertices */
    new_fvert(new_vert, 0, NULL, NULL);

    /* signal that we created a new vertex */
    return (1);
}

/******************************************************************************
Initialize the lists of fill vertices and triangles.
******************************************************************************/
void init_fill_lists()
{
    int i;

    /* either create or clear out the fill vertices */

    if (fverts == NULL) {
        fverts = (FillVertex**) malloc(sizeof(FillVertex*) * max_fverts);
        num_fverts = 0;
    } else {
        for (i = 0; i < num_fverts; i++)
            free(fverts[i]);
        num_fverts = 0;
    }

    /* create or clear out the fill triangles */

    if (ftris == NULL) {
        ftris = (FillTri**) malloc(sizeof(FillTri*) * max_ftris);
        num_ftris = 0;
    } else {
        for (i = 0; i < num_ftris; i++)
            free(ftris[i]);
        num_ftris = 0;
    }
}

/******************************************************************************
Add a new triangle to the list of fill triangles.

Entry:
  tri      - pointer to new triangle to add to list
******************************************************************************/
void new_ftri(Triangle* tri)
{
    FillTri* ftri;

    /* make sure there is room for this triangle */
    if (num_ftris >= max_ftris) {
        max_ftris *= 2;
        ftris = (FillTri**) realloc(ftris, sizeof(FillTri*) * max_ftris);
    }

    /* add the new triangle to the list */

    ftri = (FillTri*) malloc(sizeof(FillTri));
    ftri->tri = tri;
    ftri->index = num_ftris;
    if (tri->more == NULL)
        tri->more = (More_Tri_Stuff*) malloc(sizeof(More_Tri_Stuff));
    tri->more->mids[0] = NULL;
    tri->more->mids[1] = NULL;
    tri->more->mids[2] = NULL;
    tri->more->fill_tri = ftri;

    ftris[num_ftris] = ftri;
    num_ftris++;
}

/******************************************************************************
Delete a triangle being used to fill a hole.
******************************************************************************/
void delete_ftri(FillTri* ftri)
{
    if (ftri->tri->more)
        free(ftri->tri->more);

    ftris[ftri->index] = ftris[--num_ftris];
    ftris[ftri->index]->index = ftri->index;
    free(ftri);
}

/******************************************************************************
Add a vertex to the list of vertices that help fill a hole.

Entry:
  vert    - vertex to put on list
  on_edge - is this vertex on the hole's boundary?
  v1,v2   - adjacent vertices on hole's boundary, if on_edge is 1
******************************************************************************/
void new_fvert(Vertex* vert, int on_edge, Vertex* v1, Vertex* v2)
{
    FillVertex* fvert;

    /* make sure there is room for this vertex */
    if (num_fverts >= max_fverts) {
        max_fverts *= 2;
        fverts = (FillVertex**)
                 realloc(fverts, sizeof(FillVertex*) * max_fverts);
    }

    /* add the vertex to the list */

    fvert = (FillVertex*) malloc(sizeof(FillVertex));
    fvert->vert = vert;
    fvert->on_edge = on_edge;
    if (on_edge) {
        fvert->e1 = v1;
        fvert->e2 = v2;
    } else {
        fvert->e1 = NULL;
        fvert->e2 = NULL;
    }

    fverts[num_fverts] = fvert;
    num_fverts++;

    /* have vertex point back to position in this list */
    vert->move_to = (Vertex*) fvert;
}

/******************************************************************************
Split a hole-filling triangle into two triangles.

Entry:
  ftri  - triangle to split
  mesh  - mesh tri comes from
  index - index of midpoint to use to split with
******************************************************************************/
void fill_tri_split2(FillTri* ftri, Mesh* mesh, int index)
{
    Vertex* v1, *v2, *v3;
    Vertex* m1;
    Triangle* tri = ftri->tri;
    Triangle* new_tri;

    v1 = tri->verts[index];
    v2 = tri->verts[(index + 1) % 3];
    v3 = tri->verts[(index + 2) % 3];

    m1 = tri->more->mids[index];

    /* delete the original triangle from the hole triangle list */
    delete_ftri(ftri);

    /* delete triangle from the mesh */
    delete_triangle(tri, mesh, 0);

    /* create the new triangles */

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, v1, m1, v3, 1e20);
    new_ftri(new_tri);

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, m1, v2, v3, 1e20);
    new_ftri(new_tri);
}

/******************************************************************************
Split a hole-filling triangle into three triangles.

Entry:
  ftri  - triangle to split
  mesh  - mesh tri comes from
  index - index of midpoint to use to split with
******************************************************************************/
void fill_tri_split3(FillTri* ftri, Mesh* mesh, int index)
{
    Vertex* v1, *v2, *v3;
    Vertex* m1, *m2;
    Triangle* tri = ftri->tri;
    Triangle* new_tri;

    v1 = tri->verts[index];
    v2 = tri->verts[(index + 1) % 3];
    v3 = tri->verts[(index + 2) % 3];

    m1 = tri->more->mids[index];
    m2 = tri->more->mids[(index + 1) % 3];

    /* delete the original triangle from the hole triangle list */
    delete_ftri(ftri);

    /* delete triangle from the mesh */
    delete_triangle(tri, mesh, 0);

    /* create the new triangles */

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, v1, m1, m2, 1e20);
    new_ftri(new_tri);

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, v1, m2, v3, 1e20);
    new_ftri(new_tri);

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, m1, v2, m2, 1e20);
    new_ftri(new_tri);
}

/******************************************************************************
Split a hole-filling triangle into four triangles.

Entry:
  ftri  - triangle to split
  mesh  - mesh tri comes from
  index - index of midpoint to use to split with
******************************************************************************/
void fill_tri_split4(FillTri* ftri, Mesh* mesh, int index)
{
    Vertex* v1, *v2, *v3;
    Vertex* m1, *m2, *m3;
    Triangle* tri = ftri->tri;
    Triangle* new_tri;

    v1 = tri->verts[index];
    v2 = tri->verts[(index + 1) % 3];
    v3 = tri->verts[(index + 2) % 3];

    m1 = tri->more->mids[index];
    m2 = tri->more->mids[(index + 1) % 3];
    m3 = tri->more->mids[(index + 2) % 3];

    /* delete the original triangle from the hole triangle list */
    delete_ftri(ftri);

    /* delete triangle from the mesh */
    delete_triangle(tri, mesh, 0);

    /* create the new triangles */
    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, v1, m1, m3, 1e20);
    new_ftri(new_tri);

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, v2, m2, m1, 1e20);
    new_ftri(new_tri);

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, v3, m3, m2, 1e20);
    new_ftri(new_tri);

    /* Beware NULL triangles! - B. Curless 9/14/95 */
    new_tri = make_triangle(mesh, m1, m2, m3, 1e20);
    new_ftri(new_tri);
}

/******************************************************************************
Use Laplacian smoothing to put the vertices of a hole in better positions.

Entry:
  scan - scan that has a hole that is being filled
******************************************************************************/
void smooth_hole_vertices(Scan* scan)
{
    int i;
    FillVertex* fvert;
    Vertex* vert;
    Vector new_pos;

    /* find a smoothed position for all hole-filling vertices that */
    /* are not on the boundary of the hole */

    for (i = 0; i < num_fverts; i++) {
        fvert = fverts[i];
        if (fvert->on_edge)
            continue;
        vert = fvert->vert;
        compute_smoothing(vert, new_pos);
        vcopy(new_pos, vert->normal);  /* oh no! bad place to save it!*/
    }

    /* actually move the vertices */

    for (i = 0; i < num_fverts; i++) {
        fvert = fverts[i];
        if (fvert->on_edge)
            continue;
        vert = fvert->vert;
        vcopy(vert->normal, vert->coord);
    }
}

/******************************************************************************
Look for edges that should be "swapped" by deleting their common triangles and
creating two triangles with a shared edge that goes in the other direction.

Entry:
  sc - scan containing mesh
******************************************************************************/
void swap_hole_edges(Scan* sc)
{
    int i, j, k;
    Mesh* mesh;
    Vertex* vert1, *vert2;
    Vertex* v1, *v2, *v3, *v4;
    Vector dir1, dir2, dir3, dir4;
    float dot1, dot2, dot3, dot4;
    Triangle* t1, *t2;
    int hit1, hit3;
    int count;
    Triangle* new_tri;
    FillVertex* fv1, *fv2;

    mesh = sc->meshes[mesh_level];

    /* look at all edges in mesh */
    for (i = 0; i < num_fverts; i++) {
        vert1 = fverts[i]->vert;

        for (j = vert1->nverts - 1; j >= 0; j--) {
            vert2 = vert1->verts[j];

            /* make sure the edge between vertices isn't on the hole's boundary */
            fv1 = (FillVertex*) vert1->move_to;
            fv2 = (FillVertex*) vert2->move_to;
            if (fv1 == NULL || fv2 == NULL)
                continue;
            if (fv1->on_edge && fv2->on_edge) {
                if (fv1->e1 == v2 || fv1->e2 == v2 ||
                    fv2->e1 == v1 || fv2->e2 == v1)
                    continue;
            }

            /* see whether this edge is shared by exactly two triangles */
            count = edges_shared_count(vert1, vert2);
            if (count != 2)
                continue;
            t1 = shared_triangle(0);
            t2 = shared_triangle(1);

            /* make sure these triangles are both hole-filling tris */
            if (t1->more == NULL || t1->more->fill_tri == NULL ||
                t2->more == NULL || t2->more->fill_tri == NULL)
                continue;

            /* figure out how the four vertices are arranged around the triangles */
            hit1 = 0;
            hit3 = 0;
            v2 = vert1;
            v4 = vert2;
            for (k = 0; k < 3; k++) {
                if (t1->verts[k] == vert1 && t1->verts[(k + 1) % 3] == vert2) {
                    v1 = t1->verts[(k + 2) % 3];
                    hit1++;
                }
                if (t1->verts[k] == vert2 && t1->verts[(k + 1) % 3] == vert1) {
                    v3 = t1->verts[(k + 2) % 3];
                    hit3++;
                }
                if (t2->verts[k] == vert1 && t2->verts[(k + 1) % 3] == vert2) {
                    v1 = t2->verts[(k + 2) % 3];
                    hit1++;
                }
                if (t2->verts[k] == vert2 && t2->verts[(k + 1) % 3] == vert1) {
                    v3 = t2->verts[(k + 2) % 3];
                    hit3++;
                }
            }

            /* sanity check */
            if (hit1 != 1 || hit3 != 1) {
                fprintf(stderr, "swap_edges: incorrect values: %d %d\n", hit1, hit3);
                continue;
            }

            /* make sure these vertices are hole-filling vertices */
            if (v1->move_to == NULL || v3->move_to == NULL)
                continue;

            /* create normalized vectors around the boundary */
            vsub(v2->coord, v1->coord, dir1);
            vnorm(dir1);
            vsub(v3->coord, v2->coord, dir2);
            vnorm(dir2);
            vsub(v4->coord, v3->coord, dir3);
            vnorm(dir3);
            vsub(v1->coord, v4->coord, dir4);
            vnorm(dir4);

            /* compute dot products of edges leaving each vertex */
            dot1 = -vdot(dir1, dir4);
            dot2 = -vdot(dir2, dir1);
            dot3 = -vdot(dir3, dir2);
            dot4 = -vdot(dir4, dir3);

            /* create the two triangles that are most nearly equiangular */
            if (dot1 + dot3 < dot2 + dot4) {
                delete_ftri(t1->more->fill_tri);
                delete_triangle(t1, mesh, 0);
                delete_ftri(t2->more->fill_tri);
                delete_triangle(t2, mesh, 0);

                /* Beware NULL triangles! - B. Curless 9/14/95 */
                new_tri = make_triangle(mesh, v1, v2, v3, 1e20);
                new_ftri(new_tri);

                /* Beware NULL triangles! - B. Curless 9/14/95 */
                new_tri = make_triangle(mesh, v1, v3, v4, 1e20);
                new_ftri(new_tri);
            }
        }
    }
}
