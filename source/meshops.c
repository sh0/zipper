/*

Perform "fix-up" mesh operations.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "zipper.h"
#include "matrix.h"
#include "meshops.h"

typedef enum {FALSE = 0, TRUE} boolean;

/* really bad way to get triangles out of edges_shared_count() */
static Triangle* tris_shared[10];


absorb_transform(sc)
Scan* sc;
{
    int i;
    Mesh* mesh;
    Vertex* vert;

    mesh = sc->meshes[mesh_level];

    for (i = 0; i < mesh->nverts; i++) {
        vert = mesh->verts[i];
        mat_apply(sc->rotmat, vert->coord);
        vert->coord[0] += sc->xtrans;
        vert->coord[1] += sc->ytrans;
        vert->coord[2] += sc->ztrans;
        mat_apply(sc->rotmat, vert->normal);
    }

    mat_ident(sc->rotmat);
    sc->xtrans = 0;
    sc->ytrans = 0;
    sc->ztrans = 0;
}


/******************************************************************************
Fix up the bows at vertices (where a vertex has more than two unshared edges).

Entry:
  sc - scan containing mesh to fix
******************************************************************************/

fix_bows(sc)
Scan* sc;
{
    int i, j, k;
    Mesh* mesh;
    Vertex* v1;
    Triangle* tri;
#define VVMAX 20
    static Vertex* vin[VVMAX], *vout[VVMAX];
    int in_out_count;
    int found;
    int unshared;
    int result;

    mesh = sc->meshes[mesh_level];

    /* examine all vertices of this mesh */

    for (k = 0; k < mesh->nverts; k++) {
        v1 = mesh->verts[k];

        /* ignore vertices not on the mesh edge */
        if (!v1->on_edge)
            continue;

        /* see that we have enough room to handle the vertex */
        if (v1->ntris > VVMAX) {
            fprintf(stderr, "fix_bows: not enough room to handle vertex\n");
            continue;
        }

        /* look at all adjacent vertices */

        /* initialize list that counts how many times an edge has been marked */
        /* (each edge should be marked twice, once for each triangle, */
        /*  otherwise the vertex is on an edge) */

        for (i = 0; i < v1->nverts; i++)
            v1->verts[i]->count = 0;

        /* go through each triangle of vertex, counting how many times an */
        /* adjacent vertex has been used */

        for (i = 0; i < v1->ntris; i++) {
            tri = v1->tris[i];

            /* count the participation of each vertex of the triangle */
            tri->verts[0]->count++;
            tri->verts[1]->count++;
            tri->verts[2]->count++;
        }

        /* see how many edges are unshared.  if less than or equal to two are, */
        /* we don't have a bow and we can go to next vertex */
        unshared = 0;
        for (i = 0; i < v1->nverts; i++) {
            if (v1->verts[i]->count == 1) {
                unshared++;
            } else if (v1->verts[i]->count == 2) {
                /* this is okay */
            } else {
                /* this is not okay */
#ifdef DEBUG_CLIP
                fprintf(stderr, "Error! fix_bows: %d count on an edge\n",
                        v1->verts[i]->count);
#endif
            }
        }
        if (unshared <= 2)
            continue;

        /* collect together the vertices adjacent to v1, according to triangles */
        for (i = 0; i < v1->ntris; i++) {
            tri = v1->tris[i];
            if (tri->verts[0] == v1) {
                vout[i] = tri->verts[1];
                vin[i] = tri->verts[2];
            } else if (tri->verts[1] == v1) {
                vout[i] = tri->verts[2];
                vin[i] = tri->verts[0];
            } else if (tri->verts[2] == v1) {
                vout[i] = tri->verts[0];
                vin[i] = tri->verts[1];
            } else {
                fprintf(stderr, "fix_bows: triangle doesn't point to correct vertex\n");
                exit(-1);
            }
        }

        in_out_count = v1->ntris;

        do {

            found = 0;

            for (i = 0; i < in_out_count; i++) {
                for (j = 0; j < in_out_count; j++) {
                    if (i == j)
                        continue;
                    /* consolidate pairs if we get a match */
                    if (vout[j] == vin[i]) {
                        /* first update one of the pairs */
                        vin[i] = vin[j];
                        /* then delete the other pair */
                        in_out_count--;
                        vin[j] = vin[in_out_count];
                        vout[j] = vout[in_out_count];
                        /* make note that we consolidated */
                        found = 1;
                        goto here;
                    }
                }
            }

        here: ; /* easy way out of double loop */

        } while (found);

        /*
        printf ("in_out_count: %d, ntris: %d\n", in_out_count, v1->ntris);
        */

        /* fill in bow with triangles until there is only one gap left */

        while (in_out_count > 1) {
            result = fill_bow(mesh, v1, vin, vout, in_out_count);
            /* bail if there was some problem in fill_bow() */
            if (result) {
                for (i = 0; i < v1->ntris; i++)
                    v1->tris[i]->mark = 2;
                break;
            }
            in_out_count--;
        }
    }

    remove_unused_verts(mesh);

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Fill in the smallest gap at bow with a triangle.

Entry:
  mesh         - mesh in which to fill
  vert         - vertex to fill at
  vin          - incoming edges to vertex
  vout         - outgoing edges from vertex
  in_out_count - count of vin (and vout)

Exit:
  returns 1 if there was some problem, 0 if not
******************************************************************************/

int fill_bow(mesh, vert, vin, vout, in_out_count)
Mesh* mesh;
Vertex* vert;
Vertex* vin[], *vout[];
int in_out_count;
{
    int i, j;
    int ii, jj;
    static Vector din[20], dout[20];
    extern float edge_length_max();
    float max, dot;

    /* make vectors for each vin and vout */
    for (i = 0; i < in_out_count; i++) {
        vsub(vin[i]->coord, vert->coord, din[i]);
        vnorm(din[i]);
        vsub(vout[i]->coord, vert->coord, dout[i]);
        vnorm(dout[i]);
    }

    /* find minimum dot product between din and dout */
    max = -1e20;
    for (i = 0; i < in_out_count; i++)
        for (j = 0; j < in_out_count; j++) {
            if (i == j)
                continue;
            dot = vdot(din[i], dout[j]);
            if (dot > max) {
                max = dot;
                ii = i;
                jj = j;
            }
        }

    /* create triangle between vin[ii] and vout[jj] */

    /* Not used?? */
    edge_length_max(mesh_level);

    /* check unlikely case that we're not allowed to add such a new triangle */
    if (!check_proposed_tri(vert, vin[ii], vout[jj])) {
        int shared_count;
        Triangle* tri;
        int found;

        shared_count = edges_shared_count(vin[ii], vout[jj]);

        /* If we would be creating a triple-edge, see if one or the other */
        /* of the two existing polygons at that edge is ONLY connected at */
        /* that edge.  If so, we can delete it. */

        if (shared_count == 2) {

            tri = tris_shared[0];
            found = 0;
            for (i = 0; i < 3; i++)
                if (tri->verts[i]->ntris == 1) {
                    found = 1;
                    break;
                }
            if (found)
                delete_triangle(tri, mesh, 0);

            tri = tris_shared[1];
            found = 0;
            for (i = 0; i < 3; i++)
                if (tri->verts[i]->ntris == 1) {
                    found = 1;
                    break;
                }
            if (found)
                delete_triangle(tri, mesh, 0);

            if (found) {
                find_vertex_normal(vin[ii]);
                find_vertex_normal(vout[jj]);
            }
        }

        /* see if we've fixed anything */
        if (!check_proposed_tri(vert, vin[ii], vout[jj])) {
#ifdef DEBUG_CLIP
            printf("fill_bow: YUP, WE'VE GOT A PROBLEM HERE\n");
#endif
            return (1);
        }
    }

    make_triangle(mesh, vert, vin[ii], vout[jj], 1e20);

    /* re-assess whether these vertices are on an edge of the mesh */

    vertex_edge_test(vert);
    vertex_edge_test(vin[ii]);
    vertex_edge_test(vout[jj]);

    /* re-compute the normal direction at these vertices */

    find_vertex_normal(vert);
    find_vertex_normal(vin[ii]);
    find_vertex_normal(vout[jj]);

    /* combine vout[ii] and vin[jj] */

    vin[ii] = vin[jj];
    in_out_count--;
    vin[jj] = vin[in_out_count];
    vout[jj] = vout[in_out_count];

    /* say everything was okay */
    return (0);
}


/******************************************************************************
Split a triangle along an edge to form two new triangles.

Entry:
  sc     - scan containing triangle to split
  tri    - triangle to split
  index1 - index to first vertex on edge to split (0, 1 or 2)
       (second vertex is one greater, mod 3)
  t      - where on edge to split (0 = at v1, 1 = at v2)

Exit:
  tri1,tri2 - pointers to the two newly created triangles
******************************************************************************/

split_triangle(sc, tri, index1, t, tri1, tri2)
Scan* sc;
Triangle* tri;
int index1;
float t;
Triangle** tri1, ** tri2;
{
    Mesh* mesh;
    int index2, index3;
    Vertex* v1, *v2, *v3;
    Vector coord;
    int vert_index;
    Vertex* new_vert;
    Triangle* new1, *new2;

    mesh = sc->meshes[mesh_level];
    index2 = (index1 + 1) % 3;
    index3 = (index1 + 2) % 3;
    v1 = tri->verts[index1];
    v2 = tri->verts[index2];
    v3 = tri->verts[index3];

    /* create vertex on edge along which to split */
    coord[X] = v1->coord[X] + t * (v2->coord[X] - v1->coord[X]);
    coord[Y] = v1->coord[Y] + t * (v2->coord[Y] - v1->coord[Y]);
    coord[Z] = v1->coord[Z] + t * (v2->coord[Z] - v1->coord[Z]);
    vert_index = make_vertex(mesh, coord);
    new_vert = mesh->verts[vert_index];

    /* create new triangles */
    new1 = make_triangle(mesh, v1, new_vert, v3, 1e20);
    new2 = make_triangle(mesh, new_vert, v2, v3, 1e20);

    /* delete the original triangle */
    delete_triangle(tri, mesh, 1);

    /* re-assess whether these vertices are on an edge of the mesh */
    vertex_edge_test(v1);
    vertex_edge_test(v2);
    vertex_edge_test(v3);
    vertex_edge_test(new_vert);

    /* compute the normal at the new vertex */
    find_vertex_normal(new_vert);

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;

    /* return pointers to the two new triangles */
    *tri1 = new1;
    *tri2 = new2;
}


/******************************************************************************
Test the triangle splitting routine.

Entry:
  sc - scan containing triangles to split
******************************************************************************/

split_test(sc)
Scan* sc;
{
    int i, j;
    Mesh* mesh;
    Vertex* v1, *v2, *v3;
    int found;
    Triangle* tri;
    Triangle* new1, *new2;

    mesh = sc->meshes[mesh_level];

    /* split triangles at random */

    for (i = 0; i < 50; i++) {
        found = 0;
        do {
            j = floor(drand48() * mesh->ntris);
            tri = mesh->tris[j];
            v1 = tri->verts[0];
            v2 = tri->verts[1];
            v3 = tri->verts[2];
            if (v1->on_edge && v2->on_edge && edges_shared_count(v1, v2) == 1) {
                split_triangle(sc, tri, 0, 0.5, &new1, &new2);
                found = 1;
            } else if (v2->on_edge && v3->on_edge && edges_shared_count(v2, v3) == 1) {
                split_triangle(sc, tri, 1, 0.5, &new1, &new2);
                found = 1;
            } else if (v3->on_edge && v1->on_edge && edges_shared_count(v3, v1) == 1) {
                split_triangle(sc, tri, 2, 0.5, &new1, &new2);
                found = 1;
            }
        } while (found == 0);
    }
}


/******************************************************************************
Determine how many triangles share a pair of vertices.

Entry:
  v1,v2 - two vertices to check

Exit:
  returns number of edges shared
  side effect is to place shared triangles in the array tris_shared[]
******************************************************************************/

int edges_shared_count(v1, v2)
Vertex* v1, *v2;
{
    int i, j;
    int count = 0;
    Triangle* tri;
    int found1, found2;

    /* examine all triangles of the first vertex */

    for (i = 0; i < v1->ntris; i++) {

        tri = v1->tris[i];

        /* does the first vertex appear in this triangle? */
        found1 = 0;
        for (j = 0; j < 3; j++)
            if (tri->verts[j] == v1) {
                found1 = 1;
                break;
            }

        /* does the second vertex appear in this triangle? */
        found2 = 0;
        for (j = 0; j < 3; j++)
            if (tri->verts[j] == v2) {
                found2 = 1;
                break;
            }

        /* if they both appear in the triangle, add one to edge count */
        if (found1 && found2) {
            if (count < 10)
                tris_shared[count] = tri;
            count++;
        }
    }

    /* return the number of edges found that share both vertices */
    return (count);
}


/******************************************************************************
Return the nth shared triangle.
******************************************************************************/

Triangle* shared_triangle(index)
int index;
{
    return (tris_shared[index]);
}


/******************************************************************************
Fill in the small holes in a mesh with triangles.

Entry:
  sc - scan containing mesh to fix
******************************************************************************/

fill_small_holes(sc)
Scan* sc;
{
    int i;
    Mesh* mesh;
    Edge** loops;
    int nloops;
    int been_around;
    Edge* edge;
    int edge_count;
    Triangle* last_tri;
    int tri_count;

    /* fill in any "bows" in the mesh */
    fix_bows(sc);

    mesh = sc->meshes[mesh_level];

    /* get the list of edge loops in the mesh */
    if (!mesh->edges_valid)
        create_edge_list(mesh);
    nloops = mesh->looplist.nloops;
    loops = mesh->looplist.loops;

    /*
    printf ("nloops = %d\n", nloops);
    */

    /* examine each loop to see if it is small enough to fill with triangles */

    for (i = 0; i < nloops; i++) {

        /* Count the number of edges in a loop.  Also count the number */
        /* of distinct triangles around the loop. */

        been_around = 0;
        edge_count = 0;
        last_tri = loops[i]->prev->tri;
        tri_count = 0;

        for (edge = loops[i]; edge != loops[i] || !been_around; edge = edge->next) {
            if (edge->tri != last_tri) {
                last_tri = edge->tri;
                tri_count++;
            }
            edge_count++;
            been_around = 1;
        }

        if (tri_count == 0)
            tri_count = 1;

        /*
        printf ("edges = %d, tris = %d\n", edge_count, tri_count);
        */

        /* fill in the hole if it has the same number of triangles surrounding */
        /* it as the number of edges AND if it has just 3 or 4 edges */

        if (tri_count == edge_count && (edge_count == 3 || edge_count == 4)) {
            /*
                  printf ("filling %d hole\n", edge_count);
            */
            fill_hole(mesh, loops[i], edge_count);
        }
#if 1
        else if (tri_count == edge_count && (edge_count <= 10)) {
            /*
                  printf ("filling %d loop\n", edge_count);
            */
            fill_loop(i, sc);
        }
#endif
    }
}


/******************************************************************************
Fill in a small hole in a mesh.

Entry:
  mesh       - mesh in which there is a hole to fill
  edge       - one of the edges along the hole
  edge_count - number of edges around the hole
******************************************************************************/

fill_hole(mesh, edge, edge_count)
Mesh* mesh;
Edge* edge;
int edge_count;
{
    Vertex* v1, *v2, *v3, *v4;

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;

    /* handle the 3-edge case */
    if (edge_count == 3) {
        v3 = edge->v1;
        v2 = edge->v2;
        v1 = edge->next->v2;
        if (!check_proposed_tri(v1, v2, v3))
            return;
        make_triangle(mesh, v1, v2, v3, 1e20);
        vertex_edge_test(v1);
        vertex_edge_test(v2);
        vertex_edge_test(v3);
        find_vertex_normal(v1);
        find_vertex_normal(v2);
        find_vertex_normal(v3);
        return;
    }

    /* the 4-edge case */
    if (edge_count == 4) {
        v4 = edge->prev->v1;
        v3 = edge->v1;
        v2 = edge->v2;
        v1 = edge->next->v2;
        fill_four_hole(mesh, v1, v2, v3, v4);
        return;
    }
}


/******************************************************************************
Fill in a four-edged hole in a mesh.

Entry:
  mesh        - mesh in which there is a hole to fill
  v1,v2,v3,v4 - vertices around the hole
******************************************************************************/

fill_four_hole(mesh, v1, v2, v3, v4)
Mesh* mesh;
Vertex* v1, *v2, *v3, *v4;
{
    Vector dir1, dir2, dir3, dir4;
    float dot1, dot2, dot3, dot4;

    /* create normalized vectors around the edge of the hole */
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
        if (!check_proposed_tri(v1, v2, v3) || !check_proposed_tri(v1, v3, v4))
            return;
        make_triangle(mesh, v1, v2, v3, 1e20);
        make_triangle(mesh, v1, v3, v4, 1e20);
    } else {
        if (!check_proposed_tri(v1, v2, v4) || !check_proposed_tri(v2, v3, v4))
            return;
        make_triangle(mesh, v1, v2, v4, 1e20);
        make_triangle(mesh, v2, v3, v4, 1e20);
    }

    /* re-compute whether each vertex is on an edge */
    vertex_edge_test(v1);
    vertex_edge_test(v2);
    vertex_edge_test(v3);
    vertex_edge_test(v4);

    /* re-compute vertex normals */
    find_vertex_normal(v1);
    find_vertex_normal(v2);
    find_vertex_normal(v3);
    find_vertex_normal(v4);
}


/******************************************************************************
Remove the "cut" vertices from a mesh that were introduced during zippering.

Entry:
  scan - scan containing the zippered mesh
******************************************************************************/

remove_cut_vertices(scan)
Scan* scan;
{
    int i;
    Mesh* mesh;
    int count;
    Vertex* v;
    int removed;
    int on_edge;
    int result;

    mesh = scan->meshes[mesh_level];

    /* examine each vertex to see if it should be removed */

    count = 0;
    removed = 0;
    on_edge = 0;
    for (i = mesh->nverts - 1; i >= 0; i--) {

        v = mesh->verts[i];

        /* only look at vertices that have been marked as "cut" vertices */
        /* from zippering */
        if (v->moving)
            count++;
        else
            continue;

        /* un-mark the vertex */
        v->moving = 0;

        /* see if it is on an edge */
        vertex_edge_test(v);

        /* don't attempt to move vertices on an edge */
        if (v->on_edge) {
            on_edge++;
            continue;
        }

        /* remove the vertex */
        result = remove_a_vertex(scan, v);

        if (result)
            removed++;
    }

    /* remove any vertices that are used by no triangles */
    remove_unused_verts(mesh);

#ifdef DEBUG_CLIP
    printf("%d cut vertices were examined\n", count);
    printf("%d were removed\n", removed);
    printf("%d were on edge\n", on_edge);
#endif
}


/******************************************************************************
Remove a vertex and re-tile the surrounding region.

Entry:
  scan - scan containing the vertex to remove
  v    - the vertex to remove

Exit:
  returns 1 if successful at removing vertex, 0 if not
******************************************************************************/

int remove_a_vertex(scan, v)
Scan* scan;
Vertex* v;
{
    int j, k;
    Mesh* mesh;
    Triangle* tri;
    Vertex** vin, **vout;
    Vertex* tvert;
    int found;
    int self_intersect;
    Vector vec;
    int result;
    float w;
    int p1, p2, p3;
    int get_ntris();
    int ntris;

    mesh = scan->meshes[mesh_level];

    /* allocate room for vertices surrounding "v" */
    vin  = (Vertex**) malloc(sizeof(Vertex*) * v->ntris);
    vout = (Vertex**) malloc(sizeof(Vertex*) * v->ntris);

    /* collect together the vertices adjacent to v, according to triangles */
    for (j = 0; j < v->ntris; j++) {
        tri = v->tris[j];
        if (tri->verts[0] == v) {
            vout[j] = tri->verts[1];
            vin[j] = tri->verts[2];
        } else if (tri->verts[1] == v) {
            vout[j] = tri->verts[2];
            vin[j] = tri->verts[0];
        } else if (tri->verts[2] == v) {
            vout[j] = tri->verts[0];
            vin[j] = tri->verts[1];
        } else {
            fprintf(stderr,
                    "remove_a_vertex: triangle doesn't point to correct vertex\n");
            return (0);
        }
    }

    /* order the vertices around v */
    for (j = 0; j < v->ntris - 2; j++) {
        found = 0;
        for (k = 0; k < v->ntris; k++)
            if (vout[j] == vin[k]) {
                tvert = vin[k];
                vin[k] = vin[j + 1];
                vin[j + 1] = tvert;
                tvert = vout[k];
                vout[k] = vout[j + 1];
                vout[j + 1] = tvert;
                found = 1;
                break;
            }
        if (!found) {
#ifdef DEBUG_CLIP
            fprintf(stderr,
                    "remove_a_vertex: can't order triangles around vertex\n");
#endif
            return (0);
        }
    }

    /* find transformation from tangent plane at vertex to xy-plane */
    w = -vdot(v->normal, v->coord);

    /* re-tile the surrounding vertices WITHOUT vertex "v" */
    result = init_splitter(v->normal[X], v->normal[Y], v->normal[Z], w);
    if (result) {
#ifdef DEBUG_CLIP
        fprintf(stderr, "remove_a_vertex: bad return from init_splitter\n");
#endif
        return (0);
    }

    /* transform the face's vertices to the xy-plane */
    /* and send them to the splitter */

    ntris = v->ntris;
    for (j = v->ntris - 1; j >= 0; j--) {
        vcopy(vin[j]->coord, vec);
        add_boundary_point(vec[X], vec[Y], vec[Z], j);
    }

    /* call the splitter */
    self_intersect = greedy_connect();

    if (self_intersect) {
#ifdef DEBUG_CLIP
        fprintf(stderr, "self-intersection while removing cut point\n");
#endif
        return (0);
    }

    /* remove the original triangles */
    for (j = v->ntris - 1; j >= 0; j--)
        delete_triangle(v->tris[j], mesh, 0);

    /* get the newly-formed triangles from the splitting routine */
    /* and create these triangles */
    for (j = 0; j < get_ntris(); j++) {
        get_triangle(j, &p1, &p2, &p3);
        if (check_proposed_tri(vin[p1], vin[p2], vin[p3]))
            make_triangle(mesh, vin[p1], vin[p2], vin[p3], 1e20);
#ifdef DEBUG_CLIP
        else
            printf("remove_a_vertex: proposed triangle bad\n");
#endif
    }

    /* re-compute normals and whether vertices are on an edge */
    for (j = 0; j < ntris; j++) {
        find_vertex_normal(vin[j]);
        vertex_edge_test(vin[j]);
    }

    /* free space */
    free(vin);
    free(vout);

    /* say we removed vertex successfully */
    return (1);
}


/******************************************************************************
Remove long, thin polygons from a mesh.

Entry:
  scan  - mesh to remove slivers from
  fract - length of shortest allowed altitude, as a fraction of the spacing
      between range points
******************************************************************************/

remove_sliver_tris(scan, fract)
Scan* scan;
float fract;
{
    int i;
    Mesh* mesh;
    float* v1, *v2, *v3;
    float d1, d2, d3;
    float edge_length_max();
    float min;
    Triangle* tri;
    int removed;

    mesh = scan->meshes[mesh_level];

    /* minimum allowable length is a fraction of the maximum allowed */
    min = fract * level_to_inc(mesh_level) * ZIPPER_RESOLUTION;

    /* examine pairs of vertices connected by edges */

    for (i = 0; i < mesh->ntris; i++) {

        tri = mesh->tris[i];
        v1 = tri->verts[0]->coord;
        v2 = tri->verts[1]->coord;
        v3 = tri->verts[2]->coord;

        /* find altitudes of the triangle */
        d1 = v1[X] * tri->a[1] + v1[Y] * tri->b[1] + v1[Z] * tri->c[1] + tri->d[1];
        d2 = v2[X] * tri->a[2] + v2[Y] * tri->b[2] + v2[Z] * tri->c[2] + tri->d[2];
        d3 = v3[X] * tri->a[0] + v3[Y] * tri->b[0] + v3[Z] * tri->c[0] + tri->d[0];

        /* see if any altitude is too small */

        if (d1 < min || d2 < min || d3 < min) {

            removed = FALSE;

            /* remove the vertex with the shortest altitude */
            if (d1 < d2 && d1 < d3) {
                if (tri->verts[0]->on_edge == 0)
                    removed = removed || remove_a_vertex(scan, tri->verts[0]);
            } else if (d2 < d1 && d2 < d3) {
                if (tri->verts[1]->on_edge == 0)
                    removed = removed || remove_a_vertex(scan, tri->verts[1]);
            } else if (d3 < d1 && d3 < d2) {
                if (tri->verts[2]->on_edge == 0)
                    removed = removed || remove_a_vertex(scan, tri->verts[2]);
            }

            /* If one or more of the triangles vertices are successfully
               removed, then this triangle is gone and has been replaced
               by the last triangle in the list.  Ergo, we decrement "i"
               to visit the next triangle (formerly the last triangle) */
            if (removed)
                i--;
        }
    }

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Remove vertices that have very similar normals to their neighbors

Doesn't work well.  Generates large flat polygons; need to account for
the size of the contributing polygons, since the normal measure is not
meaningful when the polygon sizes differ significantly.

Entry:
  scan  - mesh to remove vertices from
  cos_max - maximum dot product allowed, else removal
******************************************************************************/

remove_flat_verts(scan, cos_max)
Scan* scan;
float cos_max;
{
    int i, j;
    Mesh* mesh;
    float* n1, *n2;
    float min_dot, dot;
    Vertex* vert;

    /*
      printf("%f\n", cos_max);
    */

    mesh = scan->meshes[mesh_level];

    /* Check dot products */
    for (i = 0; i < mesh->nverts; i++) {
        vert = mesh->verts[i];
        vert->moving = 0;

        if (vert->on_edge) {
            continue;
        }

        n1 = vert->normal;
        min_dot = 1;
        for (j = 0; j < vert->nverts; j++) {
            n2 = vert->verts[j]->normal;
            dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
            min_dot = MIN(min_dot, dot);
        }
        if (min_dot > cos_max) {
            vert->moving = 1;
        }
    }

    /* Remove the marked vertices */
    remove_cut_vertices(scan);

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Remove polygons with bad aspect ratios.  Currently no working very well.

Entry:
  scan   - mesh to remove slivers from
  max_aspect - aspect ratio above which we should consider vertices for removal
  min_cos - minimum cosine of angles between normals; intended to preserve
            regions of high curvature
  diff - the maximum difference between the total number of triangles
         that vertex belongs to, and the number of bad triangles
         a vertex belongs to
******************************************************************************/

remove_bad_aspect_tris(scan, max_aspect, min_cos, diff)
Scan* scan;
float max_aspect;
float min_cos;
int diff;
{
    int i;
    Mesh* mesh;
    float* v1, *v2, *v3, vec[3];
    float* n1, *n2, *n3;
    float len1, len2, len3;
    float aspect;
    float dot1, dot2, dot3;
    float edge_length_max();
    float min_len, max_len;
    Triangle* tri;

    printf("Max aspect = %f.\nMin cos = %f\n", max_aspect, min_cos);

    mesh = scan->meshes[mesh_level];

    /* We'll use the "moving" field to store counts of how many
       bad triangles a vertex belongs to */


    /* Clear the "moving" field */
    for (i = 0; i < mesh->nverts; i++) {
        mesh->verts[i]->moving = 0;
    }


    /* Count number of bad triangles per vertex */
    for (i = 0; i < mesh->ntris; i++) {
        tri = mesh->tris[i];
        v1 = tri->verts[0]->coord;
        v2 = tri->verts[1]->coord;
        v3 = tri->verts[2]->coord;
        n1 = tri->verts[0]->normal;
        n2 = tri->verts[1]->normal;
        n3 = tri->verts[2]->normal;


        /* Don't remove vertices in areas of high curvature */
        dot1 = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
        dot2 = n1[0] * n3[0] + n1[1] * n3[1] + n1[2] * n3[2];
        dot3 = n3[0] * n2[0] + n3[1] * n2[1] + n3[2] * n2[2];

        if (dot1 < min_cos || dot2 < min_cos || dot3 < min_cos) {
            continue;
        }

        /* */
        vsub(v1, v2, vec);
        len1 = vlength(vec);

        vsub(v1, v3, vec);
        len2 = vlength(vec);

        vsub(v2, v3, vec);
        len3 = vlength(vec);

        max_len = MAX(len1, MAX(len2, len3));
        min_len = MIN(len1, MIN(len2, len3));

        aspect = max_len / min_len;
        if (aspect > max_aspect) {
            tri->verts[0]->moving += 1;
            tri->verts[1]->moving += 1;
            tri->verts[2]->moving += 1;
        }
    }

    /* Update the moving field to indicate whether or not
       a vertex should be removed within remove_cut_vertices() */

    for (i = 0; i < mesh->nverts; i++) {
        if ((mesh->verts[i]->ntris - mesh->verts[i]->moving) <= diff) {
            mesh->verts[i]->moving = 1;
        } else {
            mesh->verts[i]->moving = 0;
        }
    }

    /* Remove the marked vetices */
    remove_cut_vertices(scan);

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Move a given vertex to a new position.

Entry:
  v    - vertex to move
  pos  - it's new position
  mesh - mesh containing vertex
******************************************************************************/

move_vertex(v, pos, mesh)
Vertex* v;
Vector pos;
Mesh* mesh;
{
    int i;

    /* move the vertex */
    remove_from_hash(v, mesh);
    vcopy(pos, v->coord);
    add_to_hash(v, mesh);

    /* adjust the triangles that use the vertex */
    for (i = 0; i < v->ntris; i++)
        set_triangle_geometry(v->tris[i]);

    /* re-calculate the affected vertice's normals */
    find_vertex_normal(v);
    for (i = 0; i < v->nverts; i++)
        find_vertex_normal(v->verts[i]);
}


/******************************************************************************
Remove any short edges between vertices of two meshes that were introduced
during zippering.

Entry:
  scan  - mesh that was created during zippering
  fract - length of shortest allowed edge, as a fraction of the spacing
      between range points
******************************************************************************/

remove_short_edges(scan, fract)
Scan* scan;
float fract;
{
    int i, j;
    Mesh* mesh;
    Vertex* v1, *v2;
    Vector diff;
    float len;
    float min;
    float edge_length_max();

    mesh = scan->meshes[mesh_level];

    /* minimum allowable length is a fraction of the maximum allowed */
    min = fract * level_to_inc(mesh_level) * ZIPPER_RESOLUTION;

    /* examine pairs of vertices connected by edges */

    for (i = 0; i < mesh->nverts; i++) {

    here_again:
        ;

        if (i < mesh->nverts)
            v1 = mesh->verts[i];
        else
            break;

        for (j = 0; j < v1->nverts; j++) {

            v2 = v1->verts[j];

            /* see if edge is short, and if it is then collapse the edge */
            vsub(v1->coord, v2->coord, diff);
            len = vlen(diff);

            /* don't collapse edge if the vertices are not on a common boundary */
            if (v1->on_edge && v2->on_edge) {
                if (edges_shared_count(v1, v2) != 1)
                    continue;
            }

            if (len < min) {
                collapse_edge(mesh, v1, v2);
                goto here_again;  /* oh wow, i'm so bad! */
            }
        }
    }

    /* remove all vertices that are not used by any triangles */
    remove_unused_verts(mesh);

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Collapse an edge between two vertices and combine the vertices in the process.

Entry:
  mesh  - mesh in which this is taking place
  v1,v2 - the two vertices to combine
******************************************************************************/

collapse_edge(mesh, v1, v2)
Mesh* mesh;
Vertex* v1, *v2;
{
    int i, j;
    int v1ntris, v2ntris;
    Triangle* tri;
    int del_count = 0;
    Vertex* w1, *w2, *w3;

    /* average the two vertex positions, and set this to be their position */
    remove_from_hash(v1, mesh);
    remove_from_hash(v2, mesh);
    v1->coord[X] = 0.5 * (v1->coord[X] + v2->coord[X]);
    v1->coord[Y] = 0.5 * (v1->coord[Y] + v2->coord[Y]);
    v1->coord[Z] = 0.5 * (v1->coord[Z] + v2->coord[Z]);
    vcopy(v1->coord, v2->coord);
    add_to_hash(v1, mesh);
    add_to_hash(v2, mesh);

    /* delete all triangles common to both vertices */
    v1ntris = v1->ntris;
    v2ntris = v2->ntris;
    for (i = v1->ntris - 1; i >= 0; i--) {
        tri = v1->tris[i];
        for (j = 0; j < v2->ntris; j++)
            if (v2->tris[j] == tri) {
                delete_triangle(tri, mesh, 1);
                v1ntris--;
                v2ntris--;
                del_count++;
                break;
            }
        /* bail if we've deleted a vertex (by deleting all its triangles) */
        if (v1ntris == 0 || v2ntris == 0) {
            /*
            printf ("early bail\n");
            */
            return;
        }
    }

    /*
    printf ("deleted %d triangles\n", del_count);
    */

    /* move all of v2's triangles into v1 */
    for (i = v2->ntris - 1; i >= 0; i--) {

        tri = v2->tris[i];

        if (tri->verts[0] == v2) {
            w1 = v1;
            w2 = tri->verts[1];
            w3 = tri->verts[2];
        } else if (tri->verts[1] == v2) {
            w1 = tri->verts[0];
            w2 = v1;
            w3 = tri->verts[2];
        } else if (tri->verts[2] == v2) {
            w1 = tri->verts[0];
            w2 = tri->verts[1];
            w3 = v1;
        } else {
            fprintf(stderr, "collapse_edge: can't find vertex reference\n");
            exit(-1);
        }

        delete_triangle(tri, mesh, 0);

        if (check_proposed_tri(w1, w2, w3))
            make_triangle(mesh, w1, w2, w3, 1e20);
#ifdef DEBUG_CLIP
        else
            printf("collape_edge: tried to create invalid triangle\n");
#endif
    }

    /* re-compute normals and whether vertices are on an edge */
    for (i = 0; i < v1->nverts; i++) {
        find_vertex_normal(v1->verts[i]);
        vertex_edge_test(v1->verts[i]);
    }
    find_vertex_normal(v1);
    vertex_edge_test(v1);
}


/******************************************************************************
Chop each triangle of a mesh into four triangles.
******************************************************************************/

quarter_mesh(scan)
Scan* scan;
{
    int i, j, k;
    Mesh* mesh;
    Vertex* v1, *v2, *v3;
    Vertex* m1, *m2, *m3;
    int count;
    Triangle* tri;
    Vertex* vert;
    Vector mid;
    int index;
    float big = 1e20;

    mesh = scan->meshes[mesh_level];

    /* allocate extra information at each triangle in the mesh */

    for (i = 0; i < mesh->ntris; i++) {
        tri = mesh->tris[i];
        tri->more = (More_Tri_Stuff*) malloc(sizeof(More_Tri_Stuff));
        tri->more->mids[0] = NULL;
        tri->more->mids[1] = NULL;
        tri->more->mids[2] = NULL;
    }

    /* examine each edge of a mesh to see if we need to add a new */
    /* vertex at its midpoint */

    for (i = mesh->nverts - 1; i >= 0; i--) {

        v1 = mesh->verts[i];

        /* look at all vertices radiating from v1 and check these edges */

        for (j = 0; j < v1->nverts; j++) {

            v2 = v1->verts[j];

            /* see how many triangles share this edge */
            count = edges_shared_count(v1, v2);

            /* error check */
            if (count < 1 || count > 2) {
                printf("quarter_mesh: bad count = %d\n", count);
                continue;
            }

            /* only examine each edge once */
            if (count == 2 && v1->index < v2->index)
                continue;

            /* create a vertex at the midpoint of this edge */
            mid[X] = 0.5 * (v1->coord[X] + v2->coord[X]);
            mid[Y] = 0.5 * (v1->coord[Y] + v2->coord[Y]);
            mid[Z] = 0.5 * (v1->coord[Z] + v2->coord[Z]);
            index = make_vertex(mesh, mid);
            vert = mesh->verts[index];
            add_to_hash(vert, mesh);

            /* store pointer to this midpoint at each triangle that uses this edge */
            for (k = 0; k < count; k++) {
                tri = tris_shared[k];
                if ((tri->verts[0] == v1 && tri->verts[1] == v2) ||
                    (tri->verts[0] == v2 && tri->verts[1] == v1)) {
                    tri->more->mids[0] = vert;
                } else if ((tri->verts[1] == v1 && tri->verts[2] == v2) ||
                           (tri->verts[1] == v2 && tri->verts[2] == v1)) {
                    tri->more->mids[1] = vert;
                } else if ((tri->verts[2] == v1 && tri->verts[0] == v2) ||
                           (tri->verts[2] == v2 && tri->verts[0] == v1)) {
                    tri->more->mids[2] = vert;
                } else {
                    fprintf(stderr, "quarter_mesh: can't find edge in triangle\n");
                    exit(-1);
                }
            }

        }
    }

    /* replace each triangle in mesh with four smaller ones */

    for (i = mesh->ntris - 1; i >= 0; i--) {

        tri = mesh->tris[i];

        m1 = tri->more->mids[0];
        m2 = tri->more->mids[1];
        m3 = tri->more->mids[2];

        /* error check */
        if (m1 == NULL || m2 == NULL || m3 == NULL) {
            fprintf(stderr, "quarter_mesh: didn't create all midpoints\n");
            exit(-1);
        }

        /* delete the old triangle */
        v1 = tri->verts[0];
        v2 = tri->verts[1];
        v3 = tri->verts[2];
        delete_triangle(tri, mesh, 0);

        /* make the new triangles */
        make_triangle(mesh, v1, m1, m3, big);
        make_triangle(mesh, v2, m2, m1, big);
        make_triangle(mesh, v3, m3, m2, big);
        make_triangle(mesh, m1, m2, m3, big);
    }

    /* recompute info about all vertices in the mesh */
    for (i = 0; i < mesh->nverts; i++) {
        find_vertex_normal(mesh->verts[i]);
        vertex_edge_test(mesh->verts[i]);
    }

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Fill in a loop of a mesh.

Entry:
  loop - index of loop to fill
  scan - scan containing loop

Exit:
  returns 0 if successful, 1 if it couldn't fill the loop
******************************************************************************/

int fill_loop(loop, scan)
int loop;
Scan* scan;
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
#if 0
            fprintf(stderr, "can't fill this loop\n");
#endif
            return (1);
        }
        vadd(norm, e->v1->normal, norm);
        vcount++;
    }

    /* find coefficients for plane equation */
    vnorm(norm);

    /* initialize the polygon splitter */
    result = init_splitter(norm[X], norm[Y], norm[Z], 0.0);
    if (result) {
#if 0
        fprintf(stderr, "can't fill this loop\n");
#endif
        return (1);
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
#if 0
        fprintf(stderr, "can't fill this loop\n");
#endif
        return (1);
    }

    /* create the new triangles */
    for (i = 0; i < get_ntris(); i++) {
        get_triangle(i, &p1, &p2, &p3);
        if (check_proposed_tri(vlist[p1], vlist[p2], vlist[p3])) {
            make_triangle(mesh, vlist[p1], vlist[p2], vlist[p3], 1e20);
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

    /* signal normal termination */
    return (0);
}


/******************************************************************************
Look for edges that should be "swapped" by deleting their common triangles and
creating two triangles with a shared edge that goes in the other direction.

Entry:
  sc - scan containing mesh
******************************************************************************/

swap_edges(sc)
Scan* sc;
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

    mesh = sc->meshes[mesh_level];

    /* look at all edges in mesh */

    for (i = 0; i < mesh->nverts; i++) {
        vert1 = mesh->verts[i];

        for (j = 0; j < vert1->nverts; j++) {
            vert2 = vert1->verts[j];

            /* see whether this edge is shared by exactly two triangles */
            count = edges_shared_count(vert1, vert2);
            if (count != 2)
                continue;
            t1 = tris_shared[0];
            t2 = tris_shared[1];

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
                delete_triangle(t1, mesh, 0);
                delete_triangle(t2, mesh, 0);
                if (!check_proposed_tri(v1, v2, v3) ||
                    !check_proposed_tri(v1, v3, v4)) {
                    fprintf(stderr, "swap_edges: error checking new tris\n");
                    continue;
                }
                make_triangle(mesh, v1, v2, v3, 1e20);
                make_triangle(mesh, v1, v3, v4, 1e20);

                /* re-compute whether each vertex is on an edge */
                vertex_edge_test(v1);
                vertex_edge_test(v2);
                vertex_edge_test(v3);
                vertex_edge_test(v4);

                /* re-compute vertex normals */
                find_vertex_normal(v1);
                find_vertex_normal(v2);
                find_vertex_normal(v3);
                find_vertex_normal(v4);
            }
        }
    }
}


/******************************************************************************
Return the average of the surrounding vertex positions adjacent to a given
vertex.

Entry:
  v - vertex to "smooth"

Exit:
  new_pos - average position of neighboring vertices
******************************************************************************/

compute_smoothing(v, new_pos)
Vertex* v;
Vector new_pos;
{
    int i;
    int count = 0;
    Vector sum;
    float k;

    vset(sum, 0.0, 0.0, 0.0);

    for (i = 0; i < v->nverts; i++) {
        vadd(v->verts[i]->coord, sum, sum);
        count++;
    }

    if (count != 0) {
        k = 1.0 / count;
        new_pos[X] = k * sum[X];
        new_pos[Y] = k * sum[Y];
        new_pos[Z] = k * sum[Z];
    } else
        vcopy(v->coord, new_pos);
}


/******************************************************************************
Smooth all the vertex positions in a mesh by averaging the nearby neighbor
positions.

Entry:
  sc - scan containing mesh
******************************************************************************/

smooth_vertices(sc)
Scan* sc;
{
    int i;
    Mesh* mesh;
    Vertex* v;
    Vector new_pos;

    mesh = sc->meshes[mesh_level];

    /* find the smoothed position of all vertices in the mesh */

    for (i = 0; i < mesh->nverts; i++) {

        v = mesh->verts[i];

        /* skip vertices on the edge of the mesh */
        if (v->on_edge) {
            vcopy(v->coord, v->normal);
        } else {
            compute_smoothing(v, new_pos);
            vcopy(new_pos, v->normal);  /* oh no! bad place to save it!*/
        }
    }

    /* change the positions */

    for (i = 0; i < mesh->nverts; i++) {
        v = mesh->verts[i];
        remove_from_hash(v, mesh);
        vcopy(v->normal, v->coord);
        add_to_hash(v, mesh);
    }

    /* fix the triangle geometry */

    for (i = 0; i < mesh->ntris; i++)
        set_triangle_geometry(mesh->tris[i]);

    /* find the correct normals */

    for (i = 0; i < mesh->nverts; i++)
        find_vertex_normal(mesh->verts[i]);
}

