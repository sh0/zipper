/*

Keep track of edges of a mesh.

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
#include "malloc.h"
#include "cyfile.h"
#include "zipper.h"
#include "zipglobal.h"
#include "matrix.h"

/* define either VERBOSE_EDGES or NO_VERBOSE_EDGES */

#define NO_VERBOSE_EDGES

extern int bound_flag;


/******************************************************************************
See that the boundary edges are drawn.  May have to create the boundary
edge list.
******************************************************************************/

void make_edge_proc()
{
    bound_flag = 1 - bound_flag;
}


/******************************************************************************
Create the edge list for a mesh.

Entry:
  mesh - mesh to create edges of
******************************************************************************/

create_edge_list(mesh)
Mesh* mesh;
{
    int i, j, k;
    Vertex* v1, *v2;
    Triangle* tri;
    int vcount[10];
    int abnormal;
    int num_adj;
    int nverts;
    static Vertex* near_verts[100];

    /*** do this cleaner !!!! ***/
    /*** do this cleaner !!!! ***/
    /*** do this cleaner !!!! ***/
    /*** do this cleaner !!!! ***/
    /* wipe out old description of edges */
    if (mesh->nedges != 0) {
        mesh->nedges = 0;
    }

    /* zero out each vertices list of edges */
    for (i = 0; i < mesh->nverts; i++) {
        mesh->verts[i]->nedges = 0;
    }
    /*** do this cleaner !!!! ***/
    /*** do this cleaner !!!! ***/
    /*** do this cleaner !!!! ***/
    /*** do this cleaner !!!! ***/


    /* delete all abnormal vertices and the associated triangles */

    for (k = mesh->nverts - 1; k >= 0; k--) {
        v1 = mesh->verts[k];

        /* ignore vertices not on the mesh edge */
        if (!v1->on_edge)
            continue;

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

        /* examine the counts of the neighboring vertices to see if they */
        /* tell if the vertex has more than 3 edges meeting or if some edges */
        /* are triple edges */

        abnormal = 0;
        num_adj = 0;
        for (i = 0; i < v1->nverts; i++) {
            v2 = v1->verts[i];
            /* see if this vertex has and edge shared by too many triangles */
            if (v2->count != 1 && v2->count != 2) {
                abnormal = 1;
                break;
            }
            /* count the number of edges meeting at the vertex */
            if (v2->count == 1) {
                num_adj++;
            }
        }

        /* delete the triangles from any abnormal vertex */
        if (abnormal || num_adj > 2) {

#ifdef VERBOSE

            if (abnormal)
                printf("(getting rid of abnormal vertex)\n");
            else
                printf("(getting rid of vertex on a triple edge)\n");

#endif

            /* save a list of adjacent vertices */
            nverts = v1->nverts;
            for (i = 0; i < nverts; i++)
                near_verts[i] = v1->verts[i];

            /* delete the triangles */
            for (i = v1->ntris - 1; i >= 0; i--)
                delete_triangle(v1->tris[i], mesh, 0);

            /* check to see if the nearby vertices are now on an edge */
            for (i = 0; i < nverts; i++)
                vertex_edge_test(near_verts[i]);
        }
    }

    /* delete any triangles that now have no associated triangles */
    remove_unused_verts(mesh);



    /* find all the edges that belong in the list */

    for (k = mesh->nverts - 1; k >= 0; k--) {
        v1 = mesh->verts[k];

        /* ignore vertices not on the mesh edge */
        if (!v1->on_edge)
            continue;

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

        /* examine the counts of the neighboring vertices to see if they */
        /* prove the vertex to be an edge.  also check for mesh consistancy */

        for (i = 0; i < v1->nverts; i++) {
            /* add a new edge if this vertex is on the mesh edge */
            /* (but take care to add edge only once) */
            v2 = v1->verts[i];
            if (v2->count == 1 && (size_t) v1 < (size_t) v2) {
                add_edge_to_mesh(mesh, v1, v2);
            }
        }
    }

    /* say the edges of this mesh are valid */
    mesh->edges_valid = 1;

#if 0

    /* examine number of edges per vertex */

    printf("\n");

    for (i = 0; i < 10; i++)
        vcount[i] = 0;

    for (i = 0; i < mesh->nverts; i++) {
        v1 = mesh->verts[i];
        vcount[v1->nedges]++;
    }

    for (i = 0; i < 10; i++)
        if (vcount[i])
            printf("%d edges: %d\n", i, vcount[i]);

#endif

    /* create edge loop list */
    make_edge_loops(mesh);
}


/******************************************************************************
Link together edges on the boundary into edge loops.

Entry:
  mesh - mesh to make loops for
******************************************************************************/

make_edge_loops(mesh)
Mesh* mesh;
{
    int i, j;
    EdgeLoop* list;
    Vertex* vert;
    Edge* e;
    Vector pos;

    list = &mesh->looplist;
    list->nloops = 0;

    /* use "count" to say how many of the edges of a vertex have been */
    /* added to a loop list */

    for (i = 0; i < mesh->nverts; i++)
        mesh->verts[i]->count = 0;

    /* go through each vertex to see if it has edges with which to create */
    /* an edge loop */

    for (i = 0; i < mesh->nverts; i++) {
        vert = mesh->verts[i];
        if (vert->count < vert->nedges) {
            int found;
            /* find an un-used edge */
            found = 0;
            for (j = 0; j < vert->nedges; j++) {
                e = vert->edges[j];
                if (!e->used && e->v1 == vert) {
                    found = 1;
                    break;
                }
            }
            /* check to make sure we found an un-used edge */
            if (!found) {
                fprintf(stderr, "make_edge_loops: can't find un-used edge\n");

#if 1
                init_extra_lines();
                mesh_to_world(mesh->parent_scan, vert->coord, pos);
                add_extra_line(pos, pos, 0x00ff00);
#endif

                printf("vert = %d, count = %d, nedges = %d\n",
                       vert->index, vert->count, vert->nedges);

                for (j = 0; j < vert->nedges; j++) {
                    e = vert->edges[j];
                    printf("v1 v2 used: %d %d %d\n", e->v1->index, e->v2->index, e->used);
                }

                for (j = 0; j < vert->nverts; j++)
                    printf("v: %d\n", vert->verts[j]->index);

                for (j = 0; j < vert->ntris; j++)
                    printf("t: %d %d %d\n", vert->tris[j]->verts[0]->index,
                           vert->tris[j]->verts[1]->index,
                           vert->tris[j]->verts[2]->index);

                for (j = 0; j < vert->ntris; j++)
                    vert->tris[j]->mark = 3;

                /*
                exit (-1);
                */
                continue;
            }
            /* say we've used another edge of both vertices of this edge */
            e->v1->count++;
            e->v2->count++;
            /* this edge will start an edge loop */
            e->used = 1;
            if (list->nloops < LOOP_MAX) {
                list->loops[list->nloops++] = e;
            } else {
                fprintf(stderr, "make_edge_loops: not enough room for another loop\n");
                exit(-1);
            }
            follow_edges(mesh, e, list->nloops - 1);
        }
    }

#ifdef VERBOSE_EDGES
    printf("%d loops\n", list->nloops);
#endif
}


/******************************************************************************
Follow the edges around to form a loop of edges.

Entry:
  mesh   - mesh to add new loop to
  e_orig - edge to start the loop with
  num    - number of this loop
******************************************************************************/

follow_edges(mesh, e_orig, num)
Mesh* mesh;
Edge* e_orig;
int num;
{
    int i, j;
    Vertex* vert;
    Edge* e;
    Edge* e_next;
    Vertex* v_next;

    e = e_orig;
    e->num = num;
    vert = e->v2;

    /* keep looking for new edges for loop until we've come full circle */

    while (vert->count < vert->nedges) {

        /* look for next edge to add to loop */

        e_next = NULL;
        for (i = 0; i < vert->nedges; i++) {
            e_next = vert->edges[i];
            if (!e_next->used && e_next->v1 == vert) {
                break;
            }
        }

        /* see if we got an edge */

        if (e_next == NULL) {
            fprintf(stderr, "follow_edges: can't find another edge\n");
            exit(-1);
        }

        /* add this edge to the loop */

        e->next = e_next;
        e = e_next;
        e->used = 1;
        e->num = num;

        e->v1->count++;
        e->v2->count++;

        vert = e->v2;
    }

    /* close up loop */
    e->next = e_orig;

    /* create the backwards pointers */

    e = e_orig;

    do {
        e_next = e->next;
        e_next->prev = e;
        e = e_next;
        i++;
    } while (e != e_orig);
}


/******************************************************************************
Swap the order of the vertices in an edge.

Entry:
  e - edge in which to swap vertices
******************************************************************************/

swap_verts_in_edge(e)
Edge* e;
{
    Vertex* vert;

    vert = e->v1;
    e->v1 = e->v2;
    e->v2 = vert;
}


/******************************************************************************
Add an edge to the list of edges of a mesh.

Entry:
  mesh  - mesh to add edges to
  v1,v2 - vertices of the edge to add
******************************************************************************/

add_edge_to_mesh(mesh, v1, v2)
Mesh* mesh;
Vertex* v1, *v2;
{
    Edge* e;
    int i, j;
    int found;
    Triangle* tri;
    Vertex* vert;

    /* see if there is room for another edge */

    if (mesh->nedges >= mesh->max_edges) {
        mesh->max_edges += 40;
        mesh->edges = (Edge**)
                      realloc(mesh->edges, sizeof(Edge*) * mesh->max_edges);
    }

    /* create edge and place it on edge list of the mesh */

    e = (Edge*) myalloc(sizeof(Edge));
    e->v1 = v1;
    e->v2 = v2;
    e->used = 0;
    e->cuts = NULL;
    mesh->edges[mesh->nedges++] = e;

    /* add edge to each vertices list of edges */

    if (v1->nedges >= v1->max_edges) {
        v1->max_edges += 2;
        if (v1->edges == NULL)
            v1->edges = (Edge**) myalloc(sizeof(Edge*) * v1->max_edges);
        else
            v1->edges = (Edge**)
                        realloc(v1->edges, sizeof(Edge*) * v1->max_edges);
    }
    v1->edges[v1->nedges++] = e;

    if (v2->nedges >= v2->max_edges) {
        v2->max_edges += 2;
        if (v2->edges == NULL)
            v2->edges = (Edge**) myalloc(sizeof(Edge*) * v2->max_edges);
        else
            v2->edges = (Edge**)
                        realloc(v2->edges, sizeof(Edge*) * v2->max_edges);
    }
    v2->edges[v2->nedges++] = e;

    /* determine which triangle this edge belongs to */

    found = 0;
    for (i = 0; i < v1->ntris; i++) {
        for (j = 0; j < 3; j++) {
            tri = v1->tris[i];
            if (tri->verts[j] == v2) {
                found = 1;
                /* record which triangle the edge belongs to */
                e->tri = tri;
                /* make sure the vertices are in the proper order */
                if (tri->verts[(j + 1) % 3] == v1) {
                    vert = e->v1;
                    e->v1 = e->v2;
                    e->v2 = vert;
                } else if (tri->verts[(j + 2) % 3] == v1) {
                    ; /* do nothing */
                } else {
                    fprintf(stderr,
                            "add_edge_to_mesh: couldn't find vertex in triangle\n");
                }
                break;
            }
        }
        if (found)
            break;
    }

    /* error check */
    if (!found) {
        fprintf(stderr, "add_edge_to_mesh: couldn't find triangle for edge\n");
    }
}


/******************************************************************************
Zipper together a model.
******************************************************************************/

void new_zipper_proc()
{
    join_loops(scans[1], scans[0]);
}


/******************************************************************************
Join boundary loops from two meshes together.

Entry:
  sc1,sc2 - scans containing the two meshes
******************************************************************************/

join_loops(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    EdgeLoop* list1, *list2;
    Edge* e, *e_next, *e_orig;
    NearPosition near_info;
    Vector pos, norm;
    int result;
    Vertex* vert;
    Vertex* vnear;
    int inc;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    inc = level_to_inc(mesh_level);

    if (!m1->edges_valid)
        create_edge_list(m1);

    if (!m2->edges_valid)
        create_edge_list(m2);

    list1 = &m1->looplist;
    list2 = &m2->looplist;

    /* mark all vertices from both meshes as untouched */

    for (i = 0; i < m1->nverts; i++)
        m1->verts[i]->count = 0;

    for (i = 0; i < m2->nverts; i++)
        m2->verts[i]->count = 0;

    /* examine each loop from mesh 1 to see what loop of mesh 2 it should */
    /* be joined to */

    for (i = 0; i < list1->nloops; i++) {

        /*
        printf ("loop number %d:\n", i);
        */

        /* look through each vertex in an edge */

        e_orig = list1->loops[i];
        for (e = e_orig; e->next != e_orig; e = e->next) {

            vert = e->v1;
            if (vert->count == 0) {

                mesh_to_world(sc1, vert->coord, pos);
                mesh_to_world_normal(sc1, vert->normal, norm);
                result = nearest_on_mesh(sc2, m2, NULL, pos, norm,
                                         inc * ZIPPER_RESOLUTION, FIND_COS, &near_info);

                /* go on to next vertex if this one has already been matched */
                if (!result) {
                    vert->count = 1;
                    continue;
                }

                /* go on if the nearest vertex isn't part of a loop */
                /* or if it has already been used in joining two loops */
                vnear = near_info.v1;
                if (vnear->nedges == 0 || vnear->count != 0) {
                    vert->count = 1;
                    continue;
                }

                /* follow the adjacent loops */
                follow_loops(sc1, sc2, m1, m2, vert, vnear, e);
            }
        }
    }
}


/******************************************************************************
Follow around two loops to find how long the portion they should be zippered
is.

Entry:
  sc1,sc2 - scans of the two meshes
  m1,m2   - two meshes that the loops belong to
  v1,v2   - vertices from first and second mesh
  e1      - edge from first mesh
******************************************************************************/

follow_loops(sc1, sc2, m1, m2, v1, v2, e1)
Scan* sc1, *sc2;
Mesh* m1, *m2;
Vertex* v1, *v2;
Edge* e1;
{
    int i, j;
    EdgeLoop* list1, *list2;
    int num1, num2;
    Edge* e;
    int still_adjacent;
    int result;
    NearPosition near_info;
    Vector pos, norm;
    Vertex* vert, *vnear;
    Edge* e_min, *e_max;
    int count;
    int inc;

    inc = level_to_inc(mesh_level);

    list1 = &m1->looplist;
    list2 = &m2->looplist;

    num1 = v1->edges[0]->num;
    num2 = v2->edges[0]->num;

    /* follow the first loop down the "next" chain */

    e = e1;
    e_max = e;
    still_adjacent = 1;

    while (still_adjacent) {
        vert = e->v1;
        mesh_to_world(sc1, vert->coord, pos);
        mesh_to_world_normal(sc1, vert->normal, norm);
        result = nearest_on_mesh(sc2, m2, NULL, pos, norm,
                                 inc * ZIPPER_RESOLUTION, FIND_COS, &near_info);
        vnear = near_info.v1;
        if (result == 0 || vnear->count != 0 || vnear->nedges == 0 ||
            vnear->edges[0]->num != num2) {
            still_adjacent = 0;
        } else {
            e_max = e;
            e = e->next;
        }
    }

    /* follow the first loop up the "prev" chain */

    e = e1;
    e_min = e;
    still_adjacent = 1;

    while (still_adjacent) {
        vert = e->v1;
        mesh_to_world(sc1, vert->coord, pos);
        mesh_to_world_normal(sc1, vert->normal, norm);
        result = nearest_on_mesh(sc2, m2, NULL, pos, norm,
                                 inc * ZIPPER_RESOLUTION, FIND_COS, &near_info);
        vnear = near_info.v1;
        if (result == 0 || vnear->count != 0 || vnear->nedges == 0 ||
            vnear->edges[0]->num != num2) {
            still_adjacent = 0;
        } else {
            e_min = e;
            e = e->prev;
        }
    }

    /* count number of links in chain */

    printf("loop: ");

    count = 0;
    for (e = e_min; e != e_max; e = e->next) {
        count++;

        printf("%d ", e->v1->index);

#if 1
        e->v1->count = 1;
#endif

    }

#if 1
    e_max->v1->count = 1;
#endif

    printf("%d ", e_max->v1->index);
    printf("\n");
    printf("chain count = %d\n", count);
}

