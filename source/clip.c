/*

Clip the triangles of one mesh to the edge of another mesh.

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

#include "zipper.h"
#include "matrix.h"

extern float point_project_line();

/* set of points near the edge of a mesh */
static Vertex** pts_near = NULL;
static int pts_near_num;
static int pts_near_max;

/* set of edges near the edge of a mesh */
static Edge** edges_near = NULL;
static int edges_near_num;
static int edges_near_max;

#define MESH_A    1
#define MESH_B    2
#define CUT       3
#define USED_CUT  4

static float CLIP_NEAR_DIST_FACTOR;
static float CLIP_NEAR_DIST;
static float CLIP_NEAR_COS;
static float CLIP_BOUNDARY_DIST_FACTOR;
static float CLIP_BOUNDARY_DIST;
static float CLIP_BOUNDARY_COS;

Clip_List* split_list(Clip_List* clist, int index1, int index2);

update_clip_resolution()
{
    CLIP_NEAR_DIST = ZIPPER_RESOLUTION * CLIP_NEAR_DIST_FACTOR;
    CLIP_BOUNDARY_DIST = ZIPPER_RESOLUTION * CLIP_BOUNDARY_DIST_FACTOR;
}


set_clip_near_dist_factor(factor)
float factor;
{
    CLIP_NEAR_DIST_FACTOR = factor;
    CLIP_NEAR_DIST = ZIPPER_RESOLUTION * CLIP_NEAR_DIST_FACTOR;
}


float
get_clip_near_dist_factor()
{
    return CLIP_NEAR_DIST_FACTOR;
}


set_clip_near_cos(cosine)
float cosine;
{
    CLIP_NEAR_COS = cosine;
}


float
get_clip_near_cos()
{
    return CLIP_NEAR_COS;
}


set_clip_boundary_dist_factor(factor)
float factor;
{
    CLIP_BOUNDARY_DIST_FACTOR = factor;
    CLIP_BOUNDARY_DIST = ZIPPER_RESOLUTION * CLIP_BOUNDARY_DIST_FACTOR;
}


float
get_clip_boundary_dist_factor()
{
    return CLIP_BOUNDARY_DIST_FACTOR;
}


set_clip_boundary_cos(cosine)
float cosine;
{
    CLIP_BOUNDARY_COS = cosine;
}


float
get_clip_boundary_cos()
{
    return CLIP_BOUNDARY_COS;
}




/******************************************************************************
Clip one set of triangles to the edges of another.  Actually, the triangles
have already been gathered into one mesh by gather_triangles().

Entry:
  sc1 - first mesh containing edge to clip to
  sc2 - second mesh, containing triangles to clip
******************************************************************************/

clip_triangles(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    Vertex* vert;
    Triangle* tri;
    Vector c1, c2, c3;
    Vector norm1, norm2, norm3;
    int r1, r2, r3;
    NearPosition n1, n2, n3;
    int count;
    float max_length;
    float time;
    extern float time_it();
    extern float edge_length_max();
    int inc;

#ifdef DEBUG_CLIP
    printf("start of clip_triangles\n");
#endif

    time = time_it();

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    inc = level_to_inc(mesh_level);

    /* mark all triangles that came from mesh 2 that may need clipping */

    for (i = 0; i < m1->ntris; i++) {

        tri = m1->tris[i];

        /* assume at first that this triangle shouldn't be marked for clipping */
        tri->mark = 0;

        /* only look at triangles that used to be in the second mesh */
        if (tri->verts[0]->old_mesh != m2)
            continue;

#if 1
        /* don't look at triangles that were marked as "don't touch" from */
        /* the intersection routines */

        if (tri->dont_touch)
            continue;
#endif

#if 0
        /* only look at triangles that are on the mesh edge */
        r1 = tri->verts[0]->on_edge;
        r2 = tri->verts[1]->on_edge;
        r3 = tri->verts[2]->on_edge;
        if (r1 + r2 + r3 == 0)
            continue;
#endif

        /* transform vertices of triangle to world coordinates */
        mesh_to_world(sc1, tri->verts[0]->coord, c1);
        mesh_to_world(sc1, tri->verts[1]->coord, c2);
        mesh_to_world(sc1, tri->verts[2]->coord, c3);

        mesh_to_world_normal(sc1, tri->verts[0]->normal, norm1);
        mesh_to_world_normal(sc1, tri->verts[1]->normal, norm2);
        mesh_to_world_normal(sc1, tri->verts[2]->normal, norm3);

        /* find nearest places on the other mesh */
        r1 = nearest_on_mesh(sc1, m1, m2, c1, norm1,
                             inc * CLIP_NEAR_DIST, CLIP_NEAR_COS, &n1);
        r2 = nearest_on_mesh(sc1, m1, m2, c2, norm2,
                             inc * CLIP_NEAR_DIST, CLIP_NEAR_COS, &n2);
        r3 = nearest_on_mesh(sc1, m1, m2, c3, norm3,
                             inc * CLIP_NEAR_DIST, CLIP_NEAR_COS, &n3);

        /* mark for clipping */
        if (r1 || r2 || r3)
            tri->mark = 1;

        /* see if these places are all valid and away from the edge of the mesh */
        r1 = r1 && (n1.on_edge == 0);
        r2 = r2 && (n2.on_edge == 0);
        r3 = r3 && (n3.on_edge == 0);

        /* save how many vertices are on the other mesh */
        tri->eat_mark = r1 + r2 + r3;
    }

    /* initialize the intersection list for edges on the boundary */
    init_cuts(sc1, m1);

    /* make triangles perpendicular to the edge that will be used for clipping */
    make_clip_triangles(sc1, m1);

    /*
    return;
    */

    /* mark all vertices as untouched */
    /* (count will be marker saying that a vertex is in pts_near) */

    for (i = 0; i < m1->nverts; i++)
        m1->verts[i]->count = 0;

    /* examine each marked triangle in turn and clip them */

    init_extra_lines();
    max_length = edge_length_max(mesh_level);

    for (i = 0; i < m1->ntris; i++) {

        tri = m1->tris[i];

        /* only examine marked triangles */
        if (tri->mark == 0)
            continue;

        /* find nearby vertices that are on the mesh edge */
        pts_near_num = 0;
        for (j = 0; j < 3; j++)
            verts_near_edges(m1, m2, tri->verts[j]->coord, tri->verts[j]->normal,
                             max_length, CLIP_NEAR_COS);

        /* un-mark all the vertices in pts_near */
        for (j = 0; j < pts_near_num; j++)
            pts_near[j]->count = 0;

#if 0
        /* (DEBUGGING) */
        /* (DEBUGGING) */
        /* mark for drawing those strange triangles that have no nearby verts */
        /* (DEBUGGING) */
        /* (DEBUGGING) */
        if (pts_near_num == 0)
            tri->mark++;
#endif

        /* find which edges are near the current triangle */
        edges_near_edges(tri, m1, m2, sc1);

        /* see which edges in edges_near clip the current triangle */
        cut_triangle(tri, m1, m2, sc1);

#if 0
        /* see which edges in pts_near cut across the current triangle */
        old_cut_triangle(tri, m1, m2, sc1);
#endif

    }

    /* create new vertices at all the intersection points */
    create_cut_vertices(m1);

    /* clip all triangles */
    perform_triangle_clipping(sc1, sc2);

    /* replace old triangles with new ones that incorporate the cut points */
    introduce_all_cuts(m1, m2);

    /* make sure all vertices now belong to the first mesh */
    for (i = 0; i < m1->nverts; i++)
        m1->verts[i]->old_mesh = m1;

    /* clean up info in triangles */
    for (i = 0; i < m1->ntris; i++)
        tri = m1->tris[i];
    if (tri->dont_touch == 0 && tri->clips) {
        free(tri->clips->cuts);
        free(tri->clips);
        tri->clips = NULL;
    }

    /* remove un-used vertices */
    remove_unused_verts(m1);

#if 0
    /* the extra lines give (roughly) how many cut vertices were introduced */
    print_extra_lines();
#endif

#ifdef DEBUG_CLIP
    printf("end of clip_triangles\n");
#endif

    time = time_it();
    printf("time: %.3f\n", time);
}


/******************************************************************************
Actually do the triangle clipping, now that we have all the information
about intersections between triangles and mesh edges.

Entry:
  sc1 - first mesh containing edge to clip to
  sc2 - second mesh, containing triangles to clip
******************************************************************************/

perform_triangle_clipping(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j, k;
    Mesh* m1, *m2;
    Triangle* tri;
    int count;
    Clip_List* potential_vertices(), *split_list();
    Clip_List* clist;
    Vertex* v1, *v2, *v3;
    int result;
    int tri_index;
    Vector tri_norm;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    /* find all triangles that have been clipped */
    for (i = m1->ntris - 1; i >= 0; i--) {

        tri = m1->tris[i];

        /* only examine marked triangles */
        if (tri->mark == 0)
            continue;

        /* triangles that haven't been cut are a special case */
        if (tri->clips == NULL) {
            /* delete un-cut triangles that are sitting on the clip mesh */
            if (tri->eat_mark > 1) {
                v1 = tri->verts[0];
                v2 = tri->verts[1];
                v3 = tri->verts[2];
                delete_triangle(tri, m1, 1);
                vertex_edge_test(v1);
                vertex_edge_test(v2);
                vertex_edge_test(v3);
                find_vertex_normal(v1);
                find_vertex_normal(v2);
                find_vertex_normal(v3);
            }
            continue;
        }

        /* count up all the places the triangle has been cut */
        count = 0;
        for (j = 0; j < 3; j++)
            count += tri->clips[j].cut_num;

        /* sort the cuts along a triangle's edges */
        sort_triangle_cuts(tri);

        /* create list of potential vertices from triangle to be clipped */
        clist = potential_vertices(tri);

        /* determine correspondence between partner cuts (two cuts that */
        /* slice across a triangle */
        result = find_partner_cuts(tri, clist, m1);

        /* clip the triangle if we successfully found partners to all cuts */

        if (result) {

            /* remember old triangle index and normal */
            tri_index = tri->index;
            tri_norm[X] = -tri->aa;
            tri_norm[Y] = -tri->bb;
            tri_norm[Z] = -tri->cc;

            /* delete the old triangle */
            delete_triangle(tri, m1, 0);

            /* process the list by either dividing it into smaller lists or */
            /* by creating triangles that are bordered by the vertices in the list */
            process_vertices(tri_norm, tri_index, clist, m1);

            /* re-compute vertex info around the old triangle */
            for (j = 0; j < clist->count; j++) {
                vertex_edge_test(clist->list[j].vert);
                find_vertex_normal(clist->list[j].vert);
            }
        } else {

#if 0
            if (count % 2 == 1)
                tri->mark++;     /* mark odd-numbered cuts */
            else
                tri->mark += 2;  /* mark other problem triangles */
#endif

#if 1
            /* if we get here, we couldn't figure out how to clip the */
            /* triangle, so we'll delete it */

            /* delete the old triangle */
            delete_triangle(tri, m1, 0);

            /* re-compute vertex info around the old triangle */
            for (j = 0; j < clist->count; j++) {
                vertex_edge_test(clist->list[j].vert);
                find_vertex_normal(clist->list[j].vert);
            }
#endif

        }

        /* free up space from the list */
        free(clist->list);
        free(clist);
    }
}


/******************************************************************************
Take a list of vertices from a triangle being clipped.  If the list requires
it, further divide the list at partnered cut vertices.  Otherwise, create
new triangles that are bordered by the list.

Entry:
  tnorm  - normal of triangle being clipped
  tindex - normal of triangle being clipped
  clist  - the list of clip vertices to process
  mesh   - mesh that triangle belongs to
******************************************************************************/

process_vertices(tnorm, tindex, clist, mesh)
Vector tnorm;
int tindex;
Clip_List* clist;
Mesh* mesh;
{
    int i, j;
    int index1, index2;
    Clip_List* clist2;

    /* find where to split the list */
    index1 = index2 = -1;
    for (j = 0; j < clist->count; j++)
        if (clist->list[j].type == CUT) {
            index1 = j;
            break;
        }

    /* if we don't find any CUT vertices, we're ready to turn the list */
    /* into a group of triangles */

#if 0
    print_list("orig", clist);
#endif

    if (index1 == -1) {

        /* create new triangles */
        if (outside_mesh(clist))
            new_list_to_tris(tnorm, tindex, clist, mesh);

        return;
    }

    /* if we get here, we've got to look for a matching vertex */

    for (j = index1 + 1; j < clist->count; j++)
        if (clist->list[j].type == CUT &&
            clist->list[j].cut->partner == clist->list[index1].cut) {
            index2 = j;
            break;
        }
    if (index2 == -1) {
#ifdef DEBUG_CLIP
        fprintf(stderr, "process_vertices: can't find matching cut\n");
#endif
        return;
    }

    /* split the list across the partnered cut vertices */
    clist2 = split_list(clist, index1, index2);

#if 0
    print_list("new1", clist);
    print_list("new2", clist2);
    printf("\n");
#endif

    /* process the two newly created vertex lists by making a recursive */
    /* call to process_vertices() */
    process_vertices(tnorm, tindex, clist, mesh);
    process_vertices(tnorm, tindex, clist2, mesh);

    /* free up memory */
    free(clist2->list);
    free(clist2);
}


/******************************************************************************
Test whether the vertices in a list lie outside of the mesh.

Entry:
  clist - vertex list to check

Exit:
  returns 1 if outside, 0 if inside
******************************************************************************/

int outside_mesh(clist)
Clip_List* clist;
{
    int j;
    int index1, index2;

    /* find a used cut */
    index1 = index2 = -1;
    for (j = 0; j < clist->count; j++)
        if (clist->list[j].type == USED_CUT) {
            index1 = j;
            break;
        }
    if (index1 == -1) {
        fprintf(stderr, "outside_mesh: cut inconsistancy\n");
        exit(-1);
    }

    /* find the first used cut's partner */
    for (j = index1 + 1; j < clist->count; j++)
        if (clist->list[j].type == USED_CUT &&
            clist->list[j].cut->partner == clist->list[index1].cut) {
            index2 = j;
            break;
        }
    if (index2 == -1) {
        fprintf(stderr, "outside_mesh: other cut inconsistancy\n");
        exit(-1);
    }


#if 0

    /* search for the first used cut */
    for (i = 0; i < ct; i++)
        if (clist->list[i].type == CUT) {
            index1 = i;
            break;
        }

    /* see if the second cut is before or after this one */
    if (clist->list[(index1 + 1) % ct].type == CUT) {
        index2 = (index1 + 1) % ct;
    } else if (clist->list[(index1 + ct - 1) % ct].type == CUT) {
        index2 = (index1 + ct - 1) % ct;
    } else {
#ifdef DEBUG_CLIP
        fprintf(stderr, "outside_mesh: can't find second cut\n");
#endif
        return (0);
    }

#endif


    /* see which side of mesh the vertex list is on */
    if (clist->list[index1].inward == 0 && clist->list[index2].inward == 1)
        return (1);
    else if (clist->list[index1].inward == 1 && clist->list[index2].inward == 0)
        return (0);
    else {
#ifdef DEBUG_CLIP
        fprintf(stderr, "outside_mesh: cuts have same direction\n");
#endif
        return (0);
    }
}


/******************************************************************************
Creat new triangles based on a list of vertices.

Entry:
  tri   - triangle that is currently being clipped
  clist - list of vertices to create triangles from
  mesh  - mesh that triangle should be added to
******************************************************************************/

list_to_tris(tri, clist, mesh)
Triangle* tri;
Clip_List* clist;
Mesh* mesh;
{
    int i, j;
    Clip_Vertex* v = clist->list;
    Triangle* new_tri;

    /* error check */
    if (clist->count < 3) {
#ifdef DEBUG_CLIP
        fprintf(stderr, "list_to_tris: too few vertices (%d)\n", clist->count);
#endif
        return;
    }

    /* make triangles that radiate from one particular vertex */

    for (i = 1; i < clist->count - 1; i++) {
        if (check_proposed_tri(v[0].vert, v[i].vert, v[i + 1].vert)) {
            new_tri = make_triangle(mesh, v[0].vert, v[i].vert, v[i + 1].vert, 1e20);

            /* check_proposed_tri() is inadequate!  - B. Curless  9/14/95 */
            if (new_tri != NULL)
                new_tri->mark = 1;

        }
    }

    /* re-compute the vertex normals and mark if they are on the mesh edge */
    for (i = 0; i < clist->count; i++) {
        vertex_edge_test(v[i].vert);
        find_vertex_normal(v[i].vert);
    }
}


/******************************************************************************
Creat new triangles based on a list of vertices.

Entry:
  tnorm - normal of triangle that is currently being clipped
  index - index of triangle being clipped
  clist - list of vertices to create triangles from
  mesh  - mesh that triangle should be added to
******************************************************************************/

new_list_to_tris(tnorm, index, clist, mesh)
Vector tnorm;
int index;
Clip_List* clist;
Mesh* mesh;
{
    int i, j;
    Clip_Vertex* v = clist->list;
    Vector vec;
    Triangle* new_tri;
    int result;
    int self_intersect;
    int p1, p2, p3;
    Vector norm;
    float n4;

    /* see if the list of vertices around the area to cover is */
    /* properly oriented */

    for (i = 0; i < clist->count; i++) {
        result = check_proposed_edge(v[i].vert, v[(i + 1) % clist->count].vert);
        if (result == 0) {
#ifdef DEBUG_CLIP
            printf("new_list_to_tris: list improperly oriented\n");
#endif
            return;
        }
    }

    /* maybe don't re-do triangle if it has too many cut points */
    if (clist->count > 10) {
#ifdef DEBUG_CLIP
        printf("clist->count = %d for tri %d\n", clist->count, index);
#endif
        return;
    }

    /* error check */
    if (clist->count < 3) {
#ifdef DEBUG_CLIP
        fprintf(stderr, "list_to_tris: too few vertices (%d)\n", clist->count);
#endif
        return;
    }

    /* initialize the polygon splitter */
    result = init_splitter(tnorm[X], tnorm[Y], tnorm[Z], 0.0);

    /* go ahead if we didn't get error signal from splitter initialization */

    if (!result) {

        /* transform the face's vertices to the xy-plane */
        /* and send them to the splitter */

        for (i = 0; i < clist->count; i++) {
            vcopy(v[i].vert->coord, vec);
            add_boundary_point(vec[X], vec[Y], vec[Z], i);
        }

        /* call the splitter */
        self_intersect = greedy_connect();
    } else {
#ifdef DEBUG_CLIP
        printf("new_list_to_tris: bad return from init_splitter\n");
#endif
    }

    if (self_intersect || result) {

#ifdef DEBUG_CLIP
        printf("self intersection\n");
#endif

        /* go to ugly set of triangles */
        for (i = 1; i < clist->count - 1; i++) {

            /* look out for flipped triangles */
            plane_thru_vectors(v[0].vert->coord, v[i].vert->coord,
                               v[i + 1].vert->coord, &norm[X], &norm[Y], &norm[Z], &n4);
            if (vdot(norm, tnorm) > 0) {
#ifdef DEBUG_CLIP
                printf("normal flipped for new tri made from tri %d.\n", index);
#endif
                continue;
            }

            if (check_proposed_tri(v[0].vert, v[i].vert, v[i + 1].vert)) {
                new_tri = make_triangle(mesh, v[0].vert, v[i].vert, v[i + 1].vert, 1e20);

                /* check_proposed_tri() is inadequate! - B. Curless  9/14/95  */
                if (new_tri != NULL)
                    new_tri->mark = 1;
            }
        }
    } else {
        /* get the newly-formed triangles from the splitting routines */
        for (i = 0; i < get_ntris(); i++) {

            get_triangle(i, &p1, &p2, &p3);

            /* look out for flipped triangles */
            plane_thru_vectors(v[p1].vert->coord, v[p2].vert->coord,
                               v[p3].vert->coord, &norm[X], &norm[Y], &norm[Z], &n4);
            if (vdot(norm, tnorm) > 0) {
#ifdef DEBUG_CLIP
                printf("normal flipped for new tri made from tri %d.\n", index);
#endif
                continue;
            }

            if (check_proposed_tri(v[p1].vert, v[p2].vert, v[p3].vert)) {
                new_tri = make_triangle(mesh, v[p1].vert, v[p2].vert, v[p3].vert, 1e20);

                /* check_proposed_tri() is inadequate! - B. Curless  9/14/95  */
                if (new_tri != NULL)
                    new_tri->mark = 1;
            }

        }
    }

    /* re-compute the vertex normals and mark if they are on the mesh edge */
    for (i = 0; i < clist->count; i++) {
        vertex_edge_test(v[i].vert);
        find_vertex_normal(v[i].vert);
    }
}


/******************************************************************************
Split a list of vertices into two lists.

Entry:
  clist - original list to spit
  index1 - index of the cut vertex at which to begin the split
  index2 - where to end the split

Exit:
  clist - one of the created lists
  returns the other new list
******************************************************************************/

Clip_List* split_list(clist, index1, index2)
Clip_List* clist;
int index1, index2;
{
    int i;
    int num;
    Clip_List* clist2, *old_list;
    Clip_List* inter_list;
    Clip_List* make_between_list();
    int forward;

    /* possibly swap the order of the indices */
    if (index1 > index2) {
        int temp;
        temp = index1;
        index1 = index2;
        index2 = temp;
    }

    /* mark the two indexed cuts as being used */
    clist->list[index1].type = USED_CUT;
    clist->list[index2].type = USED_CUT;

    /* create list of vertices along loop edge that are between the two cuts */
    inter_list = make_between_list(clist,
                                   clist->list[index1].cut, clist->list[index2].cut);

    /* determine which cut is first along the loop edge */
    if (clist->list[index1].cut->first)
        forward = 0;
    else if (clist->list[index2].cut->first)
        forward = 1;
    else {
        fprintf(stderr, "split_list: can't find which cut is first\n");
        exit(-1);
    }

    /* make a copy of the original list */
    old_list = (Clip_List*) malloc(sizeof(Clip_List));
    old_list->list = (Clip_Vertex*) malloc(sizeof(Clip_Vertex) * clist->count);
    old_list->count = clist->count;
    for (i = 0; i < old_list->count; i++)
        old_list->list[i] = clist->list[i];

    /*** create one of the new lists ***/

    clist2 = (Clip_List*) malloc(sizeof(Clip_List));
    clist2->count = index2 - index1 + 1 + inter_list->count;
    clist2->list = (Clip_Vertex*) malloc(sizeof(Clip_Vertex) * clist2->count);

    /* copy appropriate subset of old vertices */
    num = 0;
    for (i = index1; i <= index2; i++)
        clist2->list[num++] = old_list->list[i];

    /* add vertices along loop edge between cuts */
    if (forward)
        for (i = 0; i < inter_list->count; i++)
            clist2->list[num++] = inter_list->list[i];
    else
        for (i = inter_list->count - 1; i >= 0; i--)
            clist2->list[num++] = inter_list->list[i];

    /*** create the second list ***/

    /* copy appropriate subset of old vertices */
    num = 0;
    for (i = index2; i != index1; i = (i + 1) % old_list->count)
        clist->list[num++] = old_list->list[i];
    clist->list[num++] = old_list->list[index1];

    /* make sure there is enough room for inter_list vertices */
    if (num + inter_list->count > clist->count)
        clist->list = (Clip_Vertex*)
                      realloc(clist->list, sizeof(Clip_Vertex) * (num + inter_list->count));

    /* add vertices along loop edge between cuts */
    if (!forward)
        for (i = 0; i < inter_list->count; i++)
            clist->list[num++] = inter_list->list[i];
    else
        for (i = inter_list->count - 1; i >= 0; i--)
            clist->list[num++] = inter_list->list[i];
    clist->count = num;

    /* free up helper lists */
    free(old_list->list);
    free(old_list);
    free(inter_list->list);
    free(inter_list);

    /* return new list */
    return (clist2);
}


/******************************************************************************
Create a list of vertices that lie between two partner cuts along a loop edge.

Entry:
  clist - list of clip vertices around a triangle
  cut1  - one of the cuts
  cut2  - the other cut

Exit:
  returns a list of vertices that lie between cuts
******************************************************************************/

Clip_List* make_between_list(clist, cut1, cut2)
Clip_List* clist;
Cut* cut1, *cut2;
{
    Edge* edge, *edge2;
    Clip_List* inter_list;
    int inter_num, inter_max;

    /* create list of vertices along loop edge that are between the two cuts */
    inter_list = (Clip_List*) malloc(sizeof(Clip_List));
    inter_max = 4;
    inter_num = 0;
    inter_list->list = (Clip_Vertex*)
                       malloc(sizeof(Clip_Vertex) * inter_max);

    /* return an empty list if there are no vertices between the cuts */
    if (cut1->edge == cut2->edge) {
        inter_list->count = 0;
        return (inter_list);
    }

    /* see which cut comes first around the edge loop */

    if (cut1->first) {
        edge = cut1->edge;
        edge2 = cut2->edge;
    } else if (cut2->first) {
        edge = cut2->edge;
        edge2 = cut1->edge;
    } else {
        fprintf(stderr, "make_between_list: inconsistant cuts\n");
        exit(-1);
    }

    /* travel around the loop edge in one direction or the other, */
    /* recording all vertices between the two cuts */

    while (edge != edge2) {

        /* make sure there is room for new vertex in the list */
        if (inter_num == inter_max) {
            inter_max += 4;
            inter_list->list = (Clip_Vertex*)
                               realloc(inter_list->list, sizeof(Clip_Vertex) * inter_max);
        }

        /* add vertex to list */
        inter_list->list[inter_num].vert = edge->v2;
        inter_list->list[inter_num].type = MESH_A;
        inter_num++;

        /* go to next edge around the loop */
        edge = edge->next;
    }

    /* return the list */
    inter_list->count = inter_num;
    return (inter_list);
}


/******************************************************************************
Create the list of potential vertices for a clipped triangle.

Entry:
  tri - triangle being clipped

Exit:
  returns pointer to list of potential vertices
******************************************************************************/

Clip_List* potential_vertices(tri)
Triangle* tri;
{
    int i, j;
    int count;
    int num;
    Clip_List* clist;
    Clip_Vertex* list;
    Cut** cuts;

    /* create an empty list of potential vertices */

    count = tri->clips[0].cut_num + tri->clips[1].cut_num +
            tri->clips[2].cut_num + 3;
    clist = (Clip_List*) malloc(sizeof(Clip_List));
    clist->list = (Clip_Vertex*) malloc(sizeof(Clip_Vertex) * count);
    clist->count = count;
    list = clist->list;

    /* examine each edge of the triangle */

    num = 0;
    for (i = 0; i < 3; i++) {

        /* add old triangle vertex to list */
        list[num].type = MESH_B;
        list[num].vert = tri->verts[i];
        num++;

        /* add all the cuts from one triangle edge to the list */
        cuts = tri->clips[i].cuts;
        for (j = 0; j < tri->clips[i].cut_num; j++) {

            /* mark cut as having no partner */
            cuts[j]->partner = NULL;

            /* add to list */
            list[num].type = CUT;
            list[num].cut = cuts[j];
            list[num].vert = cuts[j]->new_vert;
            list[num].side_index = i;

            /* figure out if it is an inward or outward pointing cut */
            if (cuts[j]->v1 == tri->clips[i].v1) {
                if (cuts[j]->inward)
                    list[num].inward = 1;
                else
                    list[num].inward = 0;
            } else if (cuts[j]->v1 == tri->clips[i].v2) {
                if (cuts[j]->inward)
                    list[num].inward = 0;
                else
                    list[num].inward = 1;
            } else {
                fprintf(stderr, "potential vertices: wrong vertex correspondence\n");
                exit(-1);
            }
            num++;
        }
    }

    /* consistancy check */
    if (num != count) {
        fprintf(stderr, "potential vertices: wrong number of vertices\n");
        exit(-1);
    }

    /* return pointer to the list */
    return (clist);
}


/******************************************************************************
Figure out which cuts should be paired with each other because the two
of them slice across a triangle.

Entry:
  tri   - triangle whose cuts need to be partnered
  clist - the list of vertices that clip the triangle
  mesh  - mesh that triangle belongs to

Exit:
  returns 1 if successful, 0 if we got an error
******************************************************************************/

int find_partner_cuts(tri, clist, mesh)
Triangle* tri;
Clip_List* clist;
Mesh* mesh;
{
    int i, j, k;
    Cut* cut, *tcut;
    Clip_Vertex* list = clist->list;
    Clip_Vertex** loop_list;
    Clip_Vertex* clipvert;
    Clip_Vertex* cv1, *cv2;
    int nlist;
    Edge* edge, *fedge;
    int been_around;
    int num;
    int cut_index;
    Cut* next_similar_cut();

    /* allocate space for temporary cut list */
    loop_list = (Clip_Vertex**) malloc(sizeof(Clip_Vertex*) * clist->count);

    /* set all cuts to be initially un-partnered */
    for (i = 0; i < clist->count; i++)
        if (list[i].type == CUT)
            list[i].cut->partner = NULL;

    /* examine each cut of the clip list in turn to find it's partner */
    for (i = 0; i < clist->count; i++) {

        /* see if this cut is un-partnered */
        if (list[i].type == CUT && list[i].cut->partner == NULL) {

            /* collect all cuts that are around the same loop of edges */
            nlist = 0;
            cut = list[i].cut;
            for (j = 0; j < clist->count; j++)
                if (list[j].type == CUT && list[j].cut->edge->num == cut->edge->num)
                    loop_list[nlist++] = &list[j];

            /* we've got an error if there is an odd number of cuts in loop_list */
            if (nlist % 2 == 1) {
#ifdef DEBUG_CLIP
                fprintf(stderr, "find_partner_cuts: odd number for tri %d: %d\n",
                        tri->index, nlist);
#endif
                for (j = 0; j < nlist; j++) {
                    tcut = loop_list[j]->cut;
#ifdef DEBUG_CLIP
                    printf("s = %f, t = %f\n", tcut->s, tcut->t);
#endif
                    /*
                    printf ("cut %d, v1 = %d, v2 = %d\n", (int) tcut,
                      tcut->v1->index, tcut->v2->index);
                    printf ("edge: v1 = %d, v2 = %d\n",
                      (int) tcut->edge->v1->index,
                      (int) tcut->edge->v2->index);
                    printf ("location %f %f %f\n",
                               loop_list[j]->vert->coord[X],
                               loop_list[j]->vert->coord[Y],
                               loop_list[j]->vert->coord[Z]);
                    printf ("\n");
                    */
                }
                return (0);
            }

            /* travel around the loop and sort the order of the */
            /* cuts according to their position around the loop */
            fedge = mesh->looplist.loops[cut->edge->num];
            been_around = 0;
            num = 0;
            for (edge = fedge; edge != fedge || !been_around; edge = edge->next) {
                been_around = 1;
                for (j = 0; j < edge->cut_num; j++)
                    for (k = 0; k < nlist; k++)
                        if (edge->cuts[j] == loop_list[k]->cut) {
                            clipvert = loop_list[k];
                            loop_list[k] = loop_list[num];
                            loop_list[num] = clipvert;
                            num++;
                        }
            }

            /* consistancy check */
            if (num != nlist) {
#ifdef DEBUG_CLIP
                fprintf(stderr, "find_partner_cuts: couldn't sort list\n");
#endif
                return (0);
            }

            /* find two cuts that are adjacent around the triangle and that */
            /* are of the correct parity */

            cut_index = -1;
            for (j = 0; j < nlist; j++) {
                /* look at two cuts that are in order around the edge loop of mesh */
                cv1 = loop_list[j];
                cv2 = loop_list[(j + 1) % nlist];
                /* If they are correctly oriented, see if they are also adjacent */
                /* around the triangle being clipped.  (If they are adjacent on
                /* triangle, we now know that all other cuts through the same loop */
                /* can be chained together in pairs from this cut on.) */
                if (cv1->inward == 0 && cv2->inward == 1) {
                    /* check backwards around triangle */
                    tcut = next_similar_cut(cv1->cut, clist, 1);
                    if (tcut == cv2->cut) {
                        cut_index = j;
                        break;
                    } else {
                        /* if backwards didn't work, check forward */
                        tcut = next_similar_cut(cv1->cut, clist, 0);
                        if (tcut == cv2->cut) {
                            cut_index = j;
                            break;
                        }
                    }
                }
            }

            /* see if we actually got a correspondence */
            if (cut_index == -1) {

#ifdef DEBUG_CLIP

                fprintf(stderr, "find_partner_cuts: can't find correspondence\n");

                printf("triangle index = %d\n", tri->index);
                printf("nlist num loopnum: %d %d %d\n", nlist, num, cut->edge->num);

                printf("loop_list:\n");
                for (j = 0; j < nlist; j++)
                    printf("%d (%0d)\n", (int) loop_list[j]->cut, loop_list[j]->inward);

                printf("around triangle:\n");
                for (j = 0; j < clist->count; j++)
                    if (list[j].type == CUT)
                        printf("%d\n", (int) list[j].cut);

#endif

                return (0);
            }

            /* travel around the edge loop of the mesh from the cut-pair we */
            /* found above and partner the cuts in pairs */
            for (j = 0; j < nlist; j += 2) {
                k = (j + cut_index) % nlist;
                cut = loop_list[k]->cut;
                tcut = loop_list[(k + 1) % nlist]->cut;
                cut->partner = tcut;
                cut->first = 1;
                tcut->partner = cut;
                tcut->first = 0;
            }

        }
    }

    /* free the temporary list */
    free(loop_list);

    /* return value saying everything went okay */
    return (1);
}


/******************************************************************************
Find a cut that slices the same edge loop to a given cut that is the next cut
around the triangle.

Entry:
  cut   - cut to find similar cut to
  clist - the list of vertices that clip the triangle
  dir   - direction, 0 = forward, 1 = backwards

Exit:
  returns pointer to next similar cut, or NULL if there was no such cut
******************************************************************************/

Cut* next_similar_cut(cut, clist, dir)
Cut* cut;
Clip_List* clist;
int dir;
{
    int i, j;
    Clip_Vertex* list = clist->list;
    int index;

    /* go around the triangle trying to find the cut */
    index = -1;
    for (i = 0; i < clist->count; i++)
        if (list[i].type == CUT && list[i].cut == cut) {
            index = i;
            break;
        }

    /* sanity check */
    if (index == -1) {
        fprintf(stderr, "next_similar_cut: can't find the original cut\n");
        exit(-1);
    }

    /* find the next cut that cuts the same edge loop of the mesh */
    /* (actually, we travel backwards around the triangle's cuts) */

    for (i = 1; i < clist->count; i++) {
        if (dir == 1)
            j = (index - i + clist->count) % clist->count;    /* backwards */
        else
            j = (index + i) % clist->count;           /* forward */
        /* if this is a cut that cuts the same loop edge, return pointer to it */
        if (list[j].type == CUT && list[j].cut->edge->num == cut->edge->num)
            return (list[j].cut);
    }

    /* if we get here, we didn't find a correspondence */
    return (NULL);
}


/******************************************************************************
Print a list of vertices.

Entry:
  str   - string to print at head of list
  clist - list of vertices
******************************************************************************/

#if 0
print_list(str, clist)
char* str;
Clip_List* clist;
{
    int i;

    printf("%s: ", str);
    for (i = 0; i < clist->count; i++) {
        switch (clist->list[i].type) {
            case CUT:
                printf("C%0d ", clist->list[i].inward);
                break;
            case USED_CUT:
                printf("U%0d ", clist->list[i].inward);
                break;
            case MESH_A:
                printf("A ", (int) clist->list[i].vert);
                break;
            case MESH_B:
                printf("B ", (int) clist->list[i].vert);
                break;
            default:
                fprintf(stderr, "print_list: bad switch = %d\n", clist->list[i].type);
                break;
        }
    }
    printf("\n");
}
#endif

/******************************************************************************
Print a list of vertices.

Entry:
  str   - string to print at head of list
  clist - list of vertices
******************************************************************************/

#if 0
long_print_list(str, clist)
char* str;
Clip_List* clist;
{
    int i;

    printf("%s:\n", str);
    for (i = 0; i < clist->count; i++) {
        switch (clist->list[i].type) {
            case CUT:
                printf("C%0d: %f %f %f\n", clist->list[i].inward,
                       clist->list[i].vert->coord[X],
                       clist->list[i].vert->coord[Y],
                       clist->list[i].vert->coord[Z]);
                break;
            case USED_CUT:
                printf("U%0d: %f %f %f\n", clist->list[i].inward,
                       clist->list[i].vert->coord[X],
                       clist->list[i].vert->coord[Y],
                       clist->list[i].vert->coord[Z]);
                break;
            case MESH_A:
                printf("A: %f %f %f\n", (int) clist->list[i].vert,
                       clist->list[i].vert->coord[X],
                       clist->list[i].vert->coord[Y],
                       clist->list[i].vert->coord[Z]);
                break;
            case MESH_B:
                printf("B: %f %f %f\n", (int) clist->list[i].vert,
                       clist->list[i].vert->coord[X],
                       clist->list[i].vert->coord[Y],
                       clist->list[i].vert->coord[Z]);
                break;
            default:
                fprintf(stderr, "print_list: bad switch = %d\n", clist->list[i].type);
                break;
        }
    }
    printf("\n");
}
#endif

/******************************************************************************
Sort the points of intersection along each of a triangle's three edges.

Entry:
  tri - triangle whose cuts to sort
******************************************************************************/

sort_triangle_cuts(tri)
Triangle* tri;
{
    int i, j, k;
    int index;
    Vertex* vert;
    Cut* ctemp;
    Clipped_Edge* clip;
    float s_j, s_index;

    /* examine each of the triangle's edges */

    for (k = 0; k < 3; k++) {

        clip = &tri->clips[k];
        vert = tri->verts[k];

#if 0
        if (clip->cut_num > 1) {
            printf("before:  ");
            for (i = 0; i < clip->cut_num; i++)
                printf("%f ", clip->cuts[i]->s);
            printf("\n");
        }
#endif

        /* sort the cuts along this edge (selection sort) */

        for (i = clip->cut_num - 1; i > 0; i--) {

            /* find maximum for j in [0,i] */

            index = 0;
            if (clip->cuts[index]->v1 == vert)
                s_index = clip->cuts[index]->s;
            else {
                s_index = 1 - clip->cuts[index]->s;
#if 0
                printf("swapped\n");
#endif
            }

            for (j = 1; j <= i; j++) {

                if (clip->cuts[j]->v1 == vert)
                    s_j = clip->cuts[j]->s;
                else {
                    s_j = 1 - clip->cuts[j]->s;
#if 0
                    printf("swapped\n");
#endif
                }

                if (s_j > s_index) {
                    index = j;
                    s_index = s_j;
                }
            }

            /* swap maximum with position i */
            ctemp = clip->cuts[index];
            clip->cuts[index] = clip->cuts[i];
            clip->cuts[i] = ctemp;
        }

#if 0
        if (clip->cut_num > 1) {
            printf("after:  ");
            for (i = 0; i < clip->cut_num; i++)
                printf("%f ", clip->cuts[i]->s);
            printf("\n");
            printf("\n");
        }
#endif

    }

}


/******************************************************************************
See how a triangle is cut by the clipping triangles for the edges
in edges_near.

Entry:
  tri   - triangle to cut
  m1,m2 - meshes that are concerned
  scan  - scan that triangles are in
******************************************************************************/

cut_triangle(tri, m1, m2, scan)
Triangle* tri;
Mesh* m1, *m2;
Scan* scan;
{
    int i, j, k;
    Vector v1, v2;
    Edge* edge;
    Vector p1, p2;
    Vector dir;
    Vector pos;
    int result;
    Vertex** vs;
    int ct;
    int in_count, out_count;
    int inward[4];
    float s[4], t[4];
    int ind;
    float dot;
    Vector barycentric;

    /* clip each segment of the triangle against the clipping edges */

    for (i = 0; i < 3; i++) {

        /* get the two endpoints of a triangle edge */
        vcopy(tri->verts[i]->coord, v1);
        vcopy(tri->verts[(i + 1) % 3]->coord, v2);

        /* compare this segment with each clip edge */
        for (j = 0; j < edges_near_num; j++) {

            edge = edges_near[j];

#if 1
            if (edge->tri->dont_touch)
                continue;
#endif

            vcopy(edge->v1->coord, p1);
            vcopy(edge->v2->coord, p2);
            vsub(p2, p1, dir);

            /* intersect with each of the four clipping polygons of the edge */
            /* and save the information in a list of up to four items */
            ct = in_count = out_count = 0;
            for (k = 0; k < 4; k++) {

                /* see if there is an intersection */
                result = line_intersect_tri(v1, v2, edge->t[k],
                                            pos, &t[ct], &inward[ct], barycentric);

                if (result) {

                    /* see if the clipper and the clippee have opposite facing normals */

                    dot = tri->aa * edge->tri->aa +
                          tri->bb * edge->tri->bb +
                          tri->cc * edge->tri->cc;

                    if (dot < CLIP_BOUNDARY_COS) {
#ifdef DEBUG_CLIP
                        printf("opposite normals on clipper and clippee\n");
#endif
                        goto skip_cut;
                    }

                    vs = edge->t[k]->verts;
                    s[ct] = point_project_line(vs[0]->coord, vs[1]->coord,
                                               vs[2]->coord, pos);
                    if (k == 1 || k == 2)
                        s[ct] = 1 - s[ct];

                    /* count whether this was an inward or outward intersection */
                    if (inward[ct])
                        in_count++;
                    else
                        out_count++;

                    ct++;

                skip_cut:
                    ;

                }
            }

            /* Only one of the (up to four) intersections gets to represent */
            /* an actual intersection.  There is NO intersection if we get */
            /* an even number of intersections.  Otherwise, majority rules */
            /* with respect to inward or ourward. */

            if (ct > 1) {
#ifdef DEBUG_CLIP
                printf("more than one edge intersection on tri %d: %d\n", tri->index, ct);
                for (k = 0; k < ct; k++)
                    printf("s t: %f %f\n", s[k], t[k]);
#endif
            }

            if (in_count == out_count)
                continue;

            for (k = 0; k < ct; k++) {
                if (in_count > out_count && inward[k]) {
                    ind = k;
                    break;
                }
                if (in_count < out_count && inward[k] == 0) {
                    ind = k;
                    break;
                }
            }

            result = new_cut(edge, tri->verts[i], tri->verts[(i + 1) % 3],
                             s[ind], t[ind], inward[ind]);

            if (result) {
                vcopy(dir, pos);
                vscale(pos, s[ind]);
                vadd(p1, pos, pos);
                mesh_to_world(scan, pos, pos);
#if 0
                add_extra_line(pos, pos, 0x00ff00);
#endif
            }

        }
    }
}


/******************************************************************************
See which edges of triangles in pts_near cut across a particular triangle.

Entry:
  tri   - triangle to cut
  m1,m2 - meshes that are concerned
  scan  - scan that triangles are in
******************************************************************************/

old_cut_triangle(tri, m1, m2, scan)
Triangle* tri;
Mesh* m1, *m2;
Scan* scan;
{
    int i, j, k;
    Vertex* vert;
    Edge* edge;
    Vector x1, x2;
    float t1, t2;
    int result;
    int count = 0;

    /* mark each edge of the vertices as untouched */

    for (i = 0; i < pts_near_num; i++) {
        vert = pts_near[i];
        for (j = 0; j < vert->nedges; j++)
            vert->edges[j]->used = 0;
    }

    /* see which edges cut across the given triangle */

    for (i = 0; i < pts_near_num; i++) {

        vert = pts_near[i];

        for (j = 0; j < vert->nedges; j++) {

            edge = vert->edges[j];

            /* don't examine the same edge twice */
            if (edge->used)
                continue;

            /* mark this edge */
            edge->used = 1;

            /* see if edge cuts across any of the triangle's edges */
            for (j = 0; j < 3; j++) {

                /* find the closest approach between the lines of the edges */
                result = two_line_approach(tri->verts[j]->coord,
                                           tri->verts[(j + 1) % 3]->coord,
                                           edge->v1->coord, edge->v2->coord,
                                           x1, x2, &t1, &t2);
                if (result) {
                    fprintf(stderr, "cut_triangle: lines parallel\n");
                    continue;
                }

                /* see if the nearest approach between lines was within */
                /* the two line segments */
                if (t1 > 0 && t1 < 1 && t2 > 0 && t2 < 1) {
                    Vector xx1, xx2;
                    mesh_to_world(scan, x1, xx1);
                    mesh_to_world(scan, x2, xx2);
                    add_extra_line(xx1, xx2, 0x00ff00);
                    count++;
                }
            }

        }
    }

    printf("pts_near_num count: %d %d\n", pts_near_num, count);

}


/******************************************************************************
Find the closest approach of a point with a line.

Entry:
  pt    - the point
  q1,q2 - two points on the line

Exit:
  x - place of closest approach
  t - parameter saying where x is: x = q1 + t * (q2 - q1)
******************************************************************************/

point_line_approach(pt, q1, q2, x, t)
Vector pt;
Vector q1, q2;
Vector x;
float* t;
{
    Vector dir;
    float dirdot;
    float tt;

    /* find parameter of closest approach */
    vsub(q2, q1, dir);
    dirdot = vdot(dir, dir);
    tt = (vdot(pt, dir) - vdot(q1, dir)) / dirdot;

    /* find position of closest approach */
    vcopy(dir, x);
    vscale(x, tt);
    vadd(q1, x, x);

    /* set return value t */
    *t = tt;
}


/******************************************************************************
Find the closest approach of two lines, or signal that they are parallel.

Entry:
  p1,q1 - two points on the first line
  p2,q2 - points on second line

Exit:
  x1,x2 - points of closest approach on lines one and two, respectively
  t1,t2 - parameters of x1 and x2 when writing them as: x1 = p1 + t1 * (q1 - p1)
  returns 1 if lines are parallel, 0 if they are skew
******************************************************************************/

int two_line_approach(p1, q1, p2, q2, x1, x2, t1, t2)
Vector p1, q1;
Vector p2, q2;
Vector x1, x2;
float* t1, *t2;
{
    Vector v1, v2;
    Vector v12;
    Vector p21;
    Vector temp;
    float dot12;
    float tt1, tt2;
    static float epsilon = 0.0000001;

    vsub(q1, p1, v1);
    vsub(q2, p2, v2);
    vcross(v1, v2, v12);
    dot12 = vdot(v12, v12);

    /*
    printf ("p1: %f %f %f\n", p1[X], p1[Y], p1[Z]);
    printf ("q1: %f %f %f\n", q1[X], q1[Y], q1[Z]);
    printf ("p2: %f %f %f\n", p2[X], p2[Y], p2[Z]);
    printf ("q2: %f %f %f\n", q2[X], q2[Y], q2[Z]);
    printf ("v1: %f %f %f\n", v1[X], v1[Y], v1[Z]);
    printf ("v2: %f %f %f\n", v2[X], v2[Y], v2[Z]);
    printf ("v12: %g %g %g\n", v12[X], v12[Y], v12[Z]);
    printf ("dot12: %g\n", dot12);
    */

    if (dot12 == 0.0)
        return (1);

    vsub(p2, p1, p21);

    vcross(p21, v2, temp);
    tt1 = vdot(v12, temp);
    tt1 /= dot12;
    vcopy(v1, x1);
    vscale(x1, tt1);
    vadd(p1, x1, x1);

    vcross(p21, v1, temp);
    tt2 = vdot(v12, temp);
    tt2 /= dot12;
    vcopy(v2, x2);
    vscale(x2, tt2);
    vadd(p2, x2, x2);

    /*
    printf ("dot12 tt1 tt2: %g %f %f\n", dot12, tt1, tt2);
    */

    *t1 = tt1;
    *t2 = tt2;

    return (0);
}


/******************************************************************************
Find the collection of nearby vertices to a given position.  These vertices
must be on the edge of the mesh.

Entry:
  mesh     - the mesh
  not_mesh - mesh to reject matches from (old_mesh field of Vertex)
  pnt      - position (in mesh's coordinate system) to find nearest vertex to
  norm     - surface normal at pnt
  radius   - distance within which to check
  min_dot  - minimum allowed dot product between given normal and match point

Exit:
  places nearby vertices in "pts_near"
******************************************************************************/

verts_near_edges(mesh, not_mesh, pnt, norm, radius, min_dot)
Mesh* mesh;
Mesh* not_mesh;
Vector pnt;
Vector norm;
float radius;
float min_dot;
{
    int i;
    int a, b, c;
    int aa, bb, cc;
    int index;
    Hash_Table* table = mesh->table;
    Vertex* ptr;
    Vertex* min_ptr = NULL;
    float dx, dy, dz;
    float dist;
    float dot;

    /* allocate room to keep nearby points */
    if (pts_near == NULL) {
        pts_near_max = 50;
        pts_near = (Vertex**) malloc(sizeof(Vertex*) * pts_near_max);
    }

    /* use squared distance */
    radius = radius * radius;

    /* determine which cell the position lies within */
    aa = floor(table->scale * pnt[X]);
    bb = floor(table->scale * pnt[Y]);
    cc = floor(table->scale * pnt[Z]);

    /* look at nine cells, centered at cell containing location */

    for (a = aa - 1; a <= aa + 1; a++)
        for (b = bb - 1; b <= bb + 1; b++)
            for (c = cc - 1; c <= cc + 1; c++) {

                /* compute position in hash table */
                index = (a * PR1 + b * PR2 + c) % table->num_entries;
                if (index < 0)
                    index += table->num_entries;

                /* examine all points hashed to this cell */
                for (ptr = table->verts[index]; ptr != NULL; ptr = ptr->next) {

                    /* don't examine point if it's old mesh is not_mesh */
                    if (ptr->old_mesh == not_mesh)
                        continue;

                    /* go on if this vertex isn't on the mesh edge */
                    if (ptr->nedges == 0)
                        continue;

                    /* go on if this vertex has already been placed on list */
                    /* (count is a marker that says vertex is in pts_near) */
                    if (ptr->count)
                        continue;

                    /* distance (squared) to this point */
                    dx = ptr->coord[X] - pnt[X];
                    dy = ptr->coord[Y] - pnt[Y];
                    dz = ptr->coord[Z] - pnt[Z];
                    dist = dx * dx + dy * dy + dz * dz;

                    /* maybe we've found new closest point */
                    if (dist < radius) {

                        /* make sure the surface normals are roughly in the same direction */
                        dot = vdot(norm, ptr->normal);
                        if (dot < min_dot)
                            continue;

                        /* add this vertex to our list */
                        if (pts_near_num == pts_near_max) {
                            pts_near_max += 20;
                            pts_near = (Vertex**)
                                       realloc(pts_near, sizeof(Vertex*) * pts_near_max);
                        }
                        pts_near[pts_near_num] = ptr;
                        pts_near_num++;
                        ptr->count = 1;  /* this marks vertex as being in pts_near */
                    }
                }
            }
}


/******************************************************************************
Find nearby edges to a triangle, given a pre-computed list of nearby points
that are on the edge.

Entry:
  tri   - triangle to cut
  m1,m2 - meshes that are concerned
  scan  - scan that triangles are in
******************************************************************************/

edges_near_edges(tri, m1, m2, scan)
Triangle* tri;
Mesh* m1, *m2;
Scan* scan;
{
    int i, j, k;
    Vertex* vert;
    Edge* edge;
    Vector x1, x2;
    float t1, t2;
    int result;
    int count = 0;

    /* empty the list of nearby edges */
    edges_near_num = 0;

    /* allocate room to keep nearby edges */
    if (edges_near == NULL) {
        edges_near_max = 50;
        edges_near = (Edge**) malloc(sizeof(Edge*) * edges_near_max);
    }

    /* mark each edge of the vertices as untouched */

    for (i = 0; i < pts_near_num; i++) {
        vert = pts_near[i];
        for (j = 0; j < vert->nedges; j++) {
            edge = vert->edges[j];
            edge->used = 0;
            edge->prev->used = 0;
            edge->next->used = 0;
        }
    }

    /* collect together nearby edges */

    for (i = 0; i < pts_near_num; i++) {

        vert = pts_near[i];

        for (j = 0; j < vert->nedges; j++) {

            edge = vert->edges[j];

            /* place the edge in the list of nearby edges */
            /* if the edge isn't already in the list */

            if (edge->used == 0) {
                if (edges_near_num == edges_near_max) {
                    edges_near_max += 20;
                    edges_near = (Edge**)
                                 realloc(edges_near, sizeof(Edge*) * edges_near_max);
                }
                edges_near[edges_near_num] = edge;
                edges_near_num++;
                edge->used = 1;     /* mark this edge as being in edges_near */
            }

            /* see if the previous and next edges in the loop have been added */

            if (edge->prev->used == 0) {
                if (edges_near_num == edges_near_max) {
                    edges_near_max += 20;
                    edges_near = (Edge**)
                                 realloc(edges_near, sizeof(Edge*) * edges_near_max);
                }
                edges_near[edges_near_num] = edge->prev;
                edges_near_num++;
                edge->prev->used = 1;     /* mark this edge as being in edges_near */
            }

            if (edge->next->used == 0) {
                if (edges_near_num == edges_near_max) {
                    edges_near_max += 20;
                    edges_near = (Edge**)
                                 realloc(edges_near, sizeof(Edge*) * edges_near_max);
                }
                edges_near[edges_near_num] = edge->next;
                edges_near_num++;
                edge->next->used = 1;     /* mark this edge as being in edges_near */
            }
        }
    }
}


/******************************************************************************
Make a collection of triangles that border the mesh edge, to use for clipping.

Entry:
  scan   - scan to make clip triangles for
  clipto - mesh to clip to
******************************************************************************/

make_clip_triangles(scan, clipto)
Scan* scan;
Mesh* clipto;
{
    int i, j;
    Mesh* mesh;
    EdgeLoop* looplist;
    Edge* fedge, *edge, *nedge;
    int index;
    Vertex* v1, *v2;
    Vertex* v1first;
    Vertex* c1, *c2, *c3, *c4;
    Vertex* c1first, *c2first;
    Triangle* t0, *t1, *t2;
    Vector dir, pos;
    int been_around;
    Triangle* atri;
    More_Tri_Stuff* more;

    /* clear out any old mesh for the edges */
    /*** THIS SHOULD ACTUALLY FREE UP MORE STUFF !!! ***/
    /*** THIS SHOULD ACTUALLY FREE UP MORE STUFF !!! ***/
    /*** THIS SHOULD ACTUALLY FREE UP MORE STUFF !!! ***/
    if (scan->edge_mesh) {
        mesh = scan->edge_mesh;
        /* free the triangles */
        for (i = 0; i < mesh->ntris; i++) {
            if (mesh->tris[i]->more)
                free(mesh->tris[i]->more);
            free(mesh->tris[i]);
        }
        /* free the vertices */
        for (i = 0; i < mesh->nverts; i++)
            free(mesh->verts[i]);
        free(mesh->verts);
        free(mesh->tris);
        free(mesh->edges);
        scan->edge_mesh = NULL;
    }
    /*** THIS SHOULD ACTUALLY FREE UP MORE STUFF !!! ***/
    /*** THIS SHOULD ACTUALLY FREE UP MORE STUFF !!! ***/
    /*** THIS SHOULD ACTUALLY FREE UP MORE STUFF !!! ***/

    /* create new mesh for these edges */

    scan->edge_mesh = (Mesh*) malloc(sizeof(Mesh));
    mesh = scan->edge_mesh;

    /* allocate space for new triangles and vertices */

    mesh->nverts = 0;
    mesh->max_verts = 100;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    mesh->ntris = 0;
    mesh->max_tris = mesh->max_verts * 2;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 20;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 20;

    /* create triangles for the mesh */

    /* examine all edge loops in the mesh */

    looplist = &clipto->looplist;
    for (i = 0; i < looplist->nloops; i++) {

        /* look at each edge in the loop */

        fedge = looplist->loops[i];
        been_around = 0;
        for (edge = fedge; edge != fedge || !been_around; edge = nedge) {

            nedge = edge->next;

            t0 = edge->prev->tri;
            t1 = edge->tri;
            t2 = edge->next->tri;

            if (!been_around) {
                been_around = 1;

                index = make_vertex(mesh, edge->v1->coord);
                v1 = mesh->verts[index];
                v1first = v1;

                /*
                vadd (t0->normal, t1->normal, dir);
                       */
                dir[X] = -t0->aa - t1->aa;
                dir[Y] = -t0->bb - t1->bb;
                dir[Z] = -t0->cc - t1->cc;
                vnorm(dir);
                vscale(dir, CLIP_BOUNDARY_DIST);

                vadd(v1->coord, dir, pos);
                index = make_vertex(mesh, pos);
                c1 = mesh->verts[index];
                c1first = c1;

                vsub(v1->coord, dir, pos);
                index = make_vertex(mesh, pos);
                c2 = mesh->verts[index];
                c2first = c2;
            }

            if (edge->next == fedge) {
                v2 = v1first;
                c3 = c1first;
                c4 = c2first;
            } else {
                index = make_vertex(mesh, edge->v2->coord);
                v2 = mesh->verts[index];

                /*
                vadd (t1->normal, t2->normal, dir);
                       */
                dir[X] = -t1->aa - t2->aa;
                dir[Y] = -t1->bb - t2->bb;
                dir[Z] = -t1->cc - t2->cc;
                vnorm(dir);
                vscale(dir, CLIP_BOUNDARY_DIST);

                vadd(v2->coord, dir, pos);
                index = make_vertex(mesh, pos);
                c3 = mesh->verts[index];

                vsub(v2->coord, dir, pos);
                index = make_vertex(mesh, pos);
                c4 = mesh->verts[index];
            }

            /* make four triangles at the edge */

            /* What if the triangle is NULL?? - B. Curless  9/14/95 */
            edge->t[0] = make_triangle(mesh, c3, v1, c1, 1e20);
            edge->t[1] = make_triangle(mesh, v1, c3, v2, 1e20);
            edge->t[2] = make_triangle(mesh, v1, v2, c4, 1e20);
            edge->t[3] = make_triangle(mesh, c4, c2, v1, 1e20);

            /* create double-precision info for triangles */
            double_stuff(edge->t[0]);
            double_stuff(edge->t[1]);
            double_stuff(edge->t[2]);
            double_stuff(edge->t[3]);

            find_vertex_normal(v1);
            find_vertex_normal(v2);
            find_vertex_normal(c1);
            find_vertex_normal(c2);
            find_vertex_normal(c3);
            find_vertex_normal(c4);

            v1 = v2;
            c1 = c3;
            c2 = c4;
        }
    }
}


/******************************************************************************
Compute intersection between a line segment and a triangle.  This is the
single-precision version.

Entry:
  p1,p2 - endpoints of the line segment
  tri   - triangle to intersect

Exit:
  pos         - position of intersection
  tt          - parameter along line segment for where intersection occured
  inward      - whether line segment was travelling into or out of the triangle
  barycentric - barycentric coordinates of intersection (if any)
  returns 1 if they intersect, 0 if not
******************************************************************************/

int line_intersect_tri_single(p1, p2, tri, pos, tt, inward, barycentric)
Vector p1, p2;
Triangle* tri;
Vector pos;
float* tt;
int* inward;
Vector barycentric;
{
    double t;
    double pdir[3];
    double ldir[3];
    double dot1, dot2;
    double r1, r2, r3;

    /* direction of line */
    ldir[X] = p2[X] - p1[X];
    ldir[Y] = p2[Y] - p1[Y];
    ldir[Z] = p2[Z] - p1[Z];

    /* normal to plane of triangle */
    pdir[X] = tri->aa;
    pdir[Y] = tri->bb;
    pdir[Z] = tri->cc;

    /* find the intersection between the line and the plane of the tri */

    dot1 = ldir[X] * pdir[X] + ldir[Y] * pdir[Y] + ldir[Z] * pdir[Z];

    /* no intersection if line and plane are parallel */
    if (dot1 == 0.0)
        return (0);

    dot2 = p1[X] * pdir[X] + p1[Y] * pdir[Y] + p1[Z] * pdir[Z];
    t = - (dot2 + tri->dd) / dot1;

    /* no intersection if we're no longer between the segment's endpoints */
    if (t < 0 || t > 1)
        return (0);

    pos[X] = p1[X] + t * ldir[X];
    pos[Y] = p1[Y] + t * ldir[Y];
    pos[Z] = p1[Z] + t * ldir[Z];

    /* see if the intersection "pos" is on the right side of the triangle edges */

    r1 = tri->a[0] * pos[X] + tri->b[0] * pos[Y] + tri->c[0] * pos[Z] + tri->d[0];
    if (r1 < 0) return (0);

    r2 = tri->a[1] * pos[X] + tri->b[1] * pos[Y] + tri->c[1] * pos[Z] + tri->d[1];
    if (r2 < 0) return (0);

    r3 = tri->a[2] * pos[X] + tri->b[2] * pos[Y] + tri->c[2] * pos[Z] + tri->d[2];
    if (r3 < 0) return (0);

    /* if we get here then the intersection point is in the triangle */
    *tt = t;

    /* return the barycentric coordinates of the intersection */
    barycentric[X] = r2;  /* weight for tri->verts[0] */
    barycentric[Y] = r3;  /* weight for tri->verts[1] */
    barycentric[Z] = r1;  /* weight for tri->verts[2] */

    /* find which direction the line segment is travelling through the tri */

    if (dot1 < 0)
        *inward = 1;
    else
        *inward = 0;

    return (1);
}


/******************************************************************************
Compute intersection between a line segment and a triangle.  This is the
double-precision version.

Entry:
  p1,p2 - endpoints of the line segment
  tri   - triangle to intersect

Exit:
  pos         - position of intersection
  tt          - parameter along line segment for where intersection occured
  inward      - whether line segment was travelling into or out of the triangle
  barycentric - barycentric coordinates of intersection (if any)
  returns 1 if they intersect, 0 if not
******************************************************************************/

int line_intersect_tri(p1, p2, tri, pos, tt, inward, barycentric)
Vector p1, p2;
Triangle* tri;
Vector pos;
float* tt;
int* inward;
Vector barycentric;
{
    double t;
    double pdir[3];
    double ldir[3];
    double dot1, dot2;
    double r1, r2, r3;
    More_Tri_Stuff* more = tri->more;

    /* direction of line */
    ldir[X] = p2[X] - p1[X];
    ldir[Y] = p2[Y] - p1[Y];
    ldir[Z] = p2[Z] - p1[Z];

    /* normal to plane of triangle */
    pdir[X] = more->aa;
    pdir[Y] = more->bb;
    pdir[Z] = more->cc;

    /* find the intersection between the line and the plane of the tri */

    dot1 = ldir[X] * pdir[X] + ldir[Y] * pdir[Y] + ldir[Z] * pdir[Z];

    /* no intersection if line and plane are parallel */
    if (dot1 == 0.0)
        return (0);

    dot2 = p1[X] * pdir[X] + p1[Y] * pdir[Y] + p1[Z] * pdir[Z];
    t = - (dot2 + more->dd) / dot1;

    /* no intersection if we're no longer between the segment's endpoints */
    if (t < 0 || t > 1)
        return (0);

    pos[X] = p1[X] + t * ldir[X];
    pos[Y] = p1[Y] + t * ldir[Y];
    pos[Z] = p1[Z] + t * ldir[Z];

    /* see if the intersection "pos" is on the right side of the triangle edges */

    r1 = more->a[0] * pos[X] + more->b[0] * pos[Y] +
         more->c[0] * pos[Z] + more->d[0];
    if (r1 < 0) return (0);

    r2 = more->a[1] * pos[X] + more->b[1] * pos[Y] +
         more->c[1] * pos[Z] + more->d[1];
    if (r2 < 0) return (0);

    r3 = more->a[2] * pos[X] + more->b[2] * pos[Y] +
         more->c[2] * pos[Z] + more->d[2];
    if (r3 < 0) return (0);

    /* if we get here then the intersection point is in the triangle */
    *tt = t;

    /* return the barycentric coordinates of the intersection */
    barycentric[X] = r2;  /* weight for tri->verts[0] */
    barycentric[Y] = r3;  /* weight for tri->verts[1] */
    barycentric[Z] = r1;  /* weight for tri->verts[2] */

    /* find which direction the line segment is travelling through the tri */

    if (dot1 < 0)
        *inward = 1;
    else
        *inward = 0;

    return (1);
}


/******************************************************************************
Tell how far away a point in a triangle is from one of the vertices.

Entry:
  v1,v2,v3 - vertices of triangle, with v1 being the special one
  p        - point in triangle

Exit
  returns where along segment the point projects to (1 is at v1)
******************************************************************************/

float point_project_line(v1, v2, v3, p)
Vector v1, v2, v3;
Vector p;
{
    Vector d1, d2;
    Vector cross;
    float a1, a2;

    vsub(v3, v2, d1);
    vsub(p, v2, d2);
    vcross(d1, d2, cross);
    a1 = vlen(cross);

    vsub(v2, v1, d1);
    vsub(v3, v1, d2);
    vcross(d1, d2, cross);
    a2 = vlen(cross);

    return (a1 / a2);
}


/******************************************************************************
Add a new entry to the list of intersections between a triangle's line
segment and an edge.

Entry:
  edge   - edge whose list we're adding to
  v1,v2  - vertices that define line segment that intersects
  t      - parameter saying where intersection occured along edge
  s      - parameter saying where intersection occured along v1 - v2
  inward - whether line segment is passing into or out of the edge's mesh

Exit:
  returns 1 if intersection was added, 0 if the segment was added previously
******************************************************************************/

int new_cut(edge, v1, v2, t, s, inward)
Edge* edge;
Vertex* v1, *v2;
float t;
float s;
int inward;
{
    int i, j;
    int index;
    Vertex* vtemp;
    Triangle* tri1, *tri2;
    Cut* cut;

    /* see that there isn't already an intersection between this segment */
    /* and the given edge */
    for (i = 0; i < edge->cut_num; i++)
        if ((v1 == edge->cuts[i]->v1 && v2 == edge->cuts[i]->v2) ||
            (v1 == edge->cuts[i]->v2 && v2 == edge->cuts[i]->v1))
            return (0);

    /* make sure there is room for a new cut */
    if (edge->cut_num == edge->cut_max) {
        edge->cut_max += 4;
        edge->cuts = (Cut**) realloc(edge->cuts, sizeof(Cut*) * edge->cut_max);
    }

    /* add the cut to the list */
    cut = (Cut*) malloc(sizeof(Cut));
    cut->v1 = v1;
    cut->v2 = v2;
    cut->edge = edge;
    cut->t = t;
    cut->s = s;
    cut->inward = inward;
    edge->cuts[edge->cut_num++] = cut;

    /* find the one or two triangles that are bordered by the line segment */

    for (i = 0; i < v1->ntris; i++) {
        tri1 = v1->tris[i];
        for (j = 0; j < v2->ntris; j++) {
            tri2 = v2->tris[j];
            if (tri1 == tri2) {
                add_cut_to_triangle(tri1, cut);
            }
        }
    }

    /* say we've added the cut */
    return (1);
}


/******************************************************************************
Add a cut to the list of cuts for the triangle's edges.

Entry:
  tri - triangle in question
  cut - cut to add
******************************************************************************/

add_cut_to_triangle(tri, cut)
Triangle* tri;
Cut* cut;
{
    int i, j;
    int forward;
    Clipped_Edge* clips;
    Clipped_Edge* clip;
    int index;

    /* see if triangle already has a list of clipped edges */
    /* and if not, create a list of three clipped edges */

    if (tri->clips == NULL) {
        clips = (Clipped_Edge*) malloc(sizeof(Clipped_Edge) * 3);
        tri->clips = clips;
        for (i = 0; i < 3; i++) {
            clips[i].v1 = tri->verts[i];
            clips[i].v2 = tri->verts[(i + 1) % 3];
            clips[i].cut_max = 2;
            clips[i].cut_num = 0;
            clips[i].cuts = (Cut**) malloc(sizeof(Cut*) * clips[i].cut_max);
            clips[i].perp_intersect = 0;
        }
    }

    /* look for appropriate clipping edge to add the cut to */
    forward = -1;
    for (i = 0; i < 3; i++) {
        if (tri->verts[i] == cut->v1 && tri->verts[(i + 1) % 3] == cut->v2) {
            clip = &tri->clips[i];
            forward = 1;
            index = i;  /* for debug */
            break;
        }
        if (tri->verts[i] == cut->v2 && tri->verts[(i + 1) % 3] == cut->v1) {
            clip = &tri->clips[i];
            forward = 0;
            index = i;  /* for debug */
            break;
        }
    }

    /* consistancy check */
    if (forward == -1) {
        fprintf(stderr, "add_cut_to_triangle: couldn't find correct edge\n");
        exit(-1);
    }

    /* add this cut to the clipped edge's list of cuts */
    if (clip->cut_num == clip->cut_max) { /* check to see if there is room */
        clip->cut_max += 2;
        clip->cuts = (Cut**) realloc(clip->cuts, sizeof(Cut*) * clip->cut_max);
    }
    clip->cuts[clip->cut_num++] = cut;
}


/******************************************************************************
Initialize the list of intersection points for each edge on the boundary
of the clipping mesh.

Entry:
  scan   - scan containing mesh
  clipto - mesh that will be clipped to
******************************************************************************/

init_cuts(scan, clipto)
Scan* scan;
Mesh* clipto;
{
    int i, j;
    EdgeLoop* looplist;
    Edge* fedge, *edge;
    int been_around;

    /* examine all edge loops in the mesh */

    looplist = &clipto->looplist;
    for (i = 0; i < looplist->nloops; i++) {

        /* look at each edge in the loop */

        fedge = looplist->loops[i];
        been_around = 0;
        for (edge = fedge; edge != fedge || !been_around; edge = edge->next) {

            been_around = 1;

            if (edge->cuts == NULL) {
                edge->cut_max = 4;
                edge->cuts = (Cut**) malloc(sizeof(Cut*) * edge->cut_max);
            }

            edge->cut_num = 0;
        }
    }

}


/******************************************************************************
Create the vertices at intersections between triangles and the mesh edge.
Also sort the cuts in distance along each edge.

Entry:
  mesh - the mesh in which to introduce the new vertices
******************************************************************************/

create_cut_vertices(mesh)
Mesh* mesh;
{
    int i, j;
    EdgeLoop* looplist;
    Edge* fedge, *edge;
    int been_around;
    Vector c1, c2;
    Vector coord;
    float t;
    int vert_index;
    Vertex* vert;

    /* examine all edge loops in the mesh */
    looplist = &mesh->looplist;
    for (i = 0; i < looplist->nloops; i++) {

        /* look at each edge in the loop */
        fedge = looplist->loops[i];
        been_around = 0;
        for (edge = fedge; edge != fedge || !been_around; edge = edge->next) {

            been_around = 1;

            /* sort the cuts along the edge */

            if (edge->cut_num > 0)
                sort_cuts(edge);

            /* create a vertex at each cut */

            vcopy(edge->v1->coord, c1);
            vcopy(edge->v2->coord, c2);

            for (j = 0; j < edge->cut_num; j++) {
                t = edge->cuts[j]->t;
                coord[X] = c1[X] + t * (c2[X] - c1[X]);
                coord[Y] = c1[Y] + t * (c2[Y] - c1[Y]);
                coord[Z] = c1[Z] + t * (c2[Z] - c1[Z]);
                vert_index = make_vertex(mesh, coord);
                vert = mesh->verts[vert_index];
                add_to_hash(vert, mesh);
                edge->cuts[j]->new_vert = vert;
                vert->confidence = edge->v1->confidence + t * edge->v2->confidence;

                /* signal that this is a new vertex (kluge!) */
                vert->moving = 1;
            }

            /* mark all edges as not yet examined (for next phase) */
            edge->used = 0;
        }
    }

}


/******************************************************************************
Introduce the cut points (of edges) as new vertices of a mesh.

Entry:
  mesh     - the mesh to alter
  not_mesh - the other mesh that has been merged with "mesh" that we DON'T
         want to modify
******************************************************************************/

introduce_all_cuts(mesh, not_mesh)
Mesh* mesh;
Mesh* not_mesh;
{
    int i, j, k;
    EdgeLoop* looplist;
    Edge* fedge, *edge;
    Edge* edges[5];
    int edge_count;
    int been_around;
    Triangle* tri;
    Triangle* tnew1, *tnew2;
    Triangle* dummy1, *dummy2;
    Vertex* v;

    /* Examine each triangle on the boundary of the mesh and introduce */
    /* the appropriate cuts into its edges.  This loop counts backwards */
    /* through the triangles because we're creating and deleting triangles. */

    for (i = mesh->ntris - 1; i >= 0; i--) {

        tri = mesh->tris[i];

        /* don't examine triangles from the other mesh */
        if (tri->verts[0]->old_mesh == not_mesh)
            continue;

        /* see if this triangle is on the mesh boundary by examining */
        /* each vertex of the triangle */
        edge_count = 0;
        for (j = 0; j < 3; j++) {
            v = tri->verts[j];
            /* look at each edge of the vertex */
            for (k = 0; k < v->nedges; k++) {
                edge = v->edges[k];
                if (edge->cut_num > 0 && edge->used == 0 && edge->tri == tri) {
                    if (edge_count == 3) {
                        fprintf(stderr, "introduce_all_cuts: too many edges on triangle\n");
                        exit(-1);
                    }
                    edges[edge_count] = edge;
                    edge_count++;
                    edge->used = 1;
                }
            }
        }

        /*
        if (edge_count > 1)
          printf ("edge_count: %d\n", edge_count);
        */

        /* introduce cuts into the triangle */

        if (edge_count > 0)
            introduce_cuts(tri, edges[0], mesh, &tnew1, &tnew2);

        if (edge_count > 1) {
            if (edges[1]->v2 == edges[0]->v1)
                introduce_cuts(tnew1, edges[1], mesh, &dummy1, &dummy2);
            else if (edges[1]->v1 == edges[0]->v2)
                introduce_cuts(tnew2, edges[1], mesh, &dummy1, &dummy2);
            else {
                fprintf(stderr, "introduce_all_cuts: bad vertex correspondence\n");
                exit(-1);
            }
        }

        if (edge_count > 2) {
            if (edges[2]->v2 == edges[0]->v1)
                introduce_cuts(tnew1, edges[2], mesh, &dummy1, &dummy2);
            else if (edges[2]->v1 == edges[0]->v2)
                introduce_cuts(tnew2, edges[2], mesh, &dummy1, &dummy2);
            else {
                fprintf(stderr, "introduce_all_cuts: bad vertex correspondence\n");
                exit(-1);
            }
        }

    }
}


/******************************************************************************
Cut a triangle into several pieces based on the "cuts" described along
an edge.

Entry:
  tri  - triangle to cut
  edge - edge along which to cut
  mesh - the mesh containing the triangle

Exit:
  first_tri - first new triangle created
  last_tri  - last new triangle created
******************************************************************************/

introduce_cuts(tri, edge, mesh, first_tri, last_tri)
Triangle* tri;
Edge* edge;
Mesh* mesh;
Triangle** first_tri, ** last_tri;
{
    int i, j, k;
    int index;
    Vertex* v1, *v2, *v3;
    Vector coord;
    int vert_index;
    Vertex* old_vert, *new_vert;
    float t;

    /*
    for (i = 0; i < edge->cut_num; i++)
      printf ("%10f ", edge->cuts[i].t);
    printf ("\n");
    */

    /* Warning:  The code below assumes that the vertices in an */
    /* edge appear in a consistant order around a triangle.  This */
    /* should be the case when an edge was created by add_edge_to_mesh(). */

    index = -1;
    for (i = 0; i < 3; i++) {
        if (tri->verts[i] == edge->v1 && tri->verts[(i + 1) % 3] == edge->v2) {
            index = i;
            break;
        }
    }

    if (index == -1) {
        printf("introduce_cuts: vertices in wrong order along an edge\n");
        exit(-1);
    }

    v1 = tri->verts[index];
    v2 = tri->verts[(index + 1) % 3];
    v3 = tri->verts[(index + 2) % 3];

    /* create a new vertex and triangle for every cut along an edge */

    old_vert = v1;
    for (i = 0; i < edge->cut_num; i++) {

        new_vert = edge->cuts[i]->new_vert;

        /* create new triangle, saving pointer to the first one */
        if (i == 0) {
            *first_tri = make_triangle(mesh, old_vert, new_vert, v3, 1e20);
        } else {
            make_triangle(mesh, old_vert, new_vert, v3, 1e20);
        }

        old_vert = new_vert;
    }

    /* create last triangle */
    *last_tri = make_triangle(mesh, new_vert, v2, v3, 1e20);

    /* delete the original triangle */
    delete_triangle(tri, mesh, 1);

    /* re-assess whether these vertices are on an edge of the mesh */
    vertex_edge_test(v1);
    vertex_edge_test(v2);
    vertex_edge_test(v3);
    for (i = 0; i < edge->cut_num; i++)
        vertex_edge_test(edge->cuts[i]->new_vert);

    /* compute the normals at the new vertices */
    for (i = 0; i < edge->cut_num; i++)
        find_vertex_normal(edge->cuts[i]->new_vert);

    /* mark that the edges are not valid in this mesh */
    mesh->edges_valid = 0;
}


/******************************************************************************
Sort the cuts along a particular edge.

Entry:
  edge - edge whose cuts are to be sorted
******************************************************************************/

sort_cuts(edge)
Edge* edge;
{
    int i, j, k;
    int index;
    Cut* temp;

    /* sort the edge's cuts by parameter that specifies position along edge */
    /* (selection sort) */

    /*
    printf ("before: ");
    for (i = 0; i < edge->cut_num; i++)
      printf ("%10f ", edge->cuts[i]->t);
    printf ("\n");
    */

    for (i = edge->cut_num - 1; i > 0; i--) {
        /* find maximum for j in [0,i] */
        index = 0;
        for (j = 1; j <= i; j++) {
            if (edge->cuts[j]->t > edge->cuts[index]->t)
                index = j;
        }
        /* swap maximum with position i */
        temp = edge->cuts[index];
        edge->cuts[index] = edge->cuts[i];
        edge->cuts[i] = temp;
    }

    /*
    printf ("after:  ");
    for (i = 0; i < edge->cut_num; i++)
      printf ("%10f ", edge->cuts[i]->t);
    printf ("\n\n");
    */

}


/******************************************************************************
Determine equation of the plane containing three vectors.  Double-precision
version.

Entry:
  v0,v1,v2 - vectors

Exit:
  aa,bb,cc,dd - describes plane
  returns 0 if completed normally, 1 if degenerate plane
******************************************************************************/

int plane_thru_vectors_double(v0, v1, v2, aa, bb, cc, dd)
double* v0, *v1, *v2;
double* aa, *bb, *cc, *dd;
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
        fprintf(stderr,
                "plane_thru_vectors_double: degenerate point configuration\n");
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
Make the "more" field for a triangle.  This includes creating double-precision
values for a triangle's plane equation and edge planes.

Entry:
  tri - triangle to compute edge planes and plane equation for
******************************************************************************/

double_stuff(tri)
Triangle* tri;
{
    int i, j;
    double v[3];
    double v0[3], v1[3], v2[3];
    double a, b, c, d;
    More_Tri_Stuff* more;

    /* create room for new info */
    more = (More_Tri_Stuff*) malloc(sizeof(More_Tri_Stuff));
    tri->more = more;

    /* find plane equation of triangle */
    for (i = 0; i < 3; i++) {
        v0[i] = tri->verts[0]->coord[i];
        v1[i] = tri->verts[1]->coord[i];
        v2[i] = tri->verts[2]->coord[i];
    }
    plane_thru_vectors_double(v0, v1, v2,
                              &more->aa, &more->bb, &more->cc, &more->dd);

    /* find plane through two vertices and perpendicular to plane of triangle */

    for (i = 0; i < 3; i++) {

        for (j = 0; j < 3; j++) {
            v0[j] = tri->verts[i]->coord[j];
            v1[j] = tri->verts[(i + 1) % 3]->coord[j];
            v2[j] = tri->verts[(i + 2) % 3]->coord[j];
        }

        /* make a vertex out of v0 perpendicular to plane of polygon */
        v[X] = v0[X] + tri->aa;
        v[Y] = v0[Y] + tri->bb;
        v[Z] = v0[Z] + tri->cc;

        /* find plane through these three vertices */
        plane_thru_vectors_double(v, v0, v1, &a, &b, &c, &d);

        /* points in polygon are positive when plugged into plane equation */
        if (v2[X] * a + v2[Y] * b + v2[Z] * c + d > 0) {
            more->a[i] = a;
            more->b[i] = b;
            more->c[i] = c;
            more->d[i] = d;
        } else {
            more->a[i] = -a;
            more->b[i] = -b;
            more->c[i] = -c;
            more->d[i] = -d;
        }
    }
}

