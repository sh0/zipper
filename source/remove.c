/*

Remove redundant triangles.

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

int draw_during_ops = 1;        /* draw during various operations? */

extern Scan* match_from[];
extern Scan* match_to[];
extern int num_drag[];
extern int num_matches;
extern int processes_forked;
extern int parallel_procs_max;

static float global_near_dist;        /* for passing to mark_for_eating */

static float EAT_NEAR_COS;
static float EAT_NEAR_DIST_FACTOR;
static float EAT_NEAR_DIST;
static int EAT_START_ITERS;
static float EAT_START_FACTOR;
static float EAT_START_DIST;


void update_eat_resolution()
{
    EAT_NEAR_DIST = ZIPPER_RESOLUTION * EAT_NEAR_DIST_FACTOR;
}

void set_eat_near_dist_factor(float factor)
{
    EAT_NEAR_DIST_FACTOR = factor;
    EAT_NEAR_DIST = ZIPPER_RESOLUTION * EAT_NEAR_DIST_FACTOR;
}

float get_eat_near_dist_factor()
{
    return EAT_NEAR_DIST_FACTOR;
}

void set_eat_near_cos(float cosine)
{
    EAT_NEAR_COS = cosine;
}

float get_eat_near_cos()
{
    return EAT_NEAR_COS;
}

void set_eat_start_iters(int iters)
{
    EAT_START_ITERS = iters;
}

int get_eat_start_iters()
{
    return EAT_START_ITERS;
}

void set_eat_start_factor(float factor)
{
    EAT_START_FACTOR = factor;
    EAT_START_DIST = ZIPPER_RESOLUTION * EAT_START_FACTOR;
}

float get_eat_start_factor()
{
    return EAT_START_FACTOR;
}

/******************************************************************************
Zipper together everything all at once.
******************************************************************************/
do_it_all()
{
    printf("Zipper: eat_edge_pair()\n");
    eat_edge_pair(scans[0], scans[1]);
    printf("Zipper: zipper_meshes()\n");
    zipper_meshes(scans[0], scans[1]);
    printf("Zipper: fill_in_holes()\n");
    fill_in_holes(scans[0], scans[1]);
    printf("Zipper: move_vertices()\n");
    move_vertices(scans[1], scans[0]);

    if (nscans >= 4) {

        printf("starting second part\n");

        eat_edge_pair(scans[2], scans[3]);
        zipper_meshes(scans[2], scans[3]);
        fill_in_holes(scans[2], scans[3]);
        move_vertices(scans[3], scans[2]);

        printf("starting third part\n");

        eat_edge_pair(scans[0], scans[2]);
        zipper_meshes(scans[0], scans[2]);
        fill_in_holes(scans[0], scans[2]);
        move_vertices(scans[2], scans[0]);
    }

    if (nscans == 6) {

        int num;

        printf("starting third part\n");

        num = 1;
        init_eating(scans[4]);
        while (num > 0)
            num = eat_mesh_edges(scans[0], scans[4], 1, 0, EAT_NEAR_DIST);
        done_eating(scans[4]);
        zipper_meshes(scans[0], scans[4]);
        fill_in_holes(scans[0], scans[4]);
        move_vertices(scans[4], scans[0]);

        printf("starting fourth part\n");

        num = 1;
        init_eating(scans[5]);
        while (num > 0)
            num = eat_mesh_edges(scans[0], scans[5], 1, 0, EAT_NEAR_DIST);
        done_eating(scans[5]);
        zipper_meshes(scans[0], scans[5]);
        fill_in_holes(scans[0], scans[5]);
        move_vertices(scans[5], scans[0]);
    }
}


/******************************************************************************
Eat away at the edges of two meshes where they are aligned.
******************************************************************************/
void eat_edge_proc()
{
    int i;
    extern float time_it();
    float time;

    time = time_it();

    if (num_matches == 0) {
        if (scans[0] != NULL && scans[1] != NULL)
            eat_edge_pair(scans[0], scans[1]);
        else {
            printf("Cannot merge scans.\n");
        }
    } else
        for (i = 0; i < num_matches; i++)
            eat_edge_pair(match_to[i], match_from[i]);

    /*
      printf ("done eating edges\n");
    */

    time = time_it();
    printf("time: %.3f\n", time);
}

/******************************************************************************
Repeatedly eat away at a pair of meshes.
******************************************************************************/
eat_edge_pair(sc1, sc2)
Scan* sc1, *sc2;
{
    int i;
    int num1, num2;

    num1 = num2 = 1;

    /* Is this necessary???? */
#if 0
    find_mesh_edges(sc1->meshes[mesh_level]);
    find_mesh_edges(sc2->meshes[mesh_level]);
#endif

    init_eating(sc1);
    init_eating(sc2);

#if 1

    /* use some larger distances for a few cycles */

    use_large_search(1);

    for (i = 0; i < EAT_START_ITERS; i++) {
        num1 = eat_mesh_edges(sc1, sc2, 0, 2, 1, EAT_NEAR_DIST * EAT_START_FACTOR);
        num2 = eat_mesh_edges(sc2, sc1, 0, 2, 1, EAT_NEAR_DIST * EAT_START_FACTOR);
        if (num1 == 0 && num2 == 0)
            break;
    }

    use_large_search(0);

#endif

    /* first eat only those triangles with a majority of less certain */
    /* vertices than their nearby positions on the other surface */

    while (num1 || num2) {
        if (num1)
            num1 = eat_mesh_edges(sc1, sc2, 0, 2, 1, EAT_NEAR_DIST);
        if (num2)
            num2 = eat_mesh_edges(sc2, sc1, 0, 2, 1, EAT_NEAR_DIST);
    }

    num1 = num2 = 1;

    /* now eat ANY triangles on top of the other surface */

    while (num1 || num2) {
        if (num1)
            num1 = eat_mesh_edges(sc1, sc2, 0, 0, 1, EAT_NEAR_DIST);
        if (num2)
            num2 = eat_mesh_edges(sc2, sc1, 0, 0, 1, EAT_NEAR_DIST);
    }

    if (processes_forked) {
        processes_forked = 0;
    }

    done_eating(sc1);
    done_eating(sc2);
}


/******************************************************************************
Make ready for eating away at a mesh.  Must be called before calling
eat_mesh_edges().

Entry:
  scan - the scan whose mesh will be eaten
******************************************************************************/
init_eating(scan)
Scan* scan;
{
    int i;
    Mesh* mesh;
    Triangle* tri;
    int r1, r2, r3;

    mesh = scan->meshes[mesh_level];

    /* categories of how a triangle lies relative to sitting on the other mesh */

#define UNTOUCHED 0  /* un-categorized triangle */
#define ON_LIST   1  /* edge triangle on list to be examined */
#define OFF_MESH  2  /* triangle that is definitely NOT on another mesh */
#define REMOVE    3  /* triangle to be removed */

    /* mark all triangles of the mesh as un-categorized */

    for (i = 0; i < mesh->ntris; i++)
        mesh->tris[i]->eat_mark = UNTOUCHED;

    /* mark all vertices of the mesh as un-categorized */

    for (i = 0; i < mesh->nverts; i++)
        mesh->verts[i]->count = UNTOUCHED;

    /* see which triangles are currently on the mesh edge, and place */
    /* these on the list to examine */

    mesh->eat_list_num = 0;
    mesh->eat_list = (Triangle**)
                     malloc(sizeof(Triangle*) * mesh->eat_list_max);

    for (i = 0; i < mesh->ntris; i++) {

        tri = mesh->tris[i];

#if 1
        if (tri->dont_touch)
            continue;
#endif

        /* only look at triangles that are on the mesh edge */
        r1 = tri->verts[0]->on_edge;
        r2 = tri->verts[1]->on_edge;
        r3 = tri->verts[2]->on_edge;
        if (r1 + r2 + r3 == 0)
            continue;

        /* place this triangle on the list */
        tri->eat_mark = ON_LIST;
        if (mesh->eat_list_num == mesh->eat_list_max) {
            mesh->eat_list_max *= 2;
            mesh->eat_list = (Triangle**)
                             realloc(mesh->eat_list, sizeof(Triangle*) * mesh->eat_list_max);
        }
        mesh->eat_list[mesh->eat_list_num++] = tri;
    }
}


/******************************************************************************
Clean up after eating away at a mesh.  Should be called after calling
eat_mesh_edges().

Entry:
  scan - the scan whose mesh was be eaten
******************************************************************************/
done_eating(scan)
Scan* scan;
{
    int i;
    Mesh* mesh;

    mesh = scan->meshes[mesh_level];

    /* remove marks from all triangles and vertices */

    for (i = 0; i < mesh->ntris; i++)
        mesh->tris[i]->eat_mark = 0;

    for (i = 0; i < mesh->nverts; i++)
        mesh->verts[i]->count = 0;

    free(mesh->eat_list);
}

/******************************************************************************
Remove the triangles from the edge of one mesh that match the surface of
another mesh.  Must have called init_eating() beforehand.

Entry:
  sc1  - scan with mesh to match
  sc2  - scan with mesh to remove triangles from
  draw - whether or not to draw after eating
  conf - how many of the nearby positions on the other surface must be
     more confident than the triangle's vertices for the tri to be deleted
  to_edge - eat all the way to the edge?  If not, then stop short.
  near_dist - distance in which to search

Exit:
  returns how many triangles were removed
******************************************************************************/
int eat_mesh_edges(sc1, sc2, draw, conf, to_edge, near_dist)
Scan* sc1, *sc2;
int draw;
int conf, to_edge;
float near_dist;
{
    int i, j, k;
    Mesh* m2;
    Triangle* tri;
    Triangle* otri;
    Vertex* vert;
    int count;
    void mark_for_eating();

    m2 = sc2->meshes[mesh_level];

    /* Look through triangles on the list of potentially deletable ones, */
    /* marking those that should be deleted.  Maybe do this in parallel */
    /* on several processors. */

    global_near_dist = near_dist;

#if 0
    m_set_procs(parallel_procs_max);
    m_fork(mark_for_eating, sc1, sc2, draw, conf, to_edge);
    processes_forked = 1;
#else
    mark_for_eating(sc1, sc2, draw, conf, to_edge);
#endif

    /* update the list of triangles to be examined */
    for (i = m2->eat_list_num - 1; i >= 0; i--) {

        tri = m2->eat_list[i];

        /* if we're going to remove this triangle, examine its neighbors */
        /* to see if they need to be placed on the list */
        if (tri->eat_mark == REMOVE) {
            /* examine each vertex of triangle */
            for (j = 0; j < 3; j++) {
                vert = tri->verts[j];
                /* look at each triangle in these vertices */
                for (k = 0; k < vert->ntris; k++) {
                    otri = vert->tris[k];
#if 1
                    if (otri->dont_touch)
                        continue;
#endif
                    /* place this triangle on the list if it is untouched */
                    if (otri->eat_mark == UNTOUCHED) {
                        otri->eat_mark = ON_LIST;
                        if (m2->eat_list_num == m2->eat_list_max) {
                            m2->eat_list_max *= 2;
                            m2->eat_list = (Triangle**)
                                           realloc(m2->eat_list, sizeof(Triangle*)*m2->eat_list_max);
                        }
                        m2->eat_list[m2->eat_list_num++] = otri;
                    }
                }
            }
        }
    }

    /* delete the marked triangles */
    count = 0;
    for (i = m2->eat_list_num - 1; i >= 0; i--) {
        tri = m2->eat_list[i];
        if (tri->eat_mark == REMOVE) {
            delete_triangle(tri, m2, 1);
            count++;
            m2->eat_list[i] = m2->eat_list[--m2->eat_list_num];
        } else if (tri->eat_mark == OFF_MESH) {
            m2->eat_list[i] = m2->eat_list[--m2->eat_list_num];
        }
    }

    /* mark the edges of this mesh as invalid */
    m2->edges_valid = 0;

    return (count);
}

/******************************************************************************
Go through list of potentially deletable triangles, marking them if they
should be deleted.

Entry:
  sc1   - scan with mesh to match
  sc2   - scan with mesh to remove triangles from
  draw  - whether or not to draw after eating
  conf  - how many of the nearby positions on the other surface must be
          more confident than the triangle's vertices for the tri to be deleted
  to_edge - mark all the way to the edge?  If not, then stop short.

******************************************************************************/
void mark_for_eating(sc1, sc2, draw, conf, to_edge)
Scan* sc1, *sc2;
int draw;
int conf, to_edge;
{
    int i, j;
    Mesh* m1, *m2;
    Triangle* tri;
    Vertex* vert;
    Vector c1;
    Vector norm1;
    int r1;
    NearPosition n1;
    int inc;
    float confidences[3];
    int conf_count;
    float near_dist;
    int myid;
    int start, end;
    int processes;

    near_dist = global_near_dist;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    inc = level_to_inc(mesh_level);

    /* figure out which sub-range to examine, given my processor id */

    myid = 0;
    processes = 1;
    start = m2->eat_list_num * (myid / (float) processes);
    end   = m2->eat_list_num * ((myid + 1) / (float) processes) - 1;

#if 0
    printf("start end num: %d %d %d\n", start, end, m2->eat_list_num);
#endif

    /* examine all triangles in the list to examine */

    for (i = start; i <= end; i++) {

        tri = m2->eat_list[i];

        /* examine each triangle vertex in turn to see if they lie */
        /* on the other mesh */

        for (j = 0; j < 3; j++) {

            vert = tri->verts[j];

            /* if we already know that this vertex is off the mesh, */
            /* mark the triangle as such and move on */
            if (vert->count == OFF_MESH) {
                tri->eat_mark = OFF_MESH;
                goto break_loop;
            }

            /* see if this vertex is on the other mesh */
            mesh_to_world(sc2, vert->coord, c1);
            mesh_to_world_normal(sc2, vert->normal, norm1);
            r1 = nearest_on_mesh(sc1, m1, NULL, c1, norm1, inc * near_dist,
                                 EAT_NEAR_COS, &n1);
            r1 = r1 && (n1.on_edge == 0);

            if (r1 && !to_edge) {
                r1 = r1 && !is_near_edge(&n1);
                if (!r1) {
#if 0
                    fprintf(stderr, "%d\n", n1.type);
                    if (n1.type == NEAR_TRIANGLE)
                        delete_triangle(n1.tri, sc1->meshes[mesh_level], 0);
#endif
                    goto break_loop;
                }
            }

            confidences[j] = n1.confidence;

            /* if the vertex is NOT on the mesh, mark the triangle and go on */
            if (!r1) {
                tri->eat_mark = OFF_MESH;
                vert->count = OFF_MESH;
                goto break_loop;
            }
        }

        /* remove this triangle if it's confidence is greater than */
        /* the near points to it */
        conf_count  = confidences[0] > tri->verts[0]->confidence;
        conf_count += confidences[1] > tri->verts[1]->confidence;
        conf_count += confidences[2] > tri->verts[2]->confidence;
        if (conf_count >= conf) {
            tri->eat_mark = REMOVE;
        }

#if 0
        /* if we get here then all three vertices are on the mesh, so */
        /* we should mark the triangle for removal */
        tri->eat_mark = REMOVE;
#endif

    break_loop: ;  /* to break out of nested loop */
    }
}


/******************************************************************************
Zipper together a model.
******************************************************************************/
void zipper_proc()
{
    move_vertices(scans[1], scans[0]);
}


/******************************************************************************
Align edges of meshes.
******************************************************************************/
void align_proc()
{
    int i;

    if (num_matches == 0) {
        zipper_meshes(scans[0], scans[1]);
        fill_in_holes(scans[0], scans[1]);
    } else
        for (i = 0; i < num_matches; i++) {
            zipper_meshes(match_from[i], match_to[i]);
            fill_in_holes(match_from[i], match_to[i]);
        }
}

/******************************************************************************
Gather the triangles of two meshes together into one mesh.

Entry:
  sc1 - scan with first mesh. all triangles go into this mesh
  sc2 - scan with second mesh. to be emptied
******************************************************************************/
gather_triangles(sc1, sc2)
Scan* sc1, *sc2;
{
    int i;
    Mesh* m1, *m2;
    Vertex* vert;
    Triangle* tri;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    /*** move the vertices from mesh 2 to mesh 1 ***/

    /* create room for new vertices */
    m1->max_verts += m2->nverts;
    m1->verts =
        (Vertex**) realloc(m1->verts, sizeof(Vertex*) * m1->max_verts);

    /* check for error */
    if (m1->verts == NULL) {
        fprintf(stderr, "could not realloc vertices\n");
        exit(-1);
    }

    /* mark the original vertices of mesh 1 as being from mesh 1 */
    for (i = 0; i < m1->nverts; i++)
        m1->verts[i]->old_mesh = m1;

    /* actually move the vertices */
    for (i = 0; i < m2->nverts; i++) {

        vert = m2->verts[i];

        /* convert between coordinate systems */
        mesh_to_world(sc2, vert->coord, vert->coord);
        world_to_mesh(sc1, vert->coord, vert->coord);
        add_to_hash(vert, m1);

        /* place vertex in mesh 1 */
        m1->verts[m1->nverts] = vert;
        vert->index = m1->nverts;
        m1->nverts++;

        /* mark the moved vertices as coming from m2 */
        vert->old_mesh = m2;
    }

    /* free up space in mesh 2 */
    free(m2->verts);
    m2->verts = NULL;

    /*** move the triangles ***/

    /* create room for new triangles */
    m1->max_tris += m2->ntris;
    m1->tris =
        (Triangle**) realloc(m1->tris, sizeof(Triangle*) * m1->max_tris);

    /* check for error */
    if (m1->tris == NULL) {
        fprintf(stderr, "could not realloc triangles\n");
        exit(-1);
    }

    /* move the triangles */
    for (i = 0; i < m2->ntris; i++) {

        tri = m2->tris[i];

        /* add triangle to destination list */
        m1->tris[m1->ntris] = tri;
        tri->index = m1->ntris;
        m1->ntris++;

        /* compute normal vector, plane coefficients and edge planes for triangle */
        set_triangle_geometry(tri);
    }

    /* re-calculate normals at vertices */
    find_vertex_normals(m1);

    /* mark any edge information as invalid */
    m1->edges_valid = 0;
    m2->edges_valid = 0;

    /* free up hash table space in mesh 2 */
    free(m2->table->verts);
    m2->table->verts = NULL;

    /* free up triangle space in mesh 2 */
    free(m2->tris);
    m2->tris = NULL;

    /* set mesh 2 to be empty */
    m2->ntris = 0;
    m2->nverts = 0;
}

/******************************************************************************
Zipper together nearby polygons on different meshes.

Entry:
  sc1 - scan with first mesh
  sc2 - scan with second mesh
******************************************************************************/
zipper_meshes(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    NearPosition near_info;
    int result;
    Vertex* vert;
    Vector pos, norm;
    Vector near_pos;
    Vector diff;
    static Vertex* near_list[100];
    static float near_dist[100];
    int nlist_count = 0;
    int index;
    float min_dist;
    int inc;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    inc = level_to_inc(mesh_level);

    /* search the vertices along the edge of mesh 2 for nearby places */
    /* on mesh 1 */

    for (i = 0; i < m2->nverts; i++) {

        vert = m2->verts[i];

        /* only examine vertices on the edge of the mesh */
        if (vert->ntris == 0 || !vert->on_edge)
            continue;

        /* convert to global coordinates */
        mesh_to_world(sc2, vert->coord, pos);
        mesh_to_world_normal(sc2, vert->normal, norm);

        /* find nearest places on the other mesh */
        result = nearest_on_mesh(sc1, m1, NULL, pos, norm,
                                 inc * EAT_NEAR_DIST, FIND_COS, &near_info);

        /* go on to next vertex if there is no nearby place on the other */
        /* mesh or if the nearest position is on the edge of its mesh */
        if (result == 0 || near_info.on_edge)
            continue;

#if 0
        /* the near vertices returned are candidates for where "vert" */
        /* should be slid over to */
        switch (near_info.type) {
            case NEAR_VERTEX:
                near_list[0] = near_info.v1;
                nlist_count = 1;
                break;
            case NEAR_EDGE:
                near_list[0] = near_info.v1;
                near_list[1] = near_info.v2;
                nlist_count = 2;
                break;
            case NEAR_TRIANGLE:
                nlist_count = 3;
                break;
            default:
                fprintf(stderr, "zipper_meshes: bad type %d\n", near_info.type);
                exit(-1);
        }
#endif

        nlist_count = 0;
        if (near_info.v1->on_edge) {
            near_list[0] = near_info.v1;
            mesh_to_world(sc1, near_info.v1->coord, near_pos);
            vsub(near_pos, pos, diff);
            near_dist[0] = vlen(diff);
            nlist_count = 1;
        }

        /* make a list of nearby vertices on mesh 1 that "vert" might */
        /* be moved to */
        for (j = 0; j < near_info.v1->nverts; j++) {
            if (near_info.v1->verts[j]->on_edge) {
                near_list[nlist_count] = near_info.v1->verts[j];
                mesh_to_world(sc1, near_list[nlist_count]->coord, near_pos);
                vsub(near_pos, pos, diff);
                near_dist[nlist_count] = vlen(diff);
                nlist_count++;
            }
        }

        /* check to see if there are ANY close edges */
        if (nlist_count == 0) {
            printf("no close edge\n");
            continue;
        }

        /* find closest vertex in the list */
        min_dist = 1e20;
        for (j = 0; j < nlist_count; j++) {
            if (near_dist[j] < min_dist) {
                min_dist = near_dist[j];
                index = j;
            }
        }

        /* for now, just move this vertex to the first edge vertex */
        mesh_to_world(sc1, near_list[index]->coord, vert->coord);
        world_to_mesh(sc2, vert->coord, vert->coord);

        /* and mark that it is being moved */
        vert->moving = 1;
        vert->move_to = near_list[index];
    }
}

Triangle* new_tris[20000];
int new_tri_count;

/******************************************************************************
Fill in triangular holes between meshes left by zippering process.

Entry:
  sc1 - scan with first mesh
  sc2 - scan with second mesh
******************************************************************************/
fill_in_holes(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j, k;
    Mesh* m1, *m2;
    EdgeLoop* list2;
    Edge* e, *e_next, *e_orig;
    Vertex* vert1, *vert2;
    Triangle* tri;
    int found;
    int result;
    int index;
    int r1;
    Vector pos;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    if (!m1->edges_valid)
        create_edge_list(m1);

    if (!m2->edges_valid)
        create_edge_list(m2);

    list2 = &m2->looplist;

    /* examine each vertex of each edge loop to see if there are any */
    /* triangular holes that can be filled */

    new_tri_count = 0;

    for (i = 0; i < list2->nloops; i++) {

        /* look through each vertex in the loop */

        e_orig = list2->loops[i];
        for (e = e_orig; e->next != e_orig; e = e->next) {

            e_next = e->next;

            /* we want an un-moved vertex between the two edges */

            if (e->v2->moving == 1)
                continue;

            /* the un-shared vertices of these two edges must be moved */
            if (e->v1->moving == 0 || e_next->v2->moving == 0)
                continue;

            /* error check */
            if (e->v2 != e_next->v1) {
                fprintf(stderr, "fill_in_holes: vertices don't match\n");
                //printf ("e: %d %d\n", e->v1, e->v2);
                //printf ("e_next: %d %d\n", e_next->v1, e_next->v2);
            }

            /* find out if the two end vertices connect to an edge on */
            /* the other mesh.  if so, we can put a triangle here */

            vert1 = e->v1->move_to;
            vert2 = e_next->v2->move_to;

            /* these vertices must be different */

            if (vert1 == vert2)
                continue;

            /* both these vertices must be on an edge of the other mesh */

            if (vert1->nedges == 0 || vert2->nedges == 0)
                continue;

            /* see if they share an edge */

            found = 0;
            for (j = 0; j < vert1->nedges; j++)
                for (k = 0; k < vert2->nedges; k++)
                    if (vert1->edges[j] == vert2->edges[k]) {
                        found = 1;
                    }

            /* go on if we didn't find a common edge */

            if (!found)
                continue;

            /* create a new triangle */

            printf("new triangle\n");

            tri = (Triangle*) malloc(sizeof(Triangle));
            new_tris[new_tri_count++] = tri;

            result = find_edge_orientation(e->v1, e->v2);

            /*
                  mesh_to_world (sc1, vert1->coord, pos);
                  world_to_mesh (sc2, pos, pos);
                  r1 = make_vertex (m2, pos);
                  m2->verts[r1]->move_to = vert1;
                  m2->verts[r1]->moving = 1;

                  mesh_to_world (sc1, vert2->coord, pos);
                  world_to_mesh (sc2, pos, pos);
                  r2 = make_vertex (m2, pos);
                  m2->verts[r2]->move_to = vert2;
                  m2->verts[r2]->moving = 1;
            */

            if (result == 1) {
                tri->verts[0] = e->v2;
                tri->verts[1] = e->v1;
                tri->verts[2] = e_next->v2;
            } else {
                tri->verts[0] = e->v2;
                tri->verts[1] = e_next->v2;
                tri->verts[2] = e->v1;
            }
        }
    }

    /* examine loops for other kind of triangle hole */

    for (i = 0; i < list2->nloops; i++) {

        /* look through each vertex in the loop */

        e_orig = list2->loops[i];
        for (e = e_orig; e->next != e_orig; e = e->next) {

            /* we want both vertices of this edge to be moved */

            if (e->v1->moving == 0 || e->v2->moving == 0)
                continue;

            /* examine those vertices that this edge has been moved to */

            vert1 = e->v1->move_to;
            vert2 = e->v2->move_to;

            /* these vertices must be distinct */

            if (vert1 == vert2)
                continue;

            /* we want both these vertices to be on an edge */

            if (vert1->nedges == 0 || vert2->nedges == 0)
                continue;

            /* see if any of these two vertices edges meet */

            /* look through all adjacent vertices to vert1 */
            found = 0;
            for (j = 0; j < vert1->nverts; j++) {
                if (!vert1->verts[j]->on_edge)
                    continue;
                /* look through all adjacent vertices to vert1 */
                for (k = 0; k < vert2->nverts; k++) {
                    if (!vert2->verts[k]->on_edge)
                        continue;
                    /* check for match */
                    if (vert1->verts[j] == vert2->verts[k] &&
                        vert1->verts[j]->nedges > 0) {
                        found = 1;
                        index = j;
                        goto here;  /* how best to get out of a doubly nested loop */
                    }
                }
            }

            /* go on if we didn't get a vertex match above */

            if (!found)
                continue;

        here: /* oh no! */

            /* create a new triangle */

            printf("other kind of new triangle\n");

            tri = (Triangle*) malloc(sizeof(Triangle));
            new_tris[new_tri_count++] = tri;

            result = find_edge_orientation(e->v1, e->v2);

            mesh_to_world(sc1, vert1->verts[index]->coord, pos);
            world_to_mesh(sc2, pos, pos);
            r1 = make_vertex(m2, pos);
            m2->verts[r1]->moving = 1;
            m2->verts[r1]->move_to = vert1->verts[index];

            if (result == 1) {
                tri->verts[0] = m2->verts[r1];
                tri->verts[1] = e->v1;
                tri->verts[2] = e->v2;
            } else {
                tri->verts[0] = m2->verts[r1];
                tri->verts[1] = e->v2;
                tri->verts[2] = e->v1;
            }

        }
    }

    /* add in the new triangles */

    for (i = 0; i < new_tri_count; i++) {
        make_triangle(m2, new_tris[i]->verts[0],
                      new_tris[i]->verts[1],
                      new_tris[i]->verts[2],
                      100.0);

#if 0
        tri = m2->tris[m2->ntris - 1];
        vadd(tri->verts[0]->coord, tri->verts[1]->coord, vec);
        vadd(tri->verts[2]->coord, vec, vec);
        vnorm(vec);
        dot = vdot(vec, tri->normal);
        if (dot > 1) {
            tri->normal[X] *= -1;
            tri->normal[Y] *= -1;
            tri->normal[Z] *= -1;
        }
#endif

    }
}


/******************************************************************************
Find out which way an edge should be oriented, based on a triangle that
already contains the two vertices in the edge.

Entry:
  v1,v2 - vertices in edge

Exit:
  returns 0 or 1 based on orientation
******************************************************************************/

int find_edge_orientation(v1, v2)
Vertex* v1, *v2;
{
    int i;
    Triangle* tri;
    int index;

    for (i = 0; i < v1->ntris; i++) {

        tri = v1->tris[i];

        if (tri->verts[0] == v2)
            index = 0;
        else if (tri->verts[1] == v2)
            index = 1;
        else if (tri->verts[2] == v2)
            index = 2;
        else
            continue;

        if (tri->verts[(index + 1) % 3] == v1)
            return (0);
        else if (tri->verts[(index + 2) % 3] == v1)
            return (1);
        else {
            printf("find_edge_orientation: bad triangle info\n");
            exit(-1);
        }
    }

    printf("find_edge_orientation: can't find such a triangle\n");
    exit(-1);
}


/******************************************************************************
Move vertices from one mesh to another.

Entry:
  source - scan to move from
  dest   - scan to move to
******************************************************************************/

move_vertices(source, dest)
Scan* source, *dest;
{
    int i, j, k;
    Mesh* msource, *mdest;
    Vertex* vert;
    Vertex* dvert;
    Triangle* tri;
    int found;

    msource = source->meshes[mesh_level];
    mdest   = dest->meshes[mesh_level];

    /* move the vertices from the source mesh to the destination mesh */

    for (i = 0; i < msource->nverts; i++) {

        vert = msource->verts[i];

        /* only examine those vertices being moved */
        if (!vert->moving)
            continue;

        /* destination vertex on other mesh */
        dvert = vert->move_to;

        /* copy the triangles from one vertex to the other */
        for (j = 0; j < vert->ntris; j++) {
            /* make sure there is room in triangle list */
            if (dvert->ntris >= dvert->max_tris) {
                dvert->max_tris += 4;
                dvert->tris = (Triangle**)
                              realloc(dvert->tris, sizeof(Triangle*) * dvert->max_tris);
            }
            /* add to the list */
            dvert->tris[dvert->ntris++] = vert->tris[j];
        }

        /* copy the vertex list from one vertex to another */

        for (j = 0; j < vert->nverts; j++) {

            /* see if this vertex is already in the list */
            found = 0;
            for (k = 0; k < dvert->nverts; k++)
                if (dvert->verts[k] == vert->verts[j]->move_to) {
                    found = 1;
                    break;
                }

            /* if so, then skip this vertex */
            if (found)
                continue;

            /* make sure there is room in vertex list */
            if (dvert->nverts >= dvert->max_verts) {
                dvert->max_verts += 4;
                dvert->verts = (Vertex**)
                               realloc(dvert->verts, sizeof(Vertex*) * dvert->max_verts);
            }

            /* add to the vertex list, using the name of the vertex AFTER */
            /* the move */
            if (vert->verts[j]->moving)
                dvert->verts[dvert->nverts++] = vert->verts[j]->move_to;
            else
                dvert->verts[dvert->nverts++] = vert->verts[j];
        }
    }

    /* have the triangles refer to the moved vertices */

    for (i = 0; i < msource->ntris; i++) {
        tri = msource->tris[i];
        for (j = 0; j < 3; j++) {
            if (tri->verts[j]->moving)
                tri->verts[j] = tri->verts[j]->move_to;
        }
    }

    /* copy over the rest of the vertices */

    for (i = 0; i < msource->nverts; i++) {

        vert = msource->verts[i];

        /* only examine those vertices that were NOT moved earlier */
        if (vert->moving)
            continue;

        /* make sure there is room for new vertex */
        if (mdest->nverts == mdest->max_verts) {
            mdest->max_verts = (int)(mdest->max_verts * 1.5);
            mdest->verts = (Vertex**)
                           realloc(mdest->verts, sizeof(Vertex*) * mdest->max_verts);
        }

        /* add vertex to the list */
        mdest->verts[mdest->nverts] = vert;
        mdest->verts[mdest->nverts]->index = mdest->nverts;
        mdest->nverts++;

        /* convert between coordinate systems */
        mesh_to_world(source, vert->coord, vert->coord);
        world_to_mesh(dest, vert->coord, vert->coord);

        /* add this vertex to hash table */
        add_to_hash(vert, mdest);
    }

    /* delete triangles that have two or more vertices the same */

    for (i = msource->ntris - 1; i >= 0; i--) {

        Vertex* v1, *v2, *v3;
        tri = msource->tris[i];

        v1 = tri->verts[0];
        v2 = tri->verts[1];
        v3 = tri->verts[2];

        if ((v1 == v2) || (v1 == v3) || (v2 == v3))
            delete_triangle(tri, msource, 1);
    }

    /* copy over the triangles */

    for (i = 0; i < msource->ntris; i++) {

        tri = msource->tris[i];

        /* maybe make room for more triangles */
        if (mdest->ntris >= mdest->max_tris) {
            mdest->max_tris = (int)(mdest->max_tris * 1.5);
            mdest->tris = (Triangle**)
                          realloc(mdest->tris, sizeof(Triangle*) * mdest->max_tris);
        }

        /* add triangle to destination list */
        mdest->tris[mdest->ntris] = tri;
        mdest->tris[mdest->ntris]->index = mdest->ntris;
        mdest->ntris++;

        /* compute normal vector, plane coefficients and edge planes for triangle */
        set_triangle_geometry(tri);
    }

#if 0
    /* add in the new triangles */

    for (i = 0; i < new_tri_count; i++)
        make_triangle(mdest, new_tris[i]->verts[0],
                      new_tris[i]->verts[1],
                      new_tris[i]->verts[2],
                      100.0);
#endif

    /* re-compute normal vectors at vertices */
    find_vertex_normals(mdest);

    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/
    /*** clean up junk in source mesh ***/

    msource->nverts = 0;
    msource->ntris = 0;

    /* mark both mesh boundary structures as invalid */
    msource->edges_valid = 0;
    mdest->edges_valid   = 0;
}

