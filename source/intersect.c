/*

Intersect two meshes by clipping them against one another.

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

/* set of points near the edge of a mesh */
static Vertex** pts_near = NULL;
static int pts_near_num;
static int pts_near_max;

static int debug_cut_count = 0;

#define MESH_A    1
#define MESH_B    2
#define CUT       3
#define USED_CUT  4


/******************************************************************************
Intersect one mesh with another.

Entry:
  sc1 - first mesh
  sc2 - second mesh
******************************************************************************/

intersect_meshes(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    Vertex* vert;
    Triangle* tri;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    init_extra_lines();

#if 1
    /* set all vertex colors to neutral */
    for (i = 0; i < m1->nverts; i++)
        m1->verts[i]->confidence = -1;

    for (i = 0; i < m2->nverts; i++)
        m2->verts[i]->confidence = -1;
#endif

    debug_cut_count = 0;

    /* mark those tris in mesh 2 that are intersected by tris in mesh 1 */
    mark_intersected_tris(sc1, sc2);

    /* mark those tris in mesh 1 that are intersected by tris in mesh 2 */
    mark_intersected_tris(sc2, sc1);

    printf("debug_cut_count = %d\n", debug_cut_count);

    /* have triangle colors reflect the dont_touch flag */

    for (i = 0; i < m1->ntris; i++) {
        tri = m1->tris[i];
        tri->mark = tri->dont_touch;
#if 1
        tri->mark = 0;
#endif
        if (tri->dont_touch == 0 && tri->clips) {
            free(tri->clips);
            tri->clips = NULL;
        }
    }

    for (i = 0; i < m2->ntris; i++) {
        tri = m2->tris[i];
        tri->mark = tri->dont_touch;
#if 1
        tri->mark = 0;
#endif
        if (tri->dont_touch == 0 && tri->clips) {
            free(tri->clips);
            tri->clips = NULL;
        }
    }
}


/******************************************************************************
Finish the mesh intersection.

Entry:
  sc1 - first mesh
  sc2 - second mesh
******************************************************************************/

finish_intersect_meshes(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    Vertex* vert;
    Triangle* tri;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];


    /* comment this out if we're also doing a "merge" */
    /* comment this out if we're also doing a "merge" */
    /* comment this out if we're also doing a "merge" */

#if 1

    /* gather the two sets of traingles into the same mesh */
    gather_triangles(sc1, sc2);

#endif

    /* comment this out if we're also doing a "merge" */
    /* comment this out if we're also doing a "merge" */
    /* comment this out if we're also doing a "merge" */


    /* add the points of intersection to the mesh as new vertices */
    add_intersect_points(sc1, sc2);

    /* actually do the clipping of the triangles */
    perform_intersect_clipping(sc1, sc2);
}


/******************************************************************************
Mark which triangles of one mesh intersect with triangles of another.  We
will later do the actual clipping of the triangles.

Entry:
  sc1,sc2 - scans containing the meshes
******************************************************************************/

mark_intersected_tris(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    Vertex* vert;
    Triangle* tri;
    float max_length;
    float edge_length_max();
    Vector coord, normal;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    max_length = edge_length_max(mesh_level);

    /* mark all vertices in mesh 2 as untouched */
    /* (count will be marker saying that a vertex is in pts_near) */
    for (i = 0; i < m2->nverts; i++)
        m2->verts[i]->count = 0;

    /* see which triangles of mesh 1 are near triangles of mesh 2 */

    for (i = 0; i < m1->ntris; i++) {

        tri = m1->tris[i];

        /* find nearby vertices to this triangle */
        pts_near_num = 0;
        for (j = 0; j < 3; j++) {

            /* transform between coordinate systems */
            mesh_to_world(sc1, tri->verts[j]->coord, coord);
            world_to_mesh(sc2, coord, coord);
            mesh_to_world_normal(sc1, tri->verts[j]->normal, normal);
            world_to_mesh_normal(sc2, normal, normal);

            /* look nearby for vertices on other meshes */
            verts_near_vert(m2, NULL, coord, normal, max_length);
        }

        /* mark the current triangle if it had nearby vertices */
        if (pts_near_num) {
            tri->mark = 1;
        }

        /* un-mark all the vertices in pts_near */
        for (j = 0; j < pts_near_num; j++)
            pts_near[j]->count = 0;

        /* intersect triangle edges with nearby triangles */
        for (j = 0; j < 3; j++) {
            intersect_edge_with_near_tris(tri->verts[j], tri->verts[(j + 1) % 3], tri,
                                          sc1, sc2);
        }

#if 0
        if (tri->clips) {
            int n;
            n = tri->clips[0].cut_num + tri->clips[1].cut_num + tri->clips[2].cut_num;
            if (n > 0)
                printf("number of cuts: %d\n", n);
        }
#endif

    }

}


/******************************************************************************
Intersect an edge with all the triangles used by the nearby points collected
in the list pts_near.

Entry:
  v1,v2   - endpoints of the edge in question
  cut_tri - the triangle from which v1 and v2 come
  sc1     - mesh that v1 and v2 are from
  sc2     - mesh that the points in pts_near come from
******************************************************************************/

intersect_edge_with_near_tris(v1, v2, cut_tri, sc1, sc2)
Vertex* v1, *v2;
Triangle* cut_tri;
Scan* sc1, *sc2;
{
    int i, j, k;
    Vertex* vert;
    Triangle* tri;
    Vector pos, pos_world;
    float t;
    int inward;
    int result;
    float dot;
    Vector norm;
    int share_count;
    int found;
    int index;
    Triangle* shared_triangle();
    Vector coord1, coord2;
    float col;
    Vector ct_norm;
    Vector temp_norm;

    /* get cut_tri's normal into mesh 2 coordinates */
    temp_norm[X] = -tri->aa;
    temp_norm[Y] = -tri->bb;
    temp_norm[Z] = -tri->cc;
    mesh_to_world_normal(sc1, temp_norm, ct_norm);
    world_to_mesh_normal(sc2, ct_norm, ct_norm);

    /* transform the endpoints into the mesh 2 coordinate system */
    mesh_to_world(sc1, v1->coord, coord1);
    world_to_mesh(sc2, coord1, coord1);
    mesh_to_world(sc1, v2->coord, coord2);
    world_to_mesh(sc2, coord2, coord2);

    /* See which other triangles (if any) share this edge.  If there */
    /* are others, then we may have already computed the intersections. */

    share_count = edges_shared_count(v1, v2);

    /* Create the "clips" field for all triangles sharing this edge. */

    for (i = 0; i < share_count; i++) {

        tri = shared_triangle(i);

        /* make sure that each triangle has a "clips" field */

        if (tri->clips == NULL) {
            Clipped_Edge* clips;
            clips = (Clipped_Edge*) malloc(sizeof(Clipped_Edge) * 3);
            tri->clips = clips;
            for (k = 0; k < 3; k++) {
                clips[k].v1 = tri->verts[k];
                clips[k].v2 = tri->verts[(k + 1) % 3];
                clips[k].t1 = NULL;
                clips[k].t2 = NULL;
                clips[k].cut_max = 2;
                clips[k].cut_num = 0;
                clips[k].cuts = (Cut**)
                                malloc(sizeof(Cut*) * clips[k].cut_max);
                clips[k].perp_intersect = 1;
                clips[k].done_edge = 0;
            }
        }
    }

    /* find the index of v1 and v2 */

    found = 0;
    for (j = 0; j < 3; j++) {
        if (cut_tri->verts[j] == v1 && cut_tri->verts[(j + 1) % 3] == v2) {
            index = j;
            found = 1;
            break;
        }
    }

    /* sanity check */
    if (!found) {
        fprintf(stderr, "intersect_edge_with_near_tris: can't find vertices\n");
        exit(-1);
    }

    /* if the edge has already been examined, we can leave this routine */
    if (cut_tri->clips[index].done_edge == 1)
        return;

    /* jump point */
we_are_okay:
    ;

    /* mark the edge in question as "examined" in each triangle that it shares */
    for (i = 0; i < share_count; i++) {

        tri = shared_triangle(i);

        found = 0;
        for (j = 0; j < 3; j++) {
            if (tri->verts[j] == v1 && tri->verts[(j + 1) % 3] == v2) {
                tri->clips[j].done_edge = 1;
                found = 1;
                break;
            } else if (tri->verts[j] == v2 && tri->verts[(j + 1) % 3] == v1) {
                tri->clips[j].done_edge = 1;
                found = 1;
                break;
            }
        }

        /* sanity check */
        if (!found) {
            fprintf(stderr, "intersect_edge_with_near_tris: can't find vertices\n");
            exit(-1);
        }
    }

    /* un-mark all nearby triangles, using "eat_mark" as the marker */
    for (i = 0; i < pts_near_num; i++) {
        vert = pts_near[i];
        for (j = 0; j < vert->ntris; j++)
            vert->tris[j]->eat_mark = 0;
    }

    /* test all nearby triangles with the edge */
    for (i = 0; i < pts_near_num; i++) {
        vert = pts_near[i];
        for (j = 0; j < vert->ntris; j++) {

            tri = vert->tris[j];

            /* don't examine the same triangle twice */
            if (tri->eat_mark)
                continue;
            tri->eat_mark = 1;

            /* don't look at triangles that are too nearly parallel to cut_tri */
            dot = -1 * (ct_norm[X] * tri->aa +
                        ct_norm[Y] * tri->bb +
                        ct_norm[Z] * tri->cc);
            if (dot > 0.8)
                continue;

            /* perform intersection */
            result = line_intersect_tri_single(coord1, coord2, tri, pos, &t,
                                               &inward);

            /* save away info on an intersection, if there was one */
            if (result) {
                mesh_to_world(sc2, pos, pos_world);
                add_extra_line(pos_world, pos_world, 0x00ff00);

#if 0
                col = dot;
                if (col > 1)
                    col = 1;
                else if (col < -1)
                    col = -1;
                col = (col + 1) * 0.5;
                v1->confidence = col;
                v2->confidence = col;
#endif

                /* record information about this intersection */

#if 0
                printf("dot = %f\n", dot);
#endif

                new_tri_intersection(v1, v2, share_count,
                                     tri, cut_tri, pos, t, inward, dot);
            }
        }
    }

    /* set all "eat_mark" flags here back to zero */
    for (i = 0; i < pts_near_num; i++) {
        vert = pts_near[i];
        for (j = 0; j < vert->ntris; j++)
            vert->tris[j]->eat_mark = 0;
    }
}


/******************************************************************************
Find the collection of nearby vertices to a given vertex.

Entry:
  mesh     - the mesh containing the (possibly) nearby vertices
  not_mesh - mesh to reject matches from (old_mesh field of Vertex)
  pnt      - position (in mesh's coordinate system) to find nearest vertex to
  norm     - surface normal at pnt
  radius   - distance within which to check

Exit:
  places nearby vertices in "pts_near"
******************************************************************************/

verts_near_vert(mesh, not_mesh, pnt, norm, radius)
Mesh* mesh;
Mesh* not_mesh;
Vector pnt;
Vector norm;
float radius;
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

                    /* go on if this vertex has already been placed on list */
                    /* (count is a marker that says vertex is in pts_near) */
                    if (ptr->count)
                        continue;

                    /* don't look at a vertex if it is used in no triangles */
                    if (ptr->ntris == 0)
                        continue;

                    /* distance (squared) to this point */
                    dx = ptr->coord[X] - pnt[X];
                    dy = ptr->coord[Y] - pnt[Y];
                    dz = ptr->coord[Z] - pnt[Z];
                    dist = dx * dx + dy * dy + dz * dz;

                    /* maybe we've found new closest point */
                    if (dist < radius) {

#if 0
                        /* make sure the surface normals are roughly in the same direction */
                        dot = vdot(norm, ptr->normal);
                        if (dot < 0.3)
                            continue;
#endif

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
Save info about when an edge of one triangle (cut_tri) passes through another
triangle (near_tri).

Entry:
  v1,v2       - vertices that make up the edge that passed through near_tri
  share_count - how many triangles share this edge
  near_tri    - triangle that was intersected
  cut_tri     - triangle that supplied the edge
  pos         - position in 3-space of the intersection
  t           - parameter of the cut along the edge
  inward      - direction that the edge passed through (0 or 1)
  dot         - dot product between the two triangles
******************************************************************************/

new_tri_intersection(v1, v2, share_count, near_tri, cut_tri, pos, t, inward, dot)
Vertex* v1, *v2;
int share_count;
Triangle* near_tri;
Triangle* cut_tri;
Vector pos;
float t;
int inward;
float dot;
{
    int i, j, k;
    More_Tri_Stuff* more;
    Cut* cut;
    Triangle* tri;
    int index;
    int found;
    Clipped_Edge* clips, *clip;
    Triangle* t1, *t2;
    Triangle* shared_triangle();

    /* figure out which triangles share this edge */

    if (share_count == 1) {
        t1 = shared_triangle(0);
        t2 = NULL;
    } else if (share_count == 2) {
        t1 = shared_triangle(0);
        t2 = shared_triangle(1);
    } else {
        t1 = NULL;
        t2 = NULL;
    }

    /* mark all involved triangles as intersected triangles */
    near_tri->dont_touch = 1;
    if (t1 != NULL)
        t1->dont_touch = 1;
    if (t2 != NULL)
        t2->dont_touch = 1;

    /*** create a cut record ***/

    cut = (Cut*) malloc(sizeof(Cut));
    cut->v1 = v1;
    cut->v2 = v2;
    cut->tri = near_tri;
    cut->s = t;
    cut->edge = NULL;
    cut->new_vert = NULL;
    cut->inward = inward;
    cut->dot = dot;


    /*** add this to the list of each triangle that shares this edge ***/

    for (i = 0; i < share_count; i++) {

        tri = shared_triangle(i);
        clips = tri->clips;

        /* find which indices of this triangle point to v1 and v2 */
        /* so that we know in which clip edge we should save this cut */

        found = 0;
        for (j = 0; j < 3; j++) {
            if (tri->verts[j] == v1 && tri->verts[(j + 1) % 3] == v2) {
                index = j;
                found = 1;
                break;
            } else if (tri->verts[j] == v2 && tri->verts[(j + 1) % 3] == v1) {
                index = j;
                found = 1;
                break;
            }
        }

        /* sanity check */
        if (!found) {
            fprintf(stderr, "new_tri_intersection: can't find vertices\n");
            exit(-1);
        }

        /* make sure there is enough room to store cut */
        clip = &clips[index];
        if (clip->cut_num >= clip->cut_max) {
            clip->cut_max += 2;
            clip->cuts = (Cut**)
                         realloc(clip->cuts, sizeof(Cut*) * clip->cut_max);
        }

        /* save info about which triangles share this edge */
        if (tri == t1) {
            clip->t1 = t1;
            clip->t2 = t2;
        } else if (tri == t2) {
            clip->t1 = t2;
            clip->t2 = t1;
        }

        /* store the cut (finally!) */
        clip->cuts[clip->cut_num] = cut;
        clip->cut_num++;
    }


    /*** store a record of this cut in the triangle that was pierced ***/

    /* make sure there is a place to store this info */
    if (near_tri->more == NULL)
        init_tri_intersection(near_tri);

    /* make sure there is room for a new cut */
    more = near_tri->more;
    if (more->cut_num >= more->cut_max) {
        more->cut_max += 4;
        more->cuts = (Cut**) realloc(more->cuts, sizeof(Cut*) * more->cut_max);
    }

    /* add to the list of cuts */
    more->cuts[more->cut_num] = cut;
    more->cut_num++;


    debug_cut_count++;

}


/******************************************************************************
Initialize the cut list of a triangle.

Entry:
  tri - triangle to initialize
******************************************************************************/

init_tri_intersection(tri)
Triangle* tri;
{
    More_Tri_Stuff* more;

    more = (More_Tri_Stuff*) malloc(sizeof(More_Tri_Stuff));
    tri->more = more;

    more->cut_max = 4;
    more->cut_num = 0;
    more->cuts = (Cut**) malloc(sizeof(Cut*) * more->cut_max);
    more->clip_count = 0;
    more->clip_flag = 0;
}


/******************************************************************************
Add new vertices to the mesh where intersections occured.

Entry:
  sc1,sc2 - the two meshes that were intersected
******************************************************************************/

add_intersect_points(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j;
    Mesh* m1, *m2;
    Triangle* tri;
    More_Tri_Stuff* more;
    Cut* cut;
    Vector pos;
    float* p1, *p2;
    float t;
    int index;
    Vertex* vert;
    int count = 0;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    /* find all intersection points by looking at the list of cut */
    /* points that are stored at each triangle */

    for (i = 0; i < m1->ntris; i++) {

        tri = m1->tris[i];

        /* see if there might be any info on intersections at this triangle */
        more = tri->more;
        if (more == NULL)
            continue;

        /* make sure we're not examining wrong triangle */
        if (tri->dont_touch == 0 && more->cut_num > 0)
            fprintf(stderr,
                    "add_intersect_points: ooooooops! we're looking at wrong tri\n");

        /* introduce all the cuts from this triangle into the mesh */
        for (j = 0; j < more->cut_num; j++) {

            cut = more->cuts[j];

            /* compute position of intersection */
            t = cut->s;
            p1 = cut->v1->coord;
            p2 = cut->v2->coord;
            pos[X] = p1[X] + t * (p2[X] - p1[X]);
            pos[Y] = p1[Y] + t * (p2[Y] - p1[Y]);
            pos[Z] = p1[Z] + t * (p2[Z] - p1[Z]);

            /* create the new vertex */
            index = make_vertex(m1, pos);
            vert = m1->verts[index];
            add_to_hash(vert, m1);

            /* save a pointer to this vertex */
            cut->new_vert = vert;
            count++;
        }
    }

    printf("%d new vertices added\n", count);
}


/******************************************************************************
Perform the actual clipping of intersection triangles.

Entry:
  sc1,sc2 - the two meshes that were intersected
******************************************************************************/

perform_intersect_clipping(sc1, sc2)
Scan* sc1, *sc2;
{
    int i, j, k;
    Mesh* m1, *m2;
    Triangle* tri;
    int cut_count;
    Clip_List* clist;
    Clip_List* potential_vertices();
    Clip_Vertex* v;
    int in_vert, out_vert;
    int result;
    int self_intersect;
    int p1, p2, p3;
    Vector vec;
    Triangle* ntri;
    Cut* cut;
    More_Tri_Stuff* more;
    static Vertex* between_list[30];
    int between_count;
    Vertex** cv;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    /* find all intersection points by looking at the list of cut */
    /* points that are stored at each triangle */

    for (i = 0; i < m1->ntris; i++) {

        tri = m1->tris[i];
        if (tri->clips == NULL)
            continue;

        /* see how many places this triangle has been cut */
        cut_count = tri->clips[0].cut_num + tri->clips[1].cut_num +
                    tri->clips[2].cut_num;
        if (cut_count == 0)
            continue;

        /* sort the cuts along the triangle's edges */
        sort_triangle_cuts(tri);

        /* create list of potential vertices from triangle to be clipped */
        clist = potential_vertices(tri);

        /* examine the cuts to see which vertices should be retained in */
        /* the clipped polygon */

        v = clist->list;

#if 0
        print_list("tri: ", clist);
#endif

        in_vert = out_vert = -1;
        for (j = 0; j < clist->count; j++) {
            if (v[j].type != CUT)
                continue;
            if (v[j].inward) {
                in_vert = j;
            } else {
                out_vert = j;
            }
        }

        if (in_vert == -1 || out_vert == -1) {
            printf("can't find in/out match\n");
            continue;
        }

        /* determine which cut points on the *other* mesh should be */
        /* added to this clipped polygon */
        between_count = between_cuts(tri, clist, in_vert, out_vert, between_list);

        /* create a list of vertices that will eventually replace this triangle */
        if (tri->more == NULL)
            init_tri_intersection(tri);

        more = tri->more;
        more->clip_flag = 1;

        /* count how many vertices there should be */
        more->clip_count = 1;  /* so we also count in_vert */
        for (j = out_vert; j != in_vert; j = (j + 1) % clist->count)
            more->clip_count++;
        more->clip_count += between_count;

        /* create clip list */
        more->clip_verts = (Vertex**)
                           malloc(sizeof(Vertex*) * more->clip_count);
        k = 0;
        for (j = out_vert; j != in_vert; j = (j + 1) % clist->count) {
            more->clip_verts[k] = v[j].vert;
            k++;
        }
        more->clip_verts[k] = v[in_vert].vert;
        k++;
        for (j = 0; j < between_count; j++) {
            more->clip_verts[k] = between_list[j];
            k++;
        }

        /* consistancy check */
        if (k != more->clip_count) {
            fprintf(stderr,
                    "perform_intersect_clipping: number of vertices don't match\n");
            exit(-1);
        }

    }

    /* actually do the clipping by replacing the old triangle with a */
    /* clipped version */

    for (i = 0; i < m1->ntris; i++) {

        tri = m1->tris[i];
        if (tri->clips == NULL || tri->more == NULL || !tri->more->clip_flag)
            continue;

        more = tri->more;
        cv = more->clip_verts;

        /* initialize the polygon splitter */
        result = init_splitter(-tri->aa, -tri->bb, -tri->cc, 0.0);
        if (result)
            continue;

        /* send the vertices to the splitter */
        for (j = 0; j < more->clip_count; j++) {
            vcopy(cv[j]->coord, vec);
            add_boundary_point(vec[X], vec[Y], vec[Z], j);
        }

        /* call the splitter */
        self_intersect = greedy_connect();

        /* delete the original triangle */
        delete_triangle(tri, m1, 0);

        /* create the clipped triangle */
        if (!self_intersect) {

            for (j = 0; j < get_ntris(); j++) {
                get_triangle(j, &p1, &p2, &p3);
                if (check_proposed_tri(cv[p1], cv[p2], cv[p3])) {
                    ntri = make_triangle(m1, cv[p1], cv[p2], cv[p3], 1e20);
                }
            }

            /* re-compute normals and edge conditions around the new triangles */
            for (j = 0; j < more->clip_count; j++) {
                vertex_edge_test(cv[j]);
                find_vertex_normal(cv[j]);
            }
        }

        /* free up the vertex list "clip_verts" ??? */
        /* free up the vertex list "clip_verts" ??? */
        /* free up the vertex list "clip_verts" ??? */
        /* free up the vertex list "clip_verts" ??? */

    }

}


/******************************************************************************
Find which vertices should be strung between two cut points of a triangle
that is being clipped.

Entry:
  tri      - triangle being clipped
  clist    - list of vertices around clipped triangle (so far)
  in_vert  - index of one vertex that we're trying to join
  out_vert - index of other vertex

Exit:
  between_list - list of vertices at the cut points that we're finding
  returns number of elements in between_list
******************************************************************************/

int between_cuts(tri, clist, in_vert, out_vert, between_list)
Triangle* tri;
Clip_List* clist;
int in_vert;
int out_vert;
Vertex* between_list[];
{
    int i, j;
    Vertex* v1, *v2;
    Clip_Vertex* v = clist->list;
    int count;
    More_Tri_Stuff* more;
    Cut* cut;
    Cut* cut1, *cut2;
    Triangle* tri_first, *tri_last;
    Cut* cut_list[10];
    Triangle* tri_list[10];
    Triangle* tri_cur, *tri_old;
    int cut_count;
    int error_flag;
    int found;
    int index;
    int chain_count = 0;

    cut1 = v[in_vert].cut;
    cut2 = v[out_vert].cut;
    v1 = v[in_vert].vert;
    v2 = v[out_vert].vert;

    more = tri->more;
    if (more == NULL || more->cut_num == 0)
        return (0);

#if 0
    printf("%d intersections\n", more->cut_num);
    for (i = 0; i < more->cut_num; i++) {
        cut = more->cuts[i];
    }
#endif

    tri_first = cut1->tri;
    tri_last = cut2->tri;

    tri_old = NULL;
    tri_cur = tri_first;
    error_flag = 0;

    while (tri_cur != tri_last) {

        /* find where the current triangle cuts through "tri" */
        cut_count = tri_cut_by_tri(tri, tri_cur, cut_list, tri_list);

        /* we've got an error if we didn't find any such cuts */
        if (cut_count == 0) {
            printf("cut_count == 0\n");
            error_flag = 1;
            break;
        }

        found = 0;
        for (i = 0; i < cut_count; i++) {
            if (tri_list[i] != NULL && tri_list[i] != tri_old) {
                index = i;
                found++;
            }
        }

        /* we should have found exactly one such triangle */
        if (found != 1) {
            printf("found = %d\n", found);
            error_flag = 1;
            break;
        }

        /* move on to find the next cut */
        tri_old = tri_cur;
        tri_cur = tri_list[index];
        between_list[chain_count] = cut_list[index]->new_vert;
        chain_count++;
    }

    if (error_flag) {
        printf("intersect error\n");
    }

    if (chain_count != more->cut_num)
        printf("chain_count intersections: %d %d\n", chain_count, more->cut_num);

    return (chain_count);
}


/******************************************************************************
Did an edge from one triangle (tri2) intersect another triangle (tri1)?
If so, return the cut place(s).

Entry:
  tri1 - triangle that may have been intersected
  tri2 - triangle whose edge may have passed through another triangle

Exit:
  cuts - list of places of intersect (if any)
  tris - list of triangles that share the cut edges with tri2
  returns number of such intersections
******************************************************************************/

int tri_cut_by_tri(tri1, tri2, cuts, tris)
Triangle* tri1, *tri2;
Cut* cuts[];
Triangle* tris[];
{
    int i, j;
    Clipped_Edge* clips, *clip;
    Cut* cut;
    int cut_count = 0;
    Triangle* other_tri;

    clips = tri2->clips;
    if (clips == NULL)
        return (0);

    /* look along each edge of tri2 for any cut points where tri2's edge */
    /* intersected tri1 */

    for (i = 0; i < 3; i++) {

        clip = &clips[i];

        /* look at all cuts along this edge for the approprate intersection */

        for (j = 0; j < clip->cut_num; j++) {

            cut = clip->cuts[j];

            /* determine what other triangle (if any) shares this edge with tri2 */
            if (clip->t1 == tri2)
                other_tri = clip->t2;
            else if (clip->t2 == tri2)
                other_tri = clip->t1;
            else
                other_tri = NULL;

            if (cut->tri == tri1) {
                cuts[cut_count] = cut;
                tris[cut_count] = other_tri;
                cut_count++;
            }
        }
    }

    /* return the number of matching cuts found */
    return (cut_count);
}

