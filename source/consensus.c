/*

Determine consensus position of mesh surface as a weighted average of
what several scans say.

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
#include "matrix.h"
#include "zipglobal.h"


typedef struct Cvert {  /* vertex to help find consensus */
    Vector coord;
    struct Ctri** tris;
    unsigned char ntris;
    struct Cvert* next;
} Cvert;

typedef struct Ctri {   /* triangle to help find consensus */
    Cvert* verts[3];
} Ctri;

typedef struct Cinfo {
    Vector pos;       /* consensus position */
    Vector normal;    /* consensus surface normal (in global coords) */
    float intensity;      /* consensus intensity */
    float weights;
    float red, grn, blu;  /* color */
    int count;
} Cinfo;

static float CONSENSUS_POSITION_DIST_FACTOR;
static float CONSENSUS_POSITION_DIST;

static float CONSENSUS_NORMAL_DIST_FACTOR;
static float CONSENSUS_NORMAL_DIST;

static float CONSENSUS_JITTER_DIST_FACTOR;
static float CONSENSUS_JITTER_DIST;

Vertex* found_vert_near_vert(int n);

update_consensus_resolution()
{
    CONSENSUS_POSITION_DIST =
        ZIPPER_RESOLUTION * CONSENSUS_POSITION_DIST_FACTOR;

    CONSENSUS_NORMAL_DIST =
        ZIPPER_RESOLUTION * CONSENSUS_NORMAL_DIST_FACTOR;

    CONSENSUS_JITTER_DIST =
        ZIPPER_RESOLUTION * CONSENSUS_JITTER_DIST_FACTOR;
}


set_consensus_position_dist_factor(factor)
float factor;
{
    CONSENSUS_POSITION_DIST_FACTOR = factor;

    CONSENSUS_POSITION_DIST =
        ZIPPER_RESOLUTION * CONSENSUS_POSITION_DIST_FACTOR;
}


float
get_consensus_position_dist_factor()
{
    return CONSENSUS_POSITION_DIST_FACTOR;
}


set_consensus_normal_dist_factor(factor)
float factor;
{
    CONSENSUS_NORMAL_DIST_FACTOR = factor;

    CONSENSUS_NORMAL_DIST =
        ZIPPER_RESOLUTION * CONSENSUS_NORMAL_DIST_FACTOR;
}


float
get_consensus_normal_dist_factor()
{
    return CONSENSUS_NORMAL_DIST_FACTOR;
}


set_consensus_jitter_dist_factor(factor)
float factor;
{
    CONSENSUS_JITTER_DIST_FACTOR = factor;

    CONSENSUS_JITTER_DIST =
        ZIPPER_RESOLUTION * CONSENSUS_JITTER_DIST_FACTOR;
}


float
get_consensus_jitter_dist_factor()
{
    return CONSENSUS_JITTER_DIST_FACTOR;
}


/******************************************************************************
Create an "average" surface by weighted average of multiple scans.

Entry:
  scan    - scan to average out
  level   - level of mesh spacing (0 = every point, 1 = every other, 2 = every
        four, 3 = every eight)
  k_scale - scaling coefficient for search distance
******************************************************************************/

consensus_surface(scan, level, k_scale)
Scan* scan;
int level;
float k_scale;
{
    int i, j;
    Mesh* mesh;
    Cinfo* cinfo;
    Triangle* tri;
    int zeros;
    Vertex* v;

    mesh = scan->meshes[mesh_level];

    init_extra_lines();

    /* initialize the consensus geometry information at each vertex */

    for (i = 0; i < mesh->nverts; i++) {
        cinfo = (Cinfo*) malloc(sizeof(Cinfo));
        mesh->verts[i]->cinfo = cinfo;
        vset(cinfo->pos, 0.0, 0.0, 0.0);
        vset(cinfo->normal, 0.0, 0.0, 0.0);
        cinfo->weights = 0;
        cinfo->intensity = 0;
        cinfo->count = 0;
        cinfo->red = cinfo->grn = cinfo->blu = 0;

        mesh->verts[i]->old_mesh = 0;  /* kluge for mesh tags */

    }

    /* determine "average" position for vertices, based on all meshes */

    marc_find_average_positions(scan, level, k_scale);
    /*
    find_average_positions (scan, level, k_scale);
    */

    /* use the averaging information to move the vertices */

    zeros = 0;
    for (i = 0; i < mesh->nverts; i++) {
        float s;

        v = mesh->verts[i];
        cinfo = v->cinfo;

        /* go on to next vertex if this one had NO near positions */
        if (cinfo->count == 0 || cinfo->weights == 0) {
            zeros++;
            continue;
        }

        /* average the nearby positions */
        /*
        s = 1.0 / cinfo->count;
        */
        s = 1.0 / cinfo->weights;
        cinfo->pos[X] *= s;
        cinfo->pos[Y] *= s;
        cinfo->pos[Z] *= s;

        /* move vertex to the average position */
        remove_from_hash(v, mesh);
        world_to_mesh(scan, cinfo->pos, v->coord);
        add_to_hash(v, mesh);
    }

    /* re-compute triangle normals and edge planes */
    for (i = 0; i < mesh->ntris; i++) {
        tri = mesh->tris[i];
        plane_thru_vectors(tri->verts[0], tri->verts[1], tri->verts[2],
                           &tri->aa, &tri->bb, &tri->cc, &tri->dd);
        compute_edge_planes(tri);
    }

    /* re-compute vertex normals */
    find_vertex_normals(mesh);

    /* compute the weighted average of confidences and intensities */
    for (i = 0; i < mesh->nverts; i++) {
        v = mesh->verts[i];
        if (v->cinfo->count == 0 || v->cinfo->weights == 0) {
            v->confidence = 0.0;
            v->intensity = 0.0;
            v->red = v->grn = v->blu = 0;
        } else {
            v->confidence = v->cinfo->weights / (float) v->cinfo->count;
            v->intensity = v->cinfo->intensity / v->cinfo->weights;
            v->red = v->cinfo->red / v->cinfo->weights;
            v->grn = v->cinfo->grn / v->cinfo->weights;
            v->blu = v->cinfo->blu / v->cinfo->weights;
        }
    }

    /* free the vertex geometry info */
    for (i = 0; i < mesh->nverts; i++) {
        v = mesh->verts[i];
        free(v->cinfo);
        v->cinfo = NULL;
    }
}


/******************************************************************************
Create an "average" surface by weighted average of multiple scans.

Entry:
******************************************************************************/

new_consensus_surface(con_scan, scan_list, read_list, num_scans, level)
Scan* con_scan;
Scan** scan_list;
int* read_list;
int num_scans, level;
{
    int i, j;
    Mesh* mesh;
    Cinfo* cinfo;
    Triangle* tri;
    int zeros;
    Vertex* v;

    mesh = con_scan->meshes[level];

    init_extra_lines();

    /* initialize the consensus geometry information at each vertex */

    for (i = 0; i < mesh->nverts; i++) {
        cinfo = (Cinfo*) malloc(sizeof(Cinfo));
        mesh->verts[i]->cinfo = cinfo;
        vset(cinfo->pos, 0.0, 0.0, 0.0);
        vset(cinfo->normal, 0.0, 0.0, 0.0);
        cinfo->weights = 0;
        cinfo->intensity = 0;
        cinfo->count = 0;
        cinfo->red = cinfo->grn = cinfo->blu = 0;

        mesh->verts[i]->old_mesh = 0;  /* kluge for mesh tags */

    }

    /* determine "average" position for vertices, based on all meshes */

    new_find_average_positions(con_scan, scan_list, read_list, num_scans, level);
    /*
    find_average_positions (scan, level, k_scale);
    */

    /* use the averaging information to move the vertices */

    zeros = 0;
    for (i = 0; i < mesh->nverts; i++) {
        float s;

        v = mesh->verts[i];
        cinfo = v->cinfo;

        /* go on to next vertex if this one had NO near positions */
        if (cinfo->count == 0 || cinfo->weights == 0) {
            zeros++;
            continue;
        }

        /* average the nearby positions */
        /*
        s = 1.0 / cinfo->count;
        */
        s = 1.0 / cinfo->weights;
        cinfo->pos[X] *= s;
        cinfo->pos[Y] *= s;
        cinfo->pos[Z] *= s;

        /* move vertex to the average position */
        remove_from_hash(v, mesh);
        world_to_mesh(con_scan, cinfo->pos, v->coord);
        add_to_hash(v, mesh);
    }

    /* re-compute triangle normals and edge planes */
    for (i = 0; i < mesh->ntris; i++) {
        tri = mesh->tris[i];
        plane_thru_vectors(tri->verts[0], tri->verts[1], tri->verts[2],
                           &tri->aa, &tri->bb, &tri->cc, &tri->dd);
        compute_edge_planes(tri);
    }

    /* re-compute vertex normals */
    find_vertex_normals(mesh);

    /* compute the weighted average of confidences and intensities */
    for (i = 0; i < mesh->nverts; i++) {
        v = mesh->verts[i];
        if (v->cinfo->count == 0 || v->cinfo->weights == 0) {
            v->confidence = 0.0;
            v->intensity = 0.0;
            v->red = v->grn = v->blu = 0;
        } else {
            v->confidence = v->cinfo->weights / (float) v->cinfo->count;
            v->intensity = v->cinfo->intensity / v->cinfo->weights;
            v->red = v->cinfo->red / v->cinfo->weights;
            v->grn = v->cinfo->grn / v->cinfo->weights;
            v->blu = v->cinfo->blu / v->cinfo->weights;
        }
    }

    /* free the vertex geometry info */
    for (i = 0; i < mesh->nverts; i++) {
        v = mesh->verts[i];
        free(v->cinfo);
        v->cinfo = NULL;
    }
}


/******************************************************************************
Find the "average" positions over all scans to the vertices in a given mesh.
This version uses Marc's suggestion of first coming up with an average
normal.

Entry:
  scan    - scan containing the given mesh
  level   - level of mesh spacing (0 = every point, 1 = every other, 2 = every
            four, 3 = every eight)
  k_scale - scaling coefficient for search distance
******************************************************************************/

marc_find_average_positions(scan, level, k_scale)
Scan* scan;
int level;
float k_scale;
{
    int i, j, k, n;
    Mesh* mesh, *tmesh;
    NearPosition near_info;
    Vector pos, norm;
    int result;
    int zeros;
    int spacing;
    float search_dist;
    float normal_dist;
    Vertex* v;
    Mesh* make_mesh();
    Mesh* make_mesh_raw();
    Mesh* make_mesh_ply();
    Vertex* near_vert;
    Vertex* found_vert_near_vert();
    int count;
    Vector diff;
    float len;
    int mesh_index;

    mesh = scan->meshes[mesh_level];

    /* determine spacing between mesh elements */
    spacing = level_to_inc(level);

    /* use this to bound the search distance */
    search_dist = CONSENSUS_POSITION_DIST * spacing;  /* intersection search */
    normal_dist = CONSENSUS_NORMAL_DIST * spacing;    /* normal search */

    printf("using normal search distance of %f\n", normal_dist);
    printf("using intersection search distance of %f\n", search_dist);

    /* Examine each mesh to see which points on these meshes are */
    /* close to the vertices of the given mesh.  We're going to */
    /* arrive at a consensus normal first. */

    for (i = 0; i < nscans; i++) {

        /* don't bother with meshes read from a polygon file */
        if (scans[i]->file_type == POLYFILE)
            continue;

        /* create a mesh of the appropriate level */
        if (scans[i]->file_type == CYFILE)
            tmesh = make_mesh(scans[i], level, normal_dist);
        else if (scans[i]->file_type == RAWFILE)
            tmesh = make_mesh_raw(scans[i], level, normal_dist);
        else if (scans[i]->file_type == PLYRANGEFILE)
            tmesh = make_mesh_ply(scans[i], level, normal_dist);
        else {
            fprintf(stderr, "consensus_surface: bad scan type: %d\n",
                    scans[i]->file_type);
            continue;
        }

        printf("mesh has %d vertices and %d triangles\n", tmesh->nverts, tmesh->ntris);

        /* do neighbor search around vertices */

        for (j = 0; j < mesh->nverts; j++) {

            v = mesh->verts[j];

            /* find vertex position in "tmesh" coordinates */
            mesh_to_world(scan, v->coord, pos);
            world_to_mesh(scans[i], pos, pos);
            mesh_to_world_normal(scan, v->normal, norm);
            world_to_mesh_normal(scans[i], norm, norm);

            /* find nearby vertices */
            verts_near_pos(tmesh, pos, norm, normal_dist);
            count = count_near_vert();

            /* compute consensus normal */
            for (k = 0; k < count; k++) {
                near_vert = found_vert_near_vert(k);
                vsub(pos, near_vert->coord, diff);
                mesh_to_world(scans[i], near_vert->normal, norm);
                vadd(v->cinfo->normal, norm, v->cinfo->normal);
            }

#if 0
            if (j % 100 == 0) {
                printf("vertex %d: count = %d, norm = %f %f %f\n", j, count,
                       v->cinfo->normal[X], v->cinfo->normal[Y], v->cinfo->normal[Z]);
            }
#endif

        }

        /* free the mesh info */
        clear_mesh(tmesh);
    }

    /* scale all those consensus normals by the search distance */

    for (j = 0; j < mesh->nverts; j++) {
        v = mesh->verts[j];
        vnorm(v->cinfo->normal);
        v->cinfo->normal[X] *= search_dist;
        v->cinfo->normal[Y] *= search_dist;
        v->cinfo->normal[Z] *= search_dist;
    }

    /* now look in this normal direction to see where we intersect */
    /* with these other meshes */

    mesh_index = 0;
    for (i = 0; i < nscans; i++) {

        /* don't bother with meshes read from a polygon file */
        if (scans[i]->file_type == POLYFILE)
            continue;

        /* create a mesh of the appropriate level */
        if (scans[i]->file_type == CYFILE)
            tmesh = make_mesh(scans[i], level, search_dist);
        else if (scans[i]->file_type == RAWFILE)
            tmesh = make_mesh_raw(scans[i], level, search_dist);
        else if (scans[i]->file_type == PLYRANGEFILE)
            tmesh = make_mesh_ply(scans[i], level, search_dist);
        else {
            fprintf(stderr, "consensus_surface: bad scan type: %d\n",
                    scans[i]->file_type);
            continue;
        }

        printf("mesh has %d vertices and %d triangles\n", tmesh->nverts, tmesh->ntris);

        /* do neighbor search around vertices */

        for (j = 0; j < mesh->nverts; j++) {

            v = mesh->verts[j];

            /* intersect line segment through "v" with "tmesh", adding */
            /* the intersection info to the consensus record of "v" */
            intersect_segment_with_mesh(v, tmesh, scan, scans[i], search_dist,
                                        mesh_index);
        }

        /* if we get here, we need to increment the mesh index */
        mesh_index++;

        /* free the mesh info */
        clear_mesh(tmesh);
    }
}


/******************************************************************************
Find the "average" positions over all scans to the vertices in a given mesh.
This version uses Marc's suggestion of first coming up with an average
normal.

Entry:
******************************************************************************/

new_find_average_positions(con_scan, scan_list, use_old_mesh, num_scans, level)
Scan* con_scan;
Scan** scan_list;
int* use_old_mesh;
int num_scans, level;
{
    int i, j, k, n;
    Mesh* mesh, *tmesh;
    NearPosition near_info;
    Vector pos, norm;
    int result;
    int zeros;
    int spacing;
    float search_dist;
    float normal_dist;
    Vertex* v;
    Mesh* make_mesh();
    Mesh* make_mesh_raw();
    Mesh* make_mesh_ply();
    Vertex* near_vert;
    Vertex* found_vert_near_vert();
    int count;
    Vector diff;
    float len;
    int mesh_index;
    float old_size;

    mesh = con_scan->meshes[level];

    /* determine spacing between mesh elements */
    spacing = level_to_inc(level);

    /* use this to bound the search distance */
    search_dist = CONSENSUS_POSITION_DIST * spacing;  /* intersection search */
    normal_dist = CONSENSUS_NORMAL_DIST * spacing;    /* normal search */

    printf("using normal search distance of %f\n", normal_dist);
    printf("using intersection search distance of %f\n", search_dist);

    /* Examine each mesh to see which points on these meshes are */
    /* close to the vertices of the given mesh.  We're going to */
    /* arrive at a consensus normal first. */

    for (i = 0; i < num_scans; i++) {

        if (!use_old_mesh[i] && scan_list[i]->file_type != POLYFILE) {
            /* don't bother with meshes read from a polygon file */
            if (scan_list[i]->file_type == POLYFILE)
                continue;

            /* create a mesh of the appropriate level */
            if (scan_list[i]->file_type == CYFILE)
                tmesh = make_mesh(scan_list[i], level, normal_dist);
            else if (scan_list[i]->file_type == RAWFILE)
                tmesh = make_mesh_raw(scan_list[i], level, normal_dist);
            else if (scan_list[i]->file_type == PLYRANGEFILE)
                tmesh = make_mesh_ply(scan_list[i], level, normal_dist);
            else {
                fprintf(stderr, "consensus_surface: bad scan type: %d\n",
                        scan_list[i]->file_type);
                continue;
            }
        } else {
            tmesh = scan_list[i]->meshes[level];

            /* Delete old hash table, remember old cell size */
            old_size = 1 / tmesh->table->scale;
            free(tmesh->table);
            tmesh->table = NULL;

            /* Make new table */
            init_table(tmesh, normal_dist);
        }

        printf("mesh has %d vertices and %d triangles\n", tmesh->nverts, tmesh->ntris);

        /* do neighbor search around vertices */

        for (j = 0; j < mesh->nverts; j++) {

            v = mesh->verts[j];

            /* find vertex position in "tmesh" coordinates */
            mesh_to_world(con_scan, v->coord, pos);
            world_to_mesh(scan_list[i], pos, pos);
            mesh_to_world_normal(con_scan, v->normal, norm);
            world_to_mesh_normal(scan_list[i], norm, norm);

            /* find nearby vertices */
            verts_near_pos(tmesh, pos, norm, normal_dist);
            count = count_near_vert();

            /* compute consensus normal */
            for (k = 0; k < count; k++) {
                near_vert = found_vert_near_vert(k);
                vsub(pos, near_vert->coord, diff);
                mesh_to_world(scan_list[i], near_vert->normal, norm);
                vadd(v->cinfo->normal, norm, v->cinfo->normal);
            }

#if 0
            if (j % 100 == 0) {
                printf("vertex %d: count = %d, norm = %f %f %f\n", j, count,
                       v->cinfo->normal[X], v->cinfo->normal[Y], v->cinfo->normal[Z]);
            }
#endif

        }

        /* free the mesh info */
        if (!use_old_mesh[i] && scan_list[i]->file_type != POLYFILE)
            clear_mesh(tmesh);
        else {
            if (tmesh->table->verts != NULL)
                free(tmesh->table->verts);

            /* Re-build old hash table */
            init_table(tmesh, old_size);
        }
    }

    /* scale all those consensus normals by the search distance */

    for (j = 0; j < mesh->nverts; j++) {
        v = mesh->verts[j];
        vnorm(v->cinfo->normal);
        v->cinfo->normal[X] *= search_dist;
        v->cinfo->normal[Y] *= search_dist;
        v->cinfo->normal[Z] *= search_dist;
    }

    /* now look in this normal direction to see where we intersect */
    /* with these other meshes */

    mesh_index = 0;
    for (i = 0; i < num_scans; i++) {

        if (!use_old_mesh[i] && scan_list[i]->file_type != POLYFILE) {
            /* don't bother with meshes read from a polygon file */
            if (scan_list[i]->file_type == POLYFILE)
                continue;

            /* create a mesh of the appropriate level */
            if (scan_list[i]->file_type == CYFILE)
                tmesh = make_mesh(scan_list[i], level, search_dist);
            else if (scan_list[i]->file_type == RAWFILE)
                tmesh = make_mesh_raw(scan_list[i], level, search_dist);
            else if (scan_list[i]->file_type == PLYRANGEFILE)
                tmesh = make_mesh_ply(scan_list[i], level, search_dist);
            else {
                fprintf(stderr, "consensus_surface: bad scan type: %d\n",
                        scan_list[i]->file_type);
                continue;
            }
        } else {
            tmesh = scan_list[i]->meshes[level];

            /* Delete old hash table, remember old cell size */
            old_size = 1 / tmesh->table->scale;
            free(tmesh->table);
            tmesh->table = NULL;

            /* Make new table */
            init_table(tmesh, search_dist);
        }


        printf("mesh has %d vertices and %d triangles\n", tmesh->nverts, tmesh->ntris);

        /* do neighbor search around vertices */

        for (j = 0; j < mesh->nverts; j++) {

            v = mesh->verts[j];

            /* intersect line segment through "v" with "tmesh", adding */
            /* the intersection info to the consensus record of "v" */
            intersect_segment_with_mesh(v, tmesh, con_scan, scan_list[i],
                                        search_dist, mesh_index);
        }

        /* if we get here, we need to increment the mesh index */
        mesh_index++;

        /* free the mesh info */
        if (!use_old_mesh[i] && scan_list[i]->file_type != POLYFILE)
            clear_mesh(tmesh);
        else {
            if (tmesh->table->verts != NULL)
                free(tmesh->table->verts);

            /* Re-build old hash table */
            init_table(tmesh, old_size);
        }
    }
}


/******************************************************************************
Intersect the line segment at a vertex that is pointing in the surface
normal's direction with a mesh.

Entry:
  v           - vertex that we construct line segment through
  mesh        - mesh to intersect this line segment with
  vscan       - scan containing "v"
  mscan       - scan containing "mesh"
  search_dist - distance to extend the line segment
  mesh_index  - mesh index to use in mesh tags

Exit:
  adds the nearest intersection to the vertice's consensus info
******************************************************************************/

intersect_segment_with_mesh(v, mesh, vscan, mscan, search_dist, mesh_index)
Vertex* v;
Mesh* mesh;
Scan* vscan;
Scan* mscan;
float search_dist;
int mesh_index;
{
    int i, j, k;
    Vertex* near_vert;
    Vector pos, norm;
    Vector wpos;
    Vector near_pos;
    Triangle* tri;
    Cinfo* cinfo;
    Vector end1, end2;
    float t, tmin;
    int result;
    int found;
    int count;
    int in;
    float small;
    NearPosition near_info;
    float conf, nearest_conf;
    float intensity, nearest_intensity;
    Vector barycentric;
    float bary_sum;
    float red, grn, blu;
    Vertex* v1, *v2, *v3;

    cinfo = v->cinfo;

    /* transform v's position and normal into "mesh" coordinates */
    mesh_to_world(vscan, v->coord, wpos);
    world_to_mesh(mscan, wpos, pos);
    mesh_to_world_normal(vscan, v->normal, norm);
    world_to_mesh_normal(mscan, norm, norm);

#if 1
    /* jitter the position of the vertex */
    wpos[X] += (drand48() - 0.5) * CONSENSUS_JITTER_DIST;
    wpos[Y] += (drand48() - 0.5) * CONSENSUS_JITTER_DIST;
    wpos[Z] += (drand48() - 0.5) * CONSENSUS_JITTER_DIST;
#endif

    /* build a line segment in the consensus normal direction */
    /* (cinfo->pos contains global location of "v", computed above) */

    end1[X] = wpos[X] + cinfo->normal[X];
    end1[Y] = wpos[Y] + cinfo->normal[Y];
    end1[Z] = wpos[Z] + cinfo->normal[Z];

    end2[X] = wpos[X] - cinfo->normal[X];
    end2[Y] = wpos[Y] - cinfo->normal[Y];
    end2[Z] = wpos[Z] - cinfo->normal[Z];

#if 0
    if (mscan == scans[0])
        add_extra_line(end1, end2, 0x00ff00);
#endif

    /* transform the endpoints into the mesh's local coordinates */
    world_to_mesh(mscan, end1, end1);
    world_to_mesh(mscan, end2, end2);

    /* find out which points of "mesh" are near "v" */
    verts_near_pos(mesh, pos, norm, search_dist);
    count = count_near_vert();

    /* mark all nearby triangles as not looked at (mark = 0) */
    for (i = 0; i < count; i++) {
        near_vert = found_vert_near_vert(i);
        for (j = 0; j < near_vert->ntris; j++)
            near_vert->tris[j]->mark = 0;
    }

    /* intersect this segment with any nearby triangles, looking */
    /* for the nearest intersection */

    tmin = 1e20;
    found = 0;

    for (i = 0; i < count; i++) {
        near_vert = found_vert_near_vert(i);
        for (j = 0; j < near_vert->ntris; j++) {
            tri = near_vert->tris[j];
            /* don't look at a triangle again */
            if (tri->mark)
                continue;
            tri->mark = 1;

            /* perform intersection test */

#if 0
            /* single-precision intersection */
            result = line_intersect_tri_single(end1, end2, tri, pos, &t, &in,
                                               barycentric);
#endif

#if 1
            /* double-precision intersection */
            if (tri->more == NULL)
                double_stuff(tri);
            result = line_intersect_tri(end1, end2, tri, pos, &t, &in, barycentric);
#endif

            /* if we get an intersection, remember it if it is nearest yet */
            if (result) {
                t = fabs(t - 0.5);
                if (t < tmin) {
                    t = tmin;
                    found = 1;
                    mesh_to_world(mscan, pos, near_pos);
                    bary_sum = barycentric[X] + barycentric[Y] + barycentric[Z];
                    nearest_conf = (barycentric[X] * tri->verts[0]->confidence +
                                    barycentric[Y] * tri->verts[1]->confidence +
                                    barycentric[Z] * tri->verts[2]->confidence) /
                                   bary_sum;
                    nearest_intensity = (barycentric[X] * tri->verts[0]->intensity +
                                         barycentric[Y] * tri->verts[1]->intensity +
                                         barycentric[Z] * tri->verts[2]->intensity) /
                                        bary_sum;
                    red = (barycentric[X] * tri->verts[0]->red +
                           barycentric[Y] * tri->verts[1]->red +
                           barycentric[Z] * tri->verts[2]->red) / bary_sum;
                    grn = (barycentric[X] * tri->verts[0]->grn +
                           barycentric[Y] * tri->verts[1]->grn +
                           barycentric[Z] * tri->verts[2]->grn) / bary_sum;
                    blu = (barycentric[X] * tri->verts[0]->blu +
                           barycentric[Y] * tri->verts[1]->blu +
                           barycentric[Z] * tri->verts[2]->blu) / bary_sum;
                }
            }

        }
    }

    /* if we found a nearby intersection, average it into the vertex position */

    if (found) {

#if 0
        vadd(cinfo->pos, near_pos, cinfo->pos);
        cinfo->weights += 1;
#endif

        cinfo->pos[X] += nearest_conf * near_pos[X];
        cinfo->pos[Y] += nearest_conf * near_pos[Y];
        cinfo->pos[Z] += nearest_conf * near_pos[Z];
        cinfo->intensity += nearest_intensity * nearest_conf;
        cinfo->red += nearest_conf * red;
        cinfo->grn += nearest_conf * grn;
        cinfo->blu += nearest_conf * blu;

        cinfo->weights += nearest_conf;
        cinfo->count++;

        /* add bit to mesh tags */
        v->old_mesh = (Mesh*)(((size_t) v->old_mesh) | (1 << mesh_index));
    } else {
        /* otherwise, see if we can find a nearest point and use that */
        mesh_to_world_normal(vscan, v->normal, norm);
        result = nearest_on_mesh(mscan, mesh, NULL, wpos, norm,
                                 search_dist, 0.0, &near_info);
        if (result && near_info.on_edge == 0) {

#if 0
            vadd(cinfo->pos, near_info.pos, cinfo->pos);
            v->cinfo->weights += 1;
#endif

            cinfo->pos[X] += near_info.confidence * near_info.pos[X];
            cinfo->pos[Y] += near_info.confidence * near_info.pos[Y];
            cinfo->pos[Z] += near_info.confidence * near_info.pos[Z];

            v->cinfo->weights += near_info.confidence;
            v->cinfo->count++;

            /* add bit to mesh tags */
            v->old_mesh = (Mesh*)(((size_t) v->old_mesh) | (1 << mesh_index));

            switch (near_info.type) {
                case NEAR_VERTEX:
                    red = near_info.v1->red;
                    grn = near_info.v1->grn;
                    blu = near_info.v1->blu;
                    break;
                case NEAR_EDGE:
                    red = near_info.v1->red * near_info.b1 +
                          near_info.v2->red * near_info.b2;
                    grn = near_info.v1->grn * near_info.b1 +
                          near_info.v2->grn * near_info.b2;
                    blu = near_info.v1->blu * near_info.b1 +
                          near_info.v2->blu * near_info.b2;
                    break;
                case NEAR_TRIANGLE:
                    v1 = near_info.tri->verts[0];
                    v2 = near_info.tri->verts[1];
                    v3 = near_info.tri->verts[2];
                    bary_sum = near_info.b1 + near_info.b2 + near_info.b3;
                    red = (near_info.b1 * v1->red +
                           near_info.b2 * v2->red +
                           near_info.b3 * v3->red) / bary_sum;
                    grn = (near_info.b1 * v1->grn +
                           near_info.b2 * v2->grn +
                           near_info.b3 * v3->grn) / bary_sum;
                    blu = (near_info.b1 * v1->blu +
                           near_info.b2 * v2->blu +
                           near_info.b3 * v3->blu) / bary_sum;
                    break;
                default:
                    fprintf(stderr, "intersect_segment_with_mesh: bad switch = %d\n",
                            near_info.type);
                    break;
            }
            v->cinfo->red += near_info.confidence * red;
            v->cinfo->grn += near_info.confidence * grn;
            v->cinfo->blu += near_info.confidence * blu;
        }
    }
}


/******************************************************************************
Find the "average" positions over all scans to the vertices in a given mesh.

Entry:
  scan    - scan containing the given mesh
  level   - level of mesh spacing (0 = every point, 1 = every other, 2 = every
            four, 3 = every eight)
  k_scale - scaling coefficient for search distance
******************************************************************************/

find_average_positions(scan, level, k_scale)
Scan* scan;
int level;
float k_scale;
{
    int i, j;
    Mesh* mesh, *tmesh;
    NearPosition near_info;
    Vector pos, norm;
    int result;
    int zeros;
    int spacing;
    float search_dist;
    Vertex* v;
    Mesh* make_mesh();
    Mesh* make_mesh_raw();
    Mesh* make_mesh_ply();

    mesh = scan->meshes[mesh_level];

    /* determine spacing between mesh elements */
    spacing = level_to_inc(level);

    /* use this to bound the search distance */
    search_dist = CONSENSUS_POSITION_DIST * spacing;

    printf("using search distance of %f\n", search_dist);

    /* examine each mesh to see which points on these meshes are */
    /* close to the vertices of the given mesh */

    for (i = 0; i < nscans; i++) {

        /* don't bother with meshes read from a polygon file */
        if (scans[i]->file_type == POLYFILE)
            continue;

        /* create a mesh of the appropriate level */
        if (scans[i]->file_type == CYFILE)
            tmesh = make_mesh(scans[i], level, search_dist);
        else if (scans[i]->file_type == RAWFILE)
            tmesh = make_mesh_raw(scans[i], level, search_dist);
        else if (scans[i]->file_type == PLYRANGEFILE)
            tmesh = make_mesh_ply(scans[i], level, search_dist);
        else {
            fprintf(stderr, "consensus_surface: bad scan type: %d\n",
                    scans[i]->file_type);
            continue;
        }

        printf("mesh has %d vertices and %d triangles\n", tmesh->nverts, tmesh->ntris);

        /* do neighbor search around vertices */

        for (j = 0; j < mesh->nverts; j++) {
            v = mesh->verts[j];
            mesh_to_world(scan, v->coord, pos);
            mesh_to_world_normal(scan, v->normal, norm);
            result = nearest_on_mesh(scans[i], tmesh, NULL, pos, norm,
                                     search_dist, 0.0, &near_info);
            if (result && near_info.on_edge == 0) {
#if 1
                vadd(v->cinfo->pos, near_info.pos, v->cinfo->pos);
                v->cinfo->weights += 1;
#endif
#if 0
                cinfo->pos[X] += near_info.confidence * near_info.pos[X];
                cinfo->pos[Y] += near_info.confidence * near_info.pos[Y];
                cinfo->pos[Z] += near_info.confidence * near_info.pos[Z];
                v->cinfo->weights += near_info.confidence;
#endif
                v->cinfo->count++;
            }
        }

        /* free the mesh info */
        clear_mesh(tmesh);
    }
}


/******************************************************************************
Print information about sizes of data structures.
******************************************************************************/

#if 0
print_sizes()
{
    int i;
    int tsize, vsize;
    float mverts, mtris;
    Mesh* mesh;

    vsize = sizeof(Vertex);
    tsize = sizeof(Triangle);

    printf("%d bytes per vertex\n", vsize);
    printf("%d bytes per triangle\n", tsize);

    printf("%d bytes per Cvertex\n", sizeof(Cvert));
    printf("%d bytes per Ctriangle\n", sizeof(Ctri));

    for (i = 0; i < nscans; i++) {

        mesh = scans[i]->meshes[mesh_level];
        mverts = vsize * mesh->nverts / 1048576.0;
        mtris  = tsize * mesh->ntris / 1048576.0;

        printf("mesh %s:\n", scans[0]);
        printf("%.2f mbytes for vertices\n", mverts);
        printf("%.2f mbytes for triangles\n", mtris);
    }
}
#endif

/*** set of nearby points to a given position ***/

static Vertex** pts_near = NULL;
static int pts_near_num;
static int pts_near_max;


/******************************************************************************
Return how many vertices are in pts_near.
******************************************************************************/

int count_near_vert()
{
    return (pts_near_num);
}


/******************************************************************************
Return nth vertex in pts_near.
******************************************************************************/

Vertex* found_vert_near_vert(n)
int n;
{
    return (pts_near[n]);
}


/******************************************************************************
Find the collection of nearby vertices to a given position.

Entry:
  mesh     - the mesh containing the (possibly) nearby vertices
  pnt      - position (in mesh's coordinate system) to find nearest vertex to
  norm     - surface normal at pnt
  radius   - distance within which to check

Exit:
  places nearby vertices in "pts_near"
******************************************************************************/

verts_near_pos(mesh, pnt, norm, radius)
Mesh* mesh;
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
    pts_near_num = 0;

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

