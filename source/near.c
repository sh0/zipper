/*

Find nearby points of a mesh to a particular point.

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

static int large_search_flag = 0;


/******************************************************************************
Initialize a uniform spatial subdivision table.  This structure divides
3-space into cubical cells and deposits points into their appropriate
cells.  It uses hashing to make the table a one-dimensional array.

Entry:
  mesh - mesh that contains points to place in table
  size - size of a cell
******************************************************************************/

init_table(mesh, size)
Mesh* mesh;
float size;
{
    int i;
    int index;
    int a, b, c;
    Hash_Table* table;
    float scale;

    /* allocate new hash table */

    table = (Hash_Table*) myalloc(sizeof(Hash_Table));

    if (mesh->nverts < TABLE_SIZE1)
        table->num_entries = TABLE_SIZE1;
    else if (mesh->nverts < TABLE_SIZE2)
        table->num_entries = TABLE_SIZE2;
    else
        table->num_entries = TABLE_SIZE3;

    table->verts = (Vertex**) myalloc(sizeof(Vertex*) * table->num_entries);
    mesh->table = table;

    /* set all table elements to NULL */
    for (i = 0; i < table->num_entries; i++)
        table->verts[i] = NULL;

    /* place each point in table */

    scale = 1 / size;
    table->scale = scale;

    for (i = 0; i < mesh->nverts; i++) {
        a = floor(mesh->verts[i]->coord[X] * scale);
        b = floor(mesh->verts[i]->coord[Y] * scale);
        c = floor(mesh->verts[i]->coord[Z] * scale);
        index = (a * PR1 + b * PR2 + c) % table->num_entries;
        if (index < 0)
            index += table->num_entries;
        mesh->verts[i]->next = table->verts[index];
        table->verts[index] = mesh->verts[i];
    }
}


/******************************************************************************
Add a vertex to it's hash table.

Entry:
  vert - vertex to add
  mesh - mesh to add to
******************************************************************************/

add_to_hash(vert, mesh)
Vertex* vert;
Mesh* mesh;
{
    int index;
    int a, b, c;
    Hash_Table* table;
    float scale;

    table = mesh->table;
    scale = table->scale;

    a = floor(vert->coord[X] * scale);
    b = floor(vert->coord[Y] * scale);
    c = floor(vert->coord[Z] * scale);
    index = (a * PR1 + b * PR2 + c) % table->num_entries;
    if (index < 0)
        index += table->num_entries;
    vert->next = table->verts[index];
    table->verts[index] = vert;
}


/******************************************************************************
Remove a vertex from the hash table.

Entry:
  vert - vertex to remove from hash table
  mesh - mesh that vertex belongs to
******************************************************************************/

remove_from_hash(vert, mesh)
Vertex* vert;
Mesh* mesh;
{
    int index;
    int a, b, c;
    Hash_Table* table;
    float scale;
    int found;
    Vertex* ptr;
    static int been_here = 0;

    table = mesh->table;
    scale = table->scale;

    /* determine which hash cell vertex is in */
    a = floor(vert->coord[X] * scale);
    b = floor(vert->coord[Y] * scale);
    c = floor(vert->coord[Z] * scale);
    index = (a * PR1 + b * PR2 + c) % table->num_entries;
    if (index < 0)
        index += table->num_entries;

    /* see if this vertex is the first one in the hash cell */
    if (table->verts[index] == vert) {
        table->verts[index] = vert->next;
        return;
    }

    /* look for pointer to this vertex */
    found = 0;
    for (ptr = table->verts[index]; ptr != NULL; ptr = ptr->next) {
        if (ptr->next == vert) {
            found = 1;
            break;
        }
    }

    /* error check */
    if (!found) {
        fprintf(stderr, "remove_from_hash: can't find vertex\n");
        printf("index = %d, x y z = %f %f %f\n", vert->index,
               vert->coord[X], vert->coord[Y], vert->coord[Z]);
        printf("hash index: %d\n", index);
        return;
        /*
        exit (-1);
        */
    }

    /* remove the vertex from the hash table */
    ptr->next = vert->next;
}


/******************************************************************************
Find the nearest vertex in a mesh to a given position.

Entry:
  mesh     - the mesh
  not_mesh - mesh to reject matches from (old_mesh field of Vertex)
  pnt      - position (in mesh's coordinate system) to find nearest vertex to
  norm     - surface normal at pnt
  min_dot  - minimum allowed dot product between given normal and match point

Exit:
  returns pointer to the nearest vertex
******************************************************************************/

Vertex* find_nearest(mesh, not_mesh, pnt, norm, min_dot)
Mesh* mesh;
Mesh* not_mesh;
Vector pnt;
Vector norm;
float min_dot;
{
    int a, b, c;
    int aa, bb, cc;
    int index;
    Hash_Table* table = mesh->table;
    Vertex* ptr;
    Vertex* min_ptr = NULL;
    float dx, dy, dz;
    float dist;
    float min_dist = 1e20;
    float dot;

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

                    /* distance (squared) to this point */
                    dx = ptr->coord[X] - pnt[X];
                    dy = ptr->coord[Y] - pnt[Y];
                    dz = ptr->coord[Z] - pnt[Z];
                    dist = dx * dx + dy * dy + dz * dz;

                    /* maybe we've found new closest point */
                    if (dist < min_dist) {

                        /* make sure the surface normals are roughly in the same direction */
                        dot = vdot(norm, ptr->normal);
                        if (dot < min_dot)
                            continue;

                        min_dist = dist;
                        min_ptr = ptr;
                    }
                }
            }

    /* return nearest point */
    return (min_ptr);
}


/******************************************************************************
Hacked version of find_nearest() to search beyond the nearest hash cells.
Find the nearest vertex in a mesh to a given position.

Entry:
  mesh     - the mesh
  not_mesh - mesh to reject matches from (old_mesh field of Vertex)
  pnt      - position (in mesh's coordinate system) to find nearest vertex to
  norm     - surface normal at pnt
  min_dot  - minimum allowed dot product between given normal and match point

Exit:
  returns pointer to the nearest vertex
******************************************************************************/

Vertex* large_find_nearest(mesh, not_mesh, pnt, norm, min_dot)
Mesh* mesh;
Mesh* not_mesh;
Vector pnt;
Vector norm;
float min_dot;
{
    int a, b, c;
    int aa, bb, cc;
    int index;
    Hash_Table* table = mesh->table;
    Vertex* ptr;
    Vertex* min_ptr = NULL;
    float dx, dy, dz;
    float dist;
    float min_dist = 1e20;
    float dot;
    int width = 2;

    /* determine which cell the position lies within */
    aa = floor(table->scale * pnt[X]);
    bb = floor(table->scale * pnt[Y]);
    cc = floor(table->scale * pnt[Z]);

    /* look at nine cells, centered at cell containing location */

    for (a = aa - width; a <= aa + width; a++)
        for (b = bb - width; b <= bb + width; b++)
            for (c = cc - width; c <= cc + width; c++) {

                /* compute position in hash table */
                index = (a * PR1 + b * PR2 + c) % table->num_entries;
                if (index < 0)
                    index += table->num_entries;

                /* examine all points hashed to this cell */
                for (ptr = table->verts[index]; ptr != NULL; ptr = ptr->next) {

                    /* don't examine point if it's old mesh is not_mesh */
                    if (ptr->old_mesh == not_mesh)
                        continue;

                    /* distance (squared) to this point */
                    dx = ptr->coord[X] - pnt[X];
                    dy = ptr->coord[Y] - pnt[Y];
                    dz = ptr->coord[Z] - pnt[Z];
                    dist = dx * dx + dy * dy + dz * dz;

                    /* maybe we've found new closest point */
                    if (dist < min_dist) {

                        /* make sure the surface normals are roughly in the same direction */
                        dot = vdot(norm, ptr->normal);
                        if (dot < min_dot)
                            continue;

                        min_dist = dist;
                        min_ptr = ptr;
                    }
                }
            }

    /* return nearest point */
    return (min_ptr);
}


/******************************************************************************
Find the nearest vertex in a mesh to a given position.

Entry:
  mesh     - the mesh
  not_mesh - mesh to reject matches from (old_mesh field of Vertex)
  pnt      - position (in mesh's coordinate system) to find nearest vertex to
  norm     - surface normal at pnt
  max      - maximum acceptable distance
  min_dot  - minimum allowed dot product between given normal and match point

Exit:
  returns pointer to the nearest vertex
******************************************************************************/

Vertex* new_find_nearest(mesh, not_mesh, pnt, norm, max, min_dot)
Mesh* mesh;
Mesh* not_mesh;
Vector pnt;
Vector norm;
float max, min_dot;
{
    int a, b, c;
    int aa, bb, cc;
    int index, width;
    Hash_Table* table = mesh->table;
    Vertex* ptr;
    Vertex* min_ptr = NULL;
    float dx, dy, dz;
    float dist;
    float min_dist = 1e20;
    float dot;

    /* determine which cell the position lies within */
    aa = floor(table->scale * pnt[X]);
    bb = floor(table->scale * pnt[Y]);
    cc = floor(table->scale * pnt[Z]);



    width = ceil(max * table->scale - 1e-4);

    for (a = aa - width; a <= aa + width; a++)
        for (b = bb - width; b <= bb + width; b++)
            for (c = cc - width; c <= cc + width; c++) {

                /* compute position in hash table */
                index = (a * PR1 + b * PR2 + c) % table->num_entries;
                if (index < 0)
                    index += table->num_entries;

                /* examine all points hashed to this cell */
                for (ptr = table->verts[index]; ptr != NULL; ptr = ptr->next) {

                    /* don't examine point if it's old mesh is not_mesh */
                    if (ptr->old_mesh == not_mesh)
                        continue;

                    /* distance (squared) to this point */
                    dx = ptr->coord[X] - pnt[X];
                    dy = ptr->coord[Y] - pnt[Y];
                    dz = ptr->coord[Z] - pnt[Z];
                    dist = dx * dx + dy * dy + dz * dz;

                    /* maybe we've found new closest point */
                    if (dist < min_dist) {

                        /* make sure the surface normals are roughly in the same direction */
                        dot = vdot(norm, ptr->normal);
                        if (dot < min_dot)
                            continue;

                        min_dist = dist;
                        min_ptr = ptr;
                    }
                }
            }

    /* return nearest point */
    return (min_ptr);
}


/******************************************************************************
Kluge to search farther than nearest hash cellsd.
******************************************************************************/

use_large_search(boolean)
int boolean;
{
    large_search_flag = boolean;
}


/******************************************************************************
Find the nearest location on a mesh to a given position.  This nearest
position can be at a vertex, on an edge or within a triangle of the mesh.

Entry:
  sc       - the scan of the mesh
  mesh     - the mesh
  not_mesh - mesh to reject matches from (old_mesh field of Vertex)
  pos      - position (in global coordinates) to find nearest point to
  norm     - surface normal at pos
  max      - maximum acceptable distance
  min_dot  - minimum allowed dot product between given normal and match point

Exit:
  near_info - information about nearest position on mesh, including 3D position
          in global coordinates
  returns 1 if it found a near point, 0 if it couldn't find any near point
******************************************************************************/

int nearest_on_mesh(sc, mesh, not_mesh, pos, norm, max, min_dot, near_info)
Scan* sc;
Mesh* mesh;
Mesh* not_mesh;
Vector pos;
Vector norm;
float max;
float min_dot;
NearPosition* near_info;
{
    Vector v;
    Vector tnorm;
    Vector near_pos;
    Vector diff;
    Vertex* near, *near2;
    float min_dist;
    float edge_min_dist, tri_min_dist;
    float nearest_on_tris();
    float nearest_on_edges();
    int on_edge1, on_edge2;
    int near_type;
    Triangle* near_tri;
    float conf;
    float ival;
    Vector barycentric;

    /* transform position into the meshes coordinate space */
    world_to_mesh(sc, pos, v);
    world_to_mesh_normal(sc, norm, tnorm);

    /* look for nearby vertex of mesh */

#if 0
    if (large_search_flag)
        near = large_find_nearest(mesh, not_mesh, v, tnorm, min_dot);
    else
        near = find_nearest(mesh, not_mesh, v, tnorm, min_dot);
#else
    near = new_find_nearest(mesh, not_mesh, v, tnorm, max, min_dot);
#endif

    /* return now if there was no nearby vertex or if vertex has no triangles */
    if ((near == NULL) || (near->ntris == 0)) {
        return (0);
    }

    /* so far, this vertex is the closest thing to the given position */
    vsub(v, near->coord, diff);
    min_dist = vlen(diff);
    vcopy(near->coord, near_pos);
    on_edge1 = near->on_edge;  /* is this vertex on the edge of the mesh? */
    near_type = NEAR_VERTEX;   /* so far, nearest position is this vertex */
    near_info->v1 = near;      /* remember which vertex was nearest */
    near_info->confidence = near->confidence;

    /* see if position is nearer to any edges of the vertex */

    /*
      new_min_dist = nearest_on_edges (v, near, min_dist, near_pos,
                       &on_edge2, &near2, &conf, &ival);
    */

    edge_min_dist = nearest_on_edges(v, near, 1e20, near_pos,
                                     &on_edge2, &near2, &conf, &ival);
    near_info->v2 = near2;    /* remember other vertex involved */

    if (edge_min_dist < min_dist) {
        min_dist = edge_min_dist;
        near_type = NEAR_EDGE;    /* nearest position is on an edge */
        near_info->v2 = near2;    /* remember other vertex involved */
        near_info->confidence = conf;
        near_info->b1 = 1 - ival;
        near_info->b2 = ival;
    }

    /* get nearest position on triangles radiating from "near" vertex, */
    /* or leave current best position alone */

    tri_min_dist = nearest_on_tris(v, near, 1e20, near_pos, &near_tri, &conf,
                                   barycentric);
    near_info->tri = near_tri;

    if (tri_min_dist < min_dist) {
        min_dist = tri_min_dist;
        near_type = NEAR_TRIANGLE;    /* nearest position is on a triangle */
        near_info->tri = near_tri;
        near_info->confidence = conf;
        near_info->b1 = barycentric[X];
        near_info->b2 = barycentric[Y];
        near_info->b3 = barycentric[Z];
    }

    /* see if we found a near enough position that is not on the */
    /* edge of the mesh */

    near_info->on_edge = 0;  /* assume we're not on edge */
    near_info->dist = min_dist;
    near_info->type = near_type;

    if (min_dist < max) {
        if (near_type == NEAR_VERTEX) {
            mesh_to_world(sc, near_pos, near_info->pos);
            near_info->on_edge = on_edge1;
            return (1);
        } else if (near_type == NEAR_EDGE) {
            mesh_to_world(sc, near_pos, near_info->pos);
            near_info->on_edge = (on_edge1 && on_edge2);
            return (1);
        } else if (near_type == NEAR_TRIANGLE) {
            mesh_to_world(sc, near_pos, near_info->pos);
            return (1);
        }
    }

    /* if we get here, we didn't find a satisfactory near position */
    /* (too far away or on the mesh edge) */
    return (0);
}


/******************************************************************************
Find the nearest location on a group of edges of a vertex to a given
position.

Entry:
  pos  - position (in local coordinates) to find nearest point to
  near - the vertex whose edges will be checked
  max  - maximum acceptable distance

Exit:
  near_pos - nearest position on edges (in local coordinates)
  on_edge  - whether other end of the found edge is on the edge of the mesh
  near2    - other vertex of edge
  conf     - confidence in position on edge
  ival     - interpolation value between vertices
  returns distance to nearby position if it found a near point, max if it
  couldn't find any near point
******************************************************************************/

float nearest_on_edges(pos, near, max, near_pos, on_edge, near2, conf, ival)
Vector pos;
Vertex* near;
float max;
Vector near_pos;
int* on_edge;
Vertex** near2;
float* conf;
float* ival;
{
    int i;
    Vector p, q;
    Vector dir;
    float len;
    float t;

    /* one end of all edges comes from the given vertex */
    vcopy(near->coord, p);

    /* check each edge (p,q) that radiates from the given vertex */

    for (i = 0; i < near->nverts; i++) {

        /* other end of the current edge */
        vcopy(near->verts[i]->coord, q);

        /* compute unit length vector pointing from p towards q */
        vsub(q, p, dir);
        len = vnorm(dir);

        /* see if nearest point on line (p,q) is between p and q */

        t = vdot(pos, dir) - vdot(p, dir);

        if (t > 0 && t < len) {
            Vector close, diff;

            /* determine nearest point on line */
            close[X] = p[X] + t * dir[X];
            close[Y] = p[Y] + t * dir[Y];
            close[Z] = p[Z] + t * dir[Z];

            /* see how far this point is from "pos" */
            vsub(pos, close, diff);
            len = vlen(diff);

            /* if it is close enough, we have a new nearest point */
            if (len < max) {
                max = len;
                vcopy(close, near_pos);
                *on_edge = near->verts[i]->on_edge;
                *near2 = near->verts[i];
                *conf = near->confidence +
                        t * (near->verts[i]->confidence - near->confidence);
                *ival = t;
            }
        }
    }

    /* return nearest distance, or original maximum if we didn't find */
    /* any positions near enough */
    return (max);
}


/******************************************************************************
Find the nearest location on a group of triangles of a vertex to a given
position.

Entry:
  pos  - position (in local coordinates) to find nearest point to
  near - the vertex whose triangles will be checked
  max  - maximum acceptable distance

Exit:
  near_pos    - nearest position on triangles (in local coordinates)
  near_tri    - the nearest triangle
  conf        - error estimate at near position
  barycentric - barycentric coordinates of near point
  returns distance to nearby position if it found a near point, max if it
  couldn't find any near point
******************************************************************************/

float nearest_on_tris(pos, near, max, near_pos, near_tri, conf, barycentric)
Vector pos;
Vertex* near;
float max;
Vector near_pos;
Triangle** near_tri;
float* conf;
Vector barycentric;
{
    int i;
    Triangle* tri;
    float r1, r2, r3;
    float len;
    Vector on;

    /* examine each triangle of the vertex to see if there is a close */
    /* position on the triangle */

    for (i = 0; i < near->ntris; i++) {

        tri = near->tris[i];

        /* find nearest point in the triangle's plane */
        len = pos[X] * tri->aa + pos[Y] * tri->bb + pos[Z] * tri->cc + tri->dd;
        on[X] = pos[X] - len * tri->aa;
        on[Y] = pos[Y] - len * tri->bb;
        on[Z] = pos[Z] - len * tri->cc;

        /* see if this point is in the triangle */
        r1 = on[X] * tri->a[0] + on[Y] * tri->b[0] + on[Z] * tri->c[0] + tri->d[0];
        if (r1 < 0) continue;
        r2 = on[X] * tri->a[1] + on[Y] * tri->b[1] + on[Z] * tri->c[1] + tri->d[1];
        if (r2 < 0) continue;
        r3 = on[X] * tri->a[2] + on[Y] * tri->b[2] + on[Z] * tri->c[2] + tri->d[2];
        if (r3 < 0) continue;

        /* if we get here, the point is in the triangle, */
        /* so compare it's distance with the current upper bound on distance */
        len = fabs(len);
        if (len < max) {
            max = len;
            vcopy(on, near_pos);
            *near_tri = tri;
            *conf = (r2 * tri->verts[0]->confidence +
                     r3 * tri->verts[1]->confidence +
                     r1 * tri->verts[2]->confidence) / (r1 + r2 + r3);
            barycentric[X] = r2;
            barycentric[Y] = r3;
            barycentric[Z] = r1;
        }
    }

    /* return nearest distance, or original maximum if we didn't find */
    /* any positions near enough */
    return (max);
}


int
is_near_edge(near)
NearPosition* near;
{
    int i, is_near;
    Vertex* v;

    is_near = 0;
    if (near->type == NEAR_VERTEX) {
        v = near->v1;
        is_near = is_near || v->on_edge || is_vertex_near_edge(v);
    } else if (near->type == NEAR_EDGE) {
        v = near->v1;
        is_near = is_near || v->on_edge || is_vertex_near_edge(v);
        v = near->v2;
        is_near = is_near || v->on_edge || is_vertex_near_edge(v);
    } else if (near->type == NEAR_TRIANGLE) {
        v = near->tri->verts[0];
        is_near = is_near || v->on_edge || is_vertex_near_edge(v);
        v = near->tri->verts[1];
        is_near = is_near || v->on_edge || is_vertex_near_edge(v);
        v = near->tri->verts[2];
        is_near = is_near || v->on_edge || is_vertex_near_edge(v);
    }

    return is_near;
}


int
is_vertex_near_edge(v)
Vertex* v;
{
    int i, is_near;

    is_near = 0;

    for (i = 0; i < v->nverts; i++) {
        is_near = is_near || v->verts[i]->on_edge;
    }

    return is_near;
}

