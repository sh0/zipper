/*

Find a matching position for meshes.

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
#include <string.h>
#include "zipper.h"
#include "matrix.h"
#include "zipglobal.h"

/* list of how scans are to be matched */

#define  MATCH_MAX  20
Scan* match_from[MATCH_MAX];
Scan* match_to[MATCH_MAX];
Scan* match_drag[MATCH_MAX][SCAN_MAX];
int num_drag[MATCH_MAX];
int num_matches = 0;

/* list of matches between vertices on one mesh and positions on another mesh */

Match** pos_matches = NULL;
int global_num_matches = 0;
int global_max_matches = 1000;

int parallel_procs_max = 1;
int processes_forked = 0;

static float ALIGN_NEAR_DIST_FACTOR;
static float ALIGN_NEAR_DIST;
static float ALIGN_NEAR_COS;


update_align_resolution()
{
    ALIGN_NEAR_DIST = ZIPPER_RESOLUTION * ALIGN_NEAR_DIST_FACTOR;
}


set_align_near_dist_factor(factor)
float factor;
{
    ALIGN_NEAR_DIST_FACTOR = factor;
    ALIGN_NEAR_DIST = ZIPPER_RESOLUTION * ALIGN_NEAR_DIST_FACTOR;
}


float
get_align_near_dist_factor()
{
    return ALIGN_NEAR_DIST_FACTOR;
}


set_align_near_cos(cosine)
float cosine;
{
    ALIGN_NEAR_COS = cosine;
}


float
get_align_near_cos()
{
    return ALIGN_NEAR_COS;
}



/******************************************************************************
Match together multiple meshes.
******************************************************************************/

match_proc()
{
    int i;

    if (num_matches == 0) {
        new_iterative_match(&scans[0], &scans[1], NULL, 1, 1, 0);
    } else
        for (i = 0; i < num_matches; i++) {
            new_iterative_match(&match_from[i], &match_to[i], match_drag[i],
                                1, 1, num_drag[i]);
        }
}


/******************************************************************************
Try to match one group of meshes to another by iteratively moving them
small distances.  Drag other meshes along for the ride.

Entry:
  move_list - list meshes to move
  match_to  - list of meshes to match to
  drag_with - list of other scans to move along with move_list
  nmove     - number of scans to move
  nmatch    - number of scans to match
  ndrag     - number in drag list
******************************************************************************/

new_iterative_match(move_list, match_to, drag_with, nmove, nmatch, ndrag)
Scan** move_list;
Scan** match_to;
Scan** drag_with;
int nmove;
int nmatch;
int ndrag;
{
    int i, j;
    float dist_sum;
    float old_sum;
    Vector move_dir;
    Vector old_dir;
    float fact = 1.0;
    int count;
    Scan* sc1, *sc2;

    sc1 = move_list[0];
    sc2 = match_to[0];

    old_sum = 1e20;
    old_dir[X] = old_dir[Y] = old_dir[Z] = 0;

    for (i = 0; i < 10; i++) {

        /* find out which way to move mesh */
        one_match(sc1, sc2, move_dir, &count, &dist_sum);

        /* move the mesh */
        if (dist_sum > 0) {

            move_dir[X] /= (float) count;
            move_dir[Y] /= (float) count;
            move_dir[Z] /= (float) count;

            sc1->xtrans += move_dir[X] * fact;
            sc1->ytrans += move_dir[Y] * fact;
            sc1->ztrans += move_dir[Z] * fact;

            /* drag other mesh along */
            for (j = 0; j < ndrag; j++) {
                drag_with[j]->xtrans += move_dir[X] * fact;
                drag_with[j]->ytrans += move_dir[Y] * fact;
                drag_with[j]->ztrans += move_dir[Z] * fact;
            }
        }

#if 0
        /* change size of move if we switched direction of motion */
        dot = vdot(move_dir, old_dir);
        if (dot < 0) {
            fact *= 0.75;
            /*
            printf ("iter fact dot: %d %f %g\n", i, fact, dot);
            */
        }
#endif

        vcopy(move_dir, old_dir);
    }
}


/******************************************************************************
Try to match two meshes by iteratively moving them small distances.

Entry:
  sc1       - scan containing mesh to move
  sc2       - scan to match to
  drag_with - list of other scans to move along with sc1
  ndrag     - number in drag list
******************************************************************************/

iterative_match(sc1, sc2, drag_with, ndrag)
Scan* sc1, *sc2;
Scan* drag_with[];
int ndrag;
{
    int i, j;
    float dist_sum;
    float old_sum;
    Vector move_dir;
    Vector old_dir;
    float fact = 0.01;
    float dot;
    int count;

    old_sum = 1e20;
    old_dir[X] = old_dir[Y] = old_dir[Z] = 0;

    for (i = 0; i < 20; i++) {

        /* find out which way to move mesh */
        one_match(sc1, sc2, move_dir, &count, &dist_sum);

        /* move the mesh */
        if (dist_sum > 0) {

            move_dir[X] /= dist_sum;
            move_dir[Y] /= dist_sum;
            move_dir[Z] /= dist_sum;

            sc1->xtrans += move_dir[X] * fact;
            sc1->ytrans += move_dir[Y] * fact;
            sc1->ztrans += move_dir[Z] * fact;

            /* drag other mesh along */
            for (j = 0; j < ndrag; j++) {
                drag_with[j]->xtrans += move_dir[X] * fact;
                drag_with[j]->ytrans += move_dir[Y] * fact;
                drag_with[j]->ztrans += move_dir[Z] * fact;
            }
        }

        /* change size of move if we switched direction of motion */
        dot = vdot(move_dir, old_dir);
        if (dot < 0) {
            fact *= 0.75;
            /*
            printf ("iter fact dot: %d %f %f\n", i, fact, dot);
            */
        }

        vcopy(move_dir, old_dir);
    }
}


/******************************************************************************
Compare the lengths of two match entries.

Entry:
  m1,m2 - the match entries

Exit:
  returns 0 if lengths are same, +1 or -1 if they are different (for sorting)
******************************************************************************/

int compare_matches(m1, m2)
Match* m1, *m2;
{
    if (m1->len < m2->len)
        return (1);
    else if (m1->len > m2->len)
        return (-1);
    else
        return (0);
}


/******************************************************************************
Make a list of matches between vertices of one mesh and positions on a
second mesh.

Entry:
  sc1            - scan containing first mesh
  sc2            - scan of second mesh
  edges1_okay - if it is okay to return matches that are on an edge of mesh 1
  edges2_okay - return matches on mesh 2 edges?
******************************************************************************/

create_match_list(sc1, sc2, edges1_okay, edges2_okay)
Scan* sc1, *sc2;
int edges1_okay;
int edges2_okay;
{
    int i;
    void single_create_match_list();

    /* initialize the match list */

    if (pos_matches == NULL) {
        pos_matches = (Match**) malloc(sizeof(Match*) * global_max_matches);
    } else {
        for (i = 0; i < global_num_matches; i++)
            free(pos_matches[i]);
    }

    global_num_matches = 0;

    /* fork off several processes to make the list */

    //#define  PARALLEL  1

#ifdef PARALLEL
    m_set_procs(parallel_procs_max);
    m_fork(single_create_match_list, sc1, sc2, edges1_okay, edges2_okay);
    processes_forked = 1;
#else
    single_create_match_list(sc1, sc2, edges1_okay, edges2_okay);
#endif

    /*
    printf ("global_num_matches = %d\n", global_num_matches);
    */

#if 0
    /* sort the match list in order of descending length */
    qsort(pos_matches, global_num_matches, sizeof(Match*), compare_matches);
#endif

#if 0
    /* print first few match lengths */
    for (i = 0; i < 20; i++) {
        if (i % 4 == 0)
            printf("\n");
        printf("%2d) %.7f ", i, pos_matches[i]->len);
    }
    printf("\n");
#endif
}


/******************************************************************************
Make a list of matches between vertices of one mesh and positions on a
second mesh.  Do this only for a sub-range of the points in the mesh.

Entry:
  sc1         - scan containing first mesh
  sc2         - scan of second mesh
  edges1_okay - if it is okay to return matches that are on an edge of mesh 1
  edges2_okay - return matches on mesh 2 edges?
******************************************************************************/

void single_create_match_list(sc1, sc2, edges1_okay, edges2_okay)
Scan* sc1, *sc2;
int edges1_okay;
int edges2_okay;
{
    int i;
    int start, end;
    int result;
    Mesh* m1, *m2;
    Vector pos, norm;
    Vector near_pos;
    Vector diff;
    float len;
    NearPosition near_info;
    int inc;
    int processes;
    int myid;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    if (m1 == NULL || m2 == NULL) {
        fprintf(stderr, "single_create_match_list: null mesh\n");
    }

    /* find the nearest points on mesh #2 to points on mesh #1 */

#ifdef PARALLEL
    myid = m_get_myid();
    processes = m_get_numprocs();
    start = m1->nverts * (myid / (float) processes);
    end   = m1->nverts * ((myid + 1) / (float) processes) - 1;
#else
    start = 0;
    end   = m1->nverts - 1;
#endif

#if 0
    printf("id start end num: %d %d %d %d\n", myid, start, end, m1->nverts);
#endif

    for (i = start; i <= end; i++) {

        /* maybe don't bother to find nearby points to vertices that are on */
        /* the edge of the mesh */
        if (!edges1_okay && m1->verts[i]->on_edge)
            continue;

        /* convert a vertex position on the first mesh to global coordinates */
        mesh_to_world(sc1, m1->verts[i]->coord, pos);
        mesh_to_world_normal(sc1, m1->verts[i]->normal, norm);

        /* find nearest position on the second mesh */
        inc = level_to_inc(mesh_level);
        result = nearest_on_mesh(sc2, m2, NULL, pos, norm,
                                 ALIGN_NEAR_DIST * inc, ALIGN_NEAR_COS,
                                 &near_info);
        vcopy(near_info.pos, near_pos);

        /* go to next vertex if there is no position close enough */
        if (result == 0)
            continue;

        /* go on if we don't want edges and we're on an edge */
        if (!edges2_okay && near_info.on_edge)
            continue;

#ifdef PARALLEL
        m_lock();
#endif

        /* make sure there is room in the match list */
        if (global_num_matches >= global_max_matches) {
            global_max_matches *= 2;
            pos_matches = (Match**)
                          realloc(pos_matches, sizeof(Match*) * global_max_matches);
        }
        pos_matches[global_num_matches] = (Match*) malloc(sizeof(Match));

        /* save the information about this match in the match list */

        vsub(near_pos, pos, diff);
        len = vlen(diff);
        pos_matches[global_num_matches]->len = len;
        vcopy(diff, pos_matches[global_num_matches]->dir);
        vcopy(pos, pos_matches[global_num_matches]->pos);
        pos_matches[global_num_matches]->vert = m1->verts[i];
        pos_matches[global_num_matches]->confidence = near_info.confidence *
                m1->verts[i]->confidence;
        global_num_matches++;

#ifdef PARALLEL
        m_unlock();
#endif

    }
}


/******************************************************************************
Try to match two meshes.

Entry:
  sc1,sc2 - scans containing the two meshes to match

Exit:
  move_dir    - sum of vectors stretching from one mesh to the other
  n_matches - number of matches between meshes
  dist_sum    - sum of distances between corresponding points on mesh
******************************************************************************/

one_match(sc1, sc2, move_dir, n_matches, dist_sum)
Scan* sc1, *sc2;
Vector move_dir;
int* n_matches;
float* dist_sum;
{
    int i;
    int match_count = 0;

    /* find the nearest points on mesh #2 to points on mesh #1 */

    create_match_list(sc1, sc2, 1, 0);

    move_dir[X] = move_dir[Y] = move_dir[Z] = 0;
    *dist_sum = 0;

    for (i = 0; i < global_num_matches; i++) {

        /* combine the direction towards the nearest position with */
        /* previous directions */

        *dist_sum += pos_matches[i]->len;
        vadd(move_dir, pos_matches[i]->dir, move_dir);
        match_count++;
    }

    /* return the number of matches between meshes */
    *n_matches = match_count;
}


/******************************************************************************
Returns pointer to scan, given a name.

Entry:
  name - name of scan

Exit:
  returns pointer to scan, or NULL if it can't be found
******************************************************************************/

Scan* find_scan(name)
char* name;
{
    int i;

    for (i = 0; i < nscans; i++)
        if (strcmp(name, scans[i]->name) == 0)
            return (scans[i]);

    return (NULL);
}


/******************************************************************************
Returns pointer to scan, given a name.

Entry:
  name - name of scan

Exit:
  returns pointer to scan, or NULL if it can't be found
******************************************************************************/

Scan* rename_scan(old_name, new_name)
char* old_name, *new_name;
{
    int i;

    for (i = 0; i < nscans; i++) {
        if (strcmp(old_name, scans[i]->name) == 0) {
            strcpy(scans[i]->name, new_name);
            return (scans[i]);
        }
    }

    return (NULL);
}


/******************************************************************************
Reads file telling how to match up different scans.

Entry:
  filename - name of file to read from
******************************************************************************/

read_matches(filename)
char* filename;
{
    int i, j, k;
    FILE* fp;
    char* s;
    char str[200];
    char name[200];
    char name1[200];
    char name2[200];
    int num;
    int done;

    if (strlen(filename) < 6 ||
        strcmp(filename + strlen(filename) - 6, ".match") != 0)
        strcat(filename, ".match");
    printf("reading match configuration from %s\n", filename);

    /* open file for reading */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "read_matches: can't open '%s'\n", filename);
        return;
    }

    /* read info about each scan and read in the scans */

    fscanf(fp, "matches: %d\n", &num);
    num_matches = 0;

    for (i = 0; i < num; i++) {

        /* fetch line out of file */
        fgets(str, 200, fp);

        /* get name of "from" scan */
        for (j = 0, s = str; (*s != ' ') && (*s != '\012'); s++)
            name1[j++] = *s;
        name1[j] = '\0';
        if (*s != '\012')
            s++;
        match_from[num_matches] = find_scan(name1);

        /* get name of "to" scan */
        for (j = 0; (*s != ' ') && (*s != '\012'); s++)
            name2[j++] = *s;
        name2[j] = '\0';
        if (*s == '\012')
            done = 1;
        else {
            done = 0;
            s++;
        }
        match_to[num_matches] = find_scan(name2);

        /* get names of scans to be dragged along with "from" scan */
        k = 0;
        num_drag[i] = 0;
        while (!done) {
            for (j = 0; (*s != ' ') && (*s != '\012'); s++)
                name[j++] = *s;
            name[j] = '\0';
            if (*s == '\012')
                done = 1;
            else
                s++;
            match_drag[i][k] = find_scan(name);
            if (match_drag[i][k] == NULL) {
                printf("bad scan name: '%s'\n", name);
                num_matches = 0;
                return;
            }
            num_drag[i]++;
        }

        if ((match_from[num_matches] == NULL) || (match_to[num_matches] == NULL)) {
            printf("bad scan names: '%s' and '%s'\n", name1, name2);
            num_matches = 0;
            break;
        } else {
            num_matches++;
        }
    }

    fclose(fp);
}


/******************************************************************************
Cause all meshes to conform to those meshes they are matched to.
******************************************************************************/

void all_conform()
{
    int i;

    if (num_matches == 0)
        mesh_conform(scans[0], scans[1]);
    else
        for (i = 0; i < num_matches; i++)
            mesh_conform(match_from[i], match_to[i]);
}


/******************************************************************************
Cause the vertices of one mesh to lie on the surface of another.

Entry:
  sc1 - one of the meshes
  sc2 - the other mesh
******************************************************************************/

mesh_conform(sc1, sc2)
Scan* sc1, *sc2;
{
    int i;
    Mesh* m1, *m2;
    Vertex* vert;

    m1 = sc1->meshes[mesh_level];
    m2 = sc2->meshes[mesh_level];

    /* create list of matches between vertices of sc1 to positions on sc2 */
    create_match_list(sc1, sc2, 1, 0);

    /* move all those vertices in the match list */
    for (i = 0; i < global_num_matches; i++) {
        vert = pos_matches[i]->vert;
        vadd(vert->coord, pos_matches[i]->dir, vert->coord);
    }

    /* re-compute surface normals for the mesh */
    find_vertex_normals(m1);
}

