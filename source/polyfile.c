/*

Read and write ascii polygon files.

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

Scan* new_scan(char*, int);

/******************************************************************************
Write the mesh of a scan to a file.

Entry:
  scan     - which scan to write
  filename - name of file to write to
******************************************************************************/
scan_to_file(scan, filename)
Scan* scan;
char* filename;
{
    int i;
    Mesh* mesh;
    Vertex* v;
    FILE* fp;
    Triangle* tri;

    if (strlen(filename) < 5 ||
        strcmp(filename + strlen(filename) - 5, ".poly") != 0)
        strcat(filename, ".poly");

    fp = fopen(filename, "w");

    if (fp == NULL) {
        fprintf(stderr, "Can't open file '%s'\n", filename);
        return;
    }

    printf("Writing to '%s'\n", filename);

    /* write header */

    mesh = scan->meshes[mesh_level];

    fprintf(fp, "vertices: %d\n", mesh->nverts);
    fprintf(fp, "faces: %d\n", mesh->ntris);

    /* write out the vertices */

    for (i = 0; i < mesh->nverts; i++) {
        v = mesh->verts[i];
        fprintf(fp, "v %g %g %g  %g %g %g\n",
                v->coord[X], v->coord[Y], v->coord[Z],
                v->normal[X], v->normal[Y], v->normal[Z]);
    }

    /* write out the triangles */

    for (i = 0; i < mesh->ntris; i++) {

        /* write out a triangle */

        tri = mesh->tris[i];
        fprintf(fp, "f %d %d %d\n",
                tri->verts[0]->index + 1,
                tri->verts[1]->index + 1,
                tri->verts[2]->index + 1);
    }

    fclose(fp);
}


/******************************************************************************
Read a shared-vertex format polygon file into a scan.

Entry:
  filename - name of file to read in

Exit:
  returns 0 if file was read okay, 1 if not
******************************************************************************/

int file_to_scan(filename)
char* filename;
{
    int i, j;
    FILE* fp;
    Scan* sc;
    Mesh* mesh;
    char str[200];
    char name[80];
    Vector pos;
    float nx, ny, nz;
    int in1, in2, in3;
    int result;
    Vertex* v1, *v2, *v3;
    int nverts, ntris;
    int inc;

    strcpy(name, filename);

    if (strlen(name) < 5 ||
        strcmp(name + strlen(name) - 5, ".poly") != 0)
        strcat(name, ".poly");

    fp = fopen(name, "r");

    if (fp == NULL) {
        printf("Couldn't open file '%s'\n", name);
        return (1);
    }

    printf("Reading polygons from '%s'\n", name);

    sc = new_scan(name, POLYFILE);

    /* make one mesh */

    sc->meshes[mesh_level] = (Mesh*) malloc(sizeof(Mesh));
    mesh = sc->meshes[mesh_level];

    mesh->ntris = 0;
    mesh->nverts = 0;

    /* read header information */

    result = fscanf(fp, "vertices: %d\n", &nverts);
    result = fscanf(fp, "faces: %d\n", &ntris);

    /* read in the vertices */

    mesh->nverts = 0;
    mesh->max_verts = nverts + 100;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    mesh->ntris = 0;
    mesh->max_tris = ntris + 100;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 200;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 200;
    mesh->parent_scan = sc;

    for (i = 0; i < nverts; i++) {
        fgets(str, 200, fp);
        result = sscanf(str, "v %f %f %f  %f %f %f\n",
                        &pos[X], &pos[Y], &pos[Z], &nx, &ny, &nz);
        make_vertex(mesh, pos);
    }

    /* read in the faces */

    for (i = 0; i < ntris; i++) {

        /* get string describing face */

        fgets(str, 200, fp);

        result = sscanf(str, "f %d %d %d\n", &in1, &in2, &in3);
        if (result == 0)
            break;

        v1 = mesh->verts[in1 - 1];
        v2 = mesh->verts[in2 - 1];
        v3 = mesh->verts[in3 - 1];
        make_triangle(mesh, v1, v2, v3, 100.0);
    }

    fclose(fp);

    /* print info about polygons */
    printf("%d triangles\n", mesh->ntris);
    printf("%d vertices\n", mesh->nverts);

    /* compute vertex normals */
    find_vertex_normals(mesh);

    /* make guess about what resolution this mesh was created at */
    inc = guess_mesh_inc(mesh);

    /* initialize hash table for vertices in mesh */
    init_table(mesh, TABLE_DIST * inc);

    /* find the edges of the mesh */
    find_mesh_edges(mesh);

    /* replicate this mesh at all levels */
    for (j = 0; j < MAX_MESH_LEVELS; j++)
        sc->meshes[j] = mesh;

    /* say we read the file okay */
    return (0);
}


/******************************************************************************
Take a guess at what mesh size this polygon file was created at.

Entry:
  mesh - mesh to guess about

Exit:
  returns number saying how many range data points are skipped between verts
******************************************************************************/

int guess_mesh_inc(mesh)
Mesh* mesh;
{
    int i, j;
    float num = 100;
    Vector diff;
    float len;
    float len_avg;
    Triangle* tri;
    int index;
    int inc;

    /* pick a bunch of random triangles and examine their edge lengths */
    len_avg = 0;
    for (i = 0; i < num; i++) {
        index = drand48() * mesh->ntris;
        tri = mesh->tris[index];
        for (j = 0; j < 3; j++) {
            vsub(tri->verts[j]->coord, tri->verts[(j + 1) % 3]->coord, diff);
            len = vlen(diff);
            len_avg += len / 3.0;
        }
    }
    len_avg /= num;

    /* guess at the mesh level based on this length */

    if (len_avg > 0.0044)
        inc = 8;
    else if (len_avg < 0.0044 && len_avg > 0.0022)
        inc = 4;
    else if (len_avg < 0.0022 && len_avg > 0.0011)
        inc = 2;
    else
        inc = 1;

    return (inc);
}


/******************************************************************************
Write the mesh of a scan to a binary polygon file.

Entry:
  scan     - which scan to write
  name - name of file to write to
******************************************************************************/

write_bin_polyfile(scan, name)
Scan* scan;
char* name;
{
    int i;
    FILE* fp;
    char filename[80];
    Mesh* mesh;
    Vertex** verts;
    int indices[3];
    float fconf;
    unsigned char conf;
    int conf_flag = 1;

    strcpy(filename, name);
    if (strlen(filename) < 4 ||
        strcmp(filename + strlen(filename) - 4, ".ply") != 0)
        strcat(filename, ".ply");

    fp = fopen(filename, "w");

    if (fp == NULL) {
        fprintf(stderr, "Can't open file '%s'\n", filename);
        return;
    }

    printf("Writing to '%s'\n", filename);

    mesh = scan->meshes[mesh_level];

    /* write the polygon header */

    /* magic number spells "ply<cr>" */
    fprintf(fp, "ply\n");

    /* number of properties */
    if (conf_flag) {
        fprintf(fp, "properties: 4\n");
    } else {
        fprintf(fp, "properties: 3\n");
    }

    /* number of vertices */
    fprintf(fp, "vertices: %d\n", mesh->nverts);

    /* number of triangles */
    fprintf(fp, "faces: %d\n", mesh->ntris);

    /* all faces are triangles */
    fprintf(fp, "face: all_triangles\n");

    /* all faces are triangles */
    if (conf_flag)
        fprintf(fp, "vertex: confidence\n");

    /* write out the vertices */
    for (i = 0; i < mesh->nverts; i++) {
        fwrite(mesh->verts[i]->coord, sizeof(float), 3, fp);
        /* maybe write out confidences of vertices */
        if (conf_flag) {
            fconf = mesh->verts[i]->confidence;
            if (fconf < 0)
                conf = 0;
            else if (fconf > 1)
                conf = 255;
            else
                conf = (unsigned char)(255.0 * fconf);
            fwrite(&conf, sizeof(unsigned char), 1, fp);
        }
    }

    /* write out the triangles */
    for (i = 0; i < mesh->ntris; i++) {
        verts = mesh->tris[i]->verts;
        indices[0] = verts[0]->index;
        indices[1] = verts[1]->index;
        indices[2] = verts[2]->index;
        fwrite(indices, sizeof(int), 3, fp);
    }

    fclose(fp);
}


/******************************************************************************
Read a binary polygon file into a scan.

Entry:
  filename - name of file to read in

Exit:
  returns 0 if file was read okay, 1 if not
******************************************************************************/

int read_bin_polyfile(filename)
char* filename;
{
    int i, j;
    FILE* fp;
    Scan* sc;
    Mesh* mesh;
    char* sval;
    char str[200];
    char sprop[200];
    char name[80];
    Vector pos;
    int result;
    Vertex* v1, *v2, *v3;
    int nverts, ntris;
    int inc;
    int nprops;
    int error_flag = 0;
    int indices[3];
    int found;
    int index;
    int conf_flag = 0;
    unsigned char conf;

    strcpy(name, filename);

    if (strlen(name) < 4 ||
        strcmp(name + strlen(name) - 4, ".ply") != 0)
        strcat(name, ".ply");

    fp = fopen(name, "r");

    if (fp == NULL) {
        printf("Couldn't open file '%s'\n", name);
        return (1);
    }

    printf("Reading polygons from '%s'\n", name);

    sc = new_scan(filename, POLYFILE);

    /* make one mesh */

    sc->meshes[mesh_level] = (Mesh*) malloc(sizeof(Mesh));
    mesh = sc->meshes[mesh_level];

    mesh->ntris = 0;
    mesh->nverts = 0;

    /* read header information */

    fscanf(fp, "ply\n");
    fscanf(fp, "properties: %d\n", &nprops);

    /* read properties */
    /* properties have the format "<prop_name>: <prop_value><cr>" */

    for (i = 0; i < nprops; i++) {

        fgets(str, 200, fp);

        /* search for colon */
        found = 0;
        for (j = 0; j < strlen(str); j++) {
            if (str[j] == ':') {
                sval = &str[j + 1];
                sprop[j] = '\0';
                found = 1;
                break;
            } else
                sprop[j] = str[j];
        }

        /* error check */
        if (found == 0) {
            fprintf(stderr, "can't parse string: %s\n", str);
            return (1);
        }

        /* match the property name and read the property value */
        if (strcmp(sprop, "vertices") == 0) {
            nverts = atoi(sval);
        } else if (strcmp(sprop, "faces") == 0) {
            ntris = atoi(sval);
        } else if (strcmp(sprop, "face") == 0) {
            /* sval should say "all_triangles" */
        } else if (strcmp(sprop, "vertex") == 0) {
            if (strcmp(sval, " confidence\n") == 0) {
                conf_flag = 1;
            }
        } else {
            fprintf(stderr, "can't parse string: %s\n", str);
            return (1);
        }
    }

#if 0
    fscanf(fp, "vertices: %d\n", &nverts);
    fscanf(fp, "faces: %d\n", &ntris);
    fscanf(fp, "face: all_triangles\n");
#endif

    /* read in the vertices */

    mesh->nverts = 0;
    mesh->max_verts = nverts + 100;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    mesh->ntris = 0;
    mesh->max_tris = ntris + 100;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 200;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 200;
    mesh->parent_scan = sc;

    for (i = 0; i < nverts; i++) {
        result = fread(pos, sizeof(float), 3, fp);
        if (result == 0) {
            error_flag = 1;
            break;
        }
        index = make_vertex(mesh, pos);
        if (conf_flag) {
            fread(&conf, sizeof(unsigned char), 1, fp);
            mesh->verts[index]->confidence = ((int) conf) / 255.0;
        }
    }

    if (error_flag) {
        printf("error reading file\n");
        free(mesh->verts);
        free(mesh->tris);
        free(mesh->edges);
        return (1);
    }

    /* read in the faces */

    for (i = 0; i < ntris; i++) {
        result = fread(indices, sizeof(int), 3, fp);
        if (result == 0) {
            error_flag = 1;
            break;
        }
#if 0
        if (indices[0] > nverts || indices[1] > nverts || indices[2] > nverts)
            continue;
#endif
        v1 = mesh->verts[indices[0]];
        v2 = mesh->verts[indices[1]];
        v3 = mesh->verts[indices[2]];
        make_triangle(mesh, v1, v2, v3, 1e20);
    }

    if (error_flag) {
        printf("error reading file\n");
        free(mesh->verts);
        free(mesh->tris);
        free(mesh->edges);
        return (1);
    }

    fclose(fp);

    /* print info about polygons */
    printf("%d triangles\n", mesh->ntris);
    printf("%d vertices\n", mesh->nverts);

    /* compute vertex normals */
    find_vertex_normals(mesh);

    /* make guess about what resolution this mesh was created at */
    inc = guess_mesh_inc(mesh);

    /* initialize hash table for vertices in mesh */
    init_table(mesh, TABLE_DIST * inc);

    /* find the edges of the mesh */
    find_mesh_edges(mesh);

    /* replicate this mesh at all levels */
    for (j = 0; j < MAX_MESH_LEVELS; j++)
        sc->meshes[j] = mesh;

    /* say we read the file okay */
    return (0);
}


/******************************************************************************
Create a new scan.

Entry:
  name - name of scan
  type - type of scan (POLYFILE, etc.)

Exit:
  returns pointer to the new scan
******************************************************************************/

Scan* new_scan(char* name, int type)
{
    int j, k;
    Scan* sc;

    scans[nscans] = (Scan*) malloc(sizeof(Scan));
    sc = scans[nscans];
    nscans++;

    strcpy(sc->name, name);

    sc->xtrans = 0;
    sc->ytrans = 0;
    sc->ztrans = 0;
    sc->rotate = 0;

    build_rotmat(sc);

    sc->draw_flag = 1;
    sc->edge_mesh = NULL;
    sc->file_type = type;

    for (j = 0; j < 4; j++)
        for (k = 0; k < 4; k++)
            sc->rotmat[j][k] = (j == k);

    /* zero out mesh entries */

    for (j = 0; j < MAX_MESH_LEVELS; j++)
        sc->meshes[j] = NULL;

    /* return pointer to new scan */
    return (sc);
}

