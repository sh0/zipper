/*

Routines for reading and writing to PLY polygon files.

Greg Turk, April 1994

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "matrix.h"
#include "zipper.h"
#include "raw.h"
#include "ply.h"

static float RANGE_DATA_SIGMA_FACTOR;
static float RANGE_DATA_MIN_INTENSITY;
static int RANGE_DATA_HORIZONTAL_ERODE;


set_range_data_sigma_factor(factor)
float factor;
{
    RANGE_DATA_SIGMA_FACTOR = factor;
}


float
get_range_data_sigma_factor()
{
    return RANGE_DATA_SIGMA_FACTOR;
}

set_range_data_min_intensity(intensity)
float intensity;
{
    RANGE_DATA_MIN_INTENSITY = intensity;
}


float
get_range_data_min_intensity()
{
    return RANGE_DATA_MIN_INTENSITY;
}

set_range_data_horizontal_erode(erode)
int erode;
{
    RANGE_DATA_HORIZONTAL_ERODE = erode;
}


int
get_range_data_horizontal_erode()
{
    return RANGE_DATA_HORIZONTAL_ERODE;
}


enum {CONTINUOUS, JUMP_CLOSER, JUMP_FARTHER, LEFT_EDGE, RIGHT_EDGE};

int write_intensity_flag = 1;   /* whether to write out vertex intensities */
int write_normals_flag = 0;     /* write out vertex normals? */
int write_colors_flag = 0;      /* write out diffuse colors? */

Scan* new_scan(char*, int);

/* list of the elements in the PLY files */
char* elem_names[] = {
    "vertex", "face"
};


/* list of property information for a vertex */
PlyProperty vert_props[] = {
    {"x", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"y", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"z", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"confidence", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"intensity", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"nx", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"ny", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"nz", PLY_FLOAT, PLY_FLOAT, 0, 0, 0, 0, 0},
    {"mesh_tags", PLY_UINT, PLY_UINT, 0, 0, 0, 0, 0},
    {"diffuse_red", PLY_UCHAR, PLY_UCHAR, 0, 0, 0, 0, 0},
    {"diffuse_green", PLY_UCHAR, PLY_UCHAR, 0, 0, 0, 0, 0},
    {"diffuse_blue", PLY_UCHAR, PLY_UCHAR, 0, 0, 0, 0, 0},
};

typedef struct A_Vertex {
    float x, y, z;
    float confidence;
    float intensity;
    unsigned char red, grn, blu;
} A_Vertex;

PlyProperty vert_std_props[] = {
    {"x", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, x), 0, 0, 0, 0},
    {"y", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, y), 0, 0, 0, 0},
    {"z", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, z), 0, 0, 0, 0},
    {"confidence", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, confidence), 0, 0, 0, 0},
    {"intensity", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, intensity), 0, 0, 0, 0},
    {"diffuse_red", PLY_UCHAR, PLY_UCHAR, offsetof(A_Vertex, red), 0, 0, 0, 0},
    {"diffuse_green", PLY_UCHAR, PLY_UCHAR, offsetof(A_Vertex, grn), 0, 0, 0, 0},
    {"diffuse_blue", PLY_UCHAR, PLY_UCHAR, offsetof(A_Vertex, blu), 0, 0, 0, 0},
    {"std_dev", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, confidence), 0, 0, 0, 0},
};

/* list of property information for a vertex */
#if 0
PlyProperty tri_props[] = {
    {"vertex_indices", PLY_INT, PLY_INT, 0, 1, PLY_UCHAR, PLY_UCHAR, 0},
};
#else
PlyProperty tri_props[] = {
    {"vertex_index", PLY_INT, PLY_INT, 0, 1, PLY_UCHAR, PLY_UCHAR, 0},
};
#endif

typedef struct RangePnt {
    unsigned char num_pts;
    int* pts;
} RangePnt;

/* list of property information for a range data point */
PlyProperty range_props[] = {
    {
        "vertex_indices", PLY_INT, PLY_INT, offsetof(RangePnt, pts),
        1, PLY_UCHAR, PLY_UCHAR, offsetof(RangePnt, num_pts)
    },
};

typedef struct A_Triangle {
    int* verts;
    unsigned char nverts;
} A_Triangle;

static PlyProperty avert_props[] = {
    {"x", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, x), 0, 0, 0, 0},
    {"y", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, y), 0, 0, 0, 0},
    {"z", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, z), 0, 0, 0, 0},
    {"confidence", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, confidence), 0, 0, 0, 0},
    {"intensity", PLY_FLOAT, PLY_FLOAT, offsetof(A_Vertex, intensity), 0, 0, 0, 0},  {"diffuse_red", PLY_UCHAR, PLY_UCHAR, offsetof(A_Vertex, red), 0, 0, 0, 0},
    {"diffuse_green", PLY_UCHAR, PLY_UCHAR, offsetof(A_Vertex, grn), 0, 0, 0, 0},
    {"diffuse_blue", PLY_UCHAR, PLY_UCHAR, offsetof(A_Vertex, blu), 0, 0, 0, 0},
};

static Vertex* vert_dummy;
#define voffset(field) ((char *) (&vert_dummy->field) - (char *) vert_dummy)
static A_Triangle* tri_dummy;
#define toffset(field) ((char *) (&tri_dummy->field) - (char *) tri_dummy)



/******************************************************************************
Write out polygons to a PLY file.

Entry:
  sc       - scan to write from
  filename - name of PLY file to write to
******************************************************************************/

void write_ply(Scan* sc, char* filename, int writeInfo)
{
    int i;
    PlyFile* ply;
    float version;
    Mesh* mesh;
    A_Triangle atri;

    mesh = sc->meshes[mesh_level];

    atri.verts = (int*) malloc(sizeof(int) * 3);
    atri.nverts = 3;

    /* set up the offsets into the Vertex and Triangle structures */
    vert_props[0].offset = voffset(coord[0]);
    vert_props[1].offset = voffset(coord[1]);
    vert_props[2].offset = voffset(coord[2]);
    vert_props[3].offset = voffset(confidence);
    vert_props[4].offset = voffset(intensity);
    vert_props[5].offset = voffset(normal[0]);
    vert_props[6].offset = voffset(normal[1]);
    vert_props[7].offset = voffset(normal[2]);
    vert_props[8].offset = voffset(old_mesh);
    vert_props[9].offset = voffset(red);
    vert_props[10].offset = voffset(grn);
    vert_props[11].offset = voffset(blu);

    tri_props[0].offset = toffset(verts);
    tri_props[0].count_offset = toffset(nverts);

    /* set up PLY header information */
    ply = ply_open_for_writing(filename, 2, elem_names, PLY_BINARY_LE, &version);

#if 1

    ply_element_count(ply, "vertex", mesh->nverts);

    ply_describe_property(ply, "vertex", &vert_props[0]);
    ply_describe_property(ply, "vertex", &vert_props[1]);
    ply_describe_property(ply, "vertex", &vert_props[2]);
    ply_describe_property(ply, "vertex", &vert_props[3]);

    if (write_intensity_flag) {
        ply_describe_property(ply, "vertex", &vert_props[4]);
    }

    if (write_normals_flag) {
        ply_describe_property(ply, "vertex", &vert_props[5]);
        ply_describe_property(ply, "vertex", &vert_props[6]);
        ply_describe_property(ply, "vertex", &vert_props[7]);
        ply_describe_property(ply, "vertex", &vert_props[8]);
    }

    if (write_colors_flag) {
        ply_describe_property(ply, "vertex", &vert_props[9]);
        ply_describe_property(ply, "vertex", &vert_props[10]);
        ply_describe_property(ply, "vertex", &vert_props[11]);
    }

#endif

#if 0
    if (write_normals_flag)
        ply_describe_element(ply, "vertex", mesh->nverts, 9, vert_props);
    else if (write_intensity_flag)
        ply_describe_element(ply, "vertex", mesh->nverts, 5, vert_props);
    else
        ply_describe_element(ply, "vertex", mesh->nverts, 4, vert_props);
#endif

    ply_describe_element(ply, "face", mesh->ntris, 1, tri_props);

    ply_put_comment(ply, "zipper output");

    if (writeInfo) {
        for (i = 0; i < sc->num_obj_info; i++)
            ply_put_obj_info(ply, sc->obj_info[i]);

    }

    ply_header_complete(ply);

    /* write out the vertices */
    ply_put_element_setup(ply, "vertex");
    for (i = 0; i < mesh->nverts; i++)
        ply_put_element(ply, (void*) mesh->verts[i]);

    /* write out the triangles */
    ply_put_element_setup(ply, "face");
    for (i = 0; i < mesh->ntris; i++) {
        atri.verts[0] = mesh->tris[i]->verts[0]->index;
        atri.verts[1] = mesh->tris[i]->verts[1]->index;
        atri.verts[2] = mesh->tris[i]->verts[2]->index;
        ply_put_element(ply, (void*) &atri);
    }

    ply_close(ply);
}


int is_range_grid_file(char* filename)
{
    int i;
    PlyFile* ply;
    int nelems;
    char** elist;
    int file_type;
    float version;

    ply = ply_open_for_reading(filename, &nelems, &elist, &file_type, &version);
    if (ply == NULL)
        return 0;
    ply_close(ply);

    for (i = 0; i < nelems; i++) {
        if (!strcmp(elist[i], "range_grid"))
            return 1;
    }

    return 0;
}


/******************************************************************************
Read in polygons from a PLY file.

Entry:
  sc       - scan to read into
  filename - name of PLY file to read from

Exit:
  returns 0 if scan created okay, 1 if there was an error
******************************************************************************/

int read_ply(char* filename)
{
    int i, j;
    PlyFile* ply;
    int nelems;
    char** elist;
    int file_type;
    float version;
    int num_elems;
    PlyProperty** plist;
    int nprops;
    A_Vertex vert;
    Vector pos;
    A_Triangle atri;
    char* elem_name;
    Scan* sc;
    Mesh* mesh;
    int inc;
    int index;
    Vertex* v1, *v2, *v3;
    int get_intensity = 0;
    int get_confidence = 0;
    int has_red = 0;
    int has_grn = 0;
    int has_blu = 0;
    char** obj_info;
    int num_obj_info;

    tri_props[0].offset = toffset(verts);
    tri_props[0].count_offset = toffset(nverts);

    ply = ply_open_for_reading(filename, &nelems, &elist, &file_type, &version);
    if (ply == NULL)
        return (1);

    obj_info = ply_get_obj_info(ply, &num_obj_info);

    sc = new_scan(filename, POLYFILE);

    sc->num_obj_info = num_obj_info;
    for (i = 0; i < num_obj_info; i++) {
        strcpy(sc->obj_info[i], obj_info[i]);
    }


    /* make one mesh */

    sc->meshes[mesh_level] = (Mesh*) malloc(sizeof(Mesh));
    mesh = sc->meshes[mesh_level];

    /* read in the vertices */

    plist = ply_get_element_description(ply, "vertex", &num_elems, &nprops);
    mesh->nverts = 0;
    mesh->max_verts = num_elems + 100;
    mesh->verts = (Vertex**) malloc(sizeof(Vertex*) * mesh->max_verts);

    plist = ply_get_element_description(ply, "face", &num_elems, &nprops);
    mesh->ntris = 0;
    mesh->max_tris = num_elems + 100;
    mesh->tris = (Triangle**) malloc(sizeof(Triangle*) * mesh->max_tris);

    mesh->nedges = 0;
    mesh->max_edges = 200;
    mesh->edges = (Edge**) malloc(sizeof(Edge*) * mesh->max_edges);
    mesh->edges_valid = 0;
    mesh->eat_list_max = 200;
    mesh->parent_scan = sc;

    /* read in the vertices and triangles */

    for (i = 0; i < nelems; i++) {

        elem_name = elist[i];
        plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

        if (equal_strings("vertex", elem_name)) {

            /* see if the file contains intensities */
            for (j = 0; j < nprops; j++) {
                if (equal_strings("confidence", plist[j]->name))
                    get_confidence = 1;
                if (equal_strings("intensity", plist[j]->name))
                    get_intensity = 1;
                if (equal_strings("diffuse_red", plist[j]->name))
                    has_red = 1;
                if (equal_strings("diffuse_green", plist[j]->name))
                    has_grn = 1;
                if (equal_strings("diffuse_blue", plist[j]->name))
                    has_blu = 1;
            }

            ply_get_property(ply, elem_name, &avert_props[0]);
            ply_get_property(ply, elem_name, &avert_props[1]);
            ply_get_property(ply, elem_name, &avert_props[2]);
            if (get_confidence)
                ply_get_property(ply, elem_name, &avert_props[3]);
            if (get_intensity)
                ply_get_property(ply, elem_name, &avert_props[4]);
            if (has_red && has_grn && has_blu) {
                ply_get_property(ply, elem_name, &avert_props[5]);
                ply_get_property(ply, elem_name, &avert_props[6]);
                ply_get_property(ply, elem_name, &avert_props[7]);
            }

            for (j = 0; j < num_elems; j++) {
                ply_get_element(ply, (void*) &vert);
                pos[X] = vert.x;
                pos[Y] = vert.y;
                pos[Z] = vert.z;
                index = make_vertex(mesh, pos);
                mesh->verts[index]->confidence = vert.confidence;
                mesh->verts[index]->intensity = vert.intensity;
                mesh->verts[index]->red = vert.red;
                mesh->verts[index]->grn = vert.grn;
                mesh->verts[index]->blu = vert.blu;
            }
        } else if (equal_strings("face", elem_name)) {
            ply_get_element_setup(ply, elem_name, 1, tri_props);
            for (j = 0; j < num_elems; j++) {
                ply_get_element(ply, (void*) &atri);
                v1 = mesh->verts[atri.verts[0]];
                v2 = mesh->verts[atri.verts[1]];
                v3 = mesh->verts[atri.verts[2]];
                make_triangle(mesh, v1, v2, v3, 100.0);
            }
        }
    }

    /* print info about polygons */
    printf("%d triangles\n", mesh->ntris);
    printf("%d vertices\n", mesh->nverts);

    /* compute vertex normals */
    find_vertex_normals(mesh);

    /* make guess about what resolution this mesh was created at */
    inc = guess_mesh_inc(mesh);

    /* initialize hash table for vertices in mesh */
    init_table(mesh, 2.0f * get_zipper_resolution() * inc);

    /* find the edges of the mesh */
    find_mesh_edges(mesh);

    /* replicate this mesh at all levels */
    for (j = 0; j < MAX_MESH_LEVELS; j++)
        sc->meshes[j] = mesh;

    /* signal that everything went okay */
    return (0);
}


/******************************************************************************
Read range data from a PLY file.

Entry:
  name - name of PLY file to read from

Exit:
  returns pointer to data, or NULL if it couldn't read from file
******************************************************************************/

RangeData* read_ply_geom(name)
char* name;
{
    int i, j, k, index, best_index;
    int xx, yy;
    PlyFile* ply;
    RangeData* plydata = NULL;
    char** obj_info;
    int num_obj_info;
    int num_elems;
    int nprops;
    int nelems;
    char** elist;
    int file_type;
    float version;
    PlyProperty** plist;
    int num_rows, num_cols;
    RangePnt range_pnt;
    A_Vertex vert;
    char* elem_name;
    int get_std_dev = 0;
    int get_confidence = 0;
    int get_intensity = 0;
    int get_color = 0;
    int has_red = 0;
    int has_green = 0;
    int has_blue = 0;
    int is_warped = 0;
    char cont;
    int erodeMax, is_cyberware_scanner_data, count;
    char* continuity;


    float conf, std;
    float min_std, max_std, max;
    float avg_std = 0;

    ply = ply_open_for_reading(name, &nelems, &elist, &file_type, &version);
    if (ply == NULL)
        return (NULL);

    /* parse the obj_info */

    obj_info = ply_get_obj_info(ply, &num_obj_info);
#if 0
    for (i = 0; i < num_obj_info; i++) {
        Tcl_SplitList(interp, obj_info[i], &argc, &argv);
        if (equal_strings("num_cols", argv[0]))
            num_cols = atoi(argv[1]);
        if (equal_strings("num_rows", argv[0]))
            num_rows = atoi(argv[1]);
        if (equal_strings("is_warped", argv[0]))
            is_warped = atoi(argv[1]);
        if (equal_strings("optimum_std_dev", argv[0]))
            avg_std = atof(argv[1]);
    }
#endif

    min_std = avg_std / RANGE_DATA_SIGMA_FACTOR;
    max_std = avg_std * RANGE_DATA_SIGMA_FACTOR;

    /* set up the range data structure */
    plydata = (RangeData*) malloc(sizeof(RangeData));
    plydata->nlg = num_rows;
    plydata->nlt = num_cols;
    plydata->has_color = 0;
    plydata->has_intensity = 0;
    plydata->has_confidence = 0;
    plydata->mult_confidence = 0;

    plydata->num_obj_info = num_obj_info;
    for (i = 0; i < num_obj_info; i++) {
        strcpy(plydata->obj_info[i], obj_info[i]);
    }

    /* see if we've got both vertex and range_grid data */

    plist = ply_get_element_description(ply, "vertex", &num_elems, &nprops);
    if (plist == NULL) {
        fprintf(stderr, "file doesn't contain vertex data\n");
        return (NULL);
    }
    plydata->num_points = num_elems;
    plist = ply_get_element_description(ply, "range_grid", &num_elems, &nprops);
    if (plist == NULL) {
        fprintf(stderr, "file doesn't contain range_grid data\n");
        return (NULL);
    }

    plydata->points = (Vector*) malloc(sizeof(Vector) * plydata->num_points);
    plydata->confidence = (float*) malloc(sizeof(float) * plydata->num_points);
    plydata->intensity = (float*) malloc(sizeof(float) * plydata->num_points);
    plydata->red = (unsigned char*)
                   malloc(sizeof(unsigned char) * plydata->num_points);
    plydata->grn = (unsigned char*)
                   malloc(sizeof(unsigned char) * plydata->num_points);
    plydata->blu = (unsigned char*)
                   malloc(sizeof(unsigned char) * plydata->num_points);
    plydata->pnt_indices = (int*) malloc
                           (sizeof(int) * plydata->nlt * plydata->nlg);

    /* read in the range data */

    for (i = 0; i < nelems; i++) {

        elem_name = elist[i];
        plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

        if (equal_strings("vertex", elem_name)) {

            /* see if the file contains intensities */
            for (j = 0; j < nprops; j++) {
                if (equal_strings("std_dev", plist[j]->name))
                    get_std_dev = 1;
                if (equal_strings("confidence", plist[j]->name))
                    get_confidence = 1;
                if (equal_strings("intensity", plist[j]->name))
                    get_intensity = 1;
                if (equal_strings("diffuse_red", plist[j]->name))
                    has_red = 1;
                if (equal_strings("diffuse_green", plist[j]->name))
                    has_green = 1;
                if (equal_strings("diffuse_blue", plist[j]->name))
                    has_blue = 1;
            }

            if (has_red && has_green && has_blue) {
                get_color = 1;
                plydata->has_color = 1;
            }

            if (get_intensity)
                plydata->has_intensity = 1;

            if (get_std_dev && is_warped) {
                plydata->has_confidence = 1;
                plydata->mult_confidence = 1;
            } else if (get_confidence)
                plydata->has_confidence = 1;

            ply_get_property(ply, "vertex", &vert_std_props[0]);
            ply_get_property(ply, "vertex", &vert_std_props[1]);
            ply_get_property(ply, "vertex", &vert_std_props[2]);
            if (get_confidence)
                ply_get_property(ply, "vertex", &vert_std_props[3]);
            if (get_intensity)
                ply_get_property(ply, "vertex", &vert_std_props[4]);
            if (get_color) {
                ply_get_property(ply, "vertex", &vert_std_props[5]);
                ply_get_property(ply, "vertex", &vert_std_props[6]);
                ply_get_property(ply, "vertex", &vert_std_props[7]);
            }
            if (get_std_dev)
                ply_get_property(ply, "vertex", &vert_std_props[8]);

            for (j = 0; j < num_elems; j++) {
                ply_get_element(ply, (void*) &vert);
                plydata->points[j][X] = vert.x;
                plydata->points[j][Y] = vert.y;
                plydata->points[j][Z] = vert.z;

                if (get_intensity) {
                    plydata->intensity[j] = vert.intensity;
                } else
                    plydata->intensity[j] = 0.5;

                if (get_std_dev && is_warped) {

                    std = vert.confidence;

                    if (std < min_std)
                        conf = 0;
                    else if (std < avg_std)
                        conf = (std - min_std) / (avg_std - min_std);
                    else if (std > max_std)
                        conf = 0;
                    else
                        conf = (max_std - std) / (max_std - avg_std);

                    /*
                           Unsafe to use vertex intensity, as aperture settings may change
                           between scans.  Instead, use std_dev confidence * orientation.

                          conf *= vert.intensity;
                    */
                    if (vert.intensity < RANGE_DATA_MIN_INTENSITY) conf = 0.0;

                    plydata->confidence[j] = conf;
                } else if (get_confidence) {
                    plydata->confidence[j] = vert.confidence;
                } else if (plydata->has_intensity &&
                           vert.intensity < RANGE_DATA_MIN_INTENSITY) {
                    plydata->confidence[j] = 0.0;
                } else {
                    plydata->confidence[j] = 0.5;
                }

                if (get_color) {
                    plydata->red[j] = vert.red;
                    plydata->grn[j] = vert.grn;
                    plydata->blu[j] = vert.blu;
                }

            }
        }

        if (equal_strings("range_grid", elem_name)) {
            ply_get_element_setup(ply, elem_name, 1, range_props);
            for (j = 0; j < num_elems; j++) {
                ply_get_element(ply, (void*) &range_pnt);
                if (range_pnt.num_pts == 0)
                    plydata->pnt_indices[j] = -1;
                else {
                    max = -FLT_MAX;
                    for (k = 0; k < range_pnt.num_pts; k++) {
                        index = range_pnt.pts[k];
                        if (plydata->intensity[index] > max) {
                            max = plydata->intensity[index];
                            best_index = index;
                        }
                    }
                    index = best_index;
                    if (plydata->confidence[index] > 0.0)
                        plydata->pnt_indices[j] = index;
                    else
                        plydata->pnt_indices[j] = -1;
                }
            }

            /*
            Check each range grid point to see if its neighbors
             within a certain number of samples are edges; if so,
             mark to be nullified.  For the raw Cyberware data, we
             should be careful about removing "OK" data due to
             detection of only one peak per scanline (the one closer to
             the scanner)
             */

#if 1
            /* Make an auxilliary array indicating connectivity */
            continuity = (char*)malloc(plydata->nlg * plydata->nlt);
            index = 0;
            for (yy = 0; yy < plydata->nlg; yy++) {
                for (xx = 0; xx < plydata->nlt; xx++, index++) {
                    continuity[index] =
                        decide_continuity(plydata, xx, yy);
                }
            }


            /* Fill these with something meaningful !!!*/
            erodeMax = RANGE_DATA_HORIZONTAL_ERODE;
            is_cyberware_scanner_data = 0;

#if 0
            if (erodeMax > 0)
                printf("Eroding horizontally by %d triangles.\n", erodeMax);
#endif

            index = 0;
            count = 0;
            for (yy = 0; yy < plydata->nlg; yy++) {
                for (xx = 0; xx < plydata->nlt; xx++, index++) {
                    cont = continuity[index];
                    if (cont == LEFT_EDGE) {
                        count = erode_forward(plydata, continuity, xx, yy, erodeMax);
                        xx += count;
                        index += count;
                        if (xx < plydata->nlt)
                            cont = continuity[index];
                    }
                    if (cont == JUMP_CLOSER) {
                        if (!is_cyberware_scanner_data) {
                            erode_backward(plydata, continuity, xx, yy, erodeMax);
                        }
                        count = erode_forward(plydata, continuity, xx, yy, erodeMax);
                        xx += count;
                        index += count;
                        if (xx < plydata->nlt)
                            cont = continuity[index];
                    }
                    if (cont == JUMP_FARTHER) {
                        erode_backward(plydata, continuity, xx, yy, erodeMax);
                        if (!is_cyberware_scanner_data) {
                            count = erode_forward(plydata, continuity, xx, yy, erodeMax);
                            xx += count;
                            index += count;
                            if (xx < plydata->nlt)
                                cont = continuity[index];
                        }
                    }
                    if (cont == RIGHT_EDGE) {
                        erode_backward(plydata, continuity, xx, yy, erodeMax);
                    }
                }
            }
            free(continuity);
#endif
        }
    }

    ply_close(ply);

    return (plydata);
}


int
erode_forward(plydata, continuity, xx, yy, erodeMax)
RangeData* plydata;
char* continuity;
int xx, yy, erodeMax;
{
    int erode_count, index;
    char cont;

    erode_count = 0;
    index = xx + yy * plydata->nlt;
    xx++;
    index++;

    if (xx < plydata->nlt)
        cont = continuity[index];
    else
        return 0;


    while (cont == CONTINUOUS && xx < plydata->nlt &&
           erode_count < erodeMax) {
        plydata->pnt_indices[index] = -1;
        xx++;
        index++;
        if (xx < plydata->nlt)
            cont = continuity[index];
        erode_count++;
    }
    return erode_count;
}


int
erode_backward(plydata, continuity,
               xx, yy, erodeMax)
RangeData* plydata;
char* continuity;
int xx, yy, erodeMax;
{
    int erode_count, index;
    char cont;

    erode_count = 0;
    index = xx + yy * plydata->nlt;
    while (cont == CONTINUOUS && xx >= 0 &&
           erode_count < erodeMax) {
        plydata->pnt_indices[index] = -1;
        xx--;
        index--;
        if (xx >= 0)
            cont = continuity[index];
        erode_count++;
    }
    return erode_count;
}


int
decide_continuity(plydata, xx, yy)
RangeData* plydata;
int xx, yy;
{
    Vector diff;
    int index, vert_index, next_vert_index;

    index = xx + yy * plydata->nlt;
    vert_index = plydata->pnt_indices[index];

    if (xx != plydata->nlt - 1)
        next_vert_index = plydata->pnt_indices[index + 1];
    else
        next_vert_index = -1;

    if (vert_index == -1) {
        if (next_vert_index != -1)
            return LEFT_EDGE;
        else
            return CONTINUOUS;
    } else {
        if (next_vert_index == -1)
            return RIGHT_EDGE;
        else {
            my_vsub(plydata->points[vert_index],
                    plydata->points[next_vert_index], diff);
            if (vlen(diff) > edge_length_max(0)) {
                if (diff[Z] > 0)
                    return JUMP_FARTHER;
                else
                    return JUMP_CLOSER;
            } else
                return CONTINUOUS;
        }
    }
}


void
delete_ply_geom(RangeData* plydata)
{
    free(plydata->points);
    free(plydata->confidence);
    free(plydata->intensity);
    free(plydata->red);
    free(plydata->grn);
    free(plydata->blu);
    free(plydata->pnt_indices);
    free(plydata);
}



