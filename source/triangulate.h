/*
 * Copyright (c) 1995-2017, Stanford University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of Stanford University nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY STANFORD UNIVERSITY ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL STANFORD UNIVERSITY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ZIPPER_TRIANGULATE_H
#define ZIPPER_TRIANGULATE_H

// Internal
#include "zipper.h"
#include "matrix.h"

#define BORDER       1
#define INTERIOR     2
#define CREATED      3
#define OLD_CREATED  4

#define NPMAX 40

typedef struct Npoly {
    int nverts;           /* number of vertices */
    unsigned char border[NPMAX];  /* if vertex is part of border polygon */
    int index[NPMAX];     /* index of vertex */
    float x[NPMAX];
    float y[NPMAX];
    int which_tri;        /* index of which triangle it belongs to */
} Npoly;

typedef struct TriangulateEdge {
    int p1, p2;
    float len;
    float a, b, c;
    int final;
} TriangulateEdge;

typedef struct TriangulatePoint {
    Vector pos;           /* location in plane */
    Vector pos3d;         /* old location in 3-space */
    int boundary;
    struct TriangulateEdge** edges;
    int nedges;
    int max_edges;
    int index;
} TriangulatePoint;

typedef struct TriangulateTriangle {
    int p1, p2, p3;
    TriangulateEdge* e1, *e2, *e3;
} TriangulateTriangle;

// Declarations
int init_splitter(float a, float b, float c, float d);
int add_boundary_point(float xx, float yy, float zz, int index);
int add_point(float xx, float yy, float zz, int index);
int split_npoints();
int close_orig(TriangulatePoint* p, float x, float y, float z);
void set_parallel_flag(int val);
void set_shuffle_flag(int val);
void set_rescale_flag(int val, int x, int y);
void add_edge(int i, int j);
void add_final_edge(TriangulateEdge* e);
void byte_copy(char* dst, char* src, int num);
void shuffle(char* list, int num, int size);
int edge_compare(const void* p1, const void* p2);
int inside_boundary(TriangulateEdge* e);
int nearly_on_edge(TriangulateEdge* e);
int any_intersection(TriangulateEdge* e);
int greedy_connect();
int collect_triangles(int whoops_flag);
void maybe_make_tri(int p1, int p2, int p3, TriangulateEdge* e1, TriangulateEdge* e2, TriangulateEdge* e3);
void new_not_used_orient_triangles();
void orient_triangles();
int triangle_direction(TriangulateTriangle* tri);
void flip_triangle(TriangulateTriangle* tri);
int get_ntris();
int get_triangle(int num, int* p1, int* p2, int* p3);
void print_poly();
int point_in_which_triangle(float x, float y, float* b1, float* b2, float* b3);
void compute_line(float x1, float y1, float x2, float y2, float* aa, float* bb, float* cc);
void rescale_points();
int fold_in_poly_check(float x, float y);
int point_in_split_poly(float x, float y);
int point_in_poly(float x, float y, int cnt, TriangulatePoint** polypts);
int whichquad(Vector pt, float x, float y);
int face_to_xy_plane(float a, float b, float c, float d, Matrix mat, Matrix imat);
void reorder_edges();

#endif
