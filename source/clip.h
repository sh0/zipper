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

#ifndef ZIPPER_CLIP_H
#define ZIPPER_CLIP_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Parameters
void update_clip_resolution();
void set_clip_near_dist_factor(float factor);
float get_clip_near_dist_factor();
void set_clip_near_cos(float cosine);
float get_clip_near_cos();
void set_clip_boundary_dist_factor(float factor);
float get_clip_boundary_dist_factor();
void set_clip_boundary_cos(float cosine);
float get_clip_boundary_cos();

// Declarations
void clip_triangles(Scan* sc1, Scan* sc2);
void perform_triangle_clipping(Scan* sc1, Scan* sc2);
void process_vertices(Vector tnorm, int tindex, Clip_List* clist, Mesh* mesh);
int outside_mesh(Clip_List* clist);
void list_to_tris(Triangle* tri, Clip_List* clist, Mesh* mesh);
void new_list_to_tris(Vector tnorm, int index, Clip_List* clist, Mesh* mesh);
Clip_List* split_list(Clip_List* clist, int index1, int index2);
Clip_List* make_between_list(Clip_List* clist, Cut* cut1, Cut* cut2);
Clip_List* potential_vertices(Triangle* tri);
int find_partner_cuts(Triangle* tri, Clip_List* clist, Mesh* mesh);
Cut* next_similar_cut(Cut* cut, Clip_List* clist, int dir);
void sort_triangle_cuts(Triangle* tri);
void cut_triangle(Triangle* tri, Mesh* m1, Mesh* m2, Scan* scan);
void old_cut_triangle(Triangle* tri, Mesh* m1, Mesh* m2, Scan* scan);
void point_line_approach(Vector pt, Vector q1, Vector q2, Vector x, float* t);
int two_line_approach(Vector p1, Vector q1, Vector p2, Vector q2, Vector x1, Vector x2, float* t1, float* t2);
void verts_near_edges(Mesh* mesh, Mesh* not_mesh, Vector pnt, Vector norm, float radius, float min_dot);
void edges_near_edges(Triangle* tri, Mesh* m1, Mesh* m2, Scan* scan);
void make_clip_triangles(Scan* scan, Mesh* clipto);
int line_intersect_tri_single(Vector p1, Vector p2, Triangle* tri, Vector pos, float* tt, int* inward, Vector barycentric);
int line_intersect_tri(Vector p1, Vector p2, Triangle* tri, Vector pos, float* tt, int* inward, Vector barycentric);
float point_project_line(Vector v1, Vector v2, Vector v3, Vector p);
int new_cut(Edge* edge, Vertex* v1, Vertex* v2, float t, float s, int inward);
void add_cut_to_triangle(Triangle* tri, Cut* cut);
void init_cuts(Scan* scan, Mesh* clipto);
void create_cut_vertices(Mesh* mesh);
void introduce_all_cuts(Mesh* mesh, Mesh* not_mesh);
void introduce_cuts(Triangle* tri, Edge* edge, Mesh* mesh, Triangle** first_tri, Triangle** last_tri);
void sort_cuts(Edge* edge);
int plane_thru_vectors_double(double* v0, double* v1, double* v2, double* aa, double* bb, double* cc, double* dd);
void double_stuff(Triangle* tri);

#endif
