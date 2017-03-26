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

#ifndef ZIPPER_MESH_H
#define ZIPPER_MESH_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Parameters
void set_conf_edge_count_factor(float factor);
float get_conf_angle();
void set_conf_angle(float factor);
float get_conf_exponent();
void set_conf_exponent(float factor);
float get_conf_edge_count_factor();
void set_conf_edge_zero(int set);
int get_conf_edge_zero();

// Declarations
void create_scan_mesh(Scan* sc, int level);
Mesh* make_mesh_raw(Scan* sc, int level, float table_dist);
Mesh* make_mesh_ply(Scan* sc, int level, float table_dist);
void vertex_errors(Mesh* mesh, Scan* scan, int rot_flag, int mult);
void lower_edge_confidence(Mesh* mesh, int level);
void clear_mesh(Mesh* mesh);
int make_vertex(Mesh* mesh, Vector vec);
Triangle* make_triangle(Mesh* mesh, Vertex* vt1, Vertex* vt2, Vertex* vt3, float max_len);
void delete_triangle(Triangle* tri, Mesh* mesh, int dverts);
void remove_tri_from_vert(Vertex* vert, Triangle* tri, int num, Mesh* mesh, int dverts);
void delete_vertex(Vertex* vert, Mesh* mesh);
void remove_unused_verts(Mesh* mesh);
int check_proposed_tri(Vertex* v1, Vertex* v2, Vertex* v3);
int check_proposed_edge(Vertex* v1, Vertex* v2);
void add_tri_to_vert(Vertex* vert, Triangle* tri);
int set_triangle_geometry(Triangle* tri);
int plane_thru_vectors(Vector v0, Vector v1, Vector v2, float* aa, float* bb, float* cc, float* dd);
int compute_edge_planes(Triangle* tri);
void find_vertex_normals(Mesh* mesh);
void find_vertex_normal(Vertex* vert);
void find_mesh_edges(Mesh* mesh);
int vertex_edge_test(Vertex* vert);
void clean_up_mesh(Scan* scan);

#endif
