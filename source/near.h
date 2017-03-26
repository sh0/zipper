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

#ifndef ZIPPER_NEAR_H
#define ZIPPER_NEAR_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Declarations
void init_table(Mesh* mesh, float size);
void add_to_hash(Vertex* vert, Mesh* mesh);
void remove_from_hash(Vertex* vert, Mesh* mesh);
Vertex* find_nearest(Mesh* mesh, Mesh* not_mesh, Vector pnt, Vector norm, float min_dot);
Vertex* large_find_nearest(Mesh* mesh, Mesh* not_mesh, Vector pnt, Vector norm, float min_dot);
Vertex* new_find_nearest(Mesh* mesh, Mesh* not_mesh, Vector pnt, Vector norm, float max, float min_dot);
int nearest_on_mesh(
    Scan* sc, Mesh* mesh, Mesh* not_mesh,
    Vector pos, Vector norm, float max, float min_dot,
    NearPosition* near_info
);
float nearest_on_edges(
    Vector pos, Vertex* near, float max, Vector near_pos,
    int* on_edge, Vertex** near2, float* conf, float* ival
);
float nearest_on_tris(
    Vector pos, Vertex* near, float max, Vector near_pos,
    Triangle** near_tri, float* conf, Vector barycentric
);
int is_near_edge(NearPosition* near);
int is_vertex_near_edge(Vertex* v);

#endif
