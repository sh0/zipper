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

#ifndef ZIPPER_FILL_H
#define ZIPPER_FILL_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Vertex that is part of a hole that is being filled
typedef struct FillVertex {
    Vertex* vert; // vertex helping cover the hole
    unsigned char on_edge; // on the edge of the hole?
    Vertex* e1, *e2; // adjacent vertices on edge, if any
} FillVertex;

// Triangle used to fill a hole
typedef struct FillTri {
    Triangle* tri; // triangle used to fill hole
    int index; // index into list containing triangle
} FillTri;

// Parameters
void update_fill_resolution();
void set_fill_edge_length_factor(float factor);
float get_fill_edge_length_factor();

// Declarations
void better_fill_loop(int loop, Scan* scan);
int fix_fill_size(int loop, Scan* scan, float max_len, float* size);
int maybe_split_edge(Triangle* tri, Scan* scan, int idx1, int idx2);
void init_fill_lists();
void new_ftri(Triangle* tri);
void delete_ftri(FillTri* ftri);
void new_fvert(Vertex* vert, int on_edge, Vertex* v1, Vertex* v2);
void fill_tri_split2(FillTri* ftri, Mesh* mesh, int index);
void fill_tri_split3(FillTri* ftri, Mesh* mesh, int index);
void fill_tri_split4(FillTri* ftri, Mesh* mesh, int index);
void smooth_hole_vertices(Scan* scan);
void swap_hole_edges(Scan* sc);

#endif
