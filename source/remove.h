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

#ifndef ZIPPER_REMOVE_H
#define ZIPPER_REMOVE_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Parameters
void update_eat_resolution();
void set_eat_near_dist_factor(float factor);
float get_eat_near_dist_factor();
void set_eat_near_cos(float cosine);
float get_eat_near_cos();
void set_eat_start_iters(int iters);
int get_eat_start_iters();
void set_eat_start_factor(float factor);
float get_eat_start_factor();

// Declarations
void do_it_all();
void eat_edge_proc();
void eat_edge_pair(Scan* sc1, Scan* sc2);
void init_eating(Scan* scan);
void done_eating(Scan* scan);
int eat_mesh_edges(Scan* sc1, Scan* sc2, int draw, int conf, int to_edge, float near_dist);
void mark_for_eating(Scan* sc1, Scan* sc2, int draw, int conf, int to_edge);
void zipper_proc();
void align_proc();
void gather_triangles(Scan* sc1, Scan* sc2);
void zipper_meshes(Scan* sc1, Scan* sc2);
void fill_in_holes(Scan* sc1, Scan* sc2);
int find_edge_orientation(Vertex* v1, Vertex* v2);
void move_vertices(Scan* source, Scan* dest);

#endif
