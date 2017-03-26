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

#ifndef ZIPPER_INTERSECT_H
#define ZIPPER_INTERSECT_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Declarations
void intersect_meshes(Scan* sc1, Scan* sc2);
void finish_intersect_meshes(Scan* sc1, Scan* sc2);
void mark_intersected_tris(Scan* sc1, Scan* sc2);
void intersect_edge_with_near_tris(Vertex* v1, Vertex* v2, Triangle* cut_tri, Scan* sc1, Scan* sc2);
void verts_near_vert(Mesh* mesh, Mesh* not_mesh, Vector pnt, Vector norm, float radius);
void new_tri_intersection(
    Vertex* v1, Vertex* v2, int share_count,
    Triangle* near_tri, Triangle* cut_tri, Vector pos,
    float t, int inward, float dot
);
void init_tri_intersection(Triangle* tri);
void add_intersect_points(Scan* sc1, Scan* sc2);
void perform_intersect_clipping(Scan* sc1, Scan* sc2);
int between_cuts(Triangle* tri, Clip_List* clist, int in_vert, int out_vert, Vertex* between_list[]);
int tri_cut_by_tri(Triangle* tri1, Triangle* tri2, Cut* cuts[], Triangle* tris[]);

#endif
