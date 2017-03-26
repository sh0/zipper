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

#ifndef ZIPPER_CONSENSUS_H
#define ZIPPER_CONSENSUS_H

// Internal
#include "zipper.h"
#include "matrix.h"

// Vertex to help find consensus
typedef struct Cvert {
    Vector coord;
    struct Ctri** tris;
    unsigned char ntris;
    struct Cvert* next;
} Cvert;

// Triangle to help find consensus
typedef struct Ctri {
    Cvert* verts[3];
} Ctri;

typedef struct Cinfo {
    Vector pos; // consensus position
    Vector normal; // consensus surface normal (in global coords)
    float intensity; // consensus intensity
    float weights;
    float red, grn, blu; // color
    int count;
} Cinfo;

// Parameters
void update_consensus_resolution();
void set_consensus_position_dist_factor(float factor);
float get_consensus_position_dist_factor();
void set_consensus_normal_dist_factor(float factor);
float get_consensus_normal_dist_factor();
void set_consensus_jitter_dist_factor(float factor);
float get_consensus_jitter_dist_factor();

// Declarations
void consensus_surface(Scan* scan, int level, float k_scale);
void new_consensus_surface(Scan* con_scan, Scan** scan_list, int* read_list, int num_scans, int level);
void marc_find_average_positions(Scan* scan, int level, float k_scale);
void new_find_average_positions(Scan* con_scan, Scan** scan_list, int* use_old_mesh, int num_scans, int level);
void intersect_segment_with_mesh(Vertex* v, Mesh* mesh, Scan* vscan, Scan* mscan, float search_dist, int mesh_index);
void find_average_positions(Scan* scan, int level, float k_scale);
int count_near_vert();
Vertex* found_vert_near_vert(int n);
void verts_near_pos(Mesh* mesh, Vector pnt, Vector norm, float radius);

#endif
