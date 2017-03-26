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

// External
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Internal
#include "zipper.h"
#include "clip.h"
#include "consensus.h"
#include "fill.h"
#include "remove.h"
#include "mesh.h"
#include "ply_wrapper.h"

// Globals
#define SCAN_MAX 200
Scan* scans[SCAN_MAX];
int nscans = 0;
float ZIPPER_RESOLUTION = 0.0005; /* The "scale" of the system; formally SPACING. */
int mesh_level = 3; /* mesh display level */

// Parameters
static float MAX_EDGE_LENGTH_FACTOR;
static float MAX_EDGE_LENGTH;

void update_edge_length_resolution()
{
    MAX_EDGE_LENGTH = ZIPPER_RESOLUTION * MAX_EDGE_LENGTH_FACTOR;
}

void set_max_edge_length_factor(float factor)
{
    MAX_EDGE_LENGTH_FACTOR = factor;
    MAX_EDGE_LENGTH = ZIPPER_RESOLUTION * MAX_EDGE_LENGTH_FACTOR;
}

float get_max_edge_length_factor()
{
    return MAX_EDGE_LENGTH_FACTOR;
}

float get_zipper_resolution()
{
    return ZIPPER_RESOLUTION;
}

void set_zipper_resolution(float res)
{
    ZIPPER_RESOLUTION = res;

    update_edge_length_resolution();
    update_fill_resolution();
    update_eat_resolution();
    update_clip_resolution();
    update_consensus_resolution();
}

/******************************************************************************
Main routine.
******************************************************************************/
static void init_parameters()
{
    set_max_edge_length_factor(4.0);

    set_fill_edge_length_factor(2.0);

    set_conf_edge_count_factor(1.0);
    set_conf_edge_zero(0);
    set_conf_angle(0);
    set_conf_exponent(1.0);

    set_eat_near_dist_factor(2.0);
    set_eat_near_cos(-0.5);
    set_eat_start_iters(2);
    set_eat_start_factor(4.0);

    set_clip_near_dist_factor(2.0);
    set_clip_near_cos(0.3);
    set_clip_boundary_dist_factor(4.0);
    set_clip_boundary_cos(0.3);

    set_consensus_position_dist_factor(1.0);
    set_consensus_normal_dist_factor(3.0);
    set_consensus_jitter_dist_factor(0.01);

    set_range_data_sigma_factor(4.0);
    set_range_data_min_intensity(0.05);
    set_range_data_horizontal_erode(1);
}

int main(int argc, char* argv[])
{
    // Setup
    init_parameters();

    // Help
    if (argc < 4) {
        printf("Usage: zipper src1.ply src2.ply dst.ply\n");
        return 0;
    }

    // Read input
    if (read_ply(argv[1]) != 0) {
        printf("Failed to read input 1: %s\n", argv[1]);
        return 1;
    }
    if (read_ply(argv[2]) != 0) {
        printf("Failed to read input 2: %s\n", argv[2]);
        return 1;
    }

    // Process
    do_it_all();
    //clip_triangles(scans[0], scans[1]);

    // Write output
    write_ply(scans[0], argv[3], 1);
}

/******************************************************************************
Return how many range image positions are between each vertex at a given
mesh level.  E.g., mesh level 3 uses every 8th range image point.

Entry:
  level - level to find out about

Exit:
  returns number of range positions between vertices
******************************************************************************/
int level_to_inc(int level)
{
    switch (level) {
        case 0:
            return (1);
        case 1:
            return (2);
        case 2:
            return (4);
        case 3:
            return (8);
        default:
            fprintf(stderr, "level_to_inc: bad switch %d\n", level);
            exit(-1);
    }
}

/******************************************************************************
Return maximum length allowed for a triangle of a given level.
******************************************************************************/
float edge_length_max(int level)
{
    float max_length;
    int inc;

    /* pick how far apart the mesh samples are, based on the level of */
    /* detail requested */
    inc = level_to_inc(level);

    /* compute maximum okay length of a triangle edge */
    max_length = MAX_EDGE_LENGTH * inc;
    return max_length;
}

/******************************************************************************
Create all the meshes for the current level of detail.
******************************************************************************/
void create_current_level()
{
    int i;
    for (i = 0; i < nscans; i++)
        create_scan_mesh(scans[i], mesh_level);
}
