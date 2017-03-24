/*

Zipper together polygon meshes derived from multiple depth scans.

Greg Turk, December 1992

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
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "zipper.h"

int global_chew_count = 8;

// Globals
#define SCAN_MAX 200
Scan* scans[SCAN_MAX];
int nscans = 0;
float ZIPPER_RESOLUTION = 0.0005; /* The "scale" of the system; formally SPACING. */
int mesh_level = 3; /* mesh display level */

// Parameters
float MAX_EDGE_LENGTH_FACTOR;
float MAX_EDGE_LENGTH;

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

void init_resolution_parameters()
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

    /* These don't really belong under "resolution", but... */
    set_range_data_sigma_factor(4.0);
    set_range_data_min_intensity(0.05);
    set_range_data_horizontal_erode(1);
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
int main(int argc, char* argv[])
{
    init_resolution_parameters();

    if (argc < 4) {
        printf("Usage: zipper src1.ply src2.ply dst1.ply dst2.ply\n");
        return 0;
    }

    if (read_ply(argv[1]) != 0) {
        printf("Failed to read input 1: %s", argv[1]);
        return 1;
    }
    if (read_ply(argv[2]) != 0) {
        printf("Failed to read input 2: %s", argv[2]);
        return 1;
    }

    do_it_all();
    //clip_triangles(scans[0], scans[1]);

    write_ply(scans[0], argv[3], 1);
    write_ply(scans[1], argv[4], 1);
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
