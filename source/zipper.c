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
#include <fcntl.h>
#include <strings.h>
#include <malloc.h>
#include <sys/time.h>
#ifdef VOID
#undef VOID
#endif
#include "cyfile.h"
#include "zipper.h"
#include "zipglobal.h"

/* list of scans */
Scan* scans[SCAN_MAX];
int nscans = 0;

/* old zippering? */
int zipper_old = 0;

int global_dont_draw = 1;
int global_chew_count = 8;

extern int move_num;
extern Matrix rotmat;
extern Matrix transmat;


float MAX_EDGE_LENGTH_FACTOR;
float MAX_EDGE_LENGTH;

char* configuration_filename = NULL;


update_edge_length_resolution()
{
    MAX_EDGE_LENGTH = ZIPPER_RESOLUTION * MAX_EDGE_LENGTH_FACTOR;
}


set_max_edge_length_factor(factor)
float factor;
{
    MAX_EDGE_LENGTH_FACTOR = factor;
    MAX_EDGE_LENGTH = ZIPPER_RESOLUTION * MAX_EDGE_LENGTH_FACTOR;
}


float
get_max_edge_length_factor()
{
    return MAX_EDGE_LENGTH_FACTOR;
}



/******************************************************************************
Main routine.
******************************************************************************/

main(argc, argv)
int argc;
char* argv[];
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

int level_to_inc(level)
int level;
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

float edge_length_max(level)
int level;
{
    float max_length;
    int inc;

    /* pick how far apart the mesh samples are, based on the level of */
    /* detail requested */

    inc = level_to_inc(level);

    /* compute maximum okay length of a triangle edge */

    max_length = MAX_EDGE_LENGTH * inc;

    return (max_length);
}


/******************************************************************************
Pre-compute information about geometry to speed up mesh position calculations.

Entry:
  sc - geometry to augment
******************************************************************************/

setup_geometry_info(sc)
Scan* sc;
{
    int lg;
    float theta;
    float to_meters = 1.0e-6;
    GSPEC* gs = sc->gs;
    short right_hand = gs->flags & FLAG_THETARIGHT;

    /* pre-compute sines and cosines if this is rotational file */

    if (!(gs->flags & FLAG_CARTESIAN)) {

        /* allocate space for these tables */
        sc->sin_theta = (float*) myalloc(sizeof(float) * gs->nlg);
        sc->cos_theta = (float*) myalloc(sizeof(float) * gs->nlg);

        /* compute and store values */
        for (lg = 0; lg < gs->nlg; lg++) {
            theta = lg * gs->lgincr * to_meters;
            sc->sin_theta[lg] =  sin(theta) * to_meters;
            sc->cos_theta[lg] = -cos(theta) * to_meters;
            if (right_hand)
                sc->sin_theta[lg] = -sc->sin_theta[lg];
        }
    }

}


/******************************************************************************
Return the 3-space position of a mesh point in a depth file.  Distances
are measured in meters.

Entry:
  sc    - depth scan information
  lt,lg - lattitude and longitude of mesh point

Exit:
  vec - position of mesh point in 3-space
  returns 0 if everything okay, 1 if mesh point has no depth value
******************************************************************************/

int get_gs_coord(sc, lt, lg, vec)
Scan* sc;
int lt, lg;
Vector vec;
{
    GSPEC* gs = sc->gs;
    float to_meters = 1.0e-6;
    float xval, yval;
    float radius;

    if (gs->flags & FLAG_CARTESIAN) {

        xval = (lg - gs->nlg / 2) * gs->lgincr * to_meters;

        if (gs->flags & FLAG_BILATERAL)
            yval = (lt % (gs->nlt / 2) - gs->nlt) * gs->ltincr * to_meters;
        else
            yval = (lt - gs->nlt / 2) * gs->ltincr * to_meters;

        radius = GETR(gs, lt, lg);

        if (radius != CYVOID) {
            vec[X] = xval;
            vec[Y] = yval;
            vec[Z] = radius * to_meters;
        } else {
            vec[X] = 0;
            vec[Y] = 0;
            vec[Z] = 0;
            return (1);
        }

    } else {

        radius = GETR(gs, lt, lg);

        if (radius != CYVOID) {
            vec[X] = radius * sc->sin_theta[lg];
            vec[Y] = (lt - gs->nlt / 2) * gs->ltincr * to_meters;
            vec[Z] = radius * sc->cos_theta[lg];
        } else {
            vec[X] = 0;
            vec[Y] = 0;
            vec[Z] = 0;
            return (1);
        }

    }

    /* specify that we got a good value */
    return (0);
}

/******************************************************************************
Return the number of seconds since this routine was last called.
Returns junk on first call.
******************************************************************************/

float time_it()
{
    static long int last_sec;
    static long int last_usec;
    long int sec;
    long int usec;
    struct timeval tp;
    struct timezone tzp;
    float time;

    gettimeofday(&tp, &tzp);
    sec = tp.tv_sec - last_sec;
    usec = tp.tv_usec - last_usec;
    last_sec = tp.tv_sec;
    last_usec = tp.tv_usec;

    if (usec < 0)
        time = (sec - 1) + (usec + 1000000) / 1000000.0;
    else
        time = sec + usec / 1000000.0;

    return (time);
}


/******************************************************************************
Print the positions of all meshes.
******************************************************************************/

#if 0
void print_positions(fp)
FILE* fp;
{
    int i, j, k;
    Quaternion quat;
    Matrix mat;

    /* print viewing parameters */
    mat_to_quat(rotmat, quat);

    fprintf(fp, "camera %g %g %g  %g %g %g %g\n",
            transmat[3][0], transmat[3][1], transmat[3][2],
            quat[0], quat[1], quat[2], quat[3]);

    /* print quaternion version */

    for (i = 0; i < nscans; i++) {
        mat_to_quat(scans[i]->rotmat, quat);
        if (scans[i]->file_type == POLYFILE)
            fprintf(fp, "bpolygon ");
        else if (scans[i]->file_type == RAWFILE)
            fprintf(fp, "rmesh ");
        else if (scans[i]->file_type == PLYRANGEFILE)
            fprintf(fp, "bmesh ");
        else
            fprintf(fp, "mesh ");
        fprintf(fp, "%s %g %g %g %g %g %g %g\n", scans[i]->name,
                scans[i]->xtrans, scans[i]->ytrans, scans[i]->ztrans,
                quat[0], quat[1], quat[2], quat[3]);
    }
    fprintf(fp, "\n");
}
#endif

/******************************************************************************
Print the positions of all meshes in matrix form.
Modified by Afra Zomorodian to print out the number of scans
too.  7/28/95
******************************************************************************/

#if 0
void print_mat_positions(fp)
FILE* fp;
{
    int i, j, k;
    Quaternion quat;
    Matrix mat, imat;
    float det;

    /* print viewing parameters */
    mat_to_quat(rotmat, quat);
    fprintf(fp, "Number of meshes:  %d\n", nscans);
    fprintf(fp, "camera %g %g %g  %g %g %g %g\n",
            transmat[3][0], transmat[3][1], transmat[3][2],
            quat[0], quat[1], quat[2], quat[3]);

    /* print matrix version */

    fprintf(fp, "\n");
    for (i = 0; i < nscans; i++) {
        fprintf(fp, "%s:\n", scans[i]->name);

        /* compute the tranformation matrix */
        mat_translate(mat, scans[i]->xtrans, scans[i]->ytrans, scans[i]->ztrans);
        mat_mult(mat, mat, scans[i]->rotmat);

        /* inverse of mat */
        mat_copy(imat, mat);
        det = mat_invert(imat);

        fprintf(fp, "matrix:\n");
        for (j = 0; j <= 3; j++) {
            for (k = 0; k <= 3; k++)
                fprintf(fp, "%f  ", mat[k][j]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");

        fprintf(fp, "inverse matrix:\n");
        for (j = 0; j <= 3; j++) {
            for (k = 0; k <= 3; k++)
                fprintf(fp, "%f  ", imat[k][j]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}
#endif

/******************************************************************************
Create all the meshes for the current level of detail.
******************************************************************************/

create_current_level()
{
    int i;

    for (i = 0; i < nscans; i++)
        create_scan_mesh(scans[i], mesh_level);
}

