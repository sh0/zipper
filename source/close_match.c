/*

Match up meshes by translation and rotation.

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "malloc.h"
#include "cyfile.h"
#include "zipper.h"

#if 0
extern double G_Tolerance;
#endif

#define  MATCH_MAX  20
extern Scan* match_from[];
extern Scan* match_to[];
extern Scan* match_drag[MATCH_MAX][SCAN_MAX];
extern int num_drag[];
extern int num_matches;
extern int draw_during_ops;
extern int processes_forked;

typedef float StateVector[7];   /* quaternion & translation */

float measure_error();

int align_on = 0;


/******************************************************************************
Match up meshes by translation and rotation.
******************************************************************************/

void closest_proc()
{
    int i;

    if (num_matches == 0) {
        six_degree_match(&scans[0], &scans[1], NULL, 1, 1, 0, 10);
    } else
        for (i = 0; i < num_matches; i++) {
            six_degree_match(&match_from[i], &match_to[i], match_drag[i],
                             1, 1, num_drag[i], 10);
        }
}


/******************************************************************************
Six degree of freedom match.

Entry:
  move_list - list meshes to move
  match_to  - list of meshes to match to
  drag_with - list of other scans to move along with move_list
  nmove     - number of scans to move
  nmatch    - number of scans to match
  ndrag     - number in drag list
  num_iters - number of alignment iterations to perform
******************************************************************************/

six_degree_match(move_list, match_to, drag_with, nmove, nmatch, ndrag, num_iters)
Scan** move_list;
Scan** match_to;
Scan** drag_with;
int nmove;
int nmatch;
int ndrag;
int num_iters;
{
    int i, j, k;
    Matrix rotmat;
    Matrix rotmat2;
    Matrix mat;
    Scan* sc1, *sc2;
    Scan* drag;
    Quaternion quat_rot;
    Vector translation;
    float time;
    extern float time_it();
    extern float error_measure();
    float error, new_error;
    StateVector svec, svec2;
    static StateVector move_hist[3];
    static float err_hist[3];
    int pipeline;
    int jump_step;
    int new;

    time = time_it();

    sc1 = match_to[0];
    sc2 = move_list[0];

    /* perform the matching operation iteratively */

    pipeline = 0;
    jump_step = 0;

    /* turn on flag saying we're aligning */
    align_on = 1;

    /* iterate the alignment process */

    for (i = 0; i < num_iters; i++) {

        /* find a transformation that minimizes distance between pairs of */
        /* corresponding points on two meshes */

        minimize_pairs_dist(sc1, sc2, svec);

        /* translate the state vector from the match into a matrix */
        sv_to_matrix(svec, rotmat);

        /* check the error */
        error = measure_error(rotmat, sc2);

#if 1
        printf("avg dist = %.6f mm\n", error * 1000);
#endif

        if (i == 0) {
            sv_copy(svec, move_hist[2]);
            err_hist[2] = error;
            set_scan_matrix(sc2, rotmat);
            pipeline = 1;
            continue;
        }

        /* shift the old state vectors down the list and save the new one */
        /* at the top of the list */
        sv_copy(move_hist[1], move_hist[0]);
        sv_copy(move_hist[2], move_hist[1]);
        sv_copy(svec, move_hist[2]);

        err_hist[0] = err_hist[1];
        err_hist[1] = err_hist[2];
        err_hist[2] = error;

        pipeline++;

        /* perhaps compute a better move */
        if (pipeline >= 3 &&
            err_hist[0] > err_hist[1] && err_hist[1] > err_hist[2]) {
            new = improve_move(move_hist[0], move_hist[1], move_hist[2],
                               err_hist[0], err_hist[1], err_hist[2], svec2);
            if (new) {
                sv_to_matrix(svec2, rotmat);
                set_scan_matrix(sc2, rotmat);
                minimize_pairs_dist(sc1, sc2, svec2);
                sv_to_matrix(svec2, rotmat);
                new_error = measure_error(rotmat, sc2);
#if 0
                printf("new avg dist = %.6f mm\n", new_error * 1000);
#endif
                if (new_error < error) {
#if 0
                    printf("(taken)\n");
#endif

#if 1
                    printf("avg dist = %.6f mm\n", new_error * 1000);
#endif
                    sv_copy(move_hist[1], move_hist[0]);
                    sv_copy(move_hist[2], move_hist[1]);
                    sv_copy(svec2, move_hist[2]);
                    sv_copy(svec2, svec);
                    err_hist[0] = err_hist[1];
                    err_hist[1] = err_hist[2];
                    err_hist[2] = new_error;
                    pipeline++;
                } else {
#if 0
                    printf("(not taken)\n");
#endif
                    /* restore old matrix */
                    sv_to_matrix(svec, rotmat);
                    set_scan_matrix(sc2, rotmat);
                }
            }
        }

        /* use the current state vector to make the correct transformation matrix */
        sv_to_matrix(svec, rotmat);

        /* set scan 2's matrix */
        set_scan_matrix(sc2, rotmat);

#if 0

        /*** how do I do this now??? ***/
        /*** how do I do this now??? ***/
        /*** how do I do this now??? ***/
        /*** how do I do this now??? ***/

        /* drag other meshes along */
        for (j = 0; j < ndrag; j++) {

            drag = drag_with[j];

            mat_copy(mat, drag->rotmat);
            mat[3][0] = drag->xtrans;
            mat[3][1] = drag->ytrans;
            mat[3][2] = drag->ztrans;

            mat_mult(drag->rotmat, rotmat, mat);

            drag->xtrans = drag->rotmat[3][0];
            drag->ytrans = drag->rotmat[3][1];
            drag->ztrans = drag->rotmat[3][2];
            drag->rotmat[3][0] = 0;
            drag->rotmat[3][1] = 0;
            drag->rotmat[3][2] = 0;
        }

#endif

        /* check to see if we've been told to stop alignment */
        if (!align_on) {
            fprintf(stderr, "Alignment stopped.\n");
            break;
        }
    }
}


/******************************************************************************
Set the transformation of a given scan.

Entry:
  scan - scan to set the matrix of
  mat  - matrix to use
******************************************************************************/

set_scan_matrix(scan, mat)
Scan* scan;
Matrix mat;
{
    /* copy "mat" into scan */
    mat_copy(scan->rotmat, mat);

    /* extract the translation from the resulting position */
    scan->xtrans = scan->rotmat[3][0];
    scan->ytrans = scan->rotmat[3][1];
    scan->ztrans = scan->rotmat[3][2];

    scan->rotmat[3][0] = 0;
    scan->rotmat[3][1] = 0;
    scan->rotmat[3][2] = 0;
}


/******************************************************************************
Find the transformation matrix that minimizes the least squared distance
between pairs of corresponding points.

Entry:
  sc1,sc2 - meshes to match

Exit:
  svec - transformation state that gives best rotation and translation
******************************************************************************/

minimize_pairs_dist(sc1, sc2, svec)
Scan* sc1, *sc2;
StateVector svec;
{
    int j;
    Quaternion quat_rot;
    Vector translation;
    Vector p1, p2;
    Matrix mat;

    /* find list of corresponding points on the two meshes */
    create_match_list(sc1, sc2, 0, 0);

    /* give these pairs of points to the rotational matching routine */

    init_pairs();

    for (j = 0; j < global_num_matches; j++) {

        p1[X] = pos_matches[j]->pos[X];
        p1[Y] = pos_matches[j]->pos[Y];
        p1[Z] = pos_matches[j]->pos[Z];

        p2[X] = p1[X] + pos_matches[j]->dir[X];
        p2[Y] = p1[Y] + pos_matches[j]->dir[Y];
        p2[Z] = p1[Z] + pos_matches[j]->dir[Z];

        world_to_mesh(sc2, p2, p2);

        add_pair(p1, p2, pos_matches[j]->confidence);
    }

    /* find the best rotational and translational match */

    match_pairs(mat, quat_rot, translation);
    svec[0] = quat_rot[0];
    svec[1] = quat_rot[1];
    svec[2] = quat_rot[2];
    svec[3] = quat_rot[3];
    svec[4] = translation[0];
    svec[5] = translation[1];
    svec[6] = translation[2];
}


/******************************************************************************
Copy a seven-component state vector.
******************************************************************************/

sv_copy(src, dst)
StateVector src, dst;
{
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
    dst[4] = src[4];
    dst[5] = src[5];
    dst[6] = src[6];
}


/******************************************************************************
Convert a seven-element state vector into a matrix.
******************************************************************************/

sv_to_matrix(state, mat)
StateVector state;
Matrix mat;
{
    Quaternion quat;

    quat[0] = state[0];
    quat[1] = state[1];
    quat[2] = state[2];
    quat[3] = state[3];

    convert_quat_to_mat(quat, mat);

    mat[3][0] = state[4];
    mat[3][1] = state[5];
    mat[3][2] = state[6];
}


/******************************************************************************
Return the mean error between point pairs aligned by matrix "mat".
******************************************************************************/

float measure_error(mat, scan)
Matrix mat;
Scan* scan;
{
    int i, j, count;
    Vector pos, diff;
    float dist;
    float sum;
    float weight_sum, error, unweighted_sum;

    sum = 0;
    weight_sum = 0;
    unweighted_sum = 0;
    count = 0;
    for (i = 0; i < global_num_matches; i++) {

        /*
        printf ("%f %f %f  %f\n", pos_matches[i]->pos[X],
                                  pos_matches[i]->pos[Y],
                                  pos_matches[i]->pos[Z],
                                  pos_matches[i]->confidence);
        */

        vadd(pos_matches[i]->pos, pos_matches[i]->dir, pos);
        world_to_mesh(scan, pos, pos);
        vapply(mat, pos, pos);
        vsub(pos, pos_matches[i]->pos, diff);
        dist = diff[X] * diff[X] + diff[Y] * diff[Y] + diff[Z] * diff[Z];
        sum += dist * pos_matches[i]->confidence;
        weight_sum += pos_matches[i]->confidence;
        unweighted_sum += dist;
        count++;
    }

    /*
    return (sum / weight_sum);
    */

    /*
      error = sqrt(unweighted_sum/count);
    */

    error = sqrt(sum / weight_sum);

    return (error);
}


/******************************************************************************
Try to improve the new state vector by looking at a recent history.

Entry:
  s1,s2,s3 - old state vectors, with s3 being newest
  e1,e2,e3 - associated error measures

Exit:
  result - (possibly improved) state vector to use
  returns 1 if the routine produced an "improved" state vector
******************************************************************************/

int improve_move(s1, s2, s3, e1, e2, e3, result)
StateVector s1, s2, s3;
float e1, e2, e3;
StateVector result;
{
    int i, j;
    StateVector diff1, diff2;
    float dlen1, dlen2;
    float dot;
    float theta;
    static float old_theta = -180;
    float x1, x2, x3;
    float a1, b1;     /* coefficients for linear fit */
    float a2, b2, c2; /* coefficients for quadratic fit */
    float x[3], y[3];
    float max_step = 25;
    float v1, v2;
    float vmax;
    float k;

#define SMALL_THETA 10.0

    /* copy the given state vector to the result, making it the default */
    for (i = 0; i < 7; i++)
        result[i] = s3[i];

    /* compute difference vectors and their lengths */

    dot = 0;
    dlen1 = dlen2 = 0;
    for (i = 0; i < 7; i++) {
        diff1[i] = s2[i] - s1[i];
        diff2[i] = s3[i] - s2[i];
        dlen1 += diff1[i] * diff1[i];
        dlen2 += diff2[i] * diff2[i];
        dot += diff1[i] * diff2[i];
    }
    dlen1 = sqrt(dlen1);
    dlen2 = sqrt(dlen2);

    theta = acos(dot / (dlen1 * dlen2)) * 180.0 / M_PI;

    if (old_theta > SMALL_THETA) {
        old_theta = theta;
        return (0);
    }

    old_theta = theta;

    if (theta > SMALL_THETA) {
        return (0);
    }

#if 1

    /* if we get here, both old and new angle are small, so we can */
    /* try to improve the state vector */

    vmax = max_step * dlen2;

    if (vmax < 0) {
        fprintf(stderr, "improve_move: vmax less than zero: %f\n", vmax);
        return (0);
    }

    x3 = 0;
    x2 = -dlen2;
    x1 = -dlen1 + x2;

#if 0
    printf("x1 x2 x3: %f %f %f\n", x1, x2, x3);
    printf("e1 e2 e3: %.9f %.9f %.9f\n", e1, e2, e3);
#endif

    /* compute linear fit */
    x[0] = x1;
    x[1] = x2;
    x[2] = x3;
    y[0] = e1;
    y[1] = e2;
    y[2] = e3;
    linear_solve(x, y, 3, &a1, &b1);
    v1 = -b1 / a1;

    /* compute quadratic fit */
    qsolve(x1, x2, x3, e1, e2, e3, &a2, &b2, &c2);
    v2 = -b2 / (2 * a2);

#if 0
    printf("linear:    a1 = %f, b1 = %f\n", a1, b1);
    printf("quadratic: a2 = %f, b2 = %f, c2 = %f\n", a2, b2, c2);
    printf("v1 = %f, v2 = %f, vmax = %f\n", v1, v2, vmax);
#endif

    if (v1 > vmax && v2 > vmax) {
#if 0
        printf("maximum case\n");
#endif
        k = vmax / dlen2;
        for (i = 0; i < 7; i++)
            result[i] = s3[i] + k * diff2[i];
        return (1);
    } else if (0 < v2 && v2 < vmax && v1 > v2) {
#if 0
        printf("quadratic case\n");
#endif
        k = v2 / dlen2;
        for (i = 0; i < 7; i++)
            result[i] = s3[i] + k * diff2[i];
        return (1);
    } else if ((0 < v1 && v1 < vmax && v2 > v1) ||
               (v2 < 0 && 0 < v1 && v1 < vmax)) {
#if 0
        printf("linear case\n");
#endif
        k = v1 / dlen2;
        for (i = 0; i < 7; i++)
            result[i] = s3[i] + k * diff2[i];
        return (1);
    } else {
#if 0
        printf("no action taken\n");
#endif
    }

#endif

    /* if we get here then no improvement was suggested */
    return (0);
}


/******************************************************************************
Linear fit to (x,y) pairs of data points.  This is a hacked version of the
Numerical Recipes routine "fit", only taking the code that doesn't assume
you have standard deviations.

Entry:
  x,y   - pairs of data values (x,y)
  ndata - number of pairs

Exit:
  a,b - coefficients of fitting a line: y = ax + b (Note order of a and b!)
******************************************************************************/

linear_solve(x, y, ndata, a, b)
float x[], y[];
int ndata;
float* a, *b;
{
    int i;
    float sx;
    float sy;
    float t;
    float st2;
    float sxoss;

    *a = 0;

    /* accumulate sums without weights */
    sx = sy = 0;
    for (i = 0; i < ndata; i++) {
        sx += x[i];
        sy += y[i];
    }
    sxoss = sx / (float) ndata;

    st2 = 0;
    for (i = 0; i < ndata; i++) {
        t = x[i] - sxoss;
        st2 += t * t;
        *a += t * y[i];
    }

    *a /= st2;
    *b = (sy - sx * (*a)) / (float) ndata;
}


/******************************************************************************
Find coefficients for a quadratic equation that match three given data points.

Entry:
  x1,x2,x3 - location of data points
  y1,y2,y3 - values at data points

Exit:
  a,b,c - coefficients of fitting quadratic: y = ax^2 + by + c
******************************************************************************/

qsolve(x1, x2, x3, y1, y2, y3, a, b, c)
float x1, x2, x3;
float y1, y2, y3;
float* a, *b, *c;
{
    Vector y;
    Matrix mat, imat;
    float det;

    vset(y, y1, y2, y3);

    mat[0][0] = x1 * x1;
    mat[0][1] = x2 * x2;
    mat[0][2] = x3 * x3;
    mat[0][3] = 0;

    mat[1][0] = x1;
    mat[1][1] = x2;
    mat[1][2] = x3;
    mat[1][3] = 0;

    mat[2][0] = 1;
    mat[2][1] = 1;
    mat[2][2] = 1;
    mat[2][3] = 0;

    mat[3][0] = 0;
    mat[3][1] = 0;
    mat[3][2] = 0;
    mat[3][3] = 1;

    mat_copy(imat, mat);
    det = mat_invert(imat);

#if 0
    printf("mat:\n");
    mat_print(mat);

    printf("imat:\n");
    mat_print(imat);

    printf("det = %f\n", det);
#endif

    *a = imat[0][0] * y1 + imat[1][0] * y2 + imat[2][0] * y3;
    *b = imat[0][1] * y1 + imat[1][1] * y2 + imat[2][1] * y3;
    *c = imat[0][2] * y1 + imat[1][2] * y2 + imat[2][2] * y3;
}

