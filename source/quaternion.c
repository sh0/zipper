/*

Manipulation of quaternions.

Some of these routines were taken from:

  Ken Shoemake, "Quaternion Calculus for Animation", published where?

*/

#include <matrix.h>
#include <math.h>


/******************************************************************************
Convert a quaternion into a rotation matrix.  Does *not* assume a unit
quaternion.  From Ken Shoemake.

Entry:
  q - quaternion

Exit:
  mat - rotation matrix
******************************************************************************/

quat_to_mat(q, mat)
Quaternion q;
Matrix mat;
{
    float s;
    float xs, ys, zs;
    float wx, wy, wz;
    float xx, xy, xz;
    float yy, yz, zz;

    /* for unit q, just set s = 2 or set xs = q[X] + q[X], etc. */

    s = 2 / (q[X] * q[X] + q[Y] * q[Y] + q[Z] * q[Z] + q[W] * q[W]);

    xs = q[X] * s;
    ys = q[Y] * s;
    zs = q[Z] * s;

    wx = q[W] * xs;
    wy = q[W] * ys;
    wz = q[W] * zs;

    xx = q[X] * xs;
    xy = q[X] * ys;
    xz = q[X] * zs;

    yy = q[Y] * ys;
    yz = q[Y] * zs;
    zz = q[Z] * zs;

    mat[X][X] = 1 - (yy + zz);
    mat[X][Y] = xy - wz;
    mat[X][Z] = xz + wy;
    mat[X][W] = 0;

    mat[Y][X] = xy + wz;
    mat[Y][Y] = 1 - (xx + zz);
    mat[Y][Z] = yz - wx;
    mat[Y][W] = 0;

    mat[Z][X] = xz - wy;
    mat[Z][Y] = yz + wx;
    mat[Z][Z] = 1 - (xx + yy);
    mat[Z][W] = 0;

    mat[W][X] = 0;
    mat[W][Y] = 0;
    mat[W][Z] = 0;
    mat[W][W] = 1;
}


/******************************************************************************
Convert a rotation matrix into a unit quaternion.  From Ken Shoemake.

Entry:
  mat - rotation matrix

Exit:
  q - quaternion
******************************************************************************/

mat_to_quat(mat, q)
Matrix mat;
Quaternion q;
{
    int i, j, k;
    float tr, s;
    static int nxt[3] = {Y, Z, X};  /* hey, this is a new one on me! */

    tr = mat[X][X] + mat[Y][Y] + mat[Z][Z];

    if (tr > 0) {
        s = sqrt(tr + 1);
        q[W] = s * 0.5;
        s = 0.5 / s;
        q[X] = (mat[Z][Y] - mat[Y][Z]) * s;
        q[Y] = (mat[X][Z] - mat[Z][X]) * s;
        q[Z] = (mat[Y][X] - mat[X][Y]) * s;
    } else {
        i = X;
        if (mat[Y][Y] > mat[X][X])
            i = Y;
        if (mat[Z][Z] > mat[i][i])
            i = Z;
        j = nxt[i];
        k = nxt[j];
        s = sqrt(1 + (mat[i][i] - (mat[j][j] + mat[k][k])));
        q[i] = s * 0.5;
        s = 0.5 / s;
        q[W] = (mat[k][j] - mat[j][k]) * s;
        q[j] = (mat[j][i] + mat[i][j]) * s;
        q[k] = (mat[k][i] + mat[i][k]) * s;
    }
}

