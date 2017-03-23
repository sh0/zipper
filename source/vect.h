/*
 * vect:
 *  Functions to support operations on vectors and matrices.
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 */

#ifndef VECTDEF
#define VECTDEF

#include <math.h>
#include <matrix.h>

float* vnew();
float* vclone(const float*);
void vcopy(const float*, float*);
void vprint(const float*);
void vset(float*, float, float, float);
void vzero(float*);
void vnormal(float*);
float vlength(const float*);
void vscale(float*, float);
void vmult(const float*, const float*, float*);
void vadd(const float*, const float*, float*);
void vsub(const float*, const float*, float*);
void vhalf(const float*, const float*, float*);
float vdot(const float*, const float*);
void vcross(const float*, const float*, float*);
void vdirection(const float*, float*);
void vreflect(const float*, const float*, float*);
void vmultmatrix(const Matrix, const Matrix, Matrix);
void vtransform(const float*, const Matrix, float*);
void vtransform4(const float*, const Matrix, float*);

extern Matrix idmatrix;

void mcopy(const Matrix, Matrix);
void minvert(const Matrix, Matrix);
void vgetmatrix(Matrix m);
void linsolve(const float *[], int, float*);

#endif /* VECTDEF */
