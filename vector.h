#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

struct vect {
    double x;
    double y;
    double z;
};

typedef struct vect vect;

vect vnew(double x, double y, double z);

vect vmult(vect a, double k);

vect vdiv(vect a, double k);

vect vneg(vect a);

vect vadd(vect a, vect b);

vect vsub(vect a, vect b);

double vdot(vect a, vect b);

vect vcross(vect a, vect b);

double vmag2(vect a);

double vmag(vect a);

vect vunit(vect a);

vect vmin(vect a, vect b);

vect vmax(vect a, vect b);

vect vinv(vect a);

double vcomp(vect a, int comp);

void print_vect(vect a);

#endif