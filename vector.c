#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include "vector.h"
#include "misc.h"

vect vnew(double x, double y, double z) {
    vect v = {.x = x, .y = y, .z = z};
	return v;
}

vect vmult(vect a, double k) {
    return vnew(a.x * k, a.y * k, a.z * k);
}

vect vdiv(vect a, double k) {
    assert(k != 0);
    return vmult(a, (double) 1 / k);
}

vect vneg(vect a) {
    return vnew(-a.x, -a.y, -a.z);
}

vect vadd(vect a, vect b) {
    return vnew(a.x + b.x, a.y + b.y, a.z + b.z);
}

vect vsub(vect a, vect b) {
    return vadd(a, vneg(b));
}

double vdot(vect a, vect b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vect vcross(vect a, vect b) {
    return vnew(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

double vmag2(vect a) {
    return vdot(a, a);
}

double vmag(vect a) {
    return sqrtl(vmag2(a));
}

vect vunit(vect a) {
    return vdiv(a, vmag(a));
}

vect vmin(vect a, vect b) {
    return vnew(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

vect vmax(vect a, vect b) {
    return vnew(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}

vect vinv(vect a) {
    assert(a.x != 0 || a.y != 0 || a.z != 0);
    return vnew(1 / a.x, 1 / a.y, 1 / a.z);
}

double vcomp(vect a, int comp) {
    if (comp == 0) return a.x;
    if (comp == 1) return a.y;
    if (comp == 2) return a.z;
    assert(false);
}

void print_vect(vect a) {
    printf("{%.6f, %.6f, %.6f}\n", a.x, a.y, a.z);
}