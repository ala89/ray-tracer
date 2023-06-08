#include <stdlib.h>
#include <assert.h>

#include "color.h"

color cnew(double r, double g, double b) {
    color c = {.r = r, .g = g, .b = b};
	return c;
}

color cmult(color a, double k) {
    return cnew(a.r * k, a.g * k, a.b * k);
}

color cdiv(color a, double k) {
    assert(k != 0);
    return cmult(a, (double) 1 / k);
}

color cneg(color a) {
    return cnew(-a.r, -a.g, -a.b);
}

color cadd(color a, color b) {
    return cnew(a.r + b.r, a.g + b.g, a.b + b.b);
}

color csub(color a, color b) {
    return cadd(a, cneg(b));
}

void printColor(color a) {
    printf("{%.3f, %.3f, %.3f}\n", a.r, a.g, a.b);
}