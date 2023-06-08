#ifndef COLOR_H
#define COLOR_H

#include <stdio.h>

struct color {
    double r;
    double g;
    double b;
};

typedef struct color color;

color cnew(double r, double g, double b);

color cmult(color a, double k);

color cdiv(color a, double k);

color cneg(color a);

color cadd(color a, color b);

color csub(color a, color b);

void printColor(color a);

#endif