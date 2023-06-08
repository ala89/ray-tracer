#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "ray.h"
#include "aabb.h"

struct triangle {
    int v[3];
    int n[3];
    vect* vertices;
};

typedef struct triangle triangle;

void triangleRayIntersection(triangle tri, Ray* r);

AABB triangleChoppedBounds(triangle tri, int dim, double leftPlanePos, double rightPlanePos);

void triangleSplitBounds(triangle tri, int dim, double planePos, AABB* leftBounds, AABB* rightBounds);

AABB triangleBounds(triangle t);

#endif