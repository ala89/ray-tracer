#ifndef AABB_H
#define AABB_H

#include <float.h>

#include "vector.h"
#include "ray.h"

struct AABB {
    vect pMin;
    vect pMax;
};

typedef struct AABB AABB;

AABB aabbNewEmpty();

AABB aabbNewPoint(vect p) ;

AABB aabbFromPoints(vect* points, int n);

AABB aabbExpand(AABB b, vect p);

AABB aabbUnion(AABB b1, AABB b2);

AABB aabbInter(AABB b1, AABB b2);

vect aabbDiag(AABB b);

double aabbSA(AABB b);

vect aabbCentroid(AABB b);

bool aabbRayIntersects(AABB b, Ray r, vect invDir);

bool aabbIsPointInside(vect p, AABB a);

bool aabbIsWithin(AABB a, AABB b);

#endif