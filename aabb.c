#include <float.h>
#include <stdbool.h>

#include "vector.h"
#include "aabb.h"
#include "misc.h"
#include "ray.h"

AABB aabbNewEmpty() {
    AABB b;
    b.pMin = vnew(INFINITY, INFINITY, INFINITY);
    b.pMax = vnew(-INFINITY, -INFINITY, -INFINITY);
    return b;
}

AABB aabbNewPoint(vect p) {
    AABB b = {.pMin = p, .pMax = p};
    return b;
}

AABB aabbNewPoints(vect pMin, vect pMax) {
    AABB b = {.pMin = pMin, .pMax = pMax};
    return b;
}

AABB aabbFromPoints(vect* points, int n) {
    AABB b = aabbNewEmpty();
    for (int i = 0; i < n; i++) {
        b = aabbExpand(b, points[i]);
    }
    return b;
}

bool aabbIsEmpty(AABB b) {
    return b.pMin.x > b.pMax.x || b.pMin.y > b.pMax.y || b.pMin.z > b.pMax.z;
}

AABB aabbExpand(AABB b, vect p) {
   return aabbNewPoints(vmin(b.pMin, p), vmax(b.pMax, p));
}

AABB aabbUnion(AABB b1, AABB b2) {
    return aabbNewPoints(vmin(b1.pMin, b2.pMin), vmax(b1.pMax, b2.pMax));
}

AABB aabbInter(AABB b1, AABB b2) {
    return aabbNewPoints(vmax(b1.pMin, b2.pMin), vmin(b1.pMax, b2.pMax));
}

vect aabbDiag(AABB b) {
    return vsub(b.pMax, b.pMin);
}

double aabbSA(AABB b) {
    if (aabbIsEmpty(b)) return 0;

    vect d = aabbDiag(b);
    return 2 * (d.x * d.y + d.y * d.z + d.x * d.z);
}

vect aabbCentroid(AABB b) {
    return vadd(vmult(b.pMin, 0.5), vmult(b.pMax, 0.5));
}

bool aabbIsPointInside(vect p, AABB a) {
    return (p.x >= a.pMin.x && p.x <= a.pMax.x)
        && (p.y >= a.pMin.y && p.y <= a.pMax.y)
        && (p.z >= a.pMin.z && p.z <= a.pMax.z);
}

bool aabbIsWithin(AABB a, AABB b) {
    return aabbIsPointInside(a.pMin, b) && aabbIsPointInside(a.pMax, b);
}

bool aabbRayIntersects(AABB b, Ray r, vect invDir) {
    if (aabbIsEmpty(b)) return false;

    double t0 = 0;
    double t1 = r.tMax;

    for (int i = 0; i < 3; i++) {
        double tNear = (vcomp(b.pMin, i) - vcomp(r.orig, i)) * vcomp(invDir, i);
        double tFar = (vcomp(b.pMax, i) - vcomp(r.orig, i)) * vcomp(invDir, i);

        if (tNear > tFar) {
            double tmp = tNear;
            tNear = tFar;
            tFar = tmp;
        }

        t0 = max(tNear, t0);
        t1 = min(tFar, t1);

        if (t0 > t1) return false;
    }

    return true;
}