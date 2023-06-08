#include <float.h>

#include "vector.h"
#include "triangle.h"
#include "ray.h"
#include "aabb.h"

void swapVects(vect* arr, int i, int j) {
    vect tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

void triangleRayIntersection(triangle tri, Ray* r) {
    vect v0 = tri.vertices[tri.v[0]],
        v1 = tri.vertices[tri.v[1]],
        v2 = tri.vertices[tri.v[2]];

    vect v0v1 = vsub(v1, v0);
    vect v0v2 = vsub(v2, v0);
    vect pVec = vcross(r->dir, v0v2);
    double det = vdot(v0v1, pVec);
    
    if (fabsl(det) < DBL_EPSILON) return;
    
    double invDet = (double) 1 / det;
    vect tVec = vsub(r->orig, v0);
    double u = vdot(tVec, pVec) * invDet;
    
    if (u < 0 || u > 1) return;
    
    vect qVec = vcross(tVec, v0v1);
    double v = vdot(r->dir, qVec) * invDet;
    
    if (v < 0 || u + v > 1) return;

    double t = vdot(v0v2, qVec) * invDet;
    if (t < r->tMax) {
        r->hitInfo.n[0] = tri.n[0];
        r->hitInfo.n[1] = tri.n[1];
        r->hitInfo.n[2] = tri.n[2];
        r->hitInfo.u = u;
        r->hitInfo.v = v;
        r->tMax = t;
        r->hasHit = true;
    }
}

int trianglePlaneIntersection(vect vertices[3], int dim, double planePos, vect intersections[2]) {
    // If all vertices lie on the same side of the plane, no intersection
    if (vcomp(vertices[2], dim) < planePos || vcomp(vertices[0], dim) >= planePos) {
        return 0;
    }

    // Lerp vertices to find the 2 intersections
    int isxCount = 0;
    for (int i = 0; i < 2; i++) {
        vect from = vertices[i];
        for (int j = i + 1; j < 3; j++) {
            vect to = vertices[j];
            if (vcomp(from, dim) < planePos && vcomp(to, dim) >= planePos) {
                double delta = vcomp(to, dim) - vcomp(from, dim);
                double t = (planePos - vcomp(from, dim)) / delta;
                intersections[isxCount++] = vadd(vmult(from, 1 - t), vmult(to, t));
            }
        }
    }

    assert(isxCount == 2);

    return 2;
}

AABB triangleChoppedBounds(triangle tri, int dim, double leftPlanePos, double rightPlanePos) {
    vect vertices[3] = { tri.vertices[tri.v[0]], tri.vertices[tri.v[1]], tri.vertices[tri.v[2]] };

    // Sort vertices along current dim
    if (vcomp(vertices[1], dim) < vcomp(vertices[0], dim)) swapVects(vertices, 0, 1);
    if (vcomp(vertices[2], dim) < vcomp(vertices[1], dim)) swapVects(vertices, 1, 2);
    if (vcomp(vertices[1], dim) < vcomp(vertices[0], dim)) swapVects(vertices, 0, 1);

    // Compute intersections
    vect leftIsx[2];
    vect rightIsx[2];
    int leftIsxCount = trianglePlaneIntersection(vertices, dim, leftPlanePos, leftIsx);
    int rightIsxCount = trianglePlaneIntersection(vertices, dim, rightPlanePos, rightIsx);

    AABB clipped = aabbUnion(aabbFromPoints(leftIsx, leftIsxCount), aabbFromPoints(rightIsx, rightIsxCount));
    for (int i = 0; i < 3; i++) {
        if (vcomp(vertices[i], dim) >= leftPlanePos && vcomp(vertices[i], dim) < rightPlanePos) {
            clipped = aabbExpand(clipped, vertices[i]);
        }
    }

    return clipped;
}

void triangleSplitBounds(triangle tri, int dim, double planePos, AABB* leftBounds, AABB* rightBounds) {
    vect vertices[3] = { tri.vertices[tri.v[0]], tri.vertices[tri.v[1]], tri.vertices[tri.v[2]] };

    // Sort vertices along current dim
    if (vcomp(vertices[1], dim) < vcomp(vertices[0], dim)) swapVects(vertices, 0, 1);
    if (vcomp(vertices[2], dim) < vcomp(vertices[1], dim)) swapVects(vertices, 1, 2);
    if (vcomp(vertices[1], dim) < vcomp(vertices[0], dim)) swapVects(vertices, 0, 1);

    vect isx[2];
    int isxCount = trianglePlaneIntersection(vertices, dim, planePos, isx);

    *leftBounds = *rightBounds = aabbFromPoints(isx, isxCount);
    for (int i = 0; i < 3; i++) {
        if (vcomp(vertices[i], dim) >= planePos) {
            *rightBounds = aabbExpand(*rightBounds, vertices[i]);
        }
        else {
            *leftBounds = aabbExpand(*leftBounds, vertices[i]);
        }
    }
}

AABB triangleBounds(triangle t) {
    return aabbExpand(aabbExpand(aabbNewPoint(t.vertices[t.v[0]]), t.vertices[t.v[1]]), t.vertices[t.v[2]]);
}