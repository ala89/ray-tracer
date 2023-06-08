#ifndef RAY_H
#define RAY_H

#include <stdbool.h>

#include "vector.h"
#include "color.h"

struct HitInfo {
    int n[3];
    double u, v;
};

typedef struct HitInfo HitInfo;

struct Ray {
    vect orig;
    vect dir;
    
    bool hasHit;
    double tMax;
    HitInfo hitInfo;
};

typedef struct Ray Ray;

Ray pointAt(vect orig, vect to);

#endif