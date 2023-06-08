#include <stdbool.h>
#include <float.h>

#include "ray.h"
#include "vector.h"
#include "color.h"

Ray pointAt(vect orig, vect to) {
    Ray r = {
        .orig = orig,
        .dir = vunit(vsub(to, orig)),
        .hasHit = false,
        .tMax = INFINITY
    };
    return r;
}