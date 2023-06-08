#ifndef RTX_H
#define RTX_H

#include "vector.h"
#include "color.h"
#include "ray.h"
#include "misc.h"
#include "triangle.h"

struct rtx {
    bool useSpatialSplits;
    double alpha;
    
    vect vertices[MAX_VERTICES];
    vect normals[MAX_NORMALS];
    triangle triangles[MAX_TRIANGLES];

    struct BVH* bvh;
    
    int nV, nT;
};

typedef struct rtx rtx;

struct renderOptions {
    vect focus;
    double r;
    double theta;
    double phi;

    int imW;
    int imH;
    double vpW;
    double vpH;
    double focL;

    bool useBVH;
    bool heatmap;
    int redThreshold;
};

typedef struct renderOptions renderOptions;

struct renderWorkerArgs {
    rtx* rt;
    renderOptions* options;
    color** image;
    color** heatmap;
    atomic_int* nextChunk;
    vect orig, corner, horiz, vert;
};

typedef struct renderWorkerArgs renderWorkerArgs;

void render(rtx* rt, renderOptions* opt, color** image, color** heatmap);

#endif