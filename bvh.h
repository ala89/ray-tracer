#ifndef BVH_H
#define BVH_H

#include "deque.h"
#include "vector.h"
#include "aabb.h"
#include "render.h"

struct Fragment {
    int primIndex;
    AABB bounds;
    vect centroid;
};

typedef struct Fragment Fragment;

struct BVHNode {
    AABB bounds;
    int nPrims;
    int refsOffset;
    int axis;
    int leftChild;
    int rightChild;
};

typedef struct BVHNode BVHNode;

struct BVH {
    BVHNode* nodes;
    int* references;
    int nNodes;
    int nRefs;
};

typedef struct BVH BVH;

struct Split {
    int splitIndex;
    double cost;
    int dim;

    int type; // Object = 0, Spatial = 1

    int leftCount;
    int rightCount;
    AABB leftBounds;
    AABB rightBounds;
};

typedef struct Split Split;

struct ObjectBin {
    AABB bounds;
    int count;
};

typedef struct ObjectBin ObjectBin;

struct SpatialBin {
    AABB bounds;
    int entries;
    int exits;
};

typedef struct SpatialBin SpatialBin;

struct BVHOptions {
    bool useSpatialSplits;
    double alpha;
};

typedef struct BVHOptions BVHOptions;

struct buildWorkerArgs {
    triangle* triangles;
    BVHOptions* options;
    double invRootSA;
    deque** queues;
    atomic_int* nTasks;
    int index;
};

typedef struct buildWorkerArgs buildWorkerArgs;

BVH* buildBVH(triangle* triangles, int n, BVHOptions* options);

void intersectBVH(triangle* triangles, BVH* bvh, Ray* r, int* hits);

void outputBVH(BVH* bvh);

#endif