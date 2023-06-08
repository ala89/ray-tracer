#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <stdatomic.h>
#include <pthread.h>
#include <unistd.h>

#include "bvh.h"
#include "vector.h"
#include "aabb.h"
#include "render.h"
#include "deque.h"
#include "misc.h"

#define EPSILON 0.00000001

Split computeObjectSplit(int* nodeRefs, int nFrags, AABB bounds, AABB centroidBounds, Fragment* fragments) {
    Split split = {.type = 0, .cost = INFINITY, .leftBounds = aabbNewEmpty(), .rightBounds = aabbNewEmpty()};

    for (int dim = 0; dim < 3; dim++) {
        double boundsMin = vcomp(centroidBounds.pMin, dim);
        double boundsMax = vcomp(centroidBounds.pMax, dim) + EPSILON;
        double binStep = (boundsMax - boundsMin) / OBJECT_BIN_COUNT;

        if (binStep <= EPSILON) continue;

        double invBoundsExt = 1 / (boundsMax - boundsMin);

        ObjectBin bins[OBJECT_BIN_COUNT];

        // Initialize bins
        for (int i = 0; i < OBJECT_BIN_COUNT; i++) {
            bins[i].bounds = aabbNewEmpty();
            bins[i].count = 0;
        }

        // Assign primitives to a bin by centroid
        for (int i = 0; i < nFrags; i++) {
            Fragment fragment = fragments[nodeRefs[i]];
            int binIndex = invBoundsExt * OBJECT_BIN_COUNT * (vcomp(fragment.centroid, dim) - boundsMin);
            bins[binIndex].bounds = aabbUnion(bins[binIndex].bounds, fragment.bounds);
            bins[binIndex].count++;
        }

        // Compute partial costs
        double costs[OBJECT_BIN_COUNT - 1];

        // First pass computes half of the SAH
        AABB leftBounds = aabbNewEmpty();
        int leftCount = 0;
        for (int i = 0; i < OBJECT_BIN_COUNT - 1; i++) {
            leftBounds = aabbUnion(leftBounds, bins[i].bounds);
            leftCount += bins[i].count;
            costs[i] = aabbSA(leftBounds) * leftCount;
        }

        // Second pass computes the second half of the SAH and finds the best split
        double minPartialCost = INFINITY;
        int minSplitIndex = -1;

        AABB rightBounds = aabbNewEmpty();
        int rightCount = 0;
        for (int j = OBJECT_BIN_COUNT - 1; j > 0; j--) {
            rightBounds = aabbUnion(rightBounds, bins[j].bounds);
            rightCount += bins[j].count;
            
            double partialCost = costs[j - 1] + aabbSA(rightBounds) * rightCount;
            if (partialCost < minPartialCost) {
                minSplitIndex = j - 1;
                minPartialCost = partialCost;
            }
        }

        assert(minSplitIndex >= 0);

        // Calculate actual cost and update the split object
        double cost = TRAVERSAL_COST + minPartialCost / aabbSA(bounds);
        if (cost < split.cost) {
            split.splitIndex = minSplitIndex;
            split.cost = cost;
            split.dim = dim;

            // Recompute children bounds and primitive counts
            split.leftCount = 0;
            split.rightCount = 0;
            split.leftBounds = aabbNewEmpty();
            split.rightBounds = aabbNewEmpty();

            for (int i = 0; i <= minSplitIndex; i++) {
                split.leftCount += bins[i].count;
                split.leftBounds = aabbUnion(split.leftBounds, bins[i].bounds);
            }
            for (int j = minSplitIndex + 1; j < OBJECT_BIN_COUNT; j++) {
                split.rightCount += bins[j].count;
                split.rightBounds = aabbUnion(split.rightBounds, bins[j].bounds);
            }
        }
    }

    return split;
}

void partitionObject(int* nodeRefs, int nFrags, Split split, int* leftChildRefs, int* rightChildRefs, AABB centroidBounds, Fragment* fragments) {
    double boundsMin = vcomp(centroidBounds.pMin, split.dim);
    double boundsMax = vcomp(centroidBounds.pMax, split.dim) + EPSILON;
    double invBoundsExt = 1 / (boundsMax - boundsMin);

    int leftI = 0;
    int rightI = 0;
    for (int i = 0; i < nFrags; i++) {
        Fragment fragment = fragments[nodeRefs[i]];
        int binIndex = invBoundsExt * OBJECT_BIN_COUNT * (vcomp(fragment.centroid, split.dim) - boundsMin);

        if (binIndex <= split.splitIndex) {
            leftChildRefs[leftI++] = nodeRefs[i];
        }
        else {
            rightChildRefs[rightI++] = nodeRefs[i];
        }
    }
}

Split computeSpatialSplit(int* nodeRefs, int nFrags, AABB bounds, triangle* triangles, Fragment* fragments) {
    Split split = {.type = 1, .cost = INFINITY};

    for (int dim = 0; dim < 3; dim++) {
        double boundsMin = vcomp(bounds.pMin, dim);
        double boundsMax = vcomp(bounds.pMax, dim) + EPSILON;
        double binStep = (boundsMax - boundsMin) / SPATIAL_BIN_COUNT;

        if (binStep <= EPSILON) continue;

        double invBoundsExt = 1 / (boundsMax - boundsMin);

        SpatialBin bins[SPATIAL_BIN_COUNT];

        // Initialize bins
        for (int i = 0; i < OBJECT_BIN_COUNT; i++) {
            bins[i].bounds = aabbNewEmpty();
            bins[i].entries = 0;
            bins[i].exits = 0;
        }

        // Assign primitives to the bins they overlap
        for (int i = 0; i < nFrags; i++) {
            Fragment fragment = fragments[nodeRefs[i]];
            int binMin = invBoundsExt * SPATIAL_BIN_COUNT * (vcomp(fragment.bounds.pMin, dim) - boundsMin);
            int binMax = invBoundsExt * SPATIAL_BIN_COUNT * (vcomp(fragment.bounds.pMax, dim) - boundsMin);

            bins[binMin].entries++;
            bins[binMax].exits++;

            for (int j = binMin; j <= binMax; j++) {
                double leftPlanePos = boundsMin + j * binStep;
                double rightPlanePos = boundsMin + (j + 1) * binStep;
                AABB clippedTriangleBounds = triangleChoppedBounds(triangles[fragment.primIndex], dim, leftPlanePos, rightPlanePos);
                
                // Restrict the clipped triangle bounds to the actual fragment
                AABB clippedBounds = aabbInter(clippedTriangleBounds, fragment.bounds);
                bins[j].bounds = aabbUnion(bins[j].bounds, clippedBounds);
            }
        }

        // Compute partial costs
        double costs[SPATIAL_BIN_COUNT - 1];

        // First pass computes half of the SAH
        AABB leftBounds = aabbNewEmpty();
        int leftCount = 0;
        for (int i = 0; i < OBJECT_BIN_COUNT - 1; i++) {
            leftBounds = aabbUnion(leftBounds, bins[i].bounds);
            leftCount += bins[i].entries;
            costs[i] = aabbSA(leftBounds) * leftCount;
        }

        // Second pass computes the second half of the SAH and finds the best split
        double minPartialCost = INFINITY;
        int minSplitIndex = -1;

        AABB rightBounds = aabbNewEmpty();
        int rightCount = 0;
        for (int j = OBJECT_BIN_COUNT - 1; j > 0; j--) {
            rightBounds = aabbUnion(rightBounds, bins[j].bounds);
            rightCount += bins[j].exits;
            
            double partialCost = costs[j - 1] + aabbSA(rightBounds) * rightCount;
            if (partialCost < minPartialCost) {
                minSplitIndex = j - 1;
                minPartialCost = partialCost;
            }
        }

        assert(minSplitIndex >= 0);

        // Calculate actual cost and update the split object
        double cost = TRAVERSAL_COST + minPartialCost / aabbSA(bounds);
        if (cost < split.cost) {
            split.splitIndex = minSplitIndex;
            split.cost = cost;
            split.dim = dim;

            // Recompute primitive counts
            split.leftCount = 0;
            split.rightCount = 0;

            for (int i = 0; i <= minSplitIndex; i++) {
                split.leftCount += bins[i].entries;
            }
            for (int j = minSplitIndex + 1; j < SPATIAL_BIN_COUNT; j++) {
                split.rightCount += bins[j].exits;
            }
        }
    }

    return split;
}

void partitionSpatial(int* nodeRefs, int nFrags, Split split, int* leftChildRefs, int* rightChildRefs, AABB bounds, atomic_int* nextFragment, triangle* triangles, Fragment* fragments) {
    double boundsMin = vcomp(bounds.pMin, split.dim);
    double boundsMax = vcomp(bounds.pMax, split.dim) + EPSILON;
    double binStep = (boundsMax - boundsMin) / SPATIAL_BIN_COUNT;
    double invBoundsExt = 1 / (boundsMax - boundsMin);
    double planePos = boundsMin + binStep * (split.splitIndex + 1);

    int leftI = 0;
    int rightI = 0;
    for (int i = 0; i < nFrags; i++) {
        Fragment fragment = fragments[nodeRefs[i]];
        int binMin = invBoundsExt * SPATIAL_BIN_COUNT * (vcomp(fragment.bounds.pMin, split.dim) - boundsMin);
        int binMax = invBoundsExt * SPATIAL_BIN_COUNT * (vcomp(fragment.bounds.pMax, split.dim) - boundsMin);

        if (binMax <= split.splitIndex) {
            leftChildRefs[leftI++] = nodeRefs[i];
        }
        else if (binMin > split.splitIndex) {
            rightChildRefs[rightI++] = nodeRefs[i];
        }
        else { // Duplicated ref
            AABB leftTriangleBounds;
            AABB rightTriangleBounds;

            triangleSplitBounds(triangles[fragment.primIndex], split.dim, planePos, &leftTriangleBounds, &rightTriangleBounds);

            // Restrict the split triangle bounds to the actual fragment bounds
            AABB leftBounds = aabbInter(leftTriangleBounds, fragment.bounds);
            AABB rightBounds = aabbInter(rightTriangleBounds, fragment.bounds);

            // Update left prim info
            fragments[nodeRefs[i]].bounds = leftBounds;
            fragments[nodeRefs[i]].centroid = aabbCentroid(leftBounds);

            // Update right prim info
            int newFragmentI = atomic_fetch_add(nextFragment, 1);
            fragments[newFragmentI].primIndex = fragment.primIndex;
            fragments[newFragmentI].bounds = rightBounds;
            fragments[newFragmentI].centroid = aabbCentroid(rightBounds);

            leftChildRefs[leftI++] = nodeRefs[i];
            rightChildRefs[rightI++] = newFragmentI;
        }
    }
}

void outputNodeBounds(AABB bounds, FILE* f) {
    vect pMin = bounds.pMin;
    vect pMax = bounds.pMax;
    fprintf(f, "v %f %f %f\n", pMin.x, pMin.y, pMin.z);
    fprintf(f, "v %f %f %f\n", pMax.x, pMin.y, pMin.z);
    fprintf(f, "v %f %f %f\n", pMin.x, pMin.y, pMax.z);
    fprintf(f, "v %f %f %f\n", pMax.x, pMin.y, pMax.z);
    fprintf(f, "v %f %f %f\n", pMin.x, pMax.y, pMin.z);
    fprintf(f, "v %f %f %f\n", pMax.x, pMax.y, pMin.z);
    fprintf(f, "v %f %f %f\n", pMin.x, pMax.y, pMax.z);
    fprintf(f, "v %f %f %f\n", pMax.x, pMax.y, pMax.z);
    fprintf(f, "f %d %d %d %d\n", -8, -7, -5, -6);
    fprintf(f, "f %d %d %d %d\n", -8, -6, -2, -4);
    fprintf(f, "f %d %d %d %d\n", -8, -7, -3, -4);
    fprintf(f, "f %d %d %d %d\n", -4, -3, -1, -2);
    fprintf(f, "f %d %d %d %d\n", -7, -5, -1, -3);
    fprintf(f, "f %d %d %d %d\n", -6, -5, -1, -2);
}

void outputBVH(BVH* bvh) {
    FILE* f = fopen("./bvh.obj", "w");
    
    for (int i = 0; i < bvh->nNodes; i++) {
        BVHNode node = bvh->nodes[i];
        outputNodeBounds(node.bounds, f);
    }

    fclose(f);
}

void makeLeafNode(int* nodeRefs, int nFrags, AABB bounds, int nodeIndex, BVHNode* nodes, atomic_int* nextRef, Fragment* fragments, int* references) {
    int offset = atomic_fetch_add(nextRef, nFrags);
    nodes[nodeIndex].bounds = bounds;
    nodes[nodeIndex].nPrims = nFrags;
    nodes[nodeIndex].refsOffset = offset;

    for (int i = 0; i < nFrags; i++) {
        references[offset + i] = fragments[nodeRefs[i]].primIndex;
    }
}

void buildRec(buildWorkerArgs* args, taskArgs task) {
    int* nodeRefs = task.nodeRefs;
    int nFrags = task.nFrags;
    int nodeIndex = task.nodeIndex;
    BVHNode* nodes = task.nodes;
    Fragment* fragments = task.fragments;
    atomic_int* nextNode = task.nextNode;
    atomic_int* nextFragment = task.nextFragment;
    atomic_int* nextRef = task.nextRef;
    int* references = task.references;

    assert(nFrags > 0);

    // Compute node bounds and centroid bounds
    AABB bounds = aabbNewEmpty();
    AABB centroidBounds = aabbNewEmpty();

    for (int i = 0; i < nFrags; i++) {
        Fragment fragment = fragments[nodeRefs[i]];
        bounds = aabbUnion(bounds, fragment.bounds);
        centroidBounds = aabbExpand(centroidBounds, fragment.centroid);
    }

    if (nFrags <= 1) {
        makeLeafNode(nodeRefs, nFrags, bounds, nodeIndex, nodes, nextRef, fragments, references);
        (*args->nTasks)--;
        free(nodeRefs);
        return;
    }

    // Compute splits
    Split bestSplit = computeObjectSplit(nodeRefs, nFrags, bounds, centroidBounds, fragments);

    if (args->options->useSpatialSplits) {
        double lambda = aabbSA(aabbInter(bestSplit.leftBounds, bestSplit.rightBounds)) * args->invRootSA;
        assert(lambda >= 0 && lambda <= 1);

        if (lambda >= args->options->alpha) {
            Split spatialSplit = computeSpatialSplit(nodeRefs, nFrags, bounds, args->triangles, fragments);
            if (spatialSplit.cost < bestSplit.cost) bestSplit = spatialSplit;
        }
    }

    // Partition primitives according to the best split
    double leafCost = nFrags;
    if (leafCost <= bestSplit.cost) {
        makeLeafNode(nodeRefs, nFrags, bounds, nodeIndex, nodes, nextRef, fragments, references);
        (*args->nTasks)--;
        free(nodeRefs);
        return;
    }
    else {
        assert(bestSplit.leftCount + bestSplit.rightCount >= nFrags);

        int leftChildIndex = atomic_fetch_add(nextNode, 2);
        nodes[nodeIndex].bounds = bounds;
        nodes[nodeIndex].nPrims = 0;
        nodes[nodeIndex].axis = bestSplit.dim;
        nodes[nodeIndex].leftChild = leftChildIndex;
        nodes[nodeIndex].rightChild = leftChildIndex + 1; // Useless, but easier for traversal

        int* leftChildRefs = malloc(bestSplit.leftCount * sizeof(int));
        int* rightChildRefs = malloc(bestSplit.rightCount * sizeof(int));

        if (bestSplit.type == 0) partitionObject(nodeRefs, nFrags, bestSplit, leftChildRefs, rightChildRefs, centroidBounds, fragments);
        else partitionSpatial(nodeRefs, nFrags, bestSplit, leftChildRefs, rightChildRefs, bounds, nextFragment, args->triangles, fragments);

        taskArgs leftTask = {.nodeRefs = leftChildRefs, .nFrags = bestSplit.leftCount, .nodeIndex = leftChildIndex, .nextNode = nextNode, .nextFragment = nextFragment, .nextRef = nextRef, .nodes = nodes, .fragments = fragments, .references = references};
        taskArgs rightTask = {.nodeRefs = rightChildRefs, .nFrags = bestSplit.rightCount, .nodeIndex = leftChildIndex + 1, .nextNode = nextNode, .nextFragment = nextFragment, .nextRef = nextRef, .nodes = nodes, .fragments = fragments, .references = references};

        free(nodeRefs);

        (*args->nTasks)++;
        if (leftTask.nFrags > rightTask.nFrags) {
            pushRight(args->queues[args->index], leftTask);
            buildRec(args, rightTask);
        }
        else {
            pushRight(args->queues[args->index], rightTask);
            buildRec(args, leftTask);
        }
    }
}

static bool trySteal(buildWorkerArgs* data, taskArgs* t){
    int index = data->index;
    for (int i = (index + 1) % N_THREADS; i != index; i = (i + 1) % N_THREADS) {
        if (tryPopLeft(data->queues[i], t)) return true;
    }
    return false;
}

static void* buildWorker(void* args) {
    buildWorkerArgs data = *(buildWorkerArgs*) args;
    int index = data.index;
    deque* q = data.queues[index];

    while (*data.nTasks > 0) {
        taskArgs t;
        if (tryPopRight(q, &t)) {
            buildRec(&data, t);
        } else if (trySteal(&data, &t)){
            buildRec(&data, t);
        }
    }
    
    return NULL;
}

BVH* buildBVH(triangle* triangles, int n, BVHOptions* options) {
    // Initialize construction variables
    BVH* bvh = (BVH*) malloc(sizeof(BVH));
    Fragment* fragments = (Fragment*) malloc(n * MAX_REFS_FACT * sizeof(Fragment));
    int* references = (int*) malloc(n * MAX_REFS_FACT * sizeof(int));
    BVHNode* nodes = (BVHNode*) malloc(n * MAX_REFS_FACT * 2 * sizeof(BVHNode));
    atomic_int* nextNode = malloc(sizeof(atomic_int));
    atomic_int* nextFragment = malloc(sizeof(atomic_int));
    atomic_int* nextRef = malloc(sizeof(atomic_int));
    *nextNode = 1;
    *nextFragment = n;
    *nextRef = 0;

    // Compute primitives' centroids and root bounds
    int* rootRefs = malloc(n * sizeof(int));
    AABB rootBounds = aabbNewEmpty();

    for (int i = 0; i < n; i++) {
        triangle t = triangles[i];
        AABB bounds = triangleBounds(t);
        fragments[i].primIndex = i;
        fragments[i].bounds = bounds;
        fragments[i].centroid = aabbCentroid(bounds);

        rootRefs[i] = i;

        rootBounds = aabbUnion(rootBounds, bounds);
    }

    // Initialize work stealing
    pthread_t* threads = (pthread_t*) malloc(N_THREADS * sizeof(pthread_t));
    buildWorkerArgs* args = (buildWorkerArgs*) malloc(N_THREADS * sizeof(buildWorkerArgs));
    deque** queues = (deque**) malloc(N_THREADS * sizeof(deque*));
    
    for (int i = 0; i < N_THREADS; i++) {
        queues[i] = newDeque();
    }

    taskArgs rootTask = {.nodeRefs = rootRefs, .nFrags = n, .nodeIndex = 0, .nextNode = nextNode, .nextFragment = nextFragment, .nextRef = nextRef, .nodes = nodes, .fragments = fragments, .references = references};
    pushLeft(queues[0], rootTask);
    atomic_int* nTasks = malloc(sizeof(atomic_int));
    *nTasks = 1;

    for (int i = 0; i < N_THREADS; i++) {
        args[i].triangles = triangles;
        args[i].options = options;
        args[i].invRootSA = 1 / aabbSA(rootBounds);
        args[i].queues = queues;
        args[i].nTasks = nTasks;
        args[i].index = i;
        pthread_create(&threads[i], NULL, buildWorker, &args[i]);
    }
    for (int i = 0; i < N_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    // Cleanup and result
    bvh->nodes = nodes;
    bvh->references = references;
    bvh->nNodes = *nextNode;
    bvh->nRefs = *nextRef;

    free(fragments);
    free(nextNode);
    free(nextFragment);
    free(nextRef);
    free(threads);
    free(args);
    for (int i = 0; i < N_THREADS; i++) freeDeque(queues[i]);
    free(queues);
    free(nTasks);

    return bvh;
}

void intersectBVH(triangle* triangles, BVH* bvh, Ray* r, int* hits) {
    vect invDir = vinv(r->dir);
    bool dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };

    int stackSize = 1;
    int stack[64];
    stack[0] = 0;

    *hits = 0;

    while (stackSize > 0) {
        BVHNode node = bvh->nodes[stack[--stackSize]];
        (*hits)++;
        if (aabbRayIntersects(node.bounds, *r, invDir)) {
            if (node.nPrims > 0) {
                for (int i = node.refsOffset; i < node.refsOffset + node.nPrims; i++) {
                    triangle tri = triangles[bvh->references[i]];
                    triangleRayIntersection(tri, r);
                    (*hits)++;
                }
            }
            else {
                if (dirIsNeg[node.axis]) {
                    stack[stackSize++] = node.leftChild;
                    stack[stackSize++] = node.rightChild;
                }
                else {
                    stack[stackSize++] = node.rightChild;
                    stack[stackSize++] = node.leftChild;
                }
            }
        }
    }
}