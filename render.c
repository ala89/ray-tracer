#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>

#include "render.h"
#include "bvh.h"
#include "vector.h"
#include "ray.h"
#include "color.h"
#include "clock.h"
#include "misc.h"
#include "triangle.h"

const color bgColor = {.r = 0.1, .g = 0.1, .b = 0.1};

color hueToRgb(double hue) { // 0 = Turquoise, >= 180 = Red
    if (hue >= 180) return cnew(1, 0, 0);
    else if (hue >= 120) return cnew(1, (180 - hue) / 60, 0);
    else if (hue >= 60) return cnew((hue - 60) / 60, 1, 0);
    else return cnew(0, 1, (60 - hue) / 60);
}

color rayColor(rtx* rt, renderOptions* opt, Ray* r, int* hits) {
    // Intersection
    if (opt->useBVH) {
        intersectBVH(rt->triangles, rt->bvh, r, hits);
    }
    else {
        for (int i = 0; i < rt->nT; i++) {
            triangleRayIntersection(rt->triangles[i], r);
        }
    }

    if (!r->hasHit) return bgColor;

    // Shading
    vect n0 = rt->normals[r->hitInfo.n[0]];
    vect n1 = rt->normals[r->hitInfo.n[1]];
    vect n2 = rt->normals[r->hitInfo.n[2]];
    double u = r->hitInfo.u;
    double v = r->hitInfo.v;

    vect n = vadd(vmult(n0, 1 - u - v), vadd(vmult(n1, u), vmult(n2, v)));
    double shade = fabsl(vdot(n, r->dir));
    return cnew(shade, shade, shade);
}

static void* renderWorker(void* args) {
    renderWorkerArgs data = *(renderWorkerArgs*)args;
    renderOptions* opt = data.options;
    color** image = data.image;
    color** heatmap = data.heatmap;
    
    int pixelTot = opt->imW * opt->imH;

    while (true) {
        int chunkId = atomic_fetch_add(data.nextChunk, 1);
        int start = chunkId * PIXEL_CHUNCK_SIZE;
        int end = min((chunkId + 1) * PIXEL_CHUNCK_SIZE, pixelTot);

        if (start >= end) break;

        for (int k = start; k < end; k++) {
            int i = k / opt->imW;
            int j = k % opt->imW;
            double u = (double) j / opt->imW;
            double v = (double) i / opt->imH;

            int hits;
            Ray r = pointAt(data.orig, vadd(data.corner, vadd(vmult(data.horiz, u), vmult(data.vert, v))));
            image[opt->imH - i - 1][j] = rayColor(data.rt, opt, &r, &hits);

            if (opt->heatmap) {
                heatmap[opt->imH - i - 1][j] = hueToRgb((double) hits * 180 / opt->redThreshold);
            }
        }
    }
    
    return NULL;
}

void render(rtx* rt, renderOptions* opt, color** image, color** heatmap) {
    // Calculate a few useful vectors
    vect dir = vneg(vnew(sin(opt->theta) * sin(opt->phi), cos(opt->theta), sin(opt->theta) * cos(opt->phi)));
    vect orig = vadd(opt->focus, vmult(dir, -opt->r));
    vect horiz = vmult(vnew(cos(opt->phi), 0, -sin(opt->phi)), opt->vpW);
    vect vert = vneg(vmult(vnew(cos(opt->theta) * sin(opt->phi), -sin(opt->theta), cos(opt->theta) * cos(opt->phi)), opt->vpH));
    vect corner = vsub(vadd(orig, vmult(dir, opt->focL)), vadd(vdiv(horiz, 2), vdiv(vert, 2)));

    // Init render parallelization
    pthread_t* threads = malloc(N_THREADS * sizeof(pthread_t));
    
    atomic_int* nextChunk = malloc(sizeof(atomic_int));
    *nextChunk = 0;
    renderWorkerArgs args = { .rt = rt, .options = opt, .image = image, .heatmap = heatmap, .orig = orig, .corner = corner, .horiz = horiz, .vert = vert, .nextChunk = nextChunk };

    for (int i = 0; i < N_THREADS; i++) {
        pthread_create(&threads[i], NULL, renderWorker, &args);
    }
    for (int i = 0; i < N_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
    free(nextChunk);
}
