/**
 * HEADERS
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>

#include "render.h"
#include "triangle.h"
#include "vector.h"
#include "clock.h"
#include "color.h"
#include "bvh.h"

/**
 * FUNCTIONS
*/
void parse_obj(rtx* rt, char* filename) {
    FILE* f = fopen(filename, "r");
    assert(f != NULL);

    int iV = 0;
    int iN = 0;
    int iT = 0;

    char* line = NULL;
    long unsigned int len;
    while (getline(&line, &len, f) != EOF) {
        double x, y, z;
        int v0, v1, v2;
        int n0, n1, n2;
        if (sscanf(line, "v %lf %lf %lf", &x, &y, &z) == 3) {
            rt->vertices[iV] = vnew(x, y, z);
            iV++;
        }
        else if (sscanf(line, "vn %lf %lf %lf", &x, &y, &z) == 3) {
            vect v = vnew(x, y, z);
            rt->normals[iN] = vmag2(v) == 0 ? v : vunit(v);
            iN++;
        }
        else if (sscanf(line, "f %d//%d %d//%d %d//%d", &v0, &n0, &v1, &n1, &v2, &n2) == 6 || sscanf(line, "f %d/%*d/%d %d/%*d/%d %d/%*d/%d\n", &v0, &n0, &v1, &n1, &v2, &n2) == 6) {
            rt->triangles[iT].v[0] = v0 < 0 ? iV - v0 : v0 - 1;
            rt->triangles[iT].v[1] = v1 < 0 ? iV - v1 : v1 - 1;
            rt->triangles[iT].v[2] = v2 < 0 ? iV - v2 : v2 - 1;
            rt->triangles[iT].n[0] = n0 < 0 ? iN - n0 : n0 - 1;
            rt->triangles[iT].n[1] = n1 < 0 ? iN - n1 : n1 - 1;
            rt->triangles[iT].n[2] = n2 < 0 ? iN - n2 : n2 - 1;
            rt->triangles[iT].vertices = rt->vertices;
            iT++;
        }
        free(line);
        line = NULL;
    }
    free(line);

    rt->nV = iV;
    rt->nT = iT;

    fclose(f);
}

void write_image(color** image, int w, int h, char* filename) {
    FILE* f = fopen(filename, "wb");
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            char color[3];
            color[0] = image[i][j].r * 255;
            color[1] = image[i][j].g * 255;
            color[2] = image[i][j].b * 255;
            fwrite(color, 1, 3, f);
        }
    }
    fclose(f);
}

void parse_config(char* filename, renderOptions* opt) {
    FILE* f = fopen(filename, "r");
    assert(f != NULL);

    double u, v, vpH, focL, fx, fy, fz, r, phi, theta;
    int imH;
    fscanf(f, "RATIO: %lf/%lf\n", &u, &v);
    fscanf(f, "IM_H: %d\n", &imH);
    fscanf(f, "VP_H: %lf\n", &vpH);
    fscanf(f, "FOC_L: %lf\n", &focL);
    fscanf(f, "FOCUS: %lf %lf %lf\n", &fx, &fy, &fz);
    fscanf(f, "R: %lf\n", &r);
    fscanf(f, "THETA: %lf\n", &theta);
    fscanf(f, "PHI: %lf", &phi);

    double ratio = u / v;
    opt->imW = (int)(imH * ratio);
    opt->imH = imH;
    opt->vpW = vpH * ratio;
    opt->vpH = vpH;
    opt->focL = focL;
    
    opt->focus = vnew(fx, fy, fz);
    opt->r = r;
    opt->theta = theta * M_PI / 180;
    opt->phi = phi * M_PI / 180;

    fclose(f);
}

renderOptions* initRenderOptions() {
    renderOptions* opt = malloc(sizeof(renderOptions));
    
    opt->focus = vnew(0, 0, 0);
    opt->r = 1;
    opt->theta = 90;
    opt->phi = 0;

    opt->imW = 640;
    opt->imH = 480;
    opt->vpW = (double) 4/3;
    opt->vpH = 1;
    opt->focL = 1;

    opt->useBVH = false;
    opt->heatmap = false;
    opt->redThreshold = 200;
    
    return opt;
}

BVHOptions* initBVHOptions() {
    BVHOptions* opt = malloc(sizeof(BVHOptions));

    opt->useSpatialSplits = false;
    opt->alpha = DEFAULT_ALPHA;

    return opt;
}

color** initImage(int w, int h) {
    color** image = malloc(h * sizeof(color*));
    
    for (int i = 0; i < h; i++) image[i] = malloc(w * sizeof(color));

    return image;
}

void freeImage(color** image, int h) {
    for (int i = 0; i < h; i++) free(image[i]);
    free(image);
}

int main(int argc, char* argv[]) {
    rtx* rt = malloc(sizeof(rtx));
    rt->bvh = NULL;

    assert(argc >= 2);

    parse_obj(rt, argv[1]);
    printf("Parsed %d triangles\n", rt->nT);

    while (true) {
        char* line = NULL;
        long unsigned int len;
        getline(&line, &len, stdin);
        
        char* arg = strtok(line, " \n");

        if (arg == NULL) {
            free(line);
            continue;
        }
        
        if (strcmp(arg, "render") == 0) {
            renderOptions* opt = initRenderOptions();

            // Parse args
            arg = strtok(NULL, " \n");
            while (arg != NULL) {
                int red;

                if (strcmp(arg, "bvh") == 0) opt->useBVH = true;
                else if (strcmp(arg, "heatmap") == 0) opt->heatmap = true;
                else if (sscanf(arg, "red=%d", &red) == 1) opt->redThreshold = red;
                else parse_config(arg, opt);

                arg = strtok(NULL, " \n");
            }

            if (opt->useBVH) assert(rt->bvh != NULL);

            // Allocate image
            color** image = initImage(opt->imW, opt->imH);
            color** heatmap = initImage(opt->imW, opt->imH);
            
            // Render
            timestamp t1 = getTime();
            render(rt, opt, image, heatmap);
            printf("Render time: %.2fs\n", deltaT(t1, getTime()));

            // Finish render task
            write_image(image, opt->imW, opt->imH, "output.ppm");
            if (opt->heatmap) write_image(heatmap, opt->imW, opt->imH, "heatmap.ppm");

            freeImage(image, opt->imH);
            freeImage(heatmap, opt->imH);
            free(opt);
        }
        else if (strcmp(arg, "bvh") == 0) {
            BVHOptions* opt = initBVHOptions();
            bool doOutputBVH = false;

            // Parse args
            arg = strtok(NULL, " \n");
            while (arg != NULL) {
                double alpha;

                if (strcmp(arg, "spatial") == 0) opt->useSpatialSplits = true;
                else if (strcmp(arg, "output") == 0) doOutputBVH = true;
                else if (sscanf(arg, "alpha=%lf", &alpha) == 1) opt->alpha = alpha;

                arg = strtok(NULL, " \n");
            }

            timestamp t0 = getTime();
            rt->bvh = buildBVH(rt->triangles, rt->nT, opt);
            printf("BVH build time: %.2fs\n", deltaT(t0, getTime()));
            printf("Nodes: %d\nReferences: %d\nRatio: %.2f\n", rt->bvh->nNodes, rt->bvh->nRefs, (double) rt->bvh->nRefs / rt->nT);

            if (doOutputBVH) outputBVH(rt->bvh);
            
            free(opt);
        }
        else {
            printf("Didn't get it\n");
        }
        
        free(line);
    }
}