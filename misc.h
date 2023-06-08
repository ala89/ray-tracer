#ifndef MISC_H
#define MISC_H

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

#define MAX_VERTICES 10000000
#define MAX_TRIANGLES 20000000
#define MAX_NORMALS 3 * MAX_TRIANGLES
#define N_THREADS 8
#define PIXEL_CHUNCK_SIZE 1000
#define MAX_REFS_FACT 3
#define OBJECT_BIN_COUNT 16
#define SPATIAL_BIN_COUNT 16
#define TRAVERSAL_COST 0.125
#define DEFAULT_ALPHA 0.00001

#endif