#include <time.h>

#include "clock.h"

double deltaT(timestamp t0, timestamp t1) {
    return (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
}

timestamp getTime(void) {
    timestamp t;
    clock_gettime(CLOCK_REALTIME, &t);
    return t;
}