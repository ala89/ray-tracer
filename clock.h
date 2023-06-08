#ifndef CLOCK_H
#define CLOCK_H

#include <time.h>

typedef struct timespec timestamp;

timestamp getTime(void);

double deltaT(timestamp t0, timestamp t1);

#endif