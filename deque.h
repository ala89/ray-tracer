#ifndef DEQUE_H
#define DEQUE_H

#include <stdbool.h>
#include <stdatomic.h>

struct taskArgs {
    int* nodeRefs;
    int nFrags;
    int nodeIndex;
    atomic_int* nextNode;
    atomic_int* nextFragment;
    atomic_int* nextRef;
    atomic_int* nObj;
    atomic_int* nSpat;
    struct BVHNode* nodes;
    struct Fragment* fragments;
    int* references;
};

typedef struct taskArgs taskArgs;

struct node {
    taskArgs data;
    struct node* next;
    struct node* prev;
};

typedef struct node node;

struct deque {
    node* sentinel;
    pthread_mutex_t lock;
};

typedef struct deque deque;

deque* newDeque(void);

bool tryPopLeft(deque* q, taskArgs* res);

bool tryPopRight(deque* q, taskArgs* res);

void pushLeft(deque* q, taskArgs data);

void pushRight(deque* q, taskArgs data);

void freeDeque(deque* q);

#endif