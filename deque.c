#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include <stdbool.h>

#include "deque.h"

const taskArgs DEFAULT = {};

node* newNode(taskArgs data) {
    node* n = malloc(sizeof(node));
    n->next = NULL;
    n->prev = NULL;
    n->data = data;
    return n;
}

deque* newDeque(void) {
    deque* q = malloc(sizeof(deque));
    q->sentinel = newNode(DEFAULT);
    q->sentinel->next = q->sentinel;
    q->sentinel->prev = q->sentinel;
    pthread_mutex_init(&q->lock, NULL);
    return q;
}

bool tryPopLeft(deque* q, taskArgs* res) {
    pthread_mutex_lock(&q->lock);
    node* head = q->sentinel->next;
    if (head == q->sentinel) {
        pthread_mutex_unlock(&q->lock);
        return false;
    }
    head->prev->next = head->next;
    head->next->prev = head->prev;
    *res = head->data;
    free(head);
    pthread_mutex_unlock(&q->lock);
    return true;
}

bool tryPopRight(deque* q, taskArgs* res) {
    pthread_mutex_lock(&q->lock);
    node* tail = q->sentinel->prev;
    if (tail == q->sentinel) {
        pthread_mutex_unlock(&q->lock);
        return false;
    }
    tail->next->prev = tail->prev;
    tail->prev->next = tail->next;
    *res = tail->data;
    free(tail);
    pthread_mutex_unlock(&q->lock);
    return true;
}

void pushLeft(deque* q, taskArgs data) {
    pthread_mutex_lock(&q->lock);
    node* n = newNode(data);
    n->prev = q->sentinel;
    n->next = q->sentinel->next;
    n->prev->next = n;
    n->next->prev = n;
    pthread_mutex_unlock(&q->lock);
}

void pushRight(deque* q, taskArgs data) {
    pthread_mutex_lock(&q->lock);
    node* n = newNode(data);
    n->next = q->sentinel;
    n->prev = q->sentinel->prev;
    n->prev->next = n;
    n->next->prev = n;
    pthread_mutex_unlock(&q->lock);
}

void freeDeque(deque* q) {
    pthread_mutex_lock(&q->lock);
    node* n = q->sentinel->next;
    while (n != q->sentinel) {
        node* tmp = n->next;
        free(n);
        n = tmp;
    }
    free(q->sentinel);
    pthread_mutex_unlock(&q->lock);
    pthread_mutex_destroy(&q->lock);
    free(q);
}