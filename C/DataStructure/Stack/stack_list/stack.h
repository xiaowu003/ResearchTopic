#ifndef STACK_H__
#define STACK_H__

#include "llist.h"

typedef LLIST STACK;

STACK *stack_creat(int );

int stack_push(STACK *ptr, const void *data);

int stack_pop(STACK *ptr, void *data);

void stack_destroy(STACK *ptr);

#endif