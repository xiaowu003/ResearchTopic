//顺序循环队列

#include <stdio.h>
#include <stdlib.h>

#include "queue.h"

//创建队列
queue *qu_creat(){
    queue *sq;

    sq = malloc(sizeof(*sq));
    if(sq == NULL){
        return NULL;
    }

    sq->head = 0;
    sq->tail = 0;

    return sq;
}

//判断是否为空
int qu_isempty(queue *sq){
    //成立为空，否则为非空
    return (sq->head == sq->tail);
}

//入队，入在队尾
int qu_enqueue(queue *sq, datatype *x){
    //先判断队列是否已满，通过tail + 1 与 容量  取模 的结果和head的值判断
    if((sq->tail + 1) % MAXSIZE == sq->head){
        return -1;
    }

    //有新元素入队，更新tail的值，实际上空间只有5个，所以不能直接加1，有可能会发生 数组越界，通过取模运算使其加一
    sq->tail = (sq->tail + 1) % MAXSIZE;

    //入队
    sq->data[sq->tail] = *x;

    return 0;
}

//出队
int qu_dequeue(queue *sq, datatype *x){
    if(qu_isempty(sq)){
        return -1;
    }
    
    //head位置是不存储数据的，所以(sq->head + 1)%MAXSIZE所对应的数据才是要出队的数据
    sq->head = (sq->head + 1) % MAXSIZE;
    *x = sq->data[sq->head];

    return 0;
}

//遍历
void qu_travel(queue *sq){
    if(sq->head == sq->tail){
        return;
    }

    //打印数据，从head位置开始直到tail
    int i = (sq->head + 1) % MAXSIZE;
    while(i != sq->tail){
        printf("%d ", sq->data[i]);
        i = (i + 1) % MAXSIZE;
    }
    printf("%d\n", sq->data[i]);
}

//清空队列
void qu_clear(queue *sq){
    sq->head = sq->tail;
}

//销毁队列
void qu_destroy(queue *sq){
    free(sq);
}