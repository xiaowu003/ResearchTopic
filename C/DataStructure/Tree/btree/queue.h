#ifndef QUEUE_H__
#define QUEUE_H__

#define MAXSIZE 5

typedef int datatype;

typedef struct node_st{
    datatype data[MAXSIZE];
    int head, tail;
}queue;

//创刊队列
queue *qu_creat();

//判断是否为空
int qu_isempty();

//入队
int qu_enqueue(queue *sq, datatype *x);

//出队
int qu_dequeue(queue *sq, datatype *x);

//遍历
void qu_travel(queue *sq);

//清空队列
void qu_clear(queue *sq);

//销毁队列
void qu_destroy(queue *sq);
#endif