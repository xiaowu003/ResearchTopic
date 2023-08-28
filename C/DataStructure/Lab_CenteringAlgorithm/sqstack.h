#ifndef SQSTACK_H__
#define SQSTACK_H__

#define SIZE 32

typedef int type;

typedef struct node_st{
    type data[SIZE];
    int top;
}sqstack;

//创建 栈
sqstack *st_create(void);

//判断是否为空
int st_isempty(sqstack *st);

//入栈
int st_push(sqstack *st, type *data);

//出栈
int st_pop(sqstack *st, type *data);

//取栈顶元素
int st_top(sqstack *st, type *data);

//遍历
void st_travel(sqstack *st);

//销毁栈中元素
void st_destroy(sqstack *st);

#endif