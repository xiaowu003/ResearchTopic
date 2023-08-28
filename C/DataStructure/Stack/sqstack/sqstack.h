#ifndef SQSTACK_H__
#define SQSTACK_H__

#define MAXSIZE 5

typedef int datatype;

typedef struct node_st{
    datatype data[MAXSIZE];
    int top;
}sqstack;

//创建 栈
sqstack *st_create(void);

//判断是否为空
int st_isempty(sqstack *st);

//入栈
int st_push(sqstack *st, datatype *data);

//出栈
int st_pop(sqstack *st, datatype *data);

//取栈顶元素
int st_top(sqstack *st, datatype *data);

//遍历
void st_travel(sqstack *st);

//销毁栈中元素
void st_destroy(sqstack *st);

#endif