#include <stdio.h>
#include <stdlib.h>

#include "sqstack.h"

//创建 栈
sqstack *st_create(void){
    sqstack *st;

    st = malloc(sizeof(*st));
    if(st == NULL){
        return NULL;
    }

    st->top = -1;

    return st;
}

//判断是否为空
int st_isempty(sqstack *st){
    return (st->top == -1);             //0表示空，即false，非0值表示非空，即true
}

//入栈
int st_push(sqstack *st, type *data){
    //如果当前栈已满，不能进行入栈操作
    if(st->top == (SIZE -1)){
        return -1;
    }

    st->data[++st->top] = *data;
    return 0;
}

//出栈
int st_pop(sqstack *st, type *data){
    //先判断栈是否为空
    if(st_isempty(st)){
        return -1;
    }

    *data = st->data[st->top--];    //出栈，一个数据弹出后，数据总数减一，所以 st-top-- 为先取到st-top, 再将其减一
    return 0;
}

//取栈顶元素,但是不让他出栈，即只查看栈顶元素
int st_top(sqstack *st, type *data){
    if(st_isempty(st)){
        return -1;
    }

    *data = st->data[st->top];
    return 0;
}

//遍历
void st_travel(sqstack *st){
    //先判断栈是否为空
    if(st_isempty(st)){
        return ;
    }

    for(int i = 0; i <= st->top; i++){
        printf("%d ", st->data[i]);
    }

    printf("\n");
}


//销毁栈中元素
void st_destroy(sqstack *st){
    free(st);
}