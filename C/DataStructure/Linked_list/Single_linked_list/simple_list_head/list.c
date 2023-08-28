#include <stdio.h>
#include <stdlib.h>

#include "list.h"

list *list_creat(){

    list *me;

    me = malloc(sizeof(*me));
    if(me == NULL){
        return NULL;
    }

    me->next = NULL;

    return me;
}

//在指定位置插入数据
int list_insert_at(list *me, int i, datatype *data){

    int j = 0;
    list *node = me, *newnode;

    if(i < 0){
        return -1;           //输入的i不合法，不能为负值
    }

    while (j < i && node != NULL)
    {
        //找到要插入的位置
        node = node->next;
        j++;
    }

    if(node){

        newnode = malloc(sizeof(*newnode));
        if(newnode == NULL){
            return -2;       //申请新空间失败
        }
        
        //单向链表在某个位置插入新数据，只能通过获得他的前一个的next指向
        newnode->data = *data;
        newnode->next = node->next;
        node->next = newnode;

        return 0;
    }
    else{
        return -3;          //当前输入的位置i不合法，超出了链表长度
    }
    
}

//按从小到大的顺序插入数据
int list_order_insert(list *me, datatype *data){

    list *p = me, *q = NULL;

    while (p->next && p->next->data < *data)
    {
        //p的下一节点不为空，且下一节点的数据比要插入的数据小
        p = p->next;
    }

    q = malloc(sizeof(*q));

    if(q == NULL){
        return -1;
    }

    q->data = *data;
    q->next = p->next;
    p->next = q;

    return 0;
}

int list_delete_at(list *me, int i, datatype *data){
    int j = 0;
    list *p = me, *q = NULL;

    *data = NULL;

    if(i < 0){
        return -1;
    }

    while (j < i && p)
    {
        p = p->next;
        j++;
    }

    if(p){
        q = p->next;
        p->next = q->next;
        *data = q->data;
        free(q);
        q = NULL;
        return 0;
    }
    else{
        return -2;
    }
}

//删除链表中的某个数据
int list_delete(list *me, datatype *data){

    list *p = me, *q = NULL;

    while (p->next && p->next->data != *data)
    {
        p = p->next;
    }

    if(p->next == NULL){
        return -1;
    }
    else{
        q = p->next;
        p->next = q->next;

        free(q);
        q = NULL;
    }
    
}

int list_isempty(list *me){

    if(me->next == NULL){
        return 0;
    }

    return 1;
}

void list_display(list *me){

    list *node = me;

    if(list_isempty(me) == 0){
        return ;                                //如果链表为空，直接返回
    }

    while (node != NULL)
    {
        printf("%d ",node->data);              //依次打印数据
        node = node->next;
    }

    printf("\n");

    return ;
    
}

void list_destroy(list *me){

    list *node = NULL, *next = NULL;

    for(node = me->next; node != NULL; node = next){

        next = node->next;
        free(node);
    }

    free(me);
    
    return ;
}