#include <stdio.h>
#include <stdlib.h>
#include "nohead.h"

//无头结点链表 的头部插入函数
int list_insert(struct node_st **list, struct score_st *data){
    
    struct node_st *new;

    new = malloc(sizeof(*new));
    if(new == NULL){
        return -1;               //申请内存失败
    }

    new->data = *data;           //将数据存入新申请的空间
    new->next = *list;           //将新节点的指针指向原来的首节点
    *list = new;                 //更新首节点位置，变为新申请的节点

    return 0;    
    //TODO
    //为什么不返回指针时会丢失链表?
    /*
    mian.c中传入的是指针变量list的值，即：将链表的首地址传入了函数
    实际上函数是值传递的方式
    */
}

void list_show(struct node_st *list){
    struct node_st *cur;

    for(cur = list; cur != NULL; cur = cur->next){
        printf("%d %s %d %d\n", cur->data.id, cur->data.name, cur->data.math, cur->data.Chinese);
    }
}

//首部删除
int list_delete(struct node_st **list){
    if(*list == NULL){
        return -1;
    }

    struct node_st *cur = NULL;                 //创建一个新的指针来保存要被删掉的那一部分

    cur = *list;

    *list = (*list)->next;                      //让原来的链表指向要删掉部分的下一个节点

    free(cur);                                  //释放要删掉的那部分内存

    return 0;
}

struct score_st * list_find(struct node_st *list, int id){
    struct node_st *cur;

    for(cur = list; cur != NULL; cur = cur->next){
        if(cur->data.id == id){
            //printf("%d %s %d %d\n", cur->data.id, cur->data.name, cur->data.math, cur->data.Chinese);

            return &cur->data;
        }
    }

    return NULL;
}

void list_destroy(struct node_st *list){
    if(list == NULL){
        return ;
    }

    struct node_st *cur;

    for(cur = list; cur != NULL; cur = list){
        list = cur->next;
        free(cur);
    }
}

#if 0
int list_insert(struct node_st *list, struct score_st *data){
    
    struct node_st *new;

    new = malloc(sizeof(*new));
    if(new == NULL){
        return -1;              //申请内存失败
    }

    new->data = *data;          //将数据存入新申请的空间
    new->next = list;           //将新节点的指针指向原来的首节点
    list = new;                 //更新首节点位置，变为新申请的节点

    return 0;                
}
#endif