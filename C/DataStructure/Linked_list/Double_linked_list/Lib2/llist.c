#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "llist.h"

int llist_insert(LLIST *, const void *data, int mode);
void * llist_find(LLIST *, const void *key, llist_cmp *cmp);
int llist_delete(LLIST * ptr, const void *key, llist_cmp *cmp);
int llist_fetch(LLIST *ptr, const void *key, llist_cmp *cmp, void *data);
void llist_travel(LLIST *ptr, llist_op *op);

//创建一个链表, initsize 代表用户传入数据的大小
LLIST * llist_create(int initsize){
    LLIST* new;

    new = malloc(sizeof(*new));
    if(new == NULL){
        return NULL;
    }

    new->size = initsize;                               //用户传入数据的大小
    new->head.prev = &new->head;                        //
    new->head.next = &new->head;                        //

    new->insert = llist_insert;
    new->find = llist_find;
    new->delete = llist_delete;
    new->fetch = llist_fetch;
    new->travel = llist_travel;

    return new;
}

//数据插入函数
//mode代表数据插入方式，头插和尾插
int llist_insert(LLIST *ptr, const void *data, int mode){
    struct llist_node_st *newnode;

    newnode = malloc(sizeof(*newnode) + ptr->size);
    if(newnode == NULL){
        return -1;
    }

    //拷贝函数，目的内存：newnode->data; 源数据：data; 拷贝数据量：ptr->size;
    memcpy(newnode->data, data, ptr->size);

    if(mode == LLIST_FORWARD__){                        //首部插入
        newnode->prev = &ptr->head;
        newnode->next = ptr->head.next;
    }
    else if(mode == LLIST_BACKWORD__){                  //尾部插入
        newnode->prev = ptr->head.prev;
        newnode->next = &ptr->head;
    }
    else{                                               //Error!
        return -3;
    }

    newnode->prev->next = newnode;
    newnode->next->prev = newnode;

    return 0;
}

static struct llist_node_st * find_(LLIST *ptr, const void *key, llist_cmp *cmp){
    struct llist_node_st *cur;

    for(cur = ptr->head.next; cur != &ptr->head; cur = cur->next){
        if(cmp(key, cur->data) == 0){
            break;
        }
    }

    return cur;
}

//查找函数
void * llist_find(LLIST *ptr, const void *key, llist_cmp *cmp){
    struct llist_node_st *node;

    node = find_(ptr, key, cmp);
    if(node == &ptr->head){
        return NULL;
    }

    return node->data;
}

//找到符合条件的节点，并将其脱链删除
int llist_delete(LLIST * ptr, const void *key, llist_cmp *cmp){
    struct llist_node_st *node;

    node = find_(ptr, key, cmp);
    if(node == &ptr->head){
        return -1;                      //返回值是头节点的head地址，说明没有找到符合条件的地址，终止删除操作
    }

    node->prev->next = node->next;      //脱链操作
    node->next->prev = node->prev;

    free(node);

    return 0;
}

int llist_fetch(LLIST *ptr, const void *key, llist_cmp *cmp, void *data){
    struct llist_node_st *node;

    node = find_(ptr, key, cmp);
    if(node == &ptr->head){
        return -1;
    }

    node->prev->next = node->next;
    node->next->prev = node->prev;

    if(data != NULL){
        memcpy(data, node->data, ptr->size);
    }

    free(node);

    return 0;
}


void llist_travel(LLIST *ptr, llist_op *op){
    struct llist_node_st *cur;

    for(cur  = ptr->head.next; cur != &ptr->head; cur = cur->next){
        op(cur->data);
    }
}

void llist_destroy(LLIST* ptr){
    struct llist_node_st *cur, *next;

    for(cur = ptr->head.next; cur != &ptr->head; cur = next){
        next = cur->next;                               //free cur 之前，需要一个指针将 cur 的 next 保存下来

        free(cur);
    }

    free(ptr);
}