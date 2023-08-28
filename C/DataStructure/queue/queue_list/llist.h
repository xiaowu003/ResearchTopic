#ifndef LLIST_H__
#define LLIST_H__

#define LLIST_FORWARD__     1
#define LLIST_BACKWORD__    2

//将 形为 void (const void *) 的函数命名为 llist_op
typedef void llist_op(const void *);

//将形为 int (const void *, const void *) 的函数命名为 LList_cmp
typedef int llist_cmp(const void *, const void  *);

//普通节点的结构体
struct llist_node_st {
    //void *data;
    struct llist_node_st *prev;
    struct llist_node_st *next;
    char data[0];
};

//头结点的结构体
typedef struct llist_head{
    int size;
    struct llist_node_st head;

    int (*insert)(struct llist_head *, const void *, int);
    void (*find)(struct llist_head *, const void *, llist_cmp *);
    int (*delete)(struct llist_head *, const void *, llist_cmp*);
    int (*fetch)(struct llist_head *, const void *, llist_cmp *, void *);
    void (*travel)(struct llist_head *, llist_op *);
}LLIST;

//实现方法

LLIST * llist_create(int initsize);

#if 0
//数据插入函数
// mode代表数据插入方式，头插和尾插
int llist_insert(LLIST *, const void *data, int mode);

void * llist_find(LLIST *, const void *key, llist_cmp *cmp);

int llist_delete(LLIST * ptr, const void *key, llist_cmp *cmp);

int llist_fetch(LLIST *ptr, const void *key, llist_cmp *cmp, void *data);

void llist_travel(LLIST *ptr, llist_op *op);


#endif

void llist_destroy(LLIST *ptr);

#endif