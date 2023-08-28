#ifndef SQLIST_H__
#define SQLIST_H__

#define DATASIZE 1024

typedef int datatype;

typedef struct node_st {
    datatype data[DATASIZE];
    int last;
}sqlist;

sqlist *sqlist_create();

void sqlist_create1(sqlist **ptr);

int sqlist_insert(sqlist *me, int i, datatype *data);

int sqlist_delete(sqlist *me, int i);

int sqlist_find(sqlist *me, datatype *data);

int sqlist_isempty(sqlist *me);

int sqlist_setempty(sqlist *me);

int sqlist_getnum(sqlist *me);

void sqlist_display(sqlist *me);

int sqlist_destroy(sqlist *me);

int sqlist_union(sqlist *list1, sqlist *list2);

#endif