#include <stdio.h>
#include <stdlib.h>

#include "llist.h"

#define NAMESIZE 32

//用户自定义数据结构
struct score_st {
    int id;
    char name[NAMESIZE];
    int math;
    int Chinese;
};

//用户自定义的函数，打印具体数据
void print_s(const void *record){
    const struct score_st *r = record;
    printf("%d %s %d %d\n", r->id, r->name, r->math, r->Chinese);
}

//用户自定义比较函数，通过id进行比较
static int id_cmp(const void *key, const void *record){
    const int *k = key;
    const struct score_st *r = record;

    return (*k - r->id);
}

//用户自定义函数，通过name进行比较
static int name_cmp(const void *key, const void *record){
    const char *k = key;
    const struct score_st *r = record;

    return strcmp(k, r->name);
}

int main(){
    LLIST* handler;
    struct score_st tmp;
    

    handler = llist_create(sizeof(struct score_st));
    if(handler == NULL){
        exit(0);
    }

    //将用户数据赋值，并插入
    for(int i = 0; i < 7; i++){
        tmp.id = i;
        snprintf(tmp.name, NAMESIZE, "std%d", i);
        tmp.math = rand() % 100;
        tmp.Chinese = rand() % 100;

        int ret = llist_insert(handler, &tmp, LLIST_BACKWORD__);
        if(ret){
            exit(0);                                    //如果返回值不为0，说明插入异常
        }
    }

    llist_travel(handler, print_s);

    printf("\n\n"); 

    char *del_name = "std6";
    
    int ret = llist_delete(handler, del_name, name_cmp);
    if(ret){
        printf("llist_delete failed!\n");
    }

    llist_travel(handler, print_s);

#if 0
    int id = 30;

    int ret = llist_delete(handler, &id, id_cmp);
    if(ret){
        printf("llist_delete failed!\n");
    }

    llist_travel(handler, print_s);
#endif

#if 0
    
    struct score_st *data; 
    data = llist_find(handler, &id, id_cmp);
    if(data == NULL){
        printf("Can not find!\n");
    }
    else{
        print_s(data);
    }
#endif

    llist_destroy(handler);

    exit(0);
}