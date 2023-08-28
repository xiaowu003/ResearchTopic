#include <stdio.h>
#include <stdlib.h>

#include "queue.h"

#define NAMESIZE 32

struct score_st {
    int id;
    char name[NAMESIZE];
    int math;
    int chinese;
};

static void print_s(void *record){
    struct score_st *r = record;
    
    printf("%d %s %d %d\n", r->id, r->name, r->math, r->chinese);
}

int main(){
    QUEUE *qu;
    struct score_st tmp;

    qu = queue_creat(sizeof(struct score_st));
    if(qu == NULL){
        exit(1);
    }

    for(int i = 0; i < 6; i++){
        tmp.id = i;
        snprintf(tmp.name, NAMESIZE, "stu%d",i);
        tmp.math = rand()%100;
        tmp.chinese = rand()%100;

        if(queue_en(qu, &tmp) != 0){
            break;
        }
    }
    
    while(1){
        int ret = queue_de(qu, &tmp);
        if(ret == -1){
            break;
        }

        print_s(&tmp);
    }
    

    queue_destroy(qu);

    exit(0);
}