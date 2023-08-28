#include <stdio.h>
#include <stdlib.h>

#include "nohead.h"

int main(){

    struct node_st *list = NULL;
    struct score_st tmp;

    for(int i = 0; i < 7; i++){
        tmp.id = i;
        snprintf(tmp.name, NAMESIZE, "stu%d", i);
        tmp.math = rand() % 100;
        tmp.Chinese = rand() % 100;

        int res = list_insert(&list, &tmp);
    }

    //list = list_insert(list, &tmp);

    list_show(list);

    printf("\n\n");

    //list_delete(&list);
    //list_show(list);

    int id = 2;
    struct score_st *ptr = NULL;

    ptr = list_find(list, id);
    if(ptr == NULL){
        printf("Can not find!\n");
    }
    else{
        printf("%d %s %d %d\n", ptr->id, ptr->name, ptr->math, ptr->Chinese);
    }

    list_destroy(list);

    exit(0);
}