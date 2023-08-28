#include <stdio.h>
#include <stdlib.h>

#include "list.h"

//TODO
/*
A bug: 运行时，每次在第一个位置输出一个不合法数据
*/

int main(){

    list *l;
    datatype arr[] = {12, 9, 23, 2, 34, 6, 45};

    l = list_creat();
    if(l == NULL){
        exit(1);
    }

    for(int i = 0; i < sizeof(arr)/sizeof(*arr); i++){
        if(list_order_insert(l, &arr[i])){
            exit(1);
        }
    }

    list_display(l);

    int value = 12;
    int err = list_delete_at(l, 2, &value);

    if(err){
        exit(1);
    }
    
    list_display(l);
    printf("delete:%d\n", value);

#if 0
    int value = 12;
    list_delete(l, &value);
    list_display(l);
#endif


    list_destroy(l);
    
    exit(0);
}