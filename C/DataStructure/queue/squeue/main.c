#include <stdio.h>
#include <stdlib.h>

#include "queue.h"

int main(){
    queue *sq;
    datatype arr[] = {2, 34, 98, 12};

    sq = qu_creat();
    if(sq == NULL){
        exit(1);
    }
    for(int i = 0; i < sizeof(arr)/sizeof(*arr); i++){
        qu_enqueue(sq, &arr[i]);
    }

    qu_travel(sq);

/*
    datatype tmp = 100;
    int ret = qu_enqueue(sq, &tmp);
    if(ret == -1){
        printf("queue is full!\n");
    }
    else{
        qu_travel(sq);
    }
*/

    datatype tmp;
    qu_dequeue(sq,&tmp);
    printf("DEQUEUE:%d\n", tmp);

    qu_travel(sq);

    qu_destroy(sq);

    exit(0);
}