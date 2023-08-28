#include <stdio.h>
#include <stdlib.h>

#include "queue.h"
#include "sqstack.h"

#define NR_BALL 27

//检查qu队列是否回到一开始的从1到27顺序排列的样子
static int check(queue *qu){
    int i = (qu->head + 1) % MAXSIZE;

    while(i != qu->tail){
        if(qu->data[i] > qu->data[(i + 1) % MAXSIZE]){
            //从小到大排列，不会出现这样的结果，若出现了说明还未变回原来的序列
            return 0;
        }

        i = (i + 1) % MAXSIZE; 
    }
    return 1;
}

int main(){
    queue *qu;                                  //队列的指针
    sqstack *st_min, *st_fivemin, *st_hour;     //栈的指针

    int t, value;
    int time = 0;

    qu = qu_creat();
    if(qu == NULL){
        exit(1);
    }

    st_min = st_create();
    if(st_min == NULL){
        exit(1);
    }

    st_fivemin = st_create();
    if(st_fivemin == NULL){
        exit(1);
    }

    st_hour = st_create();
    if(st_hour == NULL){
        exit(1);
    }

    for(int i = 1; i <= NR_BALL; i++){
        qu_enqueue(qu, &i);                     //将27个球入队            
    }

    qu_travel(qu);

    while(1){
        qu_dequeue(qu, &t);                     //出队
        time++;

        if(st_min->top != 3){
            //st_min未满，继续入栈
            st_push(st_min, &t);                
        }
        else{
            //st_min已满，时间进一，且将st_min全部出栈，并原来的对列
            while(!st_isempty(st_min)){
                st_pop(st_min, &value);         //出栈
                qu_enqueue(qu, &value);         //入队
            }

            if(st_fivemin->top != 10){
                //st_min栈已满的情况下st_fivemin栈还未满，继续入栈
                st_push(st_fivemin, &t);
            }
            else{
                //st_fivemin栈已满
                while(!st_isempty(st_fivemin)){
                    st_pop(st_fivemin, &value); //出栈
                    qu_enqueue(qu, &value);     //入队
                }

                if(st_hour->top != 10){
                    //st_five栈已满但st_hour栈未满，
                    st_push(st_hour, &t);
                }
                else{
                    //非空时进行出栈
                    while(!st_isempty(st_hour)){
                        st_pop(st_hour, &value);//出栈
                        qu_enqueue(qu, &value); //入队
                    }

                    qu_enqueue(qu, &t);

                    if(check(qu)){
                        break;
                    }
                }
            }
        }
    }

    printf("time = %d\n",time);
    qu_travel(qu);

    qu_destroy(qu);
    st_destroy(st_min);
    st_destroy(st_fivemin);
    st_destroy(st_hour);

    exit(0);
}