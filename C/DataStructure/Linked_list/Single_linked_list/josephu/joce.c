#include <stdio.h>
#include <stdlib.h>

#define JOSE_NR 8

struct node_st {
    int data;
    struct node_st *next;
};

struct node_st *jose_create(int n){
    struct node_st *me;
    struct node_st *newnode, *cur;
    int i = 1;

    me = malloc(sizeof(*me));

    if(me == NULL){
        return NULL;
    }

    me->data = i;
    me->next = me;
    i++;

    cur = me;

    for(; i <= n; i++){
        newnode = malloc(sizeof(*newnode));

        if(newnode == NULL){
            return NULL;
        }

        newnode->data = i;
        newnode->next = me;
        cur->next = newnode;
        cur = newnode;
    }

    return me;
}

void jose_show(struct noed_st *me){
    struct node_st *list;

    for(list = me; list->next != me; list = list->next){
        printf("%d ",list->data);
    }

    printf("%d\n",list->data);
}

void jose_kill(struct node_st **me, int n){
    struct node_st *cur = *me, *node;
    int  i = 1;

    while(cur != cur->next){
        while (i < n){
            node = cur;
            cur = cur->next;
            i++;
        }

        printf("%d ", cur->data);

        node->next = cur->next;

        free(cur);

        cur = node->next;
        i = 1;
    }

    *me = cur;

    printf("\n");
}

int main(){
    struct node_st *list;

    list = jose_create(JOSE_NR);

    jose_show(list);

    jose_kill(&list, 3);

    jose_show(list);

    exit(0);
}