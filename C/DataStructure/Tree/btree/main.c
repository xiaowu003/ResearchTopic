#include <stdio.h>
#include <stdlib.h>

#include "queue.h"

#define NAMESIZE 32

struct score_st{
    int id;
    char name[NAMESIZE];
    int math;
    int chinese;
};

struct node_st_mian{
    struct score_st data;
    struct node_st_mian *l, *r;
};

void print_s(struct score_st *d){
    printf("%d %s %d %d\n", d->id, d->name, d->math, d->chinese);
}

//插入函数
int insert(struct node_st_mian **root, struct score_st *data){
    struct node_st_mian *node;

    if(*root == NULL){
        node = malloc(sizeof(*node));
        if(node == NULL){
            return -1;
        }

        node->data = *data;
        node->l = NULL;
        node->r = NULL;

        *root = node;

        return 0;
    }

    if(data->id <= (*root)->data.id){
        return insert(&(*root)->l, data);
    }
    
    return insert(&(*root)->r, data);
}

struct score_st *find(struct node_st_mian *root, int id){
    if(root == NULL){
        return NULL;
    }

    if(id == root->data.id){
        return &root->data;
    }

    if(id < root->data.id){
        return find(root->l, id);
    }
    else{
        return find(root->r, id);
    }
}

void draw_(struct  node_st_mian *root, int level){
    if(root == NULL){
        return;
    }

    draw_(root->r, level + 1);

    for(int i = 0; i < level; i++){
        printf("    ");
    }

    print_s(&root->data);

    draw_(root->l, level + 1);
}

void draw(struct node_st_mian *root){
    draw_(root, 0);
}

static int get_num(struct node_st_mian *root){
    if(root == NULL){
        return 0;
    }

    return get_num(root->l) + 1 + get_num(root->r);
}

//找到当前节点的最左侧叶子
static struct node_st_mian * find_min(struct node_st_mian *root){
    if(root->l == NULL){
        return root;
    }

    return find_min(root->l);
}

//找到当前节点的最右侧叶子
static struct node_st_mian * find_max(struct node_st_mian *root){
    if(root->r == NULL){
        return root;
    }

    return find_max(root->r);
}

void turn_left(struct node_st_mian **root){
    struct node_st_mian *cur = *root;

    *root = cur->r;
    cur->r = NULL;

    find_min(*root)->l = cur;
}

void turn_right(struct node_st_mian **root){
    struct node_st_mian *cur = *root;

    *root = cur->l;
    cur->l = NULL;

    find_max(*root)->r = cur;
}

void balance(struct node_st_mian **root){
    int sub;

    if(*root == NULL){
        return;
    }

    while(1){
        sub = get_num((*root)->l) - get_num((*root)->r);

        if(sub >= -1 && sub <= 1){
            break;
        }

        if(sub < -1){
            turn_left(root);
        }
        else{
            turn_right(root);
        }
    }

    balance(&(*root)->l);
    balance(&(*root)->r);
}

void delete(struct node_st_mian **root, int id){
    //选择用左边的 树 顶上要删除的节点位置
    struct node_st_mian **node = root;
    struct node_st_mian *cur = NULL;

    while(*node != NULL && (*node)->data.id != id){
        if(id < (*node)->data.id){
            node = &(*node)->l;
        }
        else{
            node = &(*node)->r;
        }
    }

    if(*node == NULL){
        return;
    }

    cur = *node;

    if(cur->l == NULL){
        *node = cur->r;
    }
    else{
        *node = cur->l;
        find_max(cur->l)->r = cur->r;
    }
    
    free(cur);
}

//先序遍历，根左右
void travel_xianxu(struct node_st_mian *root){
    if(root == NULL){
        return;
    }

    print_s(&root->data);

    travel_xianxu(root->l);
    travel_xianxu(root->r);
}

//后序遍历，左右根
void travel_houxu(struct node_st_mian *root){
    if(root == NULL){
        return;
    }

    travel_houxu(root->l);
    travel_houxu(root->r);
    print_s(&root->data);
}

//TODO 按层遍历还未完成，有bug，可能是依赖的头文件“queue.h”版本不对，也许需要结合链表队列的文件
/*
//按层遍历
void travel_level(struct node_st_mian *root){
    int ret;
    queue *qu;
    struct node_st_mian *cur;

    qu = qu_creat(sizeof(struct node_st_mian *));
    if(qu == NULL){
        return;
    }

    qu_enqueue(qu, &root);

    while(1){
        ret = qu_dequeue(qu, &cur);
        if(ret == -1){
            break;
        }

        print_s(&cur->data);

        if(cur->l != NULL){
            qu_enqueue(qu, &cur->l);
        }

        if(cur->r != NULL){
            qu_enqueue(qu, &cur->r);
        }
    }

    qu_destroy(qu);
}
*/

int main(){
    struct node_st_mian *tree = NULL;
    struct score_st tmp, *datap;

    int arr[] = {1, 2, 3, 7, 6, 5, 9, 8, 4};

    for(int i = 0; i < sizeof(arr)/sizeof(*arr); i++){
        tmp.id = arr[i];
        snprintf(tmp.name, NAMESIZE, "stu%d", arr[i]);
        tmp.math = rand()%100;
        tmp.chinese = rand()%100;

        insert(&tree, &tmp);
    }

    draw(tree);
    
    printf("--------------\n");

    balance(&tree);

    draw(tree);
    
#if 0
    int tmpid = 5;
    //删除节点
    delete(&tree, tmpid);

    printf("-------\n\n");
    draw(tree);
#endif

    printf("\n\n");

    travel_xianxu(tree);

    printf("\n\n");

    travel_houxu(tree);
    
#if 0
    int tmpid = 12;
    datap = find(tree, tmpid);
    if(datap == NULL){
        printf("Can not find the id %d\n", tmpid);
    }
    else{
        print_s(datap);
    }
#endif

    exit(0);
}