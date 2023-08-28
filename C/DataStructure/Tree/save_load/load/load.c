//DONE 将广义表转换为树
#include <stdio.h>
#include <stdlib.h>

#define FNAME "E:\\VScode\\C\\linux_c\\Tree\\save_load\\load\\out.txt"

struct node_st{
    char data;
    struct node_st *l, *r;
};

void draw_(struct  node_st *root, int level){
    if(root == NULL){
        return;
    }

    draw_(root->r, level + 1);

    for(int i = 0; i < level; i++){
        printf("    ");
    }

    printf("%c\n", root->data);

    draw_(root->l, level + 1);
}

void draw(struct node_st *root){
    draw_(root, 0);
}

static int get_num(struct node_st *root){
    if(root == NULL){
        return 0;
    }

    return get_num(root->l) + 1 + get_num(root->r);
}

//找到当前节点的最左侧叶子
static struct node_st * find_min(struct node_st *root){
    if(root->l == NULL){
        return root;
    }

    return find_min(root->l);
}

struct node_st *load_(FILE *fp){
    int c;
    struct node_st *root;

    c = fgetc(fp);
    if(c != '('){
        fprintf(stderr, "fgetc():error.\n");
        exit(1);
    }

    c = fgetc(fp);
    if(c == ')'){
        return NULL;
    }

    root = malloc(sizeof(*root));
    if(root == NULL){
        exit(1);
    }

    root->data = c;
    root->l = load_(fp);
    root->r = load_(fp);

    fgetc(fp);

    return root;
}

struct node_st *load(const char *path){
    FILE *fp;
    struct node_st *root;
    
    fp = fopen(path, "r");
    if(fp == NULL){
        return NULL;
    }

    root = load_(fp);

    fclose(fp);

    return root;
}

int main(){
    struct node_st *root;

    root = load(FNAME);

    draw(root);

    exit(0);
}