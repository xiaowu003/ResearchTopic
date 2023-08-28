//树转存为广义表

#include <stdio.h>
#include <stdlib.h>

#define FNAME "E:\\VScode\\C\\linux_c\\Tree\\save_load\\save\\out.txt"

struct node_st_mian{
    char data;
    struct node_st_mian *l, *r;
};

static struct node_st_mian *tree = NULL;

//插入函数
int insert(struct node_st_mian **root, int data){
    struct node_st_mian *node;

    if(*root == NULL){
        node = malloc(sizeof(*node));
        if(node == NULL){
            return -1;
        }

        node->data = data;
        node->l = NULL;
        node->r = NULL;

        *root = node;

        return 0;
    }

    if(data <= (*root)->data){
        return insert(&(*root)->l, data);
    }
    
    return insert(&(*root)->r, data);
}

void draw_(struct  node_st_mian *root, int level){
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

void draw(struct node_st_mian *root){
    draw_(root, 0);
}

int save_(struct node_st_mian *root, FILE *fp){
    fputc('(', fp);
    /*if error*/

    if(root == NULL){
        fputc(')', fp);
        /*if error*/
        return 0;
    }
    
    fputc(root->data, fp);

    save_(root->l, fp);
    save_(root->r, fp);

    fputc(')', fp);
    /*if error*/

    return 0;
}

int save(struct node_st *root, const char *path){
    FILE *fp;
    
    fp = fopen(path, "w");
    if(fp == NULL){
        return -1;
    }

    save_(tree, fp);

    fclose(fp);
    return 0;
}

int main(){
    char arr[] = "cefadjbh";

    for(int i = 0; i < sizeof(arr)/sizeof(*arr) - 1; i++){
        insert(&tree, arr[i]);
    }

    draw(tree);

    save(tree, FNAME);

    exit(0);
}