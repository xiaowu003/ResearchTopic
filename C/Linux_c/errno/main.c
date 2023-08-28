#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

int main(){
    FILE *fp;

    //BUG   
    //这里如果使用open()函数打开文件，就会输出OK！，open()函数操作成功返回一个文件描述符，失败返回-1
    //这里的文件描述符是int型,所以不会进入if()分支
    fp = fopen("E:\\VScode\\C\\Linux_c\\errno\\tmp.txt", "r");
    if(fp == NULL){
        //第一种报错方式,使用fprintf()函数手动输出错误信息
        fprintf(stderr, "fopen() failed! errno = %d\n", errno);

        //第二种报错方式，使用函数perror()
        perror("fopen()");

        //第三种方式，使用函数strerror()以及函数fprintf()
        fprintf(stderr, "fopen():%s\n", strerror(errno));
        //此处记得要包含string.h头文件

        exit(1);
    }

    puts("OK!");

    fclose(fp);

    exit(0);
}