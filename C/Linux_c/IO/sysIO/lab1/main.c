#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define BUFSIZE 1024

int main(int argc, char **argv){
    int sfd, dfd;
    char buf[BUFSIZE];
    int len, ret, pos;

    if(argc < 3){
        fprintf(stderr, "Usage...\n");
        exit(1);
    }

    sfd = open(argv[1], O_RDONLY);
    if(sfd  < 0){
        perror("open(1)");
    }

    //BUG 代码应该是在目标文件存在的时候清空，但是实际运行并没有清空，只是将前len个字符覆盖了
    dfd = open(argv[2], O_WRONLY|O_CREAT,O_TRUNC, 0600);
    if(dfd < 0){
        close(sfd);
        perror("open(2)");
        exit(1);
    }


    while(1){
        len = read(sfd, buf, BUFSIZ);
        if(len < 0){
            perror("read(1)");
            break;
        }
        if(len == 0){
            break;
        }

        pos = 0;

        while(len > 0){
            ret = write(dfd, buf + pos, len);
            if(ret < 0){
                perror("write(1)");
                exit(1);
            }

            pos += len;
            len -= ret;
        }
        
    }

    close(dfd);
    close(sfd);

    exit(1);
}