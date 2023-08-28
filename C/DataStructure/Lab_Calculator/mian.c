#include <stdio.h>
#include <stdlib.h>

#include "sqstack.h"

int compute(sqstack *snum, datatype *op){
    datatype n1, n2, n;
    st_pop(snum, &n2);
    st_pop(snum, &n1);

    //n2为先出栈的元素， n1为后出栈的元素
    switch(*op){
        case '+':
            n = n1 + n2;
            break;
        case '-':
            n = n1 - n2;
            break;
        case '*':
            n = n1 * n2;
            break;
        case '/':
            n = n1 / n2;
            break;
        default:
            exit(1);
    }

    //将小阶段的计算结果保存起来，入栈
    st_push(snum, &n);
}

void deal_bracket(sqstack *snum, sqstack *sop){
    datatype old_op;

    //取栈顶元素
    st_top(sop, &old_op);

    //判断是不是（
    while(old_op != '('){
        //出栈
        st_pop(sop, &old_op);

        //进行计算
        compute(snum, &old_op);

        //继续取栈顶元素
        st_top(sop, &old_op);
    }
    st_pop(sop, &old_op);
}

static int get_pri(int op){
    switch(op){
        case '(':
            return 0;
        case '+':
        case '-':
            return 1;
        case '*':
        case '/':
            return 2;
    }
}

static void deal_op(sqstack *snum, sqstack *sop, int op){
    datatype old_op;

    if(st_isempty(sop) || op == '('){
        st_push(sop, &op);
        return ;
    }

    st_top(sop, &old_op);

    if(get_pri(op) > get_pri(old_op)){
        st_push(sop, &op);
        return;
    }

    while(get_pri(op) <= get_pri(old_op)){
        st_pop(sop, &old_op);
        compute(snum, &old_op);

        if(st_isempty(sop)){
            break;
        }

        st_top(sop, &old_op);
    }

    st_push(sop, &op);
}

int main(){
    int i = 0, value = 0, flag = 0;
    char str[] = "(11+3)*2-5";

    //创建两个栈，一个存放数字，一个存放运算符
    sqstack *snum, *sop;

    snum = st_create();
    if(snum == NULL){
        exit(1);
    }

    sop = st_create();
    if(sop == NULL){
        exit(1);
    }

    while(str[i] != '\0'){
        if(str[i] > '0' && str[i] < '9'){
            //is a num
            value = value * 10 + (str[i] - '0');
            flag = 1;
        }
        else{
            //is a op

            //先判断是否需要入栈
            if(flag){
                st_push(snum, &value);
                flag = 0;
                value = 0;
            }

            if(str[i] == ')'){
                /*
                ')'右括号比较特别，不需要入栈，
                只要出现就表示需要计算与前面左括号之间的算式了
                */
               deal_bracket(snum, sop);
           }
           else{
            //（ + - * /
                deal_op(snum, sop, str[i]);
           }
        }

        i++;
    }

    if(flag){
        st_push(snum, &value);
    }

    datatype old_op;
    while(!st_isempty(sop)){
        st_pop(sop, &old_op);
        compute(snum, &old_op);
    }

    st_travel(snum);

    st_destroy(snum);
    st_destroy(sop);

    exit(0);
}