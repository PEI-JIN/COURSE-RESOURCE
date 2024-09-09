// Homework - Convolutional Encoder (Example)
// 作者: 電通所 Q36114221 蘇沛錦
// 日期: 2023/04/27 (星期四) 

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void)
{
    int i,j;    // Loop counter
    int L;
    int m = 6;

    FILE * fp;
    if((fp=fopen("Sim.txt","r"))==NULL){
        printf("!!! Input File Does Not Exist !!!\n");
        system("pause");
        exit(0);
    }
    else{
        fscanf(fp, "%d", &L);
        fclose(fp);
    }

    int * x = (int*)calloc(2*(L+m),sizeof(int));
    int * u = (int*)calloc(L+m,sizeof(int));
    int * s = (int*)calloc(m,sizeof(int));

    u[0] = 1;
    for(i=1;i<L;i++){
        if(i<6){
            u[i] = 0;
        }
        else{
            u[i] = u[i-5]^u[i-6];
        }
    }

    fp = fopen("x.txt","w");
    printf("Output x.txt:\n");
    for(i=0;i<(L+m);i++){
        x[2*i] = u[i]^s[1]^s[2]^s[4]^s[5];
        x[2*i+1] = u[i]^s[0]^s[1]^s[2]^s[5];

        x[2*i] = pow(-1,x[2*i]);
        x[2*i+1] = pow(-1,x[2*i+1]);

        fprintf(fp,"%d %d ",x[2*i],x[2*i+1]);
        printf("%d %d ",x[2*i],x[2*i+1]);

        for(j=5;j>0;j--){
            s[j] = s[j-1];
        }
        s[0] = u[i];
    }
    fprintf(fp,"%%2*(L+m) elements");
    printf("%%2*(L+m) elements\n");
    fclose(fp);

    free(x);
    free(u);
    free(s);

    system("pause");
    return 0;
}
