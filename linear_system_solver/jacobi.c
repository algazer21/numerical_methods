/*********************************************
Alan García Zermeño
Para el curso de métodos numéricos.
CIMAT 3/8/2022
Resuelve Ax=b por el método de Jacobi.
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double utv(double x){
    return (x*x+x);
}

void normae(int n,double (*x1)[1]){
    int i;
    double sum = 0,sum2=0,p,t;
    p= 1.0/(n+1);
    for (i = 0; i <n; i++)
    {   
        t = utv(p*i);
        sum+= (x1[i][0]-t)*(x1[i][0]-t);
        sum2+= t*t;
    }
    printf("Error calculado con norma euclidiana: %.9lf\n",sqrt(sum)/sqrt(sum2));
}

void readmatrix(size_t rows, size_t cols, double (*a)[cols], const char* filename, int flag)
{

    int n,m;
    FILE *pf;
    pf = fopen (filename, "r");
    if(flag==1){
        fscanf(pf, "%d %d", &n,&m);
    }

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            fscanf(pf, "%lf", a[i] + j);
    }
    fclose(pf);
}

/********************************************************************/
void jacobi(size_t rows,double (*vn)[1],double (*vv)[1],double (*bx)[1],double (*Ax)[rows])
{
    int i,j,k,it=0,iter = 5000;
    double num,deno,sum,tol,valor = 1e-5;
    for(i=0;i<rows;i++){
        vv[i][0] = bx[i][0]/Ax[i][i];
    }

    while(it <iter){
        num = 0;
        deno = 0;
        for (i=0;i<rows;i++){
            sum = 0;
            for (j=0;j<rows;j++){
                if(i!=j){
                    sum += Ax[i][j]*vv[j][0];
                }
            }
            vn[i][0] = (bx[i][0]-sum)/Ax[i][i];
            num += (vv[i][0]-vn[i][0])*(vv[i][0]-vn[i][0]);
            deno += (vn[i][0]*vn[i][0]);
            tol = num/deno;
            if(tol<valor){
                printf("Convergencia\n");
                it = iter;
                break;}

            for (k=0;k<rows;k++){
                vv[i][0] = vn[i][0];
            }
            it +=1;
        }
    }
    
}
/********************************************************************/
void writematrix(size_t rows,size_t cols,double (*a)[cols],const char* filename)
{
    int i,j;
    FILE    *fptr;
    fptr = fopen(filename, "w");

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            fprintf(fptr,"%lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

/*********************************************************************/

int main()
{
    int n,m,i;
    clock_t t;
    double time;
    FILE    *fptr;
    fptr = fopen("A.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    double (*A)[m] = calloc(n, sizeof *A);
    double (*L)[m] = calloc(n, sizeof *L);
    double (*b)[1] = calloc(n, sizeof *b);
    double (*Vv)[1] = calloc(n, sizeof *Vv);
    double (*x)[1] = calloc(n, sizeof *x);

    readmatrix(n, m, A, "A.txt",1);
    readmatrix(n, 1, b, "b.txt",0);
    t = clock();
    jacobi(n,x,Vv,b,A); 
    t = clock()-t;

    writematrix(n,1,Vv,"x.txt");
    normae(n,Vv);
 
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}