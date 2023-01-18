/*********************************************
Alan García Zermeño
Para el curso de métodos numéricos.
CIMAT 26/8/2022
Resuelve Ax=b por el método de Crout.
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
}

/********************************************************************/
void sol_tri_inf(size_t rows, size_t cols, double (*bx)[1], double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j;
    double su=0;
    for(i=0;i<rows;i++){
        for(j=0;j<i;j++){
            su+=Ax[i][j]*xx[j][0];
        }
        xx[i][0] = (bx[i][0]-su)/Ax[i][i];
        su = 0;
    }
}
/********************************************************************/
void sol_tri_sup(size_t rows, size_t cols, double (*bx)[1], double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j;
    double su=0;
    for(i=0;i<rows;i++){
        Ax[i][i] = 1.0;
    }

    for(i=(rows-1);i>=0;i--){
        for(j=i+1;j<cols;j++){
            su+=Ax[i][j]*xx[j][0];
        }

        xx[i][0] = (bx[i][0]-su)/Ax[i][i];
        su = 0;
    }
}
/********************************************************************/
void crout(size_t rows, size_t cols,double (*bx)[1],double (*xx)[1],double (*Ax)[cols])
{
    int i,j,p;
    double a;

    for(j=0;j<rows;j++){
        for(i=j; i<rows; i++){
            a = 0;
            for(p=0;p<j;p++){
                a = a+Ax[i][p]*Ax[p][j];   
            }
            Ax[i][j] = Ax[i][j]-a;

        }
        for(i=j+1;i<rows;i++){
            a = 0;
            for(p=0;p<j;p++){
                a= a +Ax[j][p]*Ax[p][i];
            }
            Ax[j][i] = (Ax[j][i]-a)/Ax[j][j];
        }
    }
}

/********************************************************************/
void writematrix(size_t rows, size_t cols, double (*a)[cols])
{
    int i,j;
    FILE    *fptr;
    fptr = fopen("x.txt", "w");

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            fprintf(fptr,"x%d = %.2lf ",i,a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}
/*********************************************************************/

int main()
{
    int n,m,i,j;
    clock_t t = clock();
    FILE    *fptr;
    fptr = fopen("A.txt", "r");                        //Leer tamaño de matriz
    fscanf(fptr, "%d %d", &n,&m);
    double x[m][1],time,b[m][1],y[m][1];
    double (*A)[m] = calloc(n, sizeof *A);


    readmatrix(n, m, A, "A.txt",1);
    readmatrix(n, 1, b, "b.txt",0);
    printf("\n");
    crout(n,m,b,x,A);
    sol_tri_inf(n,m,b,y,A);
    sol_tri_sup(n,m,y,x,A);
    
    writematrix(n,1,x);

    t = clock()-t;                                      //Tiempo
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}