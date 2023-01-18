/*********************************************
Alan García Zermeño
Para el curso de métodos numéricos.
CIMAT 26/8/2022
Resuelve Ax=b por el método de Gauss Jordan.
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
void gauss(size_t rows, size_t cols,double (*bx)[1]  ,double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j,k;
    double c,sum=0.0;

    for(k=0;k<(rows-1);k++){
        for(i=k; i<(rows-1); i++){
            c = Ax[i+1][k]/Ax[k][k];
            for(j=k;j<cols;j++){
                Ax[i+1][j] = Ax[i+1][j]-Ax[k][j]*c;   
            }
            bx[i+1][0] = bx[i+1][0]-bx[k][0]*c;
        }
    }
    for(i=(rows-1);i>=0;i--){
        for(j=i+1;j<cols;j++){
            sum+=Ax[i][j]*xx[j][0];
        }
        xx[i][0] = (bx[i][0]-sum)/Ax[i][i];
        sum = 0;
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
    double x[m][1],time,b[m][1];
    double (*A)[m] = calloc(n, sizeof *A);

    readmatrix(n, m, A, "A.txt",1);
    readmatrix(n, 1, b, "b.txt",0);
    gauss(n,m,b,x,A);
    
    writematrix(n,1,x);

    t = clock()-t;                                      //Tiempo
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}