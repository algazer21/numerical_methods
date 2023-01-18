/*********************************************
Alan García Zermeño
28/8/2022
Resuelve Ax=b por el método de LDU.
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
    xx[0][0] = bx[0][0]/Ax[0][0];
    for(int i = 1;i<rows;i++){
        for (int j=0;j<i;j++){
            su += Ax[i][j]*xx[j][0];
        }
        xx[i][0] = (bx[i][0]-su)/Ax[i][i];
        su = 0.0;
    }
}
/********************************************************************/
void sol_tri_sup(size_t rows, size_t cols, double (*bx)[1], double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j;
    double su=0;
    xx[rows-1][0] = bx[rows-1][0]/Ax[rows-1][rows-1]; //trivial case//remove sums
    
    for(int i =rows-2;i>=0;i--){
        //compute sum
        for (int j=i+1;j<rows;j++){
            su += Ax[i][j]*xx[j][0];
        }
        xx[i][0] = (bx[i][0]-su)/Ax[i][i];
        su = 0.0;
    }
}
/****************************************************************************/
void sol_tri_diag(size_t rows, size_t cols, double (*bx)[1], double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j;
    for(int i=0;i<rows;i++){
        xx[i][0] = bx[i][0]/Ax[i][i];
    }
}
/********************************************************************/
void ldu(size_t rows, size_t cols,double (*bx)[1],double (*xx)[1],double (*Ax)[cols])
{
    int i,j,p;
    double a;

    for (j=0;j<rows;j++) {
        a = 0.0;
        for( p=0;p<j;p++){
            a += Ax[p][p]*Ax[j][p]*Ax[p][j];
        }
        Ax[j][j] = Ax[j][j] - a;
        
        for (i =j+1; i < rows; i++) {
            a = 0;
            for ( p = 0; p < j; p++) {
                a += Ax[p][p]*Ax[i][p]*Ax[p][j];  
            }
            Ax[i][j] = (Ax[i][j] - a)/Ax[j][j];
        }
        for ( i = j+1; i < rows; i++) {
            a = 0;
            for( p = 0; p < j; p++) {
                a += Ax[p][p] * Ax[j][p] * Ax[p][i];
            }
            Ax[j][i] = (Ax[j][i] - a) / Ax[j][j];
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
    double x[m][1],time,b[m][1],y[m][1],z[m][1];
    double (*A)[m] = calloc(n, sizeof *A);
    double (*D)[m] = calloc(n, sizeof *D);


    readmatrix(n, m, A, "A.txt",1);
    readmatrix(n, 1, b, "b.txt",0);
    printf("\n");
    ldu(n,m,b,x,A);
    for(i=0;i<n;i++){
        D[i][i] = A[i][i];
        A[i][i] = 1.0;
    }

    sol_tri_inf(n,m,b,z,A);
    sol_tri_diag(n,m,z,y,D);
    sol_tri_sup(n,m,y,x,A); 
    writematrix(n,1,x);

    t = clock()-t;                                      //Tiempo
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}
