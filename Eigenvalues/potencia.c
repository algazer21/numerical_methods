/*********************************************
Alan García Zermeño
8/9/2022
Encuentra el valor y vector propio más grandes de A
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void normalize(size_t rows, double *x);
void EigenPower(size_t rows, size_t cols, double **a, double *y, double *x, double tol);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
double prod_vec_vec(double *b,double *c, size_t rows);
void freel(int n,double **listm);

int main()
{
    int n,m,i;
    clock_t t = clock();
    double time;
    FILE    *fptr;
    fptr = fopen("As.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    double **A = (double **)malloc(n*sizeof(double*));
    double *y = (double *)malloc(n*sizeof(double));
    double *x = (double *)malloc(n*sizeof(double));

    readmatrix(n, m, A, "As.txt");
    EigenPower(n,m,A,y,x,10e-7);   
    //writematrix(n,m,A,"x.txt");
    freel(n,A);
    free(x);free(y);
 
    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("\n\nTiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}

/*********************************************************************/
void normalize(size_t rows, double *x){
    double sum=0;
    for (int i=0;i<rows;i++){
        sum+= x[i]*x[i];}
    sum = sqrt(sum);
    for (int i=0;i<rows;i++){
        x[i] = x[i]/sum;}
}

/*********************************************************************/
void EigenPower(size_t rows, size_t cols, double **a, double *y, double *x ,double tol){
    int i,j=0,nmax = 1000;
    double val=0,err = 1.0,valn;
    for(i=0;i<rows;i++){
        y[i] = 1.0;
    }
    normalize(rows,y);
    while(err>tol && (j<nmax)){
        j+=1;
        prod_matr_vec(a,y,x,rows);
        valn = prod_vec_vec(y,x,rows);
        err = fabs(valn-val);
        normalize(rows,x);

        val = valn;
        for (int i = 0;i<rows;i++){
            y[i] = x[i];}
    }
    printf("Valor propio más grande: %.5lf calculado en %d iteraciones \n\n",valn,j);
    printf("Vector propio más grande normalizado:\n");
    for (i = 0;i<rows;i++){printf("%.3lf  ",x[i]);}

}

/********************************************************************/
void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

/********************************************************************/
void readmatrix(size_t rows, size_t cols, double **a, const char* filename)
{
    int n,m;
    FILE *pf;
    pf = fopen (filename, "r");
    if (pf==NULL){
        printf("No se encontro archivo\n");
        exit(1);}
    fscanf(pf, "%d %d", &n,&m);
    for(size_t i = 0; i < rows; ++i){
        a[i] = (double *)malloc(cols*sizeof(double));
        for(size_t j = 0; j < cols; ++j)
            fscanf(pf, "%lf", a[i] + j);
    }
    fclose(pf);
}

/********************************************************************/
void writematrix(size_t rows,size_t cols,double **a,const char* filename)
{
    int i,j;
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            fprintf(fptr,"%.1lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}
/******************************************************************/
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows){
    for (int i = 0;i<rows;i++){
        v1[i] = prod_vec_vec(a[i],v0,rows);
    }
}

/******************************************************************/
double prod_vec_vec(double *b,double *c, size_t rows){
    double a= 0;
    for (int i = 0;i<rows;i++){
        a += b[i]*c[i];
    }
    return a;
}
