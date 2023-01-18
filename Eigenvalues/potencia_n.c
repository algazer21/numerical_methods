/*********************************************
Alan García Zermeño
8/9/2022
Encuentra los n valores propios más grandes de A
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 10                                                //Numero de valores propios a calcular

void normalize(size_t rows, double *x);
void writevector(size_t rows,size_t cols,double *a,const char* filename);
void EigenPower_n(size_t rows, size_t cols, double **a, double *y, double *v1, double tol,double **valp,int ii,double *vv);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
double prod_vec_vec(double *b,double *c, size_t rows);
void freel(int n,double **listm);

int main()
{
    int n,m,j=0,i,k;
    clock_t t = clock();
    double time;
    FILE    *fptr;
    fptr = fopen("AA.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    double **A = (double **)malloc(n*sizeof(double*));
    double *y = (double *)malloc(n*sizeof(double));
    double **x = (double **)malloc(n*sizeof(double*));
    double *x2 = (double *)malloc(n*sizeof(double));
    double *vp = (double *)malloc(n*sizeof(double));
    for (i = 0;i<n;i++){
        x[i] = (double *)malloc(n*sizeof(double));   
        for (k = 0;k<m;k++){
            x[i][k] = 0.0;}
        }

    readmatrix(n, m, A, "AA.txt");
    for(k=0;k<N;k++){
        EigenPower_n(n,m,A,y,x2,10e-6,x,k,vp);
        printf("Eigenvalor %d = %lf... \n",k+1,vp[k]);}

    writematrix(N,n,x,"vecM.txt");
    writevector(n,1,vp,"valM.txt");
    //TeoVal(vp,x,n);
    freel(n,A);
    freel(n,x);
    free(y);free(vp);
 
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
            fprintf(fptr,"%.6lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

/********************************************************************/
void writevector(size_t rows,size_t cols,double *a,const char* filename)
{
    int i;
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    for(int i=0;i<rows;i++){
            fprintf(fptr,"%.7lf\n",a[i]);
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
/*******************************************************************/
void EigenPower_n(size_t rows, size_t cols, double **a, double *y, double *v1 ,double tol,double **valp,int ii,double *vv){
    int i,j=0,nmax = 1000,k,kk;
    double val=0,err = 1.0,valn,an = 0;
    for(i=0;i<rows;i++){
        y[i] = rand()%3+1;
        v1[i] = 0.0;
    }
    normalize(rows,y);

    while(err>tol && (j<nmax)){
        //Para quitar contribución de los eigenvectores previos.
        if(ii>0){
        for(i=0;i<ii;i++){
            an = prod_vec_vec(valp[i],y,rows);
            for(kk=0;kk<rows;kk++){
                y[kk] -= (valp[i][kk]*an);}}}
        
        j+=1;
        prod_matr_vec(a,y,v1,rows);
        valn = prod_vec_vec(y,v1,rows);
        err = fabs(valn-val);
        normalize(rows,v1);

        val = valn;
        for (int i = 0;i<rows;i++){
            y[i] = v1[i];}
    }
    vv[ii] = valn;
    //printf("\n\nValor propio número %d: %.5lf calculado en %d iteraciones \n",ii+1,valn,j);
    //printf("Vector propio número %d:  ",ii+1);
    for (i = 0;i<rows;i++){
        //printf("%.3lf  ",v1[i]);
        valp[ii][i] = v1[i];
    }
    j+=1;

}
