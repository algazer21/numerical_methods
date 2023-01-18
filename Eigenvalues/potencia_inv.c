/*********************************************
Alan García Zermeño
8/9/2022
Encuentra el valor y vector propio más pequeño de A
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void normalize(size_t rows, double *x);
void EigenPower_inv(size_t rows, size_t cols, double **a, double *y, double *x, double tol);
void crout(size_t rows, size_t cols,double **Ax);
void sol_tri_sup(size_t rows, size_t cols, double *bx, double *xx ,double **Ax);
void sol_tri_inf(size_t rows, size_t cols, double *bx, double *xx ,double **Ax);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
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
    crout(n,m,A);
    EigenPower_inv(n,m,A,y,x,10e-9);   
    //writematrix(n,m,x,"x.txt");

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
void EigenPower_inv(size_t rows, size_t cols, double **a, double *y, double *x ,double tol){
    int i,j=0,nmax = 1000;
    double lamO=10e8,err = 1.0,lam,sum,val=1;
    double *t = (double *)malloc(rows*sizeof(double));

    for(i=0;i<rows;i++){
        y[i] = 1.0;
    }
    normalize(rows,y);

    while((val>tol) && (j<nmax)){
        j+=1;
        sol_tri_inf(rows,rows,y,t,a);
        sol_tri_sup(rows,rows,t,x,a);
        sum = 0;
        for (i=0;i<rows;i++){
            sum+= y[i]*x[i];}

        lam = 1.0/sum;
        val = fabs(lam-lamO);
        normalize(rows,x);

        for (int i = 0;i<rows;i++){
            y[i] = x[i];}
        lamO = lam;
    }
    printf("Valor propio más pequeño: %.5lf calculado en %d iteraciones \n\n",lam,j);
    printf("Vector propio más pequeño normalizado:\n");
    for (i = 0;i<rows;i++){printf("%.6lf  ",x[i]);}
    free(t);
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

void crout(size_t rows, size_t cols,double **Ax){
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
/******************************************************************/
void sol_tri_inf(size_t rows, size_t cols, double *bx, double *xx ,double **Ax)
{
    int i,j;
    double su;
    for(i=0;i<rows;i++){
        su = 0;
        for(j=0;j<i;j++){
            su+=Ax[i][j]*xx[j];
        }
        xx[i] = (bx[i]-su)/Ax[i][i];
    }
}
/********************************************************************/
void sol_tri_sup(size_t rows, size_t cols, double *bx, double *xx ,double **Ax)
{
    int i,j;
    double su=0;
    for(i=(rows-1);i>=0;i--){
        for(j=i+1;j<cols;j++){
            if(i==j){su+=xx[j];}
            else{su+=Ax[i][j]*xx[j];}
        }
        xx[i] = (bx[i]-su)/1.0;
        su = 0;
    }
}
