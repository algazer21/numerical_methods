/*********************************************
Alan García Zermeño
20/9/2022
Resuelve sistemas Ax=b usando Gradiente Conjugado.
Solo matrices simétricas positivas definidas.
*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
double error(size_t n,double **a,double *x,double *b);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
void readvector(size_t rows, double *a, const char* filename);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void writevector(size_t rows,double *a,const char* filename);
void gradienteC(size_t n, double **a, double *b, double *x, double tol);
double prod_vec_vec(double *b,double *c, size_t rows);
double norm(size_t rows, double *x);
void freel(int n,double **listm);

int main(void){
	int n,m,i;
    clock_t t = clock();
    double time,er;
    FILE    *fptr;
    fptr = fopen("A.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    fclose(fptr);
    double **A = (double **)malloc(n*sizeof(double*));
    double *b = (double *)malloc(n*sizeof(double));
    double *x = (double *)malloc(n*sizeof(double));

    readmatrix(n, m, A, "A.txt");
    readvector(n, b, "b.txt");
    gradienteC(n,A,b,x, 10e-10);
    writevector(n,x,"x.txt");
    er = error(n,A,x,b);
    printf("\nError calculado con norma euclidiana: %.20lf\n",er);

    freel(n,A);
    free(b); free(x);

    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}

/************************************************************************/

void readmatrix(size_t rows, size_t cols, double **a, const char* filename){
    int n,m;
    FILE *pf;
    pf = fopen (filename, "r");
    if (pf==NULL){
        printf("No se encontro archivo\n");
        exit(1);}
    fscanf(pf, "%d %d", &n,&m);
    for(size_t i = 0; i < rows; ++i){
        a[i] = (double *)malloc(cols*sizeof(double));
        for(size_t j = 0; j < cols; ++j){
            fscanf(pf, "%lf", a[i] + j);
        }    
    }
    fclose(pf);
}

void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

void readvector(size_t rows, double *a, const char* filename){
    int n,m;
    FILE *pf;
    pf = fopen (filename, "r");
    fscanf(pf, "%d %d", &n,&m);
    if (pf==NULL){
        printf("No se encontro archivo\n");
        exit(1);}
    for(size_t i = 0; i < rows; ++i){
            fscanf(pf, "%lf", a + i);
    }
    fclose(pf);
}

/*****************************Gradiente Conjugado****************************/
void gradienteC(size_t n, double **a, double *b, double *x, double tol){
    unsigned int i,j,k=0, itmax = 1000;
    double alpha, beta, num,den,nor;
    double *r = (double *)malloc(n*sizeof(double));
    double *p = (double *)malloc(n*sizeof(double));
    double *w = (double *)malloc(n*sizeof(double));

    for (i = 0; i < n; ++i){x[i] = rand()%10;}
    prod_matr_vec(a,x,w,n);
    for (i = 0; i < n; ++i){
        r[i] = -w[i]+b[i];                 //r = b - Ax
        p[i] = r[i];
    }

    nor = norm(n,r);
    while(nor>tol && k<itmax){
        prod_matr_vec(a,p,w,n);           //w = A*p
        num = prod_vec_vec(r,r,n);
        den = prod_vec_vec(p,w,n);
        alpha = num/den;                   //alpha = r*r/(p*A*p) = r*r/(p*w)
        for (i = 0; i < n; ++i){
            x[i]+= alpha*p[i];             //x = x + alpha*p
            r[i]-= alpha*w[i];             // r = r - alpha*Ap = r- alpha*w
        }
        den = num;
        num = prod_vec_vec(r,r,n);
        beta = num/den;                    //beta = r*r/rp*rp   --> rp es el r de la iteracion k-1
        for (i = 0; i < n; ++i){
            p[i] = r[i] + beta*p[i];}
        nor = norm(n,r);
        k+=1;
    }
    if(k!=itmax){printf("Convergencia con %d iteraciones\n",k);}
    else{printf("No se consiguió convergencia con %d iteraciones\n",itmax);}
    free(r); free(p); free(w);
}
/*****************************Gradiente Conjugado****************************/


void writevector(size_t rows,double *a,const char* filename)
{
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    printf("Soluciones al sistema:\n\n");
    for(int i=0;i<rows;i++){
            fprintf(fptr,"%.7lf\n",a[i]);
            printf("%.7lf\n",a[i]);
        }
    fclose(fptr);
}

void prod_matr_vec(double **a,double *v0,double *v1, size_t rows){
    for (int i = 0;i<rows;i++){
        v1[i] = prod_vec_vec(a[i],v0,rows);
    }
}

double prod_vec_vec(double *b,double *c, size_t rows){
    double a= 0;
    for (int i = 0;i<rows;i++){
        a += b[i]*c[i];
    }
    return a;
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
            fprintf(fptr,"%.2lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

double norm(size_t rows, double *x){
    double sum=0;
    for (int i=0;i<rows;i++){
        sum+= x[i]*x[i];}
    sum = sqrt(sum);
    return(sum);
}

double error(size_t n,double **a,double *x,double *b){   //||Ax - b||
    double sum=0.0;
    double *bx = (double *)malloc(n*sizeof(double));
    prod_matr_vec(a,x,bx,n);
    for (int i = 0;i <n;i++){
        sum+=(bx[i]-b[i])*(bx[i]-b[i]);
    }
    free(bx);
    return(sqrt(sum));
}
