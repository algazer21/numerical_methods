/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 22/9/2022
Aplica método de Iteracion en subespacio para NE valores propios
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define threads 4

void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
void sol_tri_inf(size_t rows, double *bx,double *xx ,double **Ax);
void sol_tri_sup(size_t rows, double *bx,double *xx ,double **Ax);
double normaVec(double *x1, double *x2, int N);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
void writevector(size_t rows,double *a,const char* filename);
void rayleigh(double **A, double *v, double eval, int n, double tol);
double prod_vec_vec(double *b,double *c, size_t rows);
void crout(size_t rows,double **Ax);
void normalize(size_t rows, double *x);
void freel(int n, double **listm);


int main(void){
	int n,m,i;
    clock_t t = clock();
    double time,lam;
    FILE    *fptr;
    fptr = fopen("AA.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    fclose(fptr);
    omp_set_num_threads(threads);
    double **A = (double **)malloc(n*sizeof(double*));
    double *v = (double *)malloc(n*sizeof(double));
    readmatrix(n, m, A, "A2.txt");

    lam = 50.0;
    rayleigh(A,v,lam,n,10e-8);
    writevector(n,v,"evec.txt");

    freel(n,A);
    free(v);

    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("\nTiempo: %lf segundos\nSoluciones guardadas en eval.txt y evec.txt\n",time);
    return 0;
}

/***************************Rayleigh*********************************/
void rayleigh(double **A, double *v, double eval, int n, double tol){
    int i,j,it = 0, itmax = 50;
    double er = 1.1, evalt;
    double **At = (double **)malloc(n*sizeof(double*));
    double *vt = (double *)malloc(n*sizeof(double));
    double *tmpa = (double *)malloc(n*sizeof(double));
    for (int i = 0; i < n; ++i){ 
        v[i] = (double)((rand()%3)+1);
        At[i] = (double *)malloc(n*sizeof(double));
    }

    evalt = eval;
    normalize(n,v);

    while((er>tol)&&it<itmax){                      //A* = A - sigma(I)
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                if(i==j){
                At[i][j] = A[i][j] - eval;}
                else{At[i][j] = A[i][j]-0.0;}
            }
        }

        crout(n,At);
        sol_tri_inf(n,v,tmpa,At);
        sol_tri_sup(n,tmpa,vt,At);                  //(A*)vn = v
        normalize(n,vt);

        for (int i=0;i<n;i++){v[i] = vt[i];}

        prod_matr_vec(A,v,tmpa,n);               
        eval = prod_vec_vec(v,tmpa,n);              //vn(A)vn = lamb
        er = fabs(evalt-eval);
        evalt=eval;
        it+=1;
    }

    er = fabs(eval-500.0);
    if(it == itmax){printf("No se ha alcanzado convergencia en %d iteraciones\n",it);}
    else{printf("Valor propio encontrado: %.7lf en %d iteraciones\n\nVector propio: [",eval,it);}
    for (i = 0;i<(n-1);i++){ printf("%lf, ",v[i]);} printf("%lf]\n",v[n-1]);
    printf("\nError calculado con norma euclidiana: %.5lf\n",er);
    freel(n,At); free(tmpa); free(vt);
}
/***************************Rayleigh*********************************/

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

void crout(size_t rows,double **Ax)
{
    int i,j,p;
    double a;
    for(j=0;j<rows;j++){
        for(i=j; i<rows; i++){
            a = 0;
            for(p=0;p<j;p++){
                a = a+Ax[i][p]*Ax[p][j]; }
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

void normalize(size_t rows, double *x){
    double sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<rows;i++){
        sum+= x[i]*x[i];}
    sum = sqrt(sum);
    #pragma omp parallel for
    for (int i=0;i<rows;i++){
        x[i] = x[i]/sum;}
}

void sol_tri_inf(size_t rows, double *bx, double *xx ,double **Ax)
{
    int i,j;
    double su=0;
    for(i=0;i<rows;i++){
        #pragma omp parallel for reduction(+:su)
        for(j=0;j<i;j++){
            su+=Ax[i][j]*xx[j];
        }
        xx[i] = (bx[i]-su)/Ax[i][i];
        su = 0;
    }
}

void sol_tri_sup(size_t rows, double *bx,double *xx ,double **Ax)
{
    int i,j;
    double su=0;
    #pragma omp parallel for
    for(i=0;i<rows;i++){
        Ax[i][i] = 1.0;
    }

    for(i=(rows-1);i>=0;i--){
        #pragma omp parallel for reduction(+:su)
        for(j=i+1;j<rows;j++){
            su+=Ax[i][j]*xx[j];
        }
        xx[i] = (bx[i]-su)/Ax[i][i];
        su = 0;
    }
}

void prod_matr_vec(double **a,double *v0,double *v1, size_t rows){
    #pragma omp parallel for
    for (int i = 0;i<rows;i++){
        v1[i] = prod_vec_vec(a[i],v0,rows);
    }
}

double prod_vec_vec(double *b,double *c, size_t rows){
    double a= 0;
    #pragma omp parallel for reduction(+:a)
    for (int i = 0;i<rows;i++){
        a += b[i]*c[i];
    }
    return a;
}

void writevector(size_t rows,double *a,const char* filename)
{
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