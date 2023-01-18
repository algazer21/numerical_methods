/**************************************************************
                    Alan García Zermeño
                        22/9/2022
Aplica método de Iteracion en subespacio para NE valores propios
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define threads 4

void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
void eigeniter(double *eval, double **A, double **V, int N, int M);
void updtCOL(int n,double **a,double *aa, int col);
void inicializarM(double **mtx, int cols, int rows);
double normaVec(double *x1, double *x2, int N);
void gsnorm(double **V, int N, int M);
void prod_matr_matr(int n, double **matrizA, double **matrizB, double **producto);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void writevector(size_t rows,double *a,const char* filename);
void IteracionSubespacio(double **A, double **V, double *eval, int N, int M, double tol);
void extractCOL(int n,double **a,double *aa, int col);
double prod_vec_vec(double *b,double *c, size_t rows);
void normalize(size_t rows, double *x);
void freel(int n, double **listm);
void error(double *x, int n);

int main(void){
	int n,m,i, NE = 10;         //Numero de eigenvalores y vectores a calcular
    clock_t t = clock();
    double time,er;
    FILE    *fptr;
    fptr = fopen("AA.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    if(NE > n){printf("El numero de valores propios a calcular tiene que ser menor a %d\n",n);exit(-1);}
    fclose(fptr);
    omp_set_num_threads(threads);
    double **A = (double **)malloc(n*sizeof(double*));
    double **V = (double **)malloc(n*sizeof(double*));
    double *eval = (double *)malloc(n*sizeof(double));
    for(i = 0; i < n; i++){
        V[i] = (double *)malloc(n*sizeof(double));}

    readmatrix(n, m, A, "AA.txt");
    IteracionSubespacio(A,V,eval,n,n,10e-3);
    writematrix(n,NE,V,"evecM.txt");
    writevector(NE,eval,"evalM.txt");
    //error(eval,NE);

    freel(n,A); freel(n,V);
    free(eval);

    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("\nTiempo: %lf segundos\nSoluciones guardadas en eval.txt y evec.txt\n",time);
    return 0;
}

/*************************Iteracion en Subespacio*********************************/
void IteracionSubespacio(double **A, double **V, double *eval, int N, int M, double tol){
    
    int it = 0, itmax = 30,i,j;
    double er=1.0,**tmp;
    double **Vt = malloc(N*sizeof(double*));
    double **Vt0 = malloc(N*sizeof(double*));
    double *evalt = (double *)malloc(M*sizeof(double));

    for(i = 0; i < N; i++){
        Vt0[i] = (double *)malloc(M*sizeof(double)); 
        Vt[i] = (double *)malloc(M*sizeof(double));}
 
    inicializarM(V,N,M);
    gsnorm(V,N,M);  
    eigeniter(evalt,A,V,N,M);

    while((er>tol)&&it<itmax){
        printf("Iteracion %d ...\n",it+1);
        prod_matr_matr(N, A, V, Vt0);
        gsnorm(Vt0,N,M);
        tmp = V;
        V = Vt0;
        Vt0 = tmp;
        tmp = NULL;
        eigeniter(eval,A,V,N,M);
        er = normaVec(eval,evalt,M);
        for (int i = 0;i<M;i++){evalt[i] = eval[i];}
        it+=1;
    }
    Vt0 = NULL;
    printf("Convergencia con %d iteraciones.\n",it);
    free(Vt); free(Vt0); free(evalt);
    return;

}
/*************************Iteracion en Subespacio*********************************/

void prod_matr_matr(int n, double **matrizA, double **matrizB, double **producto){
    double suma = 0;
    for (int a = 0; a < n; a++) {
        for (int i = 0; i < n; i++) {
            suma = 0;
            #pragma omp parallel for reduction(+:suma)
            for (int j = 0; j < n; j++) {
                suma += matrizA[i][j] * matrizB[j][a];
            }
            producto[i][a] = suma;
        }
    }
}

void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

void writevector(size_t rows,double *a,const char* filename){
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    printf("Eigenvalores mayores calculados:\n\n");
    for(int i=0;i<rows;i++){
            fprintf(fptr,"%.7lf\n",a[i]);
            printf("%.7lf\n",a[i]);
        }
    fclose(fptr);
}

void writematrix(size_t rows,size_t cols,double **a,const char* filename)
{
    int i,j;
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    for(j=0;j<cols;j++){
            fprintf(fptr,"%9d ",j+1);}
    fprintf(fptr,"\n");
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            fprintf(fptr,"%9.6lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

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

void inicializarM(double **mtx,int rows,int cols){
    for (int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            mtx[i][j] = (double)(rand()%15);
        }
    }
}

void gsnorm(double **V, int N, int M){
    double sigma=0.0;
    double *auxu = (double *)malloc(N*sizeof(double));
    double *auxv = (double *)malloc(N*sizeof(double));

    extractCOL(N,V,auxu,0);
    normalize(N,auxu);
    updtCOL(N,V,auxu,0);
    for(int j=1;j<M;j++){
        extractCOL(N,V,auxv,j);
        for(int k=0;k<j;k++){
            for(int i=0;i<N;i++){sigma += auxv[i]*V[i][k]; }
            for(int i=0;i<N;i++){auxv[i] -= V[i][k]*sigma;}
            sigma=0.0;}
    normalize(N,auxv);
    updtCOL(N,V,auxv,j);
    }
    free(auxu); free(auxv);
}

void extractCOL(int n,double **a,double *aa, int col){
    for(int i=0;i<n;i++){aa[i] = a[i][col];}
}
void updtCOL(int n,double **a,double *aa, int col){
    for(int i=0;i<n;i++){a[i][col] = aa[i];}
}

void normalize(size_t rows, double *x){
    double sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<rows;i++){
        sum+= x[i]*x[i];}
    sum = sqrt(sum);
    for (int i=0;i<rows;i++){
        x[i] = x[i]/sum;}
}

void eigeniter(double *eval, double **A, double **V, int N, int M){
    int i;
    double *vt = (double *)malloc(N*sizeof(double));
    double *t = (double *)malloc(N*sizeof(double));
    for (i = 0; i<M; i++){
        extractCOL(N,V,vt,i);
        prod_matr_vec(A,vt,t,N);
        eval[i] = prod_vec_vec(vt,t,N);
    }
    free(vt);free(t);
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

double normaVec(double *x1, double *x2, int N){
    double sum=0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<N;i++){
        sum += (x1[i]-x2[i])*(x1[i]-x2[i]);
    }
    return sum;
}

void error(double *x, int n){
    double g,sum=0.0;
    for (int i = n; i > 0; i--){
        g = 10.0*i;
        sum+=(g-x[n-i])*(g-x[n-i]);
    }
    printf("\nError de eigenvalores calculado con la norma euclidiana: %.7lf\n",sqrt(fabs(sum)));
}
