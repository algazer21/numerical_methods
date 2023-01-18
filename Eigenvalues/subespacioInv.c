/**************************************************************
                    Alan García Zermeño
                        22/9/2022
Método Iteracion en subespacio Inverso para NE valores propios
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define threads 4

void sol_tri_sup(double **a,double *x, double *b, int n);
void error(double *x, int n);
void sol_tri_inf(double **a,double *x, double *b, int n);
void slvMultSis(double **a, double **V, double **X, int N, int M);
void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
void crout(size_t rows,double **Ax);
void eigeniter(double *eval, double **A, double **V, int N, int M);
void updtCOL(int n,double **a,double *aa, int col);
void inicializarM(double **a, int cols, int rows);
double normaVec(double *x1, double *x2, int N);
void gsnorm(double **V, int N, int M);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void writevector(size_t rows,double *a,const char* filename);
void IteracionSubespacioInv(double **A, double **V, double *eval, int N, int M, double tol);
void extractCOL(int n,double **a,double *aa, int col);
double prod_vec_vec(double *b,double *c, size_t rows);
void normalize(size_t rows, double *x);
void freel(int n, double **listm);


int main(void){
	int n,m,i, NE = 50;         //Numero de eigenvalores y vectores a calcular
    clock_t t = clock();
    double time,er;
    FILE    *fptr;
    fptr = fopen("A2.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    if(NE > n){printf("El numero de valores propios a calcular tiene que ser menor a %d\n",n);exit(-1);}
    fclose(fptr);
    omp_set_num_threads(threads);
    double **A = (double **)malloc(n*sizeof(double*));
    double **V = (double **)malloc(n*sizeof(double*));
    double *eval = (double *)malloc(n*sizeof(double));
    for(i = 0; i < n; i++){
        V[i] = (double *)malloc(n*sizeof(double));}

    readmatrix(n, m, A, "A2.txt");
    IteracionSubespacioInv(A,V,eval,n,n,10e-10);
    writematrix(n,NE,V,"evec.txt");
    writevector(NE,eval,"eval.txt");
    error(eval,NE);

    freel(n,A); freel(n,V);
    free(eval);

    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("\nTiempo: %lf segundos\nSoluciones guardadas en eval.txt y evec.txt\n",time);
    return 0;
}


/*************************Iteracion en Subespacio Inverso***************************/
void IteracionSubespacioInv(double **A, double **V, double *eval, int N, int M, double tol){
    
    int it = 0, itmax = 1000,i,j;
    double er=1.0,**tmp;
    double **Vt = malloc(N*sizeof(double*));
    double **Vt0 = malloc(N*sizeof(double*));
    double *evalt = (double *)malloc(M*sizeof(double));

    for(i = 0; i < N; i++){
        Vt0[i] = (double *)malloc(M*sizeof(double)); 
        Vt[i] = (double *)malloc(M*sizeof(double));}
    
    crout(N,A);
    inicializarM(V,N,M);
    gsnorm(V,N,M);  
    eigeniter(evalt,A,V,N,M);

    while((er>tol)&&it<itmax){
        slvMultSis(A,V,Vt0,N,M);
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
/*************************Iteracion en Subespacio Inverso***************************/

void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

void writevector(size_t rows,double *a,const char* filename)
{
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    printf("Eigenvalores menores calculados:\n\n");
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

void inicializarM(double **a,int rows,int cols){
    for (int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            a[i][j] = (double)(rand()%15);
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

void slvMultSis(double **a, double **V, double **X, int N, int M){
    int i,j;
    double *xt = (double *)malloc(N*sizeof(double));
    double *vt = (double *)malloc(N*sizeof(double));
    double *z = (double *)malloc(N*sizeof(double));
    double *d = (double *)malloc(N*sizeof(double));

    for (i=0;i<M;i++){
        extractCOL(N,V,vt,i);
        sol_tri_inf(a,z,vt,N);
        for (int j=0;j<N;j++){
            d[j] = a[j][j];
                a[j][j]=1.0;
        }
        sol_tri_sup(a,xt,z,N);
        for (j=0;j<N;j++){
             a[j][j] = d[j];
        }
        updtCOL(N,X,xt,i);
    }
    free(xt),free(vt);free(z);free(d);
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

void sol_tri_inf(double **a,double *x, double *b, int n){
    int i,j;
    double sum = 0.0;
    x[0] = b[0]/(a[0][0]);
    for(i = 1;i<n;i++){
        #pragma omp parallel for reduction(+:sum)
        for (j=0;j<i;j++){
            sum += a[i][j]*x[j];
        }
        x[i] = (b[i]-sum)/a[i][i];
        sum = 0.0;
    }
}

void sol_tri_sup(double **a,double *x, double *b, int n){
    int i,j;
    double sum = 0.0;
    x[n-1] = b[n-1]/(a[n-1][n-1]);  
    for(i=(n-2);i>=0;i--){
        #pragma omp parallel for reduction(+:sum)
        for (j=(i+1);j<n;j++){
            sum += a[i][j]*x[j];
        }
        x[i] = (b[i]-sum)/a[i][i];
        sum = 0.0;
    }
}

void error(double *x, int n){
    double g,sum=0.0;
    for (int i = 1; i <= n; i++){
        g = 10.0*i;
        sum+=(g-x[i-1])*(g-x[i-1]);
    }
    printf("\nError de eigenvalores calculado con la norma euclidiana: %lf\n",sqrt(fabs(sum)));
}
