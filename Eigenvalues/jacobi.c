/*********************************************
Alan García Zermeño
17/9/2022
Encuentra los n valores propios aplicando rotaciones a A
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI 3.141592653

void normalize(size_t rows, double *x);
void writevector(size_t rows,size_t cols,double *a,const char* filename);
void jacobi(size_t rows,double **a,double tol,double **Ap,double *vv,double **vec);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
void TeoVal(double *eva, double **eve,int n);
void freel(int n,double **listm);
int *maxim(double **a,int rows);

int main()
{
    int n,m,i,k;
    clock_t t = clock();
    double time;
    FILE    *fptr;
    fptr = fopen("At.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    fclose(fptr);
    double **A = (double **)malloc(n*sizeof(double*));
    double **Ap = (double **)malloc(n*sizeof(double*));
    double **vecp = (double **)malloc(n*sizeof(double*));
    double *vp = (double *)malloc(n*sizeof(double));
    for (i = 0;i<n;i++){
        Ap[i] = (double *)malloc(n*sizeof(double));
        vecp[i] = (double *)malloc(n*sizeof(double));   
        for (k = 0;k<m;k++){
            Ap[i][k] = 0.0;
            if (i!=k){
                vecp[i][k] = 0.0;}
            else{vecp[i][k] = 1.0;}
            }
        }

    readmatrix(n, m, A, "At.txt");
    jacobi(n,A,10e-4,Ap,vp,vecp);

    writematrix(n,m,vecp,"vecp.txt");
    writevector(n,1,vp,"vp.txt");
    TeoVal(vp,vecp,n);

    freel(n,A);
    freel(n,Ap);
    freel(n,vecp);
    free(vp);
 
    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("\nTiempo: %lf segundos\nSoluciones guardadas en vp.txt y vecp.txt\n",time);
    return 0;
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
    for(i=0;i<rows;i++){
            fprintf(fptr,"%.7lf\n",a[i]);
        }
    fclose(fptr);
}

/*******************************************************************/
void jacobi(size_t rows,double **a,double tol,double **Ap,double *vv, double **vec){
       int p,q,i,j,it = 0,itmax = 100000,*g,k;
       double thet,c,s;
       double **vecp = (double **)malloc(rows*sizeof(double*));
       for (i = 0;i<rows;i++){
        vecp[i] = (double *)malloc(rows*sizeof(double));   
        for (k = 0;k<rows;k++){
            vecp[i][k] = vec[i][k];
            }
        }

       while(it < itmax){
            //Encontrando el máximo
            g = maxim(a,rows);
            p = g[0];
            q = g[1];
            if(fabs(a[p][q])<tol){break;}

            // Inicializar A'
            for (i = 0; i < rows; i++){
                for (j = 0; j < rows; j++){
                    if ((j!=p) && (j!=q))
                    {
                        vecp[i][j] = vec[i][j];
                    }
                    if(i != p && i != q){
                        if(j != p && j != q){
                            Ap[i][j] = a[i][j];}
                    }
                }
            }

            //Calcular A'
            thet = 0.5*atan2(2*a[p][q],(a[q][q]-a[p][p]));
            c = cos(thet);
            s = sin(thet);
            for(i=0;i<rows;i++){
                if((i != p) && (i != q)){
                    Ap[p][i] = c*a[p][i]-s*a[q][i];
                    Ap[i][p] = Ap[p][i];
                }
                vecp[i][p] = c*vec[i][p]-s*vec[i][q];
            }
            for(i=0;i<rows;i++){
                if((i != p) && (i != q)){
                    Ap[q][i] = s*a[p][i]+c*a[q][i];
                    Ap[i][q] = Ap[q][i];}
                vecp[i][q] = s*vec[i][p]+c*vec[i][q];
            }
            Ap[p][p] = c*c*a[p][p] - 2*s*c*a[p][q] + s*s*a[q][q];
            Ap[q][q] = s*s*a[p][p] + 2*s*c*a[p][q] + c*c*a[q][q];
            Ap[p][q] = c*s*(a[p][p]-a[q][q])+(c*c - s*s)*a[p][q];
            Ap[q][p] = Ap[p][q];

            //A = A'
            for (i = 0; i < rows; i++){
                for (j = 0; j < rows; j++){
                    a[i][j] = Ap[i][j];
                    vec[i][j] = vecp[i][j];}
            }

            it +=1;
       }
       freel(rows,vecp);
       for (i = 0; i < rows; i++){vv[i] = a[i][i];}
        printf("Convergencia con %d iteraciones\n",it);
}

int *maxim(double **a,int rows){
    int i,j;
    static int k[2];
    double max = 0.0;
    for(i=0;i<rows;i++){
        for(j=0;j<rows;j++){
            if(i!=j){
                if(max<fabs(a[i][j])){
                    max = fabs(a[i][j]);
                    k[0] = i;k[1] = j;}
            }     
        }
    }
    return(k);
}

void TeoVal(double *eva, double **eve,int n){             //Valores teóricos
    int i,j;
    double h = 1.0/((double)(n+1)),sum=0,g;
    double *vecpr = (double *)malloc(n*sizeof(double));
    for(i=0;i<n;i++){
        g = fabs(eva[i] - (2.0/(h*h))*(cos(PI*((double)(n-i))*h)-1.0));
        //printf("%lf\n",(2.0/(h*h))*(cos(PI*((double)(n-i))*h)-1.0));
        printf("Error de eigenvalor %d: %.11lf\n",i+1,fabs(g));
    }
    printf("\n");
    for (int i = 1; i <=n; i++){
        for (j = 1; j<=n;j++){
            vecpr[n-j] = sin((double)i*PI*(double)j*h);
        }
        normalize(n,vecpr);
        for (j = 0; j<n;j++){ 
            g=vecpr[j]-eve[i-1][j];
            sum+=g*g;}
        //printf("\n");
        printf("Error de eigenvector %d: %.11lf\n",i, sqrt(sum));
        sum=0;
    }
    free(vecpr);
}

void normalize(size_t rows, double *x){
    double sum=0;
    for (int i=0;i<rows;i++){
        sum+= x[i]*x[i];}
    sum = sqrt(sum);
    for (int i=0;i<rows;i++){
        x[i] = x[i]/sum;}
}
