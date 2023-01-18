/*********************************************
Alan García Zermeño
8/9/2022
Encuentra los n valores propios más pequeños de A
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 10                                              //Numero de valores propios a calcular
#define PI 3.141592653

void normalize(size_t rows, double *x);
void crout(size_t rows, size_t cols,double **Ax);
void sol_tri_sup(size_t rows, size_t cols, double *bx, double *xx ,double **Ax);
void sol_tri_inf(size_t rows, size_t cols, double *bx, double *xx ,double **Ax);
void writevector(size_t rows,size_t cols,double *a,const char* filename);
void EigenPower_inv_n(size_t rows, size_t cols, double **a, double *y, double *v1, double tol,double **valp,int ii,double *vv);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
double prod_vec_vec(double *b,double *c, size_t rows);
void TeoVal(double *eva, double **eve,int n);
void freel(int n,double **listm);

int main()
{
    int n,m,i,k;
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
        x[i] = (double *)malloc(m*sizeof(double));   
        for (k = 0;k<m;k++){
            x[i][k] = 0.0;}
        }

    readmatrix(n, m, A, "AA.txt");
    crout(n,m,A);
    for(k=0;k<N;k++){
        EigenPower_inv_n(n,m,A,y,x2,10e-6,x,k,vp);
        printf("Eigenvalor %d = %lf... \n",n-k,vp[k]);}

    writematrix(N,n,x,"vecm.txt");
    writevector(n,1,vp,"valm.txt");
    //TeoVal(vp,x,n);
    freel(n,A);
    freel(n,x);
    free(y);free(vp);
 
    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("\n\nTiempo: %lf segundos\nSoluciones guardadas\n",time);
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
            fprintf(fptr,"%.3lf ",a[i][j]);
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
double prod_vec_vec(double *b,double *c, size_t rows){
    double a = 0.0;
    for (int i = 0;i<rows;i++){
        a += b[i]*c[i];
    }
    return a;
}
/****************************************************************************************/
void EigenPower_inv_n(size_t rows, size_t cols, double **a, double *y, double *v1 ,double tol,double **valp,int ii,double *vv){
    int i,j=0,nmax = 300,kk;
    double lamO=10e20,lam,sum,val=1.0,an=0;
    double *t = (double *)malloc(rows*sizeof(double));

    for(i=0;i<rows;i++){y[i] = rand()%3+1;}
    normalize(rows,y);

    while((val>tol) && (j<nmax)){
        //Para quitar contribución de los eigenvectores previos. 
        if(ii>0){
        for(i=0;i<ii;i++){
            an = prod_vec_vec(valp[i],y,rows);
            for(kk=0;kk<rows;kk++){
                y[kk] -= (valp[i][kk]*an);}}}

        j+=1;
        sol_tri_inf(rows,rows,y,t,a);
        sol_tri_sup(rows,rows,t,v1,a);
        sum = 0;
        for (i=0;i<rows;i++){
            sum+= y[i]*v1[i];}

        lam = 1.0/sum;
        val = fabs(lam-lamO);
        normalize(rows,v1);

        for (i = 0;i<rows;i++){
            y[i] = v1[i];}
        lamO = lam;
    }
    vv[ii] = lam;
    //printf("\n\nValor propio más pequeño: %.5lf calculado en %d iteraciones \n\n",lam,j);
    //printf("Vector propio más pequeño normalizado:\n");
    for (i = 0;i<rows;i++){
        //printf("%.6lf  ",v1[i]);
        valp[ii][i] = v1[i];}
    free(t);
}

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
            //if(i==j){su+=xx[j];}
            su+=Ax[i][j]*xx[j];
        }
        xx[i] = (bx[i]-su)/1.0;
        su = 0;
    }
}

void TeoVal(double *eva, double **eve,int n){          //Valores reales
    int i,j;
    double h = 1.0/(1000.0),sum=0,g;
    printf("%d eigenvalores calculados: \n",N);
    for(i=0;i<N;i++){
        g = fabs(eva[i] - (2.0/(h*h))*(cos(PI*((double)(i+1))*h)-1.0));
        sum+=g*g;
    }
    printf("Error de eigenvalores: %lf\n",sqrt(sum));
    sum=0;
    for (int i = 0; i <N; i++){
        for (j = 0; j<n;j++)
        {
            g = (eve[i][j]-sin((double)i*PI*(double)j*h));
            sum+=g*g;
        }
        printf("Error de eigenvector %d: %lf\n",i+1, sqrt(sum));
        sum=0;
    }
}
