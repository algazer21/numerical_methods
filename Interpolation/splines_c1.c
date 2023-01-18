/**************************************************************
                    Alan García Zermeño
                        11/10/2022
    Interpola una funcion f(x) usando splines cuadraticos 
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
   //return(fabs(x) -x/2 -x*x);
   return(x*(1 + (sin(x/2.0))/3.0));
   //return(0.3*x*x*x*x - 2*x*x + 0.4*x -0.5);               //fabs(x) -x/2 -x*x //3*x*x*x - 10*x*x + 4*x -5
}
double prod_vec_vec(double *b,double *c, size_t rows);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);

void createvalues(size_t n, double *x, double *y,int l1, int l2);
void gauss(size_t rows,double **Ax,double *xx,double *bx);
void splines3(size_t n, double *x, double *y, int m);
double errort(size_t n,double **a,double *x,double *b);
void savepoints(int m, double *x, double *y);
void error(int n, double *x, double *y);
void freel(int n,double **listm);

/*********************************************************/
int main(void){ 
   int n = 23, p = 100;                                     //n = num ptos de f(x), p = num ptos a interpolar
   double lim1 = 0.0, lim2 = 50.0;                         //limites de f(x)

   double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));

   createvalues(n,x,y,lim1,lim2);
   splines3(n,x,y,p);

   free(x); free(y);
   return 0;
}
/*********************************************************/

void createvalues(size_t n, double *x, double *y, int l1, int l2){
   double d = fabs(l1-l2)/(n-1);
   x[0] = l1;
   y[0] = fun(x[0]);
   for (int i = 1; i < n; ++i){
      x[i] = x[i-1]+d;
      y[i] = fun(x[i]);
   }
   FILE *pf;
   pf = fopen("realval.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i < n; ++i){
      fprintf(pf,"%lf %lf\n",x[i],y[i]);
   }
   fclose(pf);
}

void splines3(size_t n, double *x, double *y, int m){
   int i,j,k;
   double t,h;

   double *bx = (double *)malloc((n-1)*sizeof(double));
   double *bp = (double *)malloc(n*sizeof(double));
   double *xp = (double *)malloc(m*sizeof(double));
   double *yp = (double *)malloc(m*sizeof(double));
   double *bb = (double *)malloc((n-1)*sizeof(double));
   double **xx = (double **)malloc((n-1)*sizeof(double*));
    for(i = 0; i < (n-1); i++){
      xx[i] = (double *)malloc((n-1)*sizeof(double));}

   
   for (i = 0; i < (n-1); ++i){
      bx[i] = (y[i+1]-y[i]);
   }
   //bx[0] = 0.0;

   for (i = 0; i < (n-1); ++i){
      h = x[i+1]-x[i];
      for (j = 0; j < (n-1); ++j){
         xx[i][j] = ((i==j)||(i==(j+1)))?0.5*h:0.0;
      }

   }
   
   gauss(n-1,xx,bb,bx);
   printf("Error del sistema de ecuaciones: %.11lf\n",errort(n-1,xx,bb,bx));

   bp[0] = 0.0;
    for (i = 1; i < n; i++){
        bp[i] = bb[i-1];
    }

   double a,b,c,step = (x[n-1]-x[0])/((double)m);;
   j=0;
   a = 0.5*((bp[j+1]-bp[j])/(x[j+1]-x[j]));
   b = bp[j];
   c = y[j];
   for (i = 0; i < m; i++)
    {    
      xp[i] = x[0]+((double)i)*step;
        if((xp[i]>x[j+1])&&(j<(n-2))){
            j+=1;
            a = 0.5*((bp[j+1]-bp[j])/(x[j+1]-x[j]));
            b = bp[j];
            c = y[j];
        }
        t = (xp[i]-x[j]);
        yp[i] = a*t*t +b*t + c;
    }
   savepoints(m,xp,yp);
   error(m,xp,yp);
 
   free(bx);freel(n-1,xx);free(bb); free(xp); free(yp); free(bp);
}

void savepoints(int m, double *x, double *y){
   FILE *pf;
   pf = fopen("interpoval.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i < m; ++i){
      fprintf(pf,"%lf %lf\n",x[i],y[i]);
   }
   fclose(pf);
   printf("%d Puntos generados guardados exitosamente en 'datos.dat'\n",m);
}

void error(int n, double *x, double *y){
   double sum = 0;
   for (int i = 0; i < n; ++i){
      sum+= (y[i]-fun(x[i]))*(y[i]-fun(x[i]));
   }
   printf("Error calculado con la norma euclidiana = %lf\n",sqrt(sum));
}

void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

double errort(size_t n,double **a,double *x,double *b){   //||Ax - b||
   double sum=0.0;
   double *bx = (double *)malloc(n*sizeof(double));
   prod_matr_vec(a,x,bx,n);
   for (int i = 0;i <n;i++){
      sum+=(bx[i]-b[i])*(bx[i]-b[i]);
   }
   free(bx);
   return(sqrt(sum));
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

void gauss(size_t rows,double **Ax,double *xx,double *bx){
    int i,j,k;
    double c,sum=0.0;

    for(k=0;k<(rows-1);k++){
        for(i=k; i<(rows-1); i++){
            c = Ax[i+1][k]/Ax[k][k];
            for(j=k;j<rows;j++){
                Ax[i+1][j] = Ax[i+1][j]-Ax[k][j]*c;   
            }
            bx[i+1] = bx[i+1]-bx[k]*c;
        }
    }
    for(i=(rows-1);i>=0;i--){
        for(j=i+1;j<rows;j++){
            sum+=Ax[i][j]*xx[j];
        }
        xx[i] = (bx[i]-sum)/Ax[i][i];
        sum = 0;
    }
}
