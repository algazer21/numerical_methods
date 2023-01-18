/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos
                     CIMAT 6/10/2022
      Aplica Mínimos cuadrados para función de funciones 
***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
   return(fabs(x) -x/2 -x*x);               //fabs(x) -x/2 -x*x //3*x*x*x - 10*x*x + 4*x -5
}

void createvalues(size_t n, double *x, double *y,int l1, int l2);
void sol_tri_inf(double **a,double *x, double *b, int n);
void sol_tri_sup(double **a,double *x, double *b, int n);
void crout(size_t rows,double **Ax);
void minSqurtFun(size_t n, double *x, double *y, int m);
void savepoints(char *file,int m, double *x, double *y);
void error(int n, double *x, double *y);
double multFunc(double x, int z);
void freel(int n,double **listm);

/*********************************************************/
int main(void){ 
   int n = 5, p = 21;                                     //n = grado del polin, p = puntos
   double lim1 = -6.0, lim2 = 6.0;                         //limites de f(x)

   double *x = (double *)malloc(p*sizeof(double));
   double *y = (double *)malloc(p*sizeof(double));

   createvalues(p,x,y,lim1,lim2);
   minSqurtFun(n,x,y,p);

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
   savepoints("realval.dat",n,x,y);
}

void minSqurtFun(size_t n, double *x, double *y, int m){
   int i,j,k;
   double t;
   double **xx = (double **)malloc((n+1)*sizeof(double*));
    for(i = 0; i < (n+1); i++){
        xx[i] = (double *)malloc((n+1)*sizeof(double));}

   double *yp = (double *)malloc(m*sizeof(double));
   double *an = (double *)malloc((n+1)*sizeof(double));
   double *b = (double *)malloc((n+1)*sizeof(double));
   double *po = (double *)malloc((n+1)*sizeof(double));
   double *aux2 = (double *)malloc((n+1)*sizeof(double));

   for (i = 0; i < (n+1); i++){
      t = 0.0;
      for (j = 0; j < m; j++){           //Llenar matriz de polinomios
         t += multFunc(x[j],i)*y[j];
      }
      b[i] = t;
      for (j = i; j<(n+1); j++){
         t = 0.0;
         for (k = 0; k < m; k++){
            t+= multFunc(x[k],i)*multFunc(x[k],j);
         }
         xx[i][j] = t;
         if(j!=i){
            xx[j][i] = t;
         }
      }
   }

   crout(n+1,xx);
   sol_tri_inf(xx,aux2,b,n+1);          //Resolver matriz
   sol_tri_sup(xx,an,aux2,n+1);

   for (i = 0; i < m; ++i){
      yp[i] = 0;
      for (j = 0; j < (n+1); ++j){
         yp[i] += an[j]*multFunc(x[i],j);      //Polinomio
      }
   }
   savepoints("interpoval.dat",m,x,yp);
   printf("%d Puntos generados guardados exitosamente en 'interpoval.dat'\n",m);
   error(m,x,yp); 
   free(yp); free(an); freel(n+1,xx); free(aux2); free(b); free(po);
}

void savepoints(char *file,int m, double *x, double *y){
   FILE *pf;
   pf = fopen(file, "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i < m; ++i){
      fprintf(pf,"%lf %lf\n",x[i],y[i]);
   }
   fclose(pf);
}

void error(int n, double *x, double *y){
   double sum = 0;
   for (int i = 0; i < n; ++i){
      sum+= (y[i]-fun(x[i]))*(y[i]-fun(x[i]));
   }
   printf("Error calculado con la norma euclidiana = %.11lf\n",sqrt(sum));
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

void sol_tri_inf(double **a,double *x, double *b, int n){
    int i,j;
    double sum = 0.0;
    x[0] = b[0]/(a[0][0]);
    for(i = 1;i<n;i++){
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
     for(i=0;i<n;i++){
        a[i][i] = 1.0;
    }
    for(i=(n-1);i>=0;i--){
        for(j=i+1;j<n;j++){
            sum+=a[i][j]*x[j];
        }
        x[i] = (b[i]-sum)/a[i][i];
        sum = 0;
    }
}

double multFunc(double x, int z){
   return((x+(double)z)*(x+(double)z)*(x+(double)z));
}