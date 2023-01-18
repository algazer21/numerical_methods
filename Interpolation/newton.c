/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 6/10/2022
    Interpola una funcion f(x) usado Polinomios de Newton 
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
   return(fabs(x) -x/2 -x*x);               //fabs(x) -x/2 -x*x //3*x*x*x - 10*x*x + 4*x -5
}

void createvalues(size_t n, double *x, double *y,int l1, int l2);
void sol_tri_inf(double **a,double *x, double *b, int n);
void newtonPoly(size_t n, double *x, double *y, int m);
void savepoints(int m, double *x, double *y);
void error(int n, double *x, double *y);
void freel(int n,double **listm);

/*********************************************************/
int main(void){ 
   int n = 15, p = 100;                                     //n = num ptos de f(x), p = num ptos a interpolar
   double lim1 = -1.0, lim2 = 1.0;                         //limites de f(x)

   double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));

   createvalues(n,x,y,lim1,lim2);
   newtonPoly(n,x,y,p);

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

void newtonPoly(size_t n, double *x, double *y, int m){
   int i,j,k;
   double Nn,L,d = fabs(x[1]-x[n-2])/(m-1);

   double **xx = (double **)malloc(n*sizeof(double*));
    for(i = 0; i < n; i++){
        xx[i] = (double *)malloc(n*sizeof(double));}

   double *xp = (double *)malloc(m*sizeof(double));
   double *yp = (double *)malloc(m*sizeof(double));
   double *an = (double *)malloc(n*sizeof(double));

   xp[0] = x[1];
   for (j = 1; j < m; ++j){               //Inicializar puntos a evaluar en
      xp[j] = xp[j-1]+ d;                 //el rango definido
      yp[j] = 0.0;
   }

   for (i = 0; i < n; ++i){
      xx[i][0] = 1.0;
      for (j = 1; j <= i; ++j){           //Llenar matriz de polinomios
         Nn = 1.0;
         for (k = 0; k < j ; ++k){
            Nn *= (x[i]-x[k]);
         }
         xx[i][j] = Nn;
      }
   }

   sol_tri_inf(xx,an,y,n);

   for (i = 0; i < m; ++i){
      for (j = 0; j < n; ++j){
         Nn = 1.0;
         for (k = 0; k < j ; ++k){       //Evaluar en el polinomio
            Nn *= (xp[i]-x[k]);
         }
         yp[i] += an[j]*Nn;
      }
   }

   savepoints(m,xp,yp);
   error(m,xp,yp); 
   free(xp); free(yp); free(an); freel(n,xx);
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

void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

