/**************************************************************
                    Alan García Zermeño
                         6/10/2022
    Interpola una funcion f(x) usado Polinomios de Lagrange 
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
   return(fabs(x) -x/2 -x*x);        //1/(1+25*x*x) //fabs(x) -x/2 -x*x)
}

void createvalues(size_t n, double *x, double *y,int l1, int l2);
void lagrange(size_t n, double *x, double *y, int m);
void savepoints(int m, double *x, double *y);
void error(int n, double *x, double *y);

/*********************************************************/
int main(void){ 
   int n = 10, p = 100;                                     //n = num ptos de f(x), p = num ptos a interpolar
   double lim1 = -1.0, lim2 = 1.0;                         //limites de f(x)

   double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));

   createvalues(n,x,y,lim1,lim2);
   lagrange(n,x,y,p);

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

void lagrange(size_t n, double *x, double *y, int m){
   int i,j,k;
   double Pn,L,d = fabs(x[1]-x[n-2])/(m-1);
   double *xp = (double *)malloc(m*sizeof(double));
   double *yp = (double *)malloc(m*sizeof(double));

   xp[0] = x[1];
   for (j = 1; j < m; ++j){               //Inicializar puntos a evaluar en el rango definido
      xp[j] = xp[j-1]+ d; 
   }

   for (j = 0; j < m; ++j){
      Pn = 0;
      for(k=0;k<n;k++){
         L=1.0;
         for (i = 0; i < n; ++i){
            if(i!=k){
               L = L*(xp[j]-x[i])/(x[k]-x[i]);
            }
         }  
         Pn+=y[k]*L;
      }
      yp[j] = Pn;
   }
   savepoints(m,xp,yp);
   error(m,xp,yp);   
   free(xp); free(yp);
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
   printf("Error calculado con la norma euclidiana = %.11lf\n",sqrt(sum));
}

