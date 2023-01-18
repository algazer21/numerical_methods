/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 11/10/2022
    Interpola una funcion f(x) usando splines lineales 
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
   //return(0.3*x*x*x*x - 2*x*x + 0.4*x -0.5);
   return(x*(1 + (sin(x/2.0))/3.0));               //fabs(x) -x/2 -x*x //3*x*x*x - 10*x*x + 4*x -5
}

void createvalues(size_t n, double
 *x, double *y,int l1, int l2);
void splines(size_t n, double *x, double *y, int m);
void savepoints(int m, double *x, double *y);
void error(int n, double *x, double *y);

/*********************************************************/
int main(void){ 
   int n = 23, p = 4;                                     //n = num ptos de f(x), p = num ptos a interpolar
   double lim1 = -0.0, lim2 = 50.0;                         //limites de f(x)

   double *x = (double *)malloc((n+1)*sizeof(double));
   double *y = (double *)malloc((n+1)*sizeof(double));

   createvalues(n,x,y,lim1,lim2);
   splines(n,x,y,p);

   free(x); free(y);
   return 0;
}
/*********************************************************/

void createvalues(size_t n, double *x, double *y, int l1, int l2){
   double d = fabs(l1-l2)/(n);
   x[0] = l1;
   y[0] = fun(x[0]);
   for (int i = 1; i <= n; ++i){
      x[i] = x[i-1]+d;
      y[i] = fun(x[i]);
   }
   FILE *pf;
   pf = fopen("realval.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i <= n; ++i){
      fprintf(pf,"%lf %lf\n",x[i],y[i]);
   }
   fclose(pf);
}

void splines(size_t n, double *x, double *y, int m){
   int i,j,m2;
   double t,t2,h = fabs(x[0]-x[n])/(n);
   n+=1;
   m2 = m*(n-1) + 1;
   double *xp = (double *)malloc(m2*sizeof(double));
   double *yp = (double *)malloc(m2*sizeof(double));
   
   for (i = 0; i < m2; ++i){
      xp[i] = x[0] + (i)*(h/(m));
   }
 
    for (i = 0; i < (n-1); ++i){
      t = (y[i+1]-y[i]);
      t2 = (x[i+1]-x[i]);
      for (j = 0; j < m; ++j){
         yp[i*m+j] = t*((xp[i*m+j]-x[i])/t2) +y[i];
      }      
   }
   yp[m2-1] = y[n-1];
   savepoints(m2,xp,yp);
   error(m2,xp,yp); 
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
   printf("%d Puntos generados guardados exitosamente en 'interpoval.dat'\n",m);
}

void error(int n, double *x, double *y){
   double sum = 0;
   for (int i = 0; i < n; ++i){
      sum+= (y[i]-fun(x[i]))*(y[i]-fun(x[i]));
   }
   printf("Error calculado con la norma euclidiana = %lf\n",sqrt(sum));
}
