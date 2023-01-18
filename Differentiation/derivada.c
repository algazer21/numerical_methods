/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 17/11/2022
    Deriva una funcion interpolada con polinomios de Lagrange
                    hasta tercer orden. 
***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Uso las derivadas analiticas solo para calcular el error

double fun(double x){ return(2.0*x*x*x*x - 3.0*x*x*x);}
double dfun(double x, int cc){
   if(cc == 1){return(8.0*x*x*x - 9.0*x*x);}
   if(cc == 2){return(24.0*x*x - 18.0*x);}
   if(cc == 3){return(48.0*x - 18.0);}
   else{printf("ERROR\n");return(0);}
}
/*
double fun(double x){ return(0.2*cos(2.0*x));}
double dfun(double x, int cc){
   if(cc == 1){return(-0.4*sin(2.0*x));}
   if(cc == 2){return(-0.8*cos(2.0*x));}
   if(cc == 3){return(1.6*sin(2.0*x));}
   else{printf("ERROR\n");return(0);}
}
*/

double dfunc(size_t,double*,double *,double,int,double);
void createpoints(size_t,double *,double *,int,int);
double lagrange(size_t, double *, double *, double);
void savepoints(int, double *, double *);
void error(int, double *, double *,int);

/*********************************************************/
int main(void){ 
   int cc,m = 22,n = 8;                                     
   double d,h,lim1 = -1.0, lim2 = 2.0;
   double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));
   double *xp = (double *)malloc(m*sizeof(double));
   double *yp = (double *)malloc(m*sizeof(double));
   createpoints(n,x,y,lim1,lim2);        //Genero n puntos en func de lim1 a lim2
                           
   cc = 3;                               
   h = 0.3;

   d = fabs((lim1+2.0*h)-(lim2-2.0*h))/(m-1);
   xp[0] = lim1+2.0*h;
   for (int i = 1; i < m; ++i){xp[i] = xp[i-1]+d;}                   

   for (int i = 0; i < m; ++i){
      yp[i] = dfunc(n,x,y,xp[i],cc,h);
   }

   savepoints(m,xp,yp);
   error(m,xp,yp,cc);

   printf("Datos guardados con exito\n");
   free(x); free(y); free(xp); free(yp);
   return 0;
}
/*********************************************************/

void createpoints(size_t n, double *x, double *y, int l1, int l2){
   double d = fabs(l1-l2)/(n-1);
   x[0] = l1;
   y[0] = fun(x[0]);
   for (int i = 1; i < n; ++i){
      x[i] = x[i-1]+d;
      y[i] = fun(x[i]);
   }
   FILE *pf;
   pf = fopen("f.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i < n; ++i){
      fprintf(pf,"%lf %lf\n",x[i],y[i]);
   }
   fclose(pf);
}

double dfunc(size_t n, double *x, double *y, double xp, int ord, double h){
   double a,b,c,d,m;
   a = lagrange(n,x,y,xp+h);
   b = lagrange(n,x,y,xp-h);
   c = lagrange(n,x,y,xp+2.0*h);
   d = lagrange(n,x,y,xp-2.0*h);
   m = lagrange(n,x,y,xp);

   if(ord == 1){return(((2.0/3.0)*a - (2.0/3.0)*b - (1.0/12.0)*c + (1.0/12.0)*d)/h);}

   if(ord == 2){ return(((4.0/3.0)*a + (4.0/3.0)*b - (1.0/12.0)*c - (1.0/12.0)*d - (5.0/2.0)*m)/(h*h));}

   if(ord == 3){return((-a + b + (1.0/2.0)*c - (1.0/2.0)*d)/(h*h*h));}

   else{printf("ERROR\n");return(0);}
}

double lagrange(size_t n, double *x, double *y, double xp){
   int i,k;
   double Pn = 0.0,L;
   if(xp > x[n-1] || xp < x[0]){
      printf("El punto no esta en el intervalo\n"); return(0);}

   for(k=0;k<n;k++){
      L=1.0;
      for(i = 0; i < n; ++i){
         if(i!=k){
            L = L*(xp-x[i])/(x[k]-x[i]);
         }
      }  
      Pn+=y[k]*L;
   }
   return(Pn);   
}

void savepoints(int m, double *x, double *y){
   FILE *pf;
   pf = fopen("fp3.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i < m; ++i){
      fprintf(pf,"%lf %lf\n",x[i],y[i]);
   }
   fclose(pf);
}

void error(int n, double *x, double *y, int cc){
   double sum = 0;
   for (int i = 0; i < n; ++i){
      sum+= (y[i]-dfun(x[i],cc))*(y[i]-dfun(x[i],cc));
   }
   printf("Error calculado con la norma euclidiana => %.1E\n",sqrt(sum));
}

