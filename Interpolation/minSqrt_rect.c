/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos
                     CIMAT 6/10/2022
         Aplica Mínimos cuadrados a rectas ax + b
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fun(double x){
   return((7.3)*x-(5.2));  //Solo rectas
}

void createvalues(size_t n, double *x, double *y,int l1, int l2);
void minSqurtrect(double *x, double *y, int m);
void savepoints(char *file,int m, double *x, double *y);
void error(int n, double *x, double *y);
void expo(int n, double *x,double xx);
void freel(int n,double **listm);

/*********************************************************/
int main(void){ 
   int p = 21;                                            //p = puntos
   double lim1 = -2.0, lim2 = 3.0;                         //limites de f(x)

   double *x = (double *)malloc(p*sizeof(double));
   double *y = (double *)malloc(p*sizeof(double));

   createvalues(p,x,y,lim1,lim2);
   minSqurtrect(x,y,p);

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

void minSqurtrect(double *x, double *y, int m){
   int i,j;
   double t=0,t2 = 0,xx[2][2],an[2],b[2],au[2];
   double *yp = (double *)malloc(m*sizeof(double));

   b[0] = 0.0; b[1] = 0.0;
   for (i = 0; i < m; i++){
      t += x[i];
      t2 += x[i]*x[i];
      b[0] += x[i]*y[i];
      b[1] += y[i];
   }
      
   xx[0][0] = t2;
   xx[1][0] = t;
   xx[0][1] = t;
   xx[1][1] = m;

   //Resolver
   an[0] = (b[0]-(xx[0][1]/xx[1][1])*b[1])/(xx[0][0]-((xx[0][1]/xx[1][1])*xx[1][0]));
   an[1] = (b[1]-xx[1][0]*an[0])/xx[1][1];
   printf("a = %.3lf, b = %.3lf\n",an[0],an[1]);

   //Evaluar
   for (i = 0; i < m; ++i){
      yp[i] = an[0]*x[i] + an[1];
   }

   savepoints("interpoval.dat",m,x,yp);
   printf("%d Puntos generados guardados exitosamente en 'interpoval.dat'\n",m);
   error(m,x,yp);
   free(yp);
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
