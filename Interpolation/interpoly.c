/**************************************************************
                    Alan García Zermeño
                        6/10/2022
            Interpola funciones usando polinomios 
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
void interPoly(size_t n, double *x, double *y, int m);
void savepoints(int m, double *x, double *y);
void error(int n, double *x, double *y);
void expo(int n, double *x,double xx);
void freel(int n,double **listm);

/*********************************************************/
int main(void){ 
   int n = 15, p = 40;                                     //n = grado del polin, p = puntos
   double lim1 = -1.0, lim2 = 1.0;                         //limites de f(x)

   double *x = (double *)malloc(n*sizeof(double));
   double *y = (double *)malloc(n*sizeof(double));

   createvalues(n,x,y,lim1,lim2);
   interPoly(n,x,y,p);

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

void interPoly(size_t n, double *x, double *y, int m){
   int i,j,k;
   double t,d = fabs(x[1]-x[n-2])/(m-1);

   double **xx = (double **)malloc(n*sizeof(double*));
    for(i = 0; i < n; i++){
        xx[i] = (double *)malloc(n*sizeof(double));}

   double *yp = (double *)malloc(m*sizeof(double));
   double *xp = (double *)malloc(m*sizeof(double));
   double *an = (double *)malloc((n)*sizeof(double));
   double *po = (double *)malloc((n)*sizeof(double));
   double *aux2 = (double *)malloc((n)*sizeof(double));

   yp[0] = 0;
   xp[0] = x[1];
   for (j = 1; j < m; ++j){               //Inicializar puntos a evaluar en el rango definido
      xp[j] = xp[j-1]+ d;
      yp[j] = 0.0; 
   }

   for (i = 0; i < n; i++){
      t = 1.0;
      xx[i][0] = t;
      for (j = 1; j < n; j++){           //Llenar matriz de polinomios
         t*=x[i];
         xx[i][j] = t;
      }
   }

   crout(n,xx);
   sol_tri_inf(xx,aux2,y,n);          //Resolver matriz
   sol_tri_sup(xx,an,aux2,n);

   for (i = 0; i < m; ++i){
      yp[i] = 0;
      expo(n,po,xp[i]);
      for (j = 0; j < n; ++j){
         yp[i] += an[j]*po[j];      //Polinomio
      }
   }

   savepoints(m,xp,yp);
   printf("%d Puntos generados guardados exitosamente en 'interpoval.dat'\n",m);
   error(m,xp,yp); 
   free(yp); free(xp); free(aux2); free(an); freel(n,xx); free(po);
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

void expo(int n, double *x,double xx){
   x[0] = 1.0;
   for (int i = 1; i < n; ++i){
      x[i] = x[i-1]*xx;
   }
}
