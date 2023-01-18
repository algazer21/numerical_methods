/**************************************************************
                    Alan García Zermeño
                        30/10/2022
      Integra funciones en el intervalo dado con Roomberg
***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592

double f(double x){
    //return(4.2*x*cos(8.0*x)*exp(-x*x));                  //0.1324588139019    a = -0.6, b= 1.4;
    //return((1.0/x)*log(x*x)*sin(4.0*x)*sin(3.0*x));   //-0.1732955599801922   lix = 3.0, ldx = 10.0;
    return(2*x*x*x + x*x/10.0 - 6.0*x + 12.0);        //21.154666...          a = -0.6, b = 1.4;
}
double romberg(double a, double b, double err);
void print_row(size_t i, double *R);

int main(){
   double ans,realval = 21.154666666666666, error = 1e-8,a = -0.6, b= 1.4;
   ans = romberg(a,b,error);
   printf("La integral de la funcion de %.1lf a %.1lf aproxima a \nint = %.11lf\n",a,b,ans);
   printf("Con un error estimado de aproximadamente %.1E\n",fabs(ans-realval));
   return 0;
}


/*****************************************************************************/
void print_row(size_t i, double *R){
   printf("R[%2zu] = ", i);
   for (int j = 0; j <= i; ++j){
      printf("%.3lf ", R[j]);
   }
   printf("\n");
}

double romberg(double a, double b, double err){  
   size_t ep,maxst = 30;
   double R1[maxst],R2[maxst],*rt,*Rp = &R1[0], *Rc = &R2[0],nk,h,c;
   h = fabs(b-a);

   Rp[0] = h*0.5*(f(a) + f(b));
   print_row(0, Rp);

   for (int i = 1; i < maxst; ++i) {
      h /= 2.0;
      c = 0;
      ep = 1 << (i-1);
      for (int j = 1; j <= ep; ++j) {
         c += f(a+(2*j-1)*h);
      }
      Rc[0] = h*c + 0.5*Rp[0];

      for (int j = 1; j <= i; ++j) {
         nk = pow(4, j);
         Rc[j] = (nk*Rc[j-1] - Rp[j-1])/(nk-1);
      }
      print_row(i, Rc);
      if (i > 1 && fabs(Rp[i-1]-Rc[i]) < err) {
         return Rc[i];
      }
      rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[maxst-1];
}
