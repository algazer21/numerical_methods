/**************************************************************
                    Alan García Zermeño
                         5/11/2022
 Integra funciones en el intervalo dado con Cuadratura de Gauss
***************************************************************/
#include <stdio.h>
#include <math.h>
#define RV 317.3442466738264

void Gaussqua(double,double);
double gauss(double,double,double);

double f(double x){
    return(x*x*x*x*x*x - x*x*sin(2.0*x));           //317.3442466738264  a = 1.0, b = 3.0;
    //return(x*exp(2.0*x));                         //7.565747863545      a = -0.6, b = 1.4;
    //return(2*x*x*x + x*x/10.0 - 6.0*x + 12.0);   //21.1546666666666          a = -0.6, b = 1.4;
}

void main(){
    double a = 1.0, b = 3.0;
    Gaussqua(a,b);
}

/*******************************************/
double gauss(double a,double b,double z){
    double x=z*(b-a)/2.0+(b+a)/2.0;
    return (f(x));
}

void Gaussqua(double a,double b){
    double in;
    printf("La integral de la funcion entre %.1lf y %.1lf aproxima a:\n",a,b);
    in=(b-a)/2.0*(1.0*gauss(a,b,-1.0/sqrt(3.0))+1.0*gauss(a,b,1.0/sqrt(3.0)));
    printf("2 puntos: in=%lf   error = %.2E\n",in,fabs(in-RV));

    in=(b-a)/2.0*(5.0/9.0*gauss(a,b,-sqrt(3.0/5.0))+8.0/9.0*gauss(a,b,0)+5.0/9.0*gauss(a,b,sqrt(3.0/5.0)));
    printf("3 puntos: in=%lf   error = %.2E\n",in,fabs(in-RV));

    in=(b-a)/2.0*(0.34785*gauss(a,b,-0.86114)+0.65215*gauss(a,b,-0.33998)+0.65215*gauss(a,b,0.33998)+0.34785*gauss(a,b,0.86114));
    printf("4 puntos: in=%lf   error = %.2E\n",in,fabs(in-RV));

    in=(b-a)/2.0*(0.236926*gauss(a,b,-0.90618) + 0.478628*gauss(a,b,-0.53847) + 0.56889*gauss(a,b,0.0) 
        + 0.236926*gauss(a,b,0.90618) + 0.478628*gauss(a,b,0.53847));
    printf("5 puntos: in=%lf   error = %.2E\n",in,fabs(in-RV));
}
