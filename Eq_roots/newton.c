/*********************************************
Método de Newton-Raphson
Alan García Zermeño
Para el curso de métodos numéricos.
CIMAT 24/8/2022
*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415

void newton(double,double,double);
double dfunc(double);
double func(double);


int main()
{
    double x0 = -1.0,t = 1e-5;
    newton(x0,1.0,t);
    return 0;
}


double func(double X){
    double f = (0.99)*X*X*X-21*X*X+120*X -100;//log(X*X+1) - exp(0.4*X)*cos(PI*X);//2-log(X)/X;//
    return(f);
}

double dfunc(double X){
    double f = (0.99)*3*X*X-42*X+120;//PI*exp(0.4*X)*sin(PI*X)- (2*exp(0.4*X)*cos(PI*X))/(5.0) +(2*X)/(X*X+1);//(log(X)-1)/(X*X);//
    return(f);
}

void newton(double p0,double err,double tol){   //(p0, p real(Para medir error relativo), tolerancia)
    int i,N0 = 5000;
    double fx,dfx,p;

    for(i=1;i<N0;i++){
        dfx = dfunc(p0);
        fx = func(p0);
        p = p0 - (fx/dfx);

        if(fabs(p-p0)<tol){
            break;
        }
        p0=p;
    }
    err = fabs(err-p)/err;

    if(i==N0){
        printf("No se pudo encontrar una raiz en %d pasos\n",i);
    }
    else{
        printf("Se encontró una raíz en: %lf, con %d iteraciones y un error relativo %lf\n",p,i,err);
    }
}