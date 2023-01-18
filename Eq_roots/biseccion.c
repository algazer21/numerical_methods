/*********************************************
Método de bisección
Alan García Zermeño
3/8/2022
*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415

void biseccion(double,double,double,double);
double func(double);

int main()
{
    double x_i = 12.0, x_d = 15.0, t = 1e-5, e = 12.49;
    biseccion(x_i,x_d,t,e);

    return 0;
}

double func(double X){
    double f = log(X*X+1) - exp(0.4*X)*cos(PI*X);//X*X*X-21*X*X+120*X -100;
    return(f);
}

void biseccion(double xi,double xd,double tol,double er){
    double x,xm,fxm,fxi;
    int k=0,b=1;

    while(fabs(xd-xi)>tol){
        k+=1;
        if(func(xi)*func(xd) > 0){
            printf("Error, no se encontró una raíz\n");
            b = 0;
            break;
        }
        else{
            xm = (xd+xi)/2;
            fxm = func(xm);
            fxi = func(xi);
            if(fxm*fxi < 0){
                xd = xm;
            }
            else{
                xi = xm;
            }
        }
    }
    er = fabs(er-xm)/er;

    if(b!=0){
        printf("Se encontró una raíz en: %lf, con %d iteraciones y un error relativo %lf\n",xm,k,er);
    }

}
