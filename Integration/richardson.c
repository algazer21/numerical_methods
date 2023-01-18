/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 24/10/2022
       Aproxima exp(x) con extrapolacion de Richardson
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double richard(double Q2, double Q, int p);
double func(double x,double h){
    return(pow(1+x*h,(1.0/h)));
}

int main(void){
    clock_t t = clock();
    double tim,Qa,Qaa,realval = exp(1.0), Q, x = 1.0, h = 0.1;
    
    Q = richard(func(1.0,h/2.0),func(1.0,h),1);
    Qaa = richard(func(1.0,h/4.0),func(1.0,h/2.0),1);

    printf("Aproximacion 1 %lf con error %.3E\n",Q,fabs(Q-realval));
    printf("Aproximacion 2 %lf con error %.3E\n",Qaa,fabs(Qaa-realval));

    Qa = Qaa;
    Qaa = richard(Qaa,Q,2);
    Q = Qa;

    printf("Aproximacion 3 %lf con error %.3E\n",Qaa,fabs(Qaa-realval));
    
    t = clock()-t;
    tim = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\n",tim);
    return 0;
}

/**************************************************************************/

double richard(double Q2, double Q, int p){
    double qq;
    qq = (pow(2.0,p)*Q2-Q)/(pow(2.0,p)-1.0);
    return qq;
}
