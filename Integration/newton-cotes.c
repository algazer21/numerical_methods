/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 24/10/2022
    Integra funciones en el intervalo dado con Newton-Cotes
 Divide el intervalo en p funciones polinomiales (Bool's rule)
***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double newtonCotes(double a, double b, double error);
double func(double x){
    //if(fabs(x)<1.0e-5){return(0.0);}                 //Para evitar el cero en algunas funciones de interes
    return((1.0/x)*log(x*x)*sin(4.0*x)*sin(3.0*x));   //-0.1732955599801922   lix = 3.0, ldx = 10.0;
    //return(4.2*x*cos(8.0*x)*exp(-x*x));             //0.1324588139019    lix = -0.6, ldx = 1.4;
    //return(x*exp(2.0*x));                         //7.565747863545      lix = -0.6, ldx = 1.4;
    //return(2*x*x*x + x*x/10.0 - 6.0*x + 12.0);   //21.154666...          lix = -0.6, ldx = 1.4;
}

int main(void){
    clock_t t = clock();
    double tim,intg, realval = -0.1732955599801922;
    double lix = 3.0, ldx = 10.0;

    intg = newtonCotes(lix,ldx,1.0e-8);
     
    printf("La integral de la funcion entre %.1lf y %.1lf aproxima a: \nf(x) = %.11lf\n",lix,ldx,intg);
    printf("El error estimado es: %.3E\n",fabs(intg-realval));
    
    t = clock()-t;
    tim = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\n",tim);
    return 0;
}

/**************************************************************************/

double newtonCotes(double a, double b, double error){
    int m=1,j=0,itmax = 30;
    double ap,bp,st,S=0.0,h,t1 = 0.0,t2 = 1e8;

    while(j<itmax){
        if (fabs(t1-t2)<error){
            printf("Convergencia con %d iteraciones y %d particiones\n\n",j,m);
            break;
        }
        S = 0.0;
        t2 = t1;
        //--------------------------------------------------------------------
        st = fabs((b-a)/(double)m);
        h = st/4.0;
        for (int i = 0; i < m; ++i){
            ap = a + i*st;
            bp = a + (i+1)*st;
            S += 2.0*h*(7.0*func(ap) + 32.0*func(ap+h) + 12.0*func(ap+2.0*h) + 
                32.0*func(ap+3.0*h) + 7.0*func(bp))/45.0;
        }
        t1 = S;
        //--------------------------------------------------------------------
        m*=2;
        j+=1;
        printf("Int= %.11lf...\n",S);
    }
    return(S);
}