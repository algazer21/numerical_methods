/**************************************************************
                    Alan García Zermeño
                        20/10/2022
   Resuelve numericamente integrales dobles usando Montecarlo
   Hay que ingresar manualmente el valor del max y min de la 
   funcion en el intervalo como hM y hm respectivamente.
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void montecarlo3D(double hm, double hM, double lix, double ldx, double liy, double ldy,int pr, long int *m);
double func(double x, double y){
    return(x*x*y*y*y);    //-585/2   lix = -3.0, ldx = 3.0, liy = -3.0, ldy = 2.0;   hm = -260.0, hM = 85.0;
    //return(y*exp(2*x));    // 3*(exp(8)-1)/(4*exp(6))   lix = -3.0, ldx = 1.0, liy = -1.0, ldy = 2.0; hm = -20.0, hM = 20.0;
}

int main(void){
	int n,i,k = 400000;
    long int m[2] = {0,0};   //Contadores del area superior e inferior
    clock_t t = clock();
    double ma=0.0,mi=0.0,tim,s,sm,area,aream;
    double hm = -260.0, hM = 85.0;
    double lix = -3.0, ldx = 3.0, liy = -3.0, ldy = 2.0;
 
    srand((unsigned int)time(NULL));
    area = fabs(ldx-lix)*hM*(ldy-liy);                    //Volumen superior: rectangulo (ld-li)*max
    aream = fabs((ldx-lix)*hm*(ldy-liy));                 //Volumen inferior: rectangulo (ld-li)*min

    for (i = 0; i < 8; ++i){
        montecarlo3D(hm,hM,lix,ldx,liy,ldy,k,m);
        s = ((double)m[0]/(double)k)*area;
        sm = ((double)m[1]/(double)k)*aream;
        ma += s;
        mi += sm;
    }
    ma = ma/8.0;
    mi = mi/8.0;

    printf("La integral de la funcion entre %.1lf y %.1lf en x y entre %.1lf y %.1lf en y aproxima a: \nf(x,y) = %lf\n",lix,ldx,liy,ldy,ma-mi);
    printf("El error absoluto es: %lf\n",fabs((ma-mi)+585/2));
    
    t = clock()-t;
    tim = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\n",tim);
    return 0;
}

/**************************************************************************/

void montecarlo3D(double hm, double hM, double lix, double ldx, double liy, double ldy, int pr, long int *m){
    double t,zt, rangx = fabs(ldx-lix), rangy = fabs(ldy-liy);
    m[0] = 0;
    m[1] = 0;
    for (int i = 0; i < pr; ++i){
        t = func(((double)rand()/(double)(RAND_MAX) * rangx) + lix, ((double)rand()/(double)(RAND_MAX) * rangy) + liy);
        zt = ((double)rand()/(double)(RAND_MAX)) * hM;
        if (zt<=t && zt >0){
            m[0]+=1;
        }
    }

    if(hm < 0){
        for (int i = 0; i < pr; ++i){
            t = func(((double)rand()/(double)(RAND_MAX) * rangx) + lix, ((double)rand()/(double)(RAND_MAX) * rangy) + liy);
            zt = ((double)rand()/(double)(RAND_MAX)) * hm;
            if (zt>=t && zt < 0){
                m[1]+=1;
            }
        }
    }
}


