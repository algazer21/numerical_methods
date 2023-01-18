/**************************************************************
                    Alan García Zermeño
             Para el curso de métodos numéricos.
                     CIMAT 17/10/2022
             Integra funciones usando Montecarlo
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void montecarlo2D(double *h, double li, double ld, int g, int pr, long int *m);
void maxi(double stp, double li, double ld, double *h);
void savepoints(int m, int *n, double *x);


double func(double x){
    //return(x*x*cos(4.0*x)); //-1.7481518546
    return(3.0*exp(-3.0*x*x));   //3.0699801238
}

int main(void){
	int n,i,k = 210000, M = 100, Mm = 1000,N[100];
    long int m[2] = {0,0};
    clock_t t = clock();
    double ma=0.0,mi=0.0,tim,s,sm,area,aream,h[2], li = -3.0, ld = 3.0, XX[100];
    double realval = 3.0699801238;

    maxi((double)1000,li,ld,h);                 //Aproxima maximos y minimos de la funcion en el intervalo [li,ld]
    srand((unsigned int)time(NULL));
    area = fabs(ld-li)*h[0];                    //Area superior: rectangulo (ld-li)*max
    aream = fabs((ld-li)*h[1]);                 //Area inferior: rectangulo (ld-li)*min

    for(int j=1;j<= M;j++){
    for (i = 0; i < 4; ++i){
        montecarlo2D(h,li,ld,i,Mm*j,m);
        s = ((double)m[0]/(double)(Mm*j))*area;
        sm = ((double)m[1]/(double)(Mm*j))*aream;
        ma += s;
        mi += sm;
    }
    ma = ma/4.0;
    mi = mi/4.0;

    printf("La integral de la funcion entre %.1lf y %.1lf aproxima a: \nf(x) = %lf\n",li,ld,ma-mi);
    printf("El error absoluto es: %lf con %d numeros aleatorios\n",fabs((ma-mi)-realval),Mm*j);

    N[j-1] = Mm*j;
    XX[j-1] = ma-mi-realval;
    mi = 0.0;
    ma = 0.0;
}   
    savepoints(M,N,XX);
    
    t = clock()-t;
    tim = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\n",tim);
    return 0;
}

/**************************************************************************/

void montecarlo2D(double *h, double li, double ld, int g, int pr, long int *m){
    double t,yt, rang = fabs(ld-li);
    m[0] = 0;
    m[1] = 0;
    for (int i = 0; i < pr; ++i){
        t = func(((double)rand()/(double)(RAND_MAX) * rang) + li);
        yt = ((double)rand()/(double)(RAND_MAX)) * h[0];
        if (yt<=t && yt >0){
            m[0]+=1;
        }
    }

    if(h[1] < 0){
        for (int i = 0; i < pr; ++i){
            t = func(((double)rand()/(double)(RAND_MAX) * rang) + li);
            yt = ((double)rand()/(double)(RAND_MAX)) * h[1];
            if (yt>=t && yt < 0){
                m[1]+=1;
            }
        }
    }
}

void maxi(double stp, double li, double ld, double *h){
    double min = 0.0,max = 0.0,f, ft = fabs(ld-li)/stp;
    int gg = (int)(stp);

    for (int i = 0; i < gg; ++i){
        f = func(li + ft*i);
        if(f>max){max = f;}
        if(f<min){min = f;}
    }
    f = func(ld);
    if(f>max){max = f;}
    if(f<min){min = f;}

    max = max + max*(1.0/10.0);
    min = min + min*(1.0/10.0);
    h[0] = max;
    h[1] = min;
}


void savepoints(int m, int *n, double *x){
   FILE *pf;
   pf = fopen("error.dat", "w");
   if (pf==NULL){
      printf("No se encontro archivo\n");
      exit(1);}

   for (int i = 0; i < m; ++i){
      fprintf(pf,"%d %lf\n",n[i],x[i]);
   }
   fclose(pf);
   printf("%d Puntos generados guardados exitosamente en 'datos.dat'\n",m);
}
