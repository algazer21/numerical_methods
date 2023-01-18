/*********************************************
Alan García Zermeño
2/8/2022
Resuelve Ax=b por el método de Cholesky.
*********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


double utv(double x){
    return (x*x+x);
}

void normae(int n,double (*x1)[1]){
    int i;
    double sum = 0,sum2=0,p,t;
    p= 1.0/(n+1);
    for (i = 0; i <n; i++)
    {   
        t = utv(p*i);
        sum+= (x1[i][0]-t)*(x1[i][0]-t);
        sum2+= t*t;
    }
    printf("Error calculado con norma euclidiana: %.9lf\n",sqrt(sum)/sqrt(sum2));
}

void readmatrix(size_t rows, size_t cols, double (*a)[cols], const char* filename, int flag)
{

    int n,m;
    FILE *pf;
    pf = fopen (filename, "r");
    if(flag==1){
        fscanf(pf, "%d %d", &n,&m);
    }

    for(size_t i = 0; i < rows; ++i)
    {
        for(size_t j = 0; j < cols; ++j)
            fscanf(pf, "%lf", a[i] + j);
    }
    fclose(pf);
}

void transpuesta(size_t rows, double (*Ax)[rows]){
	int i,j;
	for (int i = 0; i <rows;i++)
	{
		for (int j = i+1; j < rows; j++)
		{
			Ax[i][j] = Ax[j][i];
			Ax[j][i] = 0.0;
		}
	}
}
void sol_tri_sup(size_t rows, size_t cols, double (*bx)[1], double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j;
    double su=0;
    xx[rows-1][0] = bx[rows-1][0]/Ax[rows-1][rows-1];
    
    for(int i =rows-2;i>=0;i--){
        //compute sum
        for (int j=i+1;j<rows;j++){
            su += Ax[i][j]*xx[j][0];
        }
        xx[i][0] = (bx[i][0]-su)/Ax[i][i];
        su = 0.0;
    }
}

void sol_tri_inf(size_t rows, size_t cols, double (*bx)[1], double (*xx)[1] ,double (*Ax)[cols])
{
    int i,j;
    double su=0;
    xx[0][0] = bx[0][0]/Ax[0][0];
    for(int i = 1;i<rows;i++){
        for (int j=0;j<i;j++){
            su += Ax[i][j]*xx[j][0];
        }
        xx[i][0] = (bx[i][0]-su)/Ax[i][i];
        su = 0.0;
    }
}
/********************************************************************/
void cholesky(size_t rows,double (*Lx)[rows],double (*Ax)[rows])
{
    int i,j,k;
    double sum=0;
    Lx[0][0] = sqrt(Ax[0][0]);
    for(j=1;j<rows;j++){
    	Lx[j][0] = Ax[j][0]/Lx[0][0];
    }
    for(i=1;i<rows-1;i++){
    	for(k=0;k<i;k++){
			sum+= Lx[i][k]*Lx[i][k];
		}
		
		Lx[i][i] = sqrt(Ax[i][i]-sum);
		sum = 0;
		for(j=i+1;j<rows;j++){
			for(k=0;k<i;k++){
				sum+= Lx[j][k]*Lx[i][k];
			}
			Lx[j][i] = (Ax[j][i]-sum)/Lx[i][i];
			sum = 0;
		}
    }
    for(k=0;k<rows;k++){
		sum+= Lx[rows-1][k]*Lx[rows-1][k];
	}
	Lx[rows-1][rows-1] = sqrt(Ax[rows-1][rows-1]-sum);  
}
/********************************************************************/
void writematrix(size_t rows,size_t cols,double (*a)[cols],const char* filename)
{
    int i,j;
    FILE    *fptr;
    fptr = fopen(filename, "w");

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            fprintf(fptr,"%lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}
/*********************************************************************/

int main()
{
    int n,m,i;
    clock_t t = clock();
    double time;
    FILE    *fptr;
    fptr = fopen("AA.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    double (*A)[m] = calloc(n, sizeof *A);
    double (*L)[m] = calloc(n, sizeof *L);
    double (*b)[1] = calloc(n, sizeof *b);
    double (*x)[1] = calloc(n, sizeof *x);
    double (*y)[1] = calloc(n, sizeof *y);

    readmatrix(n, m, A, "AA.txt",1);
    readmatrix(n, 1, b, "bb.txt",0);
    cholesky(n,L,A);  
    //writematrix(n,m,L,"L.txt");
    sol_tri_inf(n,m,b,y,L);
    transpuesta(n,L);
    sol_tri_sup(n,m,y,x,L);
    //normae(n,x);
    writematrix(n,1,x,"x.txt");
 
    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}
