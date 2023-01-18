/**************************************************************
                    Alan García Zermeño
                         20/9/2022
Aplica Factorizacion QR para resolver sistemas Ax=b De forma que 
                         QRx = b
                      Q^tQRx = Q^tb
                         Rx = b*
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void readmatrix(size_t rows, size_t cols, double **a, const char* filename);
double error(size_t n,double **a,double *x,double *b);
void prod_matr_matr(int n, double **matrizA, double **matrizB, double **producto);
void prod_matr_vec(double **a,double *v0,double *v1, size_t rows);
void readvector(size_t rows, double *a, const char* filename);
void writematrix(size_t rows,size_t cols,double **a, const char* filename);
void writevector(size_t rows,double *a,const char* filename);
void sol_tri_sup(size_t rows, double *bx, double *xx ,double **Ax);
void qr(size_t n, double **a, double **q, double **r);
void extractCOL(int n,double **a,double *aa, int col);
void transpuesta(size_t rows, double **Ax,double **a);
double prod_vec_vec(double *b,double *c, size_t rows);
double normalize(size_t rows, double *x);
void freel(int n,double **listm);


int main(void){
	int n,m,i;
    clock_t t = clock();
    double time,er;
    FILE    *fptr;
    fptr = fopen("AA.txt", "r");
    fscanf(fptr, "%d %d", &n,&m);
    fclose(fptr);
    double **A = (double **)malloc(n*sizeof(double*));
    double **Q = (double **)malloc(n*sizeof(double*));
    double **Qt = (double **)malloc(n*sizeof(double*));
    double **R = (double **)malloc(n*sizeof(double*));
    double **aux = (double **)malloc(n*sizeof(double*));
    double *b = (double *)malloc(n*sizeof(double));
    double *bb = (double *)malloc(n*sizeof(double));
    double *x = (double *)malloc(n*sizeof(double));
    for(i = 0; i < n; i++){
        Q[i] = (double *)malloc(n*sizeof(double));
        Qt[i] = (double *)malloc(n*sizeof(double));
        aux[i] = (double *)malloc(n*sizeof(double));
        R[i] = (double *)malloc(n*sizeof(double));}

    readmatrix(n, m, A, "AA.txt");
    readvector(n, b, "bb.txt");
    qr(n,A,Q,R);
    //writematrix(n,m,R,"R.txt");
    transpuesta(n,Qt,Q);
    prod_matr_vec(Qt,b,bb,n);       //Q^t b = bb
    sol_tri_sup(n,bb,x,R);           //Rx = bb
    //prod_matr_matr(n,Qt,Q,aux);
    //writematrix(n,m,Qt,"Q.txt");
    writevector(n,x,"x.txt");
    er = error(n,A,x,b);
    printf("\nError calculado con norma euclidiana: %.20lf\n",er);

    freel(n,A); freel(n,Q); freel(n,R); freel(n,Qt); freel(n,aux);
    free(b); free(bb); free(x);

    t = clock()-t;
    time = ((double)t)/CLOCKS_PER_SEC;
    printf("Tiempo: %lf segundos\nSoluciones guardadas en x.txt\n",time);
    return 0;
}

/**************************************************************************/

void readmatrix(size_t rows, size_t cols, double **a, const char* filename){
    int n,m;
    FILE *pf;
    pf = fopen (filename, "r");
    if (pf==NULL){
        printf("No se encontro archivo\n");
        exit(1);}
    fscanf(pf, "%d %d", &n,&m);
    for(size_t i = 0; i < rows; ++i){
        a[i] = (double *)malloc(cols*sizeof(double));
        for(size_t j = 0; j < cols; ++j){
            fscanf(pf, "%lf", a[i] + j);
        }    
    }
    fclose(pf);
}

void freel(int n,double **listm){
    int i;
    for(i=0;i<n;i++){
        free(listm[i]);
    }
    free(listm);
}

void readvector(size_t rows, double *a, const char* filename){
	int n=rows,m=1;
    FILE *pf;
    pf = fopen (filename, "r");
    //fscanf(pf, "%d %d", &n,&m);
    if (pf==NULL){
        printf("No se encontro archivo\n");
        exit(1);}
    for(size_t i = 0; i < rows; ++i){
            fscanf(pf, "%lf", a + i);
    }
    fclose(pf);
}


void writematrix(size_t rows,size_t cols,double **a,const char* filename)
{
    int i,j;
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            fprintf(fptr,"%.2lf ",a[i][j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}

/********************************QR************************************/
void qr(size_t n, double **a, double **q, double **r){
	unsigned int i,j,k,kk;
	double nor,dot;
	double *aa = (double *)malloc(n*sizeof(double));
	double *aa2 = (double *)malloc(n*sizeof(double));
	double *res = (double *)malloc(n*sizeof(double));

	for(j=0;j<n;j++){
		extractCOL(n,a,aa2,j);	               //Extrae la columna j de A y la guarda en aa2
		for(k=0;k<n;k++){res[k] = aa2[k];}

		for(kk=0;kk<j;kk++){
			extractCOL(n,q,aa,kk);
			dot = prod_vec_vec(aa2,aa,n);      //Producto interno entre la aa y aa2
			r[kk][j] = dot;
			for(i=0;i<n;i++){
				res[i] -= dot*aa[i];
			}
		}
		nor = normalize(n,res);                //Normaliza res y devuelve la norma
		r[j][j] = nor;

		for(i=0;i<n;i++){
			q[i][j] = res[i];
		}
	}
	free(aa); free(aa2); free(res);
}
/*********************************QR********************************/
double normalize(size_t rows, double *x){
    double sum=0;
    for (int i=0;i<rows;i++){
        sum+= x[i]*x[i];}
    sum = sqrt(sum);
    for (int i=0;i<rows;i++){
        x[i] = x[i]/sum;}
    return(sum);
}

void extractCOL(int n,double **a,double *aa, int col){
	for(int i=0;i<n;i++){aa[i] = a[i][col];}
}

double prod_vec_vec(double *b,double *c, size_t rows){
    double a= 0;
    for (int i = 0;i<rows;i++){
        a += b[i]*c[i];
    }
    return a;
}

void transpuesta(size_t rows, double **Ax,double **a){
	int i,j;
	for (i = 0; i <rows;i++){
		for (j = 0; j < rows; j++){
			Ax[i][j] = a[j][i];
		}
	}
}

void prod_matr_matr(int n, double **matrizA, double **matrizB, double **producto){
	double suma = 0;
	for (int a = 0; a < n; a++) {
    	for (int i = 0; i < n; i++) {
        	suma = 0;
        	for (int j = 0; j < n; j++) {
            	suma += matrizA[i][j] * matrizB[j][a];
        	}
        	producto[i][a] = suma;
    	}
	}
}

void prod_matr_vec(double **a,double *v0,double *v1, size_t rows){
    for (int i = 0;i<rows;i++){
        v1[i] = prod_vec_vec(a[i],v0,rows);
    }
}

void sol_tri_sup(size_t rows, double *bx, double *xx ,double **Ax)
{
    int i,j;
    double su=0;
    for(i=(rows-1);i>=0;i--){
        for(j=i+1;j<rows;j++){
            su+=Ax[i][j]*xx[j];
        }
        xx[i] = (bx[i]-su)/Ax[i][i];
        su = 0;
    }
}

void writevector(size_t rows,double *a,const char* filename)
{
    FILE    *fptr;
    fptr = fopen(filename, "w");
    if (fptr==NULL){
        printf("No fue posible escribir %s\n",filename);
        exit(1);}
    printf("Soluciones al sistema:\n\n");
    for(int i=0;i<rows;i++){
            fprintf(fptr,"%.7lf\n",a[i]);
            printf("%.7lf\n",a[i]);
        }
    fclose(fptr);
}

double error(size_t n,double **a,double *x,double *b){   //||Ax - b||
	double sum=0.0;
	double *bx = (double *)malloc(n*sizeof(double));
	prod_matr_vec(a,x,bx,n);
	for (int i = 0;i <n;i++){
		sum+=(bx[i]-b[i])*(bx[i]-b[i]);
	}
	free(bx);
	return(sqrt(sum));
}
