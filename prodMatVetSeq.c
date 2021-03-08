#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 

//Funzione che computa il prodotto matrice-vettore
void prod_mat_vett(double *a, double *v, int ROWS, int COLS, double *w)
{
    int i, j;

    for(i=0;i<ROWS;i++)
    {
        w[i]=0;
        for(j=0;j<COLS;j++)
        { 
            w[i] += a[i*COLS+j]* v[j];
        } 
    }
    
    return;   
}


int main(int argc, char **argv) {

int m,n;                  // Dimensione della matrice
int i,j;                  // Iteratori vari 
double *A,*v,*w;	  //Matrice A, vettore x, vettore y prodotto
double execTime;	  

double T_inizio,T_fine;	  //Variabili per la misurazione del tempo di esecuzione della computazione

printf("Inserire numero di righe m: \n"); 
scanf("%d",&m); 
    
printf("Inserire numero di colonne n: \n"); 
scanf("%d",&n);

// Allocazione spazio di memoria
A = malloc(m*n*sizeof(double));
v = malloc(n*sizeof(double));
w =  malloc(m*sizeof(double));

//Inizializzazione e stampa della Matrice A
printf("\nA = \n"); 
for (i=0;i<m;i++){
	for(j=0;j<n;j++){
        	if (j==0){
        		A[i*n+j]= 1.0/(i+1)-1;
		}else{
                	A[i*n+j]= 1.0/(i+1)-pow(1.0/2.0,j); 
            	}
		printf("%f ", A[i*n+j] );
        }
	printf("\n");
}

//Inizializzazione e stampa del vettore x
printf("\nv = \n");     
	for (j=0;j<n;j++){
    		v[j]=j; 
    		printf("%f ", v[j]);
	}
printf("\n");

MPI_Init(&argc, &argv);

T_inizio = MPI_Wtime();

// Computazione
prod_mat_vett(A,v,m,n,w);
T_fine = MPI_Wtime() - T_inizio;

MPI_Finalize();

//Stampa del vettore y, risultato del prodotto matrice-vettore.
printf("\nw = \n"); 
for(i = 0; i < m; i++){
	printf("%f ", w[i]);
}
printf("\n");

printf("\nTempo di esecuzione in sequenziale: %lf\n", T_fine);

return 0;

}
