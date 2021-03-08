#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 

//Funzione che computa il prodotto matrice-vettore
void prod_mat_vett(double *a, double *v, int ROWS, int COLS, double *w){
    int i, j;
    
    for(i=0;i<ROWS;i++){
        w[i]=0;
        for(j=0;j<COLS;j++){ 
            w[i] += a[i*COLS+j]* v[j]; 
        }
    }    
}


int main(int argc, char **argv) {


int nproc;			// Numero di processi totale            	
int me;        			// Il mio id         	

double *A;			// Matrice Originale m*n
int m, n;  			// Dimensione della matrice A originale
               	
double *AT;			// Matrice Trasposta n*m
int mT, nT;			// Dimensione della matrice AT trasposta 		

double *v;			// Vettore da moltiplicare di dimensione n
double *w;			// Vettore risultato di dimensione m

double *local_A, *local_AT, *local_v, *local_w;	// Variabili locali a ciascun processore

int local_m, local_n;		// Dimensione dei dati locali per la matrice local_A       
int local_mT, local_nT;  	// Dimensione dei dati locali per la matrice local_AT    

int i, j;                	// Iteratori vari     	

double T_inizio, T_fine, T_max; //Variabili per la misurazione del tempo di esecuzione


/*Attiva MPI*/
MPI_Init(&argc, &argv);

/*Trova il numero totale dei processi*/
MPI_Comm_size (MPI_COMM_WORLD, &nproc);
/*Da ad ogni processo il proprio numero identificativo*/
MPI_Comm_rank (MPI_COMM_WORLD, &me);

// Se sono il processore root
if(me == 0){
    printf("Inserire numero di righe m: \n"); 
    scanf("%d",&m); 
    
    printf("Inserire numero di colonne n: \n"); 
    scanf("%d",&n);

    //Numero di righe e colonne della matrice Trasposta AT
    mT=n;
    nT=m;

    // Numero di righe e di colonne della matrice locale local_A da inviare a ciascun processore 
    local_m = m/nproc;
    local_n = n/nproc;
    
    // Allocazione spazio di memoria
    A = malloc(m*n*sizeof(double));
    AT = malloc(mT*nT*sizeof(double));
    v = malloc(n*sizeof(double));
    w =  malloc(m*sizeof(double)); 
    
    //Inizializzazione e stampa della matrice originale A con m RIGHE e n COLONNE
    printf("\nA = \n"); 
    for (i=0;i<m;i++){
        for(j=0;j<n;j++){
            if (j==0){
                A[i*n+j]= 1.0/(i+1)-1;
            }else{
                A[i*n+j]= 1.0/(i+1)-pow(1.0/2.0,j); 
            }
		printf("%f ", A[i*n+j]);
        }
        printf("\n");
    }

    //Calcolo e stampa della matrice trasposta AT con mT=n RIGHE e nT=m COLONNE
    printf("\nAT = \n"); 
    for (i=0;i<mT;i++){
        for(j=0;j<nT;j++){
            AT[i*nT+j]= A[j*n+i];
            printf("%f ", AT[i*nT+j]);
        }
        printf("\n");
    }

    //Inizializzazione e stampa del vettore v
    printf("\nv = \n");
    for (j=0;j<n;j++){
    	v[j]=j; 
    	printf("%f ", v[j]);
    }
    printf("\n");
    
    
     
} // fine me==0

// Spedisco m, n, local_m, local_n
MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);  
MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);          
MPI_Bcast(&local_m, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);                

//Il numero di righe che ciascuna matrice Trasposta locale di un processore riceve e` pari a local_n.
local_mT = local_n;

//Il numero di colonne che ciascuna matrice Trasposta locale di un processore riceve e` pari a m.
local_nT = m;

// Allocazione del vettore local_v di dimensione (n/nproc)
local_v = malloc(local_n*sizeof(double));

// Distribuzione del vettore v a ciascun processore che riceve la propria porzione in local_v 
MPI_Scatter(v, local_n, MPI_DOUBLE, local_v, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);            

//// Alloczione del vettore local_AT di dimensione local_mT*nT
local_AT = malloc(local_mT*local_nT*sizeof(double));

// Alloczione del vettore local_A di dimensione m*local_n
local_A = malloc(m*local_n*sizeof(double));

// Allocazione del vettore dei risultati parziali local_w = r
local_w = malloc(m*sizeof(double));

// Il processore 0 invia a tutti gli altri un blocco di colonne della matrice AT di dimensione m*(n/nproc)
MPI_Scatter(AT, m*local_n, MPI_DOUBLE, local_AT, m*local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// Stampa del blocco della matrice trasposta local_AT che ha local_m Righe e n Colonne
printf("\nlocal_AT [%d] = \n", me); 
for(i=0; i<local_mT; i++){
    for(j=0; j<local_nT; j++){
        printf("%f\t", local_AT[i*local_nT+j]);
    }
    printf("\n");
}

//Calcolo e stampa del blocco della matrice local_A che ha m Righe e local_n Colonne
printf("\nlocal_A [%d] = \n", me); 
for (i=0; i<m; i++){
    for(j=0; j<local_n; j++){
        local_A[i*local_n+j] = local_AT[j*local_nT+i];
        printf("%f ", local_A[i*local_n+j]);
    }
    printf("\n");
}

MPI_Barrier(MPI_COMM_WORLD);
T_inizio=MPI_Wtime();

// Computazione del prodotto local_A*local_v raccolto nel vettore local_w
prod_mat_vett(local_A, local_v, m, local_n, local_w);

/* raccoglie i vettori dei risultati parziali local_w, li somma, ed invia al processore 0*/
MPI_Reduce(local_w, w, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);
T_fine=MPI_Wtime()-T_inizio;

/* raccoglie il tempo di esecuzione Massimo e lo invia al processore 0*/
MPI_Reduce(&T_fine, &T_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


// 0 stampa la soluzione
if(me==0){ 

    printf("\nw = \n"); 
    for(i = 0; i < m; i++){
        printf("%f ", w[i]);
    }
    printf("\n");
   
    printf("\nTempo calcolo locale: %lf\n", T_fine);
    printf("\nMPI_Reduce max time: %f\n",T_max);

    free(A);
    free(AT);
    free(v);
    free(w);
}

free(local_A);
free(local_AT);
free(local_v);
free(local_w);

MPI_Finalize (); /* Disattiva MPI */
return 0;  
}
